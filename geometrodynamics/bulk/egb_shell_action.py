"""A specified scalar material completion and a kinetic-admission diagnostic.

All matter actions/stresses use kappa5 units, chi=sqrt(kappa5)*phi.
This does not integrate out gravitational constraints or compute a physical
radial eigenmode norm. It tests whether the earlier frozen-bulk mode is even
an admissible tangent to this completion.
"""
from dataclasses import dataclass, field
import math

import numpy as np
from scipy.integrate import quad

from geometrodynamics.bulk.egb_shell_dynamics import RoundEGBShell


@dataclass(frozen=True)
class ScalarShellAction:
    shell: RoundEGBShell = field(default_factory=RoundEGBShell)
    bulk_lambda: float = 0.

    def __post_init__(self):
        if not 0 < self.shell.eta <= 1:
            raise ValueError('surface power law requires 0<eta<=1; eta=0 is not covered')
        s,p=self.shell.equilibrium
        if self.shell.bulk_C <= 0 or s+p <= 0 or not math.isfinite(self.bulk_lambda):
            raise ValueError('require C>0, s+p>0 and finite lambda')

    @property
    def q2(self):
        return 3*self.shell.bulk_C/self.shell.radius**2

    @property
    def X0(self):
        return self.q2/2

    @property
    def potential(self):
        R=self.shell.radius
        return (9/R**2+12*self.shell.alpha/R**4)/2

    @property
    def exponent(self):
        return (1+self.shell.eta)/(2*self.shell.eta)

    def bulk_P(self,X):
        return X-self.potential+self.bulk_lambda*(X-self.X0)**2/(2*self.X0)

    def bulk_PX(self,X):
        return 1+self.bulk_lambda*(X-self.X0)/self.X0

    def bulk_stress(self,X):
        p=self.bulk_P(X)
        return 2*X*self.bulk_PX(X)-p,p

    def surface_P(self,Y):
        if np.any(np.asarray(Y)<=0):
            raise ValueError('surface scalar gradient must be timelike: Y>0')
        s,p=self.shell.equilibrium
        return p+(s+p)/(2*self.exponent)*((Y/self.X0)**self.exponent-1)

    def surface_PY(self,Y):
        if np.any(np.asarray(Y)<=0):
            raise ValueError('require Y>0')
        s,p=self.shell.equilibrium
        return (s+p)/(2*self.X0)*(Y/self.X0)**(self.exponent-1)

    def surface_PYY(self,Y):
        return self.surface_PY(Y)*(self.exponent-1)/Y

    def surface_stress(self,Y):
        p=self.surface_P(Y)
        return 2*Y*self.surface_PY(Y)-p,p

    def bare_coefficients(self):
        """Local material coefficients before gravitational constraint reduction."""
        spatial=self.surface_PY(self.X0)
        temporal=spatial+2*self.X0*self.surface_PYY(self.X0)
        bulk_temporal=1+2*self.bulk_lambda
        return dict(bulk_time=bulk_temporal,bulk_space=1.,
                    bulk_sound2=1/bulk_temporal if bulk_temporal>0 else None,
                    surface_time=temporal,surface_space=spatial,
                    surface_sound2=spatial/temporal)

    def induced_Y(self,r,v):
        """Only for the deliberately frozen chi=q*t, unit-lapse exterior."""
        f,_,_,_=self.shell.geometry(r,v)
        return self.X0*(1+v*v/f)

    def frozen_density_residual(self,r,v):
        s_material=self.surface_stress(self.induced_Y(r,v))[0]
        return self.shell.density(r,v)-s_material

    def frozen_current_residual(self,r,v,acceleration):
        """(r^-3 d_tau[r^3 P_Y chi_dot]) - sum_faces n(chi).

        A nonzero residual violates the natural scalar boundary equation. The
        global charge derivative is r^3 times this residual per unit S3 area.
        """
        f,fp,delta,_=self.shell.geometry(r,v)
        q=math.sqrt(self.q2)
        Y=self.induced_Y(r,v)
        gamma=delta/math.sqrt(f)
        Ydot=self.X0*(2*v*acceleration/f-fp*v**3/f**2)
        Z=self.surface_PY(Y)+2*Y*self.surface_PYY(Y)
        return Z*Ydot/(q*gamma)+3*v/r*self.surface_PY(Y)*q*gamma-2*q*v/math.sqrt(f)

    def linear_obstruction(self):
        m=self.shell;r=m.r0;q=math.sqrt(self.q2);c=math.cos(m.mouth_angle)
        slope=3*m.linear_data()['D']/r**2
        s,p=m.equilibrium
        return dict(geometric_density_slope=slope,
                    frozen_material_density_slope=0.,
                    current_residual_per_velocity=q*(3*self.surface_PY(self.X0)/r-2/c),
                    total_charge_slope_no_S3_area=-r**3*slope/q,
                    required_scalar_trace_rate_per_q_per_displacement=m.eta*slope/(s+p))

    def frozen_total_charge(self,r,v=0.):
        """Shift charge / unit S3 area, for the trial fields, not a solution.

        Spatial bulk interval is chi in [a,pi-a]. There is ONE shell with TWO
        incident bulk faces. q, action constants and R stay fixed as r changes.
        """
        R=self.shell.radius;q=math.sqrt(self.q2)
        f,_,delta,_=self.shell.geometry(r,v)
        a=math.asin(r/R)
        volume=R**4*quad(lambda x:math.sin(x)**3,a,math.pi-a,epsabs=1e-13)[0]
        surface=r**3*self.surface_PY(self.induced_Y(r,v))*q*delta/math.sqrt(f)
        return q*volume+surface

    def material_lagrangian(self,kind,time_derivative,spatial_derivative=0.):
        """Local Lagrangian per invariant volume for independent Hessian tests."""
        invariant=(time_derivative**2-spatial_derivative**2)/2
        if kind=='bulk':return self.bulk_P(invariant)
        if kind=='surface':return self.surface_P(invariant)
        raise ValueError('kind must be bulk or surface')

    def local_shift_stationary_energy(self,kind,time_derivative,spatial_derivative=0.):
        """Fixed-geometry local H-q*J; NOT the constrained gravitational energy.

        chi=q*t is stationary under time translation combined with a scalar
        shift. Subtracting q*J removes the linear background-charge variation.
        """
        invariant=(time_derivative**2-spatial_derivative**2)/2
        if kind=='bulk':P,PX=self.bulk_P(invariant),self.bulk_PX(invariant)
        elif kind=='surface':P,PX=self.surface_P(invariant),self.surface_PY(invariant)
        else:raise ValueError('kind must be bulk or surface')
        return PX*time_derivative*(time_derivative-math.sqrt(self.q2))-P


def second_derivative(fun,x,h):
    return (-fun(x+2*h)+16*fun(x+h)-30*fun(x)+16*fun(x-h)-fun(x-2*h))/(12*h*h)
