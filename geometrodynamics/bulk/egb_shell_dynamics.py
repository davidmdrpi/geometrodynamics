"""Symmetric moving EGB shell in a supported, ultrastatic round S4 bulk.

Signature -++++; stresses below are multiplied by kappa5=8*pi*G5.
Normals point from the shell into each retained exterior. Fixed bulk and affine
shell response only: not a coupled perturbation evolution or a shell action.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import solve_ivp


def junction_tensor(k0, q, intrinsic_acceleration, intrinsic_spatial, alpha):
    """Independent orthonormal tensor contraction of Davis hep-th/0208205 (10).

    Inputs: mixed K eigenvalues (k0,q,q,q), R_0i0i=-intrinsic_acceleration,
    R_ijij=intrinsic_spatial. Return kappa5 S_ab with both indices down.
    This routine does not call the reduced shell stress formula.
    """
    h = np.diag([-1., 1., 1., 1.])
    K = np.diag([k0, q, q, q])
    trace = np.trace(K)
    Jmixed = (2*trace*K@K + np.trace(K@K)*K - 2*K@K@K - trace**2*K)/3
    J = h@Jmixed
    riemann = np.zeros((4, 4, 4, 4))
    for a in range(4):
        for b in range(a+1, 4):
            val = -intrinsic_acceleration if a == 0 else intrinsic_spatial
            riemann[a,b,a,b] = riemann[b,a,b,a] = val
            riemann[a,b,b,a] = riemann[b,a,a,b] = -val
    ricci = np.einsum('ac,abcd->bd', h, riemann)
    scalar = np.einsum('ab,ab', h, ricci)
    P = np.empty_like(riemann)
    for a,b,c,d in np.ndindex(P.shape):
        P[a,b,c,d] = (riemann[a,b,c,d] + ricci[b,c]*h[d,a]
                     - ricci[b,d]*h[c,a] - ricci[a,c]*h[d,b]
                     + ricci[a,d]*h[c,b]
                     + scalar/2*(h[a,c]*h[d,b]-h[a,d]*h[c,b]))
    PK = np.einsum('acdb,cd->ab', P, K@h)
    return -2*(h@K-trace*h + 2*alpha*(3*J-np.trace(Jmixed)*h+2*PK))


@dataclass(frozen=True)
class RoundEGBShell:
    radius: float = 1.
    mouth_angle: float = .3
    alpha: float = -.06
    eta: float = .1

    def __post_init__(self):
        if (not all(math.isfinite(x) for x in (self.radius, self.mouth_angle, self.alpha, self.eta))
                or self.radius <= 0 or not 0 < self.mouth_angle < math.pi/2):
            raise ValueError('require finite parameters, R>0 and 0<a<pi/2')

    @property
    def r0(self):
        return self.radius*math.sin(self.mouth_angle)

    @property
    def bulk_C(self):
        return 1+4*self.alpha/self.radius**2

    def geometry(self, r, v):
        if not np.isfinite(r) or not np.isfinite(v) or not 0 < r < self.radius:
            raise ValueError('shell outside open radial patch 0<r<R')
        f = 1-r*r/self.radius**2
        delta = np.sqrt(f+v*v)
        return f, -2*r/self.radius**2, delta, delta/r

    def density(self, r, v):
        _, _, _, q = self.geometry(r, v)
        B = (1+v*v)/r**2
        return -6*q*(1+4*self.alpha*(B-q*q/3))

    def acceleration_coefficient(self, r, v):
        _, _, delta, _ = self.geometry(r, v)
        return 2*self.bulk_C/delta+16*self.alpha*delta/r**2

    def pressure(self, r, v, acceleration):
        f, fp, delta, q = self.geometry(r, v)
        return (4*q+self.acceleration_coefficient(r,v)*acceleration
                -self.bulk_C*fp*v*v/(f*delta))

    @property
    def equilibrium(self):
        return self.density(self.r0, 0.), self.pressure(self.r0, 0., 0.)

    def material_pressure(self, r, v):
        s0, p0 = self.equilibrium
        return p0+self.eta*(self.density(r,v)-s0)

    def acceleration(self, r, v):
        f, fp, delta, q = self.geometry(r, v)
        D = self.acceleration_coefficient(r,v)
        if abs(D) < 1e-12:
            raise ValueError('degenerate pressure--acceleration equation')
        return (self.material_pressure(r,v)-4*q+self.bulk_C*fp*v*v/(f*delta))/D

    def flux(self, r, v):
        """2*kappa5*(rho_bulk+p_bulk)*gamma*w, into shell convention."""
        f, _, delta, _ = self.geometry(r, v)
        return 6*self.bulk_C*v*delta/(self.radius**2*f)

    def linear_data(self):
        r = self.r0
        E = 1-4*self.alpha/self.radius**2+8*self.alpha/r**2
        if abs(E) < 1e-12:
            raise ValueError('degenerate equilibrium')
        return dict(E=E, D=self.acceleration_coefficient(r,0.),
                    omega2=-(3*self.eta+2/E)/r**2,
                    eta_critical=-2/(3*E), bulk_C=self.bulk_C,
                    bulk_tensor_speed2=1/self.bulk_C if self.bulk_C != 0 else None)

    def static_bounds(self):
        upper = -3*self.r0**2/(4*(2+math.sin(self.mouth_angle)**2))
        return dict(alpha_lower=-self.radius**2/4, alpha_WEC=upper,
                    alpha_DEC=5*upper/3,
                    minimum_GB_length_over_mouth=math.sqrt(-upper)/self.r0)

    def evolve(self, fraction, duration, *, rtol=2e-10, samples=4001):
        """Proper-time evolution plus pressure work and exterior-flux integrals.

        Abort before leaving radial patch or D=0. A pole is not regularized.
        Ledger uses unit S3 area; multiply all energies by 2*pi^2/kappa5.
        """
        if not math.isfinite(duration) or duration <= 0 or samples < 3:
            raise ValueError('invalid integration horizon or sample count')
        ri = self.r0*(1+fraction)
        self.geometry(ri, 0.)
        def rhs(t, y):
            r,v,_,_ = y
            return [v, self.acceleration(r,v),
                    self.material_pressure(r,v)*3*r*r*v, r**3*self.flux(r,v)]
        def patch(t,y):
            return min(y[0]/self.radius-1e-5, 1-1e-5-y[0]/self.radius)
        def degeneracy(t,y):
            return abs(self.acceleration_coefficient(y[0],y[1]))-1e-6
        patch.terminal = degeneracy.terminal = True
        patch.direction = degeneracy.direction = -1
        sol = solve_ivp(rhs, (0.,duration), [ri,0.,0.,0.], method='DOP853',
                        t_eval=np.linspace(0.,duration,samples), rtol=rtol,
                        atol=rtol*1e-3, max_step=duration/(samples-1),
                        events=[patch,degeneracy])
        if not sol.success or sol.status != 0 or sol.t[-1] != duration:
            raise RuntimeError('unresolved integration or admission-boundary stop: '+sol.message)
        return sol.t, sol.y


def casimir_scale(D, species, loop_ratio=1.):
    """b/lP for N/(b^(D-2)*L^2) ~ 1/(G_D*b^2), L=loop_ratio*b.

    Order-one factors omitted; no charge enhancement or near-extremal balance.
    """
    if D < 3 or int(D) != D or species <= 0 or loop_ratio <= 0:
        raise ValueError('require integer D>=3 and positive species/loop ratio')
    return (species/loop_ratio**2)**(1/(D-2))
