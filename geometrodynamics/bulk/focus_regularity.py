"""Regularity does not select a unique finite-mouth event law.

Static 4+1-dimensional *controls*, on the existing finite-mouth spatial
profile. These inverse-engineered metrics are not a dynamical matter model,
a gravitational-wave evolution, or a derivation of quantum mechanics.
All stresses below are kappa_5 times orthonormal stresses, with Lambda=0.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import quad, solve_ivp

from geometrodynamics.bulk.finite_mouth import ThroatGeometry


@dataclass(frozen=True)
class RegularFocus:
    """Same spatial data and mouth clocks, different regular interior lapse.

    N = exp(epsilon (1-(s/S)^2)^3) in the neck. Its C2 extension is N=1
    outside. Restrict epsilon to a declared compact example family; this
    range is a modelling input, not inferred from observations.
    """
    radius: float = 1.0
    mouth_angle: float = 0.5
    epsilon: float = 0.0

    def __post_init__(self):
        if not all(math.isfinite(v) for v in (self.radius, self.mouth_angle, self.epsilon)):
            raise ValueError("parameters must be finite")
        if self.radius <= 0 or not 0 < self.mouth_angle < math.pi / 2:
            raise ValueError("require R>0 and 0<a<pi/2")
        if not 0 <= self.epsilon <= 0.1:
            raise ValueError("declared lapse family is 0<=epsilon<=0.1")

    @property
    def geometry(self):
        return ThroatGeometry(self.radius, self.mouth_angle)

    def fields(self, s):
        """r,r',r'',N,(log N)',(log N)'' at neck points only."""
        s = np.asarray(s, dtype=float)
        g = self.geometry
        if np.any(~np.isfinite(s)) or np.any(abs(s) > g.half_length * (1 + 1e-14)):
            raise ValueError("point must be in the finite neck")
        u = np.clip(s / g.half_length, -1, 1)
        r = np.sqrt(s*s + g.neck_radius**2)
        f = (1-u*u)**3
        fp = -6*u*(1-u*u)**2 / g.half_length
        fpp = (-6+36*u*u-30*u**4) / g.half_length**2
        return r, s/r, g.neck_radius**2/r**3, np.exp(self.epsilon*f), self.epsilon*fp, self.epsilon*fpp

    def invariants(self, s):
        """Exact warped-product curvature, effective stress and expansions.

        k_+/- = N^-1 partial_t +/- partial_s, k_+ . k_- = -2.
        A=r''/r, B=(1-r'^2)/r^2, C=N''/N, D=N'r'/(Nr).
        """
        r, rp, rpp, lapse, lp, lpp = self.fields(s)
        A, B, C, D = rpp/r, (1-rp*rp)/r**2, lpp+lp*lp, lp*rp/r
        rho = 3*(B-A)
        radial = 3*(D-B)
        tangential = C+2*D+2*A-B
        return dict(lapse=lapse, radius=r, scalar=-2*C-6*D-6*A+6*B,
                    kretschmann=4*C*C+12*D*D+12*A*A+12*B*B,
                    kappa_rho=rho, kappa_p_radial=radial, kappa_p_tangential=tangential,
                    kappa_radial_nec=rho+radial,
                    theta_plus=3*rp/r, theta_minus=-3*rp/r)

    def metric(self, coordinates):
        """Local chart (t,s,chi,theta,phi), away from angular chart poles."""
        _, s, chi, theta, _ = coordinates
        r, _, _, lapse, _, _ = self.fields(s)
        return np.diag([-lapse*lapse, 1., r*r, r*r*np.sin(chi)**2,
                        r*r*np.sin(chi)**2*np.sin(theta)**2])

    def transit_time(self):
        """Radial null transit measured in the common N=1 mouth clock."""
        S = self.geometry.half_length
        return quad(lambda s: 1/float(self.fields(s)[3]), -S, S,
                    epsabs=1e-12*self.radius, epsrel=1e-12)[0]

    def transfer(self, omega, ell=0, rtol=1e-11):
        """Scalar control: (phi, N r^3 phi_s / b^2) at the two actual mouths.

        No asymptotic scattering phases, no gravitational-wave interpretation.
        The real first-order generator has zero trace: det M=1 and the
        directed Wronskian flux Im(conj(phi) N r^3 phi_s) is conserved.
        """
        if not math.isfinite(omega) or omega < 0 or not isinstance(ell, int) or ell < 0:
            raise ValueError("require finite omega>=0 and integer ell>=0")
        if not 0 < rtol < 1:
            raise ValueError("require 0<rtol<1")
        b, S = self.geometry.neck_radius, self.geometry.half_length

        def rhs(s, y):
            r, _, _, N, _, _ = self.fields(s)
            generator = np.array([[0., b*b/(N*r**3)],
                                  [(N*r*ell*(ell+2)-omega**2*r**3/N)/(b*b), 0.]])
            return (generator @ y.reshape(2, 2)).ravel()

        sol = solve_ivp(rhs, (-S, S), np.eye(2).ravel(), method='DOP853',
                        rtol=rtol, atol=rtol*0.01)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y[:, -1].reshape(2, 2)


def round_caustic_control(spatial_dimension=4, radius=1.0):
    """A null ray caustic in a curvature-regular ultrastatic round spacetime.

    Its required background stress is explicit; it is not vacuum Einstein.
    The transverse Jacobi factor R sin(lambda/R) vanishes at pi R.
    """
    n = spatial_dimension
    if not isinstance(n, int) or n < 2 or not math.isfinite(radius) or radius <= 0:
        raise ValueError("require integer n>=2 and finite R>0")
    return dict(focal_affine_parameter=math.pi*radius,
                jacobi_factor=radius*math.sin(math.pi),
                kretschmann=2*n*(n-1)/radius**4,
                kappa_rho=n*(n-1)/(2*radius**2),
                kappa_pressure=-(n-1)*(n-2)/(2*radius**2))
