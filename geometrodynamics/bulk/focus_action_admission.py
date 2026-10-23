"""Closed ghost-scalar action: admission obstruction and Ellis control.

Signature -++++; chi=sqrt(kappa_5)*phi. The action is
integral sqrt(-g) [R+(grad chi)^2]/(2*kappa_5). This is a negative
control, not an accepted BAM focus law. See docs/focus_action_admission.md.
"""
import math

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eigh_tridiagonal

from geometrodynamics.bulk.focus_regularity import RegularFocus


def scalar_boundary_obstruction(radius=1., mouth_angle=.5):
    """Static radial P(X,chi) requires rho+p_angular=0.

    On the original neck this is (r^2 N')'+b^2 N/r^2=0.
    Positive N and N'=0 at both mouths are incompatible.
    Equal unit boundary lapses instead force N=r_m/r.
    """
    g = RegularFocus(radius, mouth_angle).geometry
    b, S = g.neck_radius, g.half_length
    rm = math.hypot(b, S)
    integral = quad(lambda s: b*b/(s*s+b*b), -S, S, epsabs=1e-13)[0]
    return dict(neck_radius=b, half_length=S, mouth_radius=rm,
                unit_lapse_positive_integral=integral,
                integral_exact=2*b*math.atan(S/b),
                forced_symmetric_lapse_at_waist=rm/b,
                forced_right_mouth_lapse_derivative=-S/rm**2,
                original_waist_rho_plus_p_angular=1/b**2,
                ellis_radius_for_same_round_mouth=radius*math.sin(mouth_angle)**1.5)


class EllisGhost:
    """Smooth two-ended five-dimensional massless-ghost solution.

    The background ODE is solved in x=s/b on [0,extent] and reflected.
    extent is a numerical evaluation domain, not a physical mouth wall.
    """
    def __init__(self, b=1., extent=16.):
        if not math.isfinite(b) or b <= 0 or not math.isfinite(extent) or extent <= 0:
            raise ValueError('require finite positive b and extent')
        self.b, self.extent = b, extent
        self.solution = solve_ivp(lambda x, q: [q[1], 2/q[0]**5, math.sqrt(6)/q[0]**3],
                                  (0, extent), [1., 0., 0.], dense_output=True,
                                  method='DOP853', rtol=5e-14, atol=1e-14, max_step=.05)
        if not self.solution.success:
            raise RuntimeError(self.solution.message)

    def fields(self, s):
        s = np.asarray(s, dtype=float)
        x = s/self.b
        if np.any(~np.isfinite(s)) or np.any(abs(x) > self.extent*(1+1e-14)):
            raise ValueError('outside numerical background domain')
        shape = s.shape
        y, v, chi = self.solution.sol(np.minimum(abs(x), self.extent).reshape(-1))
        y, v, chi = [a.reshape(shape) for a in (y, v, chi)]
        r, rp = self.b*y, np.sign(s)*v
        return r, rp, 2/(self.b*y**5), np.sign(s)*chi, math.sqrt(6)/(self.b*y**3)

    def metric(self, coordinates):
        _, s, chi, theta, _ = coordinates
        r = self.fields(s)[0]
        return np.diag([-1., 1., r*r, r*r*np.sin(chi)**2,
                        r*r*np.sin(chi)**2*np.sin(theta)**2])

    def stress(self, s):
        """kappa_5 times orthonormal stress, directly from scalar action."""
        cp = self.fields(s)[4]
        return -.5*cp**2, -.5*cp**2, .5*cp**2

    def potential(self, x):
        """Dimensionless regular radial master potential; Torii–Shinkai (3.18)."""
        y = self.fields(self.b*np.asarray(x))[0]/self.b
        return -1/(4*y*y)-27/(4*y**6)

    def radial_mode(self, half_domain=12., spacing=.02):
        """Dirichlet approximation to decay at both spatial infinities.

        Returns dimensionless lambda=(b*omega)^2 and residual in the
        Euclidean norm. Discretization is centered second order.
        """
        if not 0 < half_domain <= self.extent or not 0 < spacing < half_domain:
            raise ValueError('invalid spectral domain or spacing')
        intervals = round(2*half_domain/spacing)
        h = 2*half_domain/intervals
        if intervals < 4:
            raise ValueError('need at least four intervals')
        x = np.linspace(-half_domain, half_domain, intervals+1)[1:-1]
        diagonal = 2/h**2+self.potential(x)
        off = np.full(len(x)-1, -1/h**2)
        values, vectors = eigh_tridiagonal(diagonal, off, select='i', select_range=(0, 1))
        u = vectors[:, 0]
        residual = (diagonal-values[0])*u
        residual[:-1] += off*u[1:]
        residual[1:] += off*u[:-1]
        return dict(half_domain_over_b=half_domain, spacing_over_b=h,
                    lambda0=float(values[0]), lambda1=float(values[1]),
                    residual_l2=float(np.linalg.norm(residual)),
                    parity_error=float(np.max(abs(u-u[::-1]))),
                    growth_time_over_b=1/math.sqrt(-values[0]) if values[0]<0 else None)

    def gaussian_rayleigh(self):
        """Continuum test u=exp(-x^2), on the full line to roundoff.

        At |x|>=12 the omitted integrals are exponentially negligible.
        A strictly negative quotient independently proves a negative mode.
        """
        L = min(self.extent, 12.)
        if L < 8:
            raise ValueError('Rayleigh tail bound requires extent>=8')
        kinetic = quad(lambda x: 4*x*x*math.exp(-2*x*x), -L, L, epsabs=1e-12)[0]
        potential = quad(lambda x: float(self.potential(x))*math.exp(-2*x*x), -L, L, epsabs=1e-12)[0]
        norm = math.sqrt(math.pi/2)
        return dict(kinetic=kinetic, potential=potential, norm=norm,
                    quotient=(kinetic+potential)/norm,
                    integration_half_domain=L)
