"""Maxwell two-port and variational loads on the repository's S3 bridge.

Geometry: -N(s)^2 dt^2 + stretch^2 ds^2 + r(s,waist)^2 dOmega^2.
The neutral spatial profile is the one in physical_throat.VacuumThroat. Its
charged extension solves the time-symmetric Einstein-Maxwell constraint in
the bridge and matches the same rim radius and slope; the global exterior is
not re-solved here. charge_fraction is q/[sin(a)*(1-cos(a))].
A positive lapse is an explicit supported Lorentzian completion, NOT vacuum
Schwarzschild. Scattering and radiation loads follow from the same Maxwell
action. No support compliance, throat inertia or mouth acceleration is assumed.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss

from .gr_handle_checks import stress_components


@dataclass(frozen=True)
class Bridge:
    angle: float = .3
    lapse_depth: float = 0.
    waist: float = 0.
    stretch: float = 1.
    charge_fraction: float = 0.

    def __post_init__(self):
        if not np.isfinite([self.angle, self.lapse_depth, self.waist, self.stretch, self.charge_fraction]).all():
            raise ValueError('Finite geometry parameters required')
        if not 0 < self.angle < np.pi/4 or not 0 <= self.lapse_depth < 1:
            raise ValueError('Require 0<angle<pi/4 and 0<=lapse_depth<1')
        if abs(self.waist) > .3 or not .5 < self.stretch < 2:
            raise ValueError('Require |waist|<=.3 and .5<stretch<2')
        if not 0 <= self.charge_fraction < 1:
            raise ValueError('Require charge below the extremal matched-bridge limit')

    @property
    def geometric_charge(self):
        # q=Q/sqrt(4*pi), with Maxwell action -F^2/4 and G=1.
        return self.charge_fraction*np.sin(self.angle)*(1-np.cos(self.angle))

    @property
    def mass(self):
        rim = np.sin(self.angle)
        return (rim**3+self.geometric_charge**2/rim)/2

    @property
    def discriminant(self):
        return np.sqrt(self.mass**2-self.geometric_charge**2)

    @property
    def b(self):
        return self.mass+self.discriminant

    @property
    def limit(self):
        return .5*np.arccosh((np.sin(self.angle)-self.mass)/self.discriminant)

    @property
    def half_length(self):
        x = self.limit
        return 2*self.mass*x+self.discriminant*np.sinh(2*x)

    def profile(self, x):
        x = np.asarray(x)
        r0 = self.mass+self.discriminant*np.cosh(2*x)
        s = 2*self.mass*x+self.discriminant*np.sinh(2*x)
        y = s/self.half_length
        h = (1-y*y)**3
        return dict(radius=r0*np.exp(self.waist*h), lapse=1-self.lapse_depth*h,
                    jacobian=2*r0, shape=h, s=s)

    def support(self, x):
        """8*pi stress of the undeformed static completion; no material law."""
        if self.waist != 0 or self.stretch != 1:
            raise ValueError('Support audit is for the undeformed profile')
        p = self.profile(x)
        r, y = p['radius'], p['s']/self.half_length
        first = self.discriminant*np.sinh(2*np.asarray(x))/r
        second = self.mass/r**2-self.geometric_charge**2/r**3
        hp = -6*y*(1-y*y)**2/self.half_length
        hpp = (-6+36*y*y-30*y**4)/self.half_length**2
        return stress_components(r, first, second, p['lapse'],
                                 -self.lapse_depth*hp, -self.lapse_depth*hpp)


def solve_port(bridge, omega, ell=1, rtol=2e-10, dense=True):
    """Exact finite-port characteristic S; no fitted area impedance.

    Canonical radial variables u and j=N*d_s u satisfy
        d_s u=j/N, d_s j=(N*l(l+1)/r^2-omega^2/N)*u.
    Characteristic amplitudes at each rim use impedance Z=omega. This is a
    boundary response to be attached to an exterior, not an asymptotic
    transmission probability for the complete closed S3 spacetime.
    """
    if omega <= 0 or not np.isfinite([omega, rtol]).all() or rtol <= 0:
        raise ValueError('Finite positive frequency and tolerance required')
    if int(ell) != ell or ell < 1:
        raise ValueError('Radiative Maxwell harmonics require integer ell>=1')
    angular = ell*(ell+1)
    def rhs(x, y):
        p = bridge.profile(x)
        r, n = p['radius'], p['lapse']
        jac = bridge.stretch*p['jacobian']
        return np.array([jac*y[1]/n, jac*(n*angular/r**2-omega**2/n)*y[0]])
    sol = solve_ivp(rhs, (bridge.limit, -bridge.limit), np.array([1.+0j, 1j*omega]),
                    method='DOP853', rtol=rtol, atol=rtol*.01, dense_output=dense)
    if not sol.success:
        raise RuntimeError(sol.message)
    u, j = sol.y[:, -1]
    incoming, reflected = .5*(u+j/(1j*omega)), .5*(u-j/(1j*omega))
    t, r = 1/incoming, reflected/incoming
    scattering = np.array([[r, t], [t, r]])
    result = dict(omega=float(omega), ell=int(ell), scattering=scattering,
                  reflection=float(abs(r)**2), transmission=float(abs(t)**2),
                  unitarity_error=float(np.max(abs(scattering.conj().T@scattering-np.eye(2)))),
                  evaluations=sol.nfev, solution=sol, incoming=incoming)
    return result


def field_basis(port, x):
    """Internal fields for left/right unit characteristic-amplitude incidence."""
    left = port['solution'].sol(x)/port['incoming']
    right = port['solution'].sol(-np.asarray(x))/port['incoming']
    right[1] *= -1  # derivative reverses under reflection about the waist
    return np.stack([left[0], right[0]], axis=-1), np.stack([left[1], right[1]], axis=-1)


def force_matrices(bridge, omega, ell=1, order=192, rtol=2e-10):
    """Time-averaged generalized loads per incident energy flux.

    Incoming vector a has |a|^2 equal to incident power. For a shape q,
    <F_q>=a^dagger W_q a. Waist r->r*exp(epsilon*h) holds proper length and
    lapse fixed; stretch scales proper length at fixed reference r,N.
    W_waist is positive semidefinite, including coherent counter-incidence.
    """
    port = solve_port(bridge, omega, ell, rtol)
    nodes, weights = leggauss(order)
    x, weights = nodes*bridge.limit, weights*bridge.limit
    p = bridge.profile(x)
    u, j = field_basis(port, x)
    r, n, jac, h = p['radius'], p['lapse'], p['jacobian'], p['shape']
    angular = ell*(ell+1)
    def gram(values, w):
        return np.einsum('ni,n,nj->ij', values.conj(), w, values)
    waist = gram(u, weights*bridge.stretch*jac*n*angular*h/(omega**2*r**2))
    # j=N/stretched_metric_factor * du/d_reference_s.
    length = (gram(u, weights*jac*(omega**2/n-n*angular/r**2))
              +gram(j, weights*jac/n))/(2*omega**2)
    stored = (gram(u, weights*bridge.stretch*jac*(omega**2/n+n*angular/r**2))
              +gram(j, weights*bridge.stretch*jac/n))/(2*omega**2)
    port.update(waist_force=waist, stretch_force=length, stored_energy=stored,
                waist_eigenvalues=np.linalg.eigvalsh(waist), stretch_eigenvalues=np.linalg.eigvalsh(length))
    return port


def scattering_force_derivative(bridge, omega, ell=1, parameter='waist', step=2e-5):
    """Independent scattering-phase identity W_q=(-i S† d_q S)/omega.

    No internal stress integral is used here. Port coordinates and lapse at
    the rims stay fixed when varying either admitted parameter.
    """
    from dataclasses import replace
    if parameter not in ('waist', 'stretch') or step <= 0:
        raise ValueError('Positive step and waist or stretch variation required')
    value = getattr(bridge, parameter)
    center = solve_port(bridge, omega, ell, rtol=2e-12, dense=False)['scattering']
    plus = solve_port(replace(bridge, **{parameter: value+step}), omega, ell, rtol=2e-12, dense=False)['scattering']
    minus = solve_port(replace(bridge, **{parameter: value-step}), omega, ell, rtol=2e-12, dense=False)['scattering']
    return -1j*center.conj().T@((plus-minus)/(2*step))/omega


def boundary_dtn(bridge, omega, ell=1, rtol=2e-11):
    """Real finite-interval Dirichlet-to-Neumann map, independently integrated.

    A pole is a finite-interval Dirichlet resonance, not nonunitary scattering.
    This routine is a validator at frequencies away from its poles.
    """
    angular = ell*(ell+1)
    def rhs(x, flat):
        p = bridge.profile(x)
        n, r = p['lapse'], p['radius']
        matrix = bridge.stretch*p['jacobian']*np.array([[0, 1/n], [n*angular/r**2-omega**2/n, 0]])
        return (matrix @ flat.reshape(2, 2)).ravel()
    sol = solve_ivp(rhs, (-bridge.limit, bridge.limit), np.eye(2).ravel(),
                    method='DOP853', rtol=rtol, atol=rtol*.01)
    if not sol.success:
        raise RuntimeError(sol.message)
    a, b, c, d = sol.y[:, -1]
    if abs(b) < 1e-13:
        raise ValueError('At or too close to a Dirichlet resonance')
    return np.array([[a, -1], [-1, d]])/b


def realify(matrix):
    """Real map on (Re a, Im a), for explicit time-reversal composition."""
    matrix = np.asarray(matrix)
    return np.block([[matrix.real, -matrix.imag], [matrix.imag, matrix.real]])


def static_flux(bridge, charge=1., order=192):
    """Topological radial electric flux, separate from radiative ell>=1 modes.

    At fixed Q, E_hat_s=Q/(4*pi*r^2). Conserved Killing energy is
    U=Q^2/(8*pi)*integral N/r^2 d(proper length). Thus it opens the waist but
    pulls toward shorter intrinsic length. No relation between intrinsic neck
    length and the mouths' separation in S3 is assumed.
    """
    if not np.isfinite(charge):
        raise ValueError('Finite charge required')
    x, w = leggauss(order)
    x, w = x*bridge.limit, w*bridge.limit
    p = bridge.profile(x)
    integrand = w*p['jacobian']*p['lapse']/p['radius']**2
    base_energy = charge**2/(8*np.pi)*np.sum(integrand)
    return dict(charge=float(charge), energy=float(bridge.stretch*base_energy),
                waist_force=float(charge**2/(4*np.pi)*bridge.stretch*np.sum(integrand*p['shape'])),
                stretch_force=float(-base_energy),
                length_force=float(-base_energy/(2*bridge.half_length)),
                scope='Fixed-flux radial field; intrinsic neck tension, not a derived S3 separation force')


def seam_map(spatial, reverse_time=False):
    """A SPECIFIED orthogonal angular map; time reversal conjugates phasors.

    The radial metric does not choose a charged/polarized Klein bundle.
    This transports declared data without treating conjugation as a complex
    gain or claiming that real orthogonality selects a global sourced history.
    """
    spatial = np.asarray(spatial, float)
    if spatial.ndim != 2 or spatial.shape[0] != spatial.shape[1]:
        raise ValueError('Square real spatial map required')
    if not np.allclose(spatial.T@spatial, np.eye(len(spatial)), atol=1e-12, rtol=0):
        raise ValueError('Orthogonal spatial map required')
    z = np.zeros_like(spatial)
    return np.block([[spatial, z], [z, (-1 if reverse_time else 1)*spatial]])
