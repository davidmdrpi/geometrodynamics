"""Einstein evolution controls for the matched spherical throat, not a BAM router.

Electrovacuum develops the spatial Bridge into Reissner--Nordstrom. The central
normal geodesic contracts; it does not describe two mouth centers moving on S3.
The charged calculation stops at the inner horizon. No stable bounce, topology
change, material support law, or finite-pulse recoil is assumed.
"""
import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss
from .maxwell_throat_response import solve_port


def initial_acceleration(bridge, x, support_radial_pressure=0., static_lapse=False):
    """Proper-time areal acceleration at K_ij=0 from the angular ADM equation.

    The support density and momentum are zero on this slice. Its radial
    pressure is an additional initial datum, not a constitutive law. With
    static_lapse=False the normal observers are geodesic (N=1). With True use
    the declared positive lapse's initial spatial gradient, not staticity.
    """
    if bridge.waist != 0 or bridge.stretch != 1:
        raise ValueError('Use the constraint-satisfying undeformed initial slice')
    if not np.isfinite(support_radial_pressure).all():
        raise ValueError('Finite pressure required')
    x = np.asarray(x, float)
    if not np.isfinite(x).all() or np.any(abs(x) > bridge.limit+1e-12):
        raise ValueError('Point must lie inside the bridge')
    p = bridge.profile(x)
    r, s = p['radius'], p['s']
    second = bridge.mass/r**2-bridge.geometric_charge**2/r**3
    acceleration = -second-4*np.pi*r*support_radial_pressure
    if static_lapse:
        first = bridge.discriminant*np.sinh(2*x)/r
        y = s/bridge.half_length
        lapse_first = 6*bridge.lapse_depth*y*(1-y*y)**2/bridge.half_length
        acceleration += first*lapse_first/p['lapse']
    return acceleration


def central_history(bridge, eta):
    """Exact electrovacuum central normal geodesic up to the first inner end.

    r=M+d*cos(eta), tau=M*eta+d*sin(eta). q=0 ends at a curvature
    singularity eta=pi (excluded); q>0 ends at a Cauchy horizon eta=pi.
    The RN time-translation Killing energy of this geodesic is zero.
    """
    if bridge.waist != 0 or bridge.stretch != 1:
        raise ValueError('Use the undeformed constraint-satisfying slice')
    eta = np.asarray(eta, float)
    if not np.isfinite(eta).all() or np.any(eta < 0) or np.any(eta > np.pi):
        raise ValueError('Require 0<=eta<=pi; no continuation beyond the inner horizon')
    m, q, d = bridge.mass, bridge.geometric_charge, bridge.discriminant
    if q == 0 and np.any(eta == np.pi):
        raise ValueError('Neutral singular endpoint has no regular state')
    # Stable at eta near pi even for small nonzero charge.
    inner = q*q/(m+d)
    r = inner+2*d*np.cos(eta/2)**2
    tau = m*eta+d*np.sin(eta)
    u = -d*np.sin(eta)/r
    acceleration = -m/r**2+q*q/r**3
    curvature = 48*m*m/r**6-96*m*q*q/r**7+56*q**4/r**8
    return dict(tau=tau, radius=r, velocity=u, acceleration=acceleration,
                null_expansion_plus=2*u/r, null_expansion_minus=2*u/r,
                kretschmann=curvature)


def integrate_central(bridge, times, rtol=2e-10):
    """Independent dimensionless second-order geodesic ODE, not eta inversion."""
    if bridge.waist != 0 or bridge.stretch != 1:
        raise ValueError('Use the undeformed constraint-satisfying slice')
    times = np.asarray(times, float)
    if (times.ndim != 1 or len(times) < 2 or not np.isfinite(times).all()
            or times[0] != 0 or np.any(np.diff(times) <= 0)
            or times[-1] >= np.pi*bridge.mass or not np.isfinite(rtol) or rtol <= 0):
        raise ValueError('Ordered times starting at zero and ending before pi*M required')
    b, m, q = bridge.b, bridge.mass, bridge.geometric_charge
    def rhs(t, y):
        return [y[1], -(m/b)/y[0]**2+(q/b)**2/y[0]**3]
    sol = solve_ivp(rhs, (0., times[-1]/b), [1., 0.], t_eval=times/b,
                    method='DOP853', rtol=rtol, atol=rtol*.01)
    if not sol.success:
        raise RuntimeError(sol.message)
    return dict(radius=b*sol.y[0], velocity=sol.y[1], evaluations=sol.nfev)


def optical_length(bridge, order=384):
    nodes, weights = leggauss(order)
    p = bridge.profile(nodes*bridge.limit)
    return float(np.sum(weights*bridge.limit*bridge.stretch*p['jacobian']/p['lapse']))


def group_delays(bridge, omega, ell=1, step=2e-4, rtol=2e-11):
    """Finite-port phase derivatives for exp(-i omega t), not arrival gates.

    Conditional reflection/transmission delays require nonzero amplitudes.
    A broadband output is obtained with the full complex S(omega), not one
    energy-weighted delay. The reference planes and Z=omega are unchanged.
    """
    if not np.isfinite(step) or step <= 0 or omega <= step:
        raise ValueError('Require omega>step>0')
    center = solve_port(bridge, omega, ell, rtol, dense=False)
    s = center['scattering']
    plus = solve_port(bridge, omega+step, ell, rtol, dense=False)['scattering']
    minus = solve_port(bridge, omega-step, ell, rtol, dense=False)['scattering']
    derivative = (plus-minus)/(2*step)
    smith = -1j*s.conj().T@derivative
    def conditional(i):
        z = s[i, 0]
        return None if abs(z) < 1e-12 else float(np.imag(derivative[i, 0]/z))
    return dict(omega=float(omega), reflection=center['reflection'], transmission=center['transmission'],
                reflection_delay=conditional(0), transmission_delay=conditional(1),
                smith=smith, mean_delay=float(smith[0, 0].real),
                hermiticity_error=float(np.max(abs(smith-smith.conj().T))))
