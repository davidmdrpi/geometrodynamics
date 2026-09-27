"""Explicit l=2 polar TT packets on the breathing S3 background.

Linear evolution only. Squared Weyl curvature is not an energy density.
The three radial components encode the exact S2-integrated tensor norm.
"""
import numpy as np
from scipy.special import eval_gegenbauer, roots_legendre, gammaln
from scipy.integrate import solve_ivp
from . import esu_floquet as fl

DEGREES = np.arange(2, 81)
WINDOWS = ((12, 3), (24, 6), (40, 10))
PHASES = (0., np.pi/4)
AMPLITUDES = (1e-6, 1e-5, 1e-4)
TIMES = np.unique(np.r_[np.linspace(0, np.pi+.15, 801), np.pi,
                        np.pi+np.linspace(-.15, .15, 301)])
CAP = .3
REGIONS = {'full': (0., np.pi), 'north': (0., CAP),
           'south': (np.pi-CAP, np.pi),
           'belt': (np.pi/2-.15, np.pi/2+.15)}


def radial(chi, degrees=DEGREES):
    """sqrt(3/2)A, sqrt(12)B/s, sqrt(12)D/s^2, before normalization.

    Gegenbauer derivatives are evaluated directly, not by the eigen-equation.
    Returned shape (degree, component, point), valid away from polar axes.
    """
    chi = np.atleast_1d(chi)
    s, c = np.sin(chi), np.cos(chi)
    out = []
    for n in degrees:
        m = int(n)-2
        A = eval_gegenbauer(m, 3, c)
        Ax = 6*eval_gegenbauer(m-1, 4, c) if m >= 1 else np.zeros_like(c)
        Axx = 48*eval_gegenbauer(m-2, 5, c) if m >= 2 else np.zeros_like(c)
        Ap, App = -s*Ax, -c*Ax+s*s*Axx
        B = s*s*Ap/6+s*c*A/2
        Bp = s*s*App/6+5*s*c*Ap/6+(c*c-s*s)*A/2
        D_over_s2 = (Bp+2*c/s*B-A/2)/2
        out.append([np.sqrt(1.5)*A, np.sqrt(12)*B/s,
                    np.sqrt(12)*D_over_s2])
    return np.asarray(out)


def quadrature(region='full', nodes=512):
    a, b = REGIONS[region]
    x, w = roots_legendre(nodes)
    chi = (a+b)/2+(b-a)*x/2
    return chi, w*(b-a)/2*np.sin(chi)**2


def normalization():
    chi, w = quadrature()
    H = radial(chi)
    return np.sqrt(np.einsum('ncp,ncp,p->n', H, H, w))


def grams(norms, nodes=512):
    result = {}
    for region in REGIONS:
        chi, w = quadrature(region, nodes)
        H = radial(chi)/norms[:, None, None]
        result[region] = np.einsum('ncp,mcp,p->nm', H, H, w)
    return result


def power_factors(norms, nodes=512):
    """Uncontracted quadrature factors: square the reconstructed field first.

    A precontracted Gram quadratic form loses tiny positive packet tails to
    cancellation between modes. These factors retain positivity without clipping.
    """
    factors = {}
    for region in REGIONS:
        chi, w = quadrature(region, nodes)
        H = radial(chi)/norms[:, None, None]
        factors[region] = (H*np.sqrt(w)[None,None,:]).reshape(len(DEGREES),-1)
    return factors


def coefficients(center, width, norms, even=False):
    # C_m^3(1) = Gamma(m+6)/(Gamma(6) Gamma(m+1)).
    m = DEGREES-2
    pole = np.exp(gammaln(m+6)-gammaln(6)-gammaln(m+1))/norms
    c = np.exp(-.5*((DEGREES-center)/width)**2)*pole
    if even:
        c *= (DEGREES % 2 == 0)
    return c/np.linalg.norm(c)


def modes(method='DOP853'):
    """Integrate two unit scaled-state columns for each n, independently.

    State is (h,v=h'/(n+1)), making absolute tolerance uniform over degree.
    Shape (time, degree, component, initial-column).
    """
    values = []
    rtol, atol = (1e-11, 1e-13) if method == 'DOP853' else (1e-10, 1e-12)
    for n in DEGREES:
        w = n+1
        def rhs(t, y):
            R, _, f, fp, _ = fl.background(t)
            M = y.reshape(2, 2)
            return np.array([w*M[1], -(n*(n+2)+2*R*R/f)*M[0]/w-fp*M[1]/f]).ravel()
        sol = solve_ivp(rhs, (0, TIMES[-1]), np.eye(2).ravel(),
                        method=method, rtol=rtol, atol=atol, t_eval=TIMES)
        if not sol.success:
            raise ArithmeticError(sol.message)
        values.append(sol.y.T.reshape(-1, 2, 2))
    return np.moveaxis(np.asarray(values), 0, 1)


def analytic_modes(kind):
    w = DEGREES+1
    freq = w if kind == 'free' else np.sqrt(DEGREES*(DEGREES+2))
    phase = TIMES[:, None]*freq
    M = np.empty((len(TIMES), len(w), 2, 2))
    M[:, :, 0, 0] = M[:, :, 1, 1] = np.cos(phase)
    M[:, :, 0, 1] = np.sin(phase)*w/freq
    M[:, :, 1, 0] = -np.sin(phase)*freq/w
    return M


def packet_state(matrices, c, phase):
    return np.einsum('tnij,n,j->tni', matrices, c, [np.cos(phase), np.sin(phase)])


def observables(state, factors, kind='supported'):
    h, hp = state[:, :, 0], state[:, :, 1]*(DEGREES+1)
    k = DEGREES*(DEGREES+2)
    R, _, f, fp, _ = fl.background(TIMES[:, None])
    if kind == 'supported':
        hpp = -fp*hp/f-(k+2*R*R/f)*h
    else:
        hpp = -((DEGREES+1)**2 if kind == 'free' else k)*h
    E = (k*h-hpp)/(4*f)
    stress = -R*R*h/(f*f) if kind == 'supported' else np.zeros_like(h)
    fields = dict(metric=h, weyl=E, stress=stress)
    powers = {field: {r: np.sum((a@F)**2,axis=1)
                      for r, F in factors.items()} for field, a in fields.items()}
    return fields, powers


def volume(region):
    a, b = REGIONS[region]
    return 4*np.pi*((b-a)/2-(np.sin(2*b)-np.sin(2*a))/4)


def scalar_budget(method='DOP853'):
    S = fl.energy_scaling('S', 2)
    inv = np.linalg.inv(S)
    def rhs(t, y):
        M = inv@y.reshape(4, 4)
        return (S@fl.rhs_scalar(t, M, 2)).ravel()
    rtol, atol = (1e-11, 1e-13) if method == 'DOP853' else (1e-10, 1e-12)
    times = TIMES[TIMES <= np.pi]
    sol = solve_ivp(rhs, (0, np.pi), np.eye(4).ravel(), method=method,
                    rtol=rtol, atol=atol, t_eval=times)
    if not sol.success:
        raise ArithmeticError(sol.message)
    matrices = sol.y.T.reshape(-1, 4, 4)
    transfer = []
    for t, M in zip(times, matrices):
        Psi, Phi, _, _ = fl.scalar_metric(t, inv@M, 2)
        transfer.append([Phi, Psi])
    return matrices, np.asarray(transfer)
