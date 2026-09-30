"""Refocusing-resonance rotation number of the exact homogeneous n=2 tensor mode.

Freeze: docs/r3_refocusing_resonance_prereg.md (2e984ac), correction 6b55c5b.
Locally rotationally symmetric (Taub) subsystem of nonlinear_supported_tt:
beta = x b0, q = (q0,0,0,0). The clock section is q0 = 0, q0' < 0. Orbits are
kept on the centre-stable manifold of the Einstein-static instability by
repeated bisection of a kick to A (straddle method).
"""
import numpy as np
from scipy.integrate import solve_ivp
from . import nonlinear_supported_tt as d

B0 = np.diag([1., 1., -2.])/np.sqrt(6)
TOLERANCES = {'primary': (1e-12, 1e-14), 'secondary': (1e-10, 1e-12)}
HORIZON, ACCEPT, K = 9, 4, 48
EXIT, SAMPLES, START_IGNORE = .5, 256, 1.
BRACKETS = ((1e-2, .32), (1e-6, 1e-3))   # (initial, widening cap): first window, later windows
WIDEN, BISECT_WIDTH, BISECT_MAX = 4., 4e-16, 70


def tensor(y):
    """LRS amplitude x = tr(beta b0) and x' = tr(L b0)."""
    _, _, _, _, M, L = d.unpack(y)
    return float(np.log(np.diag(M)) @ np.diag(B0)/2), float(np.trace(L @ B0))


def resolve_clock_velocity(y):
    """Solve the Hamiltonian constraint for q0' on the negative root."""
    y = y.copy()
    E = d.constraints(y)['residual'][0]
    base = E-y[7]**2/2
    if base >= 0:
        raise ValueError('no real clock velocity')
    y[7] = -np.sqrt(-2*base)
    return y


def section_state(eps, pol, s=0.):
    x = pol*eps
    M = np.diag(np.exp(2*x*np.diag(B0)))
    y = d.pack(1.+s, 0., np.zeros(4), np.zeros(4), M, np.zeros((3, 3)))
    return resolve_clock_velocity(y)


def _flow(y0, t_end, tol, dense=False, exit_events=True):
    rtol, atol = tol

    def clock(t, y):
        return y[3]
    clock.direction = -1

    def up(t, y):
        return y[0]-1-EXIT
    up.terminal = True

    def down(t, y):
        return y[0]-1+EXIT
    down.terminal = True
    events = [clock, up, down] if exit_events else [clock]
    return solve_ivp(lambda t, y: d.conformal_rhs(y), (0., t_end), y0, method='DOP853',
                     rtol=rtol, atol=atol, events=events, dense_output=dense)


def classify(y_start, s, tol):
    """+1 expanding, -1 collapsing (or chart failure), over HORIZON clock periods."""
    y = y_start.copy()
    y[0] += s
    try:
        y = resolve_clock_velocity(y)
        sol = _flow(y, HORIZON*np.pi+.5, tol)
    except ValueError:
        return -1
    if sol.status == 1:
        if len(sol.t_events[1]):
            return 1
        return -1
    if not sol.success:
        return -1
    return 1 if sol.y[0, -1] > 1 else -1


def bisect(y_start, tol, first):
    lo_cap = BRACKETS[0] if first else BRACKETS[1]
    width, widenings = lo_cap[0], 0
    while True:
        lo, hi = -width, width
        clo, chi = classify(y_start, lo, tol), classify(y_start, hi, tol)
        if clo != chi:
            break
        if width*WIDEN > lo_cap[1]*(1+1e-12):
            raise ArithmeticError('bracket does not straddle the centre-stable manifold')
        width *= WIDEN
        widenings += 1
    for _ in range(BISECT_MAX):
        if hi-lo < BISECT_WIDTH:
            break
        mid = (lo+hi)/2
        if classify(y_start, mid, tol) == clo:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2, widenings


def track(eps, pol, tol, periods=K):
    """Collect `periods` accepted clock periods; return per-period data."""
    y = section_state(eps, pol)
    increments, residuals, A_dev, radii, windows = [], [], [], [], []
    first = True
    while len(increments) < periods:
        s, widenings = bisect(y, tol, first)
        windows.append(dict(kick=float(s), widenings=int(widenings)))
        first = False
        y0 = y.copy()
        y0[0] += s
        y0 = resolve_clock_velocity(y0)
        sol = _flow(y0, (ACCEPT+.5)*np.pi+.5, tol, dense=True, exit_events=False)
        if not sol.success:
            raise ArithmeticError('accepted segment failed: '+sol.message)
        times = [0.]+[t for t in sol.t_events[0] if t > START_IGNORE][:ACCEPT]
        if len(times) < ACCEPT+1:
            raise ArithmeticError('fewer clock sections than accepted periods')
        for a, b in zip(times, times[1:]):
            ts = np.linspace(a, b, SAMPLES+1)
            Y = sol.sol(ts)
            xs = np.array([tensor(Y[:, j]) for j in range(len(ts))])
            theta = np.unwrap(np.arctan2(-xs[:, 1]/3, xs[:, 0]))
            increments.append(float(theta[-1]-theta[0]))
            radii.append((float(np.hypot(xs[:, 0], xs[:, 1]/3).min()), float(np.hypot(xs[:, 0], xs[:, 1]/3).max())))
            A_dev.append(float(np.abs(Y[0]-1).max()))
            residuals.append(float(abs(d.constraints(Y[:, -1])['residual'][0])))
        y = sol.sol(times[ACCEPT])
    n = periods
    return dict(increments=increments[:n], residuals=residuals[:n], A_dev=A_dev[:n],
                radii=radii[:n], windows=windows)


def birkhoff(increments):
    inc = np.asarray(increments, float)
    t = (np.arange(len(inc))+.5)/len(inc)
    w = np.exp(-1/(t*(1-t)))
    return float(w @ inc/(2*np.pi*w.sum()))
