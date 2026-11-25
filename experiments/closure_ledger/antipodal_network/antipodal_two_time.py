"""Antipodal, free-midpoint variational problem with one-way scattering probes."""
import json

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize

ANGLES = np.arange(9)*np.pi/16
TIMES = (2., 4., 8.)


def unit(x):
    return x/np.linalg.norm(x, axis=-1, keepdims=True)


def hopf(x):
    """Multiply both complex coordinates by i in real-pair notation."""
    y = np.empty_like(x)
    y[..., ::2], y[..., 1::2] = -x[..., 1::2], x[..., ::2]
    return y


def rank(x):
    return np.sum(x[..., :2]**2, axis=-1)


def arc_action(nodes, lam):
    """Exact action of the constant-speed round arc between adjacent nodes."""
    n = len(nodes)-1
    x, y = nodes[:-1], nodes[1:]
    cosines = np.clip(np.sum(x*y, axis=-1), -1., 1.)
    delta = np.arccos(cosines)
    # delta / sin(delta), with its regular zero-angle limit.
    ratio = 1/np.sinc(delta/np.pi)
    fiber = np.sum(hopf(x)*y, axis=-1)
    return n*np.sum(delta**2 + (lam-1)*(ratio*fiber)**2)/2


def solve_paths(n=8, seed=367):
    """Vary all interior nodes on the sphere. No prescribed tangent or midpoint."""
    rng = np.random.default_rng(seed)
    x = unit(np.array([.5, .7, -.4, .3]))
    temporal = np.linspace(0., np.pi, n+1)[:, None]
    rows = []
    for trial in range(3):
        v = rng.normal(size=4)
        v = unit(v-v.dot(x)*x)
        guess = np.cos(temporal)*x + np.sin(temporal)*v
        interior = unit(guess[1:-1]+.07*rng.normal(size=(n-1,4)))
        def path(y):
            return np.vstack([x, unit(y.reshape((-1,4))), -x])
        for lam in (.5, 1., 2.):
            result = minimize(lambda y: arc_action(path(y), lam), interior.ravel(),
                              method="L-BFGS-B", options=dict(ftol=1e-14,gtol=1e-7,maxiter=600))
            if not result.success:
                raise RuntimeError(result.message)
            nodes = path(result.x)
            y = nodes[n//2]
            rows.append(dict(trial=trial, lam=lam, action=float(result.fun),
                             lower_bound=float(min(1.,lam)*np.pi**2/2),
                             midpoint_rank=float(rank(y)), initial_rank=float(rank(x)),
                             vertical_component=float(np.dot(hopf(x),y))))
    return rows


def flight_time(energy, barrier):
    """Proper time from -1 to +1 on a transmitting orbit of a unit mass."""
    if not energy > barrier >= 0:
        return np.inf
    # Work with the small positive difference explicitly near the separatrix.
    gap = energy-barrier
    def inverse_speed(q):
        return 1/np.sqrt(2*(gap + barrier*q*q*(2-q*q)))
    return 2*quad(inverse_speed, 0., 1., epsabs=1e-10)[0]


def accepted_rank(energy, elapsed):
    """Measure on uniform barrier rank [0,1] crossing by this elapsed time."""
    if energy <= 0 or elapsed <= 0 or elapsed <= flight_time(energy,0.):
        return 0.
    # A simple root at finite time, bracketed by the free-flight and separatrix limits.
    return brentq(lambda b: flight_time(energy,b)-elapsed, 0., energy, xtol=1e-15)


def fractions(right, left):
    """Lengths of the two transmitting intervals and their outside complement.

    Used for complementary ranks only. At finite times the intervals shrink
    from opposite ends and remain disjoint.
    """
    return dict(right_only=float(right), left_only=float(left), both=0.,
                neither=float(max(0.,1-right-left)))


def geometric_controls(count=131072, seed=236):
    rng = np.random.default_rng(seed)
    x = unit(rng.normal(size=(count,4)))
    z = rng.normal(size=x.shape)
    y = unit(z - (x*z).sum(axis=1, keepdims=True)*x)
    u, v = rank(x), rank(y)
    round_error = abs(u+v-1)
    # Haar complex-orthogonal unit tangent, with independent midpoint phase.
    h = np.stack((-x[:,2], x[:,3], x[:,0], -x[:,1]), axis=1)
    phase = rng.uniform(0,2*np.pi,size=(count,1))
    horizontal = h*np.cos(phase) + hopf(h)*np.sin(phase)
    v_horizontal = rank(horizontal)
    return dict(count=count, seed=seed,
                horizontal_max_defect=float(np.max(abs(u+v_horizontal-1))),
                round_mean_defect=float(np.mean(round_error)),
                round_balanced_doubles=float(np.mean((u<.5)&(v<.5))),
                round_balanced_misses=float(np.mean((u>=.5)&(v>=.5))),
                vertical_mean_defect=float(np.mean(abs(2*u-1))),
                uniform_basin_max_error=float(max(abs(np.mean(u<p)-p) for p in np.cos(ANGLES)**2)),
                biased_basin_max_error_from_p_squared=float(max(
                    abs(np.sum(u*(u<p))/np.sum(u)-p*p) for p in np.cos(ANGLES)**2)))


def midpoint_loading():
    """Post-primary point-load control. A trial path bounds every global minimizer."""
    rows = []
    for k in (-.3,-.1,.1,.3):
        # x=(1,1)/sqrt(2), real horizontal semicircle. Vary its speed in
        # the two halves, keeping the endpoints and total swept angle fixed.
        delta = k/4
        trial_rank = np.sin(np.pi/4+delta)**2
        gain = k*(trial_rank-.5)-2*delta**2
        rows.append(dict(load=k, trial_shift=delta, action_improvement=float(gain),
                         guaranteed_absolute_defect=float(gain/abs(k))))
    return rows


def probe():
    rows = []
    for i,theta in enumerate(ANGLES):
        p = float(np.cos(theta)**2) if i<8 else 0.
        finite = []
        for t in TIMES:
            right = accepted_rank(p,t)
            finite.append(dict(time=t, equal_elapsed=fractions(right,accepted_rank(1-p,t)),
                               common_slice=fractions(right,accepted_rank(1-p,t-.5))))
        rows.append(dict(index=i,power=p,eventual=fractions(p,1-p),biased_right=p*p,finite=finite))
    return dict(paths=solve_paths(),ensemble=geometric_controls(),rows=rows,
                midpoint_load_followup=midpoint_loading())


if __name__ == "__main__":
    print(json.dumps(probe(),indent=2,allow_nan=False))
