"""Finite-time parametric locking test; no spatial absorber is represented."""
import numpy as np
from scipy.integrate import simpson, solve_ivp
from . import nonlinear_supported_tt as dynamics
from .action_selection import U, V

FREEZE = 'aa7ff3d19d61c8555502eb34cb0021aa13411698'
VOLUME = 2*np.pi**2
J_BG = 3*VOLUME/4
FLOOR = 1e-6*J_BG
AMPLITUDES = np.r_[0., .02*1.5**np.arange(8)]
SHAPE_PHASES = np.arange(4)*np.pi/4
FIELD_PHASES = np.arange(3)*np.pi/4
DURATIONS = np.array([1., 2., 4., 8.])
TIMES = np.r_[0., np.concatenate([t+np.arange(-5,6)*.05 for t in DURATIONS])]


def preparation(epsilon, theta, phi):
    return dynamics.initial_data(np.cos(theta)*U, -np.sin(theta)*V, epsilon,
                                 departure=.15, phase=phi)


def quantities(y):
    """Readouts and conformal-time work rates derived from q''=-Omega^2 q."""
    A, _, q, p, M, L = dynamics.unpack(y)
    H, Hp, r, ell, trinv, force, inv = dynamics.ingredients(y)
    omega2 = trinv+(r+ell)/6
    if omega2 <= 0 or not np.isfinite(omega2):
        raise ArithmeticError('nonpositive or nonfinite instantaneous scalar frequency')
    omega = np.sqrt(omega2)
    Lp = -Hp/H*L+force
    trinvp = -2*np.trace(inv @ L)
    rp = -8*np.trace((inv+M @ M) @ L)
    ellp = 2*np.trace(L @ Lp)
    omegap = (trinvp+(rp+ellp)/6)/(2*omega)
    Q, P, cross = q @ q, p @ p, q @ p
    action = VOLUME*np.array([(P+4*Q)/4, (P+omega2*Q)/(2*omega)])
    conformal_work = VOLUME/2*np.array([(4-omega2)*cross, omegap*(Q-P/omega2)])
    return action, conformal_work, omega2


def rhs(t, augmented):
    y = augmented[:29]
    _, work, _ = quantities(y)
    return np.r_[dynamics.conformal_rhs(y), work]/y[0]


def evolve(epsilon, theta, phi, method='DOP853'):
    if method not in ('DOP853','RK45'):
        raise ValueError('unregistered integrator')
    y0 = np.r_[preparation(epsilon,theta,phi), 0., 0.]
    rtol,atol,step = (1e-12,1e-14,.02) if method=='DOP853' else (1e-10,1e-12,.01)
    sol = solve_ivp(rhs,(0.,TIMES[-1]),y0,t_eval=TIMES,method=method,
                    rtol=rtol,atol=atol,max_step=step)
    if not sol.success or not np.isfinite(sol.y).all():
        raise ArithmeticError('integration failed: '+sol.message)
    return sol.y.T


def summarize_history(states):
    states = np.asarray(states, dtype=float)
    if states.shape != (len(TIMES),31) or not np.isfinite(states).all():
        raise ValueError('incomplete or nonfinite frozen history')
    action = np.array([quantities(y[:29])[0] for y in states])
    change = action-action[0]
    means, coarse = [], []
    for j in range(4):
        sl = slice(1+11*j,1+11*(j+1))
        times, values = TIMES[sl], change[sl]
        means.append(simpson(values,x=times,axis=0)/.5)
        coarse.append(simpson(values[::2],x=times[::2],axis=0)/.5)
    constraints = [dynamics.constraints(y[:29]) for y in states]
    det, sym, minimum = 0., 0., np.full(4,np.inf)
    for y in states:
        A,_,q,_,M,L = dynamics.unpack(y)
        H = A*A-q @ q/6
        omega2 = quantities(y[:29])[2]
        minimum = np.minimum(minimum,[A,H,np.linalg.eigvalsh((M+M.T)/2).min(),omega2])
        det = max(det,abs(np.linalg.det(M)-1))
        sym = max(sym,np.max(abs(M-M.T)),np.max(abs(M @ L-(M @ L).T)))
    return dict(initial_readouts=action[0].tolist(), changes=change.tolist(),
                window_means=np.array(means).tolist(), coarse_window_means=np.array(coarse).tolist(),
                quadrature_error_over_Jbg=float(np.max(abs(np.array(means)-coarse))/J_BG),
                work_error_over_Jbg=float(np.max(abs(change-states[:,29:31]))/J_BG),
                constraint_max=float(max(max(c['normalized']) for c in constraints)),
                det_error=float(det),symmetry_error=float(sym),sampled_min_A_H_M_omega2=minimum.tolist(),
                energy_terms_initial=constraints[0]['energy_terms'].tolist(),
                energy_terms_final=constraints[-1]['energy_terms'].tolist())


def classify(means):
    """Frozen test on (8 amplitudes,12 phase pairs,3 late times,2 readouts).

    Invalid/missing numbers cannot become a negative physical verdict. No
    absolute values, clipping, per-run normalization or fitted action scale.
    """
    means = np.asarray(means,dtype=float)
    if means.shape != (8,12,3,2) or not np.isfinite(means).all():
        raise ValueError('all registered plateau readouts must be finite and present')
    windows = []
    for start in range(5):
        values = means[start:start+4]
        resolved = values > FLOOR
        positive = bool(resolved.all())
        median = float(np.median(values))
        spread = float(np.ptp(values)/median) if positive else None
        valid = resolved[1:] & resolved[:-1]
        slopes = np.full_like(values[1:],np.nan)
        np.log(np.divide(values[1:],values[:-1],out=np.ones_like(values[1:]),where=valid),out=slopes,where=valid)
        slopes /= np.log(1.5)
        max_slope = float(np.max(abs(slopes[valid]))) if valid.any() else None
        gates = dict(all_positive_resolved=positive, common_value=positive and spread<=.1,
                     flat_amplitude=bool(valid.all()) and max_slope<=.1)
        windows.append(dict(amplitudes=AMPLITUDES[1+start:5+start].tolist(),
                            unresolved_or_nonpositive=int(np.size(resolved)-resolved.sum()),
                            signed_min=float(values.min()),signed_max=float(values.max()),
                            common_median=median, relative_range=spread,
                            max_resolved_abs_log_slope=max_slope,gates=gates,passes=bool(all(gates.values()))))
    return windows


def selection_verdict(numerical_gates, windows):
    if not numerical_gates or not all(numerical_gates.values()):
        return 'INCONCLUSIVE_NUMERICAL_FAILURE'
    if len(windows)!=5:
        return 'INCONCLUSIVE_NUMERICAL_FAILURE'
    return ('CANDIDATE_FINITE_TIME_PLATEAU_NOT_QUANTIZATION' if any(w['passes'] for w in windows)
            else 'NO_ROBUST_PLATEAU_IN_REGISTERED_FAMILY')
