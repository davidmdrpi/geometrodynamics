"""Budgeted external section resets for the diagonal family; not autonomous GR."""
import json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import solve_ivp
from . import r3_family as rf, nonlinear_supported_tt as d

DETUNINGS = (-.001, -.0005, .0005, .001)
METHODS = ('DOP853', 'Radau')
ARMS = ('unforced', 'controlled')


def state(z, eta=0.):
    out = rf.to_state(np.r_[np.asarray(z), np.zeros(6)])
    out[2] = eta
    return out


def coords(y):
    return np.r_[rf.to_section(y)[:6], y[7]]


def horizon(twist, detuning):
    if not np.isfinite([twist, detuning]).all() or twist*detuning == 0:
        raise ValueError('invalid drift rate')
    n = int(np.ceil(1/(2*abs(twist*detuning))))
    if n > 256:
        raise ValueError('twist horizon exceeds fixed cap')
    return n


def decision(distance, d0, cost, cumulative, controlled, completed, phase_drift, target):
    vals = [distance, d0, cost, cumulative, phase_drift]
    if not np.isfinite(vals).all() or d0 <= 1e-6 or min(distance, cost, cumulative) < 0:
        raise ValueError('invalid diagnostic')
    if distance > 2*d0:
        return 'TUBE_ESCAPE'
    if controlled and (cost > .02*d0 or cumulative+cost > .25*d0):
        return 'CONTROL_BUDGET_EXCEEDED'
    if completed == target:
        return 'HORIZON_COMPLETED' if abs(phase_drift) >= .5 else 'DRIFT_HORIZON_NOT_RESOLVED'
    return 'CONTINUE'


class Family:
    def __init__(self, F, S):
        self.Z = np.array([F['start']['v'][:6]]+[p['v'][:6] for p in F['points']])
        closed = np.vstack([self.Z, self.Z[0]])
        arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))]
        self.theta = 2*np.pi*arc[:-1]/arc[-1]
        self.curve = CubicSpline(2*np.pi*arc/arc[-1], closed, bc_type='periodic')
        phases, projectors = [], []
        for smp in S['samples']:
            M = np.array(smp['M2'])[:6, :6]
            ev, V = np.linalg.eig(M); mask = (abs(ev) < .1) | (abs(ev) > 10)
            if sum(mask) != 2:
                raise ValueError('wrong hyperbolic dimension')
            H = V@np.diag(mask.astype(float))@np.linalg.inv(V)
            if np.max(abs(H.imag)) > 1e-8:
                raise ValueError('complex hyperbolic projector')
            phases.append(self.theta[smp['index']]); projectors.append(H.real)
        self.projector = CubicSpline(np.r_[phases, 2*np.pi], np.array(projectors+[projectors[0]]), bc_type='periodic')
        self.grid = 2*np.pi*np.arange(256)/256
        self.reference_grid = np.array([self.reference(t) for t in self.grid])
        c = F['circles']; poly = np.polyfit([x['omega'] for x in c], [x['action'] for x in c], 2)
        self.twist = float(1/np.polyval(np.polyder(poly), np.pi))
        self.action = float(rf.loop_action(self.Z))

    def reference(self, theta):
        return coords(state(self.curve(theta % (2*np.pi))))

    def nearest(self, w):
        w = np.asarray(w)
        err = np.linalg.norm(self.reference_grid-w, axis=1)**2
        loc = np.flatnonzero((err <= np.roll(err, 1)) & (err <= np.roll(err, -1)))
        if not len(loc):
            raise ArithmeticError('no phase-grid minimum')
        best = None; step = 2*np.pi/256
        for i in loc:
            t = self.grid[i]
            f = lambda u: float(np.sum((self.reference(u)-w)**2))
            sol = minimize_scalar(f, bounds=(t-step, t+step), method='bounded', options={'xatol':1e-12})
            if not sol.success or not np.isfinite(sol.fun):
                raise ArithmeticError('phase minimization failed')
            if best is None or sol.fun < best[0]:
                best = (sol.fun, sol.x % (2*np.pi))
        return float(best[1]), float(np.sqrt(best[0]))

    def tune(self, z, theta):
        ref = self.curve(theta % (2*np.pi))
        L = np.linalg.svd(self.projector(theta % (2*np.pi)))[2][:2]
        condition = np.linalg.cond(L[:, :2])
        if condition > 100:
            raise ValueError('ill-conditioned fixed controller')
        tuned = np.asarray(z).copy()
        tuned[:2] = ref[:2]-np.linalg.solve(L[:, :2], L[:, 2:]@(tuned[2:]-ref[2:]))
        return tuned, float(condition)

    def trial_loop(self, scale):
        out = []
        for z, t in zip(self.Z, self.theta):
            trial = z.copy(); trial[2:] *= scale
            out.append(self.tune(trial, t)[0])
        return np.array(out)

    def initial(self, detuning):
        f = lambda s: rf.loop_action(self.trial_loop(s))-self.action-detuning
        scale = brentq(f, .9, 1.1, xtol=1e-13)
        loop = self.trial_loop(scale); achieved = float(rf.loop_action(loop)-self.action)
        if abs(achieved-detuning) > 1e-10:
            raise ValueError('action proxy does not reproduce')
        y = state(loop[0]); phase, d0 = self.nearest(coords(y))
        if d0 <= 1e-6:
            raise ValueError('degenerate initial displacement')
        untuned = self.Z[0].copy(); untuned[2:] *= scale
        initial_kick = coords(y)-coords(state(untuned))
        return dict(detuning=detuning, scale=float(scale), loop=loop.tolist(), achieved_action_offset=achieved,
                    initial=y.tolist(), d0=d0, initial_phase=phase, target=horizon(self.twist,detuning),
                    tuning_vector=initial_kick.tolist(), tuning_norm=float(np.linalg.norm(initial_kick)),
                    log10_autonomous_unstable_tolerance=float(np.log10(2*d0)-horizon(self.twist,detuning)*np.log10(7242.)))


def one_return(y, method):
    fun = lambda t, u: d.conformal_rhs(u)
    opts = dict(method=method, rtol=2e-12, atol=2e-14, max_step=.025, dense_output=True)
    a = solve_ivp(fun, (0., 1e-6), y, **opts)
    if not a.success:
        raise ArithmeticError(a.message)
    def event(t, u): return u[3]
    event.direction = -1; event.terminal = True
    b = solve_ivp(fun, (1e-6, np.pi+.8), a.y[:, -1], events=event, **opts)
    if not b.success or len(b.t_events[0]) != 1:
        raise ArithmeticError('clock return not reached')
    end = float(b.t_events[0][0]); times = np.linspace(0., end, 33)
    history = np.array([a.sol(t) if t <= 1e-6 else b.sol(t) for t in times])
    for row in history:
        d.ingredients(row)
        if np.max(abs(d.constraints(row)['residual'])) > 1e-8:
            raise ArithmeticError('constraint accuracy failed')
    return history[-1], dict(times=times.tolist(), states=history.tolist())


def simulate(family, initial, method, arm):
    y = np.array(initial['initial']); d0 = initial['d0']; target = initial['target']
    phase = initial['initial_phase']; drift = 0.; cumulative = 0.
    out = dict(method=method, arm=arm, steps=[], terminal=None)
    for n in range(1, target+1):
        rec = dict(step=n, start=y.tolist(), returns=[])
        try:
            for _ in range(2):
                y, history = one_return(y, method); rec['returns'].append(history)
            rec['pre'] = y.tolist(); t, distance = family.nearest(coords(y))
            inc = float(np.angle(np.exp(1j*(t-phase))))
            if abs(inc) >= np.pi/2:
                raise ArithmeticError('phase unwrap not resolved')
            drift += inc; phase = t
            reason = decision(distance,d0,0.,cumulative,False,n,drift,target)
            proposal = y.copy(); cost = 0.; post_distance = distance
            if reason != 'TUBE_ESCAPE' and arm == 'controlled':
                z, condition = family.tune(coords(y)[:6], t)
                proposal = state(z, eta=y[2]); cost = float(np.linalg.norm(coords(proposal)-coords(y)))
                rec['controller_condition'] = condition
                _, post_distance = family.nearest(coords(proposal))
                reason = decision(max(distance,post_distance),d0,cost,cumulative,True,n,drift,target)
            applied = arm == 'controlled' and reason not in ('TUBE_ESCAPE','CONTROL_BUDGET_EXCEEDED')
            rec.update(phase=t, phase_drift=drift, pre_distance=distance, proposed=proposal.tolist(),
                       proposed_cost=cost, proposed_distance=post_distance, applied=applied)
            if applied: y=proposal; cumulative+=cost
            rec.update(after=y.tolist(), cumulative_cost=cumulative)
            if reason != 'CONTINUE': out['terminal'] = reason
        except (ArithmeticError, ValueError, np.linalg.LinAlgError) as exc:
            out['terminal'] = 'NUMERICALLY_UNRESOLVED'; rec['error'] = str(exc)
        out['steps'].append(rec)
        if out['terminal']: break
    return out
