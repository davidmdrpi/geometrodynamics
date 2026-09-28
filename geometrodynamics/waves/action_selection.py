"""Canonical circulation controls; no absorbed-action receiver is modeled.

I is an integral over a closed *preparation loop*, not an orbit integral along
one spacetime solution. Sector circulations depend on the stated canonical
chart. Use the full constrained Einstein/quartet dynamics without projection.
"""
import numpy as np
from scipy.integrate import solve_ivp
from . import nonlinear_supported_tt as dynamics

FREEZE = '9f18137441ecb4d9158792e8d4eef11079f0ddea'
TIMES = np.array([0., .25, .5, 1., 2.])
AMPLITUDES = (0., .01, .02, .04, .08)
VOLUME = 2*np.pi**2
H0 = 1.15**2-1/8
U = np.diag([1., -1., 0.])/np.sqrt(2)
W = np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., 0.]])/np.sqrt(2)
V = U+W/2
SECTORS = ('scale', 'scalar', 'shape', 'shape_plus_scalar', 'total')


def initial_loop(epsilon, nodes=64, phase_offset=0.):
    if not np.isfinite(epsilon) or epsilon < 0:
        raise ValueError('finite nonnegative preparation amplitude required')
    if not isinstance(nodes, (int, np.integer)) or nodes < 8 or nodes % 2:
        raise ValueError('even loop grid of at least 8 nodes required')
    if not np.isfinite(phase_offset):
        raise ValueError('finite loop phase offset required')
    theta = 2*np.pi*np.arange(nodes)/nodes+phase_offset
    return np.array([dynamics.initial_data(np.cos(t)*U, -np.sin(t)*V, epsilon,
                                           departure=.15, phase=0.) for t in theta])


def evolve_loop(epsilon, nodes=64, method='DOP853'):
    initial = initial_loop(epsilon, nodes)
    if method not in ('DOP853', 'RK45'):
        raise ValueError('registered integrators are DOP853 and RK45')
    trajectories = []
    # At epsilon=0 every loop member is the same physical solution.
    for y0 in initial[:1] if epsilon == 0 else initial:
        if method == 'DOP853':
            states = dynamics.evolve(y0, TIMES)
        else:
            sol = solve_ivp(dynamics.proper_rhs, (0., TIMES[-1]), y0, t_eval=TIMES,
                            method='RK45', rtol=1e-10, atol=1e-12, max_step=.01)
            if not sol.success or not np.isfinite(sol.y).all():
                raise ArithmeticError('independent integration failed: '+sol.message)
            states = sol.y.T
        trajectories.append(states)
    return np.repeat(trajectories, nodes, axis=0) if epsilon == 0 else np.array(trajectories)


def loop_derivative(values):
    """d/dtheta on the first, uniform periodic axis, with endpoint excluded."""
    values = np.asarray(values)
    count = len(values)
    if count < 8 or count % 2 or not np.isfinite(values).all():
        raise ValueError('finite periodic samples on an even grid required')
    modes = np.fft.fftfreq(count, d=1/count)
    modes[count//2] = 0.  # derivative of the real Nyquist cosine at its nodes
    shape = (count,)+(1,)*(values.ndim-1)
    return np.fft.ifft(1j*modes.reshape(shape)*np.fft.fft(values, axis=0), axis=0).real


def circulation(states):
    """I=Vol/(2pi) integral Theta, returning signed sector arrays over time."""
    states = np.asarray(states)
    if states.ndim != 3 or states.shape[-1] != 29:
        raise ValueError('states must have shape (loop nodes, times, 29)')
    dy = loop_derivative(states)
    A, Ap = states[:, :, 0], states[:, :, 1]
    q, qp = states[:, :, 3:7], states[:, :, 7:11]
    M, L = states[:, :, 11:20].reshape(*states.shape[:2], 3, 3), states[:, :, 20:29].reshape(*states.shape[:2], 3, 3)
    H = A*A-np.sum(q*q, axis=-1)/6
    if np.min(A) <= 0 or np.min(H) <= 0 or np.min(np.linalg.eigvalsh((M+np.swapaxes(M,-1,-2))/2)) <= 0:
        raise ValueError('outside positive canonical chart')
    Pi = H[:, :, None, None]/2*(L @ np.linalg.inv(M))
    dM = dy[:, :, 11:20].reshape(*states.shape[:2], 3, 3)
    integrands = {'scale': -6*Ap*dy[:, :, 0],
                  'scalar': np.sum(qp*dy[:, :, 3:7], axis=-1),
                  'shape': np.einsum('ntij,ntji->nt', Pi, dM)}
    out = {key: VOLUME*np.mean(value, axis=0) for key, value in integrands.items()}
    out['shape_plus_scalar'] = out['shape']+out['scalar']
    out['total'] = out['scale']+out['shape_plus_scalar']
    return out


def initial_prediction(epsilon):
    return np.pi**2*H0*epsilon**2


def diagnostics(states):
    values = np.asarray(states).reshape(-1, 29)
    constraint, det, symmetric = 0., 0., 0.
    minimum = np.full(3, np.inf)
    for y in values:
        A, _, q, _, M, L = dynamics.unpack(y)
        H = A*A-q @ q/6
        minimum = np.minimum(minimum, [A, H, np.linalg.eigvalsh((M+M.T)/2).min()])
        constraint = max(constraint, float(np.max(dynamics.constraints(y)['normalized'])))
        det = max(det, abs(np.linalg.det(M)-1))
        symmetric = max(symmetric, float(np.max(abs(M-M.T))), float(np.max(abs(M @ L-(M @ L).T))))
    return dict(sampled_constraint_max=constraint, sampled_det_error=float(det),
                sampled_symmetry_error=symmetric, sampled_minimum_A_H_M=minimum.tolist())


def receiver_readiness():
    missing = ['constraint-monitored localized source/receiver evolution',
               'separated seed perturbations and field/metric interaction controls',
               'operational absorbed-action observable and full transfer ledger',
               'physical preparation-pulse duration sweep',
               'receiver-definition robustness test']
    return dict(verdict='NOT_READY_FOR_RECEIVER_ACTION_SELECTION', missing=missing,
                reason='Homogeneous sector circulation does not measure absorption by a receiver.')
