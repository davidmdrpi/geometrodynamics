"""Conditional MTY wave graph with two dynamical response coordinates.

This is a reduced scalar-field action, NOT Einstein evolution or mouth
centre-of-mass recoil. Equal final clock rates; clock preparation is prescribed.
Three ports per node: antipodal bulk channel, handle, and source/receiver lead.
All energies include the leads. No confirmation waveform or force is inserted.
"""
from dataclasses import dataclass, asdict
import warnings
import numpy as np
from scipy.signal import lfilter
from scipy.optimize import anderson, NoConvergence
from .network import NetworkMouth

MASS = np.array([.1, .15])
STIFFNESS = np.array([1., 1.4])
QUARTIC = np.array([10., 15.])
BULK_TIME = 1.
HANDLE_TIME = .125


@dataclass(frozen=True)
class Config:
    aging_duration: float
    amplitude: float
    dt: float = 1/128
    start: float = -16.
    stop: float = 64.
    connected: bool = True
    seed: int = 0
    scheme: str = 'discrete_gradient'
    basis_B: int = 1
    nonlinear: bool = True


def clock_history(duration):
    """Ideal flat-clock out-and-back history at |v|=sqrt(3)/2, then rest.

    Two inertial segments, each duration D/2; ideal instantaneous accelerations.
    Per unit preparation rest mass: kinetic energy at speed is 1, positive
    actuator work is 2, recovered work is 2, net zero. For D=0 no excursion.
    This is not a solved moving-throat metric or a free preparation operation.
    """
    if not np.isfinite(duration) or duration < 0:
        raise ValueError('invalid aging duration')
    proper = duration/2
    return dict(duration=duration, speed=np.sqrt(3)/2, proper_duration=proper,
                offset=duration-proper, final_clock_rate=1.,
                positive_work_per_rest_mass=2. if duration else 0.,
                recovered_work_per_rest_mass=2. if duration else 0.)


def delays(config):
    delta = clock_history(config.aging_duration)['offset']
    A = NetworkMouth('A', 0., 'mty-packet', clock_offset=0.)
    B = NetworkMouth('B', np.pi, 'mty-packet', clock_offset=delta)
    # B is the mouth which aged less. Entries are arrival node A, B.
    return np.array([[BULK_TIME, A.global_time(B.local_time(0.)+HANDLE_TIME)],
                     [BULK_TIME, B.global_time(A.local_time(0.)+HANDLE_TIME)]])


def grid(config):
    if not np.isfinite([config.dt,config.start,config.stop,config.amplitude]).all():
        raise ValueError('nonfinite configuration')
    if config.dt <= 0 or config.stop <= config.start or config.amplitude <= 0:
        raise ValueError('invalid interval/amplitude')
    if config.scheme not in ('discrete_gradient', 'midpoint') or config.basis_B not in (-1, 1):
        raise ValueError('unknown scheme/basis')
    n = round((config.stop-config.start)/config.dt)
    if abs(n*config.dt-(config.stop-config.start)) > 1e-12:
        raise ValueError('interval not on grid')
    return config.start+(np.arange(n)+.5)*config.dt


def packet(t, amplitude):
    """Compact C3 envelope, support [-1/2,1/2], fixed two-cycle carrier."""
    t = np.asarray(t)
    return amplitude*np.where(abs(t)<.5, np.cos(np.pi*t)**4*np.cos(4*np.pi*t), 0.)


def shift(values, delay, dt):
    """Nonperiodic shift; zero outside the finite history. Never FFT-wrap."""
    k = round(delay/dt)
    if abs(k*dt-delay) > 1e-11:
        raise ValueError('delay must be an integer number of time cells')
    out = np.zeros_like(values)
    if abs(k) >= len(values):
        return out
    if k > 0: out[k:] = values[:-k]
    elif k < 0: out[:k] = values[-k:]
    else: out[:] = values
    return out


def linear_response(force, dt):
    """Trapezoidal linear oscillator, zero initial q,p; exact digital filter.

    M vdot + 3v + Kq = force. Return q,v on N+1 cell boundaries.
    Nonlinear forces are solved simultaneously with network fields below.
    """
    n = len(force); q = np.zeros((n+1, 2)); v = q.copy()
    for j in range(2):
        L = np.array([[0., 1.], [-STIFFNESS[j]/MASS[j], -3/MASS[j]]])
        inv = np.linalg.inv(np.eye(2)-dt*L/2)
        A = inv@(np.eye(2)+dt*L/2)
        B = inv@np.array([0., dt/MASS[j]])
        tr = np.trace(A); den = [1., -tr, np.linalg.det(A)]
        numer = A@B-tr*B
        q[1:,j] = lfilter([B[0], numer[0]], den, force[:,j])
        v[1:,j] = lfilter([B[1], numer[1]], den, force[:,j])
    return q, v


def nonlinear_force(q, config):
    if not config.nonlinear: return np.zeros_like(q[:-1])
    a, b = q[:-1], q[1:]
    if config.scheme == 'midpoint': return QUARTIC*((a+b)/2)**3
    return QUARTIC*(a**3+a*a*b+a*b*b+b**3)/4


def response(x, config):
    t = grid(config); n = len(t)
    x = np.asarray(x).reshape(n, 2, 3)
    if not np.isfinite(x).all(): raise ArithmeticError('nonfinite iteration')
    incoming = np.zeros((n, 2, 3)); incoming[:,:,:2] = x[:,:,:2]
    incoming[:,0,2] = packet(t, config.amplitude)
    q, v = linear_response(2*np.sum(incoming,axis=2)-x[:,:,2], config.dt)
    if not np.isfinite(q).all() or np.max(abs(q)) > 100*config.amplitude:
        raise ArithmeticError('response exceeds fixed iteration domain')
    middle_v = (v[:-1]+v[1:])/2
    outgoing = middle_v[:,:,None]-incoming
    target = np.zeros_like(x); ds = delays(config)
    for j in range(2):
        target[:,j,0] = config.basis_B*shift(outgoing[:,1-j,0], ds[j,0], config.dt)
        if config.connected:
            # Twisted scalar handle is an explicit line-bundle assumption.
            target[:,j,1] = -config.basis_B*shift(outgoing[:,1-j,1], ds[j,1], config.dt)
    target[:,:,2] = nonlinear_force(q, config)
    return target, dict(q=q, v=v, incoming=incoming, outgoing=outgoing)


def solve(config):
    t = grid(config); x = np.zeros((len(t),2,3))
    if config.seed:
        z = .05*config.amplitude*np.exp(-((t-1.)/2)**2)*np.sin(3*t)
        x[:,:,:2] = z[:,None,None]*np.array([[1.,-1.],[-.7,.7]])
    calls = 0
    def residual(flat):
        nonlocal calls
        calls += 1
        target, _ = response(flat, config)
        return target.ravel()-flat
    status = 'CONVERGED'; message = ''
    with warnings.catch_warnings(record=True) as caught:
        try:
            solution = anderson(residual, x.ravel(), alpha=.5, M=10, w0=.01,
                                maxiter=1500, f_tol=2e-9*config.amplitude,
                                line_search='armijo')
        except NoConvergence as exc:
            solution = exc.args[0]; status = 'NUMERICALLY_UNRESOLVED'
        except (ArithmeticError, ValueError, np.linalg.LinAlgError) as exc:
            return dict(config=asdict(config), status='NUMERICALLY_UNRESOLVED',
                        message=str(exc), calls=calls)
        message = '; '.join(sorted(set(str(w.message) for w in caught)))
    target, fields = response(solution, config)
    error = float(np.max(abs(target.ravel()-solution))/config.amplitude)
    if error > 1e-8: status = 'NUMERICALLY_UNRESOLVED'
    elif status != 'CONVERGED': status = 'RESIDUAL_ACCEPTED_AT_ITERATION_LIMIT'
    return dict(config=asdict(config), status=status, message=message, calls=calls,
                residual=error, iterate=np.asarray(solution).reshape(len(t),2,3).tolist(),
                **{key:value.tolist() for key,value in fields.items()})


def energy(q, v, nonlinear=True):
    return MASS*v*v/2+STIFFNESS*q*q/2+(QUARTIC*q**4/4 if nonlinear else 0.)


def diagnose(record):
    """Recompute equations, transport, energy and canonical-response momentum.

    Canonical p=Mv is NOT longitudinal radiation momentum or gravitating-mouth
    COM momentum. The on-site potential transfers generalized impulse to an
    explicit fixed anchor; its work is included as potential energy.
    """
    config = Config(**record['config']); t = grid(config); n = len(t)
    if 'iterate' not in record: return dict(valid=False, reason='no finite iterate')
    x = np.array(record['iterate']); target, expected = response(x, config)
    for key, shape in [('q',(n+1,2)),('v',(n+1,2)),('incoming',(n,2,3)),('outgoing',(n,2,3))]:
        a = np.asarray(record[key])
        if a.shape != shape or not np.isfinite(a).all() or np.max(abs(a-expected[key])) > 1e-10*config.amplitude:
            raise ValueError('altered field/evolution: '+key)
    q,v,inc,out = (expected[k] for k in ('q','v','incoming','outgoing'))
    dt = config.dt; vm = (v[:-1]+v[1:])/2; qm = (q[:-1]+q[1:])/2
    residual = float(np.max(abs(target-x))/config.amplitude)
    input_energy = float(dt*np.sum(inc[:,:,2]**2))
    E = energy(q,v,config.nonlinear)
    flux = dt*np.sum(inc**2-out**2,axis=2)
    cumulative = E-E[0]-np.vstack([np.zeros(2),np.cumsum(flux,axis=0)])
    balance = float(np.max(abs(cumulative))/input_energy)
    restoring = STIFFNESS*qm+nonlinear_force(q,config)
    forces = np.sum(inc-out,axis=2)
    momentum_error = MASS*(v-v[0])-np.vstack([np.zeros(2),np.cumsum(dt*(forces-restoring),axis=0)])
    pscale = config.amplitude
    mom_error = float(np.max(abs(momentum_error))/pscale)
    kinematic = float(np.max(abs(np.diff(q,axis=0)-dt*vm))/config.amplitude)
    ports = out[:,:,2]**2
    loss = float(dt*np.sum(ports))
    if not config.connected: loss += float(dt*np.sum(out[:,:,1]**2))
    transport_leak = float(dt*np.sum(out[:,:,:2]**2-inc[:,:,:2]**2))
    if not config.connected: transport_leak -= float(dt*np.sum(out[:,:,1]**2))
    # Derivative matching alone can hide an additive field mismatch on an
    # advanced edge. Check field values integrated from zero in the past too.
    times=config.start+np.arange(n+1)*dt;ds=delays(config)
    incoming_field=np.concatenate([np.zeros((1,2,3)),np.cumsum(inc,axis=0)*dt])
    outgoing_field=np.concatenate([np.zeros((1,2,3)),np.cumsum(out,axis=0)*dt])
    field_error=0.
    for j in range(2):
        for port in range(2 if config.connected else 1):
            sign=config.basis_B*(1 if port==0 else -1)
            expected_field=sign*np.interp(times-ds[j,port],times,outgoing_field[:,1-j,port],left=0.)
            field_error=max(field_error,float(np.max(abs(incoming_field[:,j,port]-expected_field)))/config.amplitude)
    total_balance = abs(input_energy-loss-float(np.sum(E[-1]-E[0])))/input_energy
    edge = (t < config.start+2)|(t > config.stop-2)
    edge_energy = float(dt*np.sum((inc[edge]**2+out[edge]**2)) / input_energy)
    endpoint_energy = float(np.sum(E[-1]+E[0])/input_energy)
    pre = t < -.5
    early = float(dt*np.sum(out[pre,0,2]**2)/input_energy)
    # Physical positive handle energy in its own clock tau, not on a global
    # exterior-time slice: B's clock tau=t_B-Delta.
    delta = clock_history(config.aging_duration)['offset']
    bp = shift(out[:,1,1], -delta, dt)
    intensity = out[:,0,1]**2+bp**2
    width = round(HANDLE_TIME/dt)
    cumulative_h = np.r_[0.,np.cumsum(intensity)*dt]
    storage = cumulative_h[1:]-cumulative_h[np.maximum(0,np.arange(1,n+1)-width)]
    valid = bool(residual<=1e-8 and balance<= (2e-3 if config.scheme=='midpoint' else 1e-6)
                 and mom_error<=1e-7 and kinematic<=1e-8 and total_balance<=1e-4
                 and edge_energy<=1e-4 and endpoint_energy<=1e-4 and field_error<=1e-5)
    return dict(valid=valid, residual=residual, local_energy_error=balance,
                canonical_momentum_error=mom_error, kinematic_error=kinematic,
                input_energy=input_energy, output_energy=loss,
                full_history_energy_error=float(total_balance),
                transport_window_defect=transport_leak/input_energy,
                integrated_field_matching_error=field_error,
                edge_energy_fraction=edge_energy, endpoint_energy_fraction=endpoint_energy,
                before_source_energy_fraction=early,
                peak_response_momenta=np.max(abs(v*MASS),axis=0).tolist(),
                peak_handle_energy=float(np.max(storage)) if config.connected else 0.)
