"""Exact diagonal Bianchi IX section map; prospective freeze 0c51f9b."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from . import nonlinear_supported_tt as full

FREEZE = '0c51f9b5b528432c012fc15a4765af8829bb3867'
E = np.array([[1., 1., -2.], [np.sqrt(3.), -np.sqrt(3.), 0.]])/np.sqrt(6.)
STAR = np.array([1., 0., 0., 0., 0., 0.])
RTOL, ATOL, STEP = 2e-12, 2e-14, .025


def phase_rhs(phi, state):
    """Array of states (...,7), including elapsed conformal time."""
    A, u, x, v, y, w, eta = np.moveaxis(state, -1, 0)
    beta = x[..., None]*E[0]+y[..., None]*E[1]
    inv, squared = np.exp(-2*beta), np.exp(4*beta)
    T = inv.sum(-1)
    curvature = 2*(2*T-squared.sum(-1))
    ell = v*v+w*w
    si, co = np.sin(phi), np.cos(phi)
    R2 = (3*u*u+A*A*(curvature-ell)/2-1.5*A**4)/(2*si*si+co*co*(T/2+(curvature-ell)/12))
    Q = R2*co*co
    H, Hp = A*A-Q/6, 2*A*u+2/3*R2*si*co
    clock = 2*si*si+(T+(curvature+ell)/6)*co*co/2
    if not np.isfinite(state).all() or min(np.real(a).min() for a in (A, R2, H, clock)) <= 0:
        raise ValueError('left positive diagonal phase chart')
    force = (((-4+Q/H)[..., None]*inv-4*squared) @ E.T)
    return np.stack([u, -A*(curvature+ell)/6+A**3, v, -Hp/H*v+force[..., 0],
                     w, -Hp/H*w+force[..., 1], np.ones_like(A)], -1)/clock[..., None]


def phase_map(points, jacobian=False, tol=(RTOL, ATOL), step=STEP):
    """Batched exact flow; complex-step derivatives of the same section map."""
    points = np.atleast_2d(points)
    n = len(points)
    h = 1e-25
    if jacobian:
        inputs = np.repeat(points[:, None, :], 7, axis=1).astype(complex)
        for k in range(6):
            inputs[:, k+1, k] += 1j*h
        inputs = inputs.reshape(-1, 6)
    else:
        inputs = points
    initial = np.zeros((len(inputs), 7), dtype=inputs.dtype)
    initial[:, 0] = inputs[:, 0]
    initial[:, 1] = -inputs[:, 1]/6
    initial[:, 2] = inputs[:, 2]
    initial[:, 3] = inputs[:, 3]/inputs[:, 0]**2
    initial[:, 4] = inputs[:, 4]
    initial[:, 5] = inputs[:, 5]/inputs[:, 0]**2
    shape = initial.shape
    sol = solve_ivp(lambda t, y: phase_rhs(t, y.reshape(shape)).ravel(), (np.pi/2, 5*np.pi/2),
                    initial.ravel(), method='DOP853', rtol=tol[0], atol=tol[1], max_step=step)
    if not sol.success:
        raise ArithmeticError(sol.message)
    s = sol.y[:, -1].reshape(shape)
    result = np.column_stack([s[:, 0], -6*s[:, 1], s[:, 2], s[:, 0]**2*s[:, 3],
                              s[:, 4], s[:, 0]**2*s[:, 5]])
    if jacobian:
        result = result.reshape(n, 7, 6)
        return result[:, 0].real, np.moveaxis(result[:, 1:].imag/h, 1, 2)
    return result.real, s[:, 6].real


def full_initial(z):
    A, pA, x, px, y, py = z
    M = np.diag(np.exp(2*(x*E[0]+y*E[1])))
    L = np.diag((px*E[0]+py*E[1])/A**2)
    state = full.pack(A, -pA/6, np.zeros(4), np.zeros(4), M, L)
    energy = full.constraints(state)['residual'][0]
    if energy >= 0:
        raise ValueError('no negative real clock root')
    state[7] = -np.sqrt(-2*energy)
    return state


def canonical(state):
    A, Ap, q, qp, M, L = full.unpack(np.asarray(state))
    shape = E @ np.log(np.diag(M))/2
    momenta = (A*A-q@q/6)*(E @ np.diag(L))
    return np.array([A, -6*Ap, shape[0], momenta[0], shape[1], momenta[1]])


def full_return(z, method='DOP853', initial=None):
    y0 = full_initial(z) if initial is None else np.asarray(initial).copy()
    y0[2] = 0.
    def positive(t, y): return y[3]
    positive.direction, positive.terminal = 1, True
    def negative(t, y): return y[3]
    negative.direction, negative.terminal = -1, True
    args = dict(method=method, rtol=RTOL, atol=ATOL, max_step=.01, dense_output=True)
    first = solve_ivp(lambda t, y: full.conformal_rhs(y), (0., 8.), y0, events=positive, **args)
    if not first.success or len(first.t_events[0]) != 1:
        raise ArithmeticError('positive crossing failed: '+first.message)
    tmid = first.t[-1]
    last = solve_ivp(lambda t, y: full.conformal_rhs(y), (tmid, 8.), first.y[:, -1], events=negative, **args)
    if not last.success or len(last.t_events[0]) != 1:
        raise ArithmeticError('negative crossing failed: '+last.message)
    times = np.linspace(0., last.t[-1], 257)
    samples = np.array([first.sol(t) if t <= tmid else last.sol(t) for t in times])
    return dict(initial=y0, returned=last.y[:, -1], samples=samples, times=times,
                return_time=float(last.t[-1]), method=method)


def spectral(K, theta):
    k = np.fft.fftfreq(len(K), 1/len(K))
    return (np.exp(1j*np.asarray(theta)[:, None]*k) @ (np.fft.fft(K, axis=0)/len(K))).real


def shift_matrix(n, omega):
    k = np.fft.fftfreq(n, 1/n)
    f = np.fft.fft(np.eye(n), axis=0)
    return (np.fft.ifft(np.exp(1j*k*omega)[:, None]*f, axis=0).real,
            np.fft.ifft((1j*k*np.exp(1j*k*omega))[:, None]*f, axis=0).real)


def action(K):
    n = len(K)
    k = np.fft.fftfreq(n, 1/n)
    c = np.fft.fft(K, axis=0)/n
    return float(sum(np.sum(c[:, p].conj()*1j*k*c[:, q]) for p, q in [(1, 0), (3, 2), (5, 4)]).real)


def circular_area(K):
    n = len(K)
    k = np.fft.fftfreq(n, 1/n)
    c = np.fft.fft(K, axis=0)/n
    return float(np.sum(c[:, 2].conj()*1j*k*c[:, 4]).real)


def linear_seed(amplitude, grid=63):
    _, jac = phase_map([STAR], True)
    vals, vectors = np.linalg.eig(jac[0, 2:4, 2:4])
    i = int(np.argmax(vals.imag))
    q = vectors[:, i]/vectors[0, i]
    wave = q[None, :]*np.exp(2j*np.pi*np.arange(grid)/grid)[:, None]
    K = np.tile(STAR, (grid, 1))
    K[:, 2:4] += amplitude*wave.real
    K[:, 4:6] += amplitude*(1j*wave).real
    return K, float(np.angle(vals[i]))


def circle(amplitude, K, omega, mapper=phase_map):
    K = np.array(K, copy=True)
    n, dim = K.shape
    e1 = np.exp(-2j*np.pi*np.arange(n)/n)/n
    trace = []
    def residual(curve, om, image):
        T, _ = shift_matrix(n, om)
        c1 = e1 @ curve[:, 2]
        return np.r_[(image-T@curve).ravel(), c1.real-amplitude/2, c1.imag]
    try:
        for it in range(20):
            image, blocks = mapper(K, True)
            res = residual(K, omega, image)
            err = float(np.max(abs(res)))
            trace.append(dict(iteration=it, residual=err, omega=float(omega)))
            if err <= 5e-11:
                break
            T, dT = shift_matrix(n, omega)
            J = np.zeros((dim*n+2, dim*n+1))
            for i in range(n):
                J[dim*i:dim*(i+1), dim*i:dim*(i+1)] = blocks[i]
            J[:dim*n, :dim*n] -= np.kron(T, np.eye(dim))
            J[:dim*n, -1] = -(dT@K).ravel()
            J[-2, 2:dim*n:dim], J[-1, 2:dim*n:dim] = e1.real, e1.imag
            delta = np.linalg.lstsq(J, -res, rcond=None)[0]
            accepted = False
            for h in range(11):
                f = 2.**-h
                trial, om = K+f*delta[:-1].reshape(K.shape), omega+f*delta[-1]
                try:
                    image_t, _ = mapper(trial)
                    rtrial = np.max(abs(residual(trial, om, image_t)))
                except (ValueError, ArithmeticError, FloatingPointError) as exc:
                    trace[-1]['trial_error'] = str(exc)
                    continue
                if rtrial < err:
                    K, omega, accepted = trial, float(om), True
                    trace[-1]['step_factor'] = f
                    break
            if not accepted:
                raise ArithmeticError('Newton line search failed')
        image, times = mapper(K)
        return dict(a=float(amplitude), K=K, omega=float(omega), image=image, return_times=times, trace=trace)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
        return dict(a=float(amplitude), K=K, omega=float(omega), trace=trace, error=str(exc))


def validate_circle(rec):
    if 'error' in rec:
        return rec
    K, omega = np.asarray(rec['K']), rec['omega']
    n = len(K)
    mid = 2*np.pi*(np.arange(n)+.5)/n
    points = spectral(K, mid)
    image, times = phase_map(points)
    rec.update(offgrid_points=points, offgrid_image=image, offgrid_times=times)
    angles = np.arange(4)*np.pi/2
    points = spectral(K, angles)
    image, times = phase_map(points)
    histories = [full_return(p) for p in points]
    rec.update(full_points=points, full_phase_image=image, full_phase_times=times, full_histories=histories)
    return rec


def shoot(seed_nodes):
    trace = []
    def fun(flat):
        nodes = flat.reshape(2, 6)
        try:
            image, _ = phase_map(nodes)
            res = image-nodes[::-1]
        except (ValueError, ArithmeticError):
            res = np.full((2, 6), 1e6)
        trace.append(dict(residual=float(np.max(abs(res))), nodes=nodes.copy()))
        return res.ravel()
    def jac(flat):
        _, blocks = phase_map(flat.reshape(2, 6), True)
        return np.block([[blocks[0], -np.eye(6)], [-np.eye(6), blocks[1]]])
    try:
        sol = least_squares(fun, np.asarray(seed_nodes).ravel(), jac=jac,
                            xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=80)
        nodes = sol.x.reshape(2, 6)
        image, times = phase_map(nodes)
        return dict(seed_nodes=seed_nodes, nodes=nodes, image=image, times=times, trace=trace,
                    solver_success=bool(sol.success), solver_message=sol.message)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
        return dict(seed_nodes=seed_nodes, trace=trace, error=str(exc))
