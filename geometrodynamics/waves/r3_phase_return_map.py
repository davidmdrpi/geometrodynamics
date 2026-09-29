"""Exact LRS clock reduction and cubic constrained Poincare map.

Prospective freeze 259ab9c; no state resets or long-time phase estimator.
All canonical quantities here are per unit S3 coordinate volume.
"""
import numpy as np
from scipy.integrate import solve_ivp
from . import nonlinear_supported_tt as full
from .taylor_jets import Jet, ring

FREEZE = '259ab9c08b774af5e76819789a40a379172a318d'
KAPPA = 1/np.sqrt(6)
B0 = np.diag([1., 1., -2.])*KAPPA
J2 = np.array([[0., 1.], [-1., 0.]])
J4 = np.kron(np.eye(2), J2)
METHODS = {'DOP853': (2e-12, 2e-14, .025), 'RK45': (2e-12, 2e-14, .0125)}


def exp(x):
    return x.exp() if isinstance(x, Jet) else np.exp(x)


def phase_rhs(phi, state):
    A, u, x, v, eta = state
    a, b, c = exp(-2*KAPPA*x), exp(4*KAPPA*x), exp(-8*KAPPA*x)
    trinv, curvature = 2*a+b, 8*a-2*c
    si, co = np.sin(phi), np.cos(phi)
    radius2 = (3*u*u+A*A*(curvature-v*v)/2-1.5*A**4)/(2*si*si+co*co*(trinv/2+(curvature-v*v)/12))
    Q = radius2*co*co
    H = A*A-Q/6
    Hp = 2*A*u+(2/3)*radius2*si*co
    omega2 = trinv+(curvature+v*v)/6
    clock = 2*si*si+omega2*co*co/2
    if not isinstance(A, Jet) and min(A, radius2, H, clock) <= 0:
        raise ValueError('left positive clock/constraint chart')
    force = 2*KAPPA*((-4+Q/H)*(a-b)-4*(b-c))
    return [u/clock, (-A*(curvature+v*v)/6+A**3)/clock,
            v/clock, (-Hp/H*v+force)/clock, clock**-1]


def integrate_jets(method='DOP853'):
    algebra = ring(4)
    z = [algebra.variable(j) for j in range(4)]
    A = 1+z[0]
    initial = [A, -z[1]/6, z[2], z[3]/(A*A), algebra.constant()]
    def rhs(t, coefficients):
        y = [Jet(algebra, c) for c in coefficients.reshape(5, algebra.size)]
        return np.array([p.c for p in phase_rhs(t, y)]).ravel()
    rtol, atol, step = METHODS[method]
    sol = solve_ivp(rhs, (np.pi/2, 5*np.pi/2), np.array([p.c for p in initial]).ravel(),
                    method=method, rtol=rtol, atol=atol, max_step=step)
    if not sol.success or not np.isfinite(sol.y).all():
        raise ArithmeticError('variational integration failed')
    A, u, x, v, eta = [Jet(algebra, c) for c in sol.y[:, -1].reshape(5, algebra.size)]
    return dict(method=method, coefficients=np.array([p.c for p in [A-1, -6*u, x, A*A*v]]),
                return_time=eta.c, nfev=sol.nfev)


def map_polynomials(coefficients):
    return [Jet(ring(4), np.asarray(c)) for c in coefficients]


def matrix(polynomials):
    return np.array([[p.derivative(j).c[0] for j in range(p.r.n)] for p in polynomials])


def transform(mat, vector):
    return [sum(a*p for a, p in zip(row, vector)) for row in mat]


def homological(linear, rotation, forcing, degree):
    """Solve h(Rw)-B h(w)=forcing with ordinary monomial coefficients."""
    algebra = ring(2)
    ids = np.flatnonzero(algebra.degrees == degree)
    w = [algebra.variable(i) for i in range(2)]
    Rw = transform(rotation, w)
    columns = []
    for component in range(2):
        for index in ids:
            h = [algebra.constant(), algebra.constant()]
            h[component].c[index] = 1.
            lhs = [a.compose(Rw)-b for a, b in zip(h, transform(linear, h))]
            columns.append(np.concatenate([p.c[ids] for p in lhs]))
    operator = np.array(columns).T
    right = np.concatenate([p.c[ids] for p in forcing])
    answer = np.linalg.solve(operator, right)
    h = [algebra.constant(), algebra.constant()]
    for j in range(2):
        h[j].c[ids] = answer[j*len(ids):(j+1)*len(ids)]
    return h, float(np.linalg.cond(operator)), float(np.max(abs(operator@answer-right)))


def reduce_map(coefficients):
    F = map_polynomials(coefficients)
    L = matrix(F)
    B, T = L[:2, :2], L[2:, 2:]
    cosine = np.trace(T)/2
    if not abs(cosine) < 1 or T[0, 1] == 0:
        raise ValueError('not the registered elliptic block')
    sine = np.copysign(np.sqrt(1-cosine*cosine), T[0, 1])
    G = -J2@(T-cosine*np.eye(2))/sine
    eig, vectors = np.linalg.eigh((G+G.T)/2)
    if min(eig) <= 0:
        raise ValueError('nonpositive elliptic action metric')
    S = (vectors*eig**-.5)@vectors.T
    S /= np.sqrt(np.linalg.det(S))
    Sinv = np.linalg.inv(S)
    R = Sinv@T@S
    algebra = ring(2)
    w = [algebra.variable(i) for i in range(2)]
    Sw = transform(S, w)
    h = [algebra.constant(), algebra.constant()]
    conditions, residuals = [], []
    for degree in (2, 3):
        embedded = [p.compose(h+Sw) for p in F]
        center = transform(Sinv, embedded[2:])
        # Current invariance defect supplies the right side for the next h_k.
        forcing = [(a-b.compose(center)).homogeneous(degree) for a, b in zip(embedded[:2], h)]
        hk, cond, residual = homological(B, R, forcing, degree)
        h = [a+b for a, b in zip(h, hk)]
        conditions.append(cond)
        residuals.append(residual)
    embedded = [p.compose(h+Sw) for p in F]
    center = transform(Sinv, embedded[2:])
    defect = [a-b.compose(center) for a, b in zip(embedded[:2], h)]
    graph_residual = max(float(np.max(abs(p.c))) for p in defect)
    # Pullback symplectic area: (1+det Dh_2) dQ wedge dP, through degree 2.
    h2 = [p.homogeneous(2) for p in h]
    density = h2[0].derivative(0)*h2[1].derivative(1)-h2[0].derivative(1)*h2[1].derivative(0)
    dp = algebra.constant()
    for value, (i, j) in zip(density.c, algebra.powers):
        if i+j == 2:
            dp.c[algebra.index[(i, j+1)]] = -value/(j+1)
    D = [w[0], w[1]+dp]
    canonical = [p.compose(D) for p in center]
    canonical[1] = canonical[1]-dp.compose(canonical)
    quadratic = [p.homogeneous(2) for p in canonical]
    hnf, cond, residual = homological(R, R, quadratic, 2)
    conditions.append(cond)
    residuals.append(residual)
    # Conjugate by id+h_nf; g2=0, so g3=f3+Df2 h_nf.
    cubic = [p.homogeneous(3)+sum(q.derivative(j)*hnf[j] for j in range(2))
             for p, q in zip(canonical, quadratic)]
    z, zb = w
    Q = (z+zb)/np.sqrt(2)
    P = 1j*(z-zb)/np.sqrt(2)
    complex_cubic = ((cubic[0]-1j*cubic[1])/np.sqrt(2)).compose([Q, P])
    lam = complex(cosine, sine)
    a21 = complex(complex_cubic.c[algebra.index[(2, 1)]])
    coefficient = a21/lam
    rho0 = 1+np.arctan2(sine, cosine)/(2*np.pi)
    return dict(nu=float(coefficient.imag/(2*np.pi)), rho0=float(rho0),
                radial_resonant=float(coefficient.real), a21=[a21.real, a21.imag],
                linear=L, action_metric=G, elliptic_basis=S, rotation=R,
                graph=np.array([p.c for p in h]), center=np.array([p.c for p in center]),
                darboux=np.array([p.c for p in D]), canonical_center=np.array([p.c for p in canonical]),
                quadratic_change=np.array([p.c for p in hnf]),
                homological_conditions=conditions, homological_residuals=residuals,
                graph_residual=graph_residual)


def symplectic_residual(coefficients):
    F = map_polynomials(coefficients)
    derivative = [[p.derivative(j) for j in range(4)] for p in F]
    absolute, scale = 0., 0.
    mask = ring(4).degrees <= 2
    for i in range(4):
        for j in range(4):
            value, bound = ring(4).constant(), ring(4).constant()
            for a in range(4):
                for b in range(4):
                    if J4[a, b]:
                        value = value+J4[a, b]*derivative[a][i]*derivative[b][j]
                        bound = bound+Jet(ring(4), abs(derivative[a][i].c))*Jet(ring(4), abs(derivative[b][j].c))
            absolute = max(absolute, float(np.max(abs((value-J4[i, j]).c[mask]))))
            scale = max(scale, float(np.max(bound.c[mask])))
    return absolute, absolute/(1+scale)


def initial_full(z):
    A, u, x = 1+z[0], -z[1]/6, z[2]
    v = z[3]/A**2
    M = np.diag(np.exp(2*x*np.diag(B0)))
    y = full.pack(A, u, np.zeros(4), np.zeros(4), M, v*B0)
    energy = full.constraints(y)['residual'][0]
    if energy >= 0:
        raise ValueError('no negative scalar root')
    y[7] = -np.sqrt(-2*energy)
    return y


def canonical_state(y):
    A, u, q, p, M, L = full.unpack(y)
    x = np.log(np.diag(M))@np.diag(B0)/2
    return np.array([A-1, -6*u, x, A*A*np.trace(L@B0)])


def direct_return(z):
    y0 = initial_full(z)
    def event(t, y):
        return y[3]
    event.direction = -1
    sol = solve_ivp(lambda t, y: full.conformal_rhs(y), (0., 4.), y0,
                    events=event, dense_output=True, method='DOP853',
                    rtol=2e-12, atol=2e-14, max_step=.01)
    # Retain the initial event but select the first later negative crossing.
    events = sol.t_events[0][sol.t_events[0] > 1e-8]
    if not sol.success or len(events) != 1:
        raise ArithmeticError('full return failed or ambiguous')
    time = events[0]
    state = sol.sol(time)
    samples = sol.sol(np.linspace(0., time, 257)).T
    residual = max(float(np.max(full.constraints(y)['normalized'])) for y in samples)
    chart = min(min(y[0], full.ingredients(y)[0], np.linalg.eigvalsh(full.unpack(y)[4]).min()) for y in samples)
    initial = [1+z[0], -z[1]/6, z[2], z[3]/(1+z[0])**2, 0.]
    phase = solve_ivp(phase_rhs, (np.pi/2, 5*np.pi/2), initial, method='DOP853',
                      rtol=2e-12, atol=2e-14, max_step=.01)
    if not phase.success:
        raise ArithmeticError('phase return failed')
    A, u, x, v, eta = phase.y[:, -1]
    exact = canonical_state(state)
    other = np.array([A-1, -6*u, x, A*A*v])
    return dict(initial=y0, returned=state, return_time=float(time), canonical=exact,
                sample_times=np.linspace(0., time, 257), samples=samples,
                phase_canonical=other, phase_time=float(eta),
                formulation_error=float(max(np.linalg.norm(exact-other), abs(time-eta))),
                constraint_max=residual, chart_min=float(chart))
