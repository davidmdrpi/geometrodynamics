"""Vacuum CMC data on a compact S2 mapping torus; public freeze 2b1bfca.

This solves initial constraints, not wave evolution or a measurement model.
The scalar solve uses Fourier x Legendre collocation; physical constraint
verification below constructs coordinate metric jets independently.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from numpy.polynomial import legendre as leg
from scipy.sparse.linalg import LinearOperator, gmres

PREREG = '2b1bfca7b2650db7f8dcf74ca1c54ecd8bdffa17'
BASELINE = '4cd86541d3b836b35561b0c4a3a54629d28851cd'
AMPLITUDES = (0., .02, .05, .1, .2)
GRIDS = ((24, 12), (40, 20), (64, 32))
POINTS = tuple((s, t, .37) for s in (.23, .77, 1.39, 2.17, 2.81, 4.13)
               for t in (.41, .83, 1.21, 1.87, 2.39))


def radial(s, k=.5, correction=1.):
    """p, w, p+w', and the remaining divergence coefficient."""
    s = np.asarray(s)
    p = np.sin(k*s)
    w = correction*k*np.cos(k*s)/(k*k+4)
    P = p-correction*k*k*np.sin(k*s)/(k*k+4)
    residual = (1-correction)*k*np.cos(k*s)
    return p, w, P, residual


def norm_squared(s, u, epsilon, C=.5, k=.5):
    _, w, P, _ = radial(s, k)
    return 6*C*C + epsilon**2*(18*P**2*u**2*(1-u*u)
                                +18*w**2*(1-u*u)**2)


def conformal_tensor(s, theta, epsilon, C=.5, k=.5, correction=1.):
    sn, u = np.sin(theta), np.cos(theta)
    _, w, P, _ = radial(s, k, correction)
    A = np.diag([2*C, -C, -C*sn*sn])
    A[0, 2] = A[2, 0] = epsilon*P*3*u*sn*sn
    A[1, 2] = A[2, 1] = -epsilon*w*3*sn**3
    return A


@dataclass
class Grid:
    ns: int
    nu: int

    def __post_init__(self):
        if self.ns < 8 or self.ns % 2 or self.nu < 4:
            raise ValueError('even ns>=8 and nu>=4 required')
        self.s = 2*np.pi*np.arange(self.ns)/self.ns
        self.u, self.weights = leg.leggauss(self.nu)
        self.V = leg.legvander(self.u, self.nu-1)
        self.Vi = np.linalg.inv(self.V)
        self.freq = np.fft.fftfreq(self.ns, 1/self.ns)
        self.eigen = -self.freq[:, None]**2 - np.arange(self.nu)[None, :]*(np.arange(self.nu)[None, :]+1)

    def coefficients(self, values):
        return np.fft.fft(values, axis=0)/self.ns @ self.Vi.T

    def from_coefficients(self, coeff):
        return np.fft.ifft(coeff @ self.V.T*self.ns, axis=0).real

    def laplacian(self, values):
        return self.from_coefficients(self.coefficients(values)*self.eigen)

    def inverse_reference(self, values):
        return self.from_coefficients(self.coefficients(values)/(self.eigen-1.25))

    def evaluate(self, values, s, u, ds=0, du=0):
        return evaluate_coefficients(self.coefficients(values), self.freq, s, u, ds, du)


def evaluate_coefficients(coeff, freq, s, u, ds=0, du=0):
    c = coeff*(1j*freq[:, None])**ds
    if du:
        c = leg.legder(c, m=du, axis=1)
    # Cartesian product evaluation; scalar arguments return a (1,1) matrix.
    es = np.exp(1j*np.atleast_1d(s)[:, None]*freq[None, :])
    vu = leg.legvander(np.atleast_1d(u), c.shape[1]-1)
    return (es @ c @ vu.T).real


def equation(grid, psi, epsilon, C=.5, k=.5, Lambda=.25):
    a2 = norm_squared(grid.s[:, None], grid.u[None, :], epsilon, C, k)
    return grid.laplacian(psi)-psi/4+a2*psi**-7/8+Lambda*psi**5/4


def solve(ns, nu, epsilon, C=.5, k=.5, Lambda=.25):
    if not np.isfinite([epsilon, C, k, Lambda]).all() or k <= 0:
        raise ValueError('finite parameters and positive k required')
    grid = Grid(ns, nu)
    psi = np.ones((ns, nu))
    a2 = norm_squared(grid.s[:, None], grid.u[None, :], epsilon, C, k)
    history = []
    for iteration in range(30):
        residual = equation(grid, psi, epsilon, C, k, Lambda)
        error = float(np.max(abs(residual)))
        history.append(error)
        if error < 1e-11:
            break
        diagonal = -.25-7*a2*psi**-8/8+5*Lambda*psi**4/4
        shape = psi.shape
        op = LinearOperator((psi.size, psi.size), matvec=lambda x:
                            (grid.laplacian(x.reshape(shape))+diagonal*x.reshape(shape)).ravel())
        pre = LinearOperator(op.shape, matvec=lambda x:
                             grid.inverse_reference(x.reshape(shape)).ravel())
        delta, info = gmres(op, -residual.ravel(), M=pre, rtol=1e-12,
                            atol=1e-14, maxiter=100)
        if info:
            raise ArithmeticError(f'linear solve failed: {info}')
        delta = delta.reshape(shape)
        step = 1.
        for backtrack in range(30):
            trial = psi+step*delta
            if trial.min() > 0 and np.max(abs(equation(grid, trial, epsilon, C, k, Lambda))) < error:
                psi = trial
                break
            step /= 2
        else:
            raise ArithmeticError('positive Newton line search failed')
    else:
        raise ArithmeticError('Newton did not converge in 30 iterations')
    return dict(ns=ns, nu=nu, epsilon=epsilon, C=C, k=k, Lambda=Lambda,
                start='psi=1', psi=psi.tolist(), iterations=iteration,
                residual_history=history)


class Interpolant:
    def __init__(self, record):
        self.record = record
        self.grid = Grid(record['ns'], record['nu'])
        self.psi = np.asarray(record['psi'], dtype=float)
        if self.psi.shape != (self.grid.ns, self.grid.nu) or not np.isfinite(self.psi).all() or self.psi.min() <= 0:
            raise ValueError('invalid scalar solution')
        self.coeff = self.grid.coefficients(self.psi)

    def values(self, s, u, ds=0, du=0):
        return evaluate_coefficients(self.coeff, self.grid.freq, s, u, ds, du)

    def point(self, s, theta, ds=0, du=0):
        return float(self.values(s, np.cos(theta), ds, du)[0, 0])

    def metric_tensor(self, x, wrong_weight=False, correction=1.):
        s, theta, _ = x
        psi = self.point(s, theta)
        g = psi**4*np.diag([1., 1., np.sin(theta)**2])
        r = self.record
        A = conformal_tensor(s, theta, r['epsilon'], r['C'], r['k'], correction)
        K = psi**(2 if wrong_weight else -2)*A
        return g, K

    def offgrid_residual(self, ns=96, nu=48):
        target = Grid(ns, nu)
        s, u = target.s, target.u
        psi = self.values(s, u)
        lap = self.values(s, u, ds=2)+(1-u*u)[None, :]*self.values(s, u, du=2)-2*u[None, :]*self.values(s, u, du=1)
        r = self.record
        a2 = norm_squared(s[:, None], u[None, :], r['epsilon'], r['C'], r['k'])
        return lap-psi/4+a2*psi**-7/8+r['Lambda']*psi**5/4

    def section(self, s):
        u, weights = leg.leggauss(64)
        psi = self.values(s, u)[0]
        ps = self.values(s, u, ds=1)[0]
        mean_curvature = 4*ps/psi**3
        trace_surface = -2*self.record['C']*psi**-6
        return dict(s=float(s), area=float(2*np.pi*np.sum(weights*psi**4)),
                    H=mean_curvature.tolist(), theta_plus=(mean_curvature-trace_surface).tolist(),
                    theta_minus=(-mean_curvature-trace_surface).tolist())


def coordinate_constraints(interpolant, x, h, wrong_weight=False, correction=1.):
    """Direct g/K differences -> Christoffel/Ricci/divergence, no conformal PDE."""
    x = np.asarray(x, dtype=float)
    g, K = interpolant.metric_tensor(x, wrong_weight, correction)
    inverse = np.linalg.inv(g)
    dg, dK = np.zeros((3, 3, 3)), np.zeros((3, 3, 3))
    ddg = np.zeros((3, 3, 3, 3))
    axes = np.eye(3)*h
    for a in range(3):
        gp, kp = interpolant.metric_tensor(x+axes[a], wrong_weight, correction)
        gm, km = interpolant.metric_tensor(x-axes[a], wrong_weight, correction)
        dg[a], dK[a] = (gp-gm)/(2*h), (kp-km)/(2*h)
        ddg[a, a] = (gp-2*g+gm)/h**2
        for b in range(a):
            gpp = interpolant.metric_tensor(x+axes[a]+axes[b], wrong_weight, correction)[0]
            gpm = interpolant.metric_tensor(x+axes[a]-axes[b], wrong_weight, correction)[0]
            gmp = interpolant.metric_tensor(x-axes[a]+axes[b], wrong_weight, correction)[0]
            gmm = interpolant.metric_tensor(x-axes[a]-axes[b], wrong_weight, correction)[0]
            ddg[a, b] = ddg[b, a] = (gpp-gpm-gmp+gmm)/(4*h*h)
    G = np.zeros((3, 3, 3))
    dG = np.zeros((3, 3, 3, 3))
    dinv = np.array([-inverse@dg[a]@inverse for a in range(3)])
    for a in range(3):
        for b in range(3):
            for c in range(3):
                terms = dg[b, :, c]+dg[c, :, b]-dg[:, b, c]
                G[a, b, c] = inverse[a]@terms/2
                for d in range(3):
                    dt = ddg[d, b, :, c]+ddg[d, c, :, b]-ddg[d, :, b, c]
                    dG[d, a, b, c] = (dinv[d, a]@terms+inverse[a]@dt)/2
    ricci = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            for a in range(3):
                ricci[i, j] += dG[a, a, i, j]-dG[j, a, i, a]
                for b in range(3):
                    ricci[i, j] += G[a, a, b]*G[b, i, j]-G[a, j, b]*G[b, i, a]
    R = float(np.sum(inverse*ricci))
    mixed = inverse@K
    K2 = float(np.trace(mixed@mixed))
    trK = float(np.trace(mixed))
    dmixed = np.array([dinv[a]@K+inverse@dK[a] for a in range(3)])
    terms = np.zeros((3, 3))
    for i in range(3):
        terms[0, i] = sum(dmixed[j, j, i] for j in range(3))-np.trace(dmixed[i])
        terms[1, i] = sum(G[j, j, k]*mixed[k, i] for j in range(3) for k in range(3))
        terms[2, i] = -sum(G[k, j, i]*mixed[j, k] for j in range(3) for k in range(3))
    momentum = terms.sum(axis=0)
    Lambda = interpolant.record['Lambda']
    H = R+trK*trK-K2-2*Lambda
    return dict(point=x.tolist(), h=h, R=R, K2=K2, trace_K=trK,
                metric=g.tolist(), extrinsic_curvature=K.tolist(),
                hamiltonian=H, momentum=momentum.tolist(),
                hamiltonian_normalized=abs(H)/max(1., abs(R)+trK*trK+K2+2*abs(Lambda)),
                momentum_normalized=float(np.linalg.norm(momentum)/max(1., sum(np.linalg.norm(t) for t in terms))),
                momentum_terms=terms.tolist())


def momentum_fd(n, k=.5):
    s = 4*np.pi*np.arange(n)/n
    h = 4*np.pi/n
    matrix = (np.roll(np.eye(n), 1, axis=1)-2*np.eye(n)+np.roll(np.eye(n), -1, axis=1))/h**2-4*np.eye(n)
    rhs = -k*np.cos(k*s)
    w = np.linalg.solve(matrix, rhs)
    exact = radial(s, k)[1]
    return dict(n=n, s=s.tolist(), rhs=rhs.tolist(), w=w.tolist(),
                exact=exact.tolist(), residual=(matrix@w-rhs).tolist())


def seam_check(interpolant, k=None, twisted=True):
    r = interpolant.record
    k = r['k'] if k is None else k
    jacobian = np.diag([1., -1., 1.]) if twisted else np.eye(3)
    values, derivatives, scalars = [], [], []
    h = 1e-5
    for s, theta, _ in POINTS:
        t2 = np.pi-theta if twisted else theta
        def a(v, t):
            return conformal_tensor(v, t, r['epsilon'], r['C'], k)
        values.append((a(s, theta)-jacobian@a(s+2*np.pi, t2)@jacobian).tolist())
        da = (a(s+h, theta)-a(s-h, theta))/(2*h)
        db = (a(s+2*np.pi+h, t2)-a(s+2*np.pi-h, t2))/(2*h)
        derivatives.append((da-jacobian@db@jacobian).tolist())
        scalars.append(interpolant.point(s, theta)-interpolant.point(s+2*np.pi, t2))
    return dict(value_residuals=values, derivative_residuals=derivatives, scalar_residuals=scalars)


@lru_cache(None)
def symbolic_checks():
    import sympy as S
    t, s, C, e, k = S.symbols('theta s C epsilon k', real=True)
    sphere = S.diag(1, S.sin(t)**2)
    inv = sphere.inv()
    coords = (t, S.Symbol('phi', real=True))
    connection = [[[sum(inv[a,d]*(S.diff(sphere[d,c],coords[b])+S.diff(sphere[d,b],coords[c])-S.diff(sphere[b,c],coords[d]))/2 for d in range(2)) for c in range(2)] for b in range(2)] for a in range(2)]
    v = S.Matrix([0, 3*S.cos(t)*S.sin(t)**2])
    dv = S.Matrix(2,2,lambda a,b:S.diff(v[b],coords[a])-sum(connection[c][a][b]*v[c] for c in range(2)))
    tensor = dv+dv.T
    div = S.Matrix([sum(inv[a,b]*(S.diff(tensor[b,c],coords[a])-sum(connection[d][a][b]*tensor[d,c]+connection[d][a][c]*tensor[b,d] for d in range(2))) for a in range(2) for b in range(2)) for c in range(2)])
    w = k*S.cos(k*s)/(k*k+4)
    residuals = [S.trace(inv*dv), S.trace(inv*tensor), *(div+4*v),
                 (v.T*inv*v)[0]-9*S.cos(t)**2*S.sin(t)**2,
                 S.trace(inv*tensor*inv*tensor)-18*S.sin(t)**4,
                 tensor[0,1]+3*S.sin(t)**3,
                 S.diff(w,s,2)-4*w+k*S.cos(k*s)]
    J=S.diag(-1,1)
    residuals.extend(J*v.subs(t,S.pi-t)+v)
    residuals.extend(J*tensor.subs(t,S.pi-t)*J+tensor)
    return tuple(str(S.trigsimp(S.simplify(x))) for x in residuals)
