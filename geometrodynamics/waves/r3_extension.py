"""R3 extension: all five homogeneous n=2 tensor polarisations, and large amplitude.

Part B (polarisations). Full homogeneous matrix system (nonlinear_supported_tt
equations, q = (q,0,0,0)) with shape M = exp(2 beta), beta = sum_k x_k e_k over
an orthonormal STF basis. Section q = 0, q' < 0, Hamiltonian constraint solved
for q'. Section coordinates z = (A, p_A, x_1, p_1, ..., x_5, p_5) with
p_k = H tr(beta' e_k); these are canonical at the fixed point (the leading twist
needs only that, section 3 of the specification). The five elliptic pairs are
1:1-resonant by SO(3) symmetry, so the cubic normal form is a resonant quartic
form F(w) on C^5; the phase-advance shift of polarisation w is
nu(w) = F(w)/(2 pi c |w|^4) with c = |omega(q, qbar)|.

The momentum constraint is not imposed on the tensor data: compensating it
needs a matter current of order eps^2, whose back-reaction on (A, q, M, L)
enters only through |q|^2 and q.q' at order eps^4 (specification, section 3).

Part A (large amplitude) reuses r3_return_map.esu_map with a parallel
finite-difference Jacobian.
"""
from multiprocessing import Pool
import numpy as np
from .jets import Jet, variables, linear_part
from . import r3_return_map as rm

S2, S6 = np.sqrt(2.), np.sqrt(6.)


def stf_basis():
    E = np.zeros((5, 3, 3))
    E[0] = np.diag([1., 1., -2.])/S6
    E[1] = np.diag([1., -1., 0.])/S2
    for k, (i, j) in enumerate(((0, 1), (0, 2), (1, 2))):
        E[2+k, i, j] = E[2+k, j, i] = 1/S2
    return E


E = stf_basis()


# ---------------------------------------------------------------- matrices of jets
def mm(X, Y):
    return [[sum((X[i][k]*Y[k][j] for k in range(1, 3)), X[i][0]*Y[0][j]) for j in range(3)] for i in range(3)]


def madd(X, Y, a=1., b=1.):
    return [[X[i][j]*a+Y[i][j]*b for j in range(3)] for i in range(3)]


def mscale(X, a):
    return [[X[i][j]*a for j in range(3)] for i in range(3)]


def tr(X):
    return X[0][0]+X[1][1]+X[2][2]


def ident(like):
    return [[like*0+(1. if i == j else 0.) for j in range(3)] for i in range(3)]


def inv3(X):
    a, b, c = X[0]
    d, e, f = X[1]
    g, h, i = X[2]
    C = [[e*i-f*h, c*h-b*i, b*f-c*e], [f*g-d*i, a*i-c*g, c*d-a*f], [d*h-e*g, b*g-a*h, a*e-b*d]]
    det = a*C[0][0]+b*C[1][0]+c*C[2][0]
    rd = 1/det
    return [[C[r][s]*rd for s in range(3)] for r in range(3)]


def from_coords(xs, like):
    return [[sum((xs[k]*E[k, i, j] for k in range(1, 5)), xs[0]*E[0, i, j]) for j in range(3)] for i in range(3)]


def coords(X):
    return [sum((X[i][j]*E[k, i, j] for i in range(3) for j in range(3) if E[k, i, j]), X[0][0]*0) for k in range(5)]


# ---------------------------------------------------------------- dynamics
def rhs(y):
    A, Ap, q, qp, M, L = y
    inv = inv3(M)
    Q = q*q
    H = A*A-Q/6
    Hp = 2*A*Ap-q*qp/3
    MM = mm(M, M)
    trinv = tr(inv)
    r = 2*(2*trinv-tr(MM))
    ell = tr(mm(L, L))
    raw = madd(inv, MM, -4+Q/H, -4.)
    t3 = tr(raw)/3
    force = [[raw[i][j]-(t3 if i == j else 0.) for j in range(3)] for i in range(3)]
    App = -A*(r+ell)/6+A*A*A
    qpp = -(trinv+(r+ell)/6)*q
    Lp = madd(L, force, -(Hp/H), 1.)
    Mp = mscale(mm(M, L), 2.)
    return [Ap, App, qp, qpp, Mp, Lp]


def constraint(y):
    A, Ap, q, qp, M, L = y
    inv = inv3(M)
    H = A*A-q*q/6
    r = 2*(2*tr(inv)-tr(mm(M, M)))
    ell = tr(mm(L, L))
    return -3*Ap*Ap+qp*qp/2+H*ell/2-H*r/2+q*q*tr(inv)/2+3*A**4/2


def _axpy(y, k, h):
    out = []
    for a, b in zip(y, k):
        out.append(madd(a, b, 1., h) if isinstance(a, list) else a+b*h)
    return out


def rk4(y, h, n):
    for _ in range(n):
        k1 = rhs(y)
        k2 = rhs(_axpy(y, k1, h/2))
        k3 = rhs(_axpy(y, k2, h/2))
        k4 = rhs(_axpy(y, k3, h))
        y = [madd(a, madd(madd(b, c, 1., 2.), madd(e, f, 2., 1.)), 1., h/6)
             if isinstance(a, list) else a+(b+2*c+2*e+f)*(h/6)
             for a, b, c, e, f in zip(y, k1, k2, k3, k4)]
    return y


def section_to_state(z):
    A, pA = z[0], z[1]
    xs, ps = z[2::2], z[3::2]
    H = A*A
    beta = from_coords(xs, A)
    bdot = from_coords([p/H for p in ps], A)
    I = ident(A)
    b2 = mm(beta, beta)
    M = madd(madd(I, beta, 1., 2.), madd(b2, mm(b2, beta), 2., 4/3), 1., 1.)
    # d/dt exp(2 beta) through order three (beta, beta' have zero constant part)
    Md = madd(madd(bdot, madd(mm(bdot, beta), mm(beta, bdot)), 2., 2.),
              madd(madd(mm(mm(bdot, beta), beta), mm(mm(beta, bdot), beta)), mm(b2, bdot)), 1., 4/3)
    L = mscale(mm(inv3(M), Md), .5)
    Ap = -pA/6
    y = [A, Ap, A*0, A*0, M, L]
    E0 = constraint(y)
    qp2 = 2*(y[3]*y[3]/2-E0)
    y[3] = -(qp2.sqrt() if isinstance(qp2, Jet) else np.sqrt(qp2))
    return y


def state_to_section(y):
    A, Ap, q, qp, M, L = y
    I = ident(A)
    m = madd(M, I, 1., -1.)
    m2 = mm(m, m)
    beta = mscale(madd(madd(m, m2, 1., -.5), mm(m2, m), 1., 1/3), .5)
    md = mscale(mm(M, L), 2.)
    dlog = madd(madd(md, madd(mm(m, md), mm(md, m)), 1., -.5),
                madd(madd(mm(m2, md), mm(mm(m, md), m)), mm(md, m2)), 1., 1/3)
    H = A*A-q*q/6
    xs = coords(beta)
    ps = [c*.5*H for c in coords(dlog)]
    out = [A, -6*Ap]
    for x, p in zip(xs, ps):
        out += [x, p]
    return out


def jet_return_map(order, steps, lead=16):
    z = variables(np.r_[1., np.zeros(11)], order)
    h = np.pi/steps
    y = rk4(section_to_state(z), h, steps-lead)
    s = z[0]*0+lead*h
    for _ in range(order+3):
        yc = rk4(y, s/lead, lead)
        s = s-yc[2]/yc[3]
    yc = rk4(y, s/lead, lead)
    return dict(P=state_to_section(yc), time=s+(steps-lead)*h, constraint=constraint(yc), q=yc[2])


# ---------------------------------------------------------------- multi-mode normal form
def multimode_normal_form(P, point, hyp=(0, 1), tol_block=1e-9):
    """Order-3 normal form with one hyperbolic pair and m identical elliptic pairs.

    Returns the resonant cubic tensor T[d, a, b, c] (w_d component, monomial
    w_a w_b wbar_c, symmetrised in a,b) divided by mu, and the normalisation c.
    """
    n = len(P)
    m = (n-2)//2
    order = P[0].s['order']
    F = [Jet(p.c.copy(), n, order) for p in P]
    fixed_error = max(abs(F[i].c[0]-point[i]) for i in range(n))
    for i in range(n):
        F[i].c[0] = 0.
    A = linear_part(F)
    blocks = [hyp]+[(2+2*k, 3+2*k) for k in range(m)]
    mask = np.zeros_like(A, dtype=bool)
    for b in blocks:
        mask[np.ix_(b, b)] = True
    off_block = float(np.abs(A[~mask]).max())
    lh, vh = np.linalg.eig(A[np.ix_(hyp, hyp)])
    order_h = np.argsort(-abs(lh))
    B0 = A[2:4, 2:4]
    le, ve = np.linalg.eig(B0)
    ie = int(np.argmax(le.imag))
    mu, qv = le[ie], ve[:, ie]
    block_spread = max(float(np.abs(A[2+2*k:4+2*k, 2+2*k:4+2*k]-B0).max()) for k in range(m))
    V = np.zeros((n, n), complex)
    L = np.zeros(n, complex)
    for j, ih in enumerate(order_h):
        V[list(hyp), j] = vh[:, ih]
        L[j] = lh[ih]
    for k in range(m):
        V[2+2*k:4+2*k, 2+k] = qv
        V[2+2*k:4+2*k, 2+m+k] = qv.conj()
        L[2+k], L[2+m+k] = mu, mu.conjugate()
    Vi = np.linalg.inv(V)
    xi = variables(np.zeros(n), order, complex)
    zero = xi[0]*0
    delta = [sum((xi[j]*V[i, j] for j in range(n) if V[i, j] != 0), zero) for i in range(n)]
    Fd = [f.compose(delta) for f in F]
    G = [sum((Fd[j]*Vi[i, j] for j in range(n) if Vi[i, j] != 0), zero) for i in range(n)]
    s = xi[0].s
    h, rel_div = [], []
    for i in range(n):
        c = np.zeros(s['n'], complex)
        for mo, e in enumerate(s['exps']):
            if s['deg'][mo] == 2:
                lam_a = np.prod(L**e)
                den = lam_a-L[i]
                rel_div.append(abs(den)/max(abs(lam_a), abs(L[i])))
                c[mo] = G[i].c[mo]/den
        h.append(Jet(c, n, order))
    Phi = [xi[i]+h[i] for i in range(n)]
    GPhi = [g.compose(Phi) for g in G]
    inv = list(xi)
    for _ in range(order):
        inv = [xi[i]-h[i].compose(inv) for i in range(n)]
    Gn = [g.compose(GPhi) for g in inv]
    T = np.zeros((m, m, m, m), complex)
    for d in range(m):
        for mo, e in enumerate(s['exps']):
            if s['deg'][mo] != 3:
                continue
            if e[:2].any():
                continue
            wa, wb = e[2:2+m], e[2+m:]
            if wa.sum() == 2 and wb.sum() == 1:
                a_idx = np.repeat(np.arange(m), wa)
                c = int(np.argmax(wb))
                val = Gn[2+d].c[mo]/mu
                if a_idx[0] == a_idx[1]:
                    T[d, a_idx[0], a_idx[0], c] += val
                else:
                    T[d, a_idx[0], a_idx[1], c] += val/2
                    T[d, a_idx[1], a_idx[0], c] += val/2
    Om = np.zeros((2, 2))
    Om[1, 0], Om[0, 1] = 1., -1.
    cnorm = float(abs(qv @ Om @ qv.conj()))
    return dict(T=T, mu=complex(mu), theta=float(np.angle(mu)), lam=complex(L[0]), c=cnorm,
                fixed_error=float(fixed_error), off_block=off_block, block_spread=block_spread,
                min_relative_divisor=float(min(rel_div)), resonance_distances=[float(abs(mu**k-1)) for k in range(1, 5)])


def rayleigh(T, w):
    """<w, G(w)> with G_d(w) = sum T[d,a,b,c] w_a w_b conj(w_c)."""
    return np.einsum('d,dabc,a,b,c->', w.conj(), T, w, w, w.conj())


def quartic_matrices(T):
    """K on C^(m*m) with <w,G(w)> = (w x w)^H K (w x w), Hermitian parts split."""
    m = T.shape[0]
    K = np.zeros((m*m, m*m), complex)
    for c in range(m):
        for d in range(m):
            for a in range(m):
                for b in range(m):
                    K[c*m+d, a*m+b] = T[d, a, b, c]
    # symmetric-tensor projector
    Psym = np.zeros((m*m, m*m))
    for a in range(m):
        for b in range(m):
            Psym[a*m+b, a*m+b] += .5
            Psym[a*m+b, b*m+a] += .5
    Ks = Psym @ K @ Psym
    Kim = (Ks-Ks.conj().T)/2j
    Kre = (Ks+Ks.conj().T)/2
    return Kim, Kre, Psym


def nu_of(T, w, c):
    w = np.asarray(w, complex)
    w = w/np.linalg.norm(w)
    return float(rayleigh(T, w).imag/(2*np.pi*c))


def sup_nu(T, c, restarts=64, seed=11):
    """Numerical max and min of nu(w) over the unit sphere of C^m, plus the
    Sym^2 eigenvalue bound (an upper bound on the max and lower bound on the min)."""
    from scipy.optimize import minimize
    m = T.shape[0]
    rng = np.random.default_rng(seed)

    def f(v, sgn):
        w = v[:m]+1j*v[m:]
        return -sgn*nu_of(T, w, c)
    out = {}
    for sgn, name in ((1, 'max'), (-1, 'min')):
        best = None
        for _ in range(restarts):
            v0 = rng.normal(size=2*m)
            r = minimize(f, v0, args=(sgn,), method='BFGS', options=dict(gtol=1e-12))
            val = -sgn*r.fun
            if best is None or sgn*val > sgn*best[0]:
                best = (val, r.x)
        w = best[1][:m]+1j*best[1][m:]
        out[name] = float(best[0])
        out[name+'_polarisation'] = [[float(x.real), float(x.imag)] for x in w/np.linalg.norm(w)]
    Kim, Kre, P = quartic_matrices(T)
    ev = np.linalg.eigvalsh(P @ Kim @ P)
    out['sym2_upper_bound'] = float(ev.max()/(2*np.pi*c))
    out['sym2_lower_bound'] = float(ev.min()/(2*np.pi*c))
    out['dissipative_max'] = float(np.abs(np.linalg.eigvalsh(P @ Kre @ P)).max())
    out['twist_scale'] = float(np.abs(ev).max())
    return out


def so3_rep(R):
    """5x5 matrix of W -> R W R^T in the STF basis E."""
    D = np.zeros((5, 5))
    for k in range(5):
        X = R @ E[k] @ R.T
        for l in range(5):
            D[l, k] = np.sum(X*E[l])
    return D


# ---------------------------------------------------------------- Part A: circles, parallel FD
def _P(z):
    return rm.esu_map(np.asarray(z))[0]


def _fd_point(z, h=1e-6):
    return rm.fd_jacobian(_P, np.asarray(z), h)


def invariant_circle_parallel(a, K0, omega0, M, pool, tol=1e-11, maxit=15):
    K, omega = K0.copy(), omega0
    th = 2*np.pi*np.arange(M)/M
    e1 = np.exp(-1j*th)/M
    hist = []
    for it in range(maxit):
        PK = np.array(pool.map(_P, list(K)))
        T, dT = rm.shift_matrix(M, omega)
        R = PK-T @ K
        c1 = e1 @ K[:, 2]
        res = np.r_[R.ravel(), c1.real-a/2, c1.imag]
        hist.append(float(np.max(abs(res))))
        if not np.isfinite(hist[-1]):
            raise ArithmeticError('nonfinite residual')
        if hist[-1] < tol:
            break
        J = np.zeros((4*M+2, 4*M+1))
        blocks = pool.map(_fd_point, list(K))
        for j in range(M):
            J[4*j:4*j+4, 4*j:4*j+4] = blocks[j]
        J[:4*M, :4*M] -= np.kron(T, np.eye(4))
        J[:4*M, 4*M] = -(dT @ K).ravel()
        J[4*M, 2:4*M:4] = e1.real
        J[4*M+1, 2:4*M:4] = e1.imag
        step = np.linalg.lstsq(J, -res, rcond=None)[0]
        K = K+step[:4*M].reshape(M, 4)
        omega = omega+step[4*M]
    PK = np.array(pool.map(_P, list(K)))
    T, _ = rm.shift_matrix(M, omega)
    Ck = np.abs(np.fft.fft(K-K.mean(0), axis=0))/M
    kk = np.abs(np.fft.fftfreq(M, 1/M))
    return dict(K=K, omega=float(omega), residual=float(np.max(abs(PK-T @ K))), iterations=it+1,
                history=hist, action=abs(rm.action(K)), fourier_tail=float(Ck[kk >= M//2-2].max()))
