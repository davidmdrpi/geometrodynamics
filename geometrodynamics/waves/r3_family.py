"""Two-return family of the diagonal circular mode: action and transverse stability.

Section maps of the full homogeneous Einstein-quartet system (unchanged
nonlinear_supported_tt.conformal_rhs, q = (q,0,0,0)) in 12 coordinates
z = (A, p_A, x_1, p_1, ..., x_5, p_5), beta = sum x_k E_k (r3_extension.E),
p_k = H tr(beta' E_k). Shape and its rate are mapped exactly through the
eigenbasis (matrix exp/log and their divided-difference Frechet derivatives),
valid at anisotropic base points. The diagonal subsystem is z[:6].

Linearisations use variational equations: base orbit plus directional
derivatives of the vector field (centred, step 1e-7), integrated with DOP853,
then the section-time correction delta_t = -dq/q', then the derivative of the
section-coordinate map. Multipliers are coordinate-independent.
"""
import numpy as np
from scipy.integrate import solve_ivp
from . import nonlinear_supported_tt as d
from .r3_extension import E

TOL = (1e-12, 1e-14)
H_DIR = 1e-7


def _eig_frechet(B, D, f, df):
    """Frechet derivative of the matrix function f at symmetric B in direction D."""
    w, V = np.linalg.eigh(B)
    Dt = V.T @ D @ V
    fw = f(w)
    G = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            G[i, j] = df(w[i]) if abs(w[i]-w[j]) < 1e-12 else (fw[i]-fw[j])/(w[i]-w[j])
    return V @ (G*Dt) @ V.T, V @ np.diag(fw) @ V.T


def to_state(z):
    z = np.asarray(z, float)
    A, pA = z[0], z[1]
    beta = np.einsum('k,kij->ij', z[2::2], E)
    H = A*A
    bdot = np.einsum('k,kij->ij', z[3::2]/H, E)
    Md, M = _eig_frechet(2*beta, 2*bdot, np.exp, np.exp)
    L = np.linalg.solve(M, Md)/2
    y = d.pack(A, -pA/6, np.zeros(4), np.zeros(4), M, L)
    E0 = d.constraints(y)['residual'][0]
    base = E0-y[7]**2/2
    if base >= 0:
        raise ValueError('no real clock velocity')
    y[7] = -np.sqrt(-2*base)
    return y


def to_section(y):
    A, Ap, q, qp, M, L = d.unpack(y)
    Ms = (M+M.T)/2
    Mdot = M @ L
    Mdot = Mdot+Mdot.T
    dlog, logM = _eig_frechet(Ms, Mdot, np.log, lambda w: 1/w)
    beta, bdot = logM/2, dlog/2
    H = A*A-q @ q/6
    out = [A, -6*Ap]
    for k in range(5):
        out += [np.sum(beta*E[k]), H*np.sum(bdot*E[k])]
    return np.array(out)


def _clock(t, y):
    return y[3]
_clock.direction = -1


def P(z, method='DOP853', full=False, max_step=.05):
    """One clock return of the 12-coordinate section map."""
    y0 = to_state(z)
    s = solve_ivp(lambda t, y: d.conformal_rhs(y), (0., np.pi+.8), y0, method=method, rtol=TOL[0], atol=TOL[1],
                  events=[_clock], max_step=max_step)
    if not s.success:
        raise ArithmeticError(s.message)
    hits = [k for k, t in enumerate(s.t_events[0]) if t > 1.]
    if not hits:
        raise ArithmeticError('no clock section')
    ye = s.y_events[0][hits[0]]
    if full:
        return to_section(ye), ye, float(s.t_events[0][hits[0]]), float(abs(d.constraints(ye)['residual'][0]))
    return to_section(ye)


def _fd(fun, x, h):
    cols = []
    for k in range(len(x)):
        e = np.zeros(len(x))
        e[k] = h
        cols.append((fun(x+e)-fun(x-e))/(2*h))
    return np.array(cols).T


def _reduced(y):
    A, Ap, q, qp, M, L = d.unpack(y)
    return np.r_[A, Ap, q[0], qp[0], M.ravel(), L.ravel()]


def _full(r):
    return d.pack(r[0], r[1], np.r_[r[2], 0, 0, 0], np.r_[r[3], 0, 0, 0], r[4:13].reshape(3, 3), r[13:22].reshape(3, 3))


def _state_jet_flow(r0, T, steps, lead=16):
    """Order-1 jets of the 22-component reduced state through one clock return
    (fixed-step RK4 on r3_extension.rhs, which equals conformal_rhs for q=(q,0,0,0));
    the last `lead` steps use a jet-valued step solving q = 0 (return-time correction)."""
    from .jets import variables
    from . import r3_extension as rx
    v = variables(r0, 1)
    pack = lambda u: [u[0], u[1], u[2], u[3], [u[4:7], u[7:10], u[10:13]], [u[13:16], u[16:19], u[19:22]]]
    unpack = lambda y: [y[0], y[1], y[2], y[3]]+[e for row in y[4] for e in row]+[e for row in y[5] for e in row]
    h = T/steps
    y = rx.rk4(pack(v), h, steps-lead)
    s = v[0]*0+lead*h
    for _ in range(4):
        yc = rx.rk4(y, s/lead, lead)
        s = s-yc[2]/yc[3]
    yc = unpack(rx.rk4(y, s/lead, lead))
    return np.array([j.c[0] for j in yc]), np.array([j.c[1:] for j in yc])


def DP(z, dims=12, steps=(2048, 4096)):
    """(P(z), DP(z)) on the first `dims` section coordinates. The state-space
    linearisation uses exact order-1 jets (Richardson-extrapolated RK4); the
    section-coordinate maps are differentiated by centred differences (h=1e-6)."""
    z = np.asarray(z, float)
    emb = lambda u: np.r_[u, np.zeros(12-dims)]
    zf = emb(z[:dims])
    y0 = to_state(zf)
    T = P(zf, full=True)[2]
    r0 = _reduced(y0)
    out = [_state_jet_flow(r0, T, n) for n in steps]
    base = out[1][0]+(out[1][0]-out[0][0])/15
    Phi = out[1][1]+(out[1][1]-out[0][1])/15
    Pin = _fd(lambda u: _reduced(to_state(emb(u))), zf[:dims], 1e-6)
    Pout = _fd(lambda r: to_section(_full(r))[:dims], base, 1e-6)
    return to_section(_full(base))[:dims], Pout @ Phi @ Pin


def two_return(z0, z1, dims):
    p0, J0 = DP(z0, dims)
    p1, J1 = DP(z1, dims)
    F = np.r_[p0-z1, p1-z0]
    J = np.block([[J0, -np.eye(dims)], [-np.eye(dims), J1]])
    return F, J, J1 @ J0


def correct(v, tangent=None, v_pred=None, dims=6, tol=1e-12, maxit=8):
    """Newton for the two-node closure, with an optional pseudo-arclength row."""
    for it in range(maxit):
        F, J, _ = two_return(v[:dims], v[dims:], dims)
        if tangent is not None:
            F = np.r_[F, tangent @ (v-v_pred)]
            J = np.vstack([J, tangent])
        if np.abs(F).max() < tol:
            break
        v = v+np.linalg.lstsq(J, -F, rcond=None)[0]
    F, J, M2 = two_return(v[:dims], v[dims:], dims)
    return v, float(np.abs(F).max()), J, M2, it+1


def null_tangent(J, prev=None):
    t = np.linalg.svd(J)[2][-1]
    if prev is not None and t @ prev < 0:
        t = -t
    return t


def loop_action(Z):
    """(1/2pi)|closed integral p_A dA + p_x dx + p_y dy| over ordered closed points Z (n x 6),
    by a periodic cubic spline in cumulative chord length."""
    from scipy.interpolate import CubicSpline
    Zc = np.vstack([Z, Z[:1]])
    s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(Zc, axis=0), axis=1))]
    S = CubicSpline(s, Zc, bc_type='periodic')
    dS = S.derivative()
    u = np.linspace(0, s[-1], 20001)
    X, dX = S(u), dS(u)
    integrand = X[:, 1]*dX[:, 0]+X[:, 3]*dX[:, 2]+X[:, 5]*dX[:, 4]
    from scipy.integrate import simpson
    return abs(simpson(integrand, x=u))/(2*np.pi)


def classify_pairs(mults, tol=1e-6):
    """Group multipliers into reciprocal pairs and classify each."""
    m = list(mults)
    out = []
    while m:
        a = m.pop(0)
        j = int(np.argmin([abs(a*b-1) for b in m]))
        b = m.pop(j)
        tr = (a+b).real
        if abs(abs(a)-1) < tol and abs(abs(b)-1) < tol and abs(tr) < 2-tol:
            kind = 'ELLIPTIC'
        elif abs(tr) > 2+tol:
            kind = 'HYPERBOLIC'
        else:
            kind = 'MARGINAL'
        out.append(dict(pair=[[a.real, a.imag], [b.real, b.imag]], trace=float(tr), kind=kind))
    return out
