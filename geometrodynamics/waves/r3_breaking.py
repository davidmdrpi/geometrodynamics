"""Resonance-breaking scan: continuum of period-q points or Birkhoff chain?

Freeze: docs/r3_breaking_prereg.md.

For a map P with a resonant closed curve c(phi) of rotation p/q, solve at each
node-0 phase phi the square system in X = (z_0, ..., z_{q-1}, lambda):

    P(z_i) - z_{i+1} = 0                  (i = 0..q-2)
    P(z_{q-1}) - z_0 - lambda g = 0
    c'(phi) . (z_0 - c(phi)) = 0          (phase)

with g = OMEGA^T c'(phi) / |.|, the symplectic dual of the curve tangent (it
overlaps the cokernel of the closure Jacobian at a continuum of periodic
points). lambda(phi) = 0 exactly where a period-q orbit passes through the
phase-phi slice. A continuum gives lambda == 0 for every phi; a Birkhoff
chain gives a function with 2q zeros per turn (q elliptic, q hyperbolic).
"""
import numpy as np


def symplectic(n):
    """u^T W v = sum (u_p v_x - u_x v_p) over pairs (x, p), as r3_return_map.OMEGA."""
    W = np.zeros((n, n))
    for k in range(0, n, 2):
        W[k, k+1], W[k+1, k] = -1., 1.
    return W


def trig_curve(K):
    """Trigonometric interpolant of samples K[j] at theta_j = 2 pi j/M."""
    K = np.asarray(K, float)
    M = len(K)
    C = np.fft.fft(K, axis=0)/M
    k = np.fft.fftfreq(M, 1/M)
    if M % 2 == 0:
        C[M//2] /= 2
        C = np.vstack([C, C[M//2:M//2+1]])
        k = np.r_[k, M//2]

    def c(th):
        return (np.exp(1j*np.outer([th], k)) @ C).real[0]

    def dc(th):
        return (np.exp(1j*np.outer([th], k)) @ (1j*k[:, None]*C)).real[0]
    return c, dc


def spline_curve(Z):
    """Periodic cubic spline through closed ordered points Z (last != first), param in [0, 2 pi)."""
    from scipy.interpolate import CubicSpline
    Z = np.asarray(Z, float)
    Zc = np.vstack([Z, Z[:1]])
    s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(Zc, axis=0), axis=1))]
    S = CubicSpline(2*np.pi*s/s[-1], Zc, bc_type='periodic')
    dS = S.derivative()
    return (lambda th: S(np.mod(th, 2*np.pi))), (lambda th: dS(np.mod(th, 2*np.pi)))


def fd_jacobian(P, z, h=1e-7):
    z = np.asarray(z, float)
    cols = []
    for k in range(len(z)):
        e = np.zeros(len(z))
        e[k] = h
        cols.append((np.asarray(P(z+e))-np.asarray(P(z-e)))/(2*h))
    return np.array(cols).T


def _residual(PZ, Z, lam, g, c0, t0):
    q = len(Z)
    F = [PZ[i]-Z[(i+1) % q] for i in range(q)]
    F[-1] = F[-1]-lam*g
    return np.r_[np.concatenate(F), t0 @ (Z[0]-c0)]


def _jacobian(blocks, g, t0):
    q, n = len(blocks), len(g)
    J = np.zeros((q*n+1, q*n+1))
    for i, B in enumerate(blocks):
        J[i*n:(i+1)*n, i*n:(i+1)*n] = B
        j = (i+1) % q
        J[i*n:(i+1)*n, j*n:(j+1)*n] -= np.eye(n)
    J[(q-1)*n:q*n, q*n] = -g
    J[q*n, :n] = t0
    return J


def scan_point(P, seeds, c0, t0, W=None, jac=None, tol=1e-12, maxit=12):
    """Solve the phase-phi system from node seeds. Returns the converged nodes, lambda and diagnostics."""
    Z = [np.asarray(s, float).copy() for s in seeds]
    q, n = len(Z), len(Z[0])
    W = symplectic(n) if W is None else W
    t0 = np.asarray(t0, float)/np.linalg.norm(t0)
    g = W.T @ t0
    g = g/np.linalg.norm(g)
    jac = jac or (lambda z: fd_jacobian(P, z))
    lam, hist = 0., []
    for it in range(maxit):
        PZ = [np.asarray(P(z)) for z in Z]
        F = _residual(PZ, Z, lam, g, c0, t0)
        hist.append(float(np.abs(F).max()))
        if not np.isfinite(hist[-1]):
            raise ArithmeticError('nonfinite residual')
        if hist[-1] < tol:
            break
        J = _jacobian([jac(z) for z in Z], g, t0)
        dX = np.linalg.solve(J, -F)
        Z = [Z[i]+dX[i*n:(i+1)*n] for i in range(q)]
        lam = lam+dX[-1]
    J = _jacobian([jac(z) for z in Z], g, t0)
    sv = np.linalg.svd(J, compute_uv=False)
    return dict(nodes=[z.tolist() for z in Z], lam=float(lam), residual=hist[-1], iterations=it+1,
                history=hist, smin=float(sv[-1]), cond=float(sv[0]/sv[-1]), g=g.tolist(), t0=t0.tolist(),
                c0=np.asarray(c0, float).tolist(), J=J.tolist())


def chord_resolve(P_alt, point, tol=1e-12, maxit=4):
    """Re-converge a scan point with another integrator, keeping the primary Jacobian (chord Newton)."""
    Z = [np.array(z) for z in point['nodes']]
    q, n = len(Z), len(Z[0])
    g, t0, c0 = (np.array(point[k]) for k in ('g', 't0', 'c0'))
    J = np.array(point['J'])
    lam, hist = point['lam'], []
    for _ in range(maxit):
        F = _residual([np.asarray(P_alt(z)) for z in Z], Z, lam, g, c0, t0)
        hist.append(float(np.abs(F).max()))
        if hist[-1] < tol:
            break
        dX = np.linalg.solve(J, -F)
        Z = [Z[i]+dX[i*n:(i+1)*n] for i in range(q)]
        lam = lam+dX[-1]
    return dict(lam=float(lam), residual=hist[-1], history=hist)


def periodic_orbit(P, seeds, jac=None, tol=1e-12, maxit=15):
    """Unconstrained closure P(z_i) = z_{i+1 mod q} (lambda = 0, no phase row), least-squares Newton."""
    Z = [np.asarray(s, float).copy() for s in seeds]
    q, n = len(Z), len(Z[0])
    jac = jac or (lambda z: fd_jacobian(P, z))
    hist = []
    for it in range(maxit):
        F = np.concatenate([np.asarray(P(Z[i]))-Z[(i+1) % q] for i in range(q)])
        hist.append(float(np.abs(F).max()))
        if hist[-1] < tol:
            break
        J = _jacobian([jac(z) for z in Z], np.zeros(n), np.zeros(n))[:q*n, :q*n]
        dX = np.linalg.lstsq(J, -F, rcond=None)[0]
        Z = [Z[i]+dX[i*n:(i+1)*n] for i in range(q)]
    blocks = [jac(z) for z in Z]
    J = _jacobian(blocks, np.zeros(n), np.zeros(n))[:q*n, :q*n]
    return dict(nodes=[z.tolist() for z in Z], residual=hist[-1], iterations=it+1, history=hist,
                smin=float(np.linalg.svd(J, compute_uv=False)[-1]))


def centre_block(blocks, n_unstable=1, n_stable=1, passes=8):
    """Centre block of the product blocks[-1] @ ... @ blocks[0] by periodic orthogonal iteration;
    the hyperbolic directions are never multiplied together."""
    n = blocks[0].shape[0]
    Q = np.linalg.qr(np.random.default_rng(0).normal(size=(n, n)))[0]   # generic start: no exact alignment
    for _ in range(passes):
        Q0, Rs = Q, []
        for B in blocks:
            Q, R = np.linalg.qr(B @ Q)
            sgn = np.sign(np.diag(R))
            sgn[sgn == 0] = 1
            Q, R = Q*sgn, sgn[:, None]*R
            Rs.append(R)
    # prod R = Q^T M Q0; the centre columns of Q and Q0 span the same plane but may be rotated
    mid = slice(n_unstable, n-n_stable)
    C = np.eye(n-n_unstable-n_stable)
    for R in Rs:
        C = R[mid, mid] @ C
    return (Q0.T @ Q)[mid, mid] @ C


def sign_changes(lam):
    """Cyclic sign changes of a periodic sampled function (exact zeros are skipped)."""
    s = np.sign(np.asarray(lam, float))
    s = s[s != 0]
    return int(np.sum(s != np.roll(s, 1))) if len(s) else 0


def classify(lams, noise, q, floor=1e-11, unbroken_cap=1e-9):
    """Registered rule (docs/r3_breaking_prereg.md, section 4); orbit gate applied by the caller."""
    lams = np.asarray(lams, float)
    r = max(10*noise, floor)
    big = float(np.abs(lams).max())
    nsc = sign_changes(lams)
    spec = np.abs(np.fft.rfft(lams))/len(lams)
    k = int(np.argmax(spec[1:])+1)
    out = dict(resolution=r, Lambda=big, sign_changes=nsc, dominant_harmonic=k,
               harmonics=[float(x) for x in spec[:16]])
    if big >= 10*r and nsc >= 2*q and nsc % (2*q) == 0:
        out['label'] = 'BROKEN_CHAIN'
    elif big <= r and r <= unbroken_cap:
        out['label'] = 'UNBROKEN_LOOP'
    else:
        out['label'] = 'INDETERMINATE'
    return out


def toy_map(eps, rho0=.4, nu=-1., lam_u=85.):
    """4D test map: hyperbolic pair (u, s) x area-preserving twist map in (x, p) with a
    q=5 resonant kick of strength eps. I = (x^2+p^2)/2, angle phi = atan2(-p, x)."""
    def P(z):
        u, s, x, p = z
        I, phi = (x*x+p*p)/2, np.arctan2(-p, x)
        I2 = I+eps*np.sin(5*phi)
        phi2 = phi+2*np.pi*(rho0+nu*(I2-.02))
        r = np.sqrt(2*I2)
        return np.array([lam_u*u, s/lam_u, r*np.cos(phi2), -r*np.sin(phi2)])
    return P
