"""Constrained clock-section return map of the LRS n=2 tensor sector.

Reduced minisuperspace (kappa = a = 1, conformal time):
    L = -3A'^2 + q'^2/2 + H x'^2/2 - V,   H = A^2 - q^2/6,
    V = -H r/2 + q^2 trinv/2 + 3A^4/2,
with M = diag(e^{2x/s6}, e^{2x/s6}, e^{-4x/s6}), trinv = tr M^-1,
r = 2(2 trinv - tr M^2). The Euler-Lagrange equations coincide with
nonlinear_supported_tt.conformal_rhs restricted to beta = x b0, q = (q,0,0,0);
the energy T + V = 0 is its Hamiltonian constraint. Momenta p_A = -6A',
p_q = q', p_x = H x'. On the section q = 0 (q' < 0) with the constraint solved
for q', the coordinates z = (A, p_A, x, p_x) are canonical for the reduced form
omega = dp_A ^ dA + dp_x ^ dx, and the return map P is symplectic.

Two independent evaluations of P:
- jets (Method 1): this module's scalar equations, fixed-step RK4 with
  Richardson extrapolation, Taylor coefficients by jet transport, and the return
  time solved as a jet (return-time correction included exactly);
- circles (Method 2): nonlinear_supported_tt.conformal_rhs with DOP853 and event
  location, invariant circles by a parameterisation (Newton) method.
"""
import numpy as np
from scipy.integrate import solve_ivp
from .jets import Jet, variables, linear_part, jacobian
from . import nonlinear_supported_tt as d

S6 = np.sqrt(6.)
B0 = np.diag([1., 1., -2.])/S6
ZSTAR = np.array([1., 0., 0., 0.])
OMEGA = np.array([[0., -1., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., -1.], [0., 0., 1., 0.]])
# omega(u, v) = u^T OMEGA v = u_pA v_A - u_A v_pA + u_px v_x - u_x v_px


def _exp(v):
    return v.exp() if isinstance(v, Jet) else np.exp(v)


def _sqrt(v):
    return v.sqrt() if isinstance(v, Jet) else np.sqrt(v)


def shape(x):
    e1, e2 = _exp(2*x/S6), _exp(-4*x/S6)
    trinv = 2/e1+1/e2
    trm2 = 2*e1*e1+e2*e2
    fb_inv = (2/e1-2/e2)/S6          # tr(b0 M^-1)
    fb_sq = (2*e1*e1-2*e2*e2)/S6      # tr(b0 M^2)
    return trinv, trm2, fb_inv, fb_sq


def rhs(y):
    A, Ap, q, qp, x, xp = y
    trinv, trm2, fb_inv, fb_sq = shape(x)
    r = 2*(2*trinv-trm2)
    ell = xp*xp
    Q = q*q
    H = A*A-Q/6
    Hp = 2*A*Ap-q*qp/3
    App = -A*(r+ell)/6+A*A*A
    qpp = -(trinv+(r+ell)/6)*q
    xpp = -(Hp/H)*xp+(-4+Q/H)*fb_inv-4*fb_sq
    return [Ap, App, qp, qpp, xp, xpp]


def constraint(y):
    A, Ap, q, qp, x, xp = y
    trinv, trm2, _, _ = shape(x)
    H = A*A-q*q/6
    r = 2*(2*trinv-trm2)
    return -3*Ap*Ap+qp*qp/2+H*xp*xp/2-H*r/2+q*q*trinv/2+3*A**4/2


def section_to_state(z):
    A, pA, x, px = z
    trinv, trm2, _, _ = shape(x)
    H = A*A
    Ap, xp = -pA/6, px/H
    r = 2*(2*trinv-trm2)
    qp2 = 2*(3*Ap*Ap-H*xp*xp/2+H*r/2-3*A**4/2)
    return [A, Ap, 0*A, -_sqrt(qp2), x, xp]


def state_to_section(y):
    A, Ap, q, qp, x, xp = y
    return [A, -6*Ap, x, (A*A-q*q/6)*xp]


def rk4(y, h, n):
    for _ in range(n):
        k1 = rhs(y)
        k2 = rhs([a+b*(h/2) for a, b in zip(y, k1)])
        k3 = rhs([a+b*(h/2) for a, b in zip(y, k2)])
        k4 = rhs([a+b*h for a, b in zip(y, k3)])
        y = [a+(b+2*c+2*e+f)*(h/6) for a, b, c, e, f in zip(y, k1, k2, k3, k4)]
    return y


# ------------------------------------------------------------------ Method 1
def jet_return_map(order, steps, point=ZSTAR, lead=16, fixed_time=False):
    """Jet of P about `point` with RK4 at h = pi/steps; the last `lead` steps
    are taken with a jet-valued step solving q = 0 by Newton (return time).
    fixed_time=True is the ablation that stops at eta = pi without the
    return-time correction (not a section map)."""
    z = variables(point, order)
    h = np.pi/steps
    if fixed_time:
        y = rk4(section_to_state(z), h, steps)
        return dict(P=state_to_section(y), time=z[0]*0+np.pi, constraint=constraint(y), q=y[2])
    y = rk4(section_to_state(z), h, steps-lead)
    s = z[0]*0+lead*h
    for _ in range(order+3):
        yc = rk4(y, s/lead, lead)
        s = s-yc[2]/yc[3]
    yc = rk4(y, s/lead, lead)
    return dict(P=state_to_section(yc), time=s+(steps-lead)*h, constraint=constraint(yc), q=yc[2])


def richardson(a, b):
    """RK4 extrapolation of jets computed with steps and 2*steps."""
    return [Jet(y.c+(y.c-x.c)/15, y.s['nvar'], y.s['order']) for x, y in zip(a, b)]


def symplectic_defect(P):
    """max |DP^T OMEGA DP - OMEGA| over all jet coefficients (orders 0..order-1)."""
    Jm = jacobian(P)
    worst = 0.
    for i in range(4):
        for k in range(4):
            acc = Jm[0][0]*0-OMEGA[i, k]
            for a in range(4):
                for b in range(4):
                    if OMEGA[a, b]:
                        acc = acc+Jm[a][i]*Jm[b][k]*OMEGA[a, b]
            worst = max(worst, float(np.max(abs(acc.truncate(P[0].s['order']-1).c))))
    return worst


def normal_form(P, point=ZSTAR):
    """Order-3 Birkhoff normal form at a fixed point with one elliptic pair.

    Diagonalise, remove all (nonresonant) quadratic terms by a near-identity
    change, read the resonant w^2 wbar coefficient g of the elliptic component.
    Per iterate, arg w advances by theta + Im(g/mu)|w|^2; Re(g/mu) must vanish
    for a symplectic map. Action I = |omega(q, qbar)| |w|^2 + O(|w|^4).
    """
    order = P[0].s['order']
    if order < 3:
        raise ValueError('order >= 3 required')
    F = [Jet(p.c.copy(), 4, order) for p in P]
    fixed_error = max(abs(F[i].c[0]-point[i]) for i in range(4))
    for i in range(4):
        F[i].c[0] = 0.
    A = linear_part(F)
    lam, vec = np.linalg.eig(A)
    ell = [i for i in range(4) if abs(abs(lam[i])-1) < 1e-6]
    hyp = [i for i in range(4) if i not in ell]
    if len(ell) != 2 or len(hyp) != 2:
        raise ValueError('expected one elliptic and one hyperbolic pair')
    ie = max(ell, key=lambda i: lam[i].imag)
    mu, qv = lam[ie], vec[:, ie]
    hi = sorted(hyp, key=lambda i: -abs(lam[i]))
    V = np.column_stack([vec[:, hi[0]], vec[:, hi[1]], qv, qv.conj()]).astype(complex)
    L = np.array([lam[hi[0]], lam[hi[1]], mu, mu.conjugate()])
    Vi = np.linalg.inv(V)
    xi = variables(np.zeros(4), order, complex)
    delta = [sum((xi[j]*V[i, j] for j in range(4)), xi[0]*0) for i in range(4)]
    Fd = [f.compose(delta) for f in F]
    G = [sum((Fd[j]*Vi[i, j] for j in range(4)), xi[0]*0) for i in range(4)]
    s = xi[0].s
    h = []
    divisors = []
    for i in range(4):
        c = np.zeros(s['n'], complex)
        for m, e in enumerate(s['exps']):
            if s['deg'][m] == 2:
                den = np.prod(L**e)-L[i]      # Lambda h(eta) - h(Lambda eta) = -G2(eta)
                divisors.append(abs(den))
                c[m] = G[i].c[m]/den
        h.append(Jet(c, 4, order))
    Phi = [xi[i]+h[i] for i in range(4)]
    GPhi = [g.compose(Phi) for g in G]
    inv = list(xi)
    for _ in range(order):
        inv = [xi[i]-h[i].compose(inv) for i in range(4)]
    Gn = [g.compose(GPhi) for g in inv]
    g = Gn[2].coef((0, 0, 2, 1))
    omega_qq = float(abs(qv @ OMEGA @ qv.conj()))
    theta = float(np.angle(mu))
    nu = float((g/mu).imag/(2*np.pi*omega_qq))
    return dict(mu=complex(mu), theta=theta, lam=complex(L[0]), g_over_mu=complex(g/mu),
                dissipative=float((g/mu).real), twist=float((g/mu).imag), omega_qq=omega_qq,
                nu=nu, fixed_error=float(fixed_error), min_quadratic_divisor=float(min(divisors)),
                resonance_distances=[float(abs(mu**k-1)) for k in range(1, 5)],
                linear=A.tolist())


# ------------------------------------------------------------------ Method 2
def esu_map(z, tol=(1e-12, 1e-14), method='DOP853'):
    """P(z) by the full homogeneous system; returns z' and the constraint residual."""
    A, pA, x, px = z
    M = np.diag(np.exp(2*x*np.diag(B0)))
    y = d.pack(A, -pA/6, np.zeros(4), np.zeros(4), M, (px/(A*A))*B0)
    E = d.constraints(y)['residual'][0]
    base = E-y[7]**2/2
    if base >= 0:
        raise ValueError('no real clock velocity')
    y[7] = -np.sqrt(-2*base)

    def clock(t, yy):
        return yy[3]
    clock.direction = -1
    sol = solve_ivp(lambda t, yy: d.conformal_rhs(yy), (0., np.pi+.6), y, method=method,
                    rtol=tol[0], atol=tol[1], events=[clock])
    if not sol.success:
        raise ArithmeticError('integration failed: '+sol.message)
    hits = [k for k, t in enumerate(sol.t_events[0]) if t > 1.]
    if not hits:
        raise ArithmeticError('no clock section')
    ye = sol.y_events[0][hits[0]]
    Ae, Ape, q, qp, Me, Le = d.unpack(ye)
    H = Ae*Ae-(q @ q)/6
    xe = float(np.log(np.diag(Me)) @ np.diag(B0)/2)
    xpe = float(np.trace(Le @ B0))
    return np.array([Ae, -6*Ape, xe, H*xpe]), float(abs(d.constraints(ye)['residual'][0])), float(sol.t_events[0][hits[0]])


def fd_jacobian(P, z, h=1e-6):
    J = np.zeros((4, 4))
    for k in range(4):
        e = np.zeros(4)
        e[k] = h
        J[:, k] = (P(z+e)-P(z-e))/(2*h)
    return J


def shift_matrix(M, omega):
    k = np.fft.fftfreq(M, 1/M)
    F = np.fft.fft(np.eye(M), axis=0)
    T = np.fft.ifft(np.exp(1j*k*omega)[:, None]*F, axis=0).real
    dT = np.fft.ifft((1j*k*np.exp(1j*k*omega))[:, None]*F, axis=0).real
    return T, dT


def action(K):
    """(1/2pi) closed integral of p_A dA + p_x dx over the parameterised circle."""
    M = len(K)
    k = np.fft.fftfreq(M, 1/M)
    C = np.fft.fft(K, axis=0)/M
    val = sum(np.sum(np.conj(C[:, p])*(1j*k)*C[:, q]) for p, q in ((1, 0), (3, 2)))
    return float(val.real)


def invariant_circle(P, a, K0, omega0, M, tol=1e-12, maxit=15, h=1e-6):
    """Solve P(K(theta)) = K(theta + omega) with x-harmonic c1 = a/2 (real).

    P returns z' only. K has shape (M, 4) sampled at theta_j = 2 pi j/M.
    """
    K, omega = K0.copy(), omega0
    th = 2*np.pi*np.arange(M)/M
    e1 = np.exp(-1j*th)/M
    hist = []
    for it in range(maxit):
        PK = np.array([P(k) for k in K])
        T, dT = shift_matrix(M, omega)
        R = PK-T @ K
        c1 = e1 @ K[:, 2]
        res = np.r_[R.ravel(), c1.real-a/2, c1.imag]
        hist.append(float(np.max(abs(res))))
        if hist[-1] < tol:
            break
        J = np.zeros((4*M+2, 4*M+1))
        for j in range(M):
            J[4*j:4*j+4, 4*j:4*j+4] = fd_jacobian(P, K[j], h)
        J[:4*M, :4*M] -= np.kron(T, np.eye(4))
        J[:4*M, 4*M] = -(dT @ K).ravel()
        J[4*M, 2:4*M:4] = e1.real
        J[4*M+1, 2:4*M:4] = e1.imag
        step = np.linalg.lstsq(J, -res, rcond=None)[0]
        K = K+step[:4*M].reshape(M, 4)
        omega = omega+step[4*M]
    PK = np.array([P(k) for k in K])
    T, _ = shift_matrix(M, omega)
    Ck = np.abs(np.fft.fft(K-K.mean(0), axis=0))/M
    kk = np.abs(np.fft.fftfreq(M, 1/M))
    tail = float(Ck[kk >= M//2-2].max())
    return dict(K=K, omega=float(omega), residual=float(np.max(abs(PK-T @ K))), iterations=it+1,
                history=hist, action=action(K), fourier_tail=tail)


def linear_circle(A, a, M, point=ZSTAR):
    lam, vec = np.linalg.eig(A)
    ie = max((i for i in range(4) if abs(abs(lam[i])-1) < 1e-6), key=lambda i: lam[i].imag)
    q = vec[:, ie]/vec[2, ie]
    th = 2*np.pi*np.arange(M)/M
    K = point+a*np.real(q[None, :]*np.exp(1j*th)[:, None])
    return K, float(np.angle(lam[ie]))
