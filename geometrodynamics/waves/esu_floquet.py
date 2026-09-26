"""Floquet maps of linear perturbations about the breathing four-scalar ESU.

Einstein frame, kappa=1, a=1, Lambda=3/2, conformal time eta:
g_E = f(-deta^2+gamma), phi = R x, R = q cos(2 eta), f = 1-R^2/6.
Reduced equations are those frozen in docs/esu_floquet_refocusing_prereg.md.
"""
import numpy as np
from scipy.integrate import solve_ivp

Q = np.sqrt(3)/2
SECTORS = ('T', 'V', 'S')
DIM = dict(T=2, V=2, S=4, C=2)


def background(eta):
    R = Q*np.cos(2*eta)
    Rp = -2*Q*np.sin(2*eta)
    f = 1-R*R/6
    fp = -R*Rp/3
    return R, Rp, f, fp, fp/(2*f)


def rhs_tensor(eta, y, n):
    """y=(h, p=f h'); h''+(f'/f)h'+[n(n+2)+2R^2/f]h=0."""
    R, _, f, _, _ = background(eta)
    h, p = y
    return np.array([p/f, -(n*(n+2)*f+2*R*R)*h])


def vector_s(eta, w, wp, n):
    R, Rp, f, _, _ = background(eta)
    lam = (n+1)**2-4
    return 2*(R*wp-Rp*w)/(2*R*R+lam*f)


def rhs_vector(eta, y, n, coupled=True):
    """y=(w, w'); w''+(n+1)^2 w = R s'+2R's, s from the 0i constraint."""
    R, Rp, f, fp, _ = background(eta)
    w, wp = y
    if not coupled:
        return np.array([wp, -(n+1)**2*w])
    s = vector_s(eta, w, wp, n)
    sp = (-2*R*w-fp*s)/f
    return np.array([wp, -(n+1)**2*w+R*sp+2*Rp*s])


def scalar_metric(eta, y, n):
    """Solve (00),(0i),(traceless) for Psi, Phi, Psi', Phi' given (a,a',b,b')."""
    R, Rp, f, fp, H = background(eta)
    k = n*(n+2)
    a, ap, b, bp = y
    phi0 = -2*R*b/f                     # Phi = Psi + phi0
    J = ((R*bp-Rp*b)/f+Rp*a/f**2)/2     # Psi' = -H Phi + J
    # (00): c_a1 a' + c_P1 Psi' + c_b b + c_a a + c_Psi Psi + c_Phi Phi = 0
    c_a1, c_P1, c_b = -Rp/f**2, R*Rp/f, k*R/f
    c_a, c_Psi, c_Phi = -R*(12*f-7)/f**3, -2*(f*k-12*f+9)/f, 3*(8*f-7)/f
    const = c_a1*ap+c_P1*(J-H*phi0)+c_b*b+c_a*a+c_Phi*phi0
    coef = c_Psi+c_Phi-c_P1*H
    Psi = -const/coef
    Phi = Psi+phi0
    Psi1 = -H*Phi+J
    Phi1 = Psi1-2*(Rp*b+R*bp)/f+2*R*b*fp/f**2
    return Psi, Phi, Psi1, Phi1


def rhs_scalar(eta, y, n):
    R, Rp, f, fp, H = background(eta)
    k = n*(n+2)
    a, ap, b, bp = y
    Psi, Phi, Psi1, Phi1 = scalar_metric(eta, y, n)
    app = (-(R*Rp/(3*f))*ap+3*Rp*Psi1+Rp*Phi1+2*k*b
           -(k-6+(25*f-14)/f**2)*a-6*R*Psi-8*R*Phi)
    bpp = -(k+2*R*R/f)*b+2*a/f
    return np.array([ap, app, bp, bpp])


def rhs_free(eta, y, n):
    """Free conformal test field omega=n+1 (control C1)."""
    return np.array([y[1], -(n+1)**2*y[0]])


def rhs_bare_tensor(eta, y, n):
    """Static-ESU tensor without support mass (control C4): omega^2=n(n+2)."""
    return np.array([y[1], -n*(n+2)*y[0]])


RHS = dict(T=rhs_tensor, V=rhs_vector, S=rhs_scalar, C=rhs_free, B=rhs_bare_tensor)


def _matrix_rhs(sector, n, **kw):
    d = 4 if sector == 'S' else 2
    f = RHS[sector]
    def rhs(eta, Y):
        M = Y.reshape(d, d)
        return np.column_stack([f(eta, M[:, j], n, **kw) for j in range(d)]).ravel()
    return rhs, d


def monodromy(sector, n, t0=0., t1=np.pi, rtol=1e-12, atol=1e-14, **kw):
    rhs, d = _matrix_rhs(sector, n, **kw)
    sol = solve_ivp(rhs, (t0, t1), np.eye(d).ravel(), method='DOP853', rtol=rtol, atol=atol)
    if not sol.success:
        raise ArithmeticError(f'integration failed {sector} n={n}')
    return sol.y[:, -1].reshape(d, d)


def monodromy_rk4(sector, n, steps, t0=0., t1=np.pi, **kw):
    """Independent fixed-step classical RK4 for the same fundamental matrix."""
    rhs, d = _matrix_rhs(sector, n, **kw)
    h = (t1-t0)/steps
    Y = np.eye(d).ravel()
    t = t0
    for _ in range(steps):
        k1 = rhs(t, Y); k2 = rhs(t+h/2, Y+h*k1/2)
        k3 = rhs(t+h/2, Y+h*k2/2); k4 = rhs(t+h, Y+h*k3)
        Y = Y+h*(k1+2*k2+2*k3+k4)/6
        t += h
    return Y.reshape(d, d)


def parity(sector, n):
    """Antipodal parity of the physical carrier (field for S,V; metric for T)."""
    return (-1)**(n+1) if sector == 'S' else (-1)**n


def energy_scaling(sector, n):
    """(u, u'/(n+1)) coordinates at eta=0; T uses p=f h' so u'=p/f(0)."""
    w = n+1
    f0 = background(0.)[2]
    if sector == 'T':
        return np.diag([1., 1/(w*f0)])
    if sector == 'S':
        return np.diag([1., 1/w, 1., 1/w])
    return np.diag([1., 1/w])


def observables(sector, n, M):
    d = M.shape[0]
    Rf = parity(sector, n)*M
    S = energy_scaling(sector, n)
    Rn = S@Rf@np.linalg.inv(S)
    mu = np.linalg.eigvals(M)
    out = dict(n=n, trace=float(np.trace(M)), det=float(np.linalg.det(M)),
               max_abs_multiplier=float(np.max(abs(mu))),
               multipliers=[[float(z.real), float(z.imag)] for z in mu],
               fidelity=float(-np.trace(Rf)/d),
               defect=float(np.linalg.norm(Rn+np.eye(d), 2)))
    if d == 2:
        c = -np.trace(Rf)/2
        out['theta'] = float(np.arccos(c)) if abs(c) <= 1 else None
    else:
        out['harmonic_fidelity'] = float(-np.trace((-1)**n*M)/d)
        ph = np.angle(np.linalg.eigvals(Rf))
        out['refocusing_eigenphases'] = sorted(float(v) for v in ph)
    return out


def stability(mu_max):
    if mu_max <= 1+1e-7:
        return 'ELLIPTIC'
    if mu_max >= 1+1e-5:
        return 'HYPERBOLIC'
    return 'MARGINAL'


def wkb_masses(samples=200001):
    eta = np.linspace(0, np.pi, samples)
    R, _, f, _, H = background(eta)
    mean = lambda v: float(np.trapezoid(v, eta)/np.pi) if hasattr(np, 'trapezoid') else float(np.trapz(v, eta)/np.pi)
    m2 = mean(2*R*R/f)
    return dict(mean_2R2_over_f=m2, mean_H2=mean(H*H),
                m_T2=m2-mean(H*H)-1, m_V2=m2,
                analytic_mean_2R2_over_f=12*(np.sqrt(8/7)-1))


def scalar_trace_residual(eta, y, n):
    """Unused ij-trace equation along a solution, relative to its term scale."""
    R, Rp, f, fp, H = background(eta)
    Rpp = -4*R
    k = n*(n+2)
    a, ap, b, bp = y
    Psi, Phi, Psi1, Phi1 = scalar_metric(eta, y, n)
    _, app, _, bpp = rhs_scalar(eta, y, n)
    Hp = (-(Rp*Rp+R*Rpp)/3)/(2*f)-fp*fp/(2*f*f)
    J1 = ((Rp*bp+R*bpp-Rpp*b-Rp*bp)/f-(R*bp-Rp*b)*fp/f**2
          +(Rpp*a+Rp*ap)/f**2-2*Rp*a*fp/f**3)/2
    Psi2 = -Hp*Phi-H*Phi1+J1
    terms = [2*Psi2, -Rp*ap/f**2, -(2*R*Rp/(3*f))*Psi1, -(R*Rp/(3*f))*Phi1,
             (k*R/f)*b, -R*(6*f-7)*a/f**3, -2*(4*f-3)*Psi/f, -(8*f-9)*Phi/f]
    return abs(sum(terms))/max(sum(abs(t) for t in terms), 1e-300)


def vector_constraint_residual(eta, y, n):
    """d/deta of the 0i-constraint s versus the unused ij equation (fs)'=-2Rw."""
    R, Rp, f, fp, _ = background(eta)
    Rpp = -4*R
    lam = (n+1)**2-4
    w, wp = y
    wpp = rhs_vector(eta, y, n)[1]
    num = 2*(R*wp-Rp*w); den = 2*R*R+lam*f
    numd = 2*(Rp*wp+R*wpp-Rpp*w-Rp*wp); dend = 4*R*Rp+lam*fp
    s = num/den
    sd = numd/den-num*dend/den**2
    lhs = fp*s+f*sd
    return abs(lhs+2*R*w)/max(abs(lhs)+abs(2*R*w), 1e-300)
