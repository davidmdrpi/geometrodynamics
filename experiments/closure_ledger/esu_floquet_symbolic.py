"""Gate G1: reduced equations imposed -> full linearized Einstein-sigma equations vanish.

Explicit harmonics: zonal scalar Y=sin((n+1)chi)/sin(chi) (n=2,3,4), toroidal
vector h=C^(2)_(n-1)(cos chi) d_phi (n=2,3), left-invariant TT sigma1^2-sigma2^2 (n=2).
Every function of eta gets random values; the reduced system from
geometrodynamics.waves.esu_floquet supplies the dependent ones.
"""
import json, random, sys
import numpy as np
import sympy as sp
from .esu_linearization import linearize, eta, chi, th, ph, eps, Rf
from geometrodynamics.waves import esu_floquet as fl

x = [sp.cos(chi), sp.sin(chi)*sp.sin(th)*sp.cos(ph), sp.sin(chi)*sp.sin(th)*sp.sin(ph), sp.sin(chi)*sp.cos(th)]
a2 = 1-Rf**2/6
g_bg = sp.diag(-a2, a2, a2*sp.sin(chi)**2, a2*sp.sin(chi)**2*sp.sin(th)**2)


def evaluate(exprs, funcs, eta0, rng):
    """Substitute function jets at eta0 and random angles; return floats."""
    R0 = fl.Q*np.cos(2*eta0); R1 = -2*fl.Q*np.sin(2*eta0)
    vals = {}
    for fn, jet in funcs:
        for order in (2, 1):
            vals[sp.Derivative(fn, (eta, order)) if order > 1 else sp.Derivative(fn, eta)] = jet[order]
        vals[fn] = jet[0]
    vals[sp.Derivative(Rf, (eta, 2))] = -4*R0; vals[sp.Derivative(Rf, eta)] = R1; vals[Rf] = R0
    ang = {chi: rng.uniform(.3, 2.8), th: rng.uniform(.3, 2.8), ph: rng.uniform(0, 6)}
    out = []
    for e in exprs:
        v = e
        for key in sorted(vals, key=lambda z: -sp.count_ops(z)):
            v = v.subs(key, vals[key])
        out.append(complex(sp.N(v.subs(ang).subs(eta, eta0), 30)))
    return out


def scale(exprs_abs):
    return max(1., max(abs(v) for v in exprs_abs))


def scalar_case(n, rng):
    Y = sp.sin((n+1)*chi)/sp.sin(chi)
    P, S, al, be = [sp.Function(s)(eta) for s in ('Phi', 'Psi', 'alpha', 'beta')]
    g = g_bg.copy()
    g[0, 0] = -a2*(1+2*eps*P*Y)
    for i in (1, 2, 3):
        g[i, i] = g_bg[i, i]*(1-2*eps*S*Y)
    dxc = [sp.diff(v, chi) for v in x]
    phi = [Rf*x[A]+eps*(al*Y*x[A]+be*sp.diff(Y, chi)*dxc[A]) for A in range(4)]
    E1, FE, _, _ = linearize(g, phi)
    exprs = [E1[i, j] for i in range(4) for j in range(i, 4)]+list(FE)
    worst = 0.
    for _ in range(3):
        eta0 = rng.uniform(0, np.pi)
        y = np.array([rng.uniform(-1, 1) for _ in range(4)])
        Psi, Phi, Psi1, Phi1 = fl.scalar_metric(eta0, y, n)
        dy = fl.rhs_scalar(eta0, y, n)
        # Psi'' from the exact differentiated (0i) relation, as in the frozen trace check
        R, Rp, f, fp, H = fl.background(eta0); Rpp = -4*R
        a, ap, b, bp = y; app, bpp = dy[1], dy[3]
        Hp = (-(Rp*Rp+R*Rpp)/3)/(2*f)-fp*fp/(2*f*f)
        J1 = ((Rp*bp+R*bpp-Rpp*b-Rp*bp)/f-(R*bp-Rp*b)*fp/f**2+(Rpp*a+Rp*ap)/f**2-2*Rp*a*fp/f**3)/2
        Psi2 = -Hp*Phi-H*Phi1+J1
        funcs = [(al, (a, ap, app)), (be, (b, bp, bpp)), (P, (Phi, Phi1, 0.)), (S, (Psi, Psi1, Psi2))]
        vals = evaluate(exprs, funcs, eta0, rng)
        ref = scale([a, ap, b, bp, Psi, Phi, Psi1, Phi1])
        worst = max(worst, max(abs(v) for v in vals)/ref)
    return worst


def vector_case(n, rng):
    h = sp.gegenbauer(n-1, 2, sp.cos(chi))
    w, s = sp.Function('w')(eta), sp.Function('s')(eta)
    g = g_bg.copy()
    g[0, 3] = g[3, 0] = a2*eps*s*h*sp.sin(chi)**2*sp.sin(th)**2
    dx = [sp.diff(v, ph) for v in x]
    phi = [Rf*x[A]+eps*w*h*dx[A] for A in range(4)]
    E1, FE, _, _ = linearize(g, phi)
    exprs = [E1[i, j] for i in range(4) for j in range(i, 4)]+list(FE)
    worst = 0.
    for _ in range(3):
        eta0 = rng.uniform(0, np.pi)
        y = np.array([rng.uniform(-1, 1) for _ in range(2)])
        R, Rp, f, fp, _ = fl.background(eta0)
        sv = fl.vector_s(eta0, y[0], y[1], n)
        sp1 = (-2*R*y[0]-fp*sv)/f
        wpp = fl.rhs_vector(eta0, y, n)[1]
        vals = evaluate(exprs, [(w, (y[0], y[1], wpp)), (s, (sv, sp1, 0.))], eta0, rng)
        worst = max(worst, max(abs(v) for v in vals)/scale([y[0], y[1], sv, sp1]))
    return worst


def tensor_case(rng):
    """n=2 left-invariant TT tensor sigma1^2-sigma2^2 in polar coordinates."""
    X = [eta, chi, th, ph]
    def form(a, b, c, d):
        # 1-form x_a dx_b - x_b dx_a + x_c dx_d - x_d dx_c, pulled back to (chi,theta,phi)
        return [x[a]*sp.diff(x[b], X[m])-x[b]*sp.diff(x[a], X[m])+x[c]*sp.diff(x[d], X[m])-x[d]*sp.diff(x[c], X[m]) for m in (1, 2, 3)]
    s1, s2 = form(0, 1, 2, 3), form(0, 2, 3, 1)
    hmat = sp.Matrix(3, 3, lambda i, j: s1[i]*s1[j]-s2[i]*s2[j])
    b = sp.Function('b')(eta)
    g = g_bg.copy()
    for i in range(3):
        for j in range(3):
            g[i+1, j+1] = g_bg[i+1, j+1]+eps*a2*b*hmat[i, j]
    phi = [Rf*v for v in x]
    E1, FE, _, _ = linearize(g, phi)
    exprs = [E1[i, j] for i in range(4) for j in range(i, 4)]+list(FE)
    worst = 0.
    for _ in range(3):
        eta0 = rng.uniform(0, np.pi)
        h0, p0 = rng.uniform(-1, 1), rng.uniform(-1, 1)
        R, _, f, fp, _ = fl.background(eta0)
        hp = p0/f
        hpp = -(fp/f)*hp-(8+2*R*R/f)*h0
        vals = evaluate(exprs, [(b, (h0, hp, hpp))], eta0, rng)
        worst = max(worst, max(abs(v) for v in vals)/scale([h0, hp, hpp]))
    return worst


def run(seed=2026092610):
    rng = random.Random(seed)
    out = dict(seed=seed, scalar={}, vector={}, tensor={})
    out['tensor']['2'] = tensor_case(rng); print('tensor n=2', out['tensor']['2'], flush=True)
    for n in (2, 3):
        out['vector'][str(n)] = vector_case(n, rng); print('vector', n, out['vector'][str(n)], flush=True)
    for n in (2, 3, 4):
        out['scalar'][str(n)] = scalar_case(n, rng); print('scalar', n, out['scalar'][str(n)], flush=True)
    return out


if __name__ == '__main__':
    result = run()
    json.dump(result, open(sys.argv[1], 'w'), indent=2)
