"""Independent coordinate TT and first-order electric Weyl checks.

Uses the coordinate Christoffel definition, not the radial eigen-equation.
Weyl is formed from R_0i0j=-h''H_ij/2, Ricci and the Ricci scalar
of -deta^2+gamma+h(eta)H. Its Einstein orthonormal value divides by f.
"""
import functools
import numpy as np
import sympy as sp


def check_degree(n):
    chi, th, ph = sp.symbols('chi theta phi', real=True)
    coords = (chi, th, ph)
    s, c = sp.sin(chi), sp.cos(chi)
    g = sp.diag(1, s*s, s*s*sp.sin(th)**2)
    gi = g.inv()
    # Overall Y normalization cancels in these differential identities.
    Y = (3*sp.cos(th)**2-1)/2
    A = sp.gegenbauer(n-2, 3, c)
    B = s*s*(sp.diff(A, chi)+3*c*A/s)/6
    C = -s*s*A/2
    D = s*s*(sp.diff(B, chi)+2*c*B/s-A/2)/2
    H = sp.zeros(3)
    H[0, 0] = A*Y
    H[0, 1] = H[1, 0] = B*sp.diff(Y, th)
    H[1, 1] = C*Y+D*(sp.diff(Y, th, 2)+3*Y)
    H[2, 2] = C*sp.sin(th)**2*Y+D*(sp.sin(th)*sp.cos(th)*sp.diff(Y, th)+3*sp.sin(th)**2*Y)
    H = H.applyfunc(sp.expand_trig)
    G = [[[sp.simplify(sum(gi[a,d]*(sp.diff(g[d,b],coords[k])+sp.diff(g[d,k],coords[b])-sp.diff(g[b,k],coords[d])) for d in range(3))/2)
           for k in range(3)] for b in range(3)] for a in range(3)]
    @functools.lru_cache(None)
    def dh(k, i, j):
        return sp.diff(H[i,j],coords[k])-sum(G[a][k][i]*H[a,j]+G[a][k][j]*H[i,a] for a in range(3))
    trace = sum(gi[i,i]*H[i,i] for i in range(3))
    div = [sum(gi[i,i]*dh(i,i,j) for i in range(3)) for j in range(3)]
    lap = sp.zeros(3)
    for i in range(3):
        for j in range(i,3):
            lap[i,j] = lap[j,i] = sum(gi[k,k]*(sp.diff(dh(k,i,j),coords[k])-sum(
                G[a][k][k]*dh(a,i,j)+G[a][k][i]*dh(k,a,j)+G[a][k][j]*dh(k,i,a) for a in range(3))) for k in range(3))
    # Direct coordinate variation of Gamma and Ricci, independent of lap above.
    dgi = -gi*H*gi
    dG = [[[sp.expand(sum(gi[a,d]*(sp.diff(H[d,b],coords[k])+sp.diff(H[d,k],coords[b])-sp.diff(H[b,k],coords[d]))
                 + dgi[a,d]*(sp.diff(g[d,b],coords[k])+sp.diff(g[d,k],coords[b])-sp.diff(g[b,k],coords[d])) for d in range(3))/2)
            for k in range(3)] for b in range(3)] for a in range(3)]
    ric = sp.zeros(3)
    for i in range(3):
        for j in range(i,3):
            ric[i,j] = ric[j,i] = sum(sp.diff(dG[a][i][j],coords[a])-sp.diff(dG[a][i][a],coords[j]) for a in range(3))+sum(
                dG[a][a][b]*G[b][i][j]+G[a][a][b]*dG[b][i][j]-dG[a][j][b]*G[b][i][a]-G[a][j][b]*dG[b][i][a] for a in range(3) for b in range(3))
    delta_scalar = sum(gi[i,i]*ric[i,i]+2*dgi[i,i]*g[i,i] for i in range(3))
    # Choose independent h=1, h''=7/5; their coefficients are linear.
    hpp = sp.Rational(7,5)
    ric4 = ric+hpp*H/2
    ric00 = -hpp*trace/2
    scalar4 = delta_scalar+hpp*trace
    E = -hpp*H/2+ric4/2-g*ric00/2-H-g*scalar4/6
    k = n*(n+2)
    expected = (k-hpp)*H/4
    expr = [trace]+div+list(lap+(k-2)*H)+list(E-expected)
    fn = sp.lambdify((chi,th,ph), expr, 'numpy', cse=True)
    base = sp.lambdify((chi,th,ph), H, 'numpy', cse=True)
    residuals = []
    for point in ((.61,.73,.2),(1.21,1.07,.5),(2.31,2.01,1.)):
        scale = max(1., k*np.max(np.abs(base(*point))))
        residuals.append(np.max(np.abs(np.asarray(fn(*point),float)))/scale)
    return float(max(residuals))


def run():
    out = {}
    for n in (2,3,5):
        out[str(n)] = check_degree(n)
        print('coordinate TT/Weyl', n, out[str(n)], flush=True)
    return out
