# Linearization engine: Einstein-frame sigma model, 8piG=1, first order in eps.
import sympy as sp
eta,chi,th,ph,eps=sp.symbols('eta chi theta phi epsilon')
X=[eta,chi,th,ph]
Rf=sp.Function('R')(eta)
q2=sp.Rational(3,4)
def bg_reduce(e):
    # impose R''=-4R, R'^2=4(q^2-R^2)
    e=e.subs(sp.Derivative(Rf,(eta,2)),-4*Rf)
    Rp=sp.Derivative(Rf,eta)
    e=sp.expand(e)
    # replace even powers of R'
    e=e.subs(Rp**4,16*(q2-Rf**2)**2).subs(Rp**3,4*(q2-Rf**2)*Rp).subs(Rp**2,4*(q2-Rf**2))
    return e
def first(e):
    return sp.diff(e,eps).subs(eps,0)
def zeroth(e):
    return e.subs(eps,0)
def linearize(g,phi,extra_simplify=lambda e:e):
    """g: 4x4 Matrix in eps (first order), phi: list of 4 field exprs in eps.
    Returns first-order parts of E_mn=G_mn-T_mn and field eqs E^A."""
    g0=g.applyfunc(zeroth); g1=g.applyfunc(first)
    g0i=sp.simplify(g0.inv()); g1i=-(g0i*g1*g0i)
    gi=g0i+eps*g1i
    n=4
    def d(e,k): return sp.diff(e,X[k])
    G0=[[[sp.simplify(sum(g0i[a,dd]*(d(g0[dd,b],c)+d(g0[dd,c],b)-d(g0[b,c],dd)) for dd in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]
    G1=[[[sp.expand(sum(g0i[a,dd]*(d(g1[dd,b],c)+d(g1[dd,c],b)-d(g1[b,c],dd))+g1i[a,dd]*(d(g0[dd,b],c)+d(g0[dd,c],b)-d(g0[b,c],dd)) for dd in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]
    def Ric(Ga,Gb_list):
        pass
    R0=sp.zeros(n); R1=sp.zeros(n)
    for b in range(n):
        for c in range(b,n):
            r0=sum(d(G0[a][b][c],a)-d(G0[a][b][a],c) for a in range(n))+sum(G0[a][a][dd]*G0[dd][b][c]-G0[a][c][dd]*G0[dd][b][a] for a in range(n) for dd in range(n))
            r1=sum(d(G1[a][b][c],a)-d(G1[a][b][a],c) for a in range(n))+sum(G1[a][a][dd]*G0[dd][b][c]+G0[a][a][dd]*G1[dd][b][c]-G1[a][c][dd]*G0[dd][b][a]-G0[a][c][dd]*G1[dd][b][a] for a in range(n) for dd in range(n))
            R0[b,c]=R0[c,b]=r0; R1[b,c]=R1[c,b]=r1
    S0=sum(g0i[a,b]*R0[a,b] for a in range(n) for b in range(n))
    S1=sum(g1i[a,b]*R0[a,b]+g0i[a,b]*R1[a,b] for a in range(n) for b in range(n))
    Ein1=R1-g1*S0/2-g0*S1/2
    # matter
    ph_=sp.Matrix(phi); r2=(ph_.T*ph_)[0]; f=1-r2/6
    GAB=sp.eye(4)/f+ph_*ph_.T/(6*f*f); U=sp.Rational(3,2)/f**2
    dphi=[[sp.diff(phi[A],X[m]) for A in range(4)] for m in range(n)]
    K=sp.Matrix(n,n,lambda m,v: sum(GAB[A,B]*dphi[m][A]*dphi[v][B] for A in range(4) for B in range(4)))
    trK=sum(gi[a,b]*K[a,b] for a in range(n) for b in range(n))
    T=K-g*(trK/2+U)
    T1=T.applyfunc(first)
    E1=(Ein1-T1)
    E0=(R0-g0*S0/2-T.applyfunc(zeroth))
    # field equations: box phi^A + Gamma^A_BC d phi^B.d phi^C - phi^A/f
    detg=g.det()
    sq0=sp.sqrt(-zeroth(detg)); sq1=-first(detg)/(2*sq0)
    sq=sq0+eps*sq1
    FE=[]
    for A in range(4):
        box=sum(sp.diff(sq*sum(gi[m,v]*dphi[v][A] for v in range(n)),X[m]) for m in range(n))/sq
        conn=2*sum(gi[m,v]*dphi[m][A]*sum(phi[B]*dphi[v][B] for B in range(4)) for m in range(n) for v in range(n))/(6*f)
        FE.append(first(box+conn-phi[A]/f))
    FE0=[]
    for A in range(4):
        box=sum(sp.diff(sq0*sum(g0i[m,v]*sp.diff(zeroth(phi[A]),X[v]) for v in range(n)),X[m]) for m in range(n))/sq0
        p0=[zeroth(p) for p in phi]; f0=zeroth(f)
        conn=2*sum(g0i[m,v]*sp.diff(p0[A],X[m])*sum(p0[B]*sp.diff(p0[B],X[v]) for B in range(4)) for m in range(n) for v in range(n))/(6*f0)
        FE0.append(box+conn-p0[A]/f0)
    return E1,FE,E0,FE0
