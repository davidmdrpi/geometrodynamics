from lin import *
a2=1-Rf**2/6
Y=sp.Function('Y')(chi)
P,S,al,be=[sp.Function(n_)(eta) for n_ in ('Phi','Psi','alpha','beta')]
x=[sp.cos(chi),sp.sin(chi)*sp.sin(th)*sp.cos(ph),sp.sin(chi)*sp.sin(th)*sp.sin(ph),sp.sin(chi)*sp.cos(th)]
g=sp.diag(-a2*(1+2*eps*P*Y),a2*(1-2*eps*S*Y),a2*(1-2*eps*S*Y)*sp.sin(chi)**2,a2*(1-2*eps*S*Y)*sp.sin(chi)**2*sp.sin(th)**2)
# grad Y as R^4 vector: Y'(chi) * d x/d chi
dxc=[sp.diff(xi,chi) for xi in x]
phi=[Rf*x[A]+eps*(al*Y*x[A]+be*sp.diff(Y,chi)*dxc[A]) for A in range(4)]
E1,FE,E0,FE0=linearize(g,phi)
import pickle
out={}
for (i,j) in [(0,0),(0,1),(1,1),(2,2),(3,3),(0,2),(1,2)]:
    out[(i,j)]=sp.simplify(bg_reduce(sp.simplify(E1[i,j])))
    print((i,j),out[(i,j)],flush=True)
fe=[sp.simplify(bg_reduce(sp.simplify(e))) for e in FE]
print('FE:',fe,flush=True)
pickle.dump((out,fe),open('sca.pkl','wb'))
