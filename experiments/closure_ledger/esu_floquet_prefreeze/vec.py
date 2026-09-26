from lin import *
a2=1-Rf**2/6
hh=sp.Function('h')(chi); w=sp.Function('w')(eta); sg=sp.Function('s')(eta)
x=[sp.cos(chi),sp.sin(chi)*sp.sin(th)*sp.cos(ph),sp.sin(chi)*sp.sin(th)*sp.sin(ph),sp.sin(chi)*sp.cos(th)]
gam=sp.diag(1,sp.sin(chi)**2,sp.sin(chi)**2*sp.sin(th)**2)
# toroidal vector V = h(chi) d_phi ; lower index V_phi = sin^2chi sin^2th h
g=sp.diag(-a2,a2,a2*sp.sin(chi)**2,a2*sp.sin(chi)**2*sp.sin(th)**2)
Sphi=a2*eps*sg*hh*sp.sin(chi)**2*sp.sin(th)**2
g[0,3]=g[3,0]=Sphi
dx=[sp.diff(xi,ph) for xi in x]
phi=[Rf*x[A]+eps*w*hh*dx[A] for A in range(4)]
E1,FE,E0,FE0=linearize(g,phi)
import pickle
out={}
for (i,j) in [(0,1),(0,2),(0,3),(1,3),(2,3),(1,1),(3,3),(0,0),(1,2)]:
    out[(i,j)]=sp.simplify(bg_reduce(sp.simplify(E1[i,j])))
    print((i,j),out[(i,j)],flush=True)
fe=[sp.simplify(bg_reduce(sp.simplify(e))) for e in FE]
print('FE:',fe)
pickle.dump((out,fe),open('vec.pkl','wb'))
