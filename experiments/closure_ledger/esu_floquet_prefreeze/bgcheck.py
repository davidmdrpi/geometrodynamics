from lin import *
a2=1-Rf**2/6
x=[sp.cos(chi),sp.sin(chi)*sp.sin(th)*sp.cos(ph),sp.sin(chi)*sp.sin(th)*sp.sin(ph),sp.sin(chi)*sp.cos(th)]
g=sp.diag(-a2,a2,a2*sp.sin(chi)**2,a2*sp.sin(chi)**2*sp.sin(th)**2)
phi=[Rf*xi for xi in x]
E1,FE,E0,FE0=linearize(g,phi)
print('background Einstein residuals:',[sp.simplify(bg_reduce(sp.simplify(E0[i,j]))) for i in range(4) for j in range(i,4)])
print('background field residuals:',[sp.simplify(bg_reduce(sp.simplify(e))) for e in FE0])
