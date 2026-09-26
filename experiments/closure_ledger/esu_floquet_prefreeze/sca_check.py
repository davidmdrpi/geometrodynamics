import sympy as sp
e,k=sp.symbols('eta k')
q=sp.sqrt(3)/2; R=q*sp.cos(2*e); Rp=sp.diff(R,e); f=1-R**2/6; fp=sp.diff(f,e); H=fp/(2*f)
a0,a1,b0,b1=sp.symbols('a0 a1 b0 b1')   # alpha, alpha', beta, beta'
Psi=sp.Symbol('Psi')
Phi=Psi-2*R*b0/f
Psi1=-H*Phi+((R*b1-Rp*b0)/f+Rp*a0/f**2)/2
C0=-Rp*a1/f**2+(R*Rp/f)*Psi1+(k*R/f)*b0-R*(12*f-7)*a0/f**3-2*(f*k-12*f+9)*Psi/f+3*(8*f-7)*Phi/f
PsiSol=sp.solve(C0,Psi)[0]
PhiS=Phi.subs(Psi,PsiSol); Psi1S=Psi1.subs(Psi,PsiSol)
Phi1=Psi1S-2*(Rp*b0+R*b1)/f+2*R*b0*fp/f**2
a2=-(R*Rp/(3*f))*a1+3*Rp*Psi1S+Rp*Phi1+2*k*b0-(k-6+(25*f-14)/f**2)*a0-6*R*PsiSol-8*R*PhiS
b2=-(k+2*R**2/f)*b0+2*a0/f
def dt(expr):
    return sp.diff(expr,e)+sp.diff(expr,a0)*a1+sp.diff(expr,a1)*a2+sp.diff(expr,b0)*b1+sp.diff(expr,b1)*b2
import random
# numeric spot checks (exact rational-ish evaluation) of: d/dt(PsiSol) == Psi1S ; trace equation E11
Psi2=dt(Psi1S)
E11=2*Psi2-Rp*a1/f**2-(2*R*Rp/(3*f))*Psi1S-(R*Rp/(3*f))*Phi1+(k*R/f)*b0-R*(6*f-7)*a0/f**3-2*(4*f-3)*PsiSol/f-(8*f-9)*PhiS/f
chk1=dt(PsiSol)-Psi1S
for trial in range(4):
    vals={e:sp.Rational(random.randint(1,300),97),k:random.choice([8,15,24,35,48]),a0:sp.Rational(random.randint(-9,9),7),a1:sp.Rational(random.randint(-9,9),5),b0:sp.Rational(random.randint(-9,9),3),b1:sp.Rational(random.randint(-9,9),11)}
    print('consistency dPsi/deta - C1:',sp.N(chk1.subs(vals),30),'  trace eq E11:',sp.N(E11.subs(vals),30))
