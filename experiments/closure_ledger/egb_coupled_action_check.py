"""Independent full warped-curvature action check, including two Myers faces.

Subtract the static action at each instantaneous configuration and average a
periodic trial over time. Temporal endpoint terms cancel. This checks the
quadratic gravity kinetic coefficient without using the reduced norm formula.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from geometrodynamics.bulk.egb_coupled_radial import CoupledRadial
m=CoupledRadial();a=m.a;R=m.R;alpha=m.shell.alpha
qs,ws=leggauss(80);u=a+(qs+1)*(m.mid-a)/2;ws=ws*(m.mid-a)/2
ell=lambda u:np.cos(2*u);ellp=lambda u:-2*np.sin(2*u)
p=lambda u:.3+.2*np.cos(2*u);pp=lambda u:-.4*np.sin(2*u);ppp=lambda u:-.8*np.cos(2*u)
def bulk(u,Q,Qd,Qdd):
 L=R*np.exp(ell(u)*Q);r=R*np.sin(u)*np.exp(p(u)*Q)
 Ldot=L*ell(u)*Qd;Lddot=L*(ell(u)*Qdd+ell(u)**2*Qd**2);Lp=L*ellp(u)*Q
 rd=r*p(u)*Qd;rdd=r*(p(u)*Qdd+p(u)**2*Qd**2)
 rp=r*(1/np.tan(u)+pp(u)*Q)
 rpp=r*((1/np.tan(u)+pp(u)*Q)**2-1/np.sin(u)**2+ppp(u)*Q)
 rdp=rd*(1/np.tan(u)+pp(u)*Q)+r*pp(u)*Qd
 B=(1+rd**2-rp**2/L**2)/r**2;R2=2*Lddot/L
 huu=rpp/L**2-Ldot/L*rd-Lp*rp/L**3;htu=rdp-Ldot/L*rp
 box=-rdd+huu;diff=-2*rdd*huu+2*htu**2/L**2
 ric=R2+6*B-6*box/r;gb=12*R2*B+24*diff/r**2-24*B*box/r
 return L*r**3*(ric+alpha*gb)/2

def boundary(Q,Qd,Qdd):
 u=a;L=R*np.exp(ell(u)*Q);r=R*np.sin(u)*np.exp(p(u)*Q)
 rd=r*p(u)*Qd;rdd=r*(p(u)*Qdd+p(u)**2*Qd**2)
 rp=r*(1/np.tan(u)+pp(u)*Q);k=rp/(L*r);B=(1+rd**2)/r**2;A=rdd/r
 return -2*r**3*(3*k+2*alpha*(-2*k**3+12*A*k+6*B*k))
expected=np.sum(ws*R**4*np.sin(u)**3*m.C*(-6*ell(u)*p(u)-6*p(u)**2))+24*alpha*m.r**3/(R*np.tan(a))*p(a)**2
def check_action_hessian(eps):
 vals=[]
 for t in np.arange(128)*2*np.pi/128:
  Q=eps*np.cos(t);Qd=-eps*np.sin(t);Qdd=-Q
  vals.append(2*np.sum(ws*(bulk(u,Q,Qd,Qdd)-bulk(u,Q,0.,0.)))+boundary(Q,Qd,Qdd)-boundary(Q,0.,0.))
 measured=float(np.mean(vals)/(eps**2/2))
 return dict(amplitude=eps,expected=float(expected),measured=measured,relative_error=float(abs((measured-expected)/expected)))
