# PRE-FREEZE: linear facts only (epsilon = 0 background and its linearisation)
import numpy as np
from scipy.integrate import solve_ivp
from geometrodynamics.waves import nonlinear_supported_tt as d
y0=d.pack(1.,0.,np.r_[np.sqrt(3)/2,0,0,0],np.zeros(4),np.eye(3),np.zeros((3,3)))
print('constraint residual',d.constraints(y0)['residual'])
f=lambda t,y:d.conformal_rhs(y)
s=solve_ivp(f,(0,np.pi),y0,method='DOP853',rtol=1e-13,atol=1e-15)
print('return error',np.abs(s.y[:,-1]-y0)[[0,1,3,4,7]].max())
# homogeneous (A,A',q0,q0') linear block via finite differences
idx=[0,1,3,7]; h=1e-7; Mx=np.zeros((4,4))
for j,i in enumerate(idx):
    yp=y0.copy(); yp[i]+=h; ym=y0.copy(); ym[i]-=h
    a=solve_ivp(f,(0,np.pi),yp,method='DOP853',rtol=1e-13,atol=1e-15).y[:,-1]
    b=solve_ivp(f,(0,np.pi),ym,method='DOP853',rtol=1e-13,atol=1e-15).y[:,-1]
    Mx[:,j]=(a-b)[idx]/(2*h)
print('homogeneous-block multipliers per pi:',np.linalg.eigvals(Mx))
