# PRE-FREEZE: equation agreement and linear/quadratic facts only. No order-3 ESU jet, no ESU circle.
import numpy as np, time
from geometrodynamics.waves import r3_return_map as rm, nonlinear_supported_tt as d, esu_floquet as fl
rng=np.random.default_rng(3); worst=0
for _ in range(20):
    A,Ap,q,qp,x,xp=1+.1*rng.normal(),.1*rng.normal(),.3*rng.normal(),rng.normal(),.1*rng.normal(),.3*rng.normal()
    y=d.pack(A,Ap,np.r_[q,0,0,0],np.r_[qp,0,0,0],np.diag(np.exp(2*x*np.diag(rm.B0))),xp*rm.B0)
    f=d.conformal_rhs(y); _,_,_,_,Mp,Lp=d.unpack(f)
    red=rm.rhs([A,Ap,q,qp,x,xp]); full=[f[0],f[1],f[3],f[7],np.trace(np.linalg.solve(np.diag(np.exp(2*x*np.diag(rm.B0))),Mp))/2@np.eye(1)[0] if False else xp,np.trace(Lp@rm.B0)]
    worst=max(worst,max(abs(a-b) for a,b in zip(red,full)))
    worst=max(worst,abs(rm.constraint([A,Ap,q,qp,x,xp])-d.constraints(y)['residual'][0]))
print('reduced vs full rhs/constraint max diff: %.1e'%worst)
tr310=np.trace(fl.monodromy('T',2))
for order,steps in ((1,1024),(1,2048)):
    t=time.time(); out=rm.jet_return_map(order,steps)
    A=rm.linear_part(out['P']); ev=np.linalg.eigvals(A)
    print(order,steps,'fixed pt err %.1e'%max(abs(p.c[0]-z) for p,z in zip(out['P'],rm.ZSTAR)),'eig',np.round(ev,10),
          'tensor trace-310 %.1e'%(A[2,2]+A[3,3]-tr310),'A-x coupling %.1e'%max(abs(A[0:2,2:4]).max(),abs(A[2:4,0:2]).max()),
          'symp %.1e'%np.abs(A.T@rm.OMEGA@A-rm.OMEGA).max(),'%.1fs'%(time.time()-t))
t=time.time(); out=rm.jet_return_map(2,1024); print('order-2 time %.1fs, symplectic defect (orders<=1) %.1e, return-time jet const %.12f'%(time.time()-t,rm.symplectic_defect(out['P']),out['time'].c[0]))
