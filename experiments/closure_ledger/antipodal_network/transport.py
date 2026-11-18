"""Held-fixed exterior preparation under uniform full-mesh refinement."""
import json
import numpy as np
from scipy.sparse.linalg import splu
from .refinement import domain

def run(level,dt):
    ((tm,tg,tk),other),n,nt,seed,traces=domain(level,True)
    mass,gyro,stiff=[a+b for a,b in zip((tm,tg,tk),other)]
    u=np.zeros(n)
    v=seed/np.sqrt(seed @ (mass @ seed)/2)
    assert v @ (tm @ v)==0
    step=splu((mass+dt*gyro/2+dt**2*stiff/4).tocsc())
    work=np.zeros(2)
    balance=drift=0.
    history=[]
    for i in range(1,round(3/dt)+1):
        dv=step.solve(-dt*(gyro @ v+stiff @ (u+dt*v/2)))
        midv=v+dv/2
        midu=u+dt*midv/2
        reaction=tm @ (dv/dt)+tg @ midv+tk @ midu
        work-=dt*((midv*reaction) @ traces)
        u+=dt*midv
        v+=dv
        et=(u @ (tk @ u)+v @ (tm @ v))/2
        total=(u @ (stiff @ u)+v @ (mass @ v))/2
        balance=max(balance,abs(et+sum(work)))
        drift=max(drift,abs(total-1))
        if i % round(.1/dt)==0:
            history.append((i*dt,*work,et))
    assert max(balance,drift)<1e-10
    return dict(level=level,step=dt,cover_tetrahedra=nt,quotient_nodes=n,
                max_tube_balance_error=balance,max_total_energy_error=drift,
                columns=['t','net_out_A','net_out_B','tube_energy'], history=history)
if __name__=='__main__':
    rows=[run(j,dt) for j,dt in [(0,.01),(1,.01),(2,.01),(2,.005)]]
    print(json.dumps(rows,indent=2))
