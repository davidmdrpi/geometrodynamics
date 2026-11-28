"""Double-null Einstein--Maxwell--scalar evolution relative to an exact cold bridge."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import beta, betainc

from .inner_horizon_flux import DELTA, background, events, reference

OUT = Path(__file__).resolve().parents[1]/'runs/20261127_inner_horizon_flux/coupled.json'


def coupled(step, epsilon, delta=DELTA, end_fraction=.9999):
    ref=reference(.01)
    a,b,c=ref['squared_horizons']
    q,mu,kb=ref['charge'],ref['mu'],ref['kappa'][1]
    targets,times,_=events(ref,delta)
    last=int(np.searchsorted(times/ref['turn_tau'],end_fraction-1e-10))
    extent=int(np.ceil(targets[last]/step))+4
    axis=np.arange(-extent,extent+1)*step
    K=delta*np.sinh(axis)
    Kp=delta*np.cosh(axis)
    n0=-kb*np.tanh(axis)
    grid=8e-5
    m=int(np.ceil(np.arcsinh(K[-1]**2)/grid))
    rt,lt=background(ref,grid*np.arange(-m,m+1))
    def geometry(u,v):
        index=np.arcsinh(-K[u]*K[v])/grid+m
        j=np.clip(index.astype(int),0,len(rt)-2)
        w=index-j
        rv=rt[j]+w*(rt[j+1]-rt[j])
        loge=lt[j]+w*(lt[j+1]-lt[j])+np.log(Kp[u]*Kp[v]/kb**2)
        return rv,loge,n0[u],n0[v]

    u=np.arange(len(axis)-1,-1,-1)
    v=np.arange(len(axis))
    rv,loge,np0,nm0=geometry(u,v)
    s=rv**2
    load=epsilon*ref['mass_gap']
    aa,bb,cc=np.sort(np.roots([1.,-2.,2*(mu+load),-2*q*q]).real) if load else (a,b,c)
    g=(cc-b)/(c-b)
    sl=b+g*(s-b)
    r=np.sqrt(sl)
    dr=(g-1)*(s-b)/(r+rv)
    width=2.7*(b-a)
    z=np.clip((sl-b-.3*(b-a))/width,0.,1.)
    dm=load*betainc(9,9,z)
    # Ratios of F on the two slices, cancelling bifurcation zeros before
    # evaluation. The unsimplified expressions lose both cosmological feet.
    fr=np.ones_like(s)
    inner,outer=z==0,z==1
    shell=~(inner|outer)
    fr[inner]=(sl[inner]-a)/(s[inner]-a)*g*(c-sl[inner])/(c-s[inner])*(s[inner]/sl[inner])**2
    fr[outer]=(sl[outer]-aa)*(sl[outer]-bb)/((s[outer]-a)*(s[outer]-b))*g*(s[outer]/sl[outer])**2
    ss=s[shell]
    f0=(ss-a)*(ss-b)*(c-ss)/(2*ss**2)
    ss=sl[shell]
    fr[shell]=((ss-a)*(ss-b)*(c-ss)/(2*ss**2)-dm[shell]/ss)/f0
    jac=g*rv/r
    ds=np.log(jac/np.sqrt(fr)) if load else np.zeros_like(s)
    dsx=np.gradient(ds,step,edge_order=2)
    field=np.sqrt(6*load*(z*(1-z))**8/(width*beta(9,9)*sl))
    p=field*r**1.5*np.exp(loge/2+ds)
    # Scalar values themselves do not enter the neutral scalar equations.
    state=np.array([dr,ds,dm,dsx,-dsx,np0*(jac*np.exp(-2*ds)-1),
                    nm0*(jac*np.exp(-2*ds)-1),p,-p])
    if not load:
        state[:]=0.

    def rhs(y,geom):
        dr,ds,m,sp,sm,dnp,dnm,p,f=y
        rv,loge,np0,nm0=geom
        r=rv+dr
        e0=np.exp(loge)
        e=np.exp(loge+2*ds)
        nplus,nminus=np0+dnp,nm0+dnm
        curv=.5+3*(mu+m)/r**4-10*q*q/r**6
        curv0=.5+3*mu/rv**4-10*q*q/rv**6
        wave=e0*(np.expm1(2*ds)*curv+(curv-curv0))-p*f/(2*r**3)
        scale=np.exp(-loge/2-ds)/r
        common=e0*np.expm1(2*ds)
        # First three variables can be integrated along either ray. The
        # remaining variables each have one evolution direction.
        pos=np.array([common*np0+e*dnp,sp,-nminus*p*p/3])
        neg=np.array([common*nm0+e*dnm,sm,-nplus*f*f/3])
        adv=np.array([wave,wave,-(p*scale)**2/3,-(f*scale)**2/3,
                      -1.5*e*nplus*f/r,-1.5*e*nminus*p/r])
        return pos,neg,adv

    # S+ and P travel left, S- and M travel right. N± follow their own rays.
    left_transport=np.array([False,True,True,False,False,True])[:,None]
    geom=rv,loge,np0,nm0
    valid=np.ones(len(axis),dtype=bool)
    mid=extent
    tau=mass_integral=velocity_integral=0.
    previous_acceleration=b**.5/2-mu/b**1.5+2*q*q/b**2.5
    lapse0=np.exp(loge[mid]/2+ds[mid])
    rate0=0.
    samples=[]
    ksample=0
    previous=None
    def record(t,y,geom,tau,integral):
        rv,le,ap,am=geom
        r=rv+y[0]
        lapse=np.exp(le/2+y[1])
        U=lapse*((ap+y[5])+(am+y[6]))/2
        w=lapse*((ap+y[5])-(am+y[6]))/2
        M=r*r*(1+U*U-w*w-r*r/2)+q*q/(r*r)
        return dict(t=t,tau=tau,r=r,u=U,mass_gain=y[2],geometric_gain=M-mu,
                    stress_gain=integral,u_integral=velocity_integral,theta_plus=3*lapse*(ap+y[5])/r,
                    theta_minus=3*lapse*(am+y[6])/r)
    previous=record(0.,state[:,mid],[a[mid] for a in geom],0.,0.)
    turns=[]
    for k in range(1,extent+1):
        # Never propagate invalid parents. Disconnected components may be
        # discarded while the entire past cone of the observer stays valid.
        ok=valid[:-2]&valid[1:-1]&valid[2:]
        ids=np.flatnonzero(ok)
        mid=extent-k
        if not ok[mid]:
            break
        pu,pv=u[:-2],v[2:]
        next_geom=geometry(pu,pv)
        gg=[a[ids] for a in next_geom]
        L,R=ids,ids+2
        yl,yr=state[:,L],state[:,R]
        pl,_,al=rhs(yl,[a[L] for a in geom])
        _,nr,ar=rhs(yr,[a[R] for a in geom])
        anchor=np.vstack(((yl[:3]+yr[:3])/2,np.where(left_transport,yl[3:],yr[3:])))
        old=np.vstack(((pl+nr)/2,np.where(left_transport,al,ar)))
        pred=anchor+step*old
        def resolved(y):
            return (np.isfinite(y).all(axis=0)&(gg[0]+y[0]>.05)&(abs(y[1])<8)
                    &(step*np.maximum(abs(y[3]),abs(y[4]))<.08))
        keep=resolved(pred)
        ids,gg,anchor,old,pred=(ids[keep],[g[keep] for g in gg],anchor[:,keep],
                                old[:,keep],pred[:,keep])
        pp,nn,adj=rhs(pred,gg)
        updated=anchor+step*(old+np.vstack(((pp+nn)/2,adj)))/2
        keep=resolved(updated)
        valid=np.zeros(len(pu),dtype=bool)
        valid[ids[keep]]=True
        if not valid[mid]:
            break
        state=np.zeros((9,len(pu)))
        state[:,ids[keep]]=updated[:,keep]
        u,v,geom=pu,pv,next_geom
        y=state[:,mid]
        gpt=[a[mid] for a in geom]
        lapse=np.exp(gpt[1]/2+y[1])
        rate=-((gpt[3]+y[6])*y[7]**2+(gpt[2]+y[5])*y[8]**2)/6
        tau+=step*(lapse0+lapse)/2
        mass_integral+=step*(rate0+rate)/2
        radius=gpt[0]+y[0]
        acceleration=(radius/2-(mu+y[2])/radius**3+2*q*q/radius**5
                      -(y[7]**2+y[8]**2)/(12*radius**2*lapse**2))
        velocity_integral+=.5*step*(lapse*acceleration+lapse0*previous_acceleration)
        previous_acceleration=acceleration
        curr=record(k*step,y,gpt,tau,mass_integral)
        if previous['u']<0<=curr['u']:
            turns.append(curr)
        while ksample<=last and times[ksample]<=tau:
            w=(times[ksample]-previous['tau'])/(tau-previous['tau'])
            samples.append({key:old+w*(curr[key]-old) for key,old in previous.items()})
            ksample+=1
        previous,lapse0,rate0=curr,lapse,rate
    return dict(step=step,epsilon=epsilon,delta=delta,end_fraction=end_fraction,
                stop='optical_limit' if k==extent and valid[mid] else 'valid_stencil_limit',
                end=previous,samples=samples,turns=turns)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args()
    runs=[coupled(h,e) for e in (0.,.005,.05) for h in (.012,.006,.003)]
    runs.extend(coupled(.0015,e,end_fraction=.75) for e in (.005,.05))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(runs,indent=1,allow_nan=False)+'\n')
