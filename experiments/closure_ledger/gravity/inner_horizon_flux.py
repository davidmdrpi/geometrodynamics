"""Leading scalar mass response on the exact 5D cold electrovac bridge."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import beta, betainc, expit

from .subgap_neck import reference, exact_clock

OUT = Path(__file__).resolve().parents[1]/'runs/20261127_inner_horizon_flux/probe.json'
FRACTIONS = (.25,.5,.6,.75,.9,.95,.98,.99,.995,.999,.9995,.9999)
MODES = ('standing','one_way')
DELTA = .0001


def background(ref, signed_log):
    """Invert signed H = ±exp(2 k_b r*) on each side of the outer horizon.

    Input is asinh(H). The logit radial variable retains each small horizon
    distance even after the areal radius rounds to its limiting value.
    """
    ra,rb,rc = np.sqrt(ref['squared_horizons'])
    ka,kb,kc = ref['kappa']
    exterior = signed_log>=0
    other = np.where(exterior,rc,ra)
    width = np.abs(other-rb)
    h = np.sinh(signed_log)
    log_h = np.log(np.maximum(abs(h),1e-300))
    low,high = np.full_like(h,-1400.), np.full_like(h,1400.)
    for _ in range(60):
        l = (low+high)/2
        r = rb+(other-rb)*expit(l)
        log_b = np.log(width)-np.logaddexp(0.,-l)
        far = np.log(width)-np.logaddexp(0.,l)
        log_a = np.where(exterior,np.log(np.maximum(r-ra,1e-300)),far)
        log_c = np.where(exterior,far,np.log(np.maximum(rc-r,1e-300)))
        la,lc = log_a-np.log(r+ra),log_c-np.log(r+rc)
        found = log_b-np.log(r+rb)-kb/ka*la-kb/kc*lc
        low = np.where(found<log_h,l,low)
        high = np.where(found>=log_h,l,high)
    ratio = (log_a+np.log(r+ra)+2*np.log(r+rb)+log_c+np.log(r+rc)
             -np.log(2.)-4*np.log(r)+kb/ka*la+kb/kc*lc)
    return r,ratio


def events(ref,delta):
    radii=np.sqrt(ref['squared_horizons'])
    ka,kb,kc=ref['kappa']
    tau=ref['turn_tau']*np.array(FRACTIONS)
    orbit=np.array(exact_clock(ref,tau))
    l=np.log(abs((orbit[:,0,None]-radii)/(orbit[:,0,None]+radii)))
    log_h=l[:,1]-kb/ka*l[:,0]-kb/kc*l[:,2]
    return np.arcsinh(np.exp(log_h/2)/delta),tau,orbit


def transport(step, epsilon=.005, delta=DELTA, table_step=8e-5):
    ref=reference(.01)
    target,tau,orbit=events(ref,delta)
    limit=int(np.ceil(target[-1]/step))+2
    kb=ref['kappa'][1]
    coord=step*np.arange(-limit,limit+1)
    K=delta*np.sinh(coord)
    Kp=delta*np.cosh(coord)
    # Use a uniform signed-log table to avoid root solves in every diamond.
    edge=np.arcsinh(K[-1]**2)
    m=int(np.ceil(edge/table_step))
    table=table_step*np.arange(-m,m+1)
    radial,log_ratio=background(ref,table)
    def metric(u,v):
        value=np.arcsinh(-K[u]*K[v])/table_step+m
        idx=np.clip(value.astype(int),0,len(table)-2)
        weight=value-idx
        r=radial[idx]+weight*(radial[idx+1]-radial[idx])
        ratio=np.exp(log_ratio[idx]+weight*(log_ratio[idx+1]-log_ratio[idx]))
        return r,ratio
    u=np.arange(2*limit,-1,-1)
    v=np.arange(2*limit+1)
    r,ratio=metric(u,v)
    a,b,_=ref['squared_horizons']
    width=2.7*(b-a)
    z=np.clip((r*r-b-.3*(b-a))/width,0.,1.)
    load=epsilon*ref['mass_gap']
    f=np.sqrt(6*load*(z*(1-z))**8/(width*beta(9,9)*r*r))
    norm=r**1.5*np.sqrt(ratio)*Kp/kb
    plus=np.array([norm*f,np.sqrt(2.)*norm*f])
    minus=np.array([-norm*f,np.zeros_like(f)])
    initial=load*betainc(9,9,z)
    mass_right=np.tile(initial,(2,1))
    mass_left=mass_right.copy()
    def coefficients(u,v,r,ratio):
        return (1.5*ratio*K[u]*Kp[v]/(kb*r),1.5*ratio*Kp[u]*K[v]/(kb*r))
    cp,cm=coefficients(u,v,r,ratio)
    factor=kb*np.tanh(coord)/3
    accum=np.zeros(2)
    rate=np.zeros(2)
    # First variation of the timelike radial equation at fixed proper time.
    dr=np.zeros(2)
    du=np.zeros(2)
    force=np.zeros(2)
    kinetic=np.zeros(2)
    lapse=np.sqrt(ratio[limit])*Kp[limit]/kb
    clock=0.
    samples=[]
    smallest_det=1.
    next_event=0
    for k in range(1,limit+1):
        # The two parents originate on the supplied slice. The domain shrinks
        # one cell on each side at each step and never reads a boundary value.
        pu, pv = u[:-2], v[2:]
        r,ratio=metric(pu,pv)
        new_cp,new_cm=coefficients(pu,pv,r,ratio)
        aplus,aminus=.5*step*new_cp,.5*step*new_cm
        det=1-aplus*aminus
        smallest_det=min(smallest_det,float(min(det)))
        if smallest_det<=0:
            raise ArithmeticError('unresolved characteristic coupling')
        p0=plus[:,2:]+.5*step*cp[2:]*minus[:,2:]
        m0=minus[:,:-2]+.5*step*cm[:-2]*plus[:,:-2]
        p=(p0+aplus*m0)/det
        mfield=m0+aminus*p
        # D+ mass integrates along a right-going ray (left parent), whereas
        # D- mass integrates along the opposite ray. Each includes its own
        # initial shell mass and is independent of the observer quadrature.
        right=(mass_right[:,:-2]+.5*step*(factor[v[:-2]]*plus[:,:-2]**2+factor[pv]*p**2))
        left=(mass_left[:,2:]+.5*step*(factor[u[2:]]*minus[:,2:]**2+factor[pu]*mfield**2))
        mid=limit-k
        old_rate=rate
        rate=.5*factor[limit+k]*(p[:,mid]**2+mfield[:,mid]**2)
        old_accum=accum
        accum=accum+.5*step*(rate+old_rate)
        old_lapse=lapse
        lapse=np.sqrt(ratio[mid])*Kp[limit+k]/kb
        old_clock=clock
        clock+=.5*step*(old_lapse+lapse)
        radius=r[mid]
        curvature=.5+3*ref['mu']/radius**4-10*ref['charge']**2/radius**6
        old_kinetic=kinetic
        kinetic=(p[:,mid]**2+mfield[:,mid]**2)/(2*radius**3*lapse**2)
        source=-accum/radius**3-radius*kinetic/6
        old_dr,old_du=dr,du
        half=.5*step*lapse
        d0=dr+.5*step*old_lapse*du
        d1=du+.5*step*old_lapse*force+half*source
        dr=(d0+half*d1)/(1-half*half*curvature)
        du=d1+half*curvature*dr
        force=curvature*dr+source
        while next_event<len(target) and target[next_event]<=k*step:
            j=next_event
            w=(target[j]-(k-1)*step)/step
            measures=[]
            for i,mode in enumerate(MODES):
                gain=(1-w)*old_accum[i]+w*accum[i]
                rr=(1-w)*mass_right[i,mid+1]+w*right[i,mid]
                ll=(1-w)*mass_left[i,mid+1]+w*left[i,mid]
                delta_r=(1-w)*old_dr[i]+w*dr[i]
                delta_u=(1-w)*old_du[i]+w*du[i]
                rv,uv=orbit[j]
                derivative=2*rv*(1+uv**2-rv**2)-2*ref['charge']**2/rv**3
                geometric=derivative*delta_r+2*rv**2*uv*delta_u
                measures.append(dict(mode=mode,mass_gain=gain,over_gap=gain/ref['mass_gap'],
                                     from_right_going=rr,from_left_going=ll,
                                     delta_r=delta_r,delta_u=delta_u,geometric_variation=geometric,
                                     scalar_density=((1-w)*old_kinetic[i]+w*kinetic[i])/2))
            samples.append(dict(fraction=FRACTIONS[j],tau=tau[j],vacuum_r=orbit[j,0],
                                vacuum_u=orbit[j,1],integrated_tau=(1-w)*old_clock+w*clock,values=measures))
            next_event+=1
        u,v=pu,pv
        cp,cm=new_cp,new_cm
        plus,minus=p,mfield
        mass_right,mass_left=right,left
        if next_event==len(target):
            break
    return dict(step=step,epsilon=epsilon,delta=delta,table_step=table_step,
                minimum_determinant=smallest_det,samples=samples)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args()
    rows=[transport(h) for h in (.016,.008,.004,.002)]
    rows.append(transport(.004,delta=DELTA/2,table_step=4e-5))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(reference=reference(.01),runs=rows),
                                      allow_nan=False,indent=1)+'\n')


if __name__=='__main__':
    main()
