"""Constraint-solved localized scalar loads on a cold electrovac bifurcation slice."""
import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import beta, betainc

from .charged_characteristics import diamond_step, wave_source
from .charged_nariai import radial_force
from .cold_loop_mass import cold_family
from .nariai_loop import gradients

OUT = Path(__file__).resolve().parents[1]/'runs/20261126_subgap_neck/probe.json'
PERIOD = 2*np.pi


def reference(gap):
    ref = cold_family(gap)
    a, b, c = ref['squared_horizons']
    q = ref['charge']
    s = brentq(lambda s: s*s*(1-s)-q*q, 0., 2/3)
    cold = 2*s-1.5*s*s
    turn = quad(lambda z: np.sqrt(2*(a+(b-a)*np.sin(z)**2)/(c-a-(b-a)*np.sin(z)**2)),
                0., np.pi/2)[0]
    return dict(**ref, cold_mass=cold, mass_gap=ref['mu']-cold, turn_tau=turn)


def prepare(n, gap, epsilon, mode, optical=False):
    if n % 4 or n < 16 or epsilon < 0 or mode not in ('standing','one_way','paired'):
        raise ValueError('invalid grid, load or initial propagation direction')
    ref = reference(gap)
    a,b,c = ref['squared_horizons']
    load = epsilon*ref['mass_gap']
    mu,q = ref['mu'],ref['charge']
    # The loaded cosmological turning sphere is the outermost root of F.
    roots = np.roots([1.,-2.,2*(mu+load),-2*q*q])
    if np.any(abs(roots.imag)>1e-12):
        raise ValueError('outer vacuum lost its three-horizon branch')
    aa,bb,cc = sorted(roots.real) if load else (a,b,c)
    x = np.arange(n)*PERIOD/n
    # The cropped chart resolves exponentially small distances from the
    # cosmological bifurcations while preserving the same geometric slice.
    if optical:
        coordinate = (np.arange(n)-n//2)*16./n
        argument = 1e-4*np.sinh(4*coordinate)
        y = 4*np.arcsinh(argument)
        slope = .0016*np.cosh(4*coordinate)/np.hypot(1.,argument)
        sin_x,cos_x = np.tanh(y),1/np.cosh(y)
        jac = cos_x*slope
        parity = np.where(coordinate>=0,1.,-1.)
        dx = 16./n
    else:
        sin_x,cos_x,jac = np.sin(x),np.cos(x),np.ones(n)
        parity = np.where(x<np.pi,1.,-1.)
        dx = PERIOD/n
    s = b+(cc-b)*sin_x**2
    base_rx = (cc-b)*sin_x*cos_x/np.sqrt(s)
    rx = base_rx*jac
    lo,width = b+.3*(b-a), 2.7*(b-a)
    if not b < lo < lo+width < cc:
        raise ValueError('shell extends outside static interval')
    z = np.clip((s-lo)/width,0.,1.)
    m = mu+load*betainc(9,9,z)
    dm_ds = load*(z*(1-z))**8/(width*beta(9,9))
    radius = np.sqrt(s)
    # Cancel the endpoint factors analytically at both bifurcation spheres.
    scale2 = np.empty(n)
    inside, outside = z==0, z==1
    shell = ~(inside|outside)
    scale2[inside] = 2*s[inside]*(cc-s[inside])/((s[inside]-a)*(c-s[inside]))
    scale2[outside] = 2*s[outside]*(s[outside]-b)/((s[outside]-aa)*(s[outside]-bb))
    metric = (s-a)*(s-b)*(c-s)/(2*s*s)-(m-mu)/s
    scale2[shell] = base_rx[shell]**2/metric[shell]
    if np.any(scale2 <= 0):
        raise ValueError('loaded initial slice is not spacelike')
    scale = np.sqrt(scale2)*jac
    w = rx/scale
    f = parity*np.sign(rx)*np.sqrt(6*dm_ds/s)
    p = np.zeros(n)
    if mode != 'standing':
        f /= np.sqrt(2.)
        p = f.copy() if mode=='one_way' else parity*f
    h = np.zeros(n)
    h[shell] = radius[shell]*p[shell]*f[shell]/(3*w[shell])
    sigma = np.log(scale)
    field = cumulative_simpson(f*scale,dx=dx,initial=0.)
    pos = np.array([radius,sigma,field])
    time = np.array([np.zeros(n),scale*h,scale*p])
    space = np.array([rx,gradients(sigma,dx)[0],scale*f])
    ham = (3*((1-w*w)/s-gradients(w,dx)[0]/(scale*radius))-3
           -3*q*q/s**3-(p*p+f*f)/2)
    mom = -h*w+radius*p*f/3
    geom_mass = s*(1-w*w-s/2)+q*q/s
    diagnostics = dict(load=load, recovered_load=float(np.ptp(geom_mass)),
                       mass_error=float(max(abs(geom_mass-m))),
                       energy_each_direction=[float(dx*sum(scale*radius**3*v*v/4))
                                              for v in (p+f,p-f)],
                       scalar_shift_charge=float(dx*sum(scale*radius**3*p)),
                       ham_max=float(max(abs(ham))), mom_max=float(max(abs(mom))),
                       signed_half_shift_h_error=float(max(abs(h+np.roll(h,n//2)))),
                       null_generation_rate=float(max(abs(wave_source(pos,time+space,time-space,q)[2]))))
    return pos,time+space,time-space,diagnostics


def focused_step(pos, plus, minus, valid, dt, q):
    """Integrate both null focusing constraints along their own null rays."""
    j = np.flatnonzero(valid & np.roll(valid,1) & np.roll(valid,-1))
    l,r = (j-1)%len(valid),(j+1)%len(valid)
    left,right = pos[:,l],pos[:,r]
    pl,ml,pr,mr = plus[:,l],minus[:,l],plus[:,r],minus[:,r]
    sl,sr = wave_source(left,pl,ml,q),wave_source(right,pr,mr,q)
    avg,incident = (left+right)/2,pl+mr
    zp = avg+dt*incident/2
    pp,mm = pr+dt*sr,ml+dt*sl
    # Propagate e^(-2 sigma) D_± r without cancellation between an
    # exponentially small coordinate derivative and an additive wave source.
    def radial(start,d,z,increment):
        weight = np.exp(2*(z[1]-start[1]))
        return weight*(d[0]-increment*start[0]*d[2]**2/3)
    def resolved(z,p,m):
        return (np.isfinite(z).all(axis=0) & np.isfinite(p).all(axis=0)
                & np.isfinite(m).all(axis=0) & (z[0]>.025)
                & (z[1]<np.log(150)) & (np.abs(z[1])<300)
                & (dt*np.maximum(abs(p[0]),abs(m[0]))<.08*z[0])
                & (dt*np.maximum(np.max(abs(p[1:]),axis=0),np.max(abs(m[1:]),axis=0))<.08))
    pp[0],mm[0] = radial(left,pl,zp,dt),radial(right,mr,zp,dt)
    ok = resolved(zp,pp,mm)
    j,zp,pp,mm = j[ok],zp[:,ok],pp[:,ok],mm[:,ok]
    left,right,pl,mr,pr,ml,sl,sr,avg,incident = [
        v[:,ok] for v in (left,right,pl,mr,pr,ml,sl,sr,avg,incident)]
    src = wave_source(zp,pp,mm,q)
    pp_new,mm_new = pr+dt*(sr+src)/2, ml+dt*(sl+src)/2
    # The corrector uses the predicted endpoint in both quadratures.
    new_z = avg+dt*(incident+pp+mm)/4
    for dest,start,d,pred in ((pp_new,left,pl,pp),(mm_new,right,mr,mm)):
        dest[0] = (radial(start,d,new_z,dt/2)
                   -dt*np.exp(2*(new_z[1]-zp[1]))*zp[0]*pred[2]**2/6)
    # Use corrected null derivatives for the final position and use its
    # conformal factor for the already integrated, affine-weighted expansions.
    final = avg+dt*(incident+pp_new+mm_new)/4
    weight = np.exp(2*(final[1]-new_z[1]))
    pp_new[0] *= weight
    mm_new[0] *= weight
    ok = resolved(final,pp_new,mm_new)
    return j[ok],final[:,ok],pp_new[:,ok],mm_new[:,ok]


def interpolate(array, x, dx):
    cell = math.floor(x/dx)
    j = cell % array.shape[-1]
    frac = x/dx-cell
    return array[...,j]+frac*(array[...,(j+1)%array.shape[-1]]-array[...,j])


def observer(state, pos, plus, minus, dx, q):
    x,eta,tau = state
    r,s,_ = interpolate(pos,x,dx)
    pr,ps,pf = interpolate(plus,x,dx)
    mr,ms,mf = interpolate(minus,x,dx)
    v = np.tanh(eta)
    ds = -.5*((ps-ms)+v*(ps+ms))
    clock = np.exp(s)/np.cosh(eta)
    u = .5*np.exp(-s)*(np.exp(eta)*pr+np.exp(-eta)*mr)
    mass = r*r*(1+pr*mr*np.exp(-2*s)-r*r/2)+q*q/r**2
    w = .5*np.exp(-s)*(np.exp(eta)*pr-np.exp(-eta)*mr)
    flux = ((pf*np.exp(eta-s))**2/4,(mf*np.exp(-eta-s))**2/4)
    rate = -2*r**3*((u-w)*flux[0]+(u+w)*flux[1])/3
    return np.array([v,ds,clock]), dict(x=float(x),tau=float(tau),r=float(r),u=float(u),
                                       mass=float(mass),mass_rate=float(rate),flux_plus=float(flux[0]),flux_minus=float(flux[1]))


def run(case):
    n,eps,mode,gap,*options = case
    optical = bool(options and options[0])
    ref = reference(gap)
    q = ref['charge']
    pos,plus,minus,initial = prepare(n,gap,eps,mode,optical)
    if optical:
        initial = prepare(n,gap,eps,mode)[3]
    valid = np.ones(n,dtype=bool)
    dx = 16./n if optical else PERIOD/n
    obs = np.array([8. if optical else 0.,0.,0.])
    if optical:
        valid[:3] = valid[-3:] = False
    evo,point = observer(obs,pos,plus,minus,dx,q)
    clocks = iter(ref['turn_tau']*np.array([.25,.45,.5,.54,.55,.57,.59,.60,.62,.65,.7,.75,.9,.95,.98,.99,.995,.999,1.,1.005,1.01]))
    target = next(clocks)
    samples,turns,frames = [],[],[]
    stop = 'cropped_domain' if optical else 'one_coordinate_period'
    def constraint():
        j = round(obs[0]/dx)%n
        if not valid[[(j-1)%n,j,(j+1)%n]].all():
            return None
        sr = wave_source(pos[:,j],plus[:,j],minus[:,j],q)[0]
        err = []
        for d,sign in ((plus,1),(minus,-1)):
            terms = [sr,sign*(d[0,(j+1)%n]-d[0,(j-1)%n])/dx,
                     -2*d[1,j]*d[0,j],pos[0,j]*d[2,j]**2/3]
            err.append(abs(sum(terms))/(1+sum(abs(y) for y in terms)))
        return float(max(err))
    point.update(t=0.,constraint=constraint(),stress_gain=0.)
    frames.append(point)
    for k in range(1,n+1):
        if optical:
            active,z,p,m = focused_step(pos,plus,minus,valid,dx,q)
        else:
            active,z,p,m = diamond_step(pos,plus,minus,valid,dx,q,150.)
        valid[:] = False
        valid[active] = True
        # A timelike trajectory can move at most one cell during one null step.
        j = math.floor(obs[0]/dx)
        if not valid[np.arange(j-1,j+3)%n].all():
            stop = 'observer_stencil_limit'
            break
        pos[:,active],plus[:,active],minus[:,active] = z,p,m
        predicted = obs+dx*evo
        new_evo,_ = observer(predicted,pos,plus,minus,dx,q)
        obs = obs+dx*(evo+new_evo)/2
        evo,new = observer(obs,pos,plus,minus,dx,q)
        new.update(t=k*dx,constraint=constraint(),
                   stress_gain=point['stress_gain']+.5*(new['tau']-point['tau'])*(point['mass_rate']+new['mass_rate']))
        while target is not None and new['tau']>=target:
            weight=(target-point['tau'])/(new['tau']-point['tau'])
            samples.append({key:float(value+weight*(new[key]-value)) for key,value in point.items()})
            target=next(clocks,None)
        if point['u']<0<=new['u']:
            turns.append(new)
        if k%max(1,n//128)==0:
            frames.append(new)
        point = new
    if point is not frames[-1]:
        frames.append(point)
    return dict(n=n,epsilon=eps,mode=mode,gap=gap,optical=optical,reference=ref,initial=initial,
                stop=stop,samples=samples,frames=frames,turns=turns)


def vacuum_foot(ref,r):
    """Logarithmic distance of the initial light-ray foot to the cosmological sphere.

    On the central orbit the two outer Kruskal coordinates have equal magnitude.
    On the time-symmetric exterior this fixes the same real tortoise coordinate.
    Use log distance to evaluate the exterior arbitrarily close to r_c.
    """
    ra,rb,rc=np.sqrt(ref['squared_horizons'])
    radii=np.array([ra,rb,rc])
    coeff=np.array([-1.,1.,-1.])/(2*np.array(ref['kappa']))
    interior=sum(coeff*np.log(abs((r-radii)/(r+radii))))
    def residual(log_dist):
        dist=(rc-rb)*np.exp(log_dist)
        rr=rc-dist
        logs=np.array([np.log((rr-ra)/(rr+ra)),np.log((rr-rb)/(rr+rb)),
                       log_dist+np.log(rc-rb)-np.log(rr+rc)])
        return sum(coeff*logs)-interior
    return brentq(residual,-10000.,-1e-12)/np.log(10.)


def cosmological_collar(ref, epsilon):
    """Past light feet of the untouched expanding electrovac observer at future infinity."""
    load = epsilon*ref['mass_gap']
    mu,q = ref['mu']+load,ref['charge']
    squares = np.sort(np.roots([1.,-2.,2*mu,-2*q*q]).real)
    radii = np.sqrt(squares)
    coeff = 1/(-2*radial_force(radii,mu,q))
    def tortoise(r):
        return np.dot(coeff,np.log(abs((r-radii)/(r+radii))))
    rb,rc = radii[1:]
    root = brentq(tortoise,rb+1e-10,rc-1e-10)
    a,b,_ = ref['squared_horizons']
    shell_end = np.sqrt(b+3*(b-a))
    return dict(outer_mass=mu,cosmological_radius=rc,shell_end=shell_end,
                limiting_initial_foot=root,collar_margin=root-shell_end,
                finite_radius_feet=[brentq(lambda r: tortoise(r)-tortoise(R),
                                            rb+1e-10,rc-1e-10) for R in (2.,10.,100.)])


def exact_clock(ref,times):
    r0=np.sqrt(ref['squared_horizons'][1])
    rhs=lambda _,y:(y[1],radial_force(y[0],ref['mu'],ref['charge']))
    sol=solve_ivp(rhs,(0,max(times)),(r0,0.),method='DOP853',t_eval=times,rtol=2e-13,atol=1e-15)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T.tolist()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=1)
    args=parser.parse_args()
    cases=[(n,eps,mode,.01) for eps,mode in ((0.,'standing'),(.005,'standing'),(.05,'standing'),
                                               (1.,'standing'),(.05,'one_way'),(.05,'paired'))
           for n in (2048,4096,8192)]
    cases += [(n,eps,mode,.01,True)
              for eps,mode in ((0.,'standing'),(.005,'standing'),(.05,'standing'),
                               (1.,'standing'),(.05,'one_way'),(.05,'paired'))
              for n in (4096,8192,16384)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows=list(pool.map(run,cases))
    ref=reference(.01)
    times=(ref['turn_tau']*np.array([.9,.95,.98,.99,.995,.999,1.,1.005,1.01])).tolist()
    OUT.parent.mkdir(parents=True,exist_ok=True)
    orbit=exact_clock(ref,times)
    feet=[vacuum_foot(ref,r) for time,(r,u) in zip(times,orbit) if time<ref['turn_tau']]
    OUT.write_text(json.dumps(dict(vacuum_orbit=dict(tau=times,ru=orbit,pre_turn_log10_foot=feet),
                                   collars=[dict(epsilon=e,**cosmological_collar(ref,e))
                                            for e in (0.,.005,.05,1.)],
                                   runs=rows),allow_nan=False,indent=1)+'\n')


if __name__=='__main__':
    main()
