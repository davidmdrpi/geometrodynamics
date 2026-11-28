"""J-equivariant 5D Einstein-Lambda loop with round three-sphere orbits."""
import json

import numpy as np


def vacuum(t, x, eps):
    """Gaussian development from a time-symmetric double-cover Kottler slice."""
    s = 1-eps*np.cos(2*x)
    mu = (1-eps**2)/2
    b = (s-2*mu/s)/2
    hyper = np.cosh(np.sqrt(2)*t)
    r = np.sqrt(s+b*(hyper-1))
    a = np.sqrt(2*s)/r*(1+(1+2*mu/s**2)*(hyper-1)/2)
    velocity = b*np.sinh(np.sqrt(2)*t)/(np.sqrt(2)*r)
    spatial = eps*np.sin(2*x)/np.sqrt(2*s)  # R_x/a at fixed T
    return r, a, velocity, spatial


def vacuum_probe(eps=.1, tilt=.25, size=256):
    x = np.arange(size)*(2*np.pi/size)
    width = 2.
    plus = np.exp(width*(np.cos(x)-1))
    minus = np.exp(width*(-np.cos(x)-1))
    t = tilt*(plus-minus)/(1-np.exp(-2*width))
    slope = -tilt*width*np.sin(x)*(plus+minus)/(1-np.exp(-2*width))
    r,a,v,w = vacuum(t,x,eps)
    beta = slope/a
    boost = 1/np.sqrt(1-beta**2)
    u = boost*(v+beta*w)
    rx = boost*(w+beta*v)
    mu = (1-eps**2)/2
    mass = r**2*(1+u**2-rx**2-r**2/2)
    # Longitudinal normal expansion on the tilted slice: derivative of the
    # boost rapidity along its unit tangent plus the boosted background expansion.
    def spatial_derivative(z):
        return boost/a*gradients(z,2*np.pi/size)[0]
    s0 = 1-eps*np.cos(2*x)
    factor = (1+2*mu/s0**2)/2
    hubble_a = -v/r + factor*np.sqrt(2)*np.sinh(np.sqrt(2)*t)/(1+factor*(np.cosh(np.sqrt(2)*t)-1))
    longitudinal = spatial_derivative(np.arctanh(beta))+boost*hubble_a
    ham = 3*(1+u*u-rx*rx)/r**2+3*(longitudinal*u-spatial_derivative(rx))/r-3
    mom = spatial_derivative(u)-longitudinal*rx
    pinch = np.arccosh(1/eps)/np.sqrt(2)
    return dict(eps=eps, tilt=tilt, maximum_tilt_speed=float(max(abs(beta))),
                mass_parameter=mu, max_mass_defect=float(max(abs(mass-mu))),
                max_hamiltonian=float(max(abs(ham))), max_momentum=float(max(abs(mom))),
                max_longitudinal_pair_error=float(max(abs(longitudinal+np.roll(longitudinal,size//2)))),
                max_radius_pair_error=float(max(abs(r-np.roll(r,size//2)))),
                max_velocity_pair_error=float(max(abs(u+np.roll(u,size//2)))),
                centers=[dict(x=float(x[i]),initial_radius=float(r[i]),normal_velocity=float(u[i]),
                              theta_plus=float(3*(u[i]+rx[i])/r[i]),
                              theta_minus=float(3*(u[i]-rx[i])/r[i]),
                              proper_time_to_future_pinch=float(pinch-t[i]))
                         for i in (0,size//2)],
                near_pinch=[dict(radius=cut, kretschmann=10+72*mu**2/cut**8,
                                 proper_time_left=float(pinch-np.arccosh((1-cut**2)/eps)/np.sqrt(2)))
                            for cut in (.3,.2,.1)])


def gradients(y, dx):
    l2,l1,r1,r2 = [np.roll(y,k,axis=-1) for k in (2,1,-1,-2)]
    return ((l2-8*l1+8*r1-r2)/(12*dx),
            (-l2+16*l1-30*y+16*r1-r2)/(12*dx**2))


def pulse_profile(n,length):
    phase = np.arange(n)*(2*np.pi/n)
    g = np.exp((np.cos(2*phase)-1)/4)
    f = np.cos(phase)*g
    df = (2*np.pi/length)*(-g*np.sin(phase)-f*np.sin(2*phase)/2)
    return f,df


def initial(n, amplitude, radius=.98, length=2*np.pi):
    f,df = pulse_profile(n,length)
    u = -amplitude*f**3
    p = np.sqrt(6*(1/radius**2-1))
    field = -3*u/(radius*p)
    # Smooth extension of U_x²/U at the zeros of U.
    longitudinal = -u/radius-27*amplitude*f*df**2/(2*p*p*radius)
    return np.array([np.full(n,radius),np.zeros(n),field,u,longitudinal,np.full(n,p)])


def sources(y, dx):
    r,s,phi,v,h,p = y
    first,second = gradients(y[:3],dx)
    rx,sx,px = first
    a2 = np.exp(2*s)
    wave_r = 2*a2*r-2*(a2+v*v-rx*rx)/r
    wave_s = 2*a2-1.5*wave_r/r-(p*p-px*px)/2
    wave_phi = -3*(v*p-rx*px)/r
    acc = second+np.array([wave_r,wave_s,wave_phi])
    return first, second, acc


def rhs(y,dx):
    return np.concatenate((y[3:],sources(y,dx)[2]))


def rk4(y,dx,dt,equations,clock_rate=None):
    k1 = equations(y,dx)
    stage2 = y+dt*k1/2
    k2 = equations(stage2,dx)
    stage3 = y+dt*k2/2
    k3 = equations(stage3,dx)
    stage4 = y+dt*k3
    k4 = equations(stage4,dx)
    updated = y+dt*(k1+2*k2+2*k3+k4)/6
    if clock_rate is None:
        return updated
    elapsed = dt*(clock_rate(y,0)+2*clock_rate(stage2,.5)
                  +2*clock_rate(stage3,.5)+clock_rate(stage4,1))/6
    return updated,elapsed


def periodic_sample(arr,position,dx):
    pos = (position/dx)%len(arr)
    j = int(pos)
    return float(arr[j]+(pos-j)*(arr[(j+1)%len(arr)]-arr[j]))


def measure(y,dx,t):
    r,s,_,v,h,p = y
    first,second,acc = sources(y,dx)
    rx,sx,px = first
    vx = gradients(v,dx)[0]
    errors = []
    for sign in (-1,1):
        terms = [acc[0]+second[0], 2*sign*vx, -2*(h+sign*sx)*(v+sign*rx),
                 r*(p+sign*px)**2/3]
        denom = 1+sum(abs(term) for term in terms)
        errors.append(float(max(abs(sum(terms))/denom)))
    def point(i):
        boost = np.exp(-s[i])
        return dict(r=float(r[i]), normal_velocity=float(boost*v[i]),
                    theta_plus=float(3*boost*(v[i]+rx[i])/r[i]),
                    theta_minus=float(3*boost*(v[i]-rx[i])/r[i]),
                    angular_curvature=float((1+boost**2*(v[i]**2-rx[i]**2))/r[i]**2))
    def interpolate(arr):
        return periodic_sample(arr,t,dx)
    ray = dict(r=interpolate(r), affine_velocity=interpolate(np.exp(-2*s)*(v+rx)),
               affine_rate=interpolate(np.exp(2*s)))
    return dict(t=float(t), min_r=float(min(r)), peak_lapse=float(max(np.exp(s))),
                centers=[point(0),point(len(r)//2)], null_constraint=max(errors), ray=ray)


def evolve(n, amplitude, cut=.2):
    dx = 2*np.pi/n
    y = initial(n,amplitude)
    t = 0.
    clocks = np.zeros(3)
    def frame():
        return dict(measure(y,dx,t), center_proper_times=clocks[:2].tolist(),ray_affine=float(clocks[2]))
    frames = [frame()]
    # Fixed conformal-coordinate frames for comparing refinement, plus events.
    checkpoints = iter(np.arange(.1,2.51,.1))
    target = next(checkpoints)
    turning = None
    while True:
        dt = min(dx/6, .025*np.min(y[0]/(abs(y[3])+1)),
                 .025/(1+np.max(abs(y[4]))), target-t)
        def rates(stage,fraction):
            return np.array([np.exp(stage[1,0]),np.exp(stage[1,n//2]),
                             periodic_sample(np.exp(2*stage[1]),t+fraction*dt,dx)])
        new,elapsed = rk4(y,dx,dt,rhs,rates)
        if not np.all(np.isfinite(new)):
            raise FloatingPointError("nonfinite before declared stopping surface")
        def reached(z):
            return min(z[0]) <= cut or max(z[1]) >= np.log(50)
        stop = reached(new)
        if stop:
            lo,hi = 0.,dt
            for _ in range(30):
                m = (lo+hi)/2
                if reached(rk4(y,dx,m,rhs)):
                    hi = m
                else:
                    lo = m
            dt = hi
            new,elapsed = rk4(y,dx,dt,rhs,rates)
        i = n//2
        if turning is None and y[3,i]>0>=new[3,i]:
            lo,hi = 0.,dt
            for _ in range(25):
                m = (lo+hi)/2
                if rk4(y,dx,m,rhs)[3,i]>0:
                    lo=m
                else:
                    hi=m
            turning = dict(t=float(t+hi),r=float(rk4(y,dx,hi,rhs)[0,i]))
        clocks += elapsed
        y,t = new,t+dt
        if stop or abs(t-target)<1e-12:
            frames.append(frame())
            if stop:
                reason = "radius_cutoff" if min(y[0])<=cut+1e-9 else "lapse_cutoff"
                break
            target = next(checkpoints,None)
            if target is None:
                reason = "time_limit"
                break
    return dict(n=n, amplitude=amplitude, cut=cut, stop=reason, antipode_maximum=turning,
                frames=frames)


def gaussian_rhs(y,dx):
    """Synchronous proper time; evolve radial gradient via momentum constraint."""
    r,a,u,h,w,p,f = y
    px,fx = gradients(y[5:],dx)[0]
    acc = r-(1+u*u-w*w)/r-r*(p*p+f*f)/6
    return np.array([u,a*h,acc,2-p*p-3*acc/r-h*h,-r*p*f/3,
                     fx/a + 3*w*f/r-(h+3*u/r)*p,px/a-h*f])


def gaussian(n, amplitude=.05, cut=.2, compare_at=(), length=2*np.pi):
    dx = length/n
    r,s,field,u,h,p = initial(n,amplitude,length=length)
    x = np.arange(n)*dx
    f0,df = pulse_profile(n,length)
    f = 9*amplitude*f0**2*df/(r*p)
    y = np.array([r,np.ones(n),u,h,np.zeros(n),p,f])
    t = 0.
    frames = []
    # Future light signals from the two centers headed toward one another.
    light = np.array([0.,length/2])
    arrivals = [None,None]
    points = iter(sorted(set(np.arange(.1,5.01,.1).tolist()).union(compare_at)))
    target = next(points)
    def snapshot():
        r,a,u,h,w,p,f = y
        rx,ux,wx = gradients(y[[0,2,4]],dx)[0]/a
        # H_nn: 3[(1+u²−w²)/r² + (hu−w_x/A)/r] = Λ + rho
        geom = 3*(1+u*u-w*w)/r**2
        bend = 3*(h*u-wx)/r
        rho = (p*p+f*f)/2
        residual = abs(geom+bend-3-rho)/(1+abs(geom)+abs(bend)+3+rho)
        mom = abs(ux-h*w+r*p*f/3)/(1+abs(ux)+abs(h*w)+abs(r*p*f/3))
        i = int(np.argmin(r))
        def record(j):
            return dict(radius=float(r[j]),velocity=float(u[j]),theta_plus=float(3*(u[j]+w[j])/r[j]),
                        theta_minus=float(3*(u[j]-w[j])/r[j]),
                        angular_curvature=float((1+u[j]**2-w[j]**2)/r[j]**2))
        return dict(t=float(t), light_positions=light.tolist(), min_x=float(x[i]), minimum=record(i), centers=[record(0),record(n//2)],
                    min_a=float(min(a)), max_hamiltonian=float(max(residual)),
                    max_momentum=float(max(mom)), max_gradient_defect=float(max(abs(rx-w))))
    frames.append(snapshot())
    def reached(z):
        return min(z[0])<=cut or min(z[1])<=.04
    while True:
        dt=min(dx*min(y[1])/6,.02*np.min(y[0]/(1+abs(y[2]))),.02/(1+max(abs(y[3]))),target-t)
        new=rk4(y,dx,dt,gaussian_rhs)
        if not np.all(np.isfinite(new)):
            raise FloatingPointError("nonfinite synchronous evolution")
        stop=reached(new)
        if stop:
            lo,hi=0.,dt
            for _ in range(30):
                m=(lo+hi)/2
                if reached(rk4(y,dx,m,gaussian_rhs)):
                    hi=m
                else:
                    lo=m
            dt=hi
            new=rk4(y,dx,dt,gaussian_rhs)
        def coordinate_speed(a, pos):
            return np.array([1/periodic_sample(a,pos[0],dx),-1/periodic_sample(a,pos[1],dx)])
        # Trapezoid in proper time for the diagnostic rays; separated from
        # the fourth-order evolution of metric and matter.
        speed = coordinate_speed(y[1],light)
        travel = dt*(speed+coordinate_speed(new[1],light+dt*speed))/2
        for k, destination in enumerate((length/2,0.)):
            fraction = (destination-light[k])/travel[k]
            if arrivals[k] is None and 0 <= fraction <= 1:
                arrivals[k] = float(t+dt*fraction)
        light += travel
        t,y=t+dt,new
        if stop or abs(t-target)<1e-12:
            frames.append(snapshot())
            if stop or t>=5.:
                break
            target = next(points)
    # At the reflection center w=0 and (r²)'' <= 2(r²-1).
    # Solve the zero of its homogeneous upper solution using z=exp(sqrt(2)*tau).
    y0 = r[0]**2
    linear = np.sqrt(2)*r[0]*u[0]
    c = 1-y0-linear
    upper_time = np.log((1+np.sqrt(1-c*(1-y0+linear)))/c)/np.sqrt(2)
    return dict(n=n,amplitude=amplitude,length=length,cut=cut,frames=frames,
                center_comparison_time=float(upper_time),
                center_signal_arrivals=arrivals,
                stop="radius_cutoff" if min(y[0])<=cut+1e-9 else
                     ("caustic_guard" if min(y[1])<=.04+1e-9 else "time_limit"))


def probe():
    conformal = [evolve(n,.05) for n in (128,256,512)]
    conformal += [evolve(256,a) for a in (0.,.025,.10)]
    conformal += [evolve(512,.05,cut=c) for c in (.3,.1)]
    sync = []
    checks = []
    for source in conformal[:3]:
        n = source["n"]
        sample = source["frames"][6]  # t=.6, before global lapse guard
        times = sample["center_proper_times"]
        run = gaussian(n,cut=.1,compare_at=times)
        errors = []
        for i,t in enumerate(times):
            other = min(run["frames"],key=lambda f:abs(f["t"]-t))
            if abs(other["t"]-t)>1e-10:
                raise ArithmeticError("cross-gauge event was not reached")
            errors.append(abs(other["centers"][i]["radius"]-sample["centers"][i]["r"]))
        checks.append(dict(n=n,conformal_t=sample["t"],radius_discrepancies=errors))
        sync.append(run)
    sync += [gaussian(512,cut=c) for c in (.3,.2)]
    sync += [gaussian(256,a,cut=.1) for a in (0.,.025,.1)]
    # Keep time histories for the primary fine run and event diagnostics for controls.
    for i,run in enumerate(conformal):
        if i!=2:
            run["frames"]=[run["frames"][0],run["frames"][-1]]
    for i,run in enumerate(sync):
        if i!=2:
            run["frames"]=[run["frames"][0],run["frames"][-1]]
    short_loop = []
    for n in (128,256,512):
        # Keep the A/L² scale of the longitudinal constraint term fixed.
        run = gaussian(n,.05/(2*np.pi)**2,cut=.1,length=1.)
        run["frames"] = [run["frames"][0],run["frames"][-1]]
        short_loop.append(run)
    return dict(short_loop_followup=short_loop,vacuum=[vacuum_probe(eps) for eps in (.05,.1,.2)],
                conformal=conformal, gaussian=sync, cross_gauge=checks)


if __name__ == "__main__":
    print(json.dumps(probe(),allow_nan=False,indent=2))
