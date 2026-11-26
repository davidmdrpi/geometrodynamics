"""Local null-characteristic continuation of the primary charged scalar packet."""
import numpy as np

from scipy.optimize import brentq



def null_map(x,focus,width,strength):
    """Odd, monotone map with pi-periodic derivative and two resolved peaks."""
    b=np.tanh(width)
    def primitive(s):
        k=np.floor((s+np.pi/2)/np.pi)
        angle=s-k*np.pi
        return k*np.pi+np.arctan2(np.sin(angle),b*np.cos(angle))
    def kernel(s):
        den=(b*np.cos(s))**2+np.sin(s)**2
        return b/den,-b*(1-b*b)*np.sin(2*s)/den**2
    left,dl=kernel(x-focus)
    right,dr=kernel(x+focus)
    return (x+strength*(primitive(x-focus)+primitive(x+focus)),
            1+strength*(left+right),strength*(dl+dr))


def propagate(n, charge=.14, lapse_cap=150., remap=False, *, radius=.98,
              amplitude=.05, targets=None, time_limit=3., audit_mass=False, profile_times=()):
    """Second-order characteristic trapezoid. Only evolve fully valid stencils."""
    # Adapt both null labels, preserving the initial slice, normal data and J.
    strength=.18 if remap else 0.
    focus,width=2.292,.015
    period=2*np.pi*(1+2*strength)
    if not np.isfinite(time_limit) or time_limit < 0:
        raise ValueError("time limit must be finite and nonnegative")
    dt=period/n
    xi=np.arange(n)*dt
    def mapped(x):
        return null_map(x,focus,width,strength)
    def original(label):
        if not remap:
            return label
        cycles, local = divmod(label, period)
        return cycles*2*np.pi + brentq(lambda x: mapped(x)[0]-local,0.,2*np.pi)
    x=np.array([original(z) for z in xi]) if remap else xi
    _,jac,djac=mapped(x)
    R=radius
    p_squared=6*(1/R**2-1-charge**2/R**6)
    if p_squared<=0:
        raise ValueError("packet preparation requires positive scalar kinetic energy")
    P=np.sqrt(p_squared)
    f=np.cos(x)*np.exp((np.cos(2*x)-1)/4)
    df=-np.sin(x)*np.exp((np.cos(2*x)-1)/4)-f*np.sin(2*x)/2
    u=-amplitude*f**3
    h=-u/R-27*amplitude*f*df**2/(2*P**2*R)
    pos=np.array([np.full(n,R),-np.log(jac),-3*u/(R*P)])
    vel=np.array([u,h,np.full(n,P)])/jac
    dx=np.array([np.zeros(n),-djac/jac**2,9*amplitude*f*f*df/(R*P*jac)])
    plus,minus=vel+dx,vel-dx
    valid = np.ones(n,dtype=bool)
    center_clock = 0.
    hits = []
    crossings=[]
    cuts=iter((.95,.8,.4,.23,.205,.200,.198,.197,.196) if targets is None else targets)
    next_radius=next(cuts,None)
    frames = []
    clipped = False

    def speed(z,p,m):
        r,s,_ = z
        lapse2 = np.exp(2*s)
        em = lapse2*charge**2/r**6
        sr = 2*lapse2*r-2*(lapse2+p[0]*m[0])/r+2*r*em
        ss = 2*lapse2-1.5*sr/r-p[2]*m[2]/2-4*em
        return np.array([sr,ss,-1.5*(p[0]*m[2]+m[0]*p[2])/r])

    def usable(z,p,m):
        r,s,_ = z
        return (np.isfinite(z).all(axis=0)&np.isfinite(p).all(axis=0)&np.isfinite(m).all(axis=0)
                &(r>.025)&(s<np.log(lapse_cap))&(abs(s)<20)
                &(dt*np.maximum(abs(p[0]),abs(m[0]))<.08*r)
                &(dt*np.maximum(np.max(abs(p[1:]),axis=0),np.max(abs(m[1:]),axis=0))<.08))

    def constraint():
        if not valid[[0,1,-1]].all():
            return None
        sr=speed(pos[:,0],plus[:,0],minus[:,0])[0]
        errors=[]
        for d,sign in ((plus,1),(minus,-1)):
            terms=(sr,sign*(d[0,1]-d[0,-1])/dt,-2*d[1,0]*d[0,0],pos[0,0]*d[2,0]**2/3)
            errors.append(abs(sum(terms))/(1+sum(abs(item) for item in terms)))
        return float(max(errors))

    def mass_slope():
        # At this reflection center w=0. The radial integral remains oriented
        # if the center reverses direction.
        scalar_plus, scalar_minus = plus[2,0], minus[2,0]
        norm = (scalar_plus**2 + scalar_minus**2)*np.exp(-2*pos[1,0])/2
        return -pos[0,0]**3*norm/3

    integrated_mass = R**2*(1+amplitude**2-R**2/2)+charge**2/R**2
    if audit_mass:
        slope = mass_slope()

    def record(step):
        r,s = pos[:2,0]
        u = (plus[0,0]+minus[0,0])*np.exp(-s)/2
        w = (plus[0,0]-minus[0,0])*np.exp(-s)/2
        p=(plus[2,0]+minus[2,0])*np.exp(-s)/2
        f=(plus[2,0]-minus[2,0])*np.exp(-s)/2
        result = dict(t=float(original(step*dt)),tau=float(center_clock),r=float(r),u=float(u),
                    null_constraint=constraint(),
                    y_acc=float(2*(r*r-1+w*w)+2*charge**2/r**4-r*r*(p*p+f*f)/3),
                    scalar_kinetic=float((p*p+f*f)/2),
                    mass=float(r*r*(1+u*u-w*w-r*r/2)+charge**2/r**2),
                    theta_plus=float(3*np.exp(-s)*plus[0,0]/r),
                    theta_minus=float(3*np.exp(-s)*minus[0,0]/r),
                    angular_curvature=float((1+u*u-w*w)/r**2))
        if audit_mass:
            result['integrated_mass'] = float(integrated_mass)
        return result
    profile_steps = {round(time/dt) for time in profile_times}
    profiles = []
    def sample_slice(step):
        # A fraction labels the conformal spatial grid, including after remap.
        samples = []
        for j in (0, n//8, n//4, 3*n//8, n//2):
            sample = dict(grid_fraction=j/n, valid=bool(valid[j]))
            if valid[j]:
                r, s = pos[:2,j]
                u = (plus[0,j]+minus[0,j])*np.exp(-s)/2
                w = (plus[0,j]-minus[0,j])*np.exp(-s)/2
                sample.update(r=float(r), u=float(u), w=float(w),
                              mass=float(r*r*(1+u*u-w*w-r*r/2)+charge**2/r**2))
            samples.append(sample)
        profiles.append(dict(null_grid_time=step*dt, samples=samples))

    if 0 in profile_steps:
        sample_slice(0)
    latest=record(0)
    frames.append(latest)
    last_step=0
    max_steps = int(mapped(time_limit)[0]/dt)
    for step in range(1, max_steps+1):
        # Inactive storage retains its last values. Never read an
        # inactive predecessor into a surviving cell.
        j = np.flatnonzero(valid & np.roll(valid,1) & np.roll(valid,-1))
        old_r = pos[0,0]
        old_s = pos[1,0]
        old_u = plus[0,0]+minus[0,0]
        left,right = (j-1)%n,(j+1)%n
        sl = speed(pos[:,left],plus[:,left],minus[:,left])
        sr = speed(pos[:,right],plus[:,right],minus[:,right])
        pp,mm = plus[:,right]+dt*sr, minus[:,left]+dt*sl
        avg = (pos[:,left]+pos[:,right])/2
        incident = plus[:,left]+minus[:,right]
        zp = avg+dt*incident/2
        ok = usable(zp,pp,mm)
        j,left,right,sl,sr,avg,incident = (j[ok],left[ok],right[ok],sl[:,ok],sr[:,ok],
                                          avg[:,ok],incident[:,ok])
        source = speed(zp[:,ok],pp[:,ok],mm[:,ok])
        pp = plus[:,right]+dt*(sr+source)/2
        mm = minus[:,left]+dt*(sl+source)/2
        zp = avg+dt*(incident+pp+mm)/4
        ok = usable(zp,pp,mm)
        active = j[ok]
        if len(active)==0 or active[0]!=0:
            break
        if len(active)<n:
            clipped = True
        before=None
        if next_radius is not None and old_r>next_radius>=zp[0,ok][0]:
            before=record(step-1)
        valid[:] = False
        valid[active] = True
        pos[:,active],plus[:,active],minus[:,active] = zp[:,ok],pp[:,ok],mm[:,ok]
        if audit_mass:
            new_slope = mass_slope()
            integrated_mass += (pos[0,0]-old_r)*(slope+new_slope)/2
            slope = new_slope
        center_clock += dt*(np.exp(old_s)+np.exp(pos[1,0]))/2
        if step in profile_steps:
            sample_slice(step)
        if old_u < 0 <= plus[0,0]+minus[0,0]:
            hits.append(record(step))
        if before is not None:
            after=record(step)
            while next_radius is not None and next_radius>=pos[0,0]:
                fraction=(old_r-next_radius)/(old_r-pos[0,0])
                event={key:(left+fraction*(after[key]-left)
                            if left is not None and after[key] is not None else None)
                       for key,left in before.items()}
                crossings.append(dict(target=next_radius,**event))
                next_radius=next(cuts,None)
        if step%max(1,n//64)==0:
            frames.append(record(step))
        last_step=step
    latest=record(last_step)
    if frames[-1]['t']!=latest['t']:
        frames.append(latest)
    result = dict(n=n,charge=charge,lapse_cap=lapse_cap,remapped=remap,
                mapping=dict(focus=focus,width=width,strength=strength),
                stop="valid_stencil_limit" if last_step<max_steps else "time_limit",
                last_t=float(original(last_step*dt)),crossings=crossings,
                center_minima=hits,masked_anywhere=clipped,frames=frames)
    if profile_steps:
        result['mass_profiles'] = profiles
    return result
