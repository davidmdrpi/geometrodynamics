"""Exact electrovac, constrained packet, and conserved homogeneous-flux controls."""
import json

import numpy as np
from scipy.integrate import quad, solve_ivp

from . import nariai_loop as loop
from .charged_characteristics import propagate

ROOTS = np.array([.04,.88,1.08])  # squared inner, outer, cosmological radii


def parameters():
    a,b,c = ROOTS
    return (a*b+b*c+c*a)/2,np.sqrt(a*b*c/2)


def metric(r,mu,q):
    return 1-mu/r**2+q**2/r**4-r**2/2


def radial_force(r,mu,q):
    return r/2-mu/r**3+2*q**2/r**5


def curvature(r,mu,q):
    f = metric(r,mu,q)
    fp = -2*radial_force(r,mu,q)
    fpp = -1-6*mu/r**4+20*q**2/r**6
    return fpp**2+6*(fp/r)**2+12*((1-f)/r**2)**2


def shell(initial_r, initial_u, end, events=None):
    mu,q = parameters()
    def ode(t,y):
        r,u,j,k = y
        return [u,radial_force(r,mu,q),k,(.5+3*mu/r**4-10*q*q/r**6)*j]
    sol = solve_ivp(ode,(0,end),[initial_r,initial_u,1.,0.],method="DOP853",
                    dense_output=True,events=events,rtol=2e-12,atol=1e-13,max_step=.05)
    if not sol.success:
        raise ArithmeticError(sol.message)
    return sol


def exact():
    a,b,c = ROOTS
    mu,q = parameters()
    # Regular time-symmetric chain: s=b+(c-b)sin²x, longitudinal scale
    # A0=sqrt(2s/(s-a)). After a synchronous time T, A=A0*dR(T)/dR0.
    x = np.arange(64)*2*np.pi/64
    s = b+(c-b)*np.sin(x)**2
    t = .25*(np.exp(2*(np.cos(x)-1))-np.exp(2*(-np.cos(x)-1)))/(1-np.exp(-4))
    tx = -.5*np.sin(x)*(np.exp(2*(np.cos(x)-1))+np.exp(2*(-np.cos(x)-1)))/(1-np.exp(-4))
    vals = []
    for r0,T in zip(np.sqrt(s),t):
        vals.append(shell(r0,0.,T).y[:,-1] if abs(T)>1e-14 else [r0,0.,1.,0.])
    r,u,j,_ = np.array(vals).T
    a0 = np.sqrt(2*s/(s-a))
    w0 = (c-b)*np.sin(2*x)/(2*np.sqrt(s)*a0)
    beta = tx/(a0*j)
    if np.min(j)<=0 or np.max(abs(beta))>=1:
        raise ArithmeticError("tilted slice lost regularity")
    boost = 1/np.sqrt(1-beta**2)
    normal = boost*(u+beta*w0)
    spatial = boost*(w0+beta*u)
    masses = r**2*(1+normal**2-spatial**2-r**2/2)+q**2/r**2
    def integrand(eta):
        y = a+(b-a)*np.sin(eta)**2
        return np.sqrt(2*y/(c-y))
    clock = quad(integrand,0,np.pi/2,epsabs=1e-13)[0]
    def turn(t,y):
        return y[1]
    turn.direction=1
    def caustic(t,y):
        return y[2]
    caustic.direction=-1
    trajectories = []
    for energy in (0.,.2):
        orbit = shell(np.sqrt(b),-energy,clock+.5,(turn,caustic))
        time = orbit.t_events[0][0]
        radius = orbit.y_events[0][0,0]
        trajectories.append(dict(killing_energy=energy,time_from_outer=float(time),
                                 radius=float(radius), outward_velocity_after=float(orbit.sol(time+.005)[1]),
                                 force=float(radial_force(radius,mu,q)),
                                 curvature=float(curvature(radius,mu,q)),
                                 jacobi_zero=float(orbit.t_events[1][0]) if len(orbit.t_events[1]) else None,
                                 energy_residual=float(max(abs(orbit.y[1]**2+metric(orbit.y[0],mu,q)-energy**2)))))
    return dict(mu=mu,charge=float(q),horizons=np.sqrt(ROOTS).tolist(),
                max_tilt=float(max(abs(beta))),max_mass_defect=float(max(abs(masses-mu))),
                J_velocity_error=float(max(abs(normal+np.roll(normal,len(x)//2)))),
                initial_center=float(r[0]),paired_velocities=[float(normal[0]),float(normal[len(x)//2])],
                neck_clocks=[float(clock-.25),float(clock+.25)], trajectories=trajectories,
                curvature_at_small_r={str(z):float(curvature(z,mu,q)) for z in (.15,.1,.05)})


def homogeneous(q, scale_floor):
    """Independent homogeneous IVP using the conserved scalar momentum a*r³*p.
    Stop at specified descending longitudinal scale, without extending through zero."""
    r0=.98
    p0=loop.initial(4,0.,charge=q)[5,0]
    flux=p0*r0**3
    def rhs(t,y):
        r,a,u,h=y
        p=flux/(a*r**3)
        acc=r-(1+u*u)/r + q*q/r**5-r*p*p/6
        return np.array([u,a*h,acc,2-p*p-4*q*q/r**6-3*acc/r-h*h])
    def end(t,y):
        return y[1]-scale_floor
    end.terminal=True
    end.direction=-1
    sol=solve_ivp(rhs,(0,10),[r0,1.,0.,0.],events=end,rtol=2e-11,atol=2e-13,method="DOP853")
    if not sol.success or not len(sol.t_events[0]):
        raise ArithmeticError(sol.message)
    r,a,u,h=sol.y[:,-1]
    mass=r*r*(1+u*u-r*r/2)+q*q/r**2
    return dict(scale_floor=scale_floor,time=float(sol.t[-1]),r=float(r),u=float(u),mass=float(mass),
                peak_y_acc=float(max(2*(sol.y[0]**2-1)+2*q*q/sol.y[0]**4
                                       -flux**2/(3*sol.y[1]**2*sol.y[0]**4))))


def packet(n,q,amp=.05):
    cf=loop.evolve(n,amp,cut=.06,charge=q)
    early=cf["frames"][1]
    times=early["center_proper_times"]
    run=loop.gaussian(n,amp,cut=.06,compare_at=times,charge=q)
    errors=[]
    for i,t in enumerate(times):
        match=min(run['frames'],key=lambda z:abs(z['t']-t))
        if abs(match['t']-t)>1e-10:
            raise ArithmeticError("failed to match center clock")
        errors.append(abs(match['centers'][i]['radius']-early['centers'][i]['r']))
    return dict(n=n,amplitude=amp,charge=q,initial=run['frames'][0],
                conformal_stop=cf['stop'],conformal_end=cf['frames'][-1],
                gaussian_stop=run['stop'],gaussian_end=run['frames'][-1],
                center_minima=run.get('center_minima',[]),cross_gauge_error=errors)


def characteristic_probe():
    runs=[]
    for n,remap in ((2048,False),(2048,True),(4096,True),(8192,True),(16384,True),(32768,True)):
        run=propagate(n,remap=remap)
        run["frames"]=[run["frames"][0],run["frames"][-1]]
        runs.append(run)
    return runs


def probe():
    return dict(characteristics=characteristic_probe(),electrovac=exact(),packets=[packet(n,.14) for n in (128,256,512)]
                +[packet(256,q) for q in (0.,.08,.17)],homogeneous_packet=packet(256,.14,0.),
                homogeneous=[homogeneous(.14,z) for z in (.04,.01,.001,.0001)])


if __name__=="__main__":
    print(json.dumps(probe(),allow_nan=False,indent=2))
