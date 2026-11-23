"""Right-record probability from the backwards characteristic ending at the readout cut."""
import json

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
from scipy.special import ndtr

from . import timelike_current as kg
from .elliptic_transport import THETAS, split_weights


def sample_grid(values, points, spacing):
    """Cubic interpolation of complex periodic fields before forming a current."""
    n = values.shape[-1]
    index = points/spacing+n/2
    base = np.floor(index).astype(int)
    z = index-base
    weights = [-z*(z-1)*(z-2)/6,(z+1)*(z-1)*(z-2)/2,
               -(z+1)*z*(z-2)/2,(z+1)*z*(z-1)/6]
    return sum(w*np.take(values,(base+i-1)%n,axis=-1) for i,w in enumerate(weights))


def run(c,dt,n,offset=0.):
    x = np.linspace(-56.,56.,n,endpoint=False)
    h = x[1]-x[0]
    kx = 2*np.pi*np.fft.fftfreq(n,h)
    omega = np.sqrt(c**4+c*c*kx*kx)
    g = np.exp(-x*x/(4*kg.SIGMA**2))/(2*np.pi*kg.SIGMA**2)**.25
    u0 = np.stack([g,g])
    pi0 = np.fft.ifft(-1j*omega*np.fft.fft(u0))
    def coefficients(step):
        carrier = np.exp(1j*c*c*step)
        return carrier*np.cos(omega*step),carrier*np.sin(omega*step)
    free = {s:coefficients(s) for s in (dt,-dt,-dt/2)}
    def advance(t,step,state):
        u,pi = state
        half = np.exp(-.5j*kg.K*kg.connection(t+step/2,x,offset)*step)
        uk,pk = np.fft.fft(half*u),np.fft.fft(half*pi)
        co,si = free[step]
        return (half*np.fft.ifft(co*uk+si*pk/omega),
                half*np.fft.ifft(co*pk-si*omega*uk))

    def densities(state,t):
        u,pi = state
        grad = np.fft.ifft(1j*kx*np.fft.fft(u))
        a = kg.connection(t,x,offset)
        dtfield = pi-1j*kg.K*a*u
        e = ((pi.conj()*dtfield).real/c**2
             + (-abs(pi)**2/c**2+abs(grad)**2+c*c*abs(u)**2)/2)
        return np.array([abs(u)**2,-(u.conj()*pi).imag/c**2,e])

    weights = np.array([kg.winding_weights(theta,.02)[[1,-1]] for theta in THETAS]).T
    weights /= weights.sum(axis=0)
    initial = np.einsum("lkn,ka->lan",densities((u0,pi0),0),weights)
    nstep = round(kg.END/dt)
    state = (u0,pi0)
    for j in range(nstep):
        state = advance(j*dt,dt,state)
    final = np.einsum("lkn,ka->lan",densities(state,kg.END),weights)
    final = np.concatenate([final,final[:,:,:1]],axis=-1)
    final_total = simpson(final, dx=h, axis=-1)
    raw_right = simpson(final[:,:,n//2:],dx=h,axis=-1)/final_total
    work_right = np.zeros(len(THETAS))
    xnodes = np.append(x,56.)
    cut = np.zeros((3,len(THETAS)))  # norm, KG number and energy streams
    def velocity(t,state,cut):
        u,pi = state
        grad = np.fft.ifft(1j*kx*np.fft.fft(u))
        amp = np.sqrt(weights[:,None,:])
        u = sample_grid(u,cut,h)*amp
        pi = sample_grid(pi,cut,h)*amp
        grad = sample_grid(grad,cut,h)*amp
        # currents() uses K of shape (2,1); flatten all observation points.
        shape = cut.shape
        rho,j,e,s = kg.currents(u.reshape(2,-1),pi.reshape(2,-1),grad.reshape(2,-1),
                                c,kg.connection(t,cut,offset).ravel())
        number = kg.number_density(u,pi,c)
        if np.min(number[1])<=0 or np.min(e.reshape(shape)[2])<=0:
            raise ValueError("nonpositive transport denominator at cut")
        return np.array([(j/rho).reshape(shape)[0], j.reshape(shape)[1]/number[1],
                         (s/e).reshape(shape)[2]])
    for j in range(nstep,0,-1):
        t = j*dt
        if t <= kg.TAU:
            # Integrate the work on the moving right-hand region, per winding.
            u, pi = state
            charge = (pi.conj()*(-1j*kg.K*u)).real/c**2
            drive_rate = kg.connection(kg.TAU/2,x)*np.pi/kg.TAU*np.sin(2*np.pi*t/kg.TAU)
            power = charge*drive_rate
            periodic_power = np.concatenate([power,power[:,:1]],axis=-1)
            primitive = CubicSpline(xnodes,periodic_power,axis=-1,bc_type="periodic").antiderivative()
            partial = primitive(56.)[:,None]-primitive(cut[2])
            work_right += dt*np.sum(weights*partial,axis=0)
        mid = advance(t,-dt/2,state)
        prev = advance(t,-dt,state)
        v1 = velocity(t,state,cut)
        v2 = velocity(t-dt/2,mid,cut-dt*v1/2)
        v3 = velocity(t-dt/2,mid,cut-dt*v2/2)
        v4 = velocity(t-dt,prev,cut-dt*v3)
        cut -= dt*(v1+2*v2+2*v3+v4)/6
        state = prev
    if np.max(abs(cut))>=36:
        raise ValueError("cut outside central force region")
    trans = np.empty_like(cut)
    for index in np.ndindex(cut.shape):
        density = initial[index]
        spline = CubicSpline(xnodes,np.append(density,density[0])).antiderivative()
        trans[index] = (spline(56.)-spline(cut[index]))/(spline(56.)-spline(-56.))
    initial_energy = h*initial[2].sum(axis=-1)
    budget = (final_total[2]*raw_right[2]-trans[2]*initial_energy-work_right)/initial_energy
    d,_,width,_ = kg.packet(kg.END)
    tail = ndtr(-d/width)
    result = []
    for i,theta in enumerate(THETAS):
        w358 = split_weights(theta,.02)[0]
        result.append(dict(theta=float(theta),born=float(np.cos(theta)**2),
                           nr358=float(tail+(1-2*tail)*w358),
                           nr_same_pair=float(tail+(1-2*tail)*weights[0,i]),
                           right_work_over_initial_energy=work_right[i]/initial_energy[i],
                           right_energy_budget_residual=budget[i],
                           transported=dict(zip(("norm","number","energy"),trans[:,i])),
                           counts2048=dict(zip(("norm","number","energy"),
                                               np.floor(trans[:,i]*2048+.5).astype(int).tolist())),
                           final_fraction=dict(zip(("norm","number","energy"),raw_right[:,i]))))
    born = np.cos(THETAS)**2
    comparison = dict(
        max_shift_from_same_pair=max(abs(r["transported"]["energy"]-r["nr_same_pair"])
                                     for r in result),
        max_energy_error_from_born=float(np.max(abs(trans[2]-born))),
        max_count_error_from_born=float(np.max(abs(np.floor(2048*trans[2]+.5)/2048-born))),
        max_number_error_from_final_mass=float(np.max(abs(trans[1]-raw_right[1]))),
        max_relative_right_work_residual=float(np.max(abs(budget))),
    )
    return dict(c=c,dt=dt,n=n,offset=offset,rows=result,comparison=comparison,
                field_round_trip=float(np.max(abs(state[0]-u0))))


def probe():
    return dict(series=[run(c,.02,1024) for c in (8.,12.,16.,24.)],
                refinement=run(8.,.01,2048), high_c_refinement=run(24.,.01,2048),
                generator_control=run(8.,.02,1024,.1))


if __name__=="__main__":
    print(json.dumps(probe(),indent=2,allow_nan=False))
