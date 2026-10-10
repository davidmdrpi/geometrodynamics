"""Finite distributed scalar ports on an ultrastatic Berger S3.

Exact SU(2) spectral blocks, truncated at lmax. This is a prescribed
3-spatial-dimensional geometry and aperture action, not an excised GR mouth.
The port-observable subspace is compressed exactly inside each eigenblock;
no fitted travel-time kernel is used. All quantities use R=1/pi, c=1.
"""
from dataclasses import dataclass, asdict
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh_tridiagonal
from numpy.polynomial.legendre import leggauss


@dataclass(frozen=True)
class Config:
    squash: float = 1.0
    aperture: float = .4
    carrier: float = 12.0  # angular frequency / pi
    offset: float = 0.0  # receiver displacement from -identity, round radians
    axis: str = 'horizontal'
    lmax: int = 40
    dt: float = 1/512
    stop: float = 1.75
    connected: bool = True
    source: int = 0
    xi: float = 1/6
    gamma: float = 8.0


def validate(c):
    numbers=(c.squash,c.aperture,c.carrier,c.offset,c.dt,c.stop,c.xi,c.gamma)
    if not all(np.isfinite(numbers)) or not (.8<=c.squash<=1.2 and .1<=c.aperture<1):
        raise ValueError('invalid geometry')
    if c.dt<=0 or c.stop<=.25 or c.carrier<=0 or c.gamma<=0 or c.xi<0:
        raise ValueError('invalid dynamics')
    if c.lmax<1 or int(c.lmax)!=c.lmax or c.source not in (0,1) or c.axis not in ('horizontal','fiber'):
        raise ValueError('invalid discretization or port')
    if not (0<=c.offset<=1) or abs((c.stop+.25)/c.dt-round((c.stop+.25)/c.dt))>1e-9:
        raise ValueError('invalid offset or time grid')


@lru_cache(maxsize=32)
def cap_coefficients(aperture, lmax, quadrature=256):
    """Coefficients of an L2-normalized compact C3 radial aperture.

    Haar radial measure is (2/pi) sin(theta)^2 dtheta. Characters
    sin((l+1)theta)/sin(theta) are orthonormal. Do not renormalize the
    truncated vector: its missing norm is an independently checked tail.
    """
    z,w=leggauss(quadrature);theta=(z+1)*aperture/2;w=w*aperture/2
    shape=(1-(theta/aperture)**2)**4
    norm=np.sqrt(2/np.pi*np.sum(w*np.sin(theta)**2*shape**2))
    l=np.arange(lmax+1)
    a=2/np.pi*(np.sin((l[:,None]+1)*theta)* (w*np.sin(theta)*shape)).sum(axis=1)/norm
    return a


@lru_cache(maxsize=256)
def rotation_diagonal(l, offset, axis):
    m=np.arange(-l,l+1,2,dtype=float)
    if offset==0:return np.ones(l+1)
    if axis=='fiber':return np.cos(m*offset)
    if l==0:return np.ones(1)
    j=l/2;z=m/2
    off=.5*np.sqrt((j-z[:-1])*(j+z[:-1]+1))
    eig,U=eigh_tridiagonal(np.zeros(l+1),off)
    return (U**2)@np.cos(2*offset*eig)


def operators(c):
    validate(c);a=cap_coefficients(c.aperture,c.lmax)
    freq=[];rows=[]
    for l in range(c.lmax+1):
        m=np.arange(-l,l+1,2,dtype=float)
        k=np.pi**2*(l*(l+2)+(c.squash**-2-1)*m*m+2*c.xi*(4-c.squash**2))
        d=(-1)**l*rotation_diagonal(l,c.offset,c.axis)
        weight=c.gamma*a[l]**2/(l+1)
        for sign in (1.,-1.):
            s=np.sqrt(np.maximum(0,weight*(1+sign*d)/2))
            rows.extend(np.column_stack((s,sign*s)));freq.extend(k)
    W=np.asarray(rows);k=np.asarray(freq)
    if not c.connected:W[:,1]=0
    return k,W,float(max(0,1-np.dot(a,a)))


def packet(t,carrier):
    x=np.asarray(t);f=np.zeros_like(x);keep=abs(x)<.25
    f[keep]=np.cos(2*np.pi*x[keep])**4*np.cos(np.pi*carrier*x[keep])
    return f


def simulate(c):
    """Implicit midpoint with diagonal-plus-rank-two inversion.

    q''+Kq+W W^T q'=2W a; b=W^T q'-a.
    Exactly: E[n+1]-E[n]=dt*(|a|^2-|b|^2) at midpoint.
    """
    k,W,tail=operators(c);dt=c.dt
    n=round((c.stop+.25)/dt);t=-.25+(np.arange(n)+.5)*dt
    inc=np.zeros((n,2));inc[:,c.source]=packet(t,c.carrier)
    inv=1/(1+dt*dt*k/4);B=inv[:,None]*W
    C=np.eye(2)+dt/2*(W.T@B)
    correction=dt/2*B@np.linalg.inv(C)
    q=np.zeros(len(k));v=q.copy();out=np.empty_like(inc);energy=np.zeros(n+1)
    for i in range(n):
        rhs=inv*(v-dt/2*k*q+dt*(W@inc[i]))
        vm=rhs-correction@(W.T@rhs)
        q+=dt*vm;v=2*vm-v
        out[i]=W.T@vm-inc[i]
        energy[i+1]=.5*(np.dot(v,v)+np.dot(k*q,q))
    return dict(config=asdict(c),time=t,incoming=inc,outgoing=out,energy=energy,
                final_q=q,final_v=v,cap_tail=tail)


def scattering(c,omega):
    """S(omega) with exp(-i omega t) convention; frequencies avoid poles."""
    k,W,_=operators(c);out=[]
    for w in np.atleast_1d(omega):
        if np.min(abs(k-w*w))<1e-9:raise ValueError('frequency at a free pole')
        G=W.T@(W/(k-w*w)[:,None]);I=np.eye(2)
        out.append(-np.linalg.solve(I-1j*w*G,I+1j*w*G))
    return np.asarray(out)


def shifted_return(record,clock_offset=1.5,handle_time=.125):
    """First-transit feed-forward return, not a closed feedback solution.

    Relabel every receiver sample once; no periodic wrapping, interpolation,
    gain, or discarded shifted samples. Its flux is the removed B-lead flux.
    """
    if handle_time<=0 or not np.isfinite(clock_offset+handle_time):raise ValueError('invalid clock')
    return record['time']+handle_time-clock_offset,-record['outgoing'][:,1]


def diagnose(r):
    c=Config(**r['config']);validate(c);n=round((c.stop+.25)/c.dt)
    t=-.25+(np.arange(n)+.5)*c.dt
    for key,shape in [('time',(n,)),('incoming',(n,2)),('outgoing',(n,2)),('energy',(n+1,))]:
        if np.shape(r[key])!=shape or not np.all(np.isfinite(r[key])):raise ValueError('invalid archive array')
    if not np.array_equal(t,r['time']):raise ValueError('time grid mismatch')
    inc=np.zeros((n,2));inc[:,c.source]=packet(t,c.carrier)
    if not np.array_equal(inc,r['incoming']):raise ValueError('source mismatch')
    k,W,tail=operators(c)
    q=np.asarray(r['final_q']);v=np.asarray(r['final_v'])
    if q.shape!=k.shape or v.shape!=k.shape or not np.all(np.isfinite(q+v)):raise ValueError('invalid final state')
    ein=c.dt*np.sum(inc**2);b=r['outgoing'];E=r['energy']
    ledger=E-E[0]-c.dt*np.r_[0,np.cumsum(np.sum(inc**2-b*b,axis=1))]
    receiver=b[:,1]**2;capture=c.dt*np.sum(receiver)/ein
    mean=float(np.sum(t*receiver)/np.sum(receiver)) if capture>1e-20 else None
    spread=float(np.sqrt(np.sum((t-mean)**2*receiver)/np.sum(receiver))) if mean is not None else None
    tr,ret=shifted_return(r);t0,ret0=shifted_return(r,0)
    # Earliest possible travel between compact round-coordinate caps:
    # g_b >= min(1,b)^2 g_round. Curvature coupling changes no characteristics.
    earliest=-.25+min(1,c.squash)*(np.pi-c.offset-2*c.aperture)/np.pi
    final_error=abs(E[-1]-.5*(v@v+(k*q)@q))/ein
    # Reconstruct the state from the recorded port forces a-b, using the
    # free diagonal midpoint equation rather than the coupled solver.
    qr=np.zeros(len(k));vr=qr.copy();inv=1/(1+c.dt*c.dt*k/4)
    field_error=0.;state_energy_error=0.
    for i in range(n):
        vm=inv*(vr-c.dt*k*qr/2+c.dt/2*(W@(inc[i]-b[i])))
        field_error=max(field_error,float(np.max(abs(W.T@vm-inc[i]-b[i]))))
        qr+=c.dt*vm;vr=2*vm-vr
        state_energy_error=max(state_energy_error,abs(.5*(vr@vr+(k*qr)@qr)-E[i+1])/ein)
    state_error=float(np.linalg.norm(qr-q)+np.linalg.norm(vr-v))
    out=dict(capture_fraction=float(capture),arrival_mean=mean,arrival_spread=spread,
             early_leak_fraction=float(c.dt*np.sum(receiver[t<earliest])/ein),
             energy_error=float(np.max(abs(ledger))/ein),final_energy_error=float(final_error),
             field_error=field_error,state_energy_error=state_energy_error,state_error=state_error,
             cap_tail=tail,remaining_energy=float(E[-1]/ein),
             advanced_return_fraction=float(c.dt*np.sum(ret[tr<-.25]**2)/ein),
             causal_return_fraction=float(c.dt*np.sum(ret0[t0<-.25]**2)/ein))
    out['valid']=bool(out['energy_error']<1e-9 and final_error<1e-9 and tail<1e-4
                      and field_error<1e-8 and state_energy_error<1e-8 and state_error<1e-8
                      and out['early_leak_fraction']<1e-5 and E[0]==0 and np.min(E)>=-1e-12)
    return out
