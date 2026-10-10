"""Controlled phase-scaling family; same static Berger action as PR #323.

At antipodal ports the dark Gram row vanishes and +/-m have equal frequency.
Combining those degeneracies is exact for a zero initial field. No fitted
propagation kernel or geometric-optics approximation is used in evolution.
"""
from dataclasses import dataclass, asdict
import numpy as np
from .aperture_transfer import cap_coefficients


@dataclass(frozen=True)
class Config:
    squash: float = 1.
    carrier: int = 12
    footprint: float = 4.8  # aperture * carrier; a/lambda = footprint/(2*pi)
    lmax: int = 56
    dt: float = 1/1024
    stop: float = 1.75

    @property
    def aperture(self): return self.footprint/self.carrier
    @property
    def halfwidth(self): return 3/self.carrier
    @property
    def gamma(self): return 8*self.carrier/12


def validate(c):
    if not all(np.isfinite(x) for x in asdict(c).values()):
        raise ValueError('nonfinite configuration')
    if not (.8<=c.squash<=1.2 and c.carrier>=4 and 0<c.aperture<1):
        raise ValueError('invalid geometry/source')
    if c.lmax<1 or int(c.lmax)!=c.lmax or c.dt<=0 or c.stop<=c.halfwidth:
        raise ValueError('invalid discretization')
    if abs((c.stop+c.halfwidth)/c.dt-round((c.stop+c.halfwidth)/c.dt))>1e-8:
        raise ValueError('nonintegral time grid')


def operators(c):
    validate(c)
    coeff=cap_coefficients(c.aperture,c.lmax)
    ls=np.concatenate([np.full(l//2+1,l) for l in range(c.lmax+1)])
    ms=np.concatenate([np.arange(l%2,l+1,2) for l in range(c.lmax+1)])
    multiplicity=np.where(ms==0,1.,2.)
    weight=c.gamma*coeff[ls]**2*multiplicity/(ls+1)
    s=np.sqrt(weight)
    W=np.column_stack([s,(-1.)**ls*s])
    k=np.pi**2*(ls*(ls+2)+(c.squash**-2-1)*ms**2+(4-c.squash**2)/3)
    return k,W,float(max(0,1-coeff@coeff)),ls,ms


def packet(t,c):
    x=np.asarray(t);out=np.zeros_like(x);keep=abs(x)<c.halfwidth
    out[keep]=np.cos(np.pi*x[keep]/(2*c.halfwidth))**4*np.cos(np.pi*c.carrier*x[keep])
    return out


def grid(c):
    validate(c)
    return -c.halfwidth+(np.arange(round((c.stop+c.halfwidth)/c.dt))+.5)*c.dt


def simulate(c):
    k,W,tail,_,_=operators(c);t=grid(c);dt=c.dt
    incoming=np.zeros((len(t),2));incoming[:,0]=packet(t,c)
    inv=1/(1+dt*dt*k/4);B=inv[:,None]*W
    correction=dt/2*B@np.linalg.inv(np.eye(2)+dt/2*(W.T@B))
    q=np.zeros(len(k));v=q.copy();out=np.empty_like(incoming);E=np.zeros(len(t)+1)
    for i,inc in enumerate(incoming):
        rhs=inv*(v-dt/2*k*q+dt*(W@inc))
        vm=rhs-correction@(W.T@rhs)
        q+=dt*vm;v=2*vm-v;out[i]=W.T@vm-inc
        E[i+1]=.5*(v@v+(k*q)@q)
    return dict(config=asdict(c),time=t,incoming=incoming,outgoing=out,
                energy=E,final_q=q,final_v=v,cap_tail=tail)


def source_fourier(omega,c):
    """Real Fourier integral of the compact even packet, not an FFT."""
    h=c.halfwidth
    def box(z): return 2*h*np.sinc(z*h/np.pi)
    def envelope(z):
        return (3/8*box(z)+1/4*(box(z-np.pi/h)+box(z+np.pi/h))
                +1/16*(box(z-2*np.pi/h)+box(z+2*np.pi/h)))
    return (envelope(omega-np.pi*c.carrier)+envelope(omega+np.pi*c.carrier))/2


def phase_summary(c):
    k,W,_,l,m=operators(c)
    round_omega=np.pi*(l+1)
    # Fixed round free-source energy weights; not inferred from captured flux.
    weights=W[:,0]**2*source_fourier(round_omega,c)**2
    weights/=weights.sum()
    delta=np.sqrt(k)-round_omega
    mean=float(weights@delta)
    return dict(proxy=float(np.pi*abs(c.squash**-2-1)*c.carrier/2),
                weighted_mean=float(mean),
                weighted_std=float(np.sqrt(weights@((delta-mean)**2))),
                weighted_free_coherence=float(abs(weights@np.exp(1j*delta))**2))


def diagnose(r):
    c=Config(**r['config']);t=grid(c);n=len(t)
    for name,shape in [('time',(n,)),('incoming',(n,2)),('outgoing',(n,2)),('energy',(n+1,))]:
        if np.shape(r[name])!=shape or not np.all(np.isfinite(r[name])):
            raise ValueError('invalid archive array')
    if not np.array_equal(r['time'],t):raise ValueError('time grid mismatch')
    inc=np.zeros((n,2));inc[:,0]=packet(t,c)
    stored=r['incoming']
    source_error=float(np.max(abs(inc-stored)))
    if source_error>1e-14 or np.any(stored[:,1]!=0) or np.any(stored[abs(t)>=c.halfwidth,0]!=0):
        raise ValueError('source mismatch')
    k,W,tail,_,_=operators(c)
    q=r['final_q'];v=r['final_v'];E=r['energy'];out=r['outgoing']
    if q.shape!=k.shape or v.shape!=k.shape or not np.all(np.isfinite(q+v)):
        raise ValueError('invalid final state')
    # Retain archived source in all physics ledgers after formula validation.
    ein=float(c.dt*np.sum(stored**2))
    ledger=E-E[0]-c.dt*np.r_[0,np.cumsum(np.sum(stored**2-out**2,axis=1))]
    qr=np.zeros(len(k));vr=qr.copy();inv=1/(1+c.dt**2*k/4)
    port_error=0.;energy_reconstruction=0.
    for i in range(n):
        vm=inv*(vr-c.dt*k*qr/2+c.dt/2*(W@(stored[i]-out[i])))
        port_error=max(port_error,float(np.max(abs(W.T@vm-stored[i]-out[i]))))
        qr+=c.dt*vm;vr=2*vm-vr
        energy_reconstruction=max(energy_reconstruction,abs(.5*(vr@vr+(k*qr)@qr)-E[i+1])/ein)
    early=-c.halfwidth+min(1,c.squash)*(np.pi-2*c.aperture)/np.pi
    capture=float(c.dt*np.sum(out[:,1]**2)/ein)
    d=dict(capture=capture,reflected=float(c.dt*np.sum(out[:,0]**2)/ein),
           remaining=float(E[-1]/ein),source_roundoff=source_error,
           energy_error=float(np.max(abs(ledger))/ein),port_error=port_error,
           reconstructed_energy_error=float(energy_reconstruction),
           final_energy_error=float(abs(E[-1]-.5*(v@v+(k*q)@q))/ein),
           state_error=float(np.linalg.norm(qr-q)+np.linalg.norm(vr-v)),
           early_flux=float(c.dt*np.sum(out[t<early,1]**2)/ein),cap_tail=tail,
           phase=phase_summary(c))
    d['valid']=bool(d['energy_error']<1e-9 and d['port_error']<1e-8
                    and d['reconstructed_energy_error']<1e-8 and d['final_energy_error']<1e-9
                    and d['state_error']<1e-8 and d['early_flux']<1e-5 and tail<1e-4
                    and E[0]==0 and np.min(E)>=-1e-12)
    return d
