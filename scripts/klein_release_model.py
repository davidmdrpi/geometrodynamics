"""Illustrative sphere -> intrinsic throat -> release, with one sphere mode basis.

The sphere field is solved in conformal S2 modes. Throat readout is prescribed
one-way delayed transmission; smooth withdrawal is an added irreversible law.
The release reservoir is assigned. Only its work budget, including interference,
is balanced here; this is not a coupled GR/fluid calculation.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from neck_wave_contact import ContactNeck, Parameters


def smooth(x):
    x=np.clip(x,0,1)
    return x*x*x*(10+x*(-15+6*x))


def energy(a,p,omega):
    return .5*np.sum(p*p+omega*omega*a*a,axis=-1)


class KleinRelease:
    def __init__(self,drive=1.,gain=.18,detach=True,budget=.04):
        self.drive=drive;self.gain=gain;self.detach=detach;self.budget=budget
        self.l=np.arange(25);self.omega=self.l+.5
        self.norm=np.sqrt((2*self.l+1)/(4*np.pi));self.focus=1.2
        f=self.l*(self.l+1)*np.exp(-.24**2*self.l*(self.l+1)/2)*self.norm
        self.f=f*np.sqrt(.08/energy(f,np.zeros_like(f),self.omega))
        self.peak=float(self.norm@self.f)
        self.neck=ContactNeck(Parameters(amplitude=gain))
        self.length=2*np.sin(.3)*np.cos(.3)
        self.tc=self.find_contact();self.withdrawal=.65;self.pulse_duration=.3
        self.h=self.norm*(1-(-1.)**self.l)*np.exp(-.3**2*self.l*(self.l+1)/2)
        self.z,self.w=np.polynomial.legendre.leggauss(64)
        self.scale=0.;self.self_energy=0.;self.cross_energy=0.
        if self.tc is not None and detach:
            a,p=self.incident(self.tc+self.pulse_duration)
            q,v=self.unit_release(self.pulse_duration)
            e=float(energy(q,v,self.omega));cross=float(p@v+(self.omega*a)@(self.omega*q))
            root=np.sqrt(cross*cross+4*e*budget)
            self.scale=2*budget/(root+cross) if cross>=0 else (root-cross)/(2*e)
            self.self_energy=self.scale**2*e;self.cross_energy=self.scale*cross
        self.end=(self.tc if self.tc is not None else self.focus)+self.pulse_duration+2*np.pi+.25

    def basis(self,mu):
        return np.polynomial.legendre.legvander(mu,24)*self.norm

    def incident(self,t):
        angle=self.omega*(t-self.focus)
        return self.drive*self.f*np.cos(angle),-self.drive*self.omega*self.f*np.sin(angle)

    def mouth(self,t):
        # Vectorized normalized trace at mouth A; source amplitude is NOT
        # renormalized when drive changes, so the wave-off null is meaningful.
        t=np.asarray(t)
        return self.drive*np.sum(self.norm*self.f*np.cos((t[...,None]-self.focus)*self.omega),axis=-1)/self.peak

    def transmitted(self,x,t):
        delay=(np.asarray(x)+.6)/1.2*self.length
        # A launch window provides C2 retarded startup; the target focus occurs
        # after it ends. This one-way transmission is not a neck wave PDE solve.
        source_time=t-delay
        return self.mouth(source_time)*smooth(source_time/.15)

    def radius(self,x,t):
        x=np.asarray(x)
        return self.neck.rest(x)-self.gain*self.transmitted(x,t)*(1-(x/.6)**2)**2

    def minimum(self,t):
        x=np.linspace(-.6,.6,401);r=self.radius(x,t);i=int(np.argmin(r))
        if i in (0,len(x)-1):return float(r[i])
        return float(minimize_scalar(lambda x:self.radius(x,t),bounds=(x[i-1],x[i+1]),method='bounded',options={'xatol':1e-13}).fun)

    def find_contact(self):
        f=lambda t:self.minimum(t)-self.neck.level_radius
        times=np.linspace(0,3.,301);prev=f(0)
        for lo,hi in zip(times[:-1],times[1:]):
            now=f(hi)
            if now<=0<prev:return float(brentq(f,lo,hi,xtol=1e-12))
            prev=now
        return None

    def progress(self,t):
        return 0. if self.tc is None or not self.detach else float(smooth((t-self.tc)/self.withdrawal))

    def throat_clock(self,t):
        if self.tc is None or not self.detach or t<=self.tc:return t
        u=min((t-self.tc)/.15,1.)
        # Integral of 1-smooth(u): C2 stopping of the readout, not a hard freeze.
        return self.tc+.15*(u-2.5*u**4+3*u**5-u**6)

    def density(self,x,rho,t):
        r=self.radius(x,self.throat_clock(t));eps=.025
        a=np.exp(-((np.asarray(rho)-r)/eps)**2);b=np.exp(-((np.asarray(rho)+r)/eps)**2)
        q=self.progress(t)
        return (a+b-a*b)*(1-smooth((1.25*q-(np.asarray(x)/.6)**2)/.25))

    def unit_release(self,dt):
        if dt<=0:return np.zeros_like(self.h),np.zeros_like(self.h)
        T=self.pulse_duration;end=min(dt,T);s=(self.z+1)*end/2;w=self.w*end/2
        pulse=30*(s/T)**2*(1-s/T)**2/T
        phase=self.omega[:,None]*(dt-s)
        return self.h*np.sum(w*pulse*np.sin(phase)/self.omega[:,None],axis=1),self.h*np.sum(w*pulse*np.cos(phase),axis=1)

    def state(self,t,emitted_only=False):
        a,p=self.incident(t)
        if emitted_only:a,p=0*a,0*p
        if self.tc is not None and self.detach:
            q,v=self.unit_release(t-self.tc);a=a+self.scale*q;p=p+self.scale*v
        return a,p

    def ledger(self,t):
        a,p=self.state(t);u,v=self.incident(t)
        work=float(energy(a,p,self.omega)-energy(u,v,self.omega))
        return dict(work=work,reservoir=self.budget-work)

    def cap_energy(self,t):
        a,p=self.state(t,emitted_only=True)
        z,w=np.polynomial.legendre.leggauss(80);lo=np.cos(.3)
        mu=np.r_[lo+(z+1)*(1-lo)/2,-lo-(z+1)*(1-lo)/2]
        weights=np.tile(w*(1-lo)/2,2)
        u=self.basis(mu)@a;ut=self.basis(mu)@p
        derivative=np.polynomial.legendre.legval(mu,np.polynomial.legendre.legder(self.norm*a))
        return float(np.pi*np.sum(weights*(ut*ut+(1-mu*mu)*derivative*derivative+.25*u*u)))

    def report(self):
        ts=np.linspace(0,self.end,501);ledger=[self.ledger(t) for t in ts]
        center=(self.tc or self.focus)+self.pulse_duration/2
        mid=self.cap_energy(center+np.pi/2);focus=self.cap_energy(center+np.pi)
        return dict(status='ILLUSTRATIVE_ONE_WAY_CONTACT_AND_RELEASE',contact_time=self.tc,focus_time=self.focus,
                    gain=self.gain,drive=self.drive,detach=self.detach,intrinsic_neck_length=self.length,
                    withdrawal_duration=self.withdrawal,pulse_duration=self.pulse_duration,
                    assigned_release_budget=self.budget,release_self_energy=self.self_energy,release_interference=self.cross_energy,
                    final_added_field_energy=ledger[-1]['work'],minimum_reservoir=min(v['reservoir'] for v in ledger),
                    maximum_reservoir=max(v['reservoir'] for v in ledger),minimum_added_energy=min(v['work'] for v in ledger),
                    emitted_cap_energy_midflight=mid,emitted_cap_energy_antipode=focus,cap_refocus_ratio=focus/mid if mid>0 else None,
                    maximum_budget_residual=max(abs(v['work']+v['reservoir']-self.budget) for v in ledger),
                    mode_count=25,sphere_equation='u_tt - Delta_S2 u + u/4 = F; omega_l=l+1/2',
                    limitations=['no solved neck stress or geometry feedback','readout gain and irreversible withdrawal are imposed','assigned reservoir is not derived throat energy','smooth scalar density is not a constructed regular spacetime metric'])
