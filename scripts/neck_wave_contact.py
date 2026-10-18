"""Contact-driven *kinematic* neck study. No independent waist-collapse schedule.

The imposed packet determines a finite-thickness contact event. The post-contact
geometry is a regularized illustration, not a membrane or fluid solution.
"""
from dataclasses import dataclass
import numpy as np

H, END = .6, 9.
X = np.linspace(-H, H, 801)

def smooth(u):
    u=np.clip(u,0,1)
    return u*u*(3-2*u)

@dataclass(frozen=True)
class Parameters:
    waist: float = .16
    amplitude: float = .18
    thickness: float = .025  # radial half-thickness of the display band

class ContactNeck:
    def __init__(self, parameters=Parameters()):
        self.par=parameters
        if not (0 < parameters.thickness < parameters.waist < 1 and 0 <= parameters.amplitude <= .4):
            raise ValueError('Require 0 < thickness < waist < 1 and amplitude in [0,.4]')
        self.tc, self.xc, self.bracket=self.find_contact()

    def rest(self,x):
        return self.par.waist*np.cosh(np.arccosh(1/self.par.waist)*np.asarray(x)/H)

    def wave(self,x,t):
        x=np.asarray(x);c=-.85+.20*t
        taper=np.maximum(0,1-(x/H)**2)**2
        return -self.par.amplitude*np.exp(-((x-c)/.18)**2)*np.cos(14*(x-c))*taper

    def radius(self,x,t):
        return self.rest(x)+self.wave(x,t)

    def gap(self,t):
        return 2*(float(np.min(self.radius(X,t)))-self.par.thickness)

    def find_contact(self):
        lo=0.
        for hi in np.linspace(0,END,361)[1:]:
            if self.gap(hi)<=0:
                bracket=(lo,float(hi))
                for _ in range(36):
                    mid=(lo+hi)/2
                    if self.gap(mid)>0:lo=mid
                    else:hi=mid
                tc=(lo+hi)/2
                xc=float(X[np.argmin(self.radius(X,tc))])
                return float(tc),xc,bracket
            lo=float(hi)
        return None,None,None

    def progress(self,t):
        return 0. if self.tc is None else float(smooth((t-self.tc)/2.8))

    def branches(self,t,n=101):
        if self.tc is None or t < self.tc:
            x=np.linspace(-H,H,n)
            return [(x,self.radius(x,t),self.wave(x,t),self.rest(x))]
        p=self.progress(t);g=.018+.96*p
        # Two rounded closed tips start a finite distance apart. Local blending
        # across the old/new surfaces is a display operation, not physical time integration.
        left=self.xc-(self.xc+H)*g
        right=self.xc+(H-self.xc)*g
        parts=[]
        for rim,tip in [(-H,left),(H,right)]:
            x=tip+(rim-tip)*(1-np.linspace(0,1,n))**2;distance=abs(x-tip)
            L=min(.12,abs(rim-tip)*.65)
            u=np.clip(distance/L,0,1)
            cap=np.sqrt(np.maximum(0,2*u-u*u))
            r=self.radius(x,self.tc)*cap
            # A rounded shoulder on each cap denotes the rolled rim. This is
            # attached to the cap profile, not a detached torus with assigned mass.
            r += .018*np.exp(-((distance-.045)/.022)**2)*cap*(1-p)
            r[0]=1.;r[-1]=0.
            parts.append((x,r,np.zeros_like(x),r))
        return parts

    def phase(self,t):
        if self.tc is None:return 'Passing deformation', 'No contact: the narrow resting neck remains intact'
        dt=t-self.tc
        if dt<-.8:return 'Incoming deformation', 'A fixed narrow neck; the wave supplies the displacement'
        if dt<0:return 'Opposing faces approach', 'The inward crest closes the finite-thickness gap'
        if dt<.65:return 'Finite-width contact and reconnection', 'Local geometric blend; no zero-radius intersection is rendered'
        return 'Rounded rims retract', 'Motion guides illustrate retraction, not a solved fluid velocity'
