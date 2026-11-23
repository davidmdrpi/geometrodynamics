"""Two-port dynamic stiffness on the same handled mesh as network.py."""
import json
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from .network import build

def response(mass, stiff, ports):
    """Uniform traces fixed on both rims. All other regional nodes respond."""
    traces=np.zeros((len(mass),2))
    for j,indices in enumerate(ports):
        traces[indices,j]=1.
    interior=(np.diag(mass)>0) & (traces.sum(axis=1)==0)
    kii=stiff[np.ix_(interior,interior)]
    mii=mass[np.ix_(interior,interior)]
    lam,v=eigh(kii,mii)
    kb=stiff @ traces
    coupling=v.T @ kb[interior]
    direct=traces[:,0] @ kb[:,1]
    def off(omega):
        return direct - np.sum(coupling[:,0]*coupling[:,1]/(lam-omega**2))
    def independent(omega):
        # Reconstruct both full driven steady fields from the original nodal equations.
        z=stiff-omega**2*mass
        field=traces.copy()
        field[interior]=np.linalg.solve(z[np.ix_(interior,interior)],-(z @ traces)[interior])
        return float((traces.T @ z @ field)[0,1])
    return off, independent, np.sqrt(lam)

def run(layers,window):
    mats,_,ports,length,_=build(layers,0.)
    components=[]
    for region in ('tube','exterior'):
        m,_,k=mats[1][region]
        components.append(response(m,k,ports))
    (h,h_direct,ph),(e,e_direct,pe)=components
    poles=np.r_[ph,pe]
    def delta(w):
        return h(w)-e(w)
    lo,hi=window
    assert not np.any((poles>=lo)&(poles<=hi))
    omega=brentq(delta,lo,hi,xtol=1e-14)
    h0,e0=h_direct(omega),e_direct(omega)
    clock=np.pi/omega
    def readout(tau):
        return abs(e0+h0*np.exp(1j*omega*tau))
    ans=dict(tube_layers=layers,omega=omega,half_wave_offset=clock,
             tube_length=length,nearest_regional_dirichlet_pole_gap=float(min(abs(poles-omega))),
             tube_cross_reaction=h0,exterior_cross_reaction=e0,
             bright=readout(0),dark=readout(clock),detuned_clock=readout(1.05*clock))
    assert clock < length and min(abs(h0),abs(e0))>.01
    assert ans['dark'] < 1e-9*ans['bright'] and ans['detuned_clock']>.01*ans['bright']
    return ans

if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='scan':
        mats,_,ports,length,_=build(16,0.)
        hs=[]
        for region in ('tube','exterior'):
            m,_,k=mats[1][region]
            hs.append(response(m,k,ports))
        (h,_,ph),(e,_,pe)=hs
        ps=np.r_[ph,pe]
        for lo in np.arange(3.35,14,.025):
            hi=lo+.025
            if min(abs(ps-(lo+hi)/2))>.05 and (h(lo)-e(lo))*(h(hi)-e(hi))<0:
                print(round(lo,4),round(hi,4),h(lo),e(lo))
    else:
        print(json.dumps([run(layers,(4.575,4.6)) for layers in (8,16,32)],indent=2))
