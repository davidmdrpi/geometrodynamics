"""Linear Maxwell-reference closure on a prescribed antipodal history.

This does not replace the strict BAM router or solve the GR mouth boundary.
The free field is the existing round-sphere Maxwell field. Endpoint phase-space
jumps are explicitly specified boundary data, not a derived emission law.
"""
import numpy as np
from scipy.optimize import brentq


def closure_band(angle, k, length=None):
    """Ray-model interval with captured energy still in the throat at contact.

    The upper endpoint needs capture-before-contact ordering. A lower endpoint
    set by throat exit is open; if set by aperture inclusion it is closed.
    Coincident exit/contact at a lower endpoint is excluded here.
    """
    if not 0<angle<np.pi/4 or k<0 or int(k)!=k:
        raise ValueError('Require 0<a<pi/4 and nonnegative integer k')
    length=np.sin(2*angle) if length is None else float(length)
    if length<=0 or not np.isfinite(length):raise ValueError('Positive finite neck length required')
    focus=(2*k+1)*np.pi;distance=np.pi-2*angle
    lower=(np.pi-3*angle)/focus;upper=distance/focus
    def margin(u):return focus+length/np.sqrt(1-u*u)-distance/u
    edge=lower if margin(lower)>0 else brentq(margin,lower,upper,xtol=5e-15)
    return dict(k=int(k),capture_lower=float(lower),upper=float(upper),lower=float(edge),
                lower_closed=bool(margin(lower)>0),length=length,
                exact_image_endpoint=float(upper))


def energy(omega,q,p):
    return float(.5*np.sum(abs(omega*q)**2+abs(p)**2))


def evolve(omega,q,p,time):
    c,s=np.cos(omega*time),np.sin(omega*time)
    return q*c+p*s/omega,p*c-omega*q*s


def endpoint_image(omega,parity,q,p,duration):
    """K U(T) z, K=diag(S,-S); fixed points obey U(T) z=K z."""
    qt,pt=evolve(omega,q,p,duration)
    return parity*qt,-parity*pt


def project_history(omega,parity,q,p,duration):
    """Energy-orthogonal projection onto invariant free histories.

    Replaces the prescribed emission data by z_+=(z+K U(T)z)/2. The change is
    measured, never treated as proof that the original birth flash closes.
    Homogeneous invariant histories remain free data.
    """
    omega,parity,q,p=[np.asarray(x) for x in [omega,parity,q,p]]
    if not (omega.shape==parity.shape==q.shape==p.shape) or np.any(omega<=0) or not np.isin(parity,[-1,1]).all():
        raise ValueError('Matching positive frequencies and involutive parities required')
    if not np.isfinite(duration) or duration<=0 or not all(np.isfinite(x).all() for x in [omega,parity,q,p]):
        raise ValueError('Finite field data and positive duration required')
    qi,pi=endpoint_image(omega,parity,q,p,duration)
    qp,pp=(q+qi)/2,(p+pi)/2
    qcheck,pcheck=endpoint_image(omega,parity,qp,pp,duration)
    e0=energy(omega,q,p)
    return dict(q=qp,p=pp,initial_energy=e0,projected_energy=energy(omega,qp,pp),
        correction_energy=energy(omega,qp-q,pp-p),
        raw_endpoint_defect_energy=energy(omega,qi-q,pi-p),
        projected_endpoint_defect_energy=energy(omega,qcheck-qp,pcheck-pp))


def antipodal_overlap(omega,ell,q,p,time):
    """Best complex-amplitude antipodal template overlap in the energy norm.

    This allows an arbitrary global phase and gain, so it is weaker than
    matching specified real Q/P boundary data. Residual is not an absorption
    probability or an aperture capture fraction.
    """
    z=omega*q+1j*p;weights=abs(z)**2
    if weights.sum()==0:raise ValueError('Nonzero template required')
    overlap=np.sum(weights*np.exp(-1j*omega*time)*(-1.)**ell)/weights.sum()
    return dict(overlap_abs=float(abs(overlap)),phase=float(np.angle(overlap)),
                orthogonal_energy_fraction=float(max(0.,1-abs(overlap)**2)))


def driven_zero_initial(omega,source,time,order=128):
    """Duhamel integral for Q''+omega^2 Q=source(t), zero initial data.

    Source trajectories are prescribed; their work is external in this
    reference problem. This does not infer any BAM recoil.
    """
    nodes,weights=np.polynomial.legendre.leggauss(order)
    times=(nodes+1)*time/2;weights=weights*time/2
    values=np.array([source(t) for t in times])
    phase=(time-times[:,None])*omega
    return (np.sum(weights[:,None]*np.sin(phase)*values/omega,axis=0),
            np.sum(weights[:,None]*np.cos(phase)*values,axis=0))
