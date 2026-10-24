"""Exact GR controls for clock reversal and wave admission.

These are fixed Schwarzschild/ultrastatic geometries, not an evolved BAM pair.
Units c=G=1; stress_components returns 8*pi times orthonormal stress.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import lambertw


def kruskal_radius(time, position, neck=1.):
    """r/b = 1 + W((X^2-T^2)/e), b=2M; exclude curvature singularity."""
    z=np.asarray(position)**2-np.asarray(time)**2
    if neck<=0 or not np.isfinite(z).all() or np.any(z<=-1):
        raise ValueError('Positive neck and regular Kruskal point required')
    return neck*(1+lambertw(z/np.e).real)


def radial_metric(time, position, neck=1.):
    r=kruskal_radius(time,position,neck)
    factor=4*neck**3/r*np.exp(-r/neck)
    return factor*np.diag([-1.,1.])


def radial_boost(rapidity):
    c,s=np.cosh(rapidity),np.sinh(rapidity)
    return np.array([[c,s],[s,c]])


def deck_matrix(time_orientation=1, offset=0.):
    """End-exchanging radial map B(offset) diag(epsilon,-1).

    Add the antipodal S2 map for a free action. epsilon=+1 is the usual
    time-oriented geon; epsilon=-1 reverses physical time. The latter is an
    involution only at offset=0. We return other offsets to audit their failure,
    not to admit them as Z2 quotients.
    """
    if time_orientation not in (-1,1) or not np.isfinite(offset):
        raise ValueError('Time orientation +/-1 and finite offset required')
    return radial_boost(offset)@np.diag([float(time_orientation),-1.])


def killing_time(time, position, neck=1.):
    if neck<=0 or abs(time)>=abs(position):
        raise ValueError('Exterior point required')
    return 2*neck*np.arctanh(time/position)


def stress_components(radius, first, second, lapse=1., lapse_first=0., lapse_second=0.):
    """Einstein tensor for -N(s)^2 dt^2 + ds^2 + r(s)^2 dOmega^2.

    Returns 8*pi*(rho, p_radial, p_tangent) in the static orthonormal frame.
    A static frame at N=0 is excluded, not regularized into a traversable neck.
    """
    r,p,q,n,dn,ddn=np.broadcast_arrays(radius,first,second,lapse,lapse_first,lapse_second)
    if np.any(r<=0) or np.any(n<=0) or not all(np.isfinite(v).all() for v in [r,p,q,n,dn,ddn]):
        raise ValueError('Finite profile with positive radius and lapse required')
    rho=(1-p*p-2*r*q)/(r*r)
    radial=(p*p-1)/(r*r)+2*p*dn/(r*n)
    tangent=q/r+p*dn/(r*n)+ddn/n
    return np.array([rho,radial,tangent])


def schwarzschild_profile(radius, neck=1., vacuum_lapse=True):
    """One exterior of the scalar-flat bridge, parametrized by areal radius."""
    r=np.asarray(radius,float)
    if neck<=0 or np.any(r<=neck):
        raise ValueError('Require r>b>0; static horizon frame is singular')
    p=np.sqrt(1-neck/r);q=neck/(2*r*r)
    if vacuum_lapse:
        return stress_components(r,p,q,p,q,-neck*p/(r**3))
    return stress_components(r,p,q)


def maxwell_potential(tortoise, ell=1):
    """Dimensionless b^2 V_l with x=r*/b and b=2M, for either Maxwell parity.

    Invert y+log(y)=x-1 by Newton in log(y), avoiding exp overflow and loss
    of the tiny near-horizon f. No radiative Maxwell monopole is admitted.
    """
    if int(ell)!=ell or ell<1:
        raise ValueError('Radiative Maxwell modes require integer ell>=1')
    x=np.asarray(tortoise,float)
    z=np.where(x>2,np.log(np.maximum(x-1,1.)),x-1)
    for _ in range(12):
        y=np.exp(z);z-=(y+z-(x-1))/(y+1)
    y=np.exp(z)
    return ell*(ell+1)*y/(1+y)**3


def maxwell_scattering(omega, ell=1, left=-30., right=120., rtol=2e-10):
    """Incoming exterior wave; transmitted flux enters the future horizon.

    omega means omega*b; x means r*/b. Purely ingoing e^(-i omega x)
    at the horizon. Match to both exterior plane waves, retaining the finite
    endpoint approximation. T is NOT transmission to the other exterior.
    """
    if omega<=0 or not np.isfinite([omega,left,right,rtol]).all() or left>=right or rtol<=0:
        raise ValueError('Positive frequency/tolerance and ordered finite domain required')
    maxwell_potential(0.,ell)
    def rhs(x,y):
        return [y[1],(maxwell_potential(x,ell)-omega**2)*y[0]]
    u0=np.exp(-1j*omega*left)
    sol=solve_ivp(rhs,(left,right),np.array([u0,-1j*omega*u0]),method='DOP853',rtol=rtol,atol=rtol*.01)
    if not sol.success:raise RuntimeError(sol.message)
    u,d=sol.y[:,-1]
    incoming=.5*(u+1j*d/omega)*np.exp(1j*omega*right)
    outgoing=.5*(u-1j*d/omega)*np.exp(-1j*omega*right)
    reflection=abs(outgoing/incoming)**2
    horizon=1/abs(incoming)**2
    current=np.imag(np.conj(sol.y[0])*sol.y[1])
    return dict(omega_b=float(omega),ell=int(ell),left=float(left),right=float(right),
                reflection=float(reflection),horizon_flux=float(horizon),
                flux_error=float(abs(reflection+horizon-1)),
                current_relative_error=float(np.max(abs(current+omega))/omega),
                evaluations=int(sol.nfev),rtol=rtol)
