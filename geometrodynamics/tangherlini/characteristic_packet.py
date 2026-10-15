"""2D test scalar and O(A^2) spherical Einstein response on Tangherlini.

Two bulk patches are evolved independently. Matching their scalar traces does
not imply that their perturbed induced metrics match. No nonlinear metric or
surface matter evolution is hidden in this module.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson
from scipy.interpolate import CubicSpline


def radius(U, V, b=1.):
    z = -np.asarray(U)*np.asarray(V)/b**2
    if b <= 0 or np.any(z < 0):
        raise ValueError('require b>0 and exterior UV<=0')
    x = 1 + .5*np.log1p(z)
    for _ in range(24):
        e = np.exp(2*x)
        h = e*(x-1)/(x+1)
        hp = e*(2*(x-1)/(x+1)+2/(x+1)**2)
        step = (h-z)/hp
        x = np.maximum(1., x-step)
        if np.max(np.abs(step)) < 2e-14:
            break
    return b*x


def geometry(U, V, ell=0, b=1.):
    r = radius(U,V,b)
    F = ((r+b)/r)**2*np.exp(-2*r/b)
    f = 1-b*b/r**2
    Q = F/4*(ell*(ell+2)/r**2+3*f/(4*r**2)+3*b*b/r**4)
    return r,F,Q


def compact_packet(V, b=1., amplitude=.01):
    V = np.asarray(V)
    y = (V-2*b)/b
    out = np.zeros_like(V,dtype=float)
    mask = np.abs(y)<1
    out[mask] = amplitude*np.exp(-.5*(4*y[mask])**2+1-1/(1-y[mask]**2))
    return out


def diamond(U, V, left, bottom, Q):
    """Implicit cell-midpoint diamond, second order on smooth nonuniform grids.

    U,V increase; left is data at U[0], bottom at V[0]. Q is at cell centers.
    Independent solutions may occupy an optional last axis.
    """
    U,V = np.asarray(U),np.asarray(V)
    left,bottom = np.asarray(left),np.asarray(bottom)
    if np.any(np.diff(U)<=0) or np.any(np.diff(V)<=0):
        raise ValueError('coordinates must increase')
    if left.shape[0]!=len(V) or bottom.shape[0]!=len(U) or left.shape[1:]!=bottom.shape[1:]:
        raise ValueError('boundary shape mismatch')
    if not np.allclose(left[0],bottom[0],rtol=0,atol=1e-14):
        raise ValueError('incompatible characteristic corner')
    field = np.empty((len(U),len(V))+left.shape[1:])
    field[0] = left
    field[:,0] = bottom
    weight = np.diff(U)[:,None]*np.diff(V)[None,:]*Q/4
    if left.ndim==2:
        weight=weight[...,None]
    nu,nv=len(U)-1,len(V)-1
    for diagonal in range(nu+nv-1):
        i=np.arange(max(0,diagonal-nv+1),min(nu-1,diagonal)+1)
        j=diagonal-i
        a=weight[i,j]
        field[i+1,j+1]=((1-a)*(field[i+1,j]+field[i,j+1])-(1+a)*field[i,j])/(1+a)
    return field


def symplectic_balance(U,V,fields):
    p,q=fields[...,0],fields[...,1]
    pu,pv=np.gradient(p,U,V,edge_order=2)
    qu,qv=np.gradient(q,U,V,edge_order=2)
    ju=p*qu-q*pu
    jv=p*qv-q*pv
    edges = dict(left=float(simpson(jv[0],x=V)),right=float(simpson(jv[-1],x=V)),
                 bottom=float(simpson(ju[:,0],x=U)),top=float(simpson(ju[:,-1],x=U)))
    denom=max(abs(edges['left'])+abs(edges['bottom']),1e-30)
    edges['relative_residual']=abs(edges['right']-edges['left']+edges['top']-edges['bottom'])/denom
    return edges


def mass_response(U,V,field,b=1.,kappa5=1.):
    """Integrate the spherical Einstein mass one-form along two paths.

    m=r^2(1+4 r_U r_V/F); M=3 Omega3 m/(2 kappa5).
    Fixed-background scalar supplies the full leading O(A^2) source.
    Valid only for the spherical mode. Initial corner delta_m is set to zero.
    """
    r,_,_=geometry(U[:,None],V[None,:],b=b)
    phi=field/r**1.5
    pu,pv=np.gradient(phi,U,V,edge_order=2)
    mu=2*kappa5/(3*b)*r**3*U[:,None]*pu**2
    mv=2*kappa5/(3*b)*r**3*V[None,:]*pv**2
    iu=cumulative_trapezoid(mu,x=U,axis=0,initial=0)
    iv=cumulative_trapezoid(mv,x=V,axis=1,initial=0)
    path_uv=iu[:,0,None]+iv
    path_vu=iv[0,None,:]+iu
    scale=max(float(np.max(np.abs(path_uv))),float(np.max(np.abs(path_vu))),1e-30)
    return dict(delta_m=path_uv,other_path=path_vu,
                relative_path_discrepancy=float(np.max(np.abs(path_uv-path_vu))/scale),
                corner_delta_m=float(path_uv[-1,-1]),
                opposite_corner_path=float(path_vu[-1,-1]))


def seam_response(V,trace,b=1.,kappa5=1.):
    """Two independent horizon Raychaudhuri constraints at O(A^2).

    delta r_in''=-S; delta r_out''+2/V delta r_out'=-S.
    Equal initial radius, zero initial expansion. The mismatch cannot be removed
    by an O(A^2) lapse/map correction because the background expansion vanishes.
    """
    phip=np.gradient(trace/b**1.5,V,edge_order=2)
    S=kappa5*b/3*phip**2
    I=[cumulative_trapezoid(V**k*S,x=V,initial=0) for k in range(3)]
    rin=-V*I[0]+I[1]
    rout=-I[1]+I[2]/V
    dm=2*kappa5*b*b/3*cumulative_trapezoid(V*phip**2,x=V,initial=0)
    return dict(incoming_delta_r=rin,outgoing_delta_r=rout,
                induced_radius_mismatch=rin-rout,
                incoming_delta_m=dm,outgoing_delta_m=-dm,
                mass_parameter_jump=-2*dm)


def solve_pair(N=128,ell=0,amplitude=.01,boost=1.,b=1.,c=None):
    if not isinstance(N,int) or N<8 or boost<=0 or b<=0:
        raise ValueError('require integer N>=8 and positive boost,b')
    c=b*b if c is None else c
    if c<=0 or ell<0 or int(ell)!=ell:
        raise ValueError('require c>0 and nonnegative integer ell')
    # Domain is a fixed physical set of points; boosts change its labels.
    U=np.linspace(-8*b,0,N+1)/boost
    V=np.linspace(.5*b,4*b,2*N+1)*boost
    incoming=compact_packet(V/boost,b,amplitude)
    left=np.column_stack((incoming,(V/boost-2*b)/b*incoming))
    bottom=np.zeros((len(U),2))
    _,_,Q=geometry((U[:-1,None]+U[1:,None])/2,(V[None,:-1]+V[None,1:])/2,ell,b)
    field=diamond(U,V,left,bottom,Q)
    out_U=-c/V
    out_V=np.linspace(0,80*b,2*N+1)*boost
    _,_,Qo=geometry((out_U[:-1,None]+out_U[1:,None])/2,
                    (out_V[None,:-1]+out_V[None,1:])/2,ell,b)
    outgoing=diamond(out_U,out_V,np.zeros((len(out_V),2)),(-1)**ell*field[-1],Qo)
    # r_obs=2b worldline: UV=-(b^2/3) exp(4).
    obs_V=-(b*b/3)*np.exp(4)/out_U
    valid=obs_V<=out_V[-1]
    obs=np.array([CubicSpline(out_V,outgoing[i,:,0])(vv) for i,vv in enumerate(obs_V[valid])])
    t=b/2*np.log(obs_V[valid]/(-out_U[valid]))-b*np.log(boost)
    result=dict(U=U,V=V,field=field,out_U=out_U,out_V=out_V,outgoing=outgoing,
                observer_t=t,observer_proper_t=np.sqrt(.75)*t,observer_u=obs,
                incoming_symplectic=symplectic_balance(U,V,field),
                outgoing_symplectic=symplectic_balance(out_U,out_V,outgoing))
    if ell==0:
        result['incoming_mass']=mass_response(U,V,field[...,0],b)
        result['outgoing_mass']=mass_response(out_U,out_V,outgoing[...,0],b)
        result['seam']=seam_response(V/boost,field[-1,:,0],b)
        rr,_,_=geometry(out_U[:,None],out_V[None,:],b=b)
        phi=outgoing[...,0]/rr**1.5
        pu,pv=np.gradient(phi,out_U,out_V,edge_order=2)
        pu_obs=np.array([CubicSpline(out_V,pu[i])(vv) for i,vv in enumerate(obs_V[valid])])
        pv_obs=np.array([CubicSpline(out_V,pv[i])(vv) for i,vv in enumerate(obs_V[valid])])
        dt=(-out_U[valid]*pu_obs+obs_V[valid]*pv_obs)/b
        dr=(out_U[valid]*pu_obs+obs_V[valid]*pv_obs)/b
        omega3=2*np.pi**2
        result['observer_outward_killing_flux']=-omega3*(2*b)**3*dt*dr
        dm=result['outgoing_mass']['delta_m']
        result['observer_delta_M']=3*omega3/2*np.array([
            CubicSpline(out_V,dm[i])(vv) for i,vv in enumerate(obs_V[valid])])
    return result
