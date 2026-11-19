"""Gauge-invariant linearized vacuum control on a flat toroidal-fibre comparator."""
import json
import numpy as np
from .clock import modes, L

D=np.diag([1.,-1.,1.])
P=np.diag([-1.,-1.,1.])
def tt(k,n):
    w=np.array([2.,0.,-k])/n
    return np.diag([0.,1.,0.])-np.outer(w,w)

def field(t,x,z,n,m):
    return sum(tt(k,n)*np.cos(2*z+k*x-n*t)/2 for k in (m*np.pi/L,-m*np.pi/L))

def verify():
    errors=[]
    for n,m in modes(9):
        k=m*np.pi/L
        h=tt(k,n)
        errors.extend([abs(np.trace(h)), np.linalg.norm(h @ [k,0,2])/n,
                       abs(n*n-k*k-4)/n**2])
        for t,x,z in np.random.default_rng(42).normal(size=(12,3)):
            ht=field(t,x,z,n,m)
            errors.extend([np.linalg.norm(ht-D @ field(t+np.pi,x+L,z+np.pi,n,m) @ D),
                           np.linalg.norm(ht-P @ field(t,-x,z+np.pi,n,m) @ P)])
    curvature=float(np.linalg.norm(9*field(0,0,0,3,1)/2))
    assert max(errors)<1e-12 and curvature>1
    # Spatial mean of the positive TT quadratic integrability diagnostic on
    # the ordinary flat torus cover, for (n,m)=(3,1) and no homogeneous mode.
    # The overall Einstein-action prefactor is omitted.
    n,m=3,1
    x,z=np.meshgrid(np.arange(32)*2*L/32,np.arange(32)*2*np.pi/32)
    derivatives=np.zeros((3,32,32,3,3))
    for k in (m*np.pi/L,-m*np.pi/L):
        deriv=np.sin(2*z+k*x)[...,None,None]*tt(k,n)[None,None,:,:]/2
        for j,factor in enumerate((n,-k,-2)):
            derivatives[j]+=factor*deriv
    quadratic=float((derivatives**2).sum(axis=(0,3,4)).mean())
    assert abs(quadratic-9)<1e-12
    return dict(flat_cover_quadratic_TT_mean=quadratic,max_TT_dispersion_and_seam_defect=float(max(errors)),
                example_linear_tidal_curvature_frobenius=curvature,
                tested_modes=modes(9))
if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
