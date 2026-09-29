"""Retrospectively repaired kinematic R2 illustration (not a dynamical falsifier): global degree of phi/|phi| on closed S^3 for phi = R(eta) x + eps*dphi(x).
The random-root search and quadrature are illustrative, not exhaustive or
integer-degree certificates. The Hopf control has an analytic nonzero norm.
Original pre-freeze code is preserved in git at 2e984ac."""
import numpy as np
from scipy.optimize import root
rng=np.random.default_rng(7)

def grid(n=72):
    a=(np.arange(n)+.5)*(np.pi/2)/n; b=(np.arange(2*n)+.5)*2*np.pi/(2*n)
    A,B1,B2=np.meshgrid(a,b,b,indexing='ij')
    X=np.stack([np.cos(A)*np.cos(B1),np.cos(A)*np.sin(B1),np.sin(A)*np.cos(B2),np.sin(A)*np.sin(B2)],-1)
    return X,(a[1]-a[0],b[1]-b[0])

def degree(F,h):
    norm=np.linalg.norm(F,axis=-1,keepdims=True)
    if not np.isfinite(F).all() or np.any(norm == 0):
        raise ValueError("degree undefined at a sampled zero or nonfinite field")
    Fh=F/norm
    d=[np.gradient(Fh,h[0],axis=0),np.gradient(Fh,h[1],axis=1,edge_order=2),np.gradient(Fh,h[1],axis=2,edge_order=2)]
    # periodic in xi1, xi2: use roll-based central differences
    d[1]=(np.roll(Fh,-1,1)-np.roll(Fh,1,1))/(2*h[1]); d[2]=(np.roll(Fh,-1,2)-np.roll(Fh,1,2))/(2*h[1])
    J=np.linalg.det(np.stack([Fh,d[0],d[1],d[2]],-1))
    return J.sum()*h[0]*h[1]**2/(2*np.pi**2)

# random odd perturbation: cubic odd polynomial map R^4->R^4 plus linear part
Cl=rng.normal(size=(4,4)); Cc=rng.normal(size=(4,4,4,4))*.4
def dphi(X): return X@Cl.T+np.einsum('aijk,...i,...j,...k->...a',Cc,X,X,X)
J=np.array([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]],float)
def hopf(X): return X@J.T

Rp=np.sqrt(3)  # R ~ -sqrt(3)*tau near eta=pi/4
def tangential_zeros(f,seeds=4000):
    """zeros of f_perp on S^3: solve f(x)-(x.f)x=0 with |x|=1, dedupe; return x and a=x.f"""
    S=rng.normal(size=(seeds,4)); S/=np.linalg.norm(S,axis=1,keepdims=True); out=[]
    for s in S:
        def G(z):
            x=z[:4];lam=z[4];fx=f(x[None])[0];return np.r_[fx-lam*x,x@x-1]
        r=root(G,np.r_[s,s@f(s[None])[0]],tol=1e-13)
        if r.success and np.linalg.norm(G(r.x))<1e-10:
            x=r.x[:4]/np.linalg.norm(r.x[:4])
            if all(np.linalg.norm(x-y)>1e-6 for y,_ in out): out.append((x,r.x[4]))
    return out
def main():
    # This successor deliberately uses exact static-ansatz times. The original
    # pre-freeze script remains recoverable at commit 2e984ac.
    from experiments.closure_ledger.r3_prefreeze.r2_degree_intervals import main as intervals
    intervals()


if __name__ == '__main__':
    main()
