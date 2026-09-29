"""R2 falsifier: global degree of phi/|phi| on closed S^3 for phi = R(eta) x + eps*dphi(x).
Checks: (a) odd fields have odd degree (Borsuk-Ulam); (b) zero events cluster at the
collapse R=0 and their number is eps-independent; (c) Hopf (complex-structure) perturbation
gives no events."""
import numpy as np
from scipy.optimize import root
rng=np.random.default_rng(7)

def grid(n=72):
    a=(np.arange(n)+.5)*(np.pi/2)/n; b=(np.arange(2*n)+.5)*2*np.pi/(2*n)
    A,B1,B2=np.meshgrid(a,b,b,indexing='ij')
    X=np.stack([np.cos(A)*np.cos(B1),np.cos(A)*np.sin(B1),np.sin(A)*np.cos(B2),np.sin(A)*np.sin(B2)],-1)
    return X,(a[1]-a[0],b[1]-b[0])

def degree(F,h):
    Fh=F/np.linalg.norm(F,axis=-1,keepdims=True)
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

X,h=grid()
Rp=np.sqrt(3)  # R ~ -sqrt(3)*tau near eta=pi/4
def tangential_zeros(f,seeds=4000):
    """zeros of f_perp on S^3: solve f(x)-(x.f)x=0 with |x|=1, dedupe; return x and a=x.f"""
    S=rng.normal(size=(seeds,4)); S/=np.linalg.norm(S,axis=1,keepdims=True); out=[]
    for s in S:
        g=lambda z:np.r_[(lambda x,fx:(fx-(x@fx)*x))(z,f(z[None])[0])[:4]+0*z[:4], z@z-1][[0,1,2,3,4]]
        def G(z):
            x=z[:4];lam=z[4];fx=f(x[None])[0];return np.r_[fx-lam*x,x@x-1]
        r=root(G,np.r_[s,s@f(s[None])[0]],tol=1e-13)
        if r.success and np.linalg.norm(G(r.x))<1e-10:
            x=r.x[:4]/np.linalg.norm(r.x[:4])
            if all(np.linalg.norm(x-y)>1e-6 for y,_ in out): out.append((x,r.x[4]))
    return out
Z=tangential_zeros(dphi)
print('zeros of tangential part of dphi on S^3:',len(Z))
for eps in (1e-3,1e-2,5e-2):
    taus=sorted(eps*lam/Rp for _,lam in Z)   # R(tau)x+eps*lam x=0 -> -sqrt3 tau + eps lam=0
    print(f'eps={eps}: event times tau/eps =',np.round(np.array(taus)/eps,4))
# degree across the collapse (tau grid avoiding events), eps=1e-2
eps=1e-2; D=dphi(X)
lams=sorted(l for _,l in Z); edges=[-1]+[l*eps/Rp for l in lams]+[1]
print('degree in each interval between events (eps=1e-2):')
for lo,hi in zip(edges,edges[1:]):
    t=np.clip((lo+hi)/2,-.05,.05); R=(np.sqrt(3)/2)*np.cos(2*(np.pi/4+t))
    print(f'  tau={t:+.5f}  N={degree(R[...,None]*X+eps*D if np.ndim(R) else R*X+eps*D,h):+.4f}')
for t in (-.02,0.,.02):
    R=(np.sqrt(3)/2)*np.cos(2*(np.pi/4+t))
    print(f'Hopf perturbation tau={t:+.3f}: N={degree(R*X+eps*hopf(X),h):+.4f}, min|phi|={np.linalg.norm(R*X+eps*hopf(X),axis=-1).min():.3e}')
