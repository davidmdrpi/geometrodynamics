"""Uniform refinement of both domains, independent sparse port elimination."""
import json
import itertools as it
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu
from .network import mesh, deck

EDGES=list(it.combinations(range(4),2))
# Vertex indices 4..9 correspond to edge midpoints 01,02,03,12,13,23.
CHILDREN=np.array([[0,4,5,6],[1,4,7,8],[2,5,7,9],[3,6,8,9],
                   [4,9,5,6],[4,9,6,8],[4,9,8,7],[4,9,7,5]])

def domain(level, dynamic=False, full_ports=False):
    cells, ports, _=mesh(8)
    keys=sorted({p for ids,_,_ in cells for p in ids})
    ix={p:j for j,p in enumerate(keys)}
    antipode=np.array([ix[deck(p)] for p in keys])
    port=np.zeros(len(keys),int)
    for j,ps in enumerate(ports,1):
        port[[ix[p] for p in ps]]=j
    seed=np.array([float(len(p)==4 and abs(p[0])==12 and p[1]==0 and max(abs(p[2]),abs(p[3]))==12) for p in keys])
    nodes=np.array([[ix[p] for p in ids] for ids,_,_ in cells])
    coords=np.array([np.pad(c,((0,0),(0,5-c.shape[1]))) for _,c,_ in cells])
    tube=np.array([r=='tube' for _,_,r in cells])
    for _ in range(level):
        ed=np.sort(nodes[:,EDGES],axis=2)
        unique, inverse=np.unique(ed.reshape(-1,2),axis=0,return_inverse=True)
        old=len(antipode)
        look={tuple(pair):old+j for j,pair in enumerate(unique)}
        amid=[look[tuple(sorted(antipode[pair]))] for pair in unique]
        antipode=np.r_[antipode,amid]
        seed=np.r_[seed,seed[unique].mean(axis=1)]
        port=np.r_[port,np.where(port[unique[:,0]]==port[unique[:,1]],port[unique[:,0]],0)]
        nodes=np.concatenate((nodes,old+inverse.reshape(-1,6)),axis=1)[:,CHILDREN].reshape(-1,4)
        coords=np.concatenate((coords,coords[:,EDGES].mean(axis=2)),axis=1)[:,CHILDREN].reshape(-1,4,5)
        tube=np.repeat(tube,8)
    reps,orbit=np.unique(np.minimum(np.arange(len(port)),antipode),return_inverse=True)
    nodes=orbit[nodes]
    n=len(reps)
    ports_q=port[reps]
    regional=[]
    # Keep only true trace nodes shared by the two subdomains.
    interface=np.intersect1d(nodes[tube],nodes[~tube])
    assert np.all(ports_q[interface]>0)
    traces=np.zeros((n,2))
    traces[interface,ports_q[interface]-1]=1
    for mask in (tube,~tube):
        t=nodes[mask]
        edges=coords[mask,1:]-coords[mask,:1]
        gram=edges @ edges.transpose(0,2,1)
        volume=np.sqrt(np.linalg.det(gram))/12  # quotient has half the cover volume
        d=np.vstack((-np.ones(3),np.eye(3)))
        ke=(d @ np.linalg.inv(gram) @ d.T)*volume[:,None,None]
        row=np.broadcast_to(t[:,:,None],ke.shape).ravel()
        col=np.broadcast_to(t[:,None,:],ke.shape).ravel()
        k=coo_matrix((ke.ravel(),(row,col)),shape=(n,n)).tocsc()
        mass=np.bincount(t.ravel(),weights=np.repeat(volume/4,4),minlength=n)
        if dynamic:
            beta=.25 if mask is tube else 0.
            grad=d @ np.linalg.solve(gram,edges)
            local=beta*volume[:,None,None]*np.broadcast_to(grad[:,:,4][:,None,:],ke.shape)/4
            gyro=coo_matrix(((local.transpose(0,2,1)-local).ravel(),(row,col)),shape=(n,n)).tocsc()
            regional.append((diags(mass*(1-beta**2),format="csc"),gyro,k))
            continue
        interior=(mass>0)&(traces.sum(axis=1)==0)
        ki=k[interior,:][:,interior]
        mi=diags(mass[interior],format='csc')
        kb=k @ traces
        cross=traces[:,0] @ kb[:,1]
        if full_ports:
            def port_matrix(w,ki=ki,mi=mi,load=kb[interior,:],base=traces.T@kb,
                            mb=traces.T@(mass[:,None]*traces)):
                return base-w*w*mb-load.T@splu(ki-w*w*mi).solve(load)
            regional.append(port_matrix)
            continue
        def off(w,ki=ki,mi=mi,load=kb[interior,:],cross=cross):
            return cross-load[:,0] @ splu(ki-w*w*mi).solve(load[:,1])
        regional.append(off)
    if dynamic:
        return regional,n,len(nodes),seed[reps],traces
    return regional,n,len(nodes)

WINDOWS=[(4.575,4.6),(4.85,4.90),(5.098,5.10)]

def sample(level,window=None):
    (h,e),n,nt=domain(level)
    if window is None:
        for om in np.arange(3.55,5.16,.05):
            print(round(om,3),round(h(om),3),round(e(om),3),flush=True)
        return
    w=brentq(lambda f:h(f)-e(f),*window,xtol=1e-14)
    H,E=h(w),e(w)
    assert min(abs(H),abs(E))>.1 and abs(H-E)<1e-8
    calibrated=4.599769873349991
    H0,E0=h(calibrated),e(calibrated)
    return dict(coarse_clock_cross_reaction=abs(E0-H0),coarse_zero_clock_cross_reaction=abs(E0+H0),level=level,cover_tetrahedra=nt,quotient_nodes=n,omega=w,
                offset=np.pi/w,channel_cross_reactions=[H,E],bright=abs(H+E),
                dark=abs(H-E),root_bracket=window)

if __name__=='__main__':
    import sys
    if len(sys.argv)>1:
        sample(int(sys.argv[1]))
    else:
        print(json.dumps([sample(j,win) for j,win in enumerate(WINDOWS)],indent=2))
