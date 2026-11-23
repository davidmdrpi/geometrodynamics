"""Uniform refinement of both domains, independent sparse port elimination."""
import json
import itertools as it
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import LinearOperator, eigsh, splu
from .network import mesh, deck

EDGES=list(it.combinations(range(4),2))
# Vertex indices 4..9 correspond to edge midpoints 01,02,03,12,13,23.
CHILDREN=np.array([[0,4,5,6],[1,4,7,8],[2,5,7,9],[3,6,8,9],
                   [4,9,5,6],[4,9,6,8],[4,9,8,7],[4,9,7,5]])

def _refined_mesh(level, *, quotient):
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
    if quotient:
        reps, orbit = np.unique(np.minimum(np.arange(len(port)), antipode), return_inverse=True)
        nodes = orbit[nodes]
        port, seed = port[reps], seed[reps]
    n = len(port)
    # Keep only true trace nodes shared by the two subdomains.
    interface=np.intersect1d(nodes[tube],nodes[~tube])
    assert np.all(port[interface]>0)
    traces=np.zeros((n,2))
    traces[interface,port[interface]-1]=1
    return nodes, coords, tube, seed, traces


def _region_matrices(nodes, gram, n, *, quotient):
    volume = np.sqrt(np.linalg.det(gram)) / (12 if quotient else 6)
    d = np.vstack((-np.ones(3), np.eye(3)))
    ke = (d @ np.linalg.inv(gram) @ d.T) * volume[:, None, None]
    row = np.broadcast_to(nodes[:, :, None], ke.shape).ravel()
    col = np.broadcast_to(nodes[:, None, :], ke.shape).ravel()
    k = coo_matrix((ke.ravel(), (row, col)), shape=(n, n)).tocsc()
    mass = np.bincount(nodes.ravel(), weights=np.repeat(volume / 4, 4), minlength=n)
    return k, mass, volume, row, col


def exterior_spectrum(level, *, quotient=False):
    """Clamped polyhedral rims, with row-sum lumping and consistent P1 mass."""
    nodes, coords, tube, _, traces = _refined_mesh(level, quotient=quotient)
    n = len(traces)
    edges = coords[~tube, 1:] - coords[~tube, :1]
    gram = edges @ edges.transpose(0, 2, 1)
    k, mass, volumes, row, col = _region_matrices(nodes[~tube], gram, n, quotient=quotient)
    local_mass = volumes[:, None, None] * ((1 + np.eye(4)) / 20)
    consistent = coo_matrix((local_mass.ravel(), (row, col)), shape=(n, n)).tocsc()
    keep = (mass > 0) & ~traces.any(axis=1)
    stiffness = k[keep, :][:, keep]
    # Reuse the zero-shift factorization between the two generalized eigenproblems.
    inverse = LinearOperator(stiffness.shape, matvec=splu(stiffness).solve)

    def frequencies(m):
        lam = eigsh(stiffness, k=2, M=m, sigma=0, OPinv=inverse, which="LM",
                    v0=np.random.default_rng(0).normal(size=stiffness.shape[0]),
                    return_eigenvectors=False)
        return np.sqrt(np.sort(lam)).tolist()

    return dict(level=level, space="quotient" if quotient else "cover", cover_tetrahedra=len(nodes),
                lumped=frequencies(diags(mass[keep], format="csc")),
                consistent=frequencies(consistent[keep, :][:, keep]))


def domain(level, dynamic=False, full_ports=False):
    nodes, coords, tube, seed, traces = _refined_mesh(level, quotient=True)
    n = len(traces)
    regional = []
    for mask in (tube, ~tube):
        edges=coords[mask,1:]-coords[mask,:1]
        gram=edges @ edges.transpose(0,2,1)
        k, mass, volume, row, col = _region_matrices(nodes[mask], gram, n, quotient=True)
        if dynamic:
            beta=.25 if mask is tube else 0.
            grad=np.vstack((-np.ones(3),np.eye(3))) @ np.linalg.solve(gram,edges)
            local=beta*volume[:,None,None]*np.broadcast_to(grad[:,:,4][:,None,:],(len(volume),4,4))/4
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
        return regional,n,len(nodes),seed,traces
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

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", choices=("cover", "quotient"),
                        help="clamp both rims and compute the two lowest exterior frequencies")
    parser.add_argument("levels", type=int, nargs="*", help="red-refinement levels")
    args = parser.parse_args()
    if any(level < 0 for level in args.levels):
        parser.error("refinement levels must be nonnegative")
    if args.spectrum:
        results = [exterior_spectrum(level, quotient=args.spectrum == "quotient")
                   for level in (args.levels or range(4))]
        print(json.dumps(results, indent=2))
    elif len(args.levels) == 1:
        sample(args.levels[0])
    elif args.levels:
        parser.error("the single-level frequency scan takes one level")
    else:
        print(json.dumps([sample(j, win) for j, win in enumerate(WINDOWS)], indent=2))
