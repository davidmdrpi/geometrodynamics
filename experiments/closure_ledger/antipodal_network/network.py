"""Closed tetrahedral S3 surgery with two exchanged tubes. Standalone calculation."""
import itertools as it
import json
from collections import defaultdict, deque
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import splu

def prisms(lo, hi):
    # A common ordering on a triangulated rim makes neighbouring diagonals agree.
    return [(*lo,hi[2]), (lo[0],lo[1],hi[1],hi[2]), (lo[0],*hi)]

def unit(p):
    a = np.asarray(p,dtype=float)
    return a/np.linalg.norm(a)

def rotate(p):
    x,y,z,w = p
    return (-y,x,z,w)

def negate(p):
    return tuple(-x for x in p)

def mesh(layers):
    facets = []
    for signs in it.product((-1,1),repeat=4):
        verts = np.diag(24*np.array(signs))
        for perm in it.permutations(range(4)):
            p = tuple(tuple(verts[list(perm[:j])].sum(axis=0)//j) for j in range(1,5))
            facets.append((perm[0],p))
    cells, rims = [], defaultdict(list)
    def bulk(tet):
        cells.append((list(tet),np.array([unit(p) for p in tet]),'exterior'))
    for axis,p in facets:
        if axis>1:
            bulk(p)
        else:
            rim = tuple(tuple((np.array(p[0])+v)//2) for v in p[1:])
            rims[p[0]].append(rim)
            for tet in prisms(rim,p[1:]):
                bulk(tet)
    ports = [set(),set()]
    length = .95
    for sign in (-1,1):
        for tri in rims[(24*sign,0,0,0)]:
            ports[0].update(tri)
            ports[1].update(rotate(v) for v in tri)
            def vid(p,i):
                if i==0:
                    return p
                if i==layers:
                    return rotate(p)
                return (i,*p)
            for i in range(layers):
                for prism in prisms([(i,p) for p in tri], [(i+1,p) for p in tri]):
                    ids = [vid(p,j) for j,p in prism]
                    coords = np.array([(*unit(p),length*j/layers) for j,p in prism])
                    cells.append((ids,coords,'tube'))
    return cells,ports,length

def topology(tets):
    faces = defaultdict(list)
    for i,tet in enumerate(tets):
        for j in range(4):
            face = tet[:j]+tet[j+1:]
            inversions = sum(a>b for a,b in it.combinations(face,2))
            faces[tuple(sorted(face))].append((i,(-1)**(j+inversions)))
    assert all(len(v)==2 for v in faces.values())
    adjacent=defaultdict(list)
    for (i,p),(j,q) in faces.values():
        adjacent[i].append((j,-p*q))
        adjacent[j].append((i,-p*q))
    signs={0:1}
    todo=deque([0])
    orientable=True
    while todo:
        i=todo.popleft()
        for j,rel in adjacent[i]:
            expected=signs[i]*rel
            if j in signs:
                orientable &= signs[j]==expected
            else:
                signs[j]=expected
                todo.append(j)
    assert len(signs)==len(tets)
    return orientable

def deck(p):
    return negate(p) if len(p)==4 else (p[0],*negate(p[1:]))

def build(layers, beta):
    cells,ports,length=mesh(layers)
    verts=sorted({p for ids,_,_ in cells for p in ids})
    idx={p:i for i,p in enumerate(verts)}
    tets=[tuple(idx[p] for p in ids) for ids,_,_ in cells]
    assert not topology(tets)
    simplex_set={frozenset(t) for t in tets}
    for ids,_,_ in cells:
        assert frozenset(idx[deck(p)] for p in ids) in simplex_set
        assert not (set(ids) & {deck(p) for p in ids})  # rules out fixed simplices
    assert not (ports[0] & ports[1]) and all({deck(p) for p in port}==port for port in ports)
    # Assemble once on both parity sectors using orbit representatives.
    reps=sorted({min(p,deck(p)) for p in verts})
    lookup={p:i for i,p in enumerate(reps)}
    n=len(reps)
    def orbit(p):
        return lookup[min(p,deck(p))]
    mats = {par:{region:[np.zeros((n,n)) for _ in range(3)] for region in ('tube','exterior')}
            for par in (1,-1)}
    for ids,coords,region in cells:
        edges=coords[1:]-coords[0]
        gram=edges @ edges.T
        volume=np.sqrt(np.linalg.det(gram))/6
        gradients=np.vstack((-np.ones(3),np.eye(3))) @ np.linalg.solve(gram,edges)
        drift=beta if region=='tube' else 0.
        mass=volume*(1-drift**2)*np.eye(4)/4
        coupling=np.tile(drift*volume*gradients[:,-1]/4,(4,1))
        gyro=coupling.T-coupling
        stiff=volume*gradients @ gradients.T
        sel=np.ix_([orbit(p) for p in ids],[orbit(p) for p in ids])
        for par in (1,-1):
            lift=np.array([1 if p==min(p,deck(p)) else par for p in ids])
            outer=np.outer(lift,lift)/2   # physical quotient volume
            for target,part in zip(mats[par][region],(mass,gyro,stiff)):
                target[sel] += outer*part
    return mats,reps,[sorted({orbit(p) for p in port}) for port in ports],length,len(cells)

def run(layers=8,beta=.25,dt=.004,stop=1.2, launch="tube"):
    mats,reps,ports,length,ncells=build(layers,beta)
    tm,tg,tk = mats[1]['tube']
    m,g,k = [a+b for a,b in zip(mats[1]['tube'],mats[1]['exterior'])]
    n=len(reps)
    def quadratic(u,v,mass,stiff):
        return (v @ mass @ v + u @ stiff @ u)/2
    def static_spectrum(parts):
        mass=parts['tube'][0]/(1-beta**2)+parts['exterior'][0]
        stiff=parts['tube'][2]+parts['exterior'][2]
        return np.sqrt(np.maximum(0,eigh(stiff,mass,subset_by_index=(0,5),eigvals_only=True)))
    # Smooth compact pulse fully inside both lifted tubes. Velocity travels
    # toward increasing z in the tube with local null slope dt/dz=1+beta.
    u=np.zeros(n)
    v=np.zeros(n)
    for j,p in enumerate(reps):
        if len(p)==5:
            z=p[0]/layers
            if 0 < z < .75:
                u[j]=np.sin(np.pi*z/.75)**4
                derivative=4*np.pi/(.75*length)*np.sin(np.pi*z/.75)**3*np.cos(np.pi*z/.75)
                v[j]=-derivative/(1+beta)
    if launch == "exterior":
        # Four exterior directions adjacent to A and angularly distant from B.
        # Mass lumping gives exactly zero initial energy on tube elements.
        u.fill(0)
        v.fill(0)
        for j,p in enumerate(reps):
            if len(p)==4 and abs(p[0])==12 and p[1]==0 and max(abs(p[2]),abs(p[3]))==12:
                v[j]=1
    energy=quadratic(u,v,m,k)
    u,v=u/np.sqrt(energy),v/np.sqrt(energy)
    et0=quadratic(u,v,tm,tk)
    system=sparse.csr_matrix(m + dt*g/2 + dt*dt*k/4)
    factor=splu(system.tocsc())
    drift=0
    balance=0
    work=np.zeros(2)
    max_exit=np.zeros(2)
    min_exit=np.zeros(2)
    peak_tube=0.
    peak_bulk=0.
    for _ in range(round(stop/dt)):
        dv=factor.solve(-dt*(g @ v + k @ (u+dt*v/2)))
        vm=v+dv/2
        um=u+dt*vm/2
        reaction=tm @ (dv/dt) + tg @ vm + tk @ um
        for j,ps in enumerate(ports):
            work[j] -= dt * (vm[ps] @ reaction[ps])  # work out of tube
        u,v=u+dt*vm,v+dv
        et=quadratic(u,v,tm,tk)
        peak_bulk=max(peak_bulk,1-et)
        max_exit=np.maximum(max_exit,work)
        min_exit=np.minimum(min_exit,work)
        peak_tube=max(peak_tube,et)
        balance=max(balance,abs(et-et0+sum(work)))
        drift=max(drift,abs(quadratic(u,v,m,k)-1))
    return dict(launch=launch, tube_layers=layers, cover_tetrahedra=ncells,quotient_vertices=n,
                tube_clock_offset=beta*length, step=dt,final_time=round(stop/dt)*dt,
                zero_shift_even_frequencies=static_spectrum(mats[1]).tolist(),
                zero_shift_odd_frequencies=static_spectrum(mats[-1]).tolist(),
                initial_tube_energy=et0,max_exterior_energy=peak_bulk,
                net_out_at_ports=work.tolist(),max_cumulative_out_at_ports=max_exit.tolist(),
                min_cumulative_out_at_ports=min_exit.tolist(),max_tube_energy=peak_tube,
                max_total_energy_error=drift,max_tube_port_balance_error=balance)

if __name__=='__main__':
    rows=[run(n,dt=d) for n,d in ((8,.004),(16,.004),(16,.002),(32,.002))]
    rows.append(run(16,beta=0.,dt=.004))
    for row in rows:
        assert row['max_total_energy_error']<1e-10 and row['max_tube_port_balance_error']<1e-10
        assert abs(row['initial_tube_energy']-1)<1e-12 and row['max_exterior_energy']>.1
        assert row['zero_shift_odd_frequencies'][0]>.1
    for duration in (1.2,2.,3.):
        row=run(16,dt=.004,stop=duration, launch="exterior")
        assert row['initial_tube_energy']==0 and row['max_tube_energy']>.01
        assert max(row['max_total_energy_error'],row['max_tube_port_balance_error'])<1e-10
        rows.append(row)
    print(json.dumps(rows,indent=2))
