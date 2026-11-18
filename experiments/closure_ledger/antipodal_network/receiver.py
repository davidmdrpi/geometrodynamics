"""One driven port, one passive linear load on the connected mesh."""
import json
import numpy as np
from .refinement import domain
from .network import mesh,unit

def rim_area():
    # Each physical rim is represented by one lift. Sum the actual triangles.
    from collections import Counter
    from itertools import combinations
    cells,ports,_=mesh(8)
    faces=Counter(tuple(sorted(face)) for ids,_,r in cells if r=='tube'
                  for face in combinations(ids,3))
    area=0.
    for face,count in faces.items():
        if count==1 and all(p in ports[0] and p[0]>0 for p in face):
            edges=np.array([unit(p)-unit(face[0]) for p in face[1:]])
            area+=np.sqrt(np.linalg.det(edges@edges.T))/2
    return area

def solve(h,e,w,tau):
    phase=np.diag([1,np.exp(1j*w*tau)])
    zh=phase @ h @ phase.conj().T
    z=zh+e
    velocity_impedance=1.
    ub=-z[1,0]/(z[1,1]-1j*w*velocity_impedance)
    u=np.array([1,ub])
    fh,fe=zh @ u,e @ u
    def out(f):
        return w/2*np.imag(f*u.conj())
    ht,et=out(fh),out(fe)
    received=w*w*velocity_impedance*abs(ub)**2/2
    assert max(abs(ht.sum()),abs(et.sum()),abs(ht[1]+et[1]-received))<1e-9
    return dict(omega=w,tau=tau,receiver_power=received,tube_power_to_B=ht[1],
                exterior_power_to_B=et[1],source_power=-ht[0]-et[0])

def run(roots):
    area=rim_area()
    rows=[]
    for level in range(3):
        (tube,exterior),n,nt=domain(level,full_ports=True)
        omega=1.
        h,e=tube(omega),exterior(omega)
        analytical=-area*omega/np.sin(.95*omega)
        response=[solve(h,e,omega,tau) for tau in (0,.2375,.7125)]
        w=roots[level]['omega']
        h1,e1=tube(w),exterior(w)
        response.extend(solve(h1,e1,w,tau) for tau in (0,np.pi/w))
        assert response[1]['receiver_power']>.01 and response[1]['tube_power_to_B']>0
        assert response[-2]['receiver_power']>.01 and response[-1]['receiver_power']<1e-16
        low=solve(tube(.2),exterior(.2),.2,.2375)
        assert low['tube_power_to_B']>0
        rows.append(dict(low_frequency_control=low,level=level,nodes=n,cover_tetrahedra=nt,rim_area=area,
                         static_tube_offdiag=tube(0)[0,1], exact_static=-area/.95,
                         tube_offdiag=h[0,1],exact_tube_offdiag=analytical,runs=response))
    errs=[abs(r['tube_offdiag']-r['exact_tube_offdiag']) for r in rows]
    assert all(3 < a/b < 5 for a,b in zip(errs,errs[1:]))
    return rows
if __name__=='__main__':
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("refinement_json", type=Path, help="results of the full-domain refinement")
    args = parser.parse_args()
    roots = json.loads(args.refinement_json.read_text())
    if [row["level"] for row in roots] != list(range(3)):
        parser.error("expected one calibration each for levels 0, 1, 2 in order")
    print(json.dumps(run(roots), indent=2))
