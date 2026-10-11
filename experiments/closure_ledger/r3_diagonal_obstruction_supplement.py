"""Disclosed post-inspection seed audit and prospective supplemental root checks.

This does not alter the registered score. It calls the frozen root algorithm
on the already observed sign brackets despite the original locality failure.
"""
import argparse
import hashlib
import json
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import gmpy2 as gp
from gmpy2 import mpfr as mp
from scipy.interpolate import CubicSpline
from experiments.closure_ledger import r3_diagonal_obstruction_probe as p
from geometrodynamics.waves import r3_diagonal_symmetry as ds

PATH='experiments/closure_ledger/r3_diagonal_obstruction_supplement.py'
NOTE='docs/r3_diagonal_obstruction_followup.md'


def audit(rows):
    v=np.array([r['v'] for r in json.loads((p.ROOT/'experiments/closure_ledger/runs/20261001_r3_family/stage_F.json').read_text())['points']])
    vc=np.vstack([v,v[:1]]);s=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(vc[:,:6],axis=0),axis=1))]
    curve=CubicSpline(2*np.pi*s/s[-1],vc,bc_type='periodic');tangent=curve.derivative()
    f20,_=ds.fourier_fit(v[:,:6],20);records=[]
    for r in rows:
        seeds,c,t=p.control_seeds(r['phi']);coherent=curve(r['phi']).reshape(2,6);Z=np.array(r['check']['nodes'],float)
        records.append(dict(j=r['j'],old_seed_distance=r['seed_distance'],node0_correction=float(np.max(abs(Z[0]-seeds[0]))),
                            coherent_seed_distance=float(np.max(abs(Z-coherent))),
                            slice_anchor_difference=float(np.max(abs(coherent[0]-c))),
                            slice_tangent_difference=float(np.max(abs(tangent(r['phi'])[:6]-t))),
                            old_component_fit_distance=float(np.max(abs(Z-f20(ds.angle(Z)))))))
    return dict(rows=records,original_gate_failures=sum(r['old_seed_distance']>.01 for r in records),
                maxima={k:max(r[k] for r in records) for k in records[0] if k!='j'})


def inputs():
    paths=[PATH,NOTE]+[str((p.RUN/n).relative_to(p.ROOT)) for n in ('scan.json','scan_result.json','roots.json','result.json')]
    return {n:hashlib.sha256((p.ROOT/n).read_bytes()).hexdigest() for n in paths}


def score(rows,primary):
    complete=[r['j'] for r in rows]==primary['brackets']
    all_ok=bool(rows and complete and all(r['diagnostics']['ok'] for r in rows))
    groups=dict(P_C=p.equivalence_classes(rows),H_C=p.equivalence_classes(rows,True)) if all('nodes' in r for r in rows) else {}
    distinct=bool(groups and rows and groups['P_C']['min_distinct_separation']>1e-6)
    return dict(label='SUPPLEMENTAL_ROOTS_VERIFIED' if all_ok and distinct else 'SUPPLEMENTAL_ROOTS_UNRESOLVED',
                root_count=len(rows),all_isolated=bool(all_ok and distinct),equivalence=groups,
                leading_twelfth_pattern=bool(all_ok and distinct and len(rows)==24 and not primary['ambiguous_phases']))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['audit','roots']);parser.add_argument('--workers',type=int,default=8);args=parser.parse_args()
    scan=json.loads((p.RUN/'scan.json').read_text());primary=json.loads((p.RUN/'scan_result.json').read_text())
    if scan['bindings']!=p.bindings():raise ValueError('frozen sources differ')
    if args.stage=='audit':
        result=audit(scan['rows']);p.write(p.RUN/'seed_audit.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2));return
    fingerprint=inputs();jobs=[(j,scan['rows'][j],scan['rows'][(j+1)%p.N]) for j in primary['brackets']]
    rows=[]
    with Pool(args.workers) as pool:
        for r in pool.imap_unordered(p.safe_root,jobs):rows.append(r);print('supplement root complete',r['j'],flush=True)
    rows.sort(key=lambda r:r['j']);result=score(rows,primary)
    p.write(p.RUN/'supplement.json',dict(inputs=fingerprint,registered_label=primary['label'],rows=rows,result=result))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
