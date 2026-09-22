"""Prospective refinement dbec68f; original 6/8 remains intact."""
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves import localized_mouth as m
from . import localized_mouth_probe as p
from .evidence_archive import read_bytes

FREEZE='dbec68f2b6de8f745c269167a2ac7a41a38f4654'
ORIGINAL_HASH='e510ba44aeb7d8d8d25cfc28cfb54a6a402d45ba43d373d5f6b48b6bfdffd227'
ORIGINAL_GATES={k:k not in ('hamiltonian','physical') for k in p.GATES}


def reconstruction_agrees(actual,expected):
    """Exact metadata/mesh; roundoff-only bounds on the complete polynomials.

    CPU-dispatched NumPy arithmetic can differ in its last bits. Compare
    interval-scaled coefficients, not raw power coefficients (which contain
    inverse powers of the tiny interval width). The sum of absolute scaled
    coefficient differences bounds the error everywhere in that interval.
    Check psi through its second derivative and the stored velocity through
    its first derivative, without relying on a finite set of sample points.
    """
    if not isinstance(actual,dict) or not isinstance(expected,dict):
        return False,float('inf')
    if actual.keys()!=expected.keys():return False,float('inf')
    if any(actual[k]!=expected[k] for k in expected if k!='solution'):
        return False,float('inf')
    a,e=actual['solution'],expected['solution']
    if not isinstance(a,dict) or not isinstance(e,dict):
        return False,float('inf')
    if a.keys()!=e.keys() or any(a[k]!=e[k] for k in e if k!='c'):
        return False,float('inf')
    ac,ec=np.asarray(a['c'],dtype=float),np.asarray(e['c'],dtype=float)
    h=np.diff(np.asarray(e['x'],dtype=float))
    if (ac.shape!=ec.shape or ec.shape!=(6,len(h),2)
            or not np.isfinite(ac).all() or not np.isfinite(ec).all()
            or not np.isfinite(h).all() or not np.all(h>0)):
        return False,float('inf')
    # Arithmetic agreement only: far below the unchanged 1e-12 knot and
    # 1e-7 PDE thresholds. Never use this budget for a scientific gate.
    roundoff=32*np.finfo(float).eps
    worst=0.
    for component,orders in ((0,range(3)),(1,range(2))):
        for order in orders:
            powers=np.arange(5,order-1,-1)
            factors=np.ones(len(powers))
            for j in range(order):factors*=powers-j
            weights=factors[:,None]*h[None,:]**(powers[:,None]-order)
            delta=np.sum(abs((ac[:len(powers),:,component]-ec[:len(powers),:,component])*weights),axis=0)
            scale=np.maximum(1.,np.sum(abs(ec[:len(powers),:,component]*weights),axis=0))
            worst=max(worst,float(np.max(delta/(roundoff*scale))))
    return worst<=1.,worst


def read_original(path):
    raw=read_bytes(path);raw=gzip.decompress(raw) if path.suffix=='.gz' else raw
    if hashlib.sha256(raw).hexdigest()!=ORIGINAL_HASH:raise ValueError('original archive hash mismatch')
    data=json.loads(raw)
    if p.score(data)['gates']!=ORIGINAL_GATES:raise ValueError('original verdict changed')
    return data


def run(original):
    data=copy.deepcopy(original)
    for i,r in enumerate(data['solutions']):
        if r['n_initial']==513:data['solutions'][i]=m.reconstruct(r,data['profiles'][str(r['L'])])
    fine=[r for r in data['solutions'] if r['n_initial']==513]
    data['diagnostics']=[dict(L=r['L'],eta=r['eta'],data=m.Data(r,data['profiles'][str(r['L'])]).diagnostics()) for r in fine]
    data['seams']=[dict(L=r['L'],eta=r['eta'],data=m.seam(m.Data(r,data['profiles'][str(r['L'])]))) for r in fine]
    target=next(r for r in fine if (r['L'],r['eta'])==(5.5,.3))
    d=m.Data(target,data['profiles']['5.5'])
    data['physical']=[dict(h=h,cases=[m.physical_constraints(d,(frac*5.5,t,.37),h,order=4)
            for frac in (.07,.23,.51,.79,.93) for t in (.43,.91,1.47,2.13)]) for h in (.008,.004,.002)]
    data['refinement_freeze']=FREEZE;data['original_sha256']=ORIGINAL_HASH
    return data


def score(data,original):
    result=p.score(data,derivative_order=4)
    integrity=True;knots=[];changes=[];roundoff_ratios=[]
    try:
        integrity &= data['refinement_freeze']==FREEZE and data['original_sha256']==ORIGINAL_HASH
        original_bytes=json.dumps(original,sort_keys=True,indent=2,allow_nan=False).encode()+b'\n'
        integrity &= hashlib.sha256(original_bytes).hexdigest()==ORIGINAL_HASH
        integrity &= p.score(original)['gates']==ORIGINAL_GATES
        for before,after in zip(original['solutions'],data['solutions']):
            if before['n_initial']!=513:
                integrity &= before==after
                continue
            L=before['L'];profiles=original['profiles'][str(L)]
            expected=m.reconstruct(before,profiles)
            agrees,ratio=reconstruction_agrees(after,expected)
            integrity &= agrees
            roundoff_ratios.append(ratio)
            b=m.Data(before,profiles);a=m.Data(after,profiles)
            atknots=a.sol(b.sol.x);old=b.sol(b.sol.x)
            knots.append(p.peak(atknots-old))
            s=np.linspace(0,L,1001)
            changes.append(p.peak((a.sol(s)[0]-b.sol(s)[0])/b.sol(s)[0]))
        integrity &= len(data['solutions'])==27 and len(knots)==9
        integrity &= max(knots)<1e-12 and max(changes)<1e-8
    except (KeyError,ValueError,TypeError,IndexError,ArithmeticError):
        integrity=False
    result['gates']['evidence'] &= bool(integrity)
    result['metrics']['reconstruction_knot_error']=max(knots) if knots else None
    result['metrics']['reconstruction_relative_change']=max(changes) if changes else None
    ratio=max(roundoff_ratios) if roundoff_ratios else float('inf')
    result['metrics']['reconstruction_roundoff_budget_used']=ratio if np.isfinite(ratio) else None
    result['passed']=sum(result['gates'].values())
    result['original_gates']=dict(ORIGINAL_GATES)
    result['refinement_freeze']=FREEZE
    passed=all(v for k,v in result['gates'].items() if k!='localization')
    result['verdicts']={'FOUR_SCALAR_HANDLE_CONSTRAINT_DATA':passed,'LOCALIZED_BULK_MOUTH_INITIAL_DATA':passed and result['gates']['localization']}
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--original',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--rescore',type=Path)
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    md=args.output_dir/'refinement.md';js=args.output_dir/'refinement_verdict.json'
    md.write_text('# Refinement incomplete\n\nNo affirmative verdict.\n')
    try:
        original=read_original(args.original)
        if args.rescore:
            raw=read_bytes(args.rescore);data=json.loads(gzip.decompress(raw) if args.rescore.suffix=='.gz' else raw)
        else:data=run(original)
        result=score(data,original)
        raw=json.dumps(data,sort_keys=True,indent=2,allow_nan=False).encode()+b'\n'
        (args.output_dir/'refinement.json.gz').write_bytes(gzip.compress(raw,mtime=0))
        js.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        md.write_text('# Prospective localized-mouth refinement\n\nOriginal: **6/8**, unchanged.\n\n'+f"Refinement: **{result['passed']}/8**.\n\n"+'```json\n'+json.dumps(result,indent=2)+'\n```\n')
        print(json.dumps({'passed':result['passed'],'gates':result['gates'],'metrics':result['metrics'],'raw_sha256':hashlib.sha256(raw).hexdigest()},indent=2))
        return 0 if all(result['gates'].values()) else 1
    except Exception as error:
        md.write_text(f'# Refinement failed\n\nNo affirmative verdict.\n\n{error}\n')
        js.write_text(json.dumps({'error':str(error),'verdicts':{'FOUR_SCALAR_HANDLE_CONSTRAINT_DATA':False,'LOCALIZED_BULK_MOUTH_INITIAL_DATA':False}})+'\n')
        raise


if __name__=='__main__':raise SystemExit(main())
