"""Frozen checkpointed diagonal obstruction scan and isolated-root follow-up."""
import argparse
import hashlib
import json
from multiprocessing import Pool
from pathlib import Path
import platform
import time
import numpy as np
import scipy
import gmpy2 as gp
from gmpy2 import mpfr as mp
from geometrodynamics.waves import r3_diagonal_taylor as dt
from geometrodynamics.waves import r3_high_precision_shooting as sh
from geometrodynamics.waves import r3_diagonal_symmetry as ds
from geometrodynamics.waves import r3_breaking as b
from experiments.closure_ledger.r3_breaking_probe import control_seeds

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'experiments/closure_ledger/runs/20261011_r3_diagonal_obstruction'
N=192
SOURCES=['geometrodynamics/waves/r3_diagonal_taylor.py','geometrodynamics/waves/r3_high_precision_shooting.py',
         'geometrodynamics/waves/r3_diagonal_symmetry.py','geometrodynamics/waves/r3_breaking.py',
         'experiments/closure_ledger/r3_breaking_probe.py','experiments/closure_ledger/r3_diagonal_obstruction_probe.py',
         'experiments/closure_ledger/r3_diagonal_obstruction_controls.py','docs/r3_diagonal_obstruction_prereg.md',
         'tests/test_r3_diagonal_obstruction.py','experiments/closure_ledger/runs/20261001_r3_family/stage_F.json',
         'experiments/closure_ledger/runs/20261011_r3_diagonal_obstruction/controls.json']


def bindings():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data,indent=1,allow_nan=False)+'\n');tmp.replace(path)


def numerical_map(cfg,ledger=None):
    def P(z):
        r=dt.clock_map(z,**cfg)
        if ledger is not None:ledger.append(r['constraint'])
        return r['z']
    return P


def double_point(seeds,c,t):
    P=lambda z:ds.clock_map(z,half=False)['z']
    point=b.scan_point(P,seeds,np.array(c,float),np.array(t,float),tol=5e-13,maxit=16)
    # The historical routine may record the residual before its last update.
    Z=np.array(point['nodes']);point['residual']=float(max(abs(b._residual([P(z) for z in Z],Z,point['lam'],np.array(point['g']),np.array(point['c0']),np.array(point['t0'])))))
    return point


def phase(j):
    started=time.time();fingerprint=bindings();path=RUN/'checkpoints'/f'phase_{j:03d}.json'
    if path.exists():
        saved=json.loads(path.read_text())
        if saved['bindings']!=fingerprint:raise ValueError('phase provenance mismatch')
        return saved
    pre=RUN/'checkpoints'/f'primary_{j:03d}.json'
    if pre.exists():
        rec=json.loads(pre.read_text())
        if rec['bindings']!=fingerprint:raise ValueError('primary provenance mismatch')
    else:
        phi=2*np.pi*j/N;seeds,c,t=control_seeds(phi);point=double_point(seeds,c,t)
        ledger=[];primary=sh.solve(numerical_map(dt.PRIMARY,ledger),point,dt.PRIMARY['bits'],'1e-45')
        rec=dict(j=j,phi=phi,point=point,primary=primary,bindings=fingerprint,
                 seed_distance=float(np.max(abs(np.array(point['nodes'])-seeds))),
                 primary_constraint=str(max(ledger)))
        write(pre,rec)
    ledger=[];check=sh.solve(numerical_map(dt.CHECK,ledger),rec['point'],dt.CHECK['bits'],'1e-65',start=rec['primary'])
    rec.update(check=check,check_constraint=str(max(ledger)),seconds=time.time()-started)
    with gp.context(gp.get_context(),precision=336):
        delta=abs(mp(rec['primary']['lam'])-mp(check['lam']))
        amplification=2*max(rec['primary']['inverse_row_norm'],check['inverse_row_norm'])
        err=delta+mp(amplification)*(mp(rec['primary']['residual'])+mp(check['residual']))
        rec['error_estimate']=str(err)
    write(path,rec);return rec


def scan_score(rows,controls):
    with gp.context(gp.get_context(),precision=336):
        if len(rows)!=N or [r['j'] for r in rows]!=list(range(N)):raise ValueError('incomplete phase grid')
        if any('error' in r for r in rows):
            return dict(label='NUMERICALLY_UNRESOLVED',gates=dict(execution=False),brackets=[],ambiguous_phases=[],failed_phases=[r['j'] for r in rows if 'error' in r])
        error=max(mp(r['error_estimate']) for r in rows);F=max(mp('1e-30'),error)
        lams=[mp(r['check']['lam']) for r in rows];amplitude=max(map(abs,lams))
        gates=dict(controls=controls['ok'],solves=all(r['primary']['ok'] and r['check']['ok'] and mp(r['primary']['residual'])<=mp('1e-45') and mp(r['check']['residual'])<=mp('1e-65') for r in rows),
                   error=bool(error<=mp('1e-30')),
                   conditioning=all(r['point']['cond']<=1e7 for r in rows),
                   locality=all(r['seed_distance']<=.01 for r in rows),
                   constraints=all(max(mp(r['primary_constraint']),mp(r['check_constraint']))<=mp('1e-48') for r in rows))
        brackets=[j for j in range(N) if lams[j]*lams[(j+1)%N]<0 and min(abs(lams[j]),abs(lams[(j+1)%N]))>10*F]
        ambiguous=[j for j,v in enumerate(lams) if abs(v)<=10*F]
        label='NUMERICALLY_UNRESOLVED' if not all(gates.values()) else ('RESOLVED_OBSTRUCTION' if amplitude>100*F else 'ADDITIONAL_SUPPRESSION_CANDIDATE')
        spec=abs(np.fft.rfft(np.array(lams,float)))/N
        return dict(label=label,gates=gates,error_estimate=str(error),F=str(F),detection_threshold=str(100*F),Lambda=str(amplitude),
                    sigma_effective=str(-gp.log(amplitude)/12) if amplitude else None,
                    sign_changes=b.sign_changes(np.array(lams,float)),brackets=brackets,ambiguous_phases=ambiguous,
                    dominant_raw_harmonic=int(np.argmax(spec[1:])+1),raw_spectrum=spec.tolist())


def isolated_diagnostics(P,Q,Z,step='1e-22',check_step='1e-24'):
    """Numerical nonsingularity and root correction estimates; no interval proof."""
    A=sh.closure_jacobian(P,Z,step);D=sh.closure_jacobian(Q,Z,check_step);inv=sh.inverse(D)
    difference=[[a-c for a,c in zip(row,other)] for row,other in zip(A,D)]
    bound=sh.matrix_norm(inv)*sh.matrix_norm(difference)
    residual=sh.closure(Q,Z);correction=sh.matvec(inv,residual)
    # Numerical inverse residual guards the inversion itself.
    n=len(D);identity_error=max(abs(gp.fsum(D[i][k]*inv[k][j] for k in range(n))-int(i==j)) for i in range(n) for j in range(n))
    return dict(inverse_norm=str(sh.matrix_norm(inv)),jacobian_difference=str(sh.matrix_norm(difference)),
                inverse_times_difference=str(bound),closure_residual=str(max(map(abs,residual))),
                correction_norm=str(max(map(abs,correction))),inverse_residual=str(identity_error),
                ok=bool(bound<mp('.01') and max(map(abs,residual))<mp('1e-45') and max(map(abs,correction))<mp('1e-12') and identity_error<mp('1e-40')))


def root(job):
    j,left,right=job;path=RUN/'checkpoints'/f'root_{j:03d}.json';fingerprint=bindings()
    if path.exists():
        rec=json.loads(path.read_text())
        if rec['bindings']!=fingerprint:raise ValueError('root provenance mismatch')
        return rec
    start=time.time()
    with gp.context(gp.get_context(),precision=336):
        cL,cR=([mp(v) for v in r['point']['c0']] for r in (left,right))
        tL,tR=([mp(v) for v in r['point']['t0']] for r in (left,right))
        lo,hi=mp(0),mp(1);fl,fh=mp(left['check']['lam']),mp(right['check']['lam'])
        history=[];candidate=None
        for _ in range(48):
            u=(lo*fh-hi*fl)/(fh-fl)
            if not lo<u<hi:u=(lo+hi)/2
            c=[(1-u)*a+u*v for a,v in zip(cL,cR)];t=[(1-u)*a+u*v for a,v in zip(tL,tR)]
            norm=gp.sqrt(gp.fsum(v*v for v in t));t=[v/norm for v in t]
            g=[t[k+1] if k%2==0 else -t[k-1] for k in range(6)]
            seeds=np.array(left['check']['nodes'],float)*(1-float(u))+np.array(right['check']['nodes'],float)*float(u)
            pt=double_point(seeds,c,t);pt.update(c0=list(map(str,c)),t0=list(map(str,t)),g=list(map(str,g)))
            candidate=sh.solve(numerical_map(dt.PRIMARY),pt,224,'1e-56',maxit=32)
            f=mp(candidate['lam']);history.append(dict(u=str(u),lam=str(f),residual=candidate['residual']))
            if not candidate['ok']:raise ArithmeticError('root slice solve failed')
            if abs(f)<mp('1e-51'):break
            if fl*f>0:lo,fl=u,f
            else:hi,fh=u,f
        else:raise ArithmeticError('root bracket exhausted')
        Z=[[mp(v) for v in z] for z in candidate['nodes']]
        diag=isolated_diagnostics(numerical_map(dt.PRIMARY),numerical_map(dt.CHECK),Z)
        rec=dict(j=j,bindings=fingerprint,u=str(u),nodes=candidate['nodes'],lam=candidate['lam'],history=history,
                 diagnostics=diag,seconds=time.time()-start)
    write(path,rec);return rec


def equivalence_classes(rows, include_half=False):
    if not rows:return dict(classes=[],max_match_error=None,min_distinct_separation=None)
    nodes=np.array([r['nodes'][0] for r in rows],float);parents=list(range(len(rows)));worst=0.
    def find(i):
        while parents[i]!=i:i=parents[i]
        return i
    for i,r in enumerate(rows):
        images=[ds.spatial(nodes[i]),np.array(r['nodes'][1],float)]
        if include_half:images.append(np.array(dt.clock_map(r['nodes'][0],half=True,**dt.PRIMARY)['z'],float))
        for z in images:
            distances=np.max(abs(nodes-z),axis=1);j=int(np.argmin(distances));worst=max(worst,float(distances[j]))
            if distances[j]<=1e-8:parents[find(i)]=find(j)
    groups={}
    for i in range(len(rows)):groups.setdefault(find(i),[]).append(rows[i]['j'])
    separation=min(float(max(abs(nodes[i]-nodes[j]))) for i in range(len(rows)) for j in range(i)) if len(rows)>1 else 0.
    return dict(classes=list(groups.values()),max_match_error=worst,min_distinct_separation=separation)


def roots_score(rows,scan,groups):
    if scan['label']!='RESOLVED_OBSTRUCTION':return dict(label=scan['label'],roots=0)
    complete=[r['j'] for r in rows]==scan['brackets']
    distinct=bool(rows and groups['P_C']['min_distinct_separation']>1e-6)
    resolved=complete and distinct and all(r['diagnostics']['ok'] for r in rows)
    return dict(label='RESOLVED_BREAKING' if resolved else 'RESOLVED_OBSTRUCTION_ROOTS_UNRESOLVED',
                roots=len(rows),all_isolated=resolved,leading_twelfth_pattern=bool(resolved and len(rows)==24 and not scan['ambiguous_phases']),
                equivalence=groups)


def safe_phase(j):
    try:return phase(j)
    except Exception as exc:
        rec=dict(j=j,bindings=bindings(),error=type(exc).__name__+': '+str(exc))
        write(RUN/'checkpoints'/f'phase_{j:03d}.json',rec);return rec


def safe_root(job):
    try:return root(job)
    except Exception as exc:
        rec=dict(j=job[0],bindings=bindings(),error=type(exc).__name__+': '+str(exc),diagnostics=dict(ok=False))
        write(RUN/'checkpoints'/f'root_{job[0]:03d}.json',rec);return rec


def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['scan','roots']);parser.add_argument('--workers',type=int,default=8);args=parser.parse_args()
    fingerprint=bindings()
    if args.stage=='scan':
        rows=[]
        with Pool(args.workers) as pool:
            for r in pool.imap_unordered(safe_phase,range(N)):
                rows.append(r);print('phase complete',r['j'],'count',len(rows),flush=True)
        rows.sort(key=lambda r:r['j']);controls=json.loads((RUN/'controls.json').read_text())
        archive=dict(bindings=fingerprint,environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,gmpy2=gp.version()),rows=rows)
        write(RUN/'scan.json',archive);score=scan_score(rows,controls);write(RUN/'scan_result.json',score)
        print(json.dumps({k:v for k,v in score.items() if k!='raw_spectrum'},indent=2))
    else:
        scan=json.loads((RUN/'scan.json').read_text());score=json.loads((RUN/'scan_result.json').read_text())
        if scan['bindings']!=fingerprint:raise ValueError('scan provenance differs')
        jobs=[(j,scan['rows'][j],scan['rows'][(j+1)%N]) for j in score['brackets']] if score['label']=='RESOLVED_OBSTRUCTION' else []
        rows=[]
        with Pool(args.workers) as pool:
            for r in pool.imap_unordered(safe_root,jobs):rows.append(r);print('root complete',r['j'],flush=True)
        rows.sort(key=lambda r:r['j'])
        groups=dict(P_C=equivalence_classes(rows),H_C=equivalence_classes(rows,True)) if all('nodes' in r for r in rows) else dict(P_C=dict(min_distinct_separation=0),H_C=dict(error='root failures'))
        write(RUN/'roots.json',dict(bindings=fingerprint,rows=rows,equivalence=groups));result=roots_score(rows,score,groups);write(RUN/'result.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
