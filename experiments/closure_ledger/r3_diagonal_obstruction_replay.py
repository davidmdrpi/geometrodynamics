"""Authenticate and rescore the frozen high-precision diagonal experiment.

--full independently reevaluates every check solution with the momentum
formulation at 336 bits/order 64, then checks every saved root's closure.
No result or historical input is modified.
"""
import argparse
import hashlib
import json
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import gmpy2 as gp
from gmpy2 import mpfr as mp
from experiments.closure_ledger import r3_diagonal_obstruction_probe as probe
from experiments.closure_ledger import r3_diagonal_obstruction_supplement as supplement
from geometrodynamics.waves import r3_high_precision_shooting as sh
from geometrodynamics.waves import r3_diagonal_taylor as dt

FREEZE='e95f8582f4d1537ada9765c620cae79238565469'
SUPPLEMENT_FREEZE='ce65db88e2ae329643659b4d7242f3ea65b72dda'
MANIFEST_SHA256='2cf5f94f7e93999df53480a760890500760ec0e71ff8f71dbbaf786880f2ee57'
FILES=('controls.json','scan.json','scan_result.json','roots.json','result.json','seed_audit.json','supplement.json')


def authenticate(directory):
    directory=Path(directory);raw=(directory/'manifest.json').read_bytes()
    if hashlib.sha256(raw).hexdigest()!=MANIFEST_SHA256:raise ValueError('manifest fingerprint mismatch')
    manifest=json.loads(raw)
    if manifest['freeze_commit']!=FREEZE or set(manifest['sha256'])!=set(FILES):raise ValueError('manifest inventory mismatch')
    if manifest['supplement_freeze_commit']!=SUPPLEMENT_FREEZE:raise ValueError('supplement freeze differs')
    out={}
    for name in FILES:
        raw=(directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=manifest['sha256'][name]:raise ValueError('archive fingerprint mismatch: '+name)
        out[name]=json.loads(raw)
    expected=probe.bindings()
    if out['scan.json']['bindings']!=expected or out['roots.json']['bindings']!=expected:raise ValueError('source/input fingerprint mismatch')
    if out['supplement.json']['inputs']!=supplement.inputs():raise ValueError('supplement input fingerprint mismatch')
    reconstructed={}
    for r in out['scan.json']['rows']:
        reconstructed[f"phase_{r['j']:03d}.json"]=hashlib.sha256((json.dumps(r,indent=1,allow_nan=False)+'\n').encode()).hexdigest()
        if 'primary' in r:
            pre={k:v for k,v in r.items() if k not in ('check','check_constraint','seconds','error_estimate')}
            reconstructed[f"primary_{r['j']:03d}.json"]=hashlib.sha256((json.dumps(pre,indent=1,allow_nan=False)+'\n').encode()).hexdigest()
    for r in out['roots.json']['rows']+out['supplement.json']['rows']:
        reconstructed[f"root_{r['j']:03d}.json"]=hashlib.sha256((json.dumps(r,indent=1,allow_nan=False)+'\n').encode()).hexdigest()
    if reconstructed!=manifest['checkpoint_sha256']:raise ValueError('checkpoint reconstruction differs')
    return out


def finite_array(a,shape):
    if np.asarray(a).shape!=shape:raise ValueError('array shape differs')
    if not all(gp.is_finite(mp(v)) for v in np.asarray(a).ravel()):raise ValueError('nonfinite data')


def validate_scan(rows):
    if len(rows)!=probe.N or [r['j'] for r in rows]!=list(range(probe.N)):raise ValueError('phase grid differs')
    expected=probe.bindings()
    for r in rows:
        if r['bindings']!=expected:raise ValueError('phase fingerprint differs')
        if 'error' in r:continue
        if abs(r['phi']-2*np.pi*r['j']/probe.N)>1e-15:raise ValueError('phase value differs')
        p=r['point'];finite_array(p['nodes'],(2,6));finite_array(p['J'],(13,13))
        for key in ('c0','t0','g'):finite_array(p[key],(6,))
        # Original deterministic slices, not redefined phases after measurement.
        _,c,t=probe.control_seeds(r['phi']);t=t/np.linalg.norm(t)
        if np.max(abs(np.array(p['c0'])-c))>1e-14 or np.max(abs(np.array(p['t0'])-t))>1e-14:raise ValueError('slice changed')
        if np.max(abs(np.array(p['g'])-probe.b.symplectic(6).T@t))>1e-14:raise ValueError('unfolding changed')
        for key,tol in (('primary','1e-45'),('check','1e-65')):
            a=r[key];finite_array(a['nodes'],(2,6));finite_array([a['lam'],a['residual']],(2,))
            if a['ok']!=(mp(a['residual'])<=mp(tol)):raise ValueError('convergence flag differs')
        with gp.context(gp.get_context(),precision=336):
            error=abs(mp(r['primary']['lam'])-mp(r['check']['lam']))+mp(2*max(r['primary']['inverse_row_norm'],r['check']['inverse_row_norm']))*(mp(r['primary']['residual'])+mp(r['check']['residual']))
            if error!=mp(r['error_estimate']):raise ValueError('error estimate differs')


def compare(a,b,geometry=False):
    if isinstance(a,dict):
        if not isinstance(b,dict) or a.keys()!=b.keys():raise ValueError('score keys differ')
        for k in a:compare(a[k],b[k],geometry)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise ValueError('score length differs')
        for x,y in zip(a,b):compare(x,y,geometry)
    elif isinstance(a,float):
        if not np.isfinite(b) or not np.isclose(a,b,rtol=2e-12,atol=1e-14 if geometry else 1e-100):raise ValueError('numeric score differs')
    elif a!=b:raise ValueError('score differs')


def reevaluate(row):
    cfg=dict(dt.CHECK,formulation='momentum')
    with gp.context(gp.get_context(),precision=336):
        Z=[[mp(v) for v in z] for z in row['check']['nodes']]
        p=row['point'];g,t,c=([mp(v) for v in p[k]] for k in ('g','t0','c0'))
        R=sh.residual(probe.numerical_map(cfg),Z,mp(row['check']['lam']),g,t,c)
        err=max(map(abs,R))
        if err>mp('1e-45'):raise ValueError('high-precision augmented replay failed')
        return dict(j=row['j'],residual=str(err))


def root_closure(row):
    with gp.context(gp.get_context(),precision=336):
        Z=[[mp(v) for v in z] for z in row['nodes']]
        err=max(map(abs,sh.closure(probe.numerical_map(dict(dt.CHECK,formulation='momentum')),Z)))
        if err>mp('1e-45'):raise ValueError('root closure replay failed')
        return dict(j=row['j'],residual=str(err))


def validate_roots(rr,brackets):
    if [r['j'] for r in rr]!=brackets:raise ValueError('root grid differs')
    for r in rr:
        if r['bindings']!=probe.bindings():raise ValueError('root source binding differs')
        if 'error' in r:continue
        finite_array(r['nodes'],(2,6))
        finite_array([r['lam'],r['u']],(2,))
        if not (mp(0)<=mp(r['u'])<=mp(1)) or abs(mp(r['lam']))>=mp('1e-51'):raise ValueError('root bracket/residual differs')
        # The producer reparses the 224-bit solution into its 336-bit bracket
        # context before logging f, so equal values have different strings.
        with gp.context(gp.get_context(),precision=224):
            if not r['history'] or mp(r['history'][-1]['lam'])!=mp(r['lam']) or r['history'][-1]['u']!=r['u']:raise ValueError('root history differs')
        for step in r['history']:
            finite_array([step[k] for k in ('lam','u','residual')],(3,))
            if mp(step['residual'])>mp('1e-56'):raise ValueError('root solve residual differs')
        finite_array([v for k,v in r['diagnostics'].items() if k!='ok'],(6,))
        d=r['diagnostics'];ok=mp(d['inverse_times_difference'])<mp('.01') and mp(d['closure_residual'])<mp('1e-45') and mp(d['correction_norm'])<mp('1e-12') and mp(d['inverse_residual'])<mp('1e-40')
        if d['ok']!=ok:raise ValueError('isolation flag differs')


def replay(directory=None,full=False):
    rec=authenticate(directory or probe.RUN);rows=rec['scan.json']['rows'];validate_scan(rows)
    score=probe.scan_score(rows,rec['controls.json']);compare(score,rec['scan_result.json'])
    roots=rec['roots.json'];rr=roots['rows']
    validate_roots(rr,score['brackets'] if score['label']=='RESOLVED_OBSTRUCTION' else [])
    compare(probe.roots_score(rr,score,roots['equivalence']),rec['result.json'])
    sr=rec['supplement.json'];extra=sr['rows'];validate_roots(extra,score['brackets'])
    if sr['registered_label']!=score['label']:raise ValueError('registered verdict changed')
    compare(supplement.audit(rows),rec['seed_audit.json'],geometry=True)
    compare(supplement.score(extra,score),sr['result'],geometry=True)
    result=dict(replay='VERIFIED',label=rec['result.json']['label'],phases=len(rows),roots=len(rr),
                supplemental_label=sr['result']['label'],supplemental_roots=len(extra))
    rr=rr+extra
    if full:
        if any('error' in r for r in rows+rr):raise ValueError('cannot replay failed numerical jobs')
        with Pool(8) as pool:
            phases=pool.map(reevaluate,rows);root_rows=pool.map(root_closure,rr)
        result.update(phase_residuals=phases,root_residuals=root_rows,
                      max_phase_residual=str(max(mp(r['residual']) for r in phases)),
                      max_root_residual=str(max((mp(r['residual']) for r in root_rows),default=mp(0))))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true');parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=replay(full=args.full)
    if args.output:args.output.write_text(json.dumps(result,indent=1)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('phase_residuals','root_residuals')},indent=2))
