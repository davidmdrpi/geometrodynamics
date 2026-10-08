"""Produce/replay the fixed nonlinear intervention-budget experiment."""
import argparse
import base64
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves import r3_nonlinear_budget as nb
from geometrodynamics.waves import nonlinear_supported_tt as d
from experiments.closure_ledger import r3_family_probe as old
from experiments.closure_ledger import r3_family_replay as old_replay
from experiments.closure_ledger import r3_normal_quotient_replay as nq_replay
from experiments.closure_ledger import r3_normal_quotient_probe as nq

ROOT=old.ROOT
RUN=ROOT/'experiments/closure_ledger/runs/20261005_r3_nonlinear_budget'
FREEZE='75bed15f1afc7c932462d2721c5e76d3c2546458'
SOURCES=tuple(sorted(set(nq.SOURCES+(
    'experiments/closure_ledger/r3_normal_quotient_replay.py',
    'geometrodynamics/waves/r3_nonlinear_budget.py',
    'experiments/closure_ledger/r3_nonlinear_budget_probe.py',
    'docs/r3_nonlinear_budget_prereg.md'))))


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def sources(): return {p:digest(ROOT/p) for p in SOURCES}


def family():
    for name, sha in old_replay.SHA256.items():
        if digest(old.RUN_DIR/name)!=sha: raise ValueError('parent evidence changed')
    nq_replay.replay()
    F=json.loads((old.RUN_DIR/'stage_F.json').read_text())
    S=json.loads((old.RUN_DIR/'stage_S.json').read_text())
    return nb.Family(F,S)


def save(path, obj):
    if path.exists(): raise FileExistsError(path)
    text=json.dumps(obj,allow_nan=False,separators=(',',':')).encode()
    path.write_text(base64.b64encode(gzip.compress(text,mtime=0)).decode()+'\n')


def read(path):
    return json.loads(gzip.decompress(base64.b64decode(path.read_text(),validate=False)))


def arr(x,shape):
    a=np.asarray(x,float)
    if a.shape!=shape or not np.isfinite(a).all(): raise ValueError('bad state/shape')
    return a


def close(a,b,tol=1e-9):
    if np.max(abs(np.asarray(a)-np.asarray(b)))>tol: raise ValueError('state/diagnostic mismatch')


def validate_run(f,ini,run):
    if run['arm'] not in nb.ARMS or run['method'] not in nb.METHODS: raise ValueError('unknown run')
    target=ini['target'];d0=ini['d0'];y=np.array(ini['initial']);phase=ini['initial_phase']
    drift=0.;cumulative=0.;terminal='CONTINUE';max_distance=0.;peak_cost=0.;max_constraint=0.
    if not 1<=len(run['steps'])<=target: raise ValueError('empty or oversized run')
    for n,rec in enumerate(run['steps'],1):
        if terminal!='CONTINUE' or rec['step']!=n: raise ValueError('invalid step schedule')
        close(arr(rec['start'],(29,)),y)
        if 'error' in rec:
            # A failed solver/phase fit is conservatively unresolved, never positive.
            terminal='NUMERICALLY_UNRESOLVED'
            continue
        if len(rec['returns'])!=2: raise ValueError('missing return')
        for h in rec['returns']:
            ts=arr(h['times'],(33,));ys=arr(h['states'],(33,29))
            if ts[0]!=0 or np.any(np.diff(ts)<=0): raise ValueError('invalid sample times')
            close(ys[0],y)
            close(ys[:,2],y[2]+ts)
            for row in ys:
                d.ingredients(row)
                max_constraint=max(max_constraint,float(np.max(abs(d.constraints(row)['residual']))))
            y=ys[-1]
            if abs(y[3])>1e-9 or y[7]>=0: raise ValueError('not a descending section')
        close(arr(rec['pre'],(29,)),y)
        t,dist=f.nearest(nb.coords(y));inc=float(np.angle(np.exp(1j*(t-phase))))
        if abs(inc)>=np.pi/2: raise ValueError('phase unwrap ambiguous')
        drift+=inc;phase=t
        terminal=nb.decision(dist,d0,0.,cumulative,False,n,drift,target)
        proposal=y.copy();cost=0.;postdist=dist
        if terminal!='TUBE_ESCAPE' and run['arm']=='controlled':
            z,condition=f.tune(nb.coords(y)[:6],t);proposal=nb.state(z,eta=y[2])
            cost=float(np.linalg.norm(nb.coords(proposal)-nb.coords(y)));_,postdist=f.nearest(nb.coords(proposal))
            terminal=nb.decision(max(dist,postdist),d0,cost,cumulative,True,n,drift,target)
            close(rec['controller_condition'],condition)
        applied=run['arm']=='controlled' and terminal not in ('TUBE_ESCAPE','CONTROL_BUDGET_EXCEEDED')
        close(arr(rec['proposed'],(29,)),proposal)
        if rec['applied']!=applied: raise ValueError('incorrect application of reset')
        for key,value in dict(phase_drift=drift,pre_distance=dist,proposed_cost=cost,proposed_distance=postdist).items(): close(rec[key],value)
        if abs(np.angle(np.exp(1j*(rec['phase']-phase))))>1e-8: raise ValueError('wrong fitted phase')
        max_constraint=max(max_constraint,float(np.max(abs(d.constraints(proposal)['residual']))))
        if applied: y=proposal;cumulative+=cost
        close(arr(rec['after'],(29,)),y);close(rec['cumulative_cost'],cumulative)
        max_distance=max(max_distance,dist,postdist);peak_cost=max(peak_cost,cost)
    if terminal=='CONTINUE': raise ValueError('truncated run before horizon/failure')
    if max_constraint>1e-8: terminal='NUMERICALLY_UNRESOLVED'
    if terminal!=run['terminal']: raise ValueError('terminal decision differs')
    return dict(terminal=terminal,steps=len(run['steps']),target=target,phase_drift=drift,
                max_distance_over_d0=max_distance/d0,peak_proposed_cost_over_d0=peak_cost/d0,
                applied_cumulative_cost_over_d0=cumulative/d0,max_constraint=max_constraint)


def assess(f,raw):
    if raw['freeze']!=FREEZE or raw['sources']!=sources(): raise ValueError('source provenance differs')
    if raw['parent_manifest']!=nq_replay.MANIFEST_SHA: raise ValueError('wrong parent')
    if len(raw['cases'])!=4: raise ValueError('case count')
    summaries=[]
    for case,delta in zip(raw['cases'],nb.DETUNINGS):
        if case['detuning']!=delta: raise ValueError('case schedule')
        if 'setup_error' in case:
            summaries.append(dict(detuning=delta,controlled='INVALID_SETUP',unforced='INVALID_SETUP'));continue
        ini=f.initial(delta)
        for key in ('scale','achieved_action_offset','initial','d0','target','loop','tuning_vector','tuning_norm','log10_autonomous_unstable_tolerance'):
            close(case['initial'][key],ini[key])
        runs=case['runs'];schedule=[(m,a) for m in nb.METHODS for a in nb.ARMS]
        if [(r['method'],r['arm']) for r in runs]!=schedule: raise ValueError('run schedule')
        decisions=[validate_run(f,ini,r) for r in runs];statuses={};agreements={}
        for arm in nb.ARMS:
            idx=[i for i,r in enumerate(runs) if r['arm']==arm];a,b=(runs[i] for i in idx);ra,rb=(decisions[i] for i in idx)
            agree=ra['terminal']==rb['terminal'] and ra['steps']==rb['steps'] and ra['terminal']!='NUMERICALLY_UNRESOLVED'
            discrepancy=0.
            if agree:
                for sa,sb in zip(a['steps'],b['steps']):
                    for key in ('pre','after','proposed'):
                        discrepancy=max(discrepancy,float(np.max(abs(np.asarray(sa[key])-np.asarray(sb[key])))))
                    for key in ('pre_distance','proposed_distance','proposed_cost'):
                        discrepancy=max(discrepancy,abs(sa[key]-sb[key]))
                agree=discrepancy<=1e-4*ini['d0']
            statuses[arm]=ra['terminal'] if agree else 'NUMERICALLY_UNRESOLVED'
            agreements[arm]=dict(agrees=bool(agree),max_difference=discrepancy,limit=1e-4*ini['d0'])
        summaries.append(dict(detuning=delta,d0=ini['d0'],target=ini['target'],tuning_norm=ini['tuning_norm'],
                              precision_log10=ini['log10_autonomous_unstable_tolerance'],
                              **statuses,integrator_comparison=agreements,runs=decisions))
    controlled=[s['controlled'] for s in summaries]
    if any(s in ('TUBE_ESCAPE','CONTROL_BUDGET_EXCEEDED') for s in controlled): label='CONTROLLED_NONLINEAR_BOUND_FAILED'
    elif all(s=='HORIZON_COMPLETED' for s in controlled): label='CONTROLLED_NONLINEAR_BOUND_SUPPORTED_ON_TESTED_HORIZON'
    else: label='NONLINEAR_STUDY_UNRESOLVED'
    return dict(label=label,twist=f.twist,baseline_loop_action=f.action,cases=summaries,
                autonomous_stability='NOT_ESTABLISHED',action_selection='NOT_ESTABLISHED')


def produce():
    if RUN.exists(): raise FileExistsError('append-only output exists')
    f=family();src=sources();raw=dict(freeze=FREEZE,sources=src,parent_manifest=nq_replay.MANIFEST_SHA,cases=[])
    # Durable per-case checkpoints; never overwrite a completed trajectory.
    RUN.mkdir(parents=True)
    for k,delta in enumerate(nb.DETUNINGS):
        case=dict(detuning=delta,runs=[])
        try: case['initial']=f.initial(delta)
        except (ArithmeticError,ValueError,np.linalg.LinAlgError) as e: case['setup_error']=str(e)
        if 'initial' in case:
            for method in nb.METHODS:
                for arm in nb.ARMS:
                    r=nb.simulate(f,case['initial'],method,arm);case['runs'].append(r)
                    print(delta,method,arm,r['terminal'],len(r['steps']),flush=True)
        save(RUN/f'case_{k}.json.gz.b64',case);raw['cases'].append(case)
    if sources()!=src: raise ValueError('sources changed during experiment')
    result=assess(f,raw)
    meta={k:v for k,v in raw.items() if k!='cases'}
    (RUN/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    (RUN/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    names=[f'case_{k}.json.gz.b64' for k in range(4)]+['provenance.json','result.json']
    (RUN/'manifest.json').write_text(json.dumps({n:digest(RUN/n) for n in names},indent=2)+'\n')
    print(result['label'],flush=True)


def replay(directory=RUN,manifest_sha=None):
    if not manifest_sha or digest(directory/'manifest.json')!=manifest_sha: raise ValueError('manifest fingerprint mismatch')
    manifest=json.loads((directory/'manifest.json').read_text())
    names=[f'case_{k}.json.gz.b64' for k in range(4)]+['provenance.json','result.json']
    if set(manifest)!=set(names): raise ValueError('archive inventory')
    for name in names:
        if digest(directory/name)!=manifest[name]: raise ValueError('evidence fingerprint mismatch')
    raw=json.loads((directory/'provenance.json').read_text());raw['cases']=[read(directory/f'case_{k}.json.gz.b64') for k in range(4)]
    fresh=assess(family(),raw);saved=json.loads((directory/'result.json').read_text())
    if fresh['label']!=saved['label'] or [(x['controlled'],x['unforced']) for x in fresh['cases']]!=[(x['controlled'],x['unforced']) for x in saved['cases']]: raise ValueError('decisions changed')
    return fresh


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('mode',choices=['run','replay']);ap.add_argument('--manifest-sha');a=ap.parse_args()
    if a.mode=='run': produce()
    else: print(json.dumps(replay(manifest_sha=a.manifest_sha),indent=2))
