"""Targeted prospectively frozen refinement; original verdict and files stay intact."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson,solve_ivp
from geometrodynamics.waves import selection_plateau as p
from experiments.closure_ledger.selection_plateau_probe import ROOT,RUN

FREEZE='e1abc93f0accb1e49793fb4792f0cf5780fcb161'
PARAMETERS=[float(p.AMPLITUDES[-1]),0.,float(np.pi/2)]
TIMES=np.r_[0.,np.concatenate([t+np.arange(-10,11)*.025 for t in p.DURATIONS])]
COARSE=np.r_[0,np.concatenate([1+21*j+np.arange(0,21,2) for j in range(4)])]


def run_case(method):
    rtol,atol,step=(1e-12,1e-14,.02) if method=='DOP853' else (1e-10,1e-12,.01)
    sol=solve_ivp(p.rhs,(0.,TIMES[-1]),np.r_[p.preparation(*PARAMETERS),0.,0.],
                  t_eval=TIMES,method=method,rtol=rtol,atol=atol,max_step=step)
    if not sol.success or not np.isfinite(sol.y).all():raise ArithmeticError('refinement integration failed')
    return dict(method=method,parameters=PARAMETERS,states=sol.y.T.tolist())


def summary(row):
    states=np.array(row['states'])
    if row['parameters']!=PARAMETERS or row['method'] not in ('DOP853','RK45') or states.shape!=(85,31) or not np.isfinite(states).all():
        raise ValueError('invalid refinement data')
    np.testing.assert_allclose(states[0],np.r_[p.preparation(*PARAMETERS),0.,0.],rtol=1e-12,atol=1e-14)
    d=p.summarize_history(states[COARSE])
    action=np.array([p.quantities(y)[0] for y in states]);change=action-action[0]
    means=[simpson(change[1+21*j:1+21*(j+1)],x=TIMES[1+21*j:1+21*(j+1)],axis=0)/.5 for j in range(4)]
    old=np.array(d['window_means'])
    d['coarse_window_means']=old.tolist();d['window_means']=np.array(means).tolist()
    d['quadrature_error_over_Jbg']=float(np.max(abs(np.array(means)-old))/p.J_BG)
    d['work_error_over_Jbg']=float(np.max(abs(change-states[:,29:]))/p.J_BG)
    for y in states:
        A,_,q,_,M,L=p.dynamics.unpack(y)
        d['constraint_max']=max(d['constraint_max'],float(max(p.dynamics.constraints(y)['normalized'])))
        d['det_error']=max(d['det_error'],float(abs(np.linalg.det(M)-1)))
        d['symmetry_error']=max(d['symmetry_error'],float(max(np.max(abs(M-M.T)),np.max(abs(M @ L-(M @ L).T)))))
        d['sampled_min_A_H_M_omega2']=np.minimum(d['sampled_min_A_H_M_omega2'],[A,A*A-q@q/6,np.linalg.eigvalsh((M+M.T)/2).min(),p.quantities(y)[2]]).tolist()
    d.pop('changes')  # raw finer states reproduce them; do not retain a mislabelled coarse series
    return d


def analyze(original,rows):
    if [r['method'] for r in rows]!=['DOP853','RK45']:raise ValueError('both registered refinements required')
    replaced=[summary(r) for r in rows]
    cases=[dict(row) for row in original['cases']];controls=[dict(row) for row in original['controls']]
    differences=[]
    for group,d in zip((cases,controls),replaced):
        indices=[i for i,row in enumerate(group) if row['parameters']==PARAMETERS]
        if len(indices)!=1:raise ValueError('target must occur once')
        i=indices[0]
        differences.append(float(np.max(abs(np.array(d['window_means'])-group[i]['diagnostics']['window_means']))/p.J_BG))
        group[i]=dict(parameters=PARAMETERS,diagnostics=d)
    gates=dict(original['numerical_gates'])
    allrows=cases+controls
    gates['constraints']=all(r['diagnostics']['constraint_max']<1e-8 for r in allrows)
    gates['chart']=all(r['diagnostics']['det_error']<1e-7 and r['diagnostics']['symmetry_error']<1e-7 and min(r['diagnostics']['sampled_min_A_H_M_omega2'])>0 for r in allrows)
    gates['work_balance']=all(r['diagnostics']['work_error_over_Jbg']<1e-8 for r in allrows)
    gates['window_quadrature']=all(r['diagnostics']['quadrature_error_over_Jbg']<1e-6 for r in allrows)
    lookup={tuple(r['parameters']):r['diagnostics'] for r in cases}
    integrator=max(float(np.max(abs(np.array(lookup[tuple(r['parameters'])]['window_means'])-r['diagnostics']['window_means']))/p.J_BG) for r in controls)
    gates['independent_integrator']=integrator<1e-6
    means=np.array([row['diagnostics']['window_means'][1:] for row in cases[12:]]).reshape(8,12,3,2)
    windows=p.classify(means)
    return dict(refinement_freeze=FREEZE,original_verdict=original['selection_verdict'],
                original_numerical_gates=original['numerical_gates'],numerical_gates=gates,
                replacements=[dict(parameters=PARAMETERS,method=r['method'],diagnostics=d) for r,d in zip(rows,replaced)],
                changes_from_original_means_over_Jbg=differences,independent_integrator_error_over_Jbg=integrator,
                candidate_windows=windows,selection_verdict=p.selection_verdict(gates,windows),localized_receiver_verdict='NOT_TESTED')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay',type=Path)
    parser.add_argument('--output-dir',type=Path,default=RUN)
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    raw=json.loads(args.replay.read_text()) if args.replay else dict(times=TIMES.tolist(),rows=[run_case(m) for m in ('DOP853','RK45')])
    if not np.array_equal(raw['times'],TIMES):raise ValueError('incorrect refinement time schedule')
    statepath=args.output_dir/'refinement_states.json'
    statepath.write_text(json.dumps(raw,separators=(',',':'),allow_nan=False)+'\n')
    result=analyze(json.loads((RUN/'plateau.json').read_text()),raw['rows'])
    result['source_sha256']={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in [
        'experiments/closure_ledger/selection_plateau_refinement.py','geometrodynamics/waves/selection_plateau.py',
        'docs/nonlinear_selection_plateau_refinement_prereg.md','experiments/closure_ledger/runs/20260928_selection_plateau/plateau.json']}
    result['raw_sha256']=hashlib.sha256(statepath.read_bytes()).hexdigest()
    (args.output_dir/'refinement.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))
    return int(result['selection_verdict']=='INCONCLUSIVE_NUMERICAL_FAILURE')


if __name__=='__main__':raise SystemExit(main())
