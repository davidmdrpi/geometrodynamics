"""Execute and score the frozen compact-handle constraint experiment."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves import mouth_momentum as m

GATES = ('symbolic', 'momentum_completion', 'tensor_descent', 'hamiltonian',
         'physical_constraints', 'minimal_section', 'controls', 'evidence')
UNESTABLISHED = ('round_s3_embedding', 'four_scalar_interface', 'traversability',
                 'mouth_evolution', 'radiative_momentum_transfer', 'crossing_events',
                 'discrete_action', 'quantum_statistics')


def maximum(x):
    a = np.asarray(x, dtype=float)
    if a.size == 0 or not np.isfinite(a).all():
        raise ValueError('empty or nonfinite evidence')
    return float(np.max(abs(a)))


def all_finite(obj):
    if isinstance(obj, dict):
        return all(all_finite(x) for x in obj.values())
    if isinstance(obj, (list, tuple)):
        return all(all_finite(x) for x in obj)
    if isinstance(obj, (int, float)):
        return bool(np.isfinite(obj))
    return True


def run():
    data = dict(prereg=m.PREREG, baseline=m.BASELINE, gate_schema=list(GATES),
                symbolic=list(m.symbolic_checks()), momentum_fd=[], solutions=[],
                physical=[], controls={}, failures=[])
    for n in (32,64,128):
        data['momentum_fd'].append(m.momentum_fd(n))
    s = np.linspace(0,4*np.pi,65)
    data['momentum_negative'] = dict(s=s.tolist(),
        uncorrected=(.1*m.radial(s,correction=0)[3]).tolist(),
        damaged=(.1*m.radial(s,correction=1.1)[3]).tolist(),
        completed=(.1*m.radial(s)[3]).tolist())
    for epsilon in m.AMPLITUDES:
        for ns,nu in m.GRIDS:
            try:
                r = m.solve(ns,nu,epsilon)
                data['solutions'].append(r)
            except Exception as error:
                data['failures'].append(dict(epsilon=epsilon,ns=ns,nu=nu,error=str(error)))
    fine = next(r for r in data['solutions'] if r['epsilon']==.2 and r['ns']==64)
    f = m.Interpolant(fine)
    for h in (1e-3,5e-4,2.5e-4):
        data['physical'].append(dict(h=h,cases=[m.coordinate_constraints(f,x,h) for x in m.POINTS]))
    data['wrong_weight'] = [m.coordinate_constraints(f,x,2.5e-4,wrong_weight=True) for x in m.POINTS]
    data['seam'] = m.seam_check(f)
    data['bad_seam'] = m.seam_check(f,k=1)
    data['sections'] = [dict(epsilon=r['epsilon'],values=[m.Interpolant(r).section(s) for s in (0,.05,.1,float(np.pi))])
                        for r in data['solutions'] if r['ns']==64 and r['epsilon']>0]
    for name,epsilon,C,k in [('amplitude_reversal',-.1,.5,.5),('time_reversal',-.1,-.5,.5),('untwisted',.1,.5,1.)]:
        data['controls'][name]=m.solve(40,20,epsilon,C=C,k=k)
    data['untwisted_seam']=m.seam_check(m.Interpolant(data['controls']['untwisted']),twisted=False)
    return data


def score(data):
    gates = dict.fromkeys(GATES,False)
    metrics = {}
    try:
        integrity = (all_finite(data) and data['gate_schema']==list(GATES)
                     and data['prereg']==m.PREREG and data['baseline']==m.BASELINE
                     and data['failures']==[])
        gates['symbolic'] = data['symbolic']==list(m.symbolic_checks()) and all(x=='0' for x in data['symbolic'])
        fds=data['momentum_fd']
        integrity &= [r['n'] for r in fds]==[32,64,128]
        errors=[]
        for r in fds:
            n=r['n'];s=4*np.pi*np.arange(n)/n;h=4*np.pi/n
            w=np.asarray(r['w'])
            residual=(np.roll(w,1)-2*w+np.roll(w,-1))/h**2-4*w+kcos(s)
            integrity &= maximum(np.asarray(r['s'])-s)<1e-13 and maximum(residual)<1e-10
            integrity &= maximum(np.asarray(r['rhs'])+kcos(s))<1e-13
            integrity &= maximum(np.asarray(r['exact'])-m.radial(s)[1])<1e-13
            integrity &= maximum(np.asarray(r['residual'])-residual)<1e-12
            errors.append(maximum(w-m.radial(s)[1]))
        ratios=np.array(errors[:-1])/errors[1:]
        neg=data['momentum_negative'];s=np.asarray(neg['s'])
        for name,correction in [('uncorrected',0),('damaged',1.1),('completed',1)]:
            integrity &= maximum(np.asarray(neg[name])-.1*m.radial(s,correction=correction)[3])<1e-12
        gates['momentum_completion'] = (errors[-1]<1e-3 and np.all((ratios>=3.8)&(ratios<=4.2))
            and maximum(neg['completed'])<1e-12 and maximum(neg['uncorrected'])>1e-4 and maximum(neg['damaged'])>1e-4)
        metrics.update(momentum_errors=errors,momentum_ratios=ratios.tolist())

        records=data['solutions']
        keys=[(r['epsilon'],r['ns'],r['nu']) for r in records]
        expected=[(e,ns,nu) for e in m.AMPLITUDES for ns,nu in m.GRIDS]
        integrity &= keys==expected
        lookup={key:r for key,r in zip(keys,records)}
        ongrid=[];offgrid=[];differences=[];minima=[];norms=[]
        for r in records:
            integrity &= r['C']==.5 and r['k']==.5 and r['Lambda']==.25
            f=m.Interpolant(r)
            minima.append(float(f.psi.min()))
            norms.append(float(np.sqrt(np.max(m.norm_squared(f.grid.s[:,None],f.grid.u[None,:],r['epsilon'])*f.psi**-12))))
            ongrid.append(maximum(m.equation(f.grid,f.psi,r['epsilon'])))
            integrity &= r['start']=='psi=1' and 0<=r['iterations']<30
            integrity &= len(r['residual_history'])==r['iterations']+1
            integrity &= abs(r['residual_history'][-1]-ongrid[-1])<1e-12
            if r['ns']==64:offgrid.append(maximum(f.offgrid_residual()))
        target=m.Grid(64,32)
        for e in m.AMPLITUDES:
            a=[m.Interpolant(lookup[(e,)+grid]).values(target.s,target.u) for grid in m.GRIDS]
            differences.append([maximum(a[0]-a[1]),maximum(a[1]-a[2])])
        gates['hamiltonian'] = (keys==expected and max(ongrid)<1e-10 and max(offgrid)<1e-8
                               and max(x[1] for x in differences)<1e-7 and min(minima)>0)
        metrics.update(ongrid_max=max(ongrid),offgrid_max=max(offgrid),grid_differences=differences,
                       minimum_psi=min(minima),maximum_K_norm_on_frozen_grids=max(norms))

        fine=m.Interpolant(lookup[.2,64,32])
        for key,args in [('seam',{}),('bad_seam',{'k':1})]:
            recomputed=m.seam_check(fine,**args)
            for field,values in recomputed.items():integrity &= maximum(np.asarray(data[key][field])-values)<1e-12
        seam=data['seam'];bad=data['bad_seam']
        gates['tensor_descent']=(max(maximum(seam[x]) for x in seam)<1e-10
                                and maximum(bad['value_residuals'])>1e-3)
        metrics.update(seam_max=max(maximum(seam[x]) for x in seam),bad_seam_max=maximum(bad['value_residuals']))

        # Reconstruct constraints from archived physical jets' contractions;
        # do not trust saved normalized errors. K and g must match the solution.
        Herrors=[];Merrors=[]
        def physical_error(case,wrong=False):
            nonlocal integrity
            g,K=fine.metric_tensor(case['point'],wrong_weight=wrong)
            integrity &= maximum(np.asarray(case['metric'])-g)<1e-12
            integrity &= maximum(np.asarray(case['extrinsic_curvature'])-K)<1e-12
            mixed=np.linalg.solve(g,K);K2=float(np.trace(mixed@mixed));tr=float(np.trace(mixed))
            integrity &= abs(K2-case['K2'])<1e-12 and abs(tr-case['trace_K'])<1e-12
            H=case['R']+tr*tr-K2-.5
            terms=np.asarray(case['momentum_terms']);M=terms.sum(axis=0)
            integrity &= abs(H-case['hamiltonian'])<1e-12 and maximum(M-case['momentum'])<1e-12
            return abs(H)/max(1.,abs(case['R'])+tr*tr+K2+.5),float(np.linalg.norm(M)/max(1.,sum(np.linalg.norm(t) for t in terms)))
        integrity &= [r['h'] for r in data['physical']]==[1e-3,5e-4,2.5e-4]
        for group in data['physical']:
            integrity &= [tuple(c['point']) for c in group['cases']]==list(m.POINTS)
            integrity &= all(c['h']==group['h'] for c in group['cases'])
            err=[physical_error(c) for c in group['cases']]
            Herrors.append(max(x[0] for x in err));Merrors.append(max(x[1] for x in err))
        integrity &= [tuple(c['point']) for c in data['wrong_weight']]==list(m.POINTS)
        wrong=[physical_error(c,True) for c in data['wrong_weight']]
        wrongH=max(x[0] for x in wrong);wrongM=max(x[1] for x in wrong)
        ratios_ok=all(v[-2]<=1e-8 or 2.5<=v[-2]/v[-1]<=5.5 for v in (Herrors,Merrors))
        gates['physical_constraints']=(Herrors[-1]<1e-5 and Merrors[-1]<1e-5 and ratios_ok
                                      and wrongH>10*Herrors[-1] and wrongM>10*Merrors[-1])
        metrics.update(physical_H=Herrors,physical_M=Merrors,wrong_weight_H=wrongH,wrong_weight_M=wrongM)
        metrics['physical_sample_max_K_norm']=float(np.sqrt(max(c['K2'] for c in data['physical'][-1]['cases'])))

        neck_ok=True;area_data=[]
        integrity &= [r['epsilon'] for r in data['sections']]==list(m.AMPLITUDES[1:])
        for row in data['sections']:
            ff=m.Interpolant(lookup[row['epsilon'],64,32])
            actual=[ff.section(s) for s in (0,.05,.1,float(np.pi))]
            integrity &= len(row['values'])==4
            for stored,computed in zip(row['values'],actual):
                for key in computed:integrity &= maximum(np.asarray(stored[key])-computed[key])<1e-12
            area=[x['area'] for x in actual]
            neck_ok &= maximum(actual[0]['H'])<1e-8 and all(x>area[0] for x in area[1:])
            area_data.append(dict(epsilon=row['epsilon'],areas=area,
                                 expansions_min=float(min(actual[0]['theta_plus']+actual[0]['theta_minus'])),
                                 expansions_max=float(max(actual[0]['theta_plus']+actual[0]['theta_minus']))))
        gates['minimal_section']=neck_ok
        metrics['sections']=area_data

        reference=m.Interpolant(lookup[.1,40,20]);controls=data['controls']
        integrity &= set(controls)=={'amplitude_reversal','time_reversal','untwisted'}
        control_errors={}
        expected_params={'amplitude_reversal':(-.1,.5,.5),'time_reversal':(-.1,-.5,.5),'untwisted':(.1,.5,1.)}
        control_ok=True
        for name,r in controls.items():
            integrity &= (r['epsilon'],r['C'],r['k'])==expected_params[name] and (r['ns'],r['nu'])==(40,20) and r['Lambda']==.25
            f=m.Interpolant(r)
            residual=maximum(m.equation(f.grid,f.psi,r['epsilon'],r['C'],r['k']))
            control_ok &= residual<1e-10
            control_errors[name]=dict(pde=residual)
            if name!='untwisted':
                delta=maximum(f.psi-reference.psi)
                control_ok &= delta<1e-9
                control_errors[name]['psi_difference']=delta
                for x in m.POINTS:
                    K=f.metric_tensor(x)[1];K0=reference.metric_tensor(x)[1]
                    mismatch=maximum(K+K0) if name=='time_reversal' else max(abs(K[0,2]+K0[0,2]),abs(K[1,2]+K0[1,2]))
                    control_ok &= mismatch<1e-9
                if name=='time_reversal':
                    for s in (0.,.05,.1,float(np.pi)):
                        reverse=f.section(s);forward=reference.section(s)
                        # Reversing K exchanges outgoing/ingoing expansions
                        # with a minus sign (at H=0 both simply reverse).
                        control_ok &= maximum(np.asarray(reverse['theta_plus'])+forward['theta_minus'])<1e-9
                        control_ok &= maximum(np.asarray(reverse['theta_minus'])+forward['theta_plus'])<1e-9
        untwisted=m.seam_check(m.Interpolant(controls['untwisted']),twisted=False)
        for key in untwisted:integrity &= maximum(np.asarray(data['untwisted_seam'][key])-untwisted[key])<1e-12
        control_ok &= max(maximum(x) for x in untwisted.values())<1e-10
        # Zero-wave solution stays the exact product solution.
        control_ok &= all(maximum(np.asarray(lookup[(0.,)+grid]['psi'])-1)<1e-12 for grid in m.GRIDS)
        gates['controls']=control_ok
        metrics['controls']=control_errors
        gates['evidence']=bool(integrity)
    except (KeyError,ValueError,TypeError,IndexError,ArithmeticError) as error:
        gates['evidence']=False
        metrics['evidence_error']=str(error)
    gates={key:bool(value) for key,value in gates.items()}
    compact=all(gates[k] for k in GATES if k!='minimal_section')
    return dict(gates=gates,passed=sum(gates.values()),total=len(GATES),metrics=metrics,
                verdicts=dict(COMPACT_TWISTED_VACUUM_CONSTRAINT_DATA=compact,
                    MINIMAL_SECTION_WITH_NONTRIVIAL_MOMENTUM_COMPLETION=compact and gates['minimal_section']),
                unestablished=dict.fromkeys(UNESTABLISHED,True))


def kcos(s):
    return .5*np.cos(.5*s)


def report(result):
    lines=['# Compact twisted-handle constraint experiment','',f"Public freeze: `{m.PREREG}`.",
           '',f"Frozen gates: **{result['passed']}/{result['total']}**.",'', '| Gate | Pass |','|---|---|']
    lines.extend(f'| {name} | {value} |' for name,value in result['gates'].items())
    lines.extend(['','```json',json.dumps(result,indent=2,allow_nan=False),'```','',
                  'Initial constraints only. No traversability, crossing, momentum-transfer or quantum claim.',''])
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--rescore',type=Path)
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    # Withdraw any stale result before work or parsing can fail.
    (args.output_dir/'probe.md').write_text('# Experiment incomplete\n\nNo affirmative verdict.\n')
    try:
        if args.rescore:
            raw=args.rescore.read_bytes()
            data=json.loads(gzip.decompress(raw) if args.rescore.suffix=='.gz' else raw)
        else:
            data=run()
        result=score(data)
        raw=json.dumps(data,sort_keys=True,indent=2,allow_nan=False).encode()+b'\n'
        (args.output_dir/'probe.json.gz').write_bytes(gzip.compress(raw,mtime=0))
        (args.output_dir/'probe.md').write_text(report(result))
        (args.output_dir/'verdict.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(json.dumps(dict(passed=result['passed'],gates=result['gates'],metrics=result['metrics'],raw_sha256=hashlib.sha256(raw).hexdigest()),indent=2))
        return 0 if all(result['gates'].values()) else 1
    except Exception as error:
        (args.output_dir/'probe.md').write_text(f'# Experiment failed\n\nNo affirmative verdict.\n\n{type(error).__name__}: {error}\n')
        (args.output_dir/'verdict.json').write_text(json.dumps({'error':str(error),'verdicts':{'COMPACT_TWISTED_VACUUM_CONSTRAINT_DATA':False,'MINIMAL_SECTION_WITH_NONTRIVIAL_MOMENTUM_COMPLETION':False}})+'\n')
        raise


if __name__=='__main__':
    raise SystemExit(main())
