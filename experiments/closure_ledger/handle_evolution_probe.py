"""Frozen coupled evolution with separate crossing, balance and impulse verdicts."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves import handle_evolution as e
from . import localized_mouth_error_budget_probe as parent
from .evidence_archive import read_bytes

FREEZE = 'cdb045eddfa0f01bf560559c5ea71b90887e3d0f'
ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT/'experiments/closure_ledger/runs/20260916_localized_mouth'
PRIOR = ROOT/'experiments/closure_ledger/runs/20260922_localized_mouth_error_budget/error_budget.json'
SOURCES = ('geometrodynamics/waves/handle_evolution.py',
           'experiments/closure_ledger/handle_evolution_probe.py',
           'docs/localized_mouth_evolution_prereg.md')
VERDICTS = ('EVOLVED_TEST_WORLDLINE_CROSSING', 'COVARIANT_MATTER_MOMENTUM_BALANCE',
            'FINITE_CROSSING_IMPULSE', 'DISCRETE_RECIPROCAL_MOMENTUM_EXCHANGE')
COLUMNS = ('time','H_normalized','M_normalized','H_absolute','M_absolute','min_A','min_B',
           'min_f','neck_Jordan_radius','bulk_Jordan_radius','neck_bulk_ratio','tube_momentum',
           'integrated_balance','integrated_flux_only','neck_minimum_margin','bulk_field_norm')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance():
    return dict(freeze=FREEZE, parent_raw_sha256=digest(PRIOR),
                refinement_sha256=parent.REFINEMENT_HASH,
                sources={v:digest(ROOT/v) for v in SOURCES})


def inputs():
    previous = parent.load_inputs(OLD)
    verdict = parent.score(json.loads(PRIOR.read_text()), previous)
    if verdict['passed'] != 8:
        raise ValueError('parent 8/8 certificate did not replay')
    raw = gzip.decompress(read_bytes(OLD/'refinement.json.gz'))
    if hashlib.sha256(raw).hexdigest() != parent.REFINEMENT_HASH:
        raise ValueError('initial archive mismatch')
    data = json.loads(raw)
    records = [next(r for r in data['solutions'] if (r['L'],r['eta'],r['n_initial']) == (5.5,eta,513)) for eta in (0.,.3)]
    return records, data['profiles']['5.5']


def convergence(runs):
    fields = [r['comparison_fields'][:, [0,1,4,5], :] for r in runs]
    differences = [float(np.max(abs(a-b)/np.maximum(1,abs(b)))) for a,b in zip(fields,fields[1:])]
    radii = [np.asarray(r['diagnostics'])[:,10] for r in runs]
    rd = [float(np.max(abs(a-b))) for a,b in zip(radii,radii[1:])]
    return dict(field_differences=differences, radius_ratio_differences=rd)


def run():
    records, profiles = inputs()
    checks = e.validate_equations()
    if not(checks['target_metric_compatible'] and checks['equivariant_closure']
           and checks['round_evolution_error'] < 1e-10 and checks['coordinate_geometry_error'] < 1e-10):
        raise ArithmeticError('independent equation validation failed')
    control = e.homogeneous_controls()
    if max(v['constraint_max'] for v in control) >= 1e-9:
        raise ArithmeticError('homogeneous constraint control failed')
    groups = []
    for record in records:
        runs = []
        for n in (512,1024,2048):
            print(f'Evolving eta={record["eta"]}, N={n}',flush=True)
            runs.append(e.run(record,profiles,n))
        errors = convergence(runs)
        for row in runs:
            del row['comparison_fields']
        groups.append(dict(eta=record['eta'], runs=runs, convergence=errors))
    return dict(**provenance(), columns=list(COLUMNS), validation=checks, homogeneous=control,
                groups=groups, receiver_model='test particles; no independent gravitating receiver')


def converges(values):
    return values[-1] < 1e-3 and (values[0] <= 1e-8 or (values[1] > 0 and 4 <= values[0]/values[1] <= 32))


def group_score(group):
    runs = group['runs']
    if [v['N'] for v in runs] != [512,1024,2048]:
        raise ValueError('changed evolution schedule')
    diagnostics = [np.asarray(r['diagnostics']) for r in runs]
    for r,d in zip(runs,diagnostics):
        if d.shape != (201,len(COLUMNS)) or not np.isfinite(d).all():
            raise ValueError('missing/nonfinite diagnostics')
        if r['steps'] != 200*r['N']//512 or r['dt'] != e.T/r['steps']:
            raise ValueError('changed time schedule')
        if not np.allclose(d[:,0],np.linspace(0,e.T,201),rtol=0,atol=1e-16):
            raise ValueError('changed output times')
        if np.asarray(r['final_fields']).shape != (8,r['N']):
            raise ValueError('missing final fields')
    h = [float(np.max(d[:,1])) for d in diagnostics]
    m = [float(np.max(d[:,2])) for d in diagnostics]
    finite = all(np.isfinite(np.asarray(r['final_fields'])).all() for r in runs)
    decreasing = lambda v: all(a>b or a<1e-8 for a,b in zip(v,v[1:]))
    evolution = (finite and h[-1]<1e-3 and m[-1]<1e-3 and decreasing(h) and decreasing(m)
                 and converges(group['convergence']['field_differences'])
                 and converges(group['convergence']['radius_ratio_differences'])
                 and all(np.min(d[:,5:7])>.005 and np.min(d[:,7])>.1 for d in diagnostics))
    defects, wrong = [], []
    for d in diagnostics:
        dp = d[:,11]-d[0,11]
        scale = np.maximum(1,np.maximum(abs(d[:,11]),abs(d[:,12])))
        defects.append(float(np.max(abs(dp-d[:,12])/scale)))
        wrong.append(float(np.max(abs(dp-d[:,13])/scale)))
    balance = evolution and defects[-1]<1e-5 and decreasing(defects) and wrong[-1]>max(10*defects[-1],1e-8)
    fine, medium = runs[-1]['crossings'], runs[-2]['crossings']
    crossing = (evolution and len(fine)==len(medium)==6 and all(v['crossed'] for v in fine+medium)
                and np.min(diagnostics[-1][:,14])>0)
    trajectory = np.asarray(runs[-1]['trajectories'])
    if trajectory.shape != (201,4,6) or not np.isfinite(trajectory).all():
        raise ValueError('missing/nonfinite worldlines')
    reflection = float(np.max(abs(trajectory[:,:,:3]+trajectory[:,:,3:])))
    crossing &= reflection<1e-8 and np.max(abs(trajectory[:,3]))<1
    impulses = []
    if all(v.get('crossed',False) for v in fine+medium):
        crossing &= all(abs(a['time']-b['time'])<1e-4 and abs(a['p_hat']-b['p_hat'])<1e-3 for a,b in zip(fine,medium))
        for a,b in zip(fine,medium):
            row = {}
            for name in ('p_s','p_hat'):
                jumps = [v[name] for v in a['windows']]
                ratios = [abs(x/y) if y else None for x,y in zip(jumps,jumps[1:])]
                delta = abs(jumps[-1]-b['windows'][-1][name])
                candidate = (abs(jumps[-1])>1e-6 and all(r is not None and .8<=r<=1.2 for r in ratios)
                             and delta<.1*abs(jumps[-1]))
                row[name] = dict(jumps=jumps,ratios=ratios,medium_fine_difference=delta,candidate=bool(candidate))
            impulses.append(row)
    return dict(eta=group['eta'],evolution_valid=bool(evolution),crossing=bool(crossing),balance=bool(balance),
                H_max=h,M_max=m,convergence=group['convergence'],balance_defects=defects,
                omitted_geometry_defects=wrong,reflection_error=reflection,crossings=fine,impulses=impulses,
                finite_impulse=bool(crossing and any(v['p_s']['candidate'] or v['p_hat']['candidate'] for v in impulses)),
                bulk_Jordan_change=float(diagnostics[-1][-1,9]/diagnostics[-1][0,9]-1),
                neck_Jordan_change=float(diagnostics[-1][-1,8]/diagnostics[-1][0,8]-1),
                field_norm_change=float(diagnostics[-1][-1,15]-diagnostics[-1][0,15]))


def agreement(a,b):
    if isinstance(b,dict):
        return isinstance(a,dict) and a.keys()==b.keys() and all(agreement(a[k],v) for k,v in b.items())
    if isinstance(b,list):
        return isinstance(a,list) and len(a)==len(b) and all(agreement(x,y) for x,y in zip(a,b))
    if isinstance(b,(bool,str)) or b is None:
        return type(a) is type(b) and a==b
    return isinstance(a,(int,float)) and np.isfinite(a) and abs(a-b)<=1e-11*max(1,abs(b))


def score(data, replay=True):
    evidence = all(data.get(k)==v for k,v in provenance().items()) and data.get('columns')==list(COLUMNS)
    if replay:
        evidence &= agreement(data,run())
    if [v['eta'] for v in data['groups']] != [0.,.3]:
        raise ValueError('changed amplitude schedule')
    validation = data['validation']
    evidence &= (validation['target_metric_compatible'] and validation['equivariant_closure']
                 and validation['round_evolution_error']<1e-10 and validation['coordinate_geometry_error']<1e-10
                 and max(v['constraint_max'] for v in data['homogeneous'])<1e-9)
    groups = [group_score(g) for g in data['groups']]
    crossing = bool(evidence and all(g['crossing'] for g in groups))
    balance = bool(evidence and all(g['balance'] for g in groups))
    impulse = bool(evidence and any(g['finite_impulse'] for g in groups))
    return dict(freeze=FREEZE,evidence=bool(evidence),groups=groups,
                verdicts=dict(zip(VERDICTS,(crossing,balance,impulse,False))),
                reciprocal_exchange_obstruction='No independent gravitating receiver or reciprocal transfer ledger in this test-particle experiment.',
                parent_initial_data_verdicts='8/8 retained; not revised by this evolution',
                unestablished=['global traversability','long-time stability','discrete reciprocal exchange','action quanta','quantum statistics'])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--rescore',type=Path)
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    target=args.output_dir/'evolution_verdict.json'
    target.write_text(json.dumps(dict(verdicts=dict.fromkeys(VERDICTS,False),error='incomplete'))+'\n')
    try:
        data=json.loads(args.rescore.read_text()) if args.rescore else run()
        result=score(data,replay=bool(args.rescore))
        (args.output_dir/'evolution.json').write_text(json.dumps(data,separators=(',',':'),allow_nan=False)+'\n')
        target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(json.dumps(result,indent=2),flush=True)
        return 0 if result['evidence'] and all(g['evolution_valid'] for g in result['groups']) else 1
    except Exception as error:
        target.write_text(json.dumps(dict(verdicts=dict.fromkeys(VERDICTS,False),error=str(error)))+'\n')
        raise


if __name__=='__main__':raise SystemExit(main())
