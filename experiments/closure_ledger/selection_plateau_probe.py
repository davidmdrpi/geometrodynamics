"""Run or replay the frozen finite-time nonlinear locking test."""
import argparse
import base64
import gzip
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from geometrodynamics.waves import selection_plateau as p

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'experiments/closure_ledger/runs/20260928_selection_plateau'


def schedule():
    return [(float(e),float(t),float(f)) for e in p.AMPLITUDES
            for t in p.SHAPE_PHASES for f in p.FIELD_PHASES]


def control_schedule():
    return [(float(p.AMPLITUDES[-1]),float(t),float(f)) for t in p.SHAPE_PHASES for f in p.FIELD_PHASES]+[(0.,0.,float(f)) for f in p.FIELD_PHASES]


def integrate_case(parameters, method):
    try:
        return dict(parameters=list(parameters),states=p.evolve(*parameters,method=method).tolist())
    except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
        return dict(parameters=list(parameters),failure=str(error))


def validate_raw(raw):
    if not np.array_equal(raw.get('times'),p.TIMES):
        raise ValueError('wrong time schedule')
    for group,expected in [('cases',schedule()),('controls',control_schedule())]:
        if [tuple(row['parameters']) for row in raw.get(group,[])]!=expected:
            raise ValueError('incomplete or relabelled preparation schedule')
        for row in raw[group]:
            if 'failure' in row:
                continue
            states=np.asarray(row['states'])
            if states.shape!=(len(p.TIMES),31) or not np.isfinite(states).all():
                raise ValueError('invalid states')
            if not np.allclose(states[0,:29],p.preparation(*row['parameters']),rtol=1e-12,atol=1e-14) or np.any(states[0,29:]!=0):
                raise ValueError('changed preparation or work origin')


def summarize(raw):
    validate_raw(raw)
    result=dict(freeze=p.FREEZE,readouts=['J_ref','J_inst'],durations=p.DURATIONS.tolist(),J_bg=p.J_BG,
                scope='Homogeneous finite-time parametric locking; no localized absorbed action',cases=[],controls=[])
    gates=dict(complete=True,constraints=True,chart=True,work_balance=True,background=True,
               window_quadrature=True,independent_integrator=True)
    for group in ('cases','controls'):
        for row in raw[group]:
            if 'failure' in row:
                gates['complete']=False
                result[group].append(row.copy())
                continue
            d=p.summarize_history(row['states'])
            result[group].append(dict(parameters=row['parameters'],diagnostics=d))
            gates['constraints'] &= d['constraint_max']<1e-8
            gates['chart'] &= d['det_error']<1e-7 and d['symmetry_error']<1e-7 and min(d['sampled_min_A_H_M_omega2'])>0
            gates['work_balance'] &= d['work_error_over_Jbg']<1e-8
            gates['window_quadrature'] &= d['quadrature_error_over_Jbg']<1e-6
            if row['parameters'][0]==0:
                gates['background'] &= np.max(abs(np.array(d['changes'])))/p.J_BG<1e-9
    diffs=[]
    if gates['complete']:
        primary={tuple(row['parameters']):row['diagnostics'] for row in result['cases']}
        for row in result['controls']:
            d=primary[tuple(row['parameters'])]
            error=float(np.max(abs(np.array(d['window_means'])-row['diagnostics']['window_means']))/p.J_BG)
            diffs.append(dict(parameters=row['parameters'],difference_over_Jbg=error))
        gates['independent_integrator']=max(r['difference_over_Jbg'] for r in diffs)<1e-6
        means=np.array([row['diagnostics']['window_means'][1:] for row in result['cases'][12:]]).reshape(8,12,3,2)
        windows=p.classify(means)
        # Full signed readouts are above. Resolved slopes are secondary diagnostics.
        slopes=[]
        for i in range(7):
            low,high=means[i],means[i+1]
            slopes.append([[[float(np.log(abs(h/l))/np.log(1.5)) if abs(l)>p.FLOOR and abs(h)>p.FLOOR and l*h>0 else None
                             for l,h in zip(a,b)] for a,b in zip(c,d)] for c,d in zip(low,high)])
        result['signed_response_log_slopes']=slopes
    else:
        windows=[]
        gates['independent_integrator']=False
    result['integrator_differences']=diffs
    result['candidate_windows']=windows
    result['numerical_gates']={k:bool(v) for k,v in gates.items()}
    result['selection_verdict']=p.selection_verdict(gates,windows)
    result['localized_receiver_verdict']='NOT_TESTED'
    return result


def load_raw(path):
    return json.loads(gzip.decompress(base64.b64decode(Path(path).read_text())))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=RUN)
    parser.add_argument('--replay',type=Path)
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    if args.replay:
        raw=load_raw(args.replay)
    else:
        raw=dict(times=p.TIMES.tolist(),cases=[],controls=[])
        for group,params,method in [('cases',schedule(),'DOP853'),('controls',control_schedule(),'RK45')]:
            for i,parameters in enumerate(params):
                if i%12==0:print(group,i,'of',len(params),flush=True)
                raw[group].append(integrate_case(parameters,method))
    archive=args.output_dir/'states.json.gz.b64'
    archive.write_text(base64.b64encode(gzip.compress(json.dumps(raw,separators=(',',':'),allow_nan=False).encode(),mtime=0)).decode()+'\n')
    result=summarize(raw)
    result['environment']=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    result['source_sha256']={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in [
        'geometrodynamics/waves/selection_plateau.py','geometrodynamics/waves/nonlinear_supported_tt.py',
        'geometrodynamics/waves/action_selection.py','experiments/closure_ledger/selection_plateau_probe.py',
        'docs/nonlinear_selection_plateau_prereg.md']}
    result['archive_sha256']=hashlib.sha256(archive.read_bytes()).hexdigest()
    (args.output_dir/'plateau.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('numerical_gates','selection_verdict')},indent=2))
    return 1 if result['selection_verdict']=='INCONCLUSIVE_NUMERICAL_FAILURE' else 0


if __name__=='__main__':
    raise SystemExit(main())
