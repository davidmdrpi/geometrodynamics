"""Prospective smaller-amplitude extension; the original failed gates remain intact."""
import argparse
import json
import gzip
from pathlib import Path
import numpy as np
from . import nonlinear_supported_tt_probe as p
from geometrodynamics.waves import nonlinear_supported_tt as n

PREREG='75c7ab46939c2f883edbe9d126800272157180d5'


def run_probe(original):
    pairs=p.pairs();rows=[]
    for v in original['variations']:
        ip=v['pair'];phase=v['phase'];U,V=pairs[ip];pred=n.second_variations(U,V,p.TIMES,phase=phase)
        zero=next(r['states'] for r in original['trajectories'] if r['pair']==ip and r['phase']==phase and r['epsilon']==0)
        base=np.array([n.observables(np.array(y)) for y in zero]);levels=[]
        for eps in (.01,.005,.0025):
            data=[]
            for sign in (-1,1):
                matches=[r['states'] for r in original['trajectories'] if r['pair']==ip and r['phase']==phase and r['epsilon']==sign*eps]
                states=np.array(matches[0]) if matches else n.evolve(n.initial_data(U,V,sign*eps,phase=phase),p.TIMES)
                data.append(states)
            minus,plus=[np.array([n.observables(y) for y in states]) for states in data]
            first=(plus-minus)/(2*eps);second=(plus+minus-2*base)/(2*eps*eps)
            levels.append(dict(epsilon=eps,minus=data[0].tolist(),plus=data[1].tolist(),
                first=first.tolist(),second=second.tolist(),
                first_errors=[p.relative(z,r['first']) for z,r in zip(first,pred)],
                second_errors=[p.relative(z,r['second']) for z,r in zip(second,pred)]))
        maxima={key:[max(row[key][1:5]) for row in levels] for key in ('first_errors','second_errors')}
        ratios={key:[values[i]/values[i+1] if min(values[i:i+2])>1e-8 else None for i in range(2)] for key,values in maxima.items()}
        rows.append(dict(pair=ip,phase=phase,levels=levels,maxima=maxima,ratios=ratios))
    return dict(prereg=PREREG,original_prereg=n.PUBLIC_PREREG,original_checks=original['checks'],original_verdict=original['verdict'],rows=rows)


def valid_rows(r,original):
    try:
        if r['prereg']!=PREREG or len(r['rows'])!=29:return False
        json.dumps(r,allow_nan=False)
        expected={(i,phase) for i in range(15) for phase in ((0.,np.pi/4,np.pi/2) if i<7 else (np.pi/4,))}
        if {(v['pair'],v['phase']) for v in r['rows']}!=expected:return False
        predictions={(v['pair'],v['phase']):v for v in original['variations']}
        zeros={(v['pair'],v['phase']):np.array([n.observables(np.asarray(y)) for y in v['states']]) for v in original['trajectories'] if v['epsilon']==0}
        for v in r['rows']:
            if [row['epsilon'] for row in v['levels']]!=[.01,.005,.0025]:return False
            pred=predictions[v['pair'],v['phase']];zero=zeros[v['pair'],v['phase']]
            errors=dict(first=[],second=[])
            for row in v['levels']:
                minus=np.asarray(row['minus']);plus=np.asarray(row['plus']);eps=row['epsilon']
                if minus.shape!=(7,29) or plus.shape!=(7,29):return False
                x=np.array([n.observables(y) for y in plus]);y=np.array([n.observables(y) for y in minus])
                values=dict(first=(x-y)/(2*eps),second=(x+y-2*zero)/(2*eps*eps))
                for key,val in values.items():
                    if p.relative(val,np.asarray(row[key]))>=1e-10:return False
                    err=[p.relative(z,b) for z,b in zip(val,pred[key+'_prediction'])]
                    if max(abs(np.asarray(err)-row[key+'_errors']))>=1e-10:return False
                    errors[key].append(max(err[1:5]))
                    if eps==.0025 and max(err[1:5])>=.001:return False
            for vals in errors.values():
                for i in range(2):
                    if min(vals[i:i+2])>1e-8 and not 3.5<=vals[i]/vals[i+1]<=4.5:return False
        return True
    except (KeyError,TypeError,ValueError,OverflowError,ZeroDivisionError,IndexError,np.linalg.LinAlgError):return False


def constraints_valid(r):
    """Check every refinement state, including times outside the accuracy window."""
    try:
        if len(r['rows'])!=29:return False
        for v in r['rows']:
            if [row['epsilon'] for row in v['levels']]!=[.01,.005,.0025]:return False
            for row in v['levels']:
                for sign in ('minus','plus'):
                    states=np.asarray(row[sign],dtype=float)
                    if states.shape!=(7,29) or not np.isfinite(states).all():return False
                    for y in states:
                        residual=n.constraints(y)['normalized']
                        if not np.isfinite(residual).all() or max(residual)>=1e-8:return False
        return True
    except (KeyError,TypeError,ValueError,OverflowError,IndexError,np.linalg.LinAlgError):return False


def evidence_checks(r,original):
    valid=valid_rows(r,original)
    checks=p.evidence_gates(original)
    checks['constraint_propagation']=bool(checks['constraint_propagation'] and constraints_valid(r))
    try:operator=all(max(v['linear_operator_errors']+v['zero_solution_errors'])<1e-8 for v in original['variations'])
    except (KeyError,TypeError,ValueError):operator=False
    checks['linear_recovery']=bool(valid and operator)
    checks['quadratic_response']=bool(valid)
    return checks


def verdict(checks,r,original):
    actual=evidence_checks(r,original)
    failures={t:[k for k,targets in p.DEPENDENCIES.items() if t in targets and (checks.get(k) is not True or actual[k] is not True)] for t in p.TARGETS}
    if set(checks)-set(p.DEPENDENCIES):
        for t in failures:failures[t].append('unknown_gate')
    labels=dict(D='EXACT_HOMOGENEOUS_REDUCTION_VERIFIED',C='LOCAL_CONSTRAINT_COMPLETED_FAMILIES_VERIFIED',
                N='FINITE_AMPLITUDE_RESPONSE_VERIFIED_AFTER_PROSPECTIVE_REFINEMENT',F='FUTURE_PERSISTENCE_PROVED_WITH_REFINED_TANGENT_CHECK')
    return {**{name:labels[t] if not failures[t] else 'UNRESOLVED' for t,name in p.TARGETS.items()},**p.SCOPE,'failed_checks':failures}


def finalize(r,original):
    # The original report is neither modified nor relabeled by this extension.
    checks=evidence_checks(r,original)
    r.update(checks=checks,checks_passed=all(checks.values()),verdict=verdict(checks,r,original))
    return r


def main(argv=None):
    parser=argparse.ArgumentParser();parser.add_argument('--original',type=Path,required=True)
    parser.add_argument('--output-dir','--output',dest='output',type=Path,required=True);args=parser.parse_args(argv)
    try:
        raw=gzip.decompress(args.original.read_bytes()).decode() if args.original.suffix=='.gz' else args.original.read_text()
        original=json.loads(raw);report=finalize(run_probe(original),original)
        payload=p.serialize_report(report)
    except Exception as exc:
        report=dict(checks_passed=False,error=type(exc).__name__+': '+str(exc),verdict=p.verdict({}))
        payload=p.serialize_report(report)
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'refinement.json').write_text(payload)
    (args.output/'refinement.md').write_text('# Prospective amplitude refinement\n\nOriginal frozen gates remain 11/13.\n\n```json\n'+json.dumps(report.get('verdict'),indent=2)+'\n```\n')
    return 0 if report['checks_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
