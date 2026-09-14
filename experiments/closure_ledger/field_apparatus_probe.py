"""Frozen independent controls and scoped audit; no fitted interface model."""
import argparse
import copy
import gzip
import itertools
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves import field_apparatus as f
from . import field_apparatus_audit as audit

GATES = ('source_inventory','quartet_identity','graph_control',
         'kinetic_measure_and_coordinate_control','physical_potential_control',
         'map_and_interface_equations','reciprocal_action_and_full_residuals',
         'preparation_ledger','raw_evidence_and_failure_paths','scope_and_provenance')
VERIFICATION = tuple(k for k in GATES if k not in ('map_and_interface_equations','reciprocal_action_and_full_residuals'))
PREPARATION = dict(fields='FOUR_REAL_CONFORMAL_SCALARS_CHOSEN',
                   geometry='MOUTH_NOT_IDENTIFIED_WITH_BULK',
                   drive='GRAPH_AND_POTENTIAL_PRESCRIBED_CONTROLS',
                   phase='PREPARED_CONTROL_PARAMETER', orientation='PREPARED_CONTROL_PARAMETER',
                   population='NOT_DERIVED', physical_record_law='NOT_DERIVED')
SCOPE = dict(Phi_selection='NOT_DERIVED', Born_rule='NOT_DERIVED',
             operational_causality='NOT_ESTABLISHED',
             noncommuting_classical_maps_are_quantum_commutators=False,
             intrinsic_circle_is_embedded_mouth=False,
             audit_is_universal_no_go=False)


def graph_cases():
    return [dict(m=m,phase=phase) for m,phase in itertools.product((1,2,3),f.PHASES)]


def potential_cases():
    return [dict(v=v,m=m,phase=phase,N=N) for v,m,phase,N in itertools.product((.05,.1),(1,2,3),f.PHASES,f.GRIDS)]


def run_probe():
    # Only this classical graph operator is called from the older quantum probe.
    from .sigma_z_readout_capstone_probe import fiber_H
    report = dict(prereg=f.PREREG,baseline=f.BASELINE,seed=f.SEED,
                  exact=f.certificate(),source_audit=audit.inventory(),
                  preparation=PREPARATION.copy(),scope=SCOPE.copy())
    report['quartet'] = [{**c,**f.quartet(np.array(c['q']),c['A'],c['circle'],c['N'])} for c in f.quartet_cases()]
    report['selected_component'] = [dict(N=N,coefficient=f.encode(np.mean(np.cos(f.grid(N))**2*np.exp(-2j*f.grid(N))))) for N in f.GRIDS]
    graphs = []
    for c in graph_cases():
        levels = []
        for h in (.01,.005,.0025):
            levels.append(dict(h=h,plus=f.graph(h,**c).tolist(),minus=f.graph(-h,**c).tolist(),
                               reference_plus=fiber_H(h,c['phase'],c['m']).real.tolist(),
                               reference_minus=fiber_H(-h,c['phase'],c['m']).real.tolist()))
        graphs.append({**c,'exact_derivative':f.graph(0.,**c,derivative=True).tolist(),'levels':levels})
    report['graph'] = graphs
    report['circle'] = [{**c,**f.circle(**c)} for c in f.circle_cases()]
    report['potential'] = [{**c,'block':f.encode(f.potential(**c))} for c in potential_cases()]
    return report


def finite(value):
    try:
        json.dumps(value,allow_nan=False)
        return True
    except (TypeError,ValueError,OverflowError):
        return False


def same_cases(rows, cases):
    if len(rows)!=len(cases):return False
    # Ordering is frozen, making duplication and omitted phase cases detectable.
    return all(all(row.get(k)==v for k,v in case.items()) for row,case in zip(rows,cases))


def quartet_valid(r):
    cert=f.certificate()
    if any(r['exact'][k]!=cert[k] for k in ('quaternion_gram','intensity')):return False
    if set(cert['quaternion_gram'])!={'0'} or cert['intensity']!='0':return False
    if not same_cases(r['quartet'],list(f.quartet_cases())):return False
    for row in r['quartet']:
        Q=np.dot(row['q'],row['q'])/row['A']**2
        coeff=f.decode(row['fourier'])
        if coeff.shape!=(3,) or max(abs(coeff))/max(1.,Q)>=1e-12:return False
        if abs(row['mean']-Q)/max(1.,Q)>=1e-12:return False
    if [row['N'] for row in r['selected_component']]!=list(f.GRIDS):return False
    return all(abs(f.decode(row['coefficient'])-.25)<1e-12 for row in r['selected_component'])


def graph_errors(row):
    exact=np.asarray(row['exact_derivative'])
    return [f.relative((np.asarray(v['plus'])-v['minus'])/(2*v['h']),exact) for v in row['levels']]


def graph_valid(r):
    if not same_cases(r['graph'],graph_cases()):return False
    for row in r['graph']:
        c={k:row[k] for k in ('m','phase')}
        exact=np.asarray(row['exact_derivative'])
        if exact.shape!=(8,8) or f.relative(exact,f.graph(0.,**c,derivative=True))>=1e-12:return False
        expected=(2-np.sqrt(2))*np.exp(1j*row['phase']) if row['m']==2 else 0j
        if abs(f.graph_block(exact)[0,1]-expected)>=1e-12:return False
        if [v['h'] for v in row['levels']]!=[.01,.005,.0025]:return False
        for v in row['levels']:
            for sign,s in (('plus',1),('minus',-1)):
                actual=np.asarray(v[sign])
                if actual.shape!=(8,8) or f.relative(actual,np.asarray(v['reference_'+sign]))>=1e-12:return False
                if f.relative(actual,f.graph(s*v['h'],**c))>=1e-12:return False
        err=graph_errors(row)
        if err[-1]>=1e-4:return False
        if any(min(x,y)>1e-10 and not 3.5<=x/y<=4.5 for x,y in zip(err,err[1:])):return False
    return True


def circle_valid(r):
    cert=f.certificate()
    for key in ('circle_blocks','mapped_equation'):
        if r['exact'][key]!=cert[key]:return False
    if set(cert['mapped_equation'])!={'0'}:return False
    if any(set(b['generalized'])!={'0'} for b in cert['circle_blocks']):return False
    if not same_cases(r['circle'],list(f.circle_cases())):return False
    for row in r['circle']:
        c={k:row[k] for k in ('R0','epsilon','m','phase','N')}
        lhs,rhs=f.decode(row['lhs']),f.decode(row['rhs'])
        if lhs.shape!=(4,row['N']) or rhs.shape!=lhs.shape:return False
        # Per-mode residual and the maximum pointwise scaled residual both pass.
        if any(f.relative(x,y)>=1e-10 for x,y in zip(lhs,rhs)):return False
        if np.max(abs(lhs-rhs)/np.maximum(1.,np.maximum(abs(lhs),abs(rhs))))>=1e-10:return False
        predicted=f.circle(**c)
        for key in ('lhs','rhs','K1','W1','kinetic_gram'):
            if f.relative(f.decode(row[key]),f.decode(predicted[key]))>=1e-12:return False
        K,W=f.decode(row['K1']),f.decode(row['W1'])
        if f.relative(K,W/row['R0']**2)>=1e-12:return False
        if abs(row['circumference']-2*np.pi*row['R0'])>=1e-10:return False
        if row['N']==128 and f.relative(f.decode(row['kinetic_gram']),np.eye(4))>=1e-10:return False
    return True


def potential_valid(r):
    if r['exact']['potential_selection']!=f.certificate()['potential_selection']:return False
    if not same_cases(r['potential'],potential_cases()):return False
    for row in r['potential']:
        off=row['v']*np.exp(1j*row['phase'])/2 if row['m']==2 else 0j
        expected=np.array([[0,off],[off.conjugate(),0]])
        if f.relative(f.decode(row['block']),expected)>=1e-12:return False
    return True


def provenance_valid(r):
    return (r['prereg']==f.PREREG and r['baseline']==f.BASELINE and r['seed']==f.SEED
            and r['scope']==SCOPE)


def _base_gates(r):
    functions=dict(source_inventory=lambda:audit.valid(r['source_audit']),
                   quartet_identity=lambda:quartet_valid(r), graph_control=lambda:graph_valid(r),
                   kinetic_measure_and_coordinate_control=lambda:circle_valid(r),
                   physical_potential_control=lambda:potential_valid(r),
                   # This implementation has no physical map or response certificate.
                   # No user-supplied boolean or graph success can manufacture one.
                   map_and_interface_equations=lambda:False,
                   reciprocal_action_and_full_residuals=lambda:False,
                   preparation_ledger=lambda:r['preparation']==PREPARATION,
                   scope_and_provenance=lambda:provenance_valid(r))
    sections=dict(source_inventory=['source_audit'],quartet_identity=['quartet','selected_component'],
                  graph_control=['graph'],kinetic_measure_and_coordinate_control=['circle'],
                  physical_potential_control=['potential'],preparation_ledger=['preparation'],
                  scope_and_provenance=['scope'])
    out={}
    for key,fn in functions.items():
        try:
            out[key]=bool(all(finite(r[s]) for s in sections.get(key,[])) and fn())
        except (KeyError,ValueError,TypeError,IndexError,OverflowError,np.linalg.LinAlgError):
            out[key]=False
    # Provisional value internal to the base evaluator. Both public entry
    # points replace it with exercised failure-control results before returning
    # verdicts; tamper trials use the same _verdict selector without recursion.
    out['raw_evidence_and_failure_paths']=True
    return out


def _verdict(checks, actual):
    ok=lambda *names:all(checks.get(k) is True and actual.get(k) is True for k in names)
    global_ok=ok('scope_and_provenance','raw_evidence_and_failure_paths') and not set(checks)-set(GATES)
    result=dict(bulk_mouth_map='UNRESOLVED',operator_status=dict(reference_lattice='UNRESOLVED',intrinsic_circle='UNRESOLVED',physical_mouth='UNRESOLVED'),
                quartet_intensity='UNRESOLVED',field_response='UNRESOLVED',
                preparation_status={k:'UNRESOLVED' for k in PREPARATION},physical_mouth_matrix_element=None,
                field_generated_amplitude=None,**SCOPE)
    if not global_ok:return result
    if ok('source_inventory'):
        result['bulk_mouth_map']='BULK_MOUTH_MAP_UNSPECIFIED'
        result['field_response']='BLOCKED_BY_UNSPECIFIED_MAP'
    if ok('graph_control'):result['operator_status']['reference_lattice']='REFERENCE_LATTICE_REPRODUCED'
    if ok('kinetic_measure_and_coordinate_control'):result['operator_status']['intrinsic_circle']='INTRINSIC_CIRCLE_MIXING_CANCELS'
    if ok('quartet_identity'):result['quartet_intensity']='NO_ANGULAR_SOURCE_IN_EXACT_QUARTET'
    if ok('source_inventory','preparation_ledger'):result['preparation_status']=PREPARATION.copy()
    return result


def failure_controls(r,actual):
    rows=[]
    # Mutate the real source structures but only copy the branch being changed.
    for kind,gate in [('missing_W','kinetic_measure_and_coordinate_control'),
                      ('matrix','graph_control'),('nonfinite','physical_potential_control'),
                      ('symbolic','quartet_identity'),('inventory','source_inventory')]:
        damaged=r.copy()
        section=dict(missing_W='circle',matrix='graph',nonfinite='potential',
                     symbolic='exact',inventory='source_audit')[kind]
        try:
            if kind=='missing_W':
                damaged['circle']=[dict(v) for v in r['circle']];damaged['circle'][0].pop('W1',None)
            elif kind=='matrix':
                damaged['graph']=copy.deepcopy(r['graph']);damaged['graph'][0]['levels'][0]['plus'][0][0]+=.1
            elif kind=='nonfinite':
                damaged['potential']=copy.deepcopy(r['potential']);damaged['potential'][0]['block'][0][0][0]=float('nan')
            elif kind=='symbolic':
                damaged['exact']=copy.deepcopy(r['exact']);damaged['exact']['intensity']='1'
            else:
                damaged['source_audit']=copy.deepcopy(r['source_audit']);damaged['source_audit']['rows'].pop()
        except (KeyError,ValueError,TypeError,IndexError):
            # Already missing/malformed input must fail its own gate, without
            # making the validation machinery erase unrelated evidence.
            damaged[section]=None
        changed=_base_gates(damaged)
        out=_verdict(actual,changed)
        passed=not changed[gate] and all(changed[k]==actual[k] for k in GATES if k!=gate)
        # The result itself, not just a gate flag, must withdraw its target.
        if kind=='missing_W':passed &= out['operator_status']['intrinsic_circle']=='UNRESOLVED'
        elif kind=='matrix':passed &= out['operator_status']['reference_lattice']=='UNRESOLVED'
        elif kind=='symbolic':passed &= out['quartet_intensity']=='UNRESOLVED'
        elif kind=='inventory':passed &= out['bulk_mouth_map']=='UNRESOLVED' and out['field_response']=='UNRESOLVED'
        rows.append(dict(control=kind,gate=gate,changed_gates=[k for k in GATES if changed[k]!=actual[k]],verdict=out,passed=bool(passed)))
    forced={k:True for k in GATES}
    out=_verdict(forced,actual)
    rows.append(dict(control='prescribed_graph_cannot_supply_map',passed=out['field_response'] in ('BLOCKED_BY_UNSPECIFIED_MAP','UNRESOLVED'),verdict=out))
    return rows


def evidence_gates(r):
    checks=_base_gates(r)
    try:
        checks['raw_evidence_and_failure_paths']=all(v['passed'] for v in failure_controls(r,checks))
    except (KeyError,ValueError,TypeError,IndexError,OverflowError):
        checks['raw_evidence_and_failure_paths']=False
    return checks


def verdict(checks,r):
    return _verdict(checks,evidence_gates(r))


def finalize(r):
    checks=_base_gates(r)
    controls=failure_controls(r,checks)
    checks['raw_evidence_and_failure_paths']=all(v['passed'] for v in controls)
    r.update(checks=checks,failure_controls=controls,verdict=_verdict(checks,checks),
             verification_passed=all(checks[k] for k in VERIFICATION))
    return r


def serialize(report):
    """Retain independently earned verdicts when one evidence section is NaN."""
    paths=[]
    def clean(x,path):
        if isinstance(x,float) and not np.isfinite(x):
            paths.append(path);return None
        if isinstance(x,dict):return {k:clean(v,path+'/'+k) for k,v in x.items()}
        if isinstance(x,list):return [clean(v,path+'/'+str(i)) for i,v in enumerate(x)]
        return x
    data=clean(report,'')
    if paths:data['nonfinite_evidence_paths']=paths
    return json.dumps(data,sort_keys=True,separators=(',',':'),allow_nan=False).encode()+b'\n'


def main(argv=None):
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        report=finalize(run_probe())
        payload=serialize(report)
    except Exception as exc:
        report=dict(verification_passed=False,error=type(exc).__name__+': '+str(exc),verdict=verdict({},{}))
        payload=json.dumps(report,sort_keys=True).encode()+b'\n'
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/'probe.json.gz').write_bytes(gzip.compress(payload,mtime=0))
    (args.output_dir/'probe.md').write_text('# Field-to-apparatus controls\n\nVerification passed: '+str(report['verification_passed'])+'\n\n```json\n'+json.dumps(report['verdict'],indent=2)+'\n```\n')
    print(json.dumps(dict(verification_passed=report['verification_passed'],checks=report.get('checks'),error=report.get('error'),verdict=report['verdict']),indent=2))
    return 0 if report['verification_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
