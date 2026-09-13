"""Frozen nonlinear completion experiment. Raw failed cases are retained."""
import argparse
import copy
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from geometrodynamics.waves import nonlinear_supported_tt as n
from geometrodynamics.waves import frw_supported_tt as f

TIMES=[0.,.25,.5,1.,2.,4.,8.]
AMPLITUDES=[0.,-.01,.01,-.02,.02,-.04,.04,-.08,.08]
DEPENDENCIES=dict(conventions_and_units='DCNF',exact_field_closure='DCNF',action_full_field_agreement='DCNF',
 constraint_derivation='DCNF',constraint_completion='CNF',constraint_propagation='NF',linear_recovery='NF',
 quadratic_response='N',finite_time_evolution='NF',future_continuation_bound='F',negative_controls='DCNF',scope='DCNF',failure_paths='DCNF')
TARGETS=dict(D='exact_reduction',C='constraint_completed_families',N='finite_amplitude_response',F='future_persistence')
SCOPE=dict(Phi_selection='NOT_DERIVED',quantization='NOT_DERIVED',causality_gate='OPEN',
          nonlinear_rotor='NOT_DERIVED',inhomogeneous_stability='NOT_ESTABLISHED',preparation_selection='NOT_DERIVED',
          tensor_only_symplectic_map='NOT_ASSERTED',explicit_initial_epsilon_radius='NOT_COMPUTED')


def relative(a,b):return float(np.linalg.norm(np.asarray(a)-b)/max(1.,np.linalg.norm(b)))


def pairs():
    E=n.rt.STF_BASIS;z=np.zeros((3,3))
    out=[(E[0],E[0]),(E[0],E[1]),(E[0],E[2]),(E[2],E[3]),((E[0]+E[3])/math.sqrt(2),(E[1]+E[4])/math.sqrt(2)),(z,E[0]),(E[0],z)]
    rng=np.random.default_rng(n.SEED)
    for i in range(8):
        coeff=rng.normal(size=(2,5));coeff/=np.linalg.norm(coeff,axis=1)[:,None]
        out.append(tuple(n.rt.tensor(coeff)))
    return out


def geometry_record(y,dy=None,a=1.,kappa=1.,angles=(.83,1.07,.61)):
    dy=n.conformal_rhs(y,a,kappa) if dy is None else dy
    g=n.full_geometry(y,dy,a,kappa,angles);p=n.expected_residual(y,dy,a,kappa,angles)
    proper=n.full_geometry(y,dy,a,kappa,angles,clock='proper')
    return dict(einstein_raw=float(np.linalg.norm(g['residual']-p['residual'])),
        einstein_normalized=float(a*a*np.linalg.norm(g['residual']-p['residual'])/g['scale']),
        KG_error=float(a*a*math.sqrt(kappa)*np.linalg.norm(g['KG']-p['KG'])),
        clock_error=relative(g['residual'],proper['residual']),
        full_residual=g['residual'].tolist(),expected_residual=p['residual'].tolist(),
        KG=g['KG'].tolist(),expected_KG=p['KG'].tolist(),scale=g['scale'])


def run_probe():
    ps=pairs();rng=np.random.default_rng(n.SEED+1)
    report=dict(prereg=n.PUBLIC_PREREG,provenance=n.PROVENANCE,baseline=n.BASELINE,seed=n.SEED,
                exact=n.exact_certificate(),future=n.future_certificate(),pairs=[dict(U=u.tolist(),V=v.tolist(),commutator=float(np.linalg.norm(u@v-v@u))) for u,v in ps])
    setup=[];jac=[];off=[];units=[];chirality=[]
    # All 3240 frozen completions, including both field-zero types.
    for ip,(U,V) in enumerate(ps):
        for dep,phase,eps in itertools.product((.05,.15,.3),np.arange(8)*math.pi/8,AMPLITUDES):
            y=n.initial_data(U,V,eps,dep,phase);c=n.constraints(y)
            setup.append(dict(pair=ip,departure=dep,phase=float(phase),epsilon=eps,state=y.tolist(),
                residual=c['residual'].tolist(),normalized=c['normalized'].tolist(),scale=c['scale'].tolist()))
    # Independently differentiate constraints in the four completion variables at epsilon=0.
    for dep,phase in itertools.product((.05,.15,.3),np.arange(8)*math.pi/8):
        y=n.initial_data(*ps[0],0,dep,phase);A,Ap,q,qp,M,L=n.unpack(y);qb=q[0];vb=qp[0];D=qb*qb+vb*vb
        def fun(z):return n.constraints(n.pack(A,z[0],np.r_[qb,-vb*z[1:]/D],np.r_[vb,qb*z[1:]/D],M,L))['residual']
        z=np.r_[Ap,0.,0.,0.];step=1e-5
        J=np.column_stack([(fun(z+step*e)-fun(z-step*e))/(2*step) for e in np.eye(4)])
        expected=np.diag([-6*Ap,-1,-1,-1]);jac.append(dict(departure=dep,phase=float(phase),matrix=J.tolist(),expected=expected.tolist(),error=relative(J,expected)))
    # Independent off-shell jets, including non-diagonal matrices and general accelerations.
    for i in range(20):
        U,V=ps[i%15];eps=float(rng.uniform(-.3,.3));phase=float(rng.uniform(0,math.pi));y=n.initial_data(U,V,eps,phase=phase)
        y[1]*=1.2;y[3:11]+=rng.normal(scale=.07,size=8)
        dy=n.conformal_rhs(y);dy[1]+=.3*rng.normal();dy[7:11]+=rng.normal(scale=.2,size=4)
        A,Ap,q,qp,M,L=n.unpack(y);root=expm(eps*U);W=n.rt.tensor(rng.normal(size=5))
        dy[20:29]+=(np.linalg.solve(root,W@root)).ravel()
        off.append(dict(case=i,state=y.tolist(),jet=dy.tolist(),**geometry_record(y,dy)))
    # Same generators give constant currents; opposite chirality is spatially varying.
    points=[(.5+.07*i,.6+.08*i,.4+.11*i) for i in range(8)]
    other=n.S.copy();other[:,1:,1:]*=-1
    for i in range(3):
        correct=[];wrong=[]
        for angles in points:
            x,dx,_,E,_,_=n.cm.coframe_jets(angles)
            correct.append(np.linalg.solve(E.T,x@n.S[i]@dx))
            wrong.append(np.linalg.solve(E.T,x@other[i]@dx))
        chirality.append(dict(generator=i,correct=np.array(correct).tolist(),wrong=np.array(wrong).tolist()))
    trajectories=[];comparison=[];geometry=[];variations=[];tails=[];controls=[]
    for ip,(U,V) in enumerate(ps):
        phases=(0.,math.pi/4,math.pi/2) if ip<7 else (math.pi/4,)
        epses=AMPLITUDES if ip<7 else (0.,-.02,.02)
        for phase in phases:
            data={}
            for eps in epses:
                initial=n.initial_data(U,V,eps,phase=phase)
                fine=n.evolve(initial,TIMES);coarse=n.evolve(initial,TIMES,rtol=1e-10,atol=1e-12)
                data[eps]=fine
                cons=[n.constraints(y) for y in fine]
                trajectories.append(dict(pair=ip,phase=phase,epsilon=eps,times=TIMES,states=fine.tolist(),
                    constraints=[c['residual'].tolist() for c in cons],constraint_scales=[c['scale'].tolist() for c in cons],
                    normalized_constraints=[c['normalized'].tolist() for c in cons],
                    chart=[dict(det=float(np.linalg.det(n.unpack(y)[4])),symmetry=float(np.linalg.norm(n.unpack(y)[4]-n.unpack(y)[4].T)),
                        trace_L=float(np.trace(n.unpack(y)[5])),F=float(n.ingredients(y)[0]/y[0]**2),
                        self_adjoint=float(np.linalg.norm(n.unpack(y)[4]@n.unpack(y)[5]-(n.unpack(y)[4]@n.unpack(y)[5]).T))) for y in fine]))
                comparison.append(dict(pair=ip,phase=phase,epsilon=eps,coarse_states=coarse.tolist(),errors=[relative(c,f) for c,f in zip(coarse,fine)]))
                if ip<7 and abs(eps)==.08:
                    for ti in (0,2,4):
                        for angles in ((.83,1.07,.61),(1.11,.72,1.27),(.58,1.32,.94)):
                            geometry.append(dict(pair=ip,phase=phase,epsilon=eps,time=TIMES[ti],angles=angles,**geometry_record(fine[ti])))
            # Variations evaluated independently, not differentiated from the nonlinear RHS.
            pred=n.second_variations(U,V,TIMES,phase=phase)
            values={eps:np.array([n.observables(y) for y in states]) for eps,states in data.items()}
            reference_model=f.FRWSupport(departure=.15,phase=phase)
            linear_checks=[]
            for i,p in enumerate(pred):
                sm=f.evolve(reference_model,0,p['eta']);m0=reference_model.coefficients(0)[0]
                b=sm[0,0]*U+sm[0,1]*m0*V
                linear_checks.append(relative(p['first'][5:],n.rt.components(b)))
            rows=[]
            for eps in (.08,.04,.02,.01) if ip<7 else (.02,):
                first=(values[eps]-values[-eps])/(2*eps);second=(values[eps]+values[-eps]-2*values[0])/(2*eps*eps)
                rows.append(dict(epsilon=eps,first=first.tolist(),second=second.tolist(),
                    first_errors=[relative(v,p['first']) for v,p in zip(first,pred)],
                    second_errors=[relative(v,p['second']) for v,p in zip(second,pred)]))
            variations.append(dict(pair=ip,phase=phase,rows=rows,first_prediction=[p['first'].tolist() for p in pred],
                second_prediction=[p['second'].tolist() for p in pred],linear_operator_errors=linear_checks,
                zero_solution_errors=[relative(v,p['baseline']) for v,p in zip(values[0],pred)]))
    # Unit and proper-clock controls, including all 108 prescribed initial states.
    for ip,phase,a,kap,eps in itertools.product((1,2),(0.,math.pi/4,math.pi/2),(.7,1.,2.),(.4,1.),(0.,-.02,.02)):
        U,V=ps[ip];y=n.initial_data(U,V,eps,phase=phase,a=a,kappa=kap);states=n.evolve(y,[0.,.5,2.],a,kap)
        normalized=states.copy();normalized[:,0:2]/=a;normalized[:,3:11]*=math.sqrt(kap)/a
        reference=n.evolve(n.initial_data(U,V,eps,phase=phase),[0.,.5,2.])
        units.append(dict(pair=ip,phase=phase,a=a,kappa=kap,epsilon=eps,state=y.tolist(),
            constraint=float(max(n.constraints(y,a,kap)['normalized'])),map_error=relative(normalized,reference),
            geometry=geometry_record(states[-1],a=a,kappa=kap)))
    for ip,phase,eps in itertools.product((1,2),(0.,math.pi/4,math.pi/2),(0.,-.02,.02,-.08,.08)):
        U,V=ps[ip];samples=[0.,2.,4.,8.,10.,12.,16.];states=n.evolve(n.initial_data(U,V,eps,phase=phase),samples)
        rows=[]
        for y in states:
            A,Ap,q,qp,M,L=n.unpack(y);rows.append(dict(A=float(A),shape=float(np.linalg.norm(M-n.I)),velocity=float(np.linalg.norm(L)),
                matter=float(np.linalg.norm(q)+np.linalg.norm(qp)),F=float(n.ingredients(y)[0]/A**2),
                proper_Hubble=float(Ap/A**2),shear_A2=float(A*np.linalg.norm(L)),field_A=float(np.linalg.norm(q))))
        tails.append(dict(pair=ip,phase=phase,epsilon=eps,times=samples,states=states.tolist(),diagnostics=rows))
    for ip in (1,2):
        U,V=ps[ip];y=n.initial_data(U,V,.08,phase=.4);rigid=n.initial_data(U,V,.08,phase=.4,rigid=True)
        dy=n.conformal_rhs(y);frozen=dy.copy();frozen[1]=-y[0]+y[0]**3;frozen[7:11]=-4*y[3:7]
        good=n.full_geometry(y)
        controls.append(dict(pair=ip,commutator=float(np.linalg.norm(U@V-V@U)),
            rigid_momentum=float(np.linalg.norm(n.constraints(rigid)['residual'][1:])),
            responsive_momentum=float(np.linalg.norm(n.constraints(y)['residual'][1:])),
            rigid_background_error=float(n.full_geometry(y,frozen)['normalized']),
            minimal_stress_error=float(n.full_geometry(y,minimal=True)['normalized']),
            wrong_clock_error=float(n.full_geometry(y,clock='wrong_proper')['normalized']),good_error=float(good['normalized'])))
    report.update(initial=setup,jacobians=jac,off_shell=off,chirality=chirality,trajectories=trajectories,
        comparisons=comparison,geometry=geometry,variations=variations,units=units,tails=tails,controls=controls)
    return report


def _variation_payload(r,kind):
    return [dict(pair=v['pair'],phase=v['phase'],prediction=v[kind+'_prediction'],
        rows=[dict(epsilon=row['epsilon'],values=row[kind],errors=row[kind+'_errors']) for row in v['rows']],
        operator=v['linear_operator_errors'] if kind=='first' else [],
        zero=v['zero_solution_errors'] if kind=='first' else []) for v in r['variations']]


def _initial_valid(r):
    rows=r['initial']
    expected={(ip,d,float(phase),e) for ip in range(15) for d in (.05,.15,.3) for phase in np.arange(8)*math.pi/8 for e in AMPLITUDES}
    if len(rows)!=3240 or {(v['pair'],v['departure'],v['phase'],v['epsilon']) for v in rows}!=expected:return False
    return all(np.asarray(v['state']).shape==(29,) and max(n.constraints(np.asarray(v['state']))['normalized'])<1e-11 for v in rows)


def _propagation_valid(r):
    expected={(ip,float(phase),e) for ip in range(15) for phase in ((0.,math.pi/4,math.pi/2) if ip<7 else (math.pi/4,)) for e in (AMPLITUDES if ip<7 else (0.,-.02,.02))}
    rows=r['trajectories']
    if len(rows)!=213 or {(v['pair'],v['phase'],v['epsilon']) for v in rows}!=expected:return False
    return all(np.asarray(v['states']).shape==(7,29) and max(max(n.constraints(np.asarray(y))['normalized']) for y in v['states'])<1e-8 for v in rows)


def _evolution_valid(r):
    if len(r['comparisons'])!=213:return False
    fine={(v['pair'],v['phase'],v['epsilon']):np.asarray(v['states']) for v in r['trajectories']}
    seen=set()
    for v in r['comparisons']:
        key=(v['pair'],v['phase'],v['epsilon']);seen.add(key)
        coarse=np.asarray(v['coarse_states']);reference=fine[key]
        if coarse.shape!=(7,29) or reference.shape!=(7,29):return False
        for left,right in zip(coarse,reference):
            for block in (slice(0,1),slice(1,2),slice(2,3),slice(3,7),slice(7,11),slice(11,20),slice(20,29)):
                if relative(left[block],right[block])>=1e-8:return False
        for y in reference:
            _,_,_,_,M,L=n.unpack(y)
            if abs(np.linalg.det(M)-1)>=1e-8 or abs(np.trace(L))>=1e-8 or np.linalg.norm(M-M.T)>=1e-8 or np.linalg.norm(M@L-(M@L).T)>=1e-8:return False
            n.ingredients(y)
    return seen==set(fine)


def _variation_valid(r,kind):
    rows=[v for v in r['variations'] if v['pair']<7]
    expected={(i,phase) for i in range(7) for phase in (0.,math.pi/4,math.pi/2)}
    if len(rows)!=21 or {(v['pair'],v['phase']) for v in rows}!=expected:return False
    for v in rows:
        levels=[row for row in v['rows'] if row['epsilon']==.01]
        if len(levels)!=1:return False
        values=np.asarray(levels[0][kind]);prediction=np.asarray(v[kind+'_prediction'])
        if values.shape!=(7,10) or prediction.shape!=(7,10):return False
        if max(relative(x,y) for x,y in zip(values[:5],prediction[:5]))>=1e-3:return False
    return True


def evidence_gates(r):
    # Each gate validates its own raw evidence; malformed future evidence must
    # not erase an independently verified local constraint-completion result.
    cert=lambda:r['exact']==n.exact_certificate() and r['exact']['all_zero']
    geom=lambda:r['geometry']+r['off_shell']+[v['geometry'] for v in r['units']]
    tests=dict(
        conventions_and_units=lambda:len(r['units'])==108 and all(max(v['constraint'],v['map_error'])<1e-8 for v in r['units'])
          and len(r['chirality'])==3 and all(np.max(abs(np.array(v['correct'])+np.eye(3)[v['generator']]))<1e-12 for v in r['chirality']),
        exact_field_closure=cert,
        action_full_field_agreement=lambda:len(r['geometry'])==378 and len(r['off_shell'])==20 and all(max(v['einstein_normalized'],v['KG_error'],v['clock_error'])<1e-8 for v in geom()),
        constraint_derivation=lambda:cert() and len(r['jacobians'])==24 and all(v['error']<1e-8 for v in r['jacobians']),
        constraint_completion=lambda:_initial_valid(r),
        constraint_propagation=lambda:_propagation_valid(r),
        linear_recovery=lambda:_variation_valid(r,'first') and all(max(v['linear_operator_errors']+v['zero_solution_errors'])<1e-8 for v in r['variations']),
        quadratic_response=lambda:_variation_valid(r,'second'),
        finite_time_evolution=lambda:_evolution_valid(r),
        future_continuation_bound=lambda:r['future']==n.future_certificate() and bool(r['future']['all_positive']),
        negative_controls=lambda:len(r['controls'])==2 and r['controls'][0]['rigid_momentum']<1e-12 and r['controls'][1]['rigid_momentum']>1e-3
          and all(v['responsive_momentum']<1e-12 and min(v['rigid_background_error'],v['minimal_stress_error'],v['wrong_clock_error'])>1e-6 for v in r['controls'])
          and all(np.max(np.ptp(np.array(v['wrong']),axis=0))>.1 for v in r['chirality']),scope=lambda:True,failure_paths=lambda:True)
    sources=dict(conventions_and_units=['units','chirality'],exact_field_closure=['exact'],
        action_full_field_agreement=['geometry','off_shell','units'],constraint_derivation=['exact','jacobians'],
        constraint_completion=['initial'],constraint_propagation=['trajectories'],finite_time_evolution=['trajectories','comparisons'],
        future_continuation_bound=['future'],negative_controls=['controls','chirality'],scope=[],failure_paths=[])
    output={}
    for gate,test in tests.items():
        try:
            if gate in ('linear_recovery','quadratic_response'):
                payload=_variation_payload(r,'first' if gate=='linear_recovery' else 'second')
            else:payload=[r[key] for key in sources[gate]]
            json.dumps(payload,allow_nan=False)
            output[gate]=bool(test())
        except (KeyError,ValueError,TypeError,OverflowError,IndexError,ZeroDivisionError,np.linalg.LinAlgError):output[gate]=False
    return output


def verdict(checks,evidence=None):
    failed={t:[] for t in TARGETS}
    for gate,targets in DEPENDENCIES.items():
        if checks.get(gate) is not True:
            for target in targets:failed[target].append(gate)
    if set(checks)-set(DEPENDENCIES):
        for target in TARGETS:failed[target].append('unknown_gate')
    if evidence is None:
        for target in TARGETS:failed[target].append('missing_evidence')
    else:
        try:
            actual=evidence_gates(evidence)
            for gate,targets in DEPENDENCIES.items():
                if actual[gate] is not True:
                    for target in targets:failed[target].append('evidence_'+gate)
        except (KeyError,ValueError,TypeError,OverflowError):
            for target in TARGETS:failed[target].append('malformed_evidence')
    labels=dict(D='EXACT_HOMOGENEOUS_REDUCTION_VERIFIED',C='LOCAL_CONSTRAINT_COMPLETED_FAMILIES_VERIFIED',
        N='FINITE_AMPLITUDE_RESPONSE_VERIFIED',F='FUTURE_PERSISTENCE_PROVED_ON_STATED_DOMAIN')
    return {**{name:labels[t] if not failed[t] else 'UNRESOLVED' for t,name in TARGETS.items()},**SCOPE,'failed_checks':failed}


def finalize(r):
    checks={k:bool(v) for k,v in evidence_gates(r).items()}
    good=dict.fromkeys(DEPENDENCIES,True);failures=True
    for gate,targets in DEPENDENCIES.items():
        for missing in (False,True):
            bad=good.copy()
            if missing:bad.pop(gate)
            else:bad[gate]=False
            out=verdict(bad,r)
            failures=failures and all(out[TARGETS[t]]=='UNRESOLVED' for t in targets)
    checks['failure_paths']=failures
    r.update(checks=checks,checks_passed=all(checks.values()),verdict=verdict(checks,r))
    return r


def serialize_report(report):
    # Preserve independent verdicts when one evidence block contains NaN/Inf.
    # Invalid values remain visible as null with explicit paths, never as zero.
    invalid=[]
    def clean(value,path):
        if isinstance(value,(float,np.floating)) and not np.isfinite(value):
            invalid.append(path);return None
        if isinstance(value,dict):return {key:clean(v,path+'.'+key) for key,v in value.items()}
        if isinstance(value,list):return [clean(v,path+'['+str(i)+']') for i,v in enumerate(value)]
        return value
    safe=clean(report,'report')
    if invalid:safe['nonfinite_evidence_paths']=invalid
    return json.dumps(safe,indent=1,allow_nan=False)+'\n'


def main(argv=None):
    p=argparse.ArgumentParser();p.add_argument('--output-dir','--output',dest='output',type=Path,required=True);args=p.parse_args(argv)
    try:
        r=finalize(run_probe());payload=serialize_report(r)
    except Exception as exc:
        r=dict(checks_passed=False,error=type(exc).__name__+': '+str(exc),verdict=verdict({}))
        payload=json.dumps(r,indent=1,allow_nan=False)+'\n'
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'probe.json').write_text(payload)
    lines=['# Nonlinear supported tensor completion','',f'Freeze: `{n.PUBLIC_PREREG}`.','']
    if 'checks' in r:lines+=['| Gate | Pass |','|---|---|']+[f'| {k} | {v} |' for k,v in r['checks'].items()]
    else:lines+=[r['error']]
    lines+=['','```json',json.dumps(r['verdict'],indent=2),'```','']
    (args.output/'probe.md').write_text('\n'.join(lines))
    return 0 if r['checks_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
