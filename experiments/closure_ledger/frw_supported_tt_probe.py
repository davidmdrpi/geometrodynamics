"""Frozen full-field checks for the supported homogeneous FRW tensor block."""
import argparse
from dataclasses import asdict
import itertools
import json
import math
from pathlib import Path

import numpy as np
from geometrodynamics.waves import frw_supported_tt as f
from geometrodynamics.waves import reciprocal_scalar_tt as rt

EPSILONS=(.002,.001,.0005)
ANGLES=((.83,1.07,.61),(1.21,.72,1.31),(.57,1.43,2.14))


def finite_check(model,eta,b,v,d,angles):
    expected=f.expected_response(model,eta,b,v,d)
    steps=[];ratios=[]
    for eps in EPSILONS:
        values={clock:[f.full_geometry(model,eta,b,v,d,sgn*eps,clock,angles) for sgn in (1,-1)]
                for clock in ('conformal','proper')}
        errors={}
        for clock,pair in values.items():
            for key,target in expected.items():
                actual=(pair[0][key]-pair[1][key])/(2*eps)
                errors[clock+'_'+key]=float(np.linalg.norm(actual-target)/max(1.,np.linalg.norm(target)))
        covariance=max(np.linalg.norm(values['conformal'][i][key]-values['proper'][i][key])/
                       max(1.,np.linalg.norm(values['conformal'][i][key]))
                       for i in (0,1) for key in expected)
        residual=(values['conformal'][0]['residual']-values['conformal'][1]['residual'])/(2*eps)
        constraint=float(np.linalg.norm(np.r_[residual[0],residual[1:,0],np.trace(residual[1:,1:])]))
        steps.append(dict(epsilon=eps,errors=errors,clock_error=float(covariance),constraint=constraint))
    for i in (0,1):
        for key,value in steps[i+1]['errors'].items():
            if value>1e-10:ratios.append(steps[i]['errors'][key]/value)
    return dict(steps=steps,ratios=ratios,expected_residual_norm=float(np.linalg.norm(expected['residual'])))


def run_probe():
    rng=np.random.default_rng(f.SEED);zero=np.zeros((3,3));cases=[]
    # Explicit coverage of each of the five components and three time jets.
    for j,E in enumerate(rt.STF_BASIS):
        for jet in range(3):
            index=3*j+jet
            model=f.FRWSupport(radius=(.7,1.,2.)[j%3],kappa=(.4,1.)[j%2],
                               departure=(.1,.3,-.1)[jet],sign=(-1,1)[j%2],phase=(0.,.31,math.pi/4)[j%3])
            jets=[zero.copy() for _ in range(3)];jets[jet]=.2*E
            cases.append((f'basis_{j}_jet_{jet}',model,(0.,.37,math.pi/4,.8)[index%4],*jets,ANGLES[j%3]))
    for dep,sgn,phase in itertools.product((.1,.3,-.1),(-1,1),(0.,.31,math.pi/4)):
        model=f.FRWSupport(departure=dep,sign=sgn,phase=phase)
        b,v,d=[rt.tensor(rng.normal(size=5))*.15 for _ in range(3)]
        cases.append(('noncommuting',model,.37,b,v,d,ANGLES[len(cases)%3]))
        cases.append(('on_equation',model,.8,b,v,model.acceleration(.8,b,v),ANGLES[len(cases)%3]))
    for eta in (0.,math.pi/4):
        cases.append(('conformal_velocity_zero' if eta==0 else 'field_zero',f.FRWSupport(departure=.3),eta,
                      .2*rt.STF_BASIS[0],.3*rt.STF_BASIS[1],.2*rt.STF_BASIS[2],ANGLES[0]))
    # q'=0 does not mean (q/A)'=0. Additional deterministic phases meet the
    # frozen physical-field velocity-zero check without moving the time grid.
    for dep in (.3,-.1):
        provisional=f.FRWSupport(departure=dep);A,Ap,*_=provisional.jets(.37)
        phase=math.atan(-Ap/(2*A))-.74
        model=f.FRWSupport(departure=dep,phase=phase)
        cases.append(('physical_velocity_zero',model,.37,.2*rt.STF_BASIS[0],.3*rt.STF_BASIS[1],
                      .2*rt.STF_BASIS[2],ANGLES[1]))
    geometry=[];backgrounds=[]
    for name,model,eta,b,v,d,angles in cases:
        row=finite_check(model,eta,b,v,d,angles)
        row.update(case=name,model=asdict(model),eta=eta,angles=list(angles),
                   commutator=float(np.linalg.norm(b@v-v@b)))
        geometry.append(row)
        A,Ap,App,q,qp,qpp=model.jets(eta)
        bg=f.full_geometry(model,eta,zero,zero,zero,angles=angles)
        rho=3*model.radius**2/(2*model.kappa*A**4)
        antiangles=(math.pi-angles[0],math.pi-angles[1],angles[2]+math.pi)
        anti=f.full_geometry(model,eta,zero,zero,zero,angles=antiangles)
        backgrounds.append(dict(model=asdict(model),eta=eta,A=A,
            residual=float(model.radius**2*np.linalg.norm(bg['residual'])),
            KG=float(model.radius**2*math.sqrt(model.kappa)*np.linalg.norm(bg['KG'])),
            density_error=float(abs(bg['stress'][0,0]+rho)/max(1.,rho)),
            parity_error=float(math.sqrt(model.kappa)*np.linalg.norm(bg['phi']+anti['phi'])),
            q=q,qp=qp,field_velocity=float(qp/A-q*Ap/(A*A))))
    maps=[];composition=[]
    for dep,sgn,phase,a,kap in itertools.product((.1,.3,-.1),(-1,1),(0.,.31,math.pi/4),(.7,1.,2.),(.4,1.)):
        model=f.FRWSupport(a,kap,dep,sgn,phase)
        coarse=f.evolve(model,rtol=1e-10,atol=1e-12);fine=f.evolve(model)
        normal=f.evolve(model,normal=True);proper=f.proper_evolve(model)
        maps.append(dict(model=asdict(model),maps=[x.tolist() for x in (coarse,fine,normal,proper['map'])],
                         eta_error=proper['eta_error'],scale_error=proper['scale_error'],
                         kinetic_bound=model.interval_bound(),proper_elapsed=proper['elapsed']))
        if a==1 and kap==1:
            joined=f.evolve(model,.4,.8)@f.evolve(model,0.,.4)
            inverse=f.evolve(model,.8,0.)
            composition.append(dict(model=asdict(model),
                composition_error=float(np.linalg.norm(joined-fine)/max(1.,np.linalg.norm(fine))),
                inverse_error=float(np.linalg.norm(inverse@fine-np.eye(2)))))
    static=[];limits=[]
    for phase in (0.,.31,math.pi/4):
        model=f.FRWSupport(departure=0.,phase=phase)
        reference=f.cm.CoupledSupport(phase=phase)
        _,trajectory=f.cm.evolve(reference)
        actual=f.evolve(model,end=math.pi/2)
        coefficient_errors=[]
        for eta in (0.,.37,math.pi/4,.8):
            coefficient_errors.append(float(np.linalg.norm(np.array(model.coefficients(eta))-reference.coefficients(eta))))
        static.append(dict(phase=phase,map_error=float(np.linalg.norm(actual-trajectory[-1])),
                           coefficient_error=max(coefficient_errors),trace=float(np.trace(actual))))
        endpoint=f.evolve(model);errors=[]
        for dep in (1e-2,1e-3,1e-4):
            matrix=f.evolve(f.FRWSupport(departure=dep,phase=phase))
            errors.append(float(np.linalg.norm(matrix-endpoint)))
        limits.append(dict(phase=phase,departures=[1e-2,1e-3,1e-4],errors=errors))
    # Controls use the independent full geometry on wrongly evolved tensor data.
    model=f.FRWSupport(departure=.3,phase=.31);eta=.37
    b=.2*rt.STF_BASIS[0];v=.3*rt.STF_BASIS[1]
    A,Ap,App,q,qp,qpp=model.jets(eta);m,mp,mpp,k=model.coefficients(eta)
    correct=model.acceleration(eta,b,v)
    bad=dict(bare_frw=-2*Ap*v/A-8*b,wrong_mass_derivative=(mp*v-k*b)/m,
             missing_expansion_in_Qprime=-((mp-model.kappa*(Ap/A)*q*q/(3*model.radius**2))*v+k*b)/m,
             wrong_acceleration=1.05*correct)
    controls={}
    for name,acc in bad.items():
        plus=f.full_geometry(model,eta,b,v,acc,.0005)
        minus=f.full_geometry(model,eta,b,v,acc,-.0005)
        controls[name]=float(np.linalg.norm((plus['residual']-minus['residual'])/.001))
    controls['bare_normal_potential_difference']=abs(float(model.normal_potential(eta)-(8-App/A)))
    exact=f.exact_certificate();controls['linear_volume_lambda']=exact['bad_linear_lambda_value']
    # No results at or across a pole or a zero/negative kinetic coefficient.
    rejected=0
    for operation in (lambda:f.FRWSupport(departure=-1),lambda:f.FRWSupport(departure=.3).jets(3.),
                      lambda:f.FRWSupport(departure=-.9).coefficients(0.)):
        try:operation()
        except ValueError:rejected+=1
    return dict(prereg=f.PUBLIC_PREREG,baseline=f.BASELINE,seed=f.SEED,exact=exact,
                geometry=geometry,backgrounds=backgrounds,map_evidence=maps,composition=composition,
                static=static,static_limit=limits,controls=controls,domain_rejections=rejected)


def finalize(report):
    json.dumps(report,allow_nan=False)
    geometry=report['geometry'];backgrounds=report['backgrounds'];exact=report['exact']
    def geometry_ok(keys):
        return bool(geometry) and all(max(r['steps'][-1]['errors'][k] for k in keys)<2e-4 for r in geometry)
    step_keys=[clock+'_'+key for clock in ('conformal','proper') for key in ('einstein','stress','residual','KG','R')]
    checks=dict(
        exact_background=bool(backgrounds) and all(max(r['residual'],r['KG'],r['density_error'],r['parity_error'])<1e-9 for r in backgrounds),
        action_derivation=exact['all_zero'] is True and set(exact['residuals'])==set(f.IDENTITIES)
                          and set(exact['residuals'].values())=={'0'} and exact['bad_linear_lambda_identity']=='0',
        full_stress_response=geometry_ok(['conformal_stress','proper_stress']),
        scalar_constraint_closure=geometry_ok(['conformal_KG','proper_KG','conformal_R','proper_R'])
             and all(r['steps'][-1]['constraint']<2e-4 for r in geometry)
             and any(r['case']=='on_equation' for r in geometry)
             and all(r['expected_residual_norm']<1e-10 and r['steps'][-1]['errors']['conformal_residual']<2e-4
                     for r in geometry if r['case']=='on_equation'),
        independent_geometry=geometry_ok(step_keys) and all(3.5<x<4.5 for r in geometry for x in r['ratios']),
        two_clocks=bool(geometry) and all(s['clock_error']<1e-9 for r in geometry for s in r['steps'])
                   and bool(report['map_evidence']) and all(r['eta_error']<1e-9 and r['scale_error']<1e-9 for r in report['map_evidence']),
        static_limit=bool(report['static']) and all(max(r['map_error'],r['coefficient_error'])<1e-8 for r in report['static'])
               and len(report['static_limit'])==3 and all(r['errors'][2]<r['errors'][1]<r['errors'][0] for r in report['static_limit']),
        canonical_transport=f.valid_maps(report['map_evidence']) and bool(report['composition'])
               and all(max(r['composition_error'],r['inverse_error'])<1e-8 for r in report['composition'])
               and report['domain_rejections']==3,
        negative_controls=all(report['controls'][k]>1e-5 for k in ('bare_frw','wrong_mass_derivative',
                  'missing_expansion_in_Qprime','wrong_acceleration','bare_normal_potential_difference','linear_volume_lambda')),
        failure_paths=True)
    checks={k:bool(v) for k,v in checks.items()}
    good=dict.fromkeys(f.REQUIRED_CHECKS,True)
    passed=f.verdict(good,report['map_evidence'])['supported_operator']=='FULLY_SUPPORTED_FRW_EQUATION_VERIFIED'
    for key in f.REQUIRED_CHECKS:
        for missing in (False,True):
            bad=good.copy()
            if missing:bad.pop(key)
            else:bad[key]=False
            passed=passed and f.verdict(bad,report['map_evidence'])['supported_operator']=='UNRESOLVED'
    checks['failure_paths']=passed
    report.update(checks=checks,checks_passed=all(checks.values()),verdict=f.verdict(checks,report['map_evidence']))
    return report


def render(report):
    if 'checks' not in report:return '# Supported FRW tensor response\n\nUNRESOLVED: '+report['error']+'\n'
    lines=['# Supported FRW tensor response','',f"Freeze: `{f.PUBLIC_PREREG}`.",'',
       '`(M beta\')\' + K beta = 0`, `M=A^2/kappa-q^2/6`, `K=8A^2/kappa+2q^2/3`.','',
       f"Independent full-geometry cases: {len(report['geometry'])}; transport backgrounds: {len(report['map_evidence'])}.",'',
       '| Gate | Pass |','|---|---|']
    lines.extend(f'| {k} | {v} |' for k,v in report['checks'].items())
    lines+=['','```json',json.dumps(report['verdict'],indent=2),'```','',
            'These are finite-interval classical maps on a nonperiodic background, not Floquet or quantum verdicts.','']
    return '\n'.join(lines)


def main(argv=None):
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir','--output',dest='output',type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        report=finalize(run_probe());text=render(report);payload=json.dumps(report,indent=2,allow_nan=False)+'\n'
    except Exception as exc:
        report=dict(checks_passed=False,error=type(exc).__name__+': '+str(exc),verdict=f.verdict({}))
        text=render(report);payload=json.dumps(report,indent=2,allow_nan=False)+'\n'
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'probe.json').write_text(payload);(args.output/'probe.md').write_text(text)
    return 0 if report['checks_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
