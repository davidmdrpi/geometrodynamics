"""Full-geometry verification and linear response of the degree-1 ESU support."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from geometrodynamics.waves import coupled_multiplet_response as cm
from geometrodynamics.waves import reciprocal_scalar_tt as rt

EPSILONS=(.002,.001,.0005)
ANGLES=((.83,1.07,.61),(1.21,.72,1.31),(.57,1.43,2.14))


def geometry_checks():
    rng=np.random.default_rng(cm.SEED)
    cases=[]
    # All five components, field and velocity zeros, generic spatial points.
    for phase in (0.,np.pi/4,np.pi/2,3*np.pi/4):
        for i,E in enumerate(rt.STF_BASIS):
            cases.append((f'basis_{i}',cm.CoupledSupport(),phase/2, .2*E,.31*E,.42*E,ANGLES[i%3]))
    zero=np.zeros((3,3))
    for name,b,v,d in (('pure_beta',.2*rt.STF_BASIS[2],zero,zero),
                       ('pure_velocity',zero,.31*rt.STF_BASIS[3],zero),
                       ('pure_acceleration',zero,zero,.42*rt.STF_BASIS[4])):
        cases.append((name,cm.CoupledSupport(),.31,b,v,d,ANGLES[0]))
    for a in (.7,1.,2.):
        for kap in (.4,1.):
            m=cm.CoupledSupport(a,kap,.17)
            for angles in ANGLES:
                b=rt.tensor(rng.normal(size=5))*.15
                v=rt.tensor(rng.normal(size=5))*.2/a
                d=rt.tensor(rng.normal(size=5))*.3/a**2
                cases.append(('generic',m,.31*a,b,v,d,angles))
    # An on-equation generic tensor checks every Einstein component jointly.
    m=cm.CoupledSupport();t=.23
    b=rt.tensor(rng.normal(size=5))*.15;v=rt.tensor(rng.normal(size=5))*.2
    f,fp,_,g=m.coefficients(t); d=-(fp*v+g*b)/f
    cases.append(('on_equation',m,t,b,v,d,ANGLES[1]))
    rows=[]
    for name,m,t,b,v,d,angles in cases:
        expectedG=d+8*b/m.radius**2
        expectedT=m.predicted_stress(t,b,v,d)
        expectedR=m.residual(t,b,v,d)
        scales=dict(G=max(1.,np.linalg.norm(expectedG)),T=max(1.,np.linalg.norm(expectedT)),
                    residual=max(1.,np.linalg.norm(expectedR)))
        steps=[]
        for eps in EPSILONS:
            response=cm.finite_response(m,t,b,v,d,eps,angles)
            G=response['einstein'];T=response['stress'];R=response['residual']
            constraints=np.r_[R[0,:],R[1:,0],np.trace(R[1:,1:])]
            steps.append(dict(epsilon=eps,
                G_absolute=float(np.linalg.norm(cm.stf(G[1:,1:])-expectedG)),
                T_absolute=float(np.linalg.norm(cm.stf(T[1:,1:])-expectedT)),
                residual_absolute=float(np.linalg.norm(cm.stf(R[1:,1:])-expectedR)),
                constraints_absolute=float(np.linalg.norm(constraints)),
                KG_absolute=float(np.linalg.norm(response['KG'])),
                scalar_curvature_absolute=float(abs(response['R'])),
                complete_Einstein_response=float(np.linalg.norm(R))))
        convergence=[]
        for field in ('G_absolute','T_absolute','residual_absolute','constraints_absolute','KG_absolute'):
            for i in (0,1):
                coarse,fine=steps[i][field],steps[i+1][field]
                if coarse>1e-10 and fine>1e-10:
                    convergence.append(dict(field=field,ratio=coarse/fine))
        rows.append(dict(case=name,radius=m.radius,kappa=m.kappa,time=t,phase=m.phase,
            angles=angles,scales=scales,steps=steps,convergence=convergence))
    backgrounds=[]
    for a in (.7,1.,2.):
        for kap in (.4,1.):
            m=cm.CoupledSupport(a,kap)
            for phase in (0.,np.pi/4,np.pi/2):
                t=phase*a/2
                r=cm.full_geometry(m,t,zero,zero,zero,0.,ANGLES[0])
                antipode=(np.pi-ANGLES[0][0],np.pi-ANGLES[0][1],ANGLES[0][2]+np.pi)
                opposite=cm.full_geometry(m,t,zero,zero,zero,0.,antipode)
                backgrounds.append(dict(radius=a,kappa=kap,phase=phase,
                    einstein_scaled=float(np.linalg.norm(r['residual'])*a*a),
                    KG_absolute=float(np.linalg.norm(r['KG'])),
                    odd_absolute=float(np.linalg.norm(r['phi']+opposite['phi'])),
                    density_error=float(abs(r['stress'][0,0]+3/(2*kap*a*a)))))
    return rows,backgrounds


def control_checks():
    m=cm.CoupledSupport();E=rt.STF_BASIS[0];zero=np.zeros((3,3))
    bare=m.residual(0.,E,zero,-8*E)
    # Pure velocity at a phase of nonzero f' isolates the sign, not stiffness.
    tau=np.pi/8;f,fp,_,g=m.coefficients(tau)
    wrong_sign_acc=fp*E/f
    sign=m.residual(tau,zero,E,wrong_sign_acc)
    frozen=m.residual(tau,E,zero,-8*E)
    times=np.linspace(0,np.pi,81)
    bare_history=[np.linalg.norm(m.residual(t,np.cos(np.sqrt(8)*t)*E,
                    -np.sqrt(8)*np.sin(np.sqrt(8)*t)*E,-8*np.cos(np.sqrt(8)*t)*E)) for t in times]
    return dict(bare_turning_residual=float(np.linalg.norm(bare)),expected_bare=1.5,
                reversed_fprime_residual=float(np.linalg.norm(sign)),
                dropped_metric_response_residual=float(np.linalg.norm(frozen)),
                bare_history_max=float(max(bare_history)))


def evolution_checks():
    model=cm.CoupledSupport();T=np.pi/2
    times,coarse=cm.evolve(model,rtol=1e-10,atol=1e-12)
    _,fine=cm.evolve(model,max_step=np.pi/200)
    _,normal=cm.evolve(model,normal_form=True,max_step=np.pi/200)
    _,long=cm.evolve(model,periods=20)
    J=np.array([[0.,1.],[-1.,0.]])
    mon=fine[-1]
    period_error=max(np.linalg.norm(long[i*100]-np.linalg.matrix_power(mon,i)) for i in range(21))
    scales=[]
    for a in (.7,1.,2.):
        for kap in (.4,1.):
            m=cm.CoupledSupport(a,kap)
            physical=cm.physical_evolve(m,a*T)
            scales.append(dict(radius=a,kappa=kap,error=float(np.linalg.norm(physical-mon))))
    eta0=cm.evolve(model,eta=0.)[1][-1]
    bare_exact=expm(np.array([[0.,1.],[-8.,0.]])*T)
    rng=np.random.default_rng(cm.SEED+1)
    initial=rng.normal(size=(2,5))
    five=np.einsum('tij,ja->tia',fine,initial)
    # Independent normal-form evolution for the same five initial components.
    normal_five=np.einsum('tij,ja->tia',normal,initial)
    phases=[]
    for phase in (.31,.77):
        p=cm.evolve(cm.CoupledSupport(phase=phase))[1][-1]
        phases.append(dict(phase=phase,trace=float(np.trace(p)),trace_difference=float(abs(np.trace(p)-np.trace(mon)))))
    traces=[float(np.trace(coarse[-1])),float(np.trace(mon))]
    eigen=np.linalg.eigvals(mon)
    return dict(evidence=dict(maps=[coarse[-1].tolist(),mon.tolist()],traces=traces),
        period_tau=T,step_caps=[np.pi/100,np.pi/200],monodromy=mon.tolist(),trace=traces[-1],
        eigenvalues=[dict(real=float(z.real),imag=float(z.imag),modulus=float(abs(z))) for z in eigen],
        classification=cm.classify_period(traces),
        determinant=float(np.linalg.det(mon)),
        refinement=float(np.max(np.linalg.norm(coarse-fine,axis=(1,2)))),
        normal_form_error=float(np.max(np.linalg.norm(normal-fine,axis=(1,2)))),
        symplectic_error=float(max(np.linalg.norm(M.T@J@M-J) for M in long)),
        twenty_period_error=float(period_error),
        twenty_period_map_norm=float(np.linalg.norm(long[-1])),
        canonical_map_norm_max=float(np.max(np.linalg.norm(long,axis=(1,2)))),
        five_component_error=float(np.max(abs(five-normal_five))),physical_time=scales,
        bare_control_error=float(np.linalg.norm(eta0-bare_exact)),
        bare_trace=float(np.trace(bare_exact)),phase_controls=phases,
        quasifrequency='NOT_ASSIGNED: modulo-4 ambiguity; trace is compared directly')


def failure_controls(evidence):
    good=dict.fromkeys(cm.REQUIRED_CHECKS,True);rows=[]
    for key in cm.REQUIRED_CHECKS:
        for mode in ('missing','failed'):
            checks=good.copy()
            if mode=='missing':del checks[key]
            else:checks[key]=False
            v=cm.verdict(checks,evidence)
            rows.append(dict(key=key,mode=mode,passed=all(v[k]=='UNRESOLVED' for k in cm.PHYSICAL)
                             and key in v['failed_checks']))
    return rows


def run_probe():
    exact=cm.exact_certificate()
    geometry,backgrounds=geometry_checks()
    controls=control_checks()
    # Only integrate after independent curvature and stress have agreed.
    if not all(r['steps'][-1]['G_absolute']/r['scales']['G']<2e-4
               and r['steps'][-1]['T_absolute']/r['scales']['T']<2e-4 for r in geometry):
        raise ArithmeticError('independent response variation failed before evolution')
    evolution=evolution_checks();failures=failure_controls(evolution['evidence'])
    checks=dict(
        background_and_parity=exact['parity_degree']==1 and exact['components']==4
            and all(r['einstein_scaled']<1e-10 and r['odd_absolute']<1e-10 for r in backgrounds),
        scalar_sector_closure=all(exact['residuals'][f'KG_{i}']=='0' for i in range(4))
            and all(r['steps'][-1]['KG_absolute']<2e-4 for r in geometry),
        quadratic_action=exact['all_zero'],
        independent_curvature_stress=all(r['steps'][-1]['G_absolute']/r['scales']['G']<2e-4
            and r['steps'][-1]['T_absolute']/r['scales']['T']<2e-4
            and r['steps'][-1]['residual_absolute']/r['scales']['residual']<2e-4
            and all(3.5<c['ratio']<4.5 for c in r['convergence']) for r in geometry),
        constraint_completion=all(r['steps'][-1]['constraints_absolute']<2e-4
            and r['steps'][-1]['scalar_curvature_absolute']<2e-4 for r in geometry)
            and next(r for r in geometry if r['case']=='on_equation')['steps'][-1]['complete_Einstein_response']<2e-4,
        bare_frequency_control=abs(controls['bare_turning_residual']-1.5)<1e-12
            and controls['reversed_fprime_residual']>.1 and controls['dropped_metric_response_residual']>.1,
        canonical_propagation=max(evolution['normal_form_error'],evolution['symplectic_error'],
            evolution['five_component_error'],evolution['bare_control_error'],
            max(r['error'] for r in evolution['physical_time']),
            max(r['trace_difference'] for r in evolution['phase_controls']))<1e-8,
        period_map_convergence=max(evolution['refinement'],abs(evolution['determinant']-1),
            evolution['twenty_period_error'])<1e-8 and evolution['classification']!='UNRESOLVED',
        scope_and_order=cm.SCOPE==dict(general_scalar_vector_response='NOT_DERIVED',
            full_dynamical_stability='NOT_ESTABLISHED',nonlinear_persistence='NOT_ESTABLISHED',
            preparation_selection='NOT_DERIVED',Phi_selection='NOT_DERIVED',causality_gate='OPEN'),
        failure_paths=all(r['passed'] for r in failures),
    )
    return dict(prereg=cm.PUBLIC_PREREG,baseline=cm.BASELINE,seed=cm.SEED,
                exact=exact,geometry=geometry,backgrounds=backgrounds,controls=controls,
                evolution=evolution,failure_controls=failures,
                checks={k:bool(v) for k,v in checks.items()},period_evidence=evolution['evidence'])


def finalize(report):
    report['verdict']=cm.verdict(report.get('checks',{}),report.get('period_evidence'))
    report['checks_passed']=not report['verdict']['failed_checks']
    return report


def render(report):
    lines=['# Coupled multiplet–metric tensor response','',
           'Four degree-1 support fields, linear homogeneous TT sector, proper time.','',
           '| Gate | Passed |','|---|---|']
    lines += [f"| {k} | {report.get('checks',{}).get(k,False)} |" for k in cm.REQUIRED_CHECKS]
    lines += ['','| Verdict | Result |','|---|---|']
    lines += [f'| {k} | {v} |' for k,v in report['verdict'].items()]
    if report['checks_passed']:
        e=report['evolution']
        lines += ['',f"Period-map trace: {e['trace']:.12f}; determinant: {e['determinant']:.12f}.",
                  f"Bare-oscillator trace: {e['bare_trace']:.12f}.",
                  'Numerical ellipticity is a result for this linear tensor sector, not full-system stability.']
    else:lines += ['','Failed or missing verification leaves the physical conclusions UNRESOLVED.']
    return '\n'.join(lines)+'\n'


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args(argv)
    try:
        report=run_probe()
        required={'exact','geometry','backgrounds','controls','evolution','failure_controls','checks','period_evidence'}
        if not isinstance(report,dict) or not required.issubset(report):
            raise ValueError('missing required probe evidence')
        json.dumps(report,allow_nan=False)  # reject nonfinite evidence before writing
        report=finalize(report)
        summary=render(report)
    except Exception as exc:
        report=finalize(dict(checks={},error=f'{type(exc).__name__}: {exc}'))
        summary=render(report)
    if args.output_dir:
        args.output_dir.mkdir(parents=True,exist_ok=True)
        (args.output_dir/'probe.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        (args.output_dir/'probe.md').write_text(summary)
    print(summary)
    return 0 if report['checks_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
