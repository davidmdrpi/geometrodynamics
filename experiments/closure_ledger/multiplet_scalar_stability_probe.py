"""Reproduce low scalar mode stability and independent full-field controls."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, null_space

from geometrodynamics.waves import multiplet_scalar_stability as m

EPSILONS=(.002,.001,.0005)
TOLERANCES=((1e-10,1e-12),(1e-12,1e-14))
ANGLES=((.83,1.07,.61),(1.13,.79,1.29),(.61,1.31,.94))


def variation(fn,expected):
    steps=[];ratios=[]
    for eps in EPSILONS:
        plus=fn(eps);minus=fn(-eps)
        errors={k:float(np.linalg.norm((plus[k]-minus[k])/(2*eps)-v)/max(1.,np.linalg.norm(v)))
                for k,v in expected.items()}
        steps.append(dict(epsilon=eps,errors=errors))
    for i in range(2):
        for key,value in steps[i+1]['errors'].items():
            if value>1e-10:ratios.append(steps[i]['errors'][key]/value)
    return dict(steps=steps,ratios=ratios)


def run_probe():
    rng=np.random.default_rng(m.SEED)
    exact=m.exact_certificate()
    homogeneous=[];dipole=[];continuations=[];gauges=[];parity=[];cover=[]
    # Independent full curvature: all coefficients, not only solutions.
    for i in range(18):
        a=(.7,1.,2.)[i%3];kap=(.4,1.)[i%2];phase=(0.,.31,math.pi/4)[(i//3)%3]
        tau=(0.,math.pi/4,.37)[i%3];angles=ANGLES[i%3]
        r=rng.normal(size=3)*.15;z=rng.normal(size=3)*.15
        if i<3:r=np.eye(3)[i]*.2;z=np.zeros(3)
        if 3<=i<6:z=np.eye(3)[i-3]*.2;r=np.zeros(3)
        if i>=12:r[2]=2*r[0];z=.12*m.field_jets(tau,phase)[1:4]
        q0=a/math.sqrt(kap)*m.field_jets(tau,phase)[:3]
        fn=lambda eps:m.frw_geometry(np.array([a,0.,0.])+eps*a*r,q0+eps*a/math.sqrt(kap)*z,a,kap,angles)
        row=variation(fn,m.expected_homogeneous(tau,r,z,a,kap,phase,angles))
        row.update(case=i,a=a,kappa=kap,phase=phase,tau=tau,on_equation=i>=12)
        homogeneous.append(row)
        # Arbitrary off-shell lapse, matter and accelerations, including zeros.
        u=rng.normal(size=3)*.15;v=rng.normal(size=3)*.15;lapse=rng.normal(size=3)*.12
        direction=rng.normal(size=4);direction/=np.linalg.norm(direction)
        if i<6:
            u=np.zeros(3);v=np.zeros(3)
            (u if i<3 else v)[i%3]=.2;lapse=np.zeros(3)
        kwargs=dict(a=a,kap=kap,phase=phase,angles=angles,direction=direction)
        expected=m.expected_dipole(tau,u,v,lapse,**kwargs)
        row=variation(lambda eps:m.dipole_geometry(tau,u,v,lapse,epsilon=eps,**kwargs),expected)
        row.update(case=i,a=a,kappa=kap,phase=phase,tau=tau)
        dipole.append(row)
        gauge=m.gauge_jets(tau,phase)
        zero={k:np.zeros_like(value) for k,value in expected.items() if k in ('residual','KG')}
        gauge_row=variation(lambda eps:m.dipole_geometry(tau,epsilon=eps,**gauge,**kwargs),zero)
        gauges.append(gauge_row)
        x=m.cm.coframe_jets(angles)[0]
        # Test the perturbation itself, including nonzero norm: parity is not
        # inferred from a background-only cancellation or a zero field phase.
        Y=direction@x;pert=(u[0]-v[0])*x*Y+v[0]*direction
        anti=(u[0]-v[0])*(-x)*(-Y)+v[0]*direction
        parity.append(dict(even_error=float(np.linalg.norm(pert-anti)),norm=float(np.linalg.norm(pert)),
                           metric_even_defect=2*abs(float(Y))))
        # Fixed-energy, finite, expanding AND contracting, under/over-density.
        for sign in (-1.,1.):
            for eps in (-.01,.01):
                Aj=m.exact_scale(tau,eps,a,sign)
                result=m.frw_geometry(Aj,q0,a,kap,angles)
                prediction=m.expected_frw(Aj,q0,a,kap,angles)
                errors={k:float(np.linalg.norm(result[k]-prediction[k])/max(1.,np.linalg.norm(prediction[k])))
                        for k in prediction}
                energy=(q0[1]**2+4*q0[0]**2)/2
                constraint=Aj[1]**2+Aj[0]**2-Aj[0]**4/(2*a*a)-kap*energy/3
                continuations.append(dict(a=a,kappa=kap,phase=phase,tau=tau,epsilon=eps,sign=sign,
                    full_residual=float(a*a*np.linalg.norm(result['residual'])),
                    KG=float(a*a*math.sqrt(kap)*np.linalg.norm(result['KG'])),
                    constraint=float(constraint/(a*a)),
                    f=float(1-kap*(q0[0]/Aj[0])**2/6),errors=errors))
    # Check arbitrary finite off-shell FRW jets separately from exact histories.
    frw_off=[]
    for _ in range(10):
        Aj=np.r_[rng.uniform(.8,1.2),rng.normal(size=2)*.2]
        qj=rng.normal(size=3)
        full=m.frw_geometry(Aj,qj);pred=m.expected_frw(Aj,qj)
        frw_off.append({k:float(np.linalg.norm(full[k]-v)/max(1.,np.linalg.norm(v))) for k,v in pred.items()})
    # Closed constrained dipole subspace, all four Y directions and both ICs.
    for phase in (0.,.31,math.pi/4):
        basis=null_space(m.constraints(0.,phase))
        for tau in (0.,math.pi/4,.37,math.pi):
            for j in range(2):
                y=expm(m.DIPOLE*tau)@basis[:,j];yd=m.DIPOLE@y
                u=[y[0],y[1],yd[1]];v=[y[2],y[3],yd[3]]
                for d in np.eye(4):
                    row=variation(lambda eps:m.dipole_geometry(tau,u,v,phase=phase,direction=d,epsilon=eps),
                                  dict(residual=np.zeros((4,4)),KG=np.zeros(4)))
                    row.update(phase=phase,tau=tau,initial=j,direction=d.tolist(),
                        constraint=float(np.linalg.norm(m.constraints(tau,phase)@y)),
                        stress_coefficients=m.dipole_coefficients(tau,u,v,phase=phase))
                    cover.append(row)
    maps={key:[] for key in ('homogeneous','dipole')};map_checks=[]
    for rtol,atol in TOLERANCES:
        h=m.period_map(m.HOMOGENEOUS,rtol=rtol,atol=atol)
        d,leak=m.constrained_dipole_map(rtol=rtol,atol=atol)
        maps['homogeneous'].append(h.tolist());maps['dipole'].append(d.tolist())
        for phase in (0.,.31,math.pi/4):
            dm,leak=m.constrained_dipole_map(phase,rtol=rtol,atol=atol)
            map_checks.append(dict(kind='dipole',phase=phase,error=float(np.linalg.norm(dm+np.eye(2))),leakage=leak))
        for a in (.7,1.,2.):
            # Independent proper-time ODE (r, dr/dt); convert at both ends.
            B=np.array([[0.,1.],[2/(a*a),0.]])
            hp=m.period_map(B,period=a*math.pi,rtol=rtol,atol=atol)
            conversion=np.diag([1.,a]);hp=conversion@hp@np.linalg.inv(conversion)
            map_checks.append(dict(kind='proper_clock',a=a,error=float(np.linalg.norm(hp-h)/np.linalg.norm(h)),leakage=0.))
    half=m.period_map(m.HOMOGENEOUS,period=math.pi/2)
    # Two clocks on the same exact finite family; conformal scalar phase shifts
    # in proper time. Check field as well as scale and volume.
    clocks=[]
    for phase in (0.,.31,math.pi/4):
        tau=.4;growth=math.exp(math.sqrt(2)*tau);integral=(growth-1)/math.sqrt(2)
        p,pd,*_=m.field_jets(tau,phase)
        expected=np.array([growth,-p*growth-pd*integral,3*growth])
        for eps in EPSILONS:
            pair=[]
            for signed in (eps,-eps):
                sol=solve_ivp(lambda t,y:np.array([(y[0]**2-1)/(math.sqrt(2)*y[0]),1/y[0]]),
                              (0,tau),[1+signed,0.],rtol=1e-12,atol=1e-14,method='DOP853',max_step=.01)
                if not sol.success:raise ArithmeticError('proper clock continuation failed')
                A,eta=sol.y[:,-1];pactual=m.field_jets(eta,phase)[0]/A
                pair.append(np.array([A,pactual,A**3]))
            estimate=(pair[0]-pair[1])/(2*eps)
            clocks.append(dict(phase=phase,epsilon=eps,error=float(np.linalg.norm(estimate-expected)/np.linalg.norm(expected)),
                               omitted_clock_term=abs(float(pd*integral))))
    # Controls come from independent residuals, not verdict strings.
    q=m.field_jets(.37)[:3]
    wrong=m.frw_geometry([1.,0.,.05],q)
    energy_control=m.expected_homogeneous(.37,[0.,0.,0.],q)
    gauge=m.gauge_jets(.37)
    plus=m.dipole_geometry(.37,gauge['u'],gauge['v'],epsilon=.0005)
    minus=m.dipole_geometry(.37,gauge['u'],gauge['v'],epsilon=-.0005)
    controls=dict(wrong_acceleration=float(np.linalg.norm(wrong['residual'])),
                  nonzero_linear_energy=float(abs(energy_control['residual'][0,0])),
                  gauge_without_metric=float(np.linalg.norm((plus['residual']-minus['residual'])/.001)),
                  wrong_dipole_map_rejected=not m.valid_period_evidence({**maps,'dipole':[np.eye(2).tolist()]*2}))
    return dict(prereg=m.PUBLIC_PREREG,baseline=m.BASELINE,seed=m.SEED,exact=exact,
                homogeneous=homogeneous,dipole=dipole,continuations=continuations,frw_off_shell=frw_off,
                gauges=gauges,parity=parity,cover=cover,period_evidence=maps,map_checks=map_checks,
                half_period=half.tolist(),clocks=clocks,controls=controls)


def finalize(report):
    # Reject malformed/nonfinite input before deriving any affirmative gate.
    json.dumps(report,allow_nan=False)
    def variations_ok(rows):
        return bool(rows) and all(len(r['steps'])==3 and max(r['steps'][-1]['errors'].values())<2e-4
                 and all(3.5<x<4.5 for x in r['ratios']) for r in rows)
    checks=dict(
        field_reduction=report['exact']['all_zero'] is True and set(report['exact']['residuals'])==set(m.EXACT_IDENTITIES)
                        and all(v=='0' for v in report['exact']['residuals'].values())
                        and bool(report['frw_off_shell']) and max(max(r.values()) for r in report['frw_off_shell'])<1e-9,
        homogeneous_constraints=bool(report['homogeneous']) and variations_ok(report['homogeneous'])
                        and abs(report['controls']['nonzero_linear_energy'])>1e-3,
        exact_continuation=bool(report['continuations']) and all(max(r['full_residual'],r['KG'],abs(r['constraint']))<1e-9
                               and r['f']>0 for r in report['continuations']),
        clock_and_gauge=variations_ok(report['gauges']) and bool(report['clocks'])
                        and max(r['error'] for r in report['clocks'] if r['epsilon']==EPSILONS[-1])<2e-4,
        dipole_reduction=variations_ok(report['cover']) and all(r['constraint']<1e-9 and
                          max(abs(v) for v in r['stress_coefficients'].values())<1e-9 for r in report['cover']),
        dipole_parity=bool(report['parity']) and all(r['even_error']<1e-12 for r in report['parity'])
                        and max(r['norm'] for r in report['parity'])>1e-3,
        independent_geometry=variations_ok(report['dipole']) and variations_ok(report['homogeneous']),
        period_maps=m.valid_period_evidence(report['period_evidence']) and bool(report['map_checks'])
                      and max(max(r['error'],r['leakage']) for r in report['map_checks'])<1e-8,
        negative_controls=all(report['controls'][k]>1e-3 for k in ('wrong_acceleration','nonzero_linear_energy','gauge_without_metric'))
                          and report['controls']['wrong_dipole_map_rejected'] is True,
        failure_paths=True)
    checks={k:bool(v) for k,v in checks.items()}
    passing=dict.fromkeys(m.REQUIRED_CHECKS,True)
    affirmative=m.verdict(passing,report['period_evidence'])['homogeneous_physical_block']=='HYPERBOLIC_GROWING_MODE'
    failures=[]
    for key in m.REQUIRED_CHECKS:
        for missing in (False,True):
            bad=passing.copy()
            if missing:bad.pop(key)
            else:bad[key]=False
            failures.append(m.verdict(bad,report['period_evidence'])['homogeneous_physical_block']=='UNRESOLVED')
    checks['failure_paths']=affirmative and all(failures)
    report['checks']=checks;report['checks_passed']=all(checks.values())
    report['verdict']=m.verdict(checks,report['period_evidence'])
    # Standalone homogeneous evidence remains visible if a separate dipole gate fails.
    report['homogeneous_evidence_passed']=(all(checks[k] for k in
        ('field_reduction','homogeneous_constraints','exact_continuation'))
        and bool(report['clocks']) and max(r['error'] for r in report['clocks'] if r['epsilon']==EPSILONS[-1])<2e-4
        and m.valid_map_pair(report['period_evidence']['homogeneous'],expm(m.HOMOGENEOUS*math.pi)))
    return report


def render(report):
    if 'checks' not in report:return '# Low scalar modes\n\nUNRESOLVED: '+report['error']+'\n'
    h=np.asarray(report['period_evidence']['homogeneous'][-1])
    d=np.asarray(report['period_evidence']['dipole'][-1])
    lines=['# Low scalar modes of the four-field ESU','',f"Freeze: `{m.PUBLIC_PREREG}`.",'',
           '| Quantity | Value |','|---|---|',
           f'| Homogeneous full-period trace | {np.trace(h):.12f} |',
           f'| Homogeneous multipliers | {math.exp(math.sqrt(2)*math.pi):.12f}, {math.exp(-math.sqrt(2)*math.pi):.12f} |',
           f'| Homogeneous determinant | {np.linalg.det(h):.12f} |',
           f'| Constrained cover dipole: norm(M + I) | {np.linalg.norm(d+np.eye(2)):.3e} |',
           f"| Exact FRW: maximum Einstein residual | {max(r['full_residual'] for r in report['continuations']):.3e} |",'',
           '| Gate | Pass |','|---|---|']
    lines += [f'| {k} | {v} |' for k,v in report['checks'].items()]
    lines += ['', '```json',json.dumps(report['verdict'],indent=2),'```','',
              'The cover dipole is excluded only under the stated antipodal restrictions.',
              'This linear instability does not specify a measure on viable preparations.','']
    return '\n'.join(lines)


def main(argv=None):
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        report=finalize(run_probe());text=render(report);payload=json.dumps(report,indent=2,allow_nan=False)+'\n'
    except Exception as exc:
        report=dict(checks_passed=False,error=type(exc).__name__+': '+str(exc),verdict=m.verdict({}))
        text=render(report);payload=json.dumps(report,indent=2,allow_nan=False)+'\n'
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'probe.json').write_text(payload);(args.output/'probe.md').write_text(text)
    return 0 if report['checks_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
