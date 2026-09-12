"""All-phase certificate, named action tests and complete asymptotic data."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from geometrodynamics.waves import tt_turning_transport as t
from geometrodynamics.waves import frw_supported_tt as f


def relative(a,b):return float(np.linalg.norm(np.asarray(a)-b)/max(1.,np.linalg.norm(b)))


def action_identity():
    import sympy as s
    b,p,m,k,mp,kp=s.symbols('b p m k mp kp',positive=True)
    c=s.sqrt(m*k);a=(p*p/c+c*b*b)/2
    derivative=s.diff(a,b)*p/m-s.diff(a,p)*k*b+s.diff(a,m)*mp+s.diff(a,k)*kp
    cp=c*(mp/m+kp/k)/2
    return str(s.simplify(derivative-cp*(b*b-p*p/(c*c))/2))


def run_probe():
    rng=np.random.default_rng(t.SEED)
    phases=list(np.arange(64)*math.pi/64)+list(rng.uniform(0,math.pi,16))
    roots=[];diagnostics=[]
    for alpha in phases:
        root=t.turning_point(alpha);roots.append(root);samples=[]
        for distance in (1e-2,1e-3,1e-4,1e-5):
            W,Wp=t.potential(root['x']+distance,alpha)
            measure=abs(Wp)/(2*W**1.5)
            samples.append(dict(distance=distance,W=W,WKB=measure,
                scaled=measure*distance**1.5,expected=1/(2*math.sqrt(abs(root['Wprime'])))))
        diagnostics.append(dict(alpha=float(alpha),samples=samples,
            positive_control_min=min(t.potential(t.x_from_radius(R),alpha)[0] for R in (1.0001,1.2,2.,2.97)),
            negative_control_max=max(t.potential(t.x_from_radius(R),alpha)[0] for R in (3.03,4.,10.,1000.))))
    reductions=[];clocks=[]
    for alpha,dep,a,kap in itertools.product((0.,.31,math.pi/2),(.05,.15,.3),(.7,1.,2.),(.4,1.)):
        star=math.log((2+dep)/dep)/math.sqrt(2)
        model=f.FRWSupport(a,kap,dep,1,alpha-2*star)
        for R in (1.1,2.,3.,4.,10.):
            x=t.x_from_radius(R);eta=star-x
            ours=t.coefficients(x,alpha)[:4];reference=model.coefficients(eta)
            # coefficients returns m,mp,mpp,k in both routes.
            reductions.append(dict(alpha=alpha,departure=dep,a=a,kappa=kap,R=R,
                coefficients=relative(ours,reference),potential=abs(t.potential(x,alpha)[0]-model.normal_potential(eta))))
        if dep==.15:
            start=star-2.;end=star-1.
            proper=f.proper_evolve(model,start,end)
            reference=t.transport(2.,1.,alpha)
            clocks.append(dict(alpha=alpha,a=a,kappa=kap,map_error=relative(proper['map'],reference),
                eta_error=proper['eta_error'],scale_error=proper['scale_error']))
    transports=[];series=[];past=[];action_rows=[];endpoint=[];tuned=[];controls=[]
    xinitial=t.x_from_radius(1.05)
    for alpha in np.arange(16)*math.pi/16:
        input_maps={L:t.input_basis(L,alpha) for L in (8,12,16,20)}
        output_maps={}
        for order,x0 in itertools.product((8,10,12),(.04,.02,.01)):
            data=t.series_basis(x0,alpha,order)
            output_maps[order,x0]=t.output_basis(alpha,x0,order)
            series.append(dict(alpha=float(alpha),order=order,x0=x0,
                normalized=data['normalized'].tolist(),residual=data['residual'].tolist(),
                resonance=data['resonance'].tolist(),determinant=float(np.linalg.det(data['matrix'])),
                coefficients=data['coefficients'].tolist()))
        out=output_maps[12,.01];incoming=input_maps[20]
        main=np.linalg.solve(out,incoming)
        coarse=np.linalg.solve(t.output_basis(alpha,rtol=1e-10,atol=1e-12),
                               t.input_basis(20,alpha,rtol=1e-10,atol=1e-12))
        different_matches=[np.linalg.solve(t.output_basis(alpha,match=x),t.input_basis(20,alpha,match=x)) for x in (.8,1.2)]
        comparisons=[main,np.linalg.solve(out,input_maps[16]),np.linalg.solve(output_maps[10,.01],incoming),
                     np.linalg.solve(output_maps[12,.02],incoming),coarse,*different_matches]
        transports.append(dict(alpha=float(alpha),comparison_maps=[M.tolist() for M in comparisons],
            input_maps={str(L):M.tolist() for L,M in input_maps.items()},
            output_maps=[dict(order=n,x0=x,M=matrix.tolist()) for (n,x),matrix in output_maps.items()],
            input_convergence=[relative(input_maps[L],input_maps[20]) for L in (8,12,16)],
            absolute_map_errors=[float(np.linalg.norm(M-main)) for M in comparisons],
            relative_map_errors=[relative(M,main) for M in comparisons],
            symplectic_errors=[float(np.linalg.norm(M.T@t.J@M-t.J)) for M in comparisons],
            determinant_errors=[abs(float(np.linalg.det(M))-1) for M in comparisons],
            determinant=float(np.linalg.det(main)),rank=int(np.linalg.matrix_rank(main)),
            C_row_norm=float(np.linalg.norm(main[0]))))
        # Supported-ESU real Floquet representation, without a quantum phase choice.
        period=math.pi/2;monodromy=t.transport(0.,period,alpha,esu=True)
        cos=float(np.trace(monodromy)/2);angle=math.acos(cos)
        B=angle/(period*math.sqrt(1-cos*cos))*(monodromy-cos*np.eye(2))
        periodic_errors=[];norms=[]
        for x in (.2,.6,1.):
            P=t.transport(0.,x,alpha,esu=True)@expm(-B*x)
            Pnext=t.transport(0.,x+period,alpha,esu=True)@expm(-B*(x+period))
            periodic_errors.append(relative(Pnext,P));norms.append(float(np.linalg.norm(P)))
        past.append(dict(alpha=float(alpha),trace=float(np.trace(monodromy)),determinant=float(np.linalg.det(monodromy)),
            log_error=relative(expm(B*period),monodromy),periodic_errors=periodic_errors,periodic_norms=norms))
        # Instantaneous action over all 16 frozen oscillator phases, at common R.
        m,_,_,k,_=t.coefficients(xinitial,alpha);omega=math.sqrt(k/m)
        initial=np.array([[math.sqrt(2/(m*omega))*math.cos(chi),math.sqrt(2*m*omega)*math.sin(chi)]
                          for chi in np.arange(16)*math.pi/8]).T
        root=t.turning_point(alpha);xturn=root['x']
        points=[t.x_from_radius(R) for R in (1.05,1.2,2.,2.97)]+[xturn+1e-3,.1,.02]
        for x in points:
            S=t.transport(xinitial,x,alpha);states=S@initial
            values=[t.actions(x,alpha,z) for z in states.T]
            restored=np.linalg.solve(S,states)
            G=np.array([[2.,.3],[.3,1.]])
            invariant=np.einsum('ij,ij->j',restored,G@restored)/2
            original=np.einsum('ij,ij->j',initial,G@initial)/2
            action_rows.append(dict(alpha=float(alpha),x=x,J_beta=[float(v['J_beta']) for v in values],
                J_beta_prime=[float(v['J_beta_prime']) for v in values],
                J_y=[None if v['J_y'] is None else float(v['J_y']) for v in values],
                pulled_back_invariant_error=relative(invariant,original)))
        # Independent fourth-order derivative check using neighboring exact ODE flows.
        x=1.;state=initial[:,3];h=1e-3
        samples=[]
        for shift in (-2,-1,1,2):
            z=t.transport(x,x+shift*h,alpha)@state
            samples.append(t.actions(x+shift*h,alpha,z)['J_beta'])
        derivative_x=(samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*h)
        derivative_eta=t.actions(x,alpha,state)['J_beta_prime']
        # Generic and exceptional normal-form turning data, via regular beta variables.
        m,mp,_,_,_=t.coefficients(xturn,alpha)
        for yp in (0.,1.):
            initial_turn=np.array([1/math.sqrt(m),math.sqrt(m)*yp-mp/(2*math.sqrt(m))])
            values=[]
            for distance in (1e-2,1e-3,1e-4,1e-5):
                x=xturn+distance;z=t.transport(xturn,x,alpha)@initial_turn
                result=t.actions(x,alpha,z)
                scaled=result['J_y']/math.sqrt(result['W']) if yp==0 else result['J_y']*math.sqrt(result['W'])
                values.append(dict(distance=distance,J_y=float(result['J_y']),scaled=float(scaled)))
            tuned.append(dict(alpha=float(alpha),turn_yprime=yp,samples=values))
        # Both pure future modes: integrate from x=.01 to test endpoint powers.
        data=t.series_basis(.01,alpha,12)
        for x in (.01,.005,.0025,.00125):
            Bx=t.transport(.01,x,alpha,initial=data['matrix'])
            Js=[t.actions(x,alpha,Bx[:,j])['J_beta'] for j in (0,1)]
            endpoint.append(dict(alpha=float(alpha),x=x,beta=Bx[0].tolist(),momentum=Bx[1].tolist(),
                J_beta=[float(v) for v in Js],ratios=[float(Js[0]*x*x/math.sqrt(8)),float(Js[1]/(x*x/(4*math.sqrt(8))))]))
        wrong_out=out.copy();wrong_out[:,1]*=-6
        omitted_in_16=np.linalg.solve(out,t.transport(16,1.,alpha))
        omitted_in_20=np.linalg.solve(out,t.transport(20,1.,alpha))
        wrong_sign=data['matrix'].copy();wrong_sign[1]*=-1
        controls.append(dict(alpha=float(alpha),action_derivative_error=abs(float(derivative_eta+derivative_x))/max(1.,abs(derivative_eta)),
            eta_x_sign_control=abs(float(derivative_eta-derivative_x)),
            wrong_out_determinant=float(np.linalg.det(np.linalg.solve(wrong_out,incoming))),
            wrong_sign_wronskian=float(np.linalg.det(wrong_sign)),
            omitted_input_convergence=relative(omitted_in_16,omitted_in_20)))
    G=np.array([[0.,-1.],[8.,0.]]);state=np.array([.3,.7]);j0=(state[1]**2+8*state[0]**2)/(2*math.sqrt(8))
    const_errors=[]
    for time in (0.,.3,1.,3.):
        z=expm(G*time)@state;const_errors.append(abs((z[1]**2+8*z[0]**2)/(2*math.sqrt(8))-j0))
    return dict(prereg=t.PUBLIC_PREREG,baseline=t.BASELINE,seed=t.SEED,certificate=t.turning_certificate(),
        frobenius=t.frobenius_certificate(),action_identity=action_identity(),roots=roots,diagnostics=diagnostics,
        reductions=reductions,clocks=clocks,transports=transports,series=series,past=past,
        actions=action_rows,endpoint=endpoint,tuned=tuned,controls=controls,constant_oscillator_errors=const_errors)


def finalize(report):
    json.dumps(report,allow_nan=False)
    roots=report['roots'];diag=report['diagnostics'];series=report['series'];transports=report['transports']
    checks=dict(
        normalization=bool(report['clocks']) and all(max(r['map_error'],r['eta_error'],r['scale_error'])<1e-8 for r in report['clocks']),
        phase_reduction=bool(report['reductions']) and all(max(r['coefficients'],r['potential'])<1e-8 for r in report['reductions']),
        turning_location=t.verify_certificate(report['certificate']) and len(roots)==80
              and all(2.97<=r['R']<=3.03 and abs(r['W'])<1e-9 and r['Wprime']<-20 for r in roots)
              and all(r['positive_control_min']>0 and r['negative_control_max']<0 for r in diag),
        adiabatic_diagnostic=len(diag)==80 and all(abs(r['samples'][-1]['scaled']/r['samples'][-1]['expected']-1)<1e-3
              and all(r['samples'][i+1]['WKB']>r['samples'][i]['WKB'] for i in range(3)) for r in diag),
        instantaneous_actions=t.valid_actions(report) and report['action_identity']=='0' and len(report['tuned'])==32
              and all(abs(r['samples'][-1]['scaled']-.5)<1e-3 for r in report['tuned'])
              and all(abs(v-1)<2e-4 for r in report['endpoint'] if r['x']==.00125 for v in r['ratios'])
              and max(max(r['J_beta']) for r in report['actions'])>2,
        frobenius_basis=t.valid_future(report) and report['frobenius']==t.frobenius_certificate() and len(series)==144
              and all(max(r['normalized'])<1e-8 and max(abs(x) for x in r['resonance'])==0 for r in series),
        past_basis=len(report['past'])==16 and all(abs(r['trace'])<2-1e-3 and abs(r['determinant']-1)<1e-8
              and max([r['log_error'],*r['periodic_errors']])<1e-8 for r in report['past']),
        transport_convergence=t.valid_transport(transports),
        symplectic_completion=t.valid_transport(transports) and all(r['rank']==2 and r['C_row_norm']>0 for r in transports),
        negative_controls=bool(report['controls']) and all(r['action_derivative_error']<1e-8
              and abs(r['wrong_out_determinant']+1/6)<1e-8 and abs(r['wrong_sign_wronskian']+1)<1e-8
              and r['omitted_input_convergence']>.01 for r in report['controls'])
              and max(report['constant_oscillator_errors'])<1e-12
              and max(r['pulled_back_invariant_error'] for r in report['actions'])<1e-10,
        scope=t.SCOPE==dict(all_classical_invariants_excluded=False,quantization='NOT_DERIVED',
            preferred_complex_structure='NOT_SPECIFIED',coupled_rotor='NOT_DERIVED',
            general_support_independence='NOT_ESTABLISHED',Phi_selection='NOT_DERIVED',causality_gate='OPEN'),failure_paths=True)
    checks={k:bool(v) for k,v in checks.items()}
    good=dict.fromkeys(t.DEPENDENCIES,True)
    affirmative=t.verdict(good,report['certificate'],transports,report,report)
    failures=all(affirmative[t.TARGETS[target]]!='UNRESOLVED' for target in t.TARGETS)
    for gate,targets in t.DEPENDENCIES.items():
        for missing in (False,True):
            bad=good.copy()
            if missing:bad.pop(gate)
            else:bad[gate]=False
            result=t.verdict(bad,report['certificate'],transports,report,report)
            failures=failures and all(result[t.TARGETS[target]]=='UNRESOLVED' for target in targets)
    checks['failure_paths']=bool(failures)
    report.update(checks=checks,checks_passed=all(checks.values()),
                  verdict=t.verdict(checks,report['certificate'],transports,report,report))
    return report


def render(report):
    if 'checks' not in report:return '# TT turning and transport\n\nUNRESOLVED: '+report['error']+'\n'
    roots=report['roots']
    lines=['# TT turning and complete asymptotic transport','',f'Freeze: `{t.PUBLIC_PREREG}`.','',
        'The global turning bound is certified by exact quadratic-form identities and rational sign bounds.',
        f"Sampled R_turn range: {min(r['R'] for r in roots):.12f} to {max(r['R'] for r in roots):.12f}.",
        'The sampled extrema do not replace the frozen all-phase window [2.97,3.03].','',
        '| Gate | Pass |','|---|---|']
    lines.extend(f'| {k} | {v} |' for k,v in report['checks'].items())
    lines+=['','```json',json.dumps(report['verdict'],indent=2),'```','',
            'The named instantaneous actions fail to stay invariant; chosen exact quadratic invariants still exist.',
            'The complete map retains both frozen and decaying coefficients. No Phi or quantization is derived.','']
    return '\n'.join(lines)


def main(argv=None):
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir','--output',dest='output',type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        report=finalize(run_probe());text=render(report);payload=json.dumps(report,indent=2,allow_nan=False)+'\n'
    except Exception as exc:
        report=dict(checks_passed=False,error=type(exc).__name__+': '+str(exc),verdict=t.verdict({}))
        text=render(report);payload=json.dumps(report,indent=2,allow_nan=False)+'\n'
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'probe.json').write_text(payload);(args.output/'probe.md').write_text(text)
    return 0 if report['checks_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
