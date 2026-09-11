"""Round 12 exact support and preparation audit; no coupled evolution."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from geometrodynamics.waves import odd_multiplet_support as ms
from geometrodynamics.waves import reciprocal_scalar_tt as rt
from geometrodynamics.waves.scalar_esu_support import improved_stress, HomogeneousSupport
from geometrodynamics.waves.backreaction import stress_series

ARCHIVE = Path(__file__).parent/'runs'/'20260910_odd_multiplet_support'
PHASES = (0.,np.pi/4,np.pi/2,3*np.pi/4,np.pi)
EPSILONS = (.01,-.01,.005,-.005,.0025,-.0025)


def failure_controls():
    good=dict.fromkeys(ms.REQUIRED_CHECKS,True)
    evidence={k:dict(certified=True,nullity=0) for k in (1,3,5)}
    rows=[]
    for key in ms.REQUIRED_CHECKS:
        for mode in ('missing','failed'):
            checks=good.copy()
            if mode=='missing': del checks[key]
            else: checks[key]=False
            v=ms.verdict(checks,evidence)
            rows.append(dict(key=key,mode=mode,passed=all(v[t]=='UNRESOLVED'
                for t in ms.CHECK_TARGETS[key]) and key in v['failed_checks']))
    return rows


def background_checks(k):
    model=ms.MultipletSupport(k)
    n=model.size; lam=k*(k+2)
    rng=np.random.default_rng(ms.SEED+k)
    off=rng.normal(size=(37,4)); off/=np.linalg.norm(off,axis=1)[:,None]
    results=[]; parity=0.; addition=0.; direct_error=0.; inherited_error=0.
    individual_min=float('inf'); coherent_min=float('inf')
    rules=[]
    for radial,angular in ((8,16),(12,24)):
        points,weights=rt.sphere_quadrature(radial,angular); rules.append((points,f'{radial},{angular}'))
    rules.append((off,'off_grid'))
    for points,label in rules:
        jets=ms.basis_jets(k,points); Y,dY,hY=jets
        plus=ms.basis_jets(k,-points)[0]+Y
        parity=max(parity,float(np.max(abs(plus))))
        addition=max(addition,float(np.max(abs(np.sum(Y*Y,axis=1)-n))/n),
            float(np.max(abs(np.einsum('pni,pnj->pij',dY,dY)-np.eye(3)*n*lam/3))/(n*lam)),
            float(np.max(abs(np.einsum('pn,pni->pi',Y,dY)))/(n*math.sqrt(lam))),
            float(np.max(abs(np.trace(hY,axis1=2,axis2=3)+lam*Y))/(n*lam)))
        error=0.
        for phase in PHASES:
            stresses,j=ms.component_stress(model,phase,points,jets=jets)
            T=stresses.sum(axis=1); target=np.diag([1.,1/3,1/3,1/3])*model.density
            error=max(error,float(np.max(abs(T-target))/model.density))
            # Inherited unit-radius routine: component index plays the independent batch axis.
            sol=dict(phi=j['phi'],dt=j['derivative'][:,:,0],dtt=j['hessian'][:,:,0,0],
                     grad=j['derivative'][:,:,1:],dtgrad=j['hessian'][:,:,0,1:],
                     hess=j['hessian'][:,:,1:,1:],laplacian=np.trace(j['hessian'][:,:,1:,1:],axis1=2,axis2=3))
            inherited_error=max(inherited_error,float(np.max(abs(stress_series(sol)-stresses))/model.density))
            p=0; m=0
            scalar=improved_stress(j['metric'],j['metric'],j['einstein'],j['phi'][p,m],
                                   j['derivative'][p,m],j['hessian'][p,m])
            direct_error=max(direct_error,float(np.max(abs(scalar-stresses[p,m]))/model.density))
            if phase==0:
                Q=stresses[:,:,1:,1:].copy()
                Q-=np.trace(Q,axis1=2,axis2=3)[:,:,None,None]*np.eye(3)/3
                individual_min=min(individual_min,float(np.min(np.max(np.linalg.norm(Q,axis=(2,3)),axis=0))/model.density))
                coherent,_=ms.component_stress(model,phase,points,np.ones((1,n))*math.sqrt(model.amplitude2),jets)
                coherent_min=min(coherent_min,float(np.max(abs(coherent[:,0]-target))/model.density))
        results.append(dict(rule=label,full_stress_relative=error))
    radius_rows=[]; kinetic_rows=[]
    for a in (.7,1.,2.):
        for kap in (.4,1.):
            m=ms.MultipletSupport(k,a,kap); jets=ms.basis_jets(k,off)
            error=0.; kg=0.; kinetic_error=0.; Fmin=float('inf'); Kmin=float('inf')
            for phase in np.linspace(0,2*np.pi,65):
                T,j=ms.component_stress(m,phase,off,jets=jets); T=T.sum(axis=1)
                target=np.diag([1.,1/3,1/3,1/3])*m.density
                error=max(error,float(np.max(abs(T-target))/m.density))
                kg=max(kg,float(np.max(abs(j['box']-j['phi']/a**2))))
                F,K=m.kinetic(j['phi']); vals=np.linalg.eigvalsh(K)
                predicted=np.column_stack([np.repeat((1/(kap*F))[:,None],n-1,axis=1),1/(kap*F)**2])
                kinetic_error=max(kinetic_error,float(np.max(abs(vals-predicted))))
                Fmin=min(Fmin,float(np.min(kap*F))); Kmin=min(Kmin,float(vals.min()))
            T,j=ms.component_stress(m,0.,off,jets=jets); T=T.sum(axis=1)
            eta=j['metric']; E=j['einstein']
            einstein=float(np.max(abs(E+m.cosmological_constant*eta-kap*T))*a*a)
            wrong_lambda=float(np.max(abs(E-m.cosmological_constant*eta-kap*T))*a*a)
            wrong_amp=float(np.max(abs(E+m.cosmological_constant*eta-kap*1.1**2*T))*a*a)
            same_even=HomogeneousSupport(a,kap)
            radius_rows.append(dict(radius=a,kappa=kap,stress_relative=error,KG_absolute=kg,
                einstein_scaled=einstein,wrong_lambda=wrong_lambda,wrong_amplitude=wrong_amp,
                even_control_density_difference=m.density-same_even.density,
                even_control_Lambda_difference=m.cosmological_constant-same_even.cosmological_constant,
                physical_coefficient_squared=m.volume*m.amplitude2))
            kinetic_rows.append(dict(radius=a,kappa=kap,f_min=Fmin,expected_f_min=1-1/(2*n),
                                    K_min=Kmin,eigenvalue_error=kinetic_error))
    # Independent ambient-polynomial differentiation through the existing
    # harmonic basis, not the invariant generators used by the new module.
    euclidean_error=0.
    old=rt.harmonic_multiplet(k)
    coeff=rng.normal(size=n)*math.sqrt(model.amplitude2)
    transfer=old.B.T@old.moments@ms.numerical_basis(k)[0]
    for a in (.7,1.,2.):
        m=ms.MultipletSupport(k,a)
        inherited_model=rt.ReciprocalModel(k,a)
        for phase in (0.,.41,np.pi/2):
            q=math.sqrt(m.volume)*math.cos(phase)*(transfer@coeff)
            velocity=-math.sqrt(m.volume)*m.omega*math.sin(phase)*(transfer@coeff)
            legacy,_=rt.scalar_jets(inherited_model,q,velocity,-m.omega**2*q,off)
            actual,j=ms.component_stress(m,phase,off,coeff[None,:])
            for i in range(len(off)):
                derivative=np.r_[legacy['dt'][i,0],legacy['grad'][i,0]]
                hessian=np.empty((4,4)); hessian[0,0]=legacy['dtt'][i,0]
                hessian[0,1:]=hessian[1:,0]=legacy['dtgrad'][i,0]
                hessian[1:,1:]=legacy['hess'][i,0]
                reference=improved_stress(j['metric'],j['metric'],j['einstein'],
                                          legacy['phi'][i,0],derivative,hessian)
                euclidean_error=max(euclidean_error,float(np.max(abs(reference-actual[i,0]))/m.density))
    even=ms.basis_jets(2,off)[0]+ms.basis_jets(2,-off)[0]
    return dict(degree=k,quadrature=results,parity_absolute=parity,addition_relative=addition,
        inherited_stress_relative=inherited_error,scalar_stress_relative=direct_error,
        euclidean_stress_relative=euclidean_error,
        individual_anisotropy_min=individual_min,coherent_failure_min=coherent_min,
        even_parity_control=float(np.max(abs(even))),radii=radius_rows,kinetic=kinetic_rows)


def time_norms(model,G,points,weights,jets):
    squares=[]
    for phase in (0.,np.pi/4,np.pi/2,3*np.pi/4):
        T=ms.gram_stress(model,phase,points,G,jets)
        squares.append([ms.source_norm(T,weights,model.density,anisotropic=a)**2 for a in (True,False)])
    return np.sqrt(np.mean(squares,axis=0))


def preparation_checks(k,cert):
    model=ms.MultipletSupport(k); n=model.size
    spectral,operators=ms.singular_analysis(k,cert)
    directions=operators['directions']
    # Degree 4k squared-stress norms require angular order >4k. The freeze's
    # illustrative (8,16) grid suffices for P1 but not k=5 sensitivity norms.
    points,weights=rt.sphere_quadrature(12,24); weights/=weights.sum()
    jets=ms.basis_jets(k,points)
    rng=np.random.default_rng(ms.SEED+100*k)
    off=rng.normal(size=(13,4)); off/=np.linalg.norm(off,axis=1)[:,None]
    offjets=ms.basis_jets(k,off)
    candidates=[]
    for name in ('A','L'):
        for side in ('min','max'):
            candidates.append((name+'_'+side,spectral[name][side+'_direction'],False))
    for i in range(20):
        z=rng.normal(size=len(directions)); z/=np.linalg.norm(z)
        candidates.append(('random_'+str(i),np.einsum('d,dij->ij',z,directions),False))
    # Every exact L-kernel basis witness, not one selected example.
    for i,col in enumerate(cert['operators']['L']['kernel_columns']):
        candidates.append(('kernel_'+str(i),ms.raw_kernel_to_haar(k,col),True))
    rows=[]; reconstruction=0.; linearity=0.; normalization=0.; quad_error=0.
    for name,H,certified_zero in candidates:
        normalization=max(normalization,abs(np.trace(H))/n,abs(np.linalg.norm(H)/math.sqrt(n)-1))
        z=np.einsum('dij,ij->d',directions,H)/n
        squared=np.array([z@operators['matrices'][q]@z for q in ('A','L')])
        # An exact kernel has exact zero norm. Taking sqrt of a roundoff-sized
        # quadratic-form cancellation would artificially produce O(sqrt(eps)).
        predicted=np.zeros(2) if certified_zero else np.sqrt(np.maximum(0,squared))
        observed=time_norms(model,model.amplitude2*H,points,weights,jets)
        quad_error=max(quad_error,float(np.max(abs(predicted-observed))))
        row=dict(direction=name,certified_zero=certified_zero,raw_squared_moment_norms=squared.tolist(),
                 anisotropy_per_gram=float(observed[0]),full_per_gram=float(observed[1]),
                 predicted_anisotropy=float(predicted[0]),predicted_full=float(predicted[1]),steps=[])
        for eps in EPSILONS:
            G=model.amplitude2*(np.eye(n)+eps*H); C=ms.factor_gram(G)
            actual=np.linalg.norm(G-model.amplitude2*np.eye(n))/(model.amplitude2*math.sqrt(n))
            step_error=0.; lin_error=0.
            for phase in (0.,.41,np.pi/2):
                full,_=ms.component_stress(model,phase,off,C,offjets); full=full.sum(axis=1)
                reduced=ms.gram_stress(model,phase,off,G,offjets)
                base=ms.gram_stress(model,phase,off,model.amplitude2*np.eye(n),offjets)
                response=ms.gram_stress(model,phase,off,model.amplitude2*H,offjets)
                step_error=max(step_error,float(np.max(abs(full-reduced))/model.density))
                lin_error=max(lin_error,float(np.max(abs(full-base-eps*response))/model.density))
            reconstruction=max(reconstruction,step_error); linearity=max(linearity,lin_error)
            row['steps'].append(dict(epsilon=eps,gram_mismatch=float(actual),
                full_mismatch=float(abs(eps)*observed[1]),anisotropy=float(abs(eps)*observed[0]),
                psd_min=float(np.linalg.eigvalsh(G).min()/model.amplitude2),
                field_reconstruction_relative=step_error,linearity_relative=lin_error))
        rows.append(row)
    # A separate refined quadrature for the extreme direction, not a fitted normalization.
    p2,w2=rt.sphere_quadrature(16,32); w2/=w2.sum()
    fine=time_norms(model,model.amplitude2*spectral['L']['max_direction'],p2,w2,ms.basis_jets(k,p2))
    coarse=time_norms(model,model.amplitude2*spectral['L']['max_direction'],points,weights,jets)
    refinement=float(np.max(abs(fine-coarse)))
    diag=rng.normal(size=n); diag-=diag.mean(); diag/=np.sqrt(np.mean(diag**2))
    diagonal=[]
    for eps in EPSILONS:
        amps=1+eps*diag; amps*=math.sqrt(n/(amps@amps))
        G=model.amplitude2*np.diag(amps*amps)
        norms=time_norms(model,G-model.amplitude2*np.eye(n),points,weights,jets)
        gram=np.linalg.norm(G/model.amplitude2-np.eye(n))/math.sqrt(n)
        diagonal.append(dict(epsilon=eps,amplitude_mismatch=float(np.linalg.norm(amps-1)/math.sqrt(n)),
            gram_mismatch=float(gram),anisotropy=float(norms[0]),full_mismatch=float(norms[1]),
            anisotropy_per_gram=float(norms[0]/gram),full_per_gram=float(norms[1]/gram)))
    trace=time_norms(model,model.amplitude2*np.eye(n),points,weights,jets)
    summary={q:{key:value for key,value in r.items() if not key.endswith('_direction')}
             for q,r in spectral.items()}
    return dict(degree=k,spectrum=summary,directions=rows,diagonal=diagonal,
        normalization_error=float(normalization),quadrature_norm_error=quad_error,
        refined_norm_error=refinement,reconstruction_relative=reconstruction,linearity_relative=linearity,
        trace_direction=dict(anisotropy=float(trace[0]),full=float(trace[1]),expected_full=math.sqrt(4/3)),
        norm_quadratures=[[12,24],[16,32]],certified_kernel_witnesses=len(cert['operators']['L']['kernel_columns']))


def run_probe(recompute=False):
    exact_background=ms.exact_background_certificate()
    exact_addition=[]; backgrounds=[]; certificates=[]; preparations=[]
    for k in (1,3,5):
        print(f'Round 12 degree {k}: exact identities, support, preparation',flush=True)
        exact_addition.append(ms.exact_addition_certificate(k))
        backgrounds.append(background_checks(k))
        path=ARCHIVE/f'kernel_{k}.json'
        if path.exists() and not recompute:
            cert=json.loads(path.read_text())
            op=ms.polynomial_operators(k)
            for name,M in (('A',op['A']),('L',np.vstack([op['Z'],op['A']]))):
                cert['operators'][name]['certified']=ms.verify_certificate(np.vstack([op['trace'],M]),cert['operators'][name])
            if not all(r['certified'] for r in cert['operators'].values()):
                raise ArithmeticError('archived certificate verification failed')
        else: cert=ms.exact_kernel_certificate(k)
        certificates.append(cert)
        preparations.append(preparation_checks(k,cert))
    failures=failure_controls()
    divergence=[ms.exact_divergence_certificate(k) for k in (1,3,5)]
    checks=dict(
        harmonic_parity=all(b['parity_absolute']<1e-10 and b['even_parity_control']>.1 for b in backgrounds),
        addition_identities=all(e['all_zero'] for e in exact_addition) and all(b['addition_relative']<1e-10 for b in backgrounds),
        full_stress_background=all(max(r['full_stress_relative'] for r in b['quadrature'])<1e-10
            and max(r['stress_relative'] for r in b['radii'])<1e-10
            and b['inherited_stress_relative']<1e-10 and b['scalar_stress_relative']<1e-10
            and b['euclidean_stress_relative']<1e-10 for b in backgrounds),
        einstein_normalization=exact_background['all_zero'] and all(all(r['einstein_scaled']<1e-10
            and r['wrong_lambda']>1 and r['wrong_amplitude']>.1
            and r['even_control_density_difference']==0 and r['even_control_Lambda_difference']==0
            for r in b['radii']) for b in backgrounds),
        kinetic_matrix=all(all(r['f_min']>0 and abs(r['f_min']-r['expected_f_min'])<1e-10
            and r['K_min']>0 and r['eigenvalue_error']<1e-10 for r in b['kinetic']) for b in backgrounds),
        component_bound=exact_background['minimum_components']==4 and exact_addition[0]['all_zero'],
        independent_field_controls=all(b['individual_anisotropy_min']>1e-8 and b['coherent_failure_min']>.01 for b in backgrounds),
        exact_preparation_kernel=all(all(c['certified'] for c in cert['operators'].values()) for cert in certificates)
            and all(c['certified'] for c in divergence),
        normalized_sensitivity=all(p['normalization_error']<1e-9 and p['quadrature_norm_error']<1e-9
            and p['refined_norm_error']<1e-9 for p in preparations),
        gram_reconstruction=all(p['reconstruction_relative']<1e-9 and p['linearity_relative']<1e-9 for p in preparations),
        scope_and_order=ms.SCOPE==dict(coupled_dynamical_stability='NOT_ESTABLISHED',
            field_content_preparation_selection='NOT_DERIVED',TT_frequency_transfer='NOT_ESTABLISHED',
            coupled_support_response='NOT_DERIVED',Phi_selection='NOT_DERIVED',causality_gate='OPEN'),
        failure_paths=all(r['passed'] for r in failures),
    )
    evidence={c['degree']:{key:c['operators']['L'][key] for key in ('certified','nullity')} for c in certificates}
    return dict(prereg=ms.PUBLIC_PREREG,review_amendment=ms.REVIEW_AMENDMENT,seed=ms.SEED,
        exact_background=exact_background,exact_addition=exact_addition,backgrounds=backgrounds,
        post_freeze_divergence=divergence,
        certificates=certificates,preparations=preparations,failure_controls=failures,
        checks={key:bool(value) for key,value in checks.items()},kernel_evidence=evidence,
        provenance='P1/P2 supplied before freeze; P3 computed after published review amendment',
        norm_scope='fixed round ESU, not coupled evolution; Haar and one-period averages define norms only')


def finalize(report):
    checks=report.get('checks',{})
    evidence=report.get('kernel_evidence')
    if evidence is not None: evidence={int(k):v for k,v in evidence.items()}
    report['verdict']=ms.verdict(checks,evidence)
    report['checks_passed']=not report['verdict']['failed_checks']
    return report


def render(report):
    v=report['verdict']
    lines=['# Round 12: odd multiplet support and preparation sensitivity','',
           'Prior candidate P1/P2 independently re-derived. P3 follows the published freeze amendment.','',
           '| Required check | Passed |','|---|---|']
    lines += [f"| {key} | {report.get('checks',{}).get(key,False)} |" for key in ms.REQUIRED_CHECKS]
    lines += ['', '| Verdict | Result |','|---|---|']
    lines += [f'| {key} | {value} |' for key,value in v.items()]
    if report['checks_passed']:
        lines += ['','| k | dim fixed trace | rank A | null A | rank L | null L |','|---|---|---|---|---|---|']
        for c in report['certificates']:
            a,l=c['operators']['A'],c['operators']['L']
            lines += [f"| {c['degree']} | {c['input_dimension']} | {a['fixed_trace_rank']} | {a['nullity']} | {l['fixed_trace_rank']} | {l['nullity']} |"]
    else: lines += ['','Physical claims affected by failed/missing checks are UNRESOLVED.']
    lines += ['','No preparation selection, dynamical stability, tensor-frequency transfer, Phi selection or causality resolution is inferred.','']
    return '\n'.join(lines)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--recompute-certificates',action='store_true')
    args=parser.parse_args(argv)
    try:
        report=finalize(run_probe(recompute=args.recompute_certificates))
    except Exception as exc:
        # Computation failure must overwrite a previous affirmative archive too.
        report=finalize(dict(checks={},kernel_evidence=None,error=f'{type(exc).__name__}: {exc}'))
    if args.output_dir:
        args.output_dir.mkdir(parents=True,exist_ok=True)
        certificates=report.pop('certificates',[])
        for cert in certificates:
            (args.output_dir/f"kernel_{cert['degree']}.json").write_text(json.dumps(cert,separators=(',',':'),allow_nan=False)+'\n')
        # Keep compact rank summaries in the report; full integer witnesses are adjacent.
        report['certificates']=[{**c,'operators':{name:{key:value for key,value in r.items()
            if key not in ('kernel_columns','pivots')} for name,r in c['operators'].items()}} for c in certificates]
        (args.output_dir/'probe.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        (args.output_dir/'probe.md').write_text(render(report))
    print(render(report))
    return 0 if report['checks_passed'] else 1


if __name__=='__main__': raise SystemExit(main())
