"""Pointwise scalar ESU support audit; published freeze de55f3f."""

import argparse
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves import scalar_esu_support as ss
from geometrodynamics.waves import reciprocal_scalar_tt as rt
from geometrodynamics.waves.backreaction import stress_series


def spatial_checks():
    rng = np.random.default_rng(ss.SEED)
    rows = []
    for degree in (1, 3, 5):
        model = rt.ReciprocalModel(degree=degree)
        q = rng.normal(size=model.multiplet.dimension)
        q /= np.linalg.norm(q)
        p = rng.normal(size=q.shape)
        for radial, angular in ((4, 8), (6, 12)):
            points, _ = rt.sphere_quadrature(radial, angular)
            jets, frames = rt.scalar_jets(model, q, p, -model.omega_scalar2*q, points)
            T = stress_series(jets)[:, 0]
            phi, grad, hess = jets["phi"][:, 0], jets["grad"][:, 0], jets["hess"][:, 0]
            predicted = ss.stf((2*np.einsum("pi,pj->pij", grad, grad)-phi[:,None,None]*hess)/3)
            error = np.max(np.linalg.norm(predicted-ss.stf(T[:,1:,1:]), axis=(1,2)))
            scale = max(1., float(np.max(np.linalg.norm(T, axis=(1,2)))))
            negative, _ = rt.scalar_jets(model, q, p, -model.omega_scalar2*q, -points)
            rows.append(dict(degree=degree, points=len(points), absolute=float(error),
                             scale=scale, relative=float(error/scale),
                             odd_residual=float(np.max(np.abs(phi+negative["phi"][:,0]))),
                             max_anisotropic_stress=float(np.max(np.linalg.norm(predicted,axis=(1,2))))))
    # Spatially isotropic reciprocal controls do not imply full Einstein solutions.
    points = rng.normal(size=(200,4))
    points /= np.linalg.norm(points,axis=1)[:,None]
    frames = np.einsum("iab,pb->pia", rt.QUATERNION_DERIVATIVES, points)
    B = np.array([.2,-.3,.4,.1])
    z = points@B
    h = 2+z
    dh = np.einsum("pia,a->pi", frames,B)
    phi = 1/h
    dphi = -dh/h[:,None]**2
    ddphi = (2*np.einsum("pi,pj->pij",dh,dh)/h[:,None,None]**3
              +z[:,None,None]*np.eye(3)/h[:,None,None]**2)
    Ttf = ss.stf((2*np.einsum("pi,pj->pij",dphi,dphi)-phi[:,None,None]*ddphi)/3)
    reciprocal = dict(max_absolute=float(np.max(np.linalg.norm(Ttf,axis=(1,2)))),
                      h_lower_bound=float(2-np.linalg.norm(B)),
                      phi_lower_bound=float(1/(2+np.linalg.norm(B))),
                      nonconstant_range=float(np.ptp(phi)),
                      full_Einstein_solution="NOT_ASSERTED",
                      pole_control=dict(A=.2, norm_B=float(np.linalg.norm(B)),
                                        smooth=False, reason="abs(A)<=norm(B): denominator vanishes on S3"))
    return dict(inherited=rows, reciprocal=reciprocal)


def background_checks():
    zero=(0.,0.,0.)
    rows=[]
    for a in (.7,1.,2.):
        for k in (.4,1.):
            m=ss.HomogeneousSupport(a,k)
            worst=wave=trace=0.
            for theta in np.linspace(0,2*np.pi,13):
                d=ss.direct_geometry(m,theta*a,2,zero,zero,zero)
                residual=d["einstein"]+m.cosmological_constant*np.diag([-1,1,1,1])-k*d["stress"]
                worst=max(worst,float(np.linalg.norm(residual)))
                wave=max(wave,float(abs(d["wave"])))
                trace=max(trace,float(abs(d["trace"])))
            F,K=m.kinetic(np.linspace(0,2*np.pi*a,49))
            # Genuine changes of the same full geometry/stress calculation.
            d=ss.direct_geometry(m,.31*a,0,zero,zero,zero)
            wrong_lambda=d["einstein"]-m.cosmological_constant*np.diag([-1,1,1,1])-k*d["stress"]
            P,Pd,Pdd=m.jets(.31*a)
            scaled=ss.direct_geometry(m,.31*a,0,(.1*P,.1*Pd,.1*Pdd),zero,zero,epsilon=1.)
            wrong_amplitude=scaled["einstein"]+m.cosmological_constant*np.diag([-1,1,1,1])-k*scaled["stress"]
            empty=ss.direct_geometry(m,.31*a,0,(-P,-Pd,-Pdd),zero,zero,epsilon=1.)
            # Lambda cancels from the sum of the 00 and 11 equations.
            zero_enthalpy=(empty["einstein"][0,0]+empty["einstein"][1,1]
                           -k*(empty["stress"][0,0]+empty["stress"][1,1]))
            rows.append(dict(radius=a,kappa=k,Einstein_absolute=worst,Einstein_scale=1/a**2,
                             Einstein_relative=worst*a**2,wave_absolute=wave,trace_absolute=trace,
                             F_min=float(np.min(F)),K_min=float(np.min(K)),K_max=float(np.max(K)),
                             reversed_Lambda_residual=float(np.linalg.norm(wrong_lambda)*a*a),
                             amplitude_1p1_residual=float(np.linalg.norm(wrong_amplitude)*a*a),
                             zero_history_enthalpy_residual=float(zero_enthalpy*a*a)))
    return rows


def variation_checks():
    rng=np.random.default_rng(ss.SEED+1)
    rows=[]
    for a,k in ((1.,1.),(.7,.4)):
        m=ss.HomogeneousSupport(a,k)
        for degree in (2,3):
            for theta in (0.,.37,np.pi/2,np.pi):
                field=rng.normal(0,.08,3)/np.array([1.,a,a*a])/np.sqrt(k)
                alpha=rng.normal(0,.06,3)/np.array([1.,a,a*a])
                psi=rng.normal(0,.06,3)/np.array([1.,a,a*a])
                r=ss.stress_response(m,theta*a,degree,field,alpha,psi)
                target=ss.response_tensor(m,degree,r)
                lam=degree*(degree+2)/a**2
                geometric=ss.response_tensor(m,degree,dict(
                    rho=2*(3/a**2-lam)*psi[0], J=-2*psi[1],
                    p=2*psi[2]-2*psi[0]/a**2+2*lam*(psi[0]-alpha[0])/3,
                    Pi=psi[0]-alpha[0]))
                scale=max(1.,float(np.linalg.norm(target)))
                errors=[]; einstein_errors=[]; trace_error=0.; wave_errors=[]
                for epsilon in (1e-3,5e-4):
                    p=ss.direct_geometry(m,theta*a,degree,field,alpha,psi,epsilon)
                    n=ss.direct_geometry(m,theta*a,degree,field,alpha,psi,-epsilon)
                    errors.append(float(np.linalg.norm((p["stress"]-n["stress"])/(2*epsilon)-target)))
                    einstein_errors.append(float(np.linalg.norm((p["einstein"]-n["einstein"])/(2*epsilon)-geometric)))
                    Y=ss.zonal(degree,.83)[0]
                    wave_errors.append(float(abs((p["wave"]-n["wave"])/(2*epsilon)-r["KG"]*Y)))
                    trace_error=max(trace_error,abs(p["trace"]-p["phi"]*p["wave"]),abs(n["trace"]-n["phi"]*n["wave"]))
                rows.append(dict(radius=a,kappa=k,degree=degree,phase=theta,
                                 field_jets=field.tolist(),lapse_jets=alpha.tolist(),psi_jets=psi.tolist(),
                                 epsilon=[1e-3,5e-4],stress_absolute=errors,stress_scale=scale,
                                 stress_relative=[e/scale for e in errors],
                                 stress_error_ratio=errors[0]/max(errors[1],1e-30),
                                 Einstein_absolute=einstein_errors,
                                 Einstein_scale=max(1.,float(np.linalg.norm(geometric))),
                                 wave_absolute=wave_errors,wave_scale=max(1.,abs(r["KG"]*Y)),
                                 off_shell_trace_absolute=float(trace_error),
                                 response_trace_absolute=float(abs(-r["rho"]+3*r["p"]-m.jets(theta*a)[0]*r["KG"]))))
    return rows


def evolution_check(model,degree):
    times=np.linspace(0,2*np.pi*model.radius,401)
    initial=ss.initial_response(model,degree)
    coarse=ss.integrate_response(model,degree,times,initial)
    fine=ss.integrate_response(model,degree,times,initial,rtol=1e-12,atol=1e-14)
    names=("hamiltonian","momentum","spatial","anisotropic","KG","trace")
    scale=max(1.,degree*(degree+2)/model.radius**2*np.linalg.norm(initial))
    reports=[]
    for values in (coarse,fine):
        residuals=[]; pis=[]; rhos=[]; slips=[]
        for t,y in zip(times,values):
            dy=ss.response_rhs(model,degree,t,y)
            r=ss.response_residuals(model,t,degree,y,(dy[1],dy[3]))
            residuals.append([r[name] for name in names])
            pis.append(r["Pi"]);rhos.append(r["rho"]);slips.append(y[2]-r["alpha"])
        reports.append(dict(absolute=dict(zip(names,np.max(np.abs(residuals),axis=0).tolist())),
                            scale=scale,relative=dict(zip(names,(np.max(np.abs(residuals),axis=0)/scale).tolist()))))
    constraint_scale=scale
    initial_res=ss.response_residuals(model,0.,degree,initial)
    actual_amplitude=1e-4
    return dict(radius=model.radius,kappa=model.kappa,degree=degree,
                initial=initial.tolist(),normalization="unit response; physical perturbation epsilon=1e-4",
                initial_constraints={n:float(initial_res[n]) for n in ("hamiltonian","momentum")},
                initial_constraint_scale=constraint_scale,runs=reports,
                refinement_absolute=float(np.max(np.abs(fine-coarse))),
                refinement_scale=max(1.,float(np.max(np.abs(fine)))),
                max_abs_Pi=float(np.max(np.abs(pis))),
                omitted_Pi_Einstein_potential_residual=float(np.max(np.abs(slips))),
                max_density_cross_term=float(np.max(np.abs(rhos))),
                physical_max_field_correction=float(actual_amplitude*np.max(np.abs(fine[:,0]))),
                times=times.tolist(),trajectory=fine.tolist())


def fail_closed_controls():
    rows=[]
    good={name:True for name in ss.REQUIRED_CHECKS}
    for name in ss.REQUIRED_CHECKS:
        for mode in ("missing","failed"):
            checks=good.copy()
            if mode=="missing":del checks[name]
            else:checks[name]=False
            v=ss.verdict(checks)
            rows.append(dict(check=name,mode=mode,passed=all(v[n]=="UNRESOLVED" for n in ss.VERDICT_FIELDS)
                             and name in v["failed_checks"]))
    return rows


def run_probe(progress=lambda message:None):
    progress("exact identities and pointwise spatial controls")
    certificate=ss.exact_certificate()
    spatial=spatial_checks()
    progress("homogeneous background and direct metric variation")
    background=background_checks()
    variation=variation_checks()
    progress("independent scalar/spatial evolution with unused constraints")
    evolutions=[evolution_check(ss.HomogeneousSupport(a,k),l)
                for a,k,l in ((1.,1.,2),(1.,1.,3),(.7,.4,2))]
    controls=fail_closed_controls()
    scope=dict(single_real_scalar=True,odd_condition="IMPOSED_NOT_DERIVED",
               pointwise_not_averaged=True,homogeneous_control_parity="EVEN_EXCLUDED_FROM_PRIMARY_CLASS",
               homogeneous_control_phase="CHOSEN",Lambda="3/(2a^2)_REQUIRED_FOR_CONTROL",
               signal_stress_leading_order="LINEAR_CROSS_TERM_ABOUT_NONZERO_BACKGROUND",
               nonlinear_stability="NOT_ESTABLISHED",throat_support="NOT_TESTED",
               triangle_map="NOT_DERIVED",preparation_law="NOT_DERIVED",Phi_selection="NOT_DERIVED",
               physical_readout="NOT_DERIVED")
    # The symbolic local certificates supplement the written global proof;
    # this flag is not a machine proof of the topology argument.
    proof=Path(__file__).resolve().parents[2]/"docs/scalar_esu_support.md"
    checks=dict(
        stress_identity=certificate["all_zero"] and all(r["relative"]<1e-10 for r in spatial["inherited"]),
        global_obstruction=certificate["all_zero"] and proof.is_file()
            and spatial["reciprocal"]["max_absolute"]<1e-10 and spatial["reciprocal"]["h_lower_bound"]>0,
        homogeneous_background=certificate["all_zero"] and all(r["Einstein_relative"]<1e-10
            and r["wave_absolute"]<1e-10 and r["reversed_Lambda_residual"]>.1
            and r["amplitude_1p1_residual"]>.1 for r in background),
        admissibility=all(r["F_min"]>=1/(2*r["kappa"])-1e-12 and r["K_min"]>0 for r in background)
            and all(abs(r["zero_history_enthalpy_residual"]-2)<1e-10 for r in background)
            and not spatial["reciprocal"]["pole_control"]["smooth"]
            and all(r["odd_residual"]<1e-12 for r in spatial["inherited"]),
        response_variation=all(max(r["stress_relative"])<1e-6
            and max(r["Einstein_absolute"])/r["Einstein_scale"]<1e-6
            and max(r["wave_absolute"])/r["wave_scale"]<1e-6
            and r["off_shell_trace_absolute"]<1e-10 and r["response_trace_absolute"]<1e-10 for r in variation),
        constraint_propagation=all(max(abs(v) for v in r["initial_constraints"].values())/r["initial_constraint_scale"]<1e-8
            and all(max(run["relative"].values())<1e-8 for run in r["runs"])
            and r["refinement_absolute"]/r["refinement_scale"]<1e-7 for r in evolutions),
        nonfluid_control=any(r["max_abs_Pi"]>1e-3 and r["omitted_Pi_Einstein_potential_residual"]>1e-3 for r in evolutions),
        scope_and_order=any(r["max_density_cross_term"]>1e-3 for r in evolutions)
            and scope["homogeneous_control_parity"]=="EVEN_EXCLUDED_FROM_PRIMARY_CLASS",
        fail_closed=all(r["passed"] for r in controls)
            and bool(ss.failed_checks({"unrelated":True})),
    )
    checks={name:bool(value) for name,value in checks.items()}
    return dict(public_preregistration=ss.PUBLIC_PREREG,baseline=ss.BASELINE,seed=ss.SEED,
                exact=certificate,spatial=spatial,background=background,variation=variation,
                evolution=evolutions,scope=scope,fail_closed_controls=controls,checks=checks,
                checks_passed=not ss.failed_checks(checks),verdict=ss.verdict(checks),
                analytic_prediction_corrections=[],
                global_proof_status="WRITTEN_ANALYTIC_PROOF_WITH_SYMBOLIC_LOCAL_CHECKS_NOT_FORMAL_VERIFICATION")


def render(report):
    lines=["# Can one real conformal scalar supply the ESU support?","",
           "The odd-sector result is a pointwise global obstruction in the stated class.",
           "The positive homogeneous control is even and has a generically anisotropic response.","",
           f"Public freeze: `{ss.PUBLIC_PREREG}`.","",
           "| Verdict | Result |","|---|---|"]
    lines += [f"| {name} | {report['verdict'][name]} |" for name in ss.VERDICT_FIELDS]
    if report["verdict"].get("failed_checks"):
        lines += ["","Failed or missing gates: "+", ".join(report["verdict"]["failed_checks"])]
    lines += ["","| Required check | Pass |","|---|---|"]
    lines += [f"| {name} | {report['checks'].get(name, False)} |" for name in ss.REQUIRED_CHECKS]
    lines += ["","The numerical controls do not establish the global exclusion by search.",
              "See docs/scalar_esu_support.md for the component-boundary proof.","",
              "| a | kappa | degree | Fine constraint max (absolute) | Refinement (absolute) | max abs(Pi) |",
              "|---|---|---|---|---|---|"]
    for r in report.get("evolution",[]):
        fine=r["runs"][1]["absolute"]
        lines.append(f"| {r['radius']} | {r['kappa']} | {r['degree']} | {max(fine['hamiltonian'],fine['momentum']):.3e} | {r['refinement_absolute']:.3e} | {r['max_abs_Pi']:.6g} |")
    lines += ["","The even control is not a BAM support selection. The triangle map, history",
              "preparation/weight and physical readout remain unconstructed.",""]
    return "\n".join(lines)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path)
    args=parser.parse_args(argv)
    report=run_probe(progress=lambda message:print(message,flush=True))
    # Recompute at the publication boundary; stale cached success cannot pass.
    report["checks_passed"]=not ss.failed_checks(report.get("checks",{}))
    report["verdict"]=ss.verdict(report.get("checks",{}))
    summary=render(report)
    if args.output_dir:
        args.output_dir.mkdir(parents=True,exist_ok=True)
        (args.output_dir/"probe.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
        (args.output_dir/"probe.md").write_text(summary)
    print(summary)
    return 0 if report["checks_passed"] else 1


if __name__=="__main__":
    raise SystemExit(main())
