"""Action-level kinetic admission: fail an inadmissible tangent, do not label it a ghost."""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from geometrodynamics.bulk.egb_shell_action import ScalarShellAction,second_derivative
from geometrodynamics.bulk import negative_egb as ne
from experiments.closure_ledger.egb_shell_dynamics_probe import derivative

ROOT=Path(__file__).resolve().parents[2]
ORIGINAL=ROOT/'experiments/closure_ledger/runs/20261012_egb_shell_dynamics/samples.npz'
RUN=ROOT/'experiments/closure_ledger/runs/20261013_egb_shell_kinetic/result.json'


def measure():
    a=ScalarShellAction();m=a.shell;r=m.r0;q=np.sqrt(a.q2)
    coefficients=a.bare_coefficients();linear=a.linear_obstruction()
    rho,p=a.bulk_stress(a.X0)
    target=np.array([ne.exterior_density(m.alpha,m.radius),ne.exterior_pressure(m.alpha,m.radius)])
    bulk_error=float(np.max(abs(np.array([rho,p])-target)/(1+abs(target))))
    shell_error=float(np.max(abs(np.array(a.surface_stress(a.X0))-np.array(m.equilibrium))))
    eos_errors=[]
    for ratio in [.9,1.,1.1]:
        s,p=a.surface_stress(ratio*a.X0)
        s0,p0=m.equilibrium
        eos_errors.append(abs(p-p0-m.eta*(s-s0))/(1+abs(p)))
    hessian_rows=[]
    family=[]
    for lam in [-1.,0.,1.]:
        model=ScalarShellAction(bulk_lambda=lam)
        c=model.bare_coefficients()
        family.append(dict(bulk_lambda=lam,background_stress=list(model.bulk_stress(model.X0)),
                           coefficients=c,bulk_material_admissible=c['bulk_time']>0))
        for kind in ['bulk','surface']:
            for relative_step in [1e-3,5e-4]:
                h=q*relative_step
                ht=second_derivative(lambda z:model.material_lagrangian(kind,z),q,h)
                hs=-second_derivative(lambda z:model.material_lagrangian(kind,q,z),0.,h)
                he=second_derivative(lambda z:model.local_shift_stationary_energy(kind,z),q,h)
                exact_t=c[kind+'_time'];exact_s=c[kind+'_space']
                hessian_rows.append(dict(bulk_lambda=lam,kind=kind,relative_step=relative_step,
                                         measured_time=float(ht),measured_space=float(hs),
                                         measured_shift_stationary_energy_time=float(he),
                                         normalized_error=float(max(abs(ht-exact_t)/(1+abs(exact_t)),abs(hs-exact_s)/(1+abs(exact_s)),abs(he-exact_t)/(1+abs(exact_t))))))
    slope_fd=derivative(lambda x:m.density(x,0.),r,1e-5*r)
    charge_fd=derivative(a.frozen_total_charge,r,1e-5*r)
    current_fd=a.frozen_current_residual(r,1e-7,0.)/1e-7
    checks=dict(bulk_stress_error=bulk_error,surface_stress_error=shell_error,
                affine_eos_error=float(max(eos_errors)),
                maximum_hessian_error=max(x['normalized_error'] for x in hessian_rows),
                density_slope_error=float(abs(slope_fd-linear['geometric_density_slope'])/(1+abs(slope_fd))),
                charge_slope_error=float(abs(charge_fd-linear['total_charge_slope_no_S3_area'])/(1+abs(charge_fd))),
                current_linear_error=float(abs(current_fd-linear['current_residual_per_velocity'])/(1+abs(current_fd))),
                equilibrium_scalar_current_residual=a.frozen_current_residual(r,0.,0.))
    tangent=[]
    scale=1+abs(m.equilibrium[0])
    for eps in [1e-3,1e-4,1e-5,1e-6]:
        residual=a.frozen_density_residual(r*(1+eps),0.)
        tangent.append(dict(fraction=eps,density_residual=float(residual),
                            dimensionless_residual_per_fraction=float(abs(residual)/(scale*eps)),
                            ratio_to_linear_prediction=float(residual/(r*eps*linear['geometric_density_slope']))))
    with np.load(ORIGINAL,allow_pickle=False) as archive:
        t=archive['stable_t'];y=archive['stable_y'];rr,v=y[:2]
        residual=np.array([a.frozen_density_residual(x,z) for x,z in zip(rr,v)])
        current=np.array([a.frozen_current_residual(x,z,m.acceleration(x,z)) for x,z in zip(rr,v)])
        # Turning points are not selected retrospectively: initial v is exactly zero.
        old=dict(archive_sha256=hashlib.sha256(ORIGINAL.read_bytes()).hexdigest(),
                 sample_count=len(t),max_abs_density_mismatch=float(np.max(abs(residual))),
                 initial_density_mismatch=float(residual[0]),
                 max_abs_scalar_current_residual=float(np.max(abs(current))),
                 min_scalar_current_residual=float(current.min()),max_scalar_current_residual=float(current.max()),
                 interpretation='Residual audit of #330 trajectory; NOT an evolution of this action')
    gates=dict(background_stress=bulk_error<1e-10 and shell_error<1e-10,
               affine_eos=max(eos_errors)<1e-10,
               equilibrium_scalar_junction=checks['equilibrium_scalar_current_residual']==0,
               independent_material_hessians=checks['maximum_hessian_error']<1e-5,
               density_slope=checks['density_slope_error']<1e-7,
               global_charge=checks['charge_slope_error']<1e-7,
               scalar_current=checks['current_linear_error']<1e-7,
               wrong_sign_countercontrol=family[0]['coefficients']['bulk_time']<0,
               canonical_material_coefficients=coefficients['bulk_time']>0 and coefficients['surface_time']>0)
    gates={k:bool(v) for k,v in gates.items()}
    tangent_pass=tangent[-1]['dimensionless_residual_per_fraction']<1e-6
    label=('NUMERICALLY_UNRESOLVED' if not all(gates.values()) else
           'STATIC_ACTION_MATCH_FROZEN_MODE_REJECTED' if not tangent_pass else
           'CONSTRAINED_KINETIC_REDUCTION_REQUIRED')
    sources=['geometrodynamics/bulk/egb_shell_action.py','geometrodynamics/bulk/egb_shell_dynamics.py',
             'geometrodynamics/bulk/negative_egb.py','experiments/closure_ledger/egb_shell_kinetic_probe.py',
             'experiments/closure_ledger/egb_shell_dynamics_probe.py','tests/test_egb_shell_action.py',
             'docs/egb_shell_kinetic_protocol.md']
    return dict(label=label,numerical_gates=gates,checks=checks,
                parameters=dict(R=m.radius,a=m.mouth_angle,alpha=m.alpha,eta=m.eta,bulk_lambda=0.,
                                q=q,X0=a.X0,V=a.potential,power=a.exponent,C=m.bulk_C,
                                surface_power_amplitude=sum(m.equilibrium)/(2*a.exponent),
                                surface_constant_offset=m.equilibrium[1]-sum(m.equilibrium)/(2*a.exponent)),
                bare_coefficients=coefficients,linear_obstruction=linear,
                hessian_checks=hessian_rows,bulk_completion_countercontrols=family,
                tangent_refinement=tangent,original_trajectory_residual=old,
                physical_radial_kinetic_norm=None,
                admission=dict(static_action_matching=gates['background_stress'] and gates['affine_eos'] and gates['equilibrium_scalar_junction'],
                               nominated_bare_material_kinetic_signs_positive=gates['canonical_material_coefficients'],
                               original_frozen_radial_mode_satisfies_constraints=bool(tangent_pass),
                               coupled_metric_scalar_constraint_reduction='UNCOMPUTED',
                               physical_radial_canonical_energy='NOT_ASSIGNED: trial mode violates linearized density junction',
                               smooth_regular_focus='UNCOMPUTED'),
                corrective_kick_budget=0.,
                sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=RUN)
    args=parser.parse_args();result=measure();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['label','numerical_gates','checks','bare_coefficients','linear_obstruction','tangent_refinement','original_trajectory_residual']},indent=2))
    if result['label']=='NUMERICALLY_UNRESOLVED':raise SystemExit(1)
