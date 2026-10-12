"""Reproduce a known unstable action control; reject, do not tune it."""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from geometrodynamics.bulk.focus_action_admission import EllisGhost, scalar_boundary_obstruction

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'experiments/closure_ledger/runs/20261012_focus_action_admission/result.json'


def measure():
    m = EllisGhost()
    mesh = [m.radial_mode(spacing=h) for h in [.04, .02, .01]]
    domains = [m.radial_mode(half_domain=L) for L in [8., 12., 16.]]
    a,b,c = [p['lambda0'] for p in mesh]
    richardson = (4*c-b)/3
    rayleigh = m.gaussian_rayleigh()
    original = scalar_boundary_obstruction()
    numerical = dict(mesh_refines=abs(c-b)<.3*abs(b-a),
                     richardson_agreement=abs(richardson-(4*b-a)/3)<1e-5,
                     domain_convergence=abs(domains[0]['lambda0']-domains[-1]['lambda0'])<1e-5,
                     eigen_residuals=all(v['residual_l2']<1e-8 for v in mesh+domains),
                     negative_variational_bound=rayleigh['quotient']<-.01,
                     reproduced_reference=abs(richardson-(-2.98495893027790))<1e-5)
    admission = dict(original_mouth_contract=abs(original['forced_right_mouth_lapse_derivative'])<1e-10,
                     static_two_ended_stability=richardson>=-.01,
                     bounded_local_kinetic_hamiltonian=False,
                     closed_bulk_scalar_matching_supplied=False,
                     focal_activation_law_supplied=False)
    sources = ['geometrodynamics/bulk/focus_action_admission.py',
               'geometrodynamics/bulk/focus_regularity.py','geometrodynamics/bulk/finite_mouth.py',
               'experiments/closure_ledger/focus_action_admission_probe.py',
               'tests/test_focus_action_admission.py']
    return dict(label='CANDIDATE_REJECTED_STATIC_AND_BOUNDARY_GATES' if all(numerical.values()) else 'NUMERICALLY_UNRESOLVED',
                scope='Explicit massless ghost action, analytic original-neck obstruction, two-ended Ellis radial stability; no nonlinear event',
                known_result='Torii-Shinkai 2013 eqs 3.17-3.18, table I; replication, not discovery',
                numerical_gates=numerical, admission_gates=admission,
                original_neck=original, mesh_refinement=mesh, domain_refinement=domains,
                richardson_lambda=richardson, rayleigh=rayleigh,
                growth_time_over_b=1/np.sqrt(-richardson),
                ghost_kinetic_hessian=-1., corrective_kick_budget=0.,
                sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=RUN)
    args=p.parse_args();r=measure();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:r[k] for k in ['label','numerical_gates','admission_gates','richardson_lambda','growth_time_over_b','rayleigh']},indent=2))
    if r['label']=='NUMERICALLY_UNRESOLVED':raise SystemExit(1)
