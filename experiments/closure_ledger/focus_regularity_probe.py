"""Reproducible, static regularity/transport audit; not a caustic event solve."""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from geometrodynamics.bulk.focus_regularity import RegularFocus, round_caustic_control

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'experiments/closure_ledger/runs/20261011_focus_regularity/result.json'


def measure():
    models = [RegularFocus(epsilon=e) for e in [0., .1]]
    profiles = []
    for m in models:
        S = m.geometry.half_length
        d = m.invariants(np.linspace(-S, S, 1001))
        records = []
        for ell in [0, 1]:
            for wS in [.5, 2., 5.]:
                matrix = m.transfer(wS/S, ell)
                check = m.transfer(wS/S, ell, rtol=2e-13)
                records.append(dict(ell=ell,omega_S=wS,matrix=matrix.tolist(),
                                    determinant_error=abs(float(np.linalg.det(matrix))-1),
                                    refinement_error=float(np.max(abs(matrix-check)))))
        profiles.append(dict(epsilon=m.epsilon,transit_time=m.transit_time(),
                             neck_radius=m.geometry.neck_radius,half_length=S,
                             max_sampled_kretschmann=float(np.max(d['kretschmann'])),
                             max_constraint_density= float(np.max(abs(d['kappa_rho']))),
                             neck_radial_nec=float(m.invariants(0)['kappa_radial_nec']),
                             future_trapped_samples=int(np.sum((d['theta_plus']<0)&(d['theta_minus']<0))),
                             transfers=records))
    differences=[dict(ell=a['ell'],omega_S=a['omega_S'],matrix_max_difference=float(np.max(abs(np.array(a['matrix'])-np.array(b['matrix'])))))
                 for a,b in zip(profiles[0]['transfers'],profiles[1]['transfers'])]
    gates=dict(regular_necks=all(p['neck_radius']>0 and np.isfinite(p['max_sampled_kretschmann']) for p in profiles),
               constraints=all(p['max_constraint_density']<1e-10 for p in profiles),
               transport_accuracy=all(r['refinement_error']<1e-8 and r['determinant_error']<1e-8 for p in profiles for r in p['transfers']),
               different_transport=all(d['matrix_max_difference']>1e-4 for d in differences),
               different_transit=abs(profiles[0]['transit_time']-profiles[1]['transit_time'])>1e-3)
    paths=['geometrodynamics/bulk/focus_regularity.py','geometrodynamics/bulk/finite_mouth.py',
           'experiments/closure_ledger/focus_regularity_probe.py']
    return dict(label='REGULARITY_DOES_NOT_SELECT_EVENT_LAW' if all(gates.values()) else 'CONTROL_FAILED',
                scope='Static prescribed-metric scalar transport control; no nonlinear Einstein event evolution',
                registration='Analytic counterexample and numerical verification, not a held-out preregistered experiment',
                gates=gates,profiles=profiles,transfer_differences=differences,
                round_caustic=round_caustic_control(),
                sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=RUN)
    args=parser.parse_args();result=measure();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['label','gates','transfer_differences']},indent=2))
    if result['label']=='CONTROL_FAILED':raise SystemExit(1)
