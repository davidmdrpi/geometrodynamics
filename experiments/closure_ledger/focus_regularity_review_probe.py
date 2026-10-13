"""Supplement to the frozen original regularity control: no event solve."""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from geometrodynamics.bulk.focus_regularity_review import SmoothRegularFocus, vacuum_neck_control

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'experiments/closure_ledger/runs/20261011_focus_regularity_review/result.json'


def measure():
    families = []
    for collar in [.1, .2]:
        profiles = []
        for epsilon in [0., .1]:
            m = SmoothRegularFocus(epsilon=epsilon, collar_fraction=collar)
            S = m.geometry.half_length
            d = m.invariants(np.linspace(-S, S, 2001))
            records = []
            for ell in [0, 1]:
                for wS in [.5, 2., 5.]:
                    M = m.transfer(wS/S, ell, rtol=2e-13)
                    fine = m.transfer(wS/S, ell, rtol=5e-14)
                    records.append(dict(ell=ell, omega_S=wS, matrix=M.tolist(),
                                        determinant_error=abs(float(np.linalg.det(M))-1),
                                        refinement_error=float(np.max(abs(M-fine)))))
            profiles.append(dict(epsilon=epsilon, transit_time=m.transit_time(),
                                 max_sampled_kretschmann=float(np.max(d['kretschmann'])),
                                 min_sampled_radius=float(np.min(d['radius'])),
                                 min_sampled_kappa_rho=float(np.min(d['kappa_rho'])),
                                 max_sampled_kappa_rho=float(np.max(d['kappa_rho'])),
                                 waist_radial_nec=float(m.invariants(0)['kappa_radial_nec']),
                                 future_trapped_samples=int(np.sum((d['theta_plus']<0)&(d['theta_minus']<0))),
                                 transfers=records))
        differences = [float(np.max(abs(np.array(a['matrix'])-np.array(b['matrix']))))
                       for a, b in zip(profiles[0]['transfers'], profiles[1]['transfers'])]
        families.append(dict(collar_fraction=collar, profiles=profiles,
                             transit_reduction_percent=100*(1-profiles[1]['transit_time']/profiles[0]['transit_time']),
                             transfer_max_entry_differences=differences))
    vacuum = [dict(mouth_angle=a, **vacuum_neck_control(mouth_angle=a)) for a in [.5, 1.3]]
    gates = dict(vacuum_control_shielded=vacuum[0]['exterior_causally_excluded'],
                 causal_countercontrol_not_shielded=not vacuum[1]['exterior_causally_excluded'],
                 finite_positive_necks=all(p['min_sampled_radius']>0 and np.isfinite(p['max_sampled_kretschmann']) for f in families for p in f['profiles']),
                 different_transport=all(min(f['transfer_max_entry_differences'])>1e-4 for f in families),
                 different_transit=all(f['transit_reduction_percent']>1 for f in families),
                 integration_accuracy=all(r['determinant_error']<1e-8 and r['refinement_error']<1e-8 for f in families for p in f['profiles'] for r in p['transfers']))
    paths = ['geometrodynamics/bulk/focus_regularity_review.py',
             'geometrodynamics/bulk/focus_regularity.py', 'geometrodynamics/bulk/finite_mouth.py',
             'experiments/closure_ledger/focus_regularity_review_probe.py',
             'tests/test_focus_regularity_review.py']
    return dict(label='VACUUM_COLLAPSES_SMOOTH_REGULARITY_UNDERDETERMINED' if all(gates.values()) else 'CONTROL_FAILED',
                scope='Analytic vacuum baseline plus inverse-engineered smooth static controls; no matter law or nonlinear event evolution',
                registration='Review supplement with disclosed collar/bump inputs; not a preregistered held-out test',
                radius=1., mouth_angle=.5, production_rtol=2e-13, refinement_rtol=5e-14,
                pilot_precision_failure=dict(rtol=1e-11, refinement_rtol=2e-13,
                                             collar_fraction=.2, epsilon=.1, ell=1, omega_S=.5,
                                             max_entry_error=1.008708636618394e-7, gate=1e-8,
                                             response='Tightened integration; absolute error gate unchanged'),
                gates=gates, vacuum_controls=vacuum, smooth_families=families,
                sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=RUN)
    args = parser.parse_args()
    result = measure()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(label=result['label'], gates=result['gates'], vacuum=result['vacuum_controls'],
                          transit_reductions=[f['transit_reduction_percent'] for f in result['smooth_families']]), indent=2))
    if result['label'] == 'CONTROL_FAILED':
        raise SystemExit(1)
