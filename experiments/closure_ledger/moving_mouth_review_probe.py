#!/usr/bin/env python
"""Post-review single-packet replication and geometric comparison for #340.

No new preregistered success is claimed. The original model/archive is unchanged.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs, harmonics
from moving_mouth_pairs_probe import save_trace


def components(model, eta):
    profiles = [model.ell*(model.ell+1)
                * np.exp(-model.par.sigma**2*model.ell*(model.ell+1)/2)
                * harmonics(model.ell, model.order, point).conj()
                for point in [[1., 0., 0.], [0., 0., 1.]]]
    f = [profiles[0], .65*profiles[1]]
    v = [-eta*(model.generators[2]@f[0]), -eta*(model.generators[0]@f[1])]
    scale = np.sqrt(model.par.energy/model.wave_energy(sum(f), sum(v)))
    states = []
    for fi, vi in zip(f, v):
        s = model.initial(eta)
        s.a = -(-1.)**model.ell*vi*scale/model.omega
        s.p = (-1.)**model.ell*model.omega*fi*scale
        states.append(s)
    original = model.initial(eta)
    for key in ['a', 'p']:
        np.testing.assert_allclose(sum(getattr(s, key) for s in states),
                                   getattr(original, key), atol=3e-17, rtol=3e-15)
    return states


def run(output, original):
    output.mkdir(parents=True, exist_ok=True)
    model = MovingMouthPairs()
    summary = dict(status='POST_REVIEW_REPLICATION_NOT_PREREGISTERED', cases={}, comparisons={})
    for name, eta in [('directional', .7), ('symmetric', 0.)]:
        both = np.genfromtxt(original/f'{name}.csv', names=True, delimiter=',')
        runs = []
        for packet, state in zip(['x_only', 'z_only'], components(model, eta)):
            data, stats = model.evolve(state)
            np.testing.assert_allclose(data['t'], both['t'], atol=5e-12, rtol=0)
            save_trace(output/f'{name}_{packet}.csv', data)
            summary['cases'][f'{name}_{packet}'] = stats
            runs.append(data)
            print(name, packet, json.dumps(stats), flush=True)
        x_run, z_run = runs
        # Compare independently evolved trajectories, not scalar angle changes.
        # This is a geometric diagnostic, NOT a solution with a common field.
        composite_angle = np.arccos(np.clip(np.sum(z_run['x'][:, 0]*x_run['x'][:, 1], axis=1), -1, 1))
        both_angular = np.array([[both[f'angular_{h}_{j}'][-1] for j in range(3)] for h in range(2)])
        residual = both_angular-x_run['angular'][-1]-z_run['angular'][-1]
        comparison = np.column_stack([x_run['t'], np.degrees(both['separation'])-90,
                                      np.degrees(composite_angle)-90,
                                      np.degrees(both['separation']-composite_angle)])
        np.savetxt(output/f'{name}_separation.csv', comparison, delimiter=',',
                   header='t,both_change_deg,composite_change_deg,both_minus_composite_deg',
                   comments='', fmt='%.12e')
        summary['comparisons'][name] = dict(
            both_final_change_deg=float(comparison[-1, 1]),
            composite_final_change_deg=float(comparison[-1, 2]),
            difference_final_deg=float(comparison[-1, 3]),
            final_momentum_nonadditivity=residual.tolist(),
            explanation='Both minus separately evolved single packets includes changed trajectories and backreaction; it is not an isolated bilinear interference torque.')
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    files = sorted(output.glob('*.csv'))+[output/'summary.json']
    (output/'sha256.json').write_text(json.dumps({p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                                for p in files}, indent=2)+'\n')
    print(json.dumps(summary['comparisons'], indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('experiments/closure_ledger/runs/20261017_moving_mouth_review'))
    parser.add_argument('--original', type=Path,
                        default=Path('experiments/closure_ledger/runs/20261016_moving_mouth_pairs'))
    args = parser.parse_args()
    run(args.output, args.original)
