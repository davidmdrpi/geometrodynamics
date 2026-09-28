"""Run the frozen nonlinear initial-data experiment, never a time evolution."""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from geometrodynamics.waves import two_regions as tr

CASES = {'round': (0., 0.), 'A_only': (.02, 0.), 'B_only': (0., .03),
         'pair': (.02, .03), 'A_perturbed': (.021, .03), 'B_perturbed': (.02, .031)}
SCHEDULE = [(12, 20, 64), (20, 28, 88), (28, 36, 112)]


def run():
    xy, weights = tr.disk_grid(64, 192)
    record = dict(freeze=tr.FREEZE, scope='nonlinear initial data and first time jet; no evolution',
                  python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                  cases={}, field_convergence={}, source_sha256={})
    root = Path(__file__).resolve().parents[2]
    for path in ['geometrodynamics/waves/two_regions.py',
                 'experiments/closure_ledger/two_regions_probe.py',
                 'docs/two_gravitating_regions_prereg.md']:
        record['source_sha256'][path] = hashlib.sha256((root/path).read_bytes()).hexdigest()
    solutions, fine_fields = {}, {}
    for name, amplitudes in CASES.items():
        rows, fields = [], []
        for degree, radial, angular in SCHEDULE:
            data = tr.solve(degree, radial, angular, amplitudes)
            diagnostic = tr.diagnostics(data)
            fields.append(data.fields(xy)['psi'])
            rows.append(dict(degree=degree, radial=radial, angular=angular,
                             coefficients=data.coefficients.tolist(), iterations=data.iterations,
                             projected_residual=data.projected_residual, diagnostics=diagnostic))
        solutions[name], fine_fields[name] = data, fields[-1]
        record['cases'][name] = rows
        record['field_convergence'][name] = [float(np.max(abs(b-a))) for a, b in zip(fields, fields[1:])]
        print(name, 'Hn', diagnostic['H_normalized_max'], 'field differences', record['field_convergence'][name], flush=True)
    pair = solutions['pair']
    record['quadrature_refinement'] = [dict(grid=[r, a], diagnostics=tr.diagnostics(pair, r, a))
                                       for r, a in [(32, 96), (64, 192), (96, 288)]]
    response, metric_response = [], []
    for region, center in enumerate(tr.CENTERS):
        base = record['cases']['pair'][-1]['diagnostics']['regions'][region]['energy']
        response.append([(record['cases'][case][-1]['diagnostics']['regions'][region]['energy']-base)/.001
                         for case in ('A_perturbed', 'B_perturbed')])
        points = np.c_[xy, np.sqrt(1-np.sum(xy*xy, axis=1)), np.zeros(len(xy))]
        w, _ = tr.window(points, center)
        metric_response.append([float((weights*w) @ (fine_fields[case]-fine_fields['pair'])/(weights @ w)/.001)
                                for case in ('A_perturbed', 'B_perturbed')])
    record['energy_response_rows_AB_columns_seeds_AB'] = response
    record['mean_psi_response_rows_AB_columns_seeds_AB'] = metric_response
    record['nonadditive_psi_max'] = float(np.max(abs(fine_fields['pair']-fine_fields['A_only']
                                                      -fine_fields['B_only']+fine_fields['round'])))
    curvature = []
    for r, theta in [(.3, .3), (.6, 1.4), (.8, 2.1)]:
        point = np.array([[r*np.cos(theta), r*np.sin(theta)]])
        rho = pair.fields(point)['rho'][0]
        curvature.append(dict(coordinate=[r, theta, .37], twice_rho=float(2*rho),
                              steps=[1e-3, 5e-4, 2.5e-4],
                              residuals=[tr.coordinate_curvature(pair, [r, theta, .37], h)-2*rho
                                         for h in (1e-3, 5e-4, 2.5e-4)]))
    record['independent_coordinate_curvature'] = curvature
    finest = [rows[-1]['diagnostics'] for rows in record['cases'].values()]
    record['gates'] = dict(
        finest_hamiltonian=all(d['H_normalized_max'] < 1e-7 for d in finest),
        finest_field_convergence=all(v[-1] < 1e-7 for v in record['field_convergence'].values()),
        all_successive_field_differences=all(max(v) < 1e-7 for v in record['field_convergence'].values()),
        stress_balance=all(abs(row['balance_error']) < 1e-5 for d in finest for row in d['ledger']),
        independent_metric_responses=all(abs(metric_response[i][1-i]) > 1e-8 for i in range(2)),
        coordinate_curvature_refines=all(abs(v['residuals'][-1]) < abs(v['residuals'][0])/8 for v in curvature))
    # The freeze's wording includes successive differences, so retain its strict
    # all-level gate even when the final pair of grids has converged.
    record['verdict'] = 'PASS_INITIAL_DATA_ONLY' if all(record['gates'].values()) else 'REGISTERED_GATE_FAILURE'
    record['not_established'] = ['finite-time independent motion', 'causal momentum transfer',
                                'gravitational recoil balance', 'self-bound objects or mouths',
                                'emergent quantum mechanics']
    record['future_trapped_spheres'] = 'excluded on this K=0 slice: theta_plus=H, theta_minus=-H'
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('experiments/closure_ledger/runs/20260927_two_regions/initial_data.json'))
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'verdict': result['verdict'], 'gates': result['gates']}))
    return 0 if result['verdict'] == 'PASS_INITIAL_DATA_ONLY' else 1


if __name__ == '__main__':
    raise SystemExit(main())
