"""Post-review source/metric controls for #312; not new constrained solutions.

Reuse archived coefficients without another elliptic solve. Mixed source/metric
controls are off shell and do not define a unique gravitational-force split.
The preregistration, original archive and its failed verdict remain unchanged.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves import two_regions as tr
from experiments.closure_ledger.two_regions_probe import CASES

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = 'experiments/closure_ledger/runs/20260927_two_regions/initial_data.json'


def load_control(record, field_case, metric_case):
    """Keep one archived metric fixed while changing its supporting field."""
    row = record['cases'][metric_case][-1]
    return tr.InitialData(row['degree'], CASES[field_case], np.array(row['coefficients']),
                          row['iterations'], row['projected_residual'])


def active_density(kinetic, spatial, potential):
    """rho+tr(S), with kinetic=G_AB Pi^A Pi^B, no factor 1/2.

    This contraction controls Ricci focusing; it is not a two-body force law.
    Spatial stress need not be isotropic, so 3p means its trace only.
    """
    rho = .5*(kinetic+spatial)+potential
    trace = 1.5*kinetic-.5*spatial-3*potential
    return rho+trace


def background_active_density(eta):
    """Exact breathing background, Einstein normal derivative Pi=R'/sqrt(f)."""
    R = tr.Q0*np.cos(2*eta)
    Rp = -2*tr.Q0*np.sin(2*eta)
    f = 1-R*R/6
    kinetic, potential = Rp*Rp/f**3, 1.5/f**2
    return 2*kinetic-2*potential


def run(record):
    out = dict(kind='post-review diagnostic controls, not preregistered gates',
               original_verdict=record['verdict'], grids=[], source_sha256={})
    for path in [ARCHIVE, 'geometrodynamics/waves/two_regions.py',
                 'experiments/closure_ledger/two_regions_controls.py',
                 'docs/two_gravitating_regions_prereg.md']:
        out['source_sha256'][path] = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    configurations = [('pair', 'pair'), ('pair', 'round'), ('A_only', 'pair'), ('B_only', 'pair'),
                      ('A_only', 'round'), ('B_only', 'round')]
    for radial, angular in [(64, 192), (96, 288)]:
        xy, _ = tr.disk_grid(radial, angular)
        cache, rows = {}, []
        def diagnostic(field, metric):
            key = field, metric
            if key not in cache:
                cache[key] = tr.diagnostics(load_control(record, field, metric), radial, angular)
            return cache[key]
        for field, metric in configurations:
            d = diagnostic(field, metric)
            rows.append(dict(field=field, metric=metric, constraint_solved=field == metric,
                             H_normalized_max=d['H_normalized_max'],
                             rates_01=[r['direct_rate_01'] for r in d['ledger']],
                             energies=[r['energy'] for r in d['regions']]))
        responses = {}
        for mode in ('solved', 'fixed_round'):
            energies = []
            for case in ('pair', 'A_perturbed', 'B_perturbed'):
                metric = case if mode == 'solved' else 'round'
                energies.append(np.array([r['energy'] for r in diagnostic(case, metric)['regions']]))
            responses[mode] = np.stack([(e-energies[0])/.001 for e in energies[1:]], axis=1).tolist()
        responses['solved_minus_fixed_round'] = (np.array(responses['solved'])-np.array(responses['fixed_round'])).tolist()
        pair = load_control(record, 'pair', 'pair').fields(xy)
        spatial = pair['psi']**-4*pair['S']
        active = active_density(0., spatial, pair['U'])
        # Exact cap supremum is independent of sampled grid alignment.
        separation = float(np.arccos(abs(tr.CENTERS[0] @ tr.CENTERS[1])))
        overlap = float(np.cosh(8*np.cos(separation-.45))/np.cosh(8))
        out['grids'].append(dict(grid=[radial, angular], controls=rows, energy_response=responses,
                                active_density_range=[float(active.min()), float(active.max())],
                                active_identity_max_error=float(np.max(abs(active+2*pair['U']))),
                                other_seed_supremum_in_either_cap=overlap))
    out['background_phase_prediction'] = {
        'formula': '2*(dR/deta)^2/f^3 - 3/f^2; R=sqrt(3)/2 cos(2 eta)',
        'eta_0': float(background_active_density(0.)),
        'eta_pi_over_4': float(background_active_density(np.pi/4)),
        'limitation': 'Ricci focusing changes sign; region recoil sign is not determined by this scalar alone'}
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'experiments/closure_ledger/runs/20260927_two_regions/review_controls.json')
    args = parser.parse_args()
    result = run(json.loads((ROOT/ARCHIVE).read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result['grids'][-1], indent=2))


if __name__ == '__main__':
    main()
