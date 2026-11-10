"""Reproduce qualifications for the frozen boundary-pairing protocol.

Run with --output PATH to save JSON. These are static/algebraic controls, not
a dynamic experiment or new quantum-foundations results.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves.boundary_pairing import (
    conditional_spinor_checks, context_gram_counterexample, direct_waist_load,
    interval_dtn, join_dtn, output_gram, polarization_matrix,
    scattering_from_dtn, third_order_response,
)
from geometrodynamics.waves.maxwell_throat_response import (
    Bridge, boundary_dtn, force_matrices, scattering_force_derivative, solve_port,
)

BASELINE = "3293d5bad3c8e103c8514ad85afc0128050004e4"
PROTOCOL = "836bb278074030262819c4bbfd91c5b791068be1"
LIMITS = {
    "dtn_reciprocity": 2e-7, "dtn_whole": 2e-7, "dtn_joining": 2e-7,
    "dtn_scattering": 2e-7, "joining_associativity": 2e-7,
    "load_polarization": 2e-6, "load_phase_derivative": 2e-6,
    "gram_hermiticity": 1e-12, "gram_negative_part": 1e-12,
    "fixed_background_i3": 1e-12, "context_gram_pr_box": 1e-12,
    "conditional_spinor": 1e-12,
}


def relative_error(actual, expected):
    return float(np.linalg.norm(actual - expected) / max(1., np.linalg.norm(expected)))


def evaluate_gates(metrics):
    """Missing, nonfinite and negative error metrics cannot pass a gate."""
    checks = {}
    for name, limit in LIMITS.items():
        value = metrics.get(name)
        valid = (isinstance(value, (int, float, np.integer, np.floating))
                 and not isinstance(value, (bool, np.bool_))
                 and np.isfinite(value) and value >= 0)
        checks[name] = dict(value=float(value) if valid else None, limit=limit,
                            passed=bool(valid and value <= limit))
    passed = all(check['passed'] for check in checks.values())
    return dict(checks=checks,
                verdict="BASELINE_QUALIFIED_DYNAMIC_EXPERIMENT_PENDING" if passed
                else "UNRESOLVED_NUMERICALLY")


def run_qualification():
    metrics = {key: 0. for key in LIMITS}
    rows, loads, errors = [], [], []
    algebra = counterexample = None
    try:
        for depth, omega in [(0., 20.), (0., 60.), (.4, 20.), (.4, 60.), (.8, 20.), (.8, 60.)]:
            bridge = Bridge(lapse_depth=depth)
            x = bridge.limit
            full = interval_dtn(bridge, omega, -x, x)
            inherited = boundary_dtn(bridge, omega)
            scattering = solve_port(bridge, omega, rtol=2e-12)['scattering']
            metrics['dtn_reciprocity'] = max(metrics['dtn_reciprocity'], relative_error(full.T, full))
            metrics['dtn_whole'] = max(metrics['dtn_whole'], relative_error(full, inherited))
            metrics['dtn_scattering'] = max(metrics['dtn_scattering'], relative_error(scattering_from_dtn(full, omega), scattering))
            for cuts in [[0.], [-.37, .21]]:
                points = [-x] + [c*x for c in cuts] + [x]
                pieces = [interval_dtn(bridge, omega, a, b) for a, b in zip(points[:-1], points[1:])]
                joined = join_dtn(pieces)
                error = relative_error(joined, full)
                metrics['dtn_joining'] = max(metrics['dtn_joining'], error)
                if len(pieces) == 3:
                    left = join_dtn([join_dtn(pieces[:2]), pieces[2]])
                    right = join_dtn([pieces[0], join_dtn(pieces[1:])])
                    metrics['joining_associativity'] = max(metrics['joining_associativity'], relative_error(left, right))
                rows.append(dict(lapse_depth=depth, omega=omega, cuts_in_x_over_limit=cuts,
                                 relative_join_error=error))
            inputs = np.array([[1., 0.], [0., 1.], [.3+.1j, .2-.4j]])
            gram = output_gram(scattering, inputs)
            scale = max(1., np.linalg.norm(gram))
            metrics['gram_hermiticity'] = max(metrics['gram_hermiticity'], relative_error(gram.conj().T, gram))
            metrics['gram_negative_part'] = max(metrics['gram_negative_part'], max(0., float(-np.linalg.eigvalsh(gram).min()) / scale))
            metrics['fixed_background_i3'] = max(metrics['fixed_background_i3'], abs(third_order_response(gram)) / scale)

        for depth, omega in [(0., 60.), (.4, 20.)]:
            bridge = Bridge(lapse_depth=depth)
            integrated = force_matrices(bridge, omega, rtol=2e-12)['waist_force']
            polarized = polarization_matrix(lambda a: direct_waist_load(bridge, omega, a))
            phase = scattering_force_derivative(bridge, omega, step=3e-6)
            p_error = relative_error(polarized, integrated)
            s_error = relative_error(phase, integrated)
            metrics['load_polarization'] = max(metrics['load_polarization'], p_error)
            metrics['load_phase_derivative'] = max(metrics['load_phase_derivative'], s_error)
            loads.append(dict(lapse_depth=depth, omega=omega, polarization_error=p_error,
                              phase_derivative_error=s_error))

        counterexample = context_gram_counterexample()
        metrics['context_gram_pr_box'] = max(
            abs(counterexample['chsh'] - 4),
            *(max(row['normalization_error'], row['marginal_error'],
                  max(0., -row['min_eigenvalue'])) for row in counterexample['contexts']))
        algebra = conditional_spinor_checks()
        metrics['conditional_spinor'] = max(algebra['residuals'].values())
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        # An incomplete run must not pass because its maxima started at zero.
        errors.append(f"{type(exc).__name__}: {exc}")
        metrics = {}
    result = evaluate_gates(metrics)
    result.update(schema_version=1, baseline_commit=BASELINE, protocol_commit=PROTOCOL,
                  scope="Static boundary and conditional algebra qualification only",
                  dynamic_mouth_response="NOT_IMPLEMENTED",
                  physical_spinor_matching="NOT_DERIVED",
                  detector_frequency_law="NOT_DERIVED",
                  common_all_setting_history_kernel="NOT_DERIVED",
                  cases=rows, loads=loads, conditional_algebra=algebra,
                  context_gram_counterexample=counterexample, errors=errors)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    result = run_qualification()
    rendered = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')
    return 0 if result['verdict'] == 'BASELINE_QUALIFIED_DYNAMIC_EXPERIMENT_PENDING' else 1


if __name__ == '__main__':
    raise SystemExit(main())
