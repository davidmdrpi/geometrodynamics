"""Check the public 11e625f driven-normal-stress freeze; fail closed."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from geometrodynamics.bulk import tt_triangle_rotor as rotor
from geometrodynamics.waves import driven_normal_stress as d
from geometrodynamics.waves import reciprocal_scalar_tt as rt


def balance_checks(candidate, times):
    errors = dict.fromkeys(("constant_source", "radial", "angular", "normal", "full"), 0.)
    constant = candidate.k*rotor.axis_tensor(np.array([0., 0., 1.]))
    negative = dict.fromkeys(("zero_source", "reversed_source", "doubled_speed", "unequal_quadratures"), 0.)
    for t in times:
        z, source = candidate.orbit(t), candidate.source(t)
        actual, required = rotor.tensor_components(source, z["n"]), rotor.tensor_components(z["required_source"], z["n"])
        for key in ("radial", "angular", "normal"):
            errors[key] = max(errors[key], float(np.linalg.norm(actual[key]-required[key])/candidate.k))
        errors["full"] = max(errors["full"], float(np.linalg.norm(source-z["required_source"])/candidate.k))
        errors["constant_source"] = max(errors["constant_source"], float(np.linalg.norm(source-constant)/candidate.k))
        for key, residual in (("zero_source", z["required_source"]),
                              ("reversed_source", z["required_source"]+source),
                              ("doubled_speed", candidate.orbit(t, 2.)["required_source"]-source),
                              ("unequal_quadratures", z["required_source"]-candidate.source(t, 1.1))):
            negative[key] = max(negative[key], float(np.linalg.norm(residual)/candidate.k))
    return {"normalized_residuals": errors, "necessity_control_residuals": negative}


def evolution_checks(candidate, times):
    fine = candidate.integrate(times)
    coarse = candidate.integrate(times, rtol=1e-10, atol=1e-12)
    tensors, coarse_tensors = rt.tensor(fine[:, :5]), rt.tensor(coarse[:, :5])
    expected = np.array([candidate.orbit(t)["beta"] for t in times])
    cone = [rotor.nearest_uniaxial(beta) for beta in tensors]
    # Negative amplitude: select the isolated most-negative eigenvalue. Do
    # not track a director through a degenerate eigenspace or the cone apex.
    axes = np.linalg.eigh(tensors)[1][:, :, 0]
    predicted_axes = np.array([candidate.orbit(t)["n"] for t in times])
    projector_error = np.max(np.linalg.norm(
        np.einsum("ti,tj->tij", axes, axes)-np.einsum("ti,tj->tij", predicted_axes, predicted_axes), axis=(1, 2)))
    half_period_overlap = float(abs(axes[0] @ axes[len(times)//2]))
    q_expected = np.array([candidate.scalar(t)[0] for t in times])
    return {"trajectory_error": float(np.max(np.linalg.norm(tensors-expected, axis=(1, 2)))/abs(candidate.A)),
            "refinement_error": float(np.max(np.linalg.norm(tensors-coarse_tensors, axis=(1, 2)))/abs(candidate.A)),
            "max_cone_distance": max(row["distance"] for row in cone)/abs(candidate.A),
            "max_eigenvalue_distance": max(row["eigenvalue_distance"] for row in cone)/abs(candidate.A),
            "director_projector_error": float(projector_error),
            "initial_vs_half_period_axis_overlap": half_period_overlap,
            "scalar_free_evolution_error": float(np.max(np.linalg.norm(fine[:, 10:26]-q_expected, axis=1))/candidate.amplitude),
            "max_tensor_norm": float(np.max(np.linalg.norm(tensors, axis=(1, 2)))),
            "times": times.tolist(), "tensor_coefficients_over_abs_A": (fine[:, :5]/abs(candidate.A)).tolist(),
            "cone_distances_over_abs_A": [row["distance"]/abs(candidate.A) for row in cone]}


def run_probe(progress=lambda message: None):
    progress("exact polynomial, charge and frequency certificates")
    exact, frequency = d.exact_certificate(), d.frequency_certificate()
    candidate = d.DrivenRotor()
    times = np.linspace(0., candidate.period, 401)
    progress("full improved stress and all ten momentum charges on two grids")
    spatial = d.spatial_checks(candidate, 8, 16)+d.spatial_checks(candidate, 12, 24)
    negative_constraints = d.negative_constraint_controls()
    progress("normal, radial, angular and full tensor balance")
    balance = balance_checks(candidate, times)
    progress("unrestricted five-component tensor and free scalar evolution")
    evolution = evolution_checks(candidate, times)
    progress("amplitude and radius controls; fail-closed verdict controls")
    amplitudes = []
    for s in (.02, .01, .005):
        model = d.DrivenRotor(amplitude=s)
        amplitudes.append({"s": s, "field_source_norm": float(np.linalg.norm(model.source(0.))),
                           "tensor_norm": float(np.linalg.norm(model.orbit(0.)["beta"]))})
    scaling_error = max(abs(amplitudes[i][key]/amplitudes[i+1][key]-4)
                        for i in (0, 1) for key in ("field_source_norm", "tensor_norm"))
    radii = []
    for radius in (.7, 2.):
        model = d.DrivenRotor(radius=radius, kappa=.4)
        radii.append({"radius": radius, "kappa": .4, "A": model.A, "k": model.k,
                      "max_tensor_norm": abs(model.A)*math.sqrt(2/3),
                      **balance_checks(model, np.linspace(0., model.period, 401))})
    rot, grad = negative_constraints["missed_rotations"], negative_constraints["missed_gradients"]
    max_abs = lambda values: float(np.max(np.abs(values)))
    charges_max = max(abs(x) for row in spatial for values in row["charges_over_s2"].values() for x in values)
    checks = {
        "exact_polynomial_certificate": exact["all_exact_zero"],
        "action_and_improved_stress": max(row["source_error"] for row in spatial) < 1e-9,
        "complete_constraint_compatibility": charges_max < 1e-9 and all(row["compatible"] for row in spatial),
        "omitted_charge_controls": rot["rejected"] and grad["rejected"]
            and max_abs(rot["three_invariant_charges"]) < 1e-9
            and max_abs(rot["charges"]["rotations"]) > 1e-3
            and max_abs(grad["charges"]["Hamiltonian"]) < 1e-9
            and max_abs(grad["charges"]["rotations"]) < 1e-9
            and max_abs(grad["charges"]["gradient_dipoles"]) > 1e-3,
        "all_tensor_equations": max(balance["normalized_residuals"].values()) < 1e-9,
        "unrestricted_tensor_evolution": max(evolution[key] for key in (
            "trajectory_error", "refinement_error", "max_cone_distance", "max_eigenvalue_distance",
            "director_projector_error", "initial_vs_half_period_axis_overlap", "scalar_free_evolution_error")) < 1e-8,
        "necessity_controls": min(balance["necessity_control_residuals"].values()) > 1e-3,
        "amplitude_radius_and_smallness": scaling_error < 1e-9
            and max(max(row["normalized_residuals"].values()) for row in radii) < 1e-9
            and max([evolution["max_tensor_norm"]]+[row["max_tensor_norm"] for row in radii]) < .05,
        "frequency_restriction": frequency["all_exact"] and frequency["original_speed_excluded"],
        "fail_closed_verdicts": d.verdict_controls(),
    }
    checks = {key: bool(value) for key, value in checks.items()}
    return {"baseline": d.BASELINE, "public_preregistration": d.PUBLIC_PREREG, "seed": d.SEED,
            "preparation": {"s": candidate.amplitude, "a": candidate.radius, "kappa": candidate.kappa,
                            "tensor_A": candidate.A, "Omega": candidate.omega, "source_k": candidate.k,
                            "tensor_period": candidate.period, "scalar_frequency": candidate.omega_scalar,
                            "nonzero_free_tensor_data": "CHOSEN", "support_class": "ASSUMED_PERFECT_FLUID_NO_ANISOTROPIC_STRESS"},
            "exact_certificate": exact, "frequency_certificate": frequency,
            "spatial_checks": spatial, "negative_constraint_controls": negative_constraints,
            "balance": balance, "evolution": evolution, "amplitude_controls": amplitudes,
            "scaling_error": scaling_error, "radius_controls": radii,
            "checks": checks, "checks_passed": not d.failed_checks(checks), "verdict": d.verdict(checks),
            "scope": "O(s^2) homogeneous TT projection with chosen constraint-compatible initial data; not a full Einstein-matter or triangle history",
            "implementation_corrections": []}


def render(report):
    lines = ["# Field-derived driven normal stress", "", report["scope"], "",
             "The nonzero tensor initial data are chosen; existence does not select this preparation.", "",
             f"Public freeze: `{report['public_preregistration']}`", "",
             "| Verdict | Result |", "|---|---|"]
    lines += [f"| {key} | {report['verdict'][key]} |" for key in d.VERDICT_FIELDS]
    if report["verdict"]["failed_checks"]:
        lines += ["", "Failed or missing: "+", ".join(report["verdict"]["failed_checks"])]
    lines += ["", "The frozen analytic candidate has S = 6 s^2 Q_z/(C a^2), A = -3 s^2/(2 C), Omega = sqrt(2)/a.",
              "Its scalar and tensor frequencies have irrational ratio; no periodic joint history is inferred.", "",
              "| Required gate | Pass |", "|---|---|"]
    lines += [f"| {key} | {value} |" for key, value in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/{len(d.REQUIRED_CHECKS)} required gates.", "",
              "| Normalized residual | Maximum |", "|---|---:|"]
    lines += [f"| {key} balance | {value:.6g} |" for key, value in report["balance"]["normalized_residuals"].items()]
    lines += [f"| {key} | {report['evolution'][key]:.6g} |" for key in
              ("trajectory_error", "refinement_error", "max_cone_distance", "director_projector_error")]
    return "\n".join(lines)+"\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    report = run_probe(progress=lambda message: print(message, flush=True))
    # Recompute from the canonical gates at the CLI boundary as well.
    report["checks_passed"] = not d.failed_checks(report["checks"])
    report["verdict"] = d.verdict(report["checks"])
    rendered = render(report)
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir/"probe.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        (args.output_dir/"probe.md").write_text(rendered)
    print(rendered)
    return 0 if report["checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
