"""Round 10: what data does joint closure retain for two independent pairs?

Pre-registered at ``b78157a`` (``docs/joint_closure_composition_prereg.md``,
amendment A1), committed before the implementation. Reports three separate
verdict fields and exits nonzero when any required check fails.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from geometrodynamics.bulk import joint_closure as jc


ANGLE_PAIRS = ((0.35, 0.70), (1.0, 1.3), (math.pi / 2, math.pi / 2), (2.4, 2.7))
BASE_A = (0.0, 0.0, 1.0)


def _settings(gamma, rotation=None):
    a = np.array(BASE_A, dtype=float)
    b = np.array([math.sin(gamma), 0.0, math.cos(gamma)], dtype=float)
    if rotation is not None:
        a, b = rotation @ a, rotation @ b
    return tuple(a), tuple(b)


def run_probe(progress=lambda s: None):
    geometry, windows, punctures, grids = [], [], [], []
    for gamma1, gamma2 in ANGLE_PAIRS:
        for label, R in (("identity", None), ("R_x(0.61)R_z(0.37)", jc.ROTATION_PAIR2)):
            progress(f"angles {gamma1:.3f},{gamma2:.3f} rotation {label}")
            a1, b1 = _settings(gamma1)
            a2, b2 = _settings(gamma2, R)
            tri1, tri2 = jc.Triangle(a1, b1, 1, 1), jc.Triangle(a2, b2, 1, -1)
            masses = jc.reference_sector_masses(a1, b1, a2, b2)
            geometry.append({
                "angles": [gamma1, gamma2], "rotation": label,
                "t": [tri1.t, tri2.t],
                "sector_probabilities": masses["probabilities"],
                "product_marginal_error": masses["product_marginal_error"],
                "associativity_regression": jc.three_factor_associativity(a1, b1, a2, b2),
                "gram": jc.gram_and_metric_checks(tri1, tri2),
                "analytic_route": jc.analytic_gram_route(tri1, tri2),
                "covariance_exchange": jc.covariance_and_exchange(a1, b1, a2, b2),
                "sector_grids": jc.sector_probability_grids(a1, b1, a2, b2),
            })
            windows.append({"angles": [gamma1, gamma2], "rotation": label,
                            **jc.window_convergence(a1, b1, a2, b2)})
            punctures.append({
                "angles": [gamma1, gamma2], "rotation": label,
                "factor1": jc.puncture_geometry(tri1),
                "factor2": jc.puncture_geometry(tri2),
                "excision1": jc.excision_masses(tri1),
                "excision2": jc.excision_masses(tri2),
                "joint_excluded": jc.joint_excision_bound(tri1, tri2)})
            grids.append(jc.grid_refinement(tri1.t))

    progress("level-set controls")
    controls = jc.level_set_controls()
    progress("repository audit")
    audit = jc.repository_joint_rule_audit()

    checks = {
        "Q1 finite-difference Gram matches the closed Jacobian": all(
            g["gram"]["finest_relative_error"] < 1e-6 and g["gram"]["improves"]
            for g in geometry),
        "Q1 independent analytic route agrees to 1e-10": all(
            max(g["analytic_route"]["single_factor_relative_error"],
                g["analytic_route"]["joint_jacobian_relative_error"]) < 1e-10
            for g in geometry),
        "Q1 density is invariant under a common SO(3) frame change": all(
            g["covariance_exchange"]["common_frame_covariance"] < 1e-10
            for g in geometry),
        "Q3 level sets contain genuinely distinct histories": all(
            controls[k]["chord_separation"] > 1e-3
            for k in ("reflection", "signed_product", "absolute_product")),
        "Q1 closed form matches split quadrature": all(
            r["closed_form_error"] < 1e-12 for r in grids),
        "Q1 uniform sector-probability grids agree at 2048": all(
            g["sector_grids"]["final_grid_error"] < 1e-6 for g in geometry),
        "Q1 finite windows converge to the coarea limit": all(
            w["final_window_discrepancy"] < 2e-3 and w["monotone"] for w in windows),
        "Q1 punctures are exactly -u and -w with |dD/dpsi| = |q|": all(
            (p[f]["has_punctures"] == 0.0
             or max(p[f]["slope_minus_q"], p[f]["puncture_is_minus_u"],
                    p[f]["puncture_is_minus_w"]) < 1e-12)
            for p in punctures for f in ("factor1", "factor2")),
        "Q1 excised mass follows the analytic 2 eta^2 law": all(
            p[f"excision{i}"]["relative_error_vs_2eta2"] < 1e-3
            for p in punctures for i in (1, 2)),
        "Q1 joint excluded fraction is bounded and vanishing": all(
            row["joint_excluded_fraction"] < 1e-3
            for p in punctures for row in p["joint_excluded"]),
        "Q3 reflection preserves the pair statistic and the density":
            controls["reflection"]["statistic_match"] < 1e-12
            and controls["reflection"]["reference_density_gap"] < 1e-10,
        "Q3 same signed product separates the cubic weight":
            controls["signed_product"]["statistic_match"] < 1e-12
            and controls["signed_product"]["reference_density_gap"] < 1e-10
            and controls["signed_product"]["cubic_gap"] > 1e-6,
        "Q3 same absolute product, opposite sign separates the cubic":
            controls["absolute_product"]["statistic_match"] < 1e-12
            and controls["absolute_product"]["reference_density_gap"] < 1e-10
            and controls["absolute_product"]["cubic_gap"] > 1e-6,
        "Q2 the generic closure rule is rank one on a union":
            bool(audit["rank_one_demonstration"]["rank_one_cancellation_demonstrated"]),
        "Q2 based-loop additivity needs a common base point":
            audit["based_loop_scope"]["same_base_additivity_residual"] < 1e-12
            and audit["based_loop_scope"]["distinct_base_noncommutativity"] > 1e-3,
    }
    regressions = {
        "product marginals (structural)": max(
            g["product_marginal_error"] for g in geometry) < 1e-10,
        "three-factor associativity (structural)": max(
            g["associativity_regression"] for g in geometry) < 1e-10,
        "joint window factorisation (structural)": max(
            w["joint_factorisation_regression"] for w in windows) < 1e-10,
        "joint Gram off-diagonal (structural)": max(
            g["gram"]["max_offdiagonal_gram"] for g in geometry) == 0.0,
        "copy exchange (structural)": max(
            g["covariance_exchange"]["copy_exchange"] for g in geometry) < 1e-14,
    }
    checks = {k: bool(v) for k, v in checks.items()}
    regressions = {k: bool(v) for k, v in regressions.items()}

    report = {
        "public_preregistration": jc.PUBLIC_PREREG,
        "seed": jc.SEED,
        "conditioning": "CHOSEN: product Haar with separate phase windows; "
                        "round 8 established that phase conditioning is not "
                        "forced by the zero set",
        "geometry": geometry, "windows": windows, "punctures": punctures,
        "grids": grids, "level_set_controls": controls,
        "repository_audit": audit,
        "checks": checks, "structural_regressions": regressions,
        "checks_passed": all(checks.values()) and all(regressions.values()),
        # a failed structural regression signals a code fault, so it blocks
        # the verdict too rather than only the checks_passed flag
        "verdict": jc.verdict({**checks, **regressions}, audit, controls),
        "scope": "Two independently prepared triangles under the inherited "
                 "chosen phase conditioning. Not a Born rule, not an "
                 "operational readout, not a Hilbert tensor product.",
        "weight_selection": "NOT_DERIVED",
        "operational_source_readout": "NOT_DERIVED",
    }
    return report


def render(report):
    v = report["verdict"]
    lines = [
        "# Round 10: joint closure and sufficient data", "",
        report["scope"], "",
        f"Public freeze: `{report['public_preregistration']}`. Seed `{report['seed']}`.", "",
        "Conditioning inherited from round 8 is **chosen**, not derived.", "",
        "| verdict field | value |", "|---|---|",
        f"| reference composition | **{v['reference_composition']}** |",
        f"| reduction of specified rules | **{v['reduction_of_specified_rules']}** |",
        f"| additional physical rule | **{v['additional_physical_rule']}** |",
        f"| consequence for selection | **{v['consequence_for_selection']}** |",
        "", f"Reduction scope: {v['reduction_scope']}.", "",
        "## Level-set controls at `t1 = t2 = 1`", "",
        "| control | statistic match | reference density gap | cubic weight gap |",
        "|---|---:|---:|---:|",
    ]
    for key, name in (("reflection", "reflection `psi -> -psi` (pair)"),
                      ("signed_product", "`(1,2)` vs `(sqrt2,sqrt2)` (product)"),
                      ("absolute_product", "`(1/4,1)` vs `(-1/4,1)` (absolute product)")):
        c = report["level_set_controls"][key]
        lines.append(f"| {name} | {c['statistic_match']:.3e} | "
                     f"{c['reference_density_gap']:.3e} | {c['cubic_gap']:.12g} |")
    lines += ["", "## Repository audit for a joint weight rule", "",
              "| module | applies to disconnected pairs | supplies a joint weight rule |",
              "|---|---|---|"]
    for e in report["repository_audit"]["entries"]:
        lines.append(f"| `{e['module']}` | {e['applies_to_disconnected_pairs']} | "
                     f"{e['supplies_joint_weight_rule']} |")
    lines += ["", f"Search scope: {report['repository_audit']['search_scope']}.", "",
              "| criterion | pass |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in report["checks"].items()]
    lines += ["", "Structural regressions (guaranteed by the product "
              "construction; not independent evidence):", ""]
    lines += [f"- {k}: {v}" for k, v in report["structural_regressions"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/{len(report['checks'])} "
              f"required checks.", "",
              "No `Phi` is selected, no Born rule is derived and no operational "
              "source-local readout is constructed.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    report = run_probe(progress=lambda s: print(s, flush=True))
    summary = render(report)
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "probe.json").write_text(
            json.dumps(report, indent=2, allow_nan=False, default=float) + "\n")
        (args.output_dir / "probe.md").write_text(summary)
    print(summary)
    return 0 if report["checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
