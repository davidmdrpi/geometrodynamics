"""Klein time attachment of a classical field. Freeze ``bb6ff2c``.

Reports separate verdict fields and exits nonzero when any required check
fails.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from geometrodynamics.bulk import klein_time as kt

PREREG = "bb6ff2c"
SCOPE = ("One linear zonal conformal scalar on the Einstein static universe, "
         "global quotients with caps placed by hand. Classical throughout; "
         "no hbar, probability rule or counting function is derived.")


def run_probe(progress=lambda s: None):
    lat = kt.Lattice(kt.N_FROZEN)
    rng = np.random.default_rng(kt.SEED)
    size = 2 * lat.dim
    zero, firm = kt.TOL_ZERO, kt.TOL_NONDEGENERATE

    progress("P1 halving")
    halving = [kt.halving_report(lat, m, eta)
               for m in kt.CAPS_P1 for eta in (1, -1)]
    control = kt.orientable_control(lat)
    direct = [kt.direct_admissibility(lat, m, eta, rng)
              for m in (7, 24) for eta in (1, -1)]
    period_error = float(np.abs(
        np.linalg.matrix_power(lat.step, 2 * lat.n) - np.eye(size)).max())

    progress("P2 conjugation")
    conj = [kt.conjugation_report(lat, m, eta)
            for m in kt.CAPS_P1 for eta in (1, -1)]

    progress("P3 two caps")
    caps = [kt.two_cap_report(lat, m) for m in kt.TABLE_P3]

    progress("P4 closed-history boundary-value problem")
    kernels = [kt.kernel_report(lat, m) for m in kt.TABLE_P3]
    slaved = kt.slaved_field_report(lat, 7)
    fredholm = kt.fredholm_report(lat, 24, rng)

    checks = {
        "lattice propagation is exactly periodic": period_error < zero,
        "deck map is an involution and reverses the Wronskian": all(
            r["involution_error"] < zero and r["antisymplectic_error"] < zero
            for r in halving),
        "each eigenspace has half the dimension": all(
            r["dim_plus"] == lat.dim and r["dim_minus"] == lat.dim
            for r in halving),
        "Wronskian vanishes on each eigenspace": all(
            r["omega_on_plus"] < zero and r["omega_on_minus"] < zero
            for r in halving),
        "the two eigenspaces pair nondegenerately": all(
            r["pairing_min_singular"] > firm for r in halving),
        "orientable control keeps a symplectic subspace": (
            control["symplectic_error"] < zero
            and control["dim_fixed"] == 24
            and control["restricted_rank"] == 24
            and control["restricted_min_singular"] > firm),
        "explicit marching confirms admissibility": all(
            r["admissible_defect"] < zero and r["generic_defect"] > 0.1
            for r in direct),
        "j is a complex structure commuting with evolution": all(
            r["j_squared_error"] < zero and r["commutes_with_evolution"] < zero
            for r in conj),
        "deck map anticommutes with j": all(
            r["antilinear_error"] < zero for r in conj),
        "j maps admissible onto twisted": all(
            r["j_plus_in_minus"] < zero for r in conj),
        "metric is symmetric definite and h is real on admissible": all(
            r["metric_symmetry_error"] < zero
            and r["metric_definite_min"] > firm
            and r["hermitian_imaginary_on_plus"] < zero for r in conj),
        "frozen table equals the closure formula": all(
            kt.TABLE_P3[r["m"]] == (r["predicted_untwisted"],
                                    r["predicted_twisted"]) for r in caps),
        "untwisted two-cap modes are exactly the predicted set": all(
            r["pp"]["modes"] == kt.TABLE_P3[r["m"]][0]
            and r["pp"]["dim"] == len(kt.TABLE_P3[r["m"]][0]) for r in caps),
        "twisted two-cap modes are exactly the predicted set": all(
            r["pm"]["modes"] == kt.TABLE_P3[r["m"]][1]
            and r["pm"]["dim"] == len(kt.TABLE_P3[r["m"]][1]) for r in caps),
        "only the product of twists matters": all(
            r["mm"] == r["pp"] and r["mp"] == r["pm"] for r in caps),
        "periodic operator reproduces the two-cap count": all(
            k["nullity"] == 2 * len(kt.TABLE_P3[k["m"]][0])
            and k["invariant_nullity"] == len(kt.TABLE_P3[k["m"]][0])
            and k["deck_commutator"] < zero and k["operator_symmetry"] < zero
            for k in kernels),
        "generic closed history is uniquely solvable": (
            slaved["smallest_singular_value"] > firm),
        "its Green function is reciprocal and deck invariant": (
            slaved["green_symmetry_error"] < zero
            and slaved["green_deck_commutator"] < zero
            and slaved["field_invariance_error"] < zero),
        "its response is not retarded": (
            slaved["before_after_ratio_on_cap"] > 0.1
            and slaved["before_after_ratio_off_cap"] > 0.1),
        "resonant history obeys the Fredholm alternative": (
            fredholm["residual_orthogonal_source"] < 1e-8
            and abs(fredholm["residual_generic_source"]
                    - fredholm["kernel_component_norm"]) < 1e-8
            and fredholm["kernel_component_norm"] > 0.1),
    }
    names = list(checks)
    groups = {
        "phase_space": names[0:7],
        "complex_structure": names[7:11],
        "phase_closure": names[11:15],
        "field_freedom": names[15:20],
    }
    return {
        "public_preregistration": PREREG, "scope": SCOPE, "seed": kt.SEED,
        "lattice_n": lat.n, "period_error": period_error,
        "halving": halving, "orientable_control": control,
        "direct_admissibility": direct, "conjugation": conj,
        "two_caps": caps, "periodic_kernels": kernels,
        "slaved_field": slaved, "fredholm": fredholm,
        "checks": checks, "check_groups": groups,
        "checks_passed": all(checks.values()),
        "verdict": kt.verdict(checks, groups),
    }


def render(report):
    v = report["verdict"]
    lines = ["# Klein time attachment of a classical field", "",
             report["scope"], "",
             f"Public freeze: `{report['public_preregistration']}`. "
             f"Lattice `N = {report['lattice_n']}`. Seed `{report['seed']}`.",
             "", "| verdict field | value |", "|---|---|"]
    lines += [f"| {k} | **{v[k]}** |" for k in kt.VERDICT_FIELDS]
    if v["failed_checks"]:
        lines += ["", "**Required checks failed:** "
                  + ", ".join(f"`{c}`" for c in v["failed_checks"])]
    lines += ["", "## P1 Phase space", "",
              "| cap | eta | dim + | dim - | max Omega on + | pairing min s.v. |",
              "|---:|---:|---:|---:|---:|---:|"]
    lines += [f"| {r['cap']} | {r['eta']:+d} | {r['dim_plus']} | "
              f"{r['dim_minus']} | {r['omega_on_plus']:.1e} | "
              f"{r['pairing_min_singular']:.4f} |" for r in report["halving"]]
    c = report["orientable_control"]
    lines += ["", f"Orientable control: fixed dimension `{c['dim_fixed']}`, "
              f"restricted Wronskian rank `{c['restricted_rank']}`.",
              "", "## P3 Two caps", "",
              "| m | tau/pi | untwisted modes | twisted modes |",
              "|---:|---:|---|---|"]
    lines += [f"| {r['m']} | {r['tau_over_pi']:.4f} | "
              f"`{list(r['pp']['modes'])}` | `{list(r['pm']['modes'])}` |"
              for r in report["two_caps"]]
    s, f = report["slaved_field"], report["fredholm"]
    lines += ["", "## P4 Closed history", "",
              f"Generic `m = {s['m']}`: smallest singular value "
              f"`{s['smallest_singular_value']:.4f}`, Green symmetry error "
              f"`{s['green_symmetry_error']:.1e}`, before/after ratio "
              f"`{s['before_after_ratio_on_cap']:.6f}` on the cap and "
              f"`{s['before_after_ratio_off_cap']:.6f}` off it.", "",
              f"Resonant `m = {f['m']}`: nullity `{f['nullity']}`, residual "
              f"`{f['residual_orthogonal_source']:.1e}` for an orthogonal "
              f"source and `{f['residual_generic_source']:.4f}` for a generic "
              f"one (kernel component `{f['kernel_component_norm']:.4f}`).",
              "", "| Required check | Pass |", "|---|---|"]
    lines += [f"| {k} | {val} |" for k, val in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/"
              f"{len(report['checks'])} required checks.", "",
              "These are classical kinematic structures. No quantum rule, "
              "lepton index or counting function follows from them.", ""]
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
