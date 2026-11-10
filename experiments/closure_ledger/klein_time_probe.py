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
REVIEWED = "5b79975"
SCOPE = ("One linear zonal conformal scalar on the Einstein static universe, "
         "global quotients with caps placed by hand. Classical throughout; "
         "no hbar, probability rule or counting function is derived. "
         "Corrected after review: the frozen 'no free field' reading holds "
         "only in the frozen truncation.")
COUNTEREXAMPLE_TEXT = (
    "u = sin(48 chi) sin(48 t) satisfies both untwisted caps at "
    "tau = 7 pi / 48; every rational tau / pi = p / q admits the untwisted "
    "mode k = q. This refutes the frozen label "
    "GENERIC_CLOSED_HISTORY_HAS_NO_FREE_FIELD as a continuum statement.")


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
    sector = [kt.sector_action_report(lat, m, eta)
              for m in kt.CAPS_P1 for eta in (1, -1)]

    progress("P3 two caps")
    caps = [kt.two_cap_report(lat, m) for m in kt.TABLE_P3]
    closure_control = [kt.closure_control_report(lat, m) for m in kt.TABLE_P3]

    progress("P4 closed-history boundary-value problem, frozen truncation")
    kernels = [kt.kernel_report(lat, m) for m in kt.TABLE_P3]
    slaved = kt.slaved_field_report(lat, 7)
    fredholm = kt.fredholm_report(lat, 24, rng)

    progress("review: fixed-geometry refinement and continuum counterexample")
    refinement = kt.refinement_report()
    counter = [kt.continuum_counterexample(7, 48, False),
               kt.continuum_counterexample(7, 48, True)]
    rational_scan = [kt.continuum_counterexample(p, q, False)
                     for p, q in ((1, 3), (2, 5), (3, 7), (5, 11), (4, 9))]

    progress("review: what the attachment adds to temporal closure")
    response = kt.response_control_report(lat)

    progress("diagnostic: continuum divisors (not a required check)")
    golden = (5 ** 0.5 - 1) / 2
    divisors = {name: kt.continuum_gap(alpha) for name, alpha in {
        "7/48": 7 / 48, "golden": golden, "sqrt2 - 1": 2 ** 0.5 - 1,
        "pi - 3": float(np.pi - 3), "0.110001": 0.110001}.items()}

    groups = {
        "phase_space": {
            "lattice propagation is exactly periodic": period_error < zero,
            "deck map is an involution and reverses the Wronskian": all(
                r["involution_error"] < zero
                and r["antisymplectic_error"] < zero for r in halving),
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
        },
        "complex_structure": {
            "j is a complex structure commuting with evolution": all(
                r["j_squared_error"] < zero
                and r["commutes_with_evolution"] < zero for r in conj),
            "deck map anticommutes with j": all(
                r["antilinear_error"] < zero for r in conj),
            "j maps admissible onto twisted": all(
                r["j_plus_in_minus"] < zero for r in conj),
            "metric is symmetric definite and h is real on admissible": all(
                r["metric_symmetry_error"] < zero
                and r["metric_definite_min"] > firm
                and r["hermitian_imaginary_on_plus"] < zero for r in conj),
            "j preserves no admissible field of a fixed twist": all(
                r["dim_plus_meet_j_plus"] == 0 and r["dim"] == lat.dim
                for r in sector),
        },
        "phase_closure": {
            "frozen table equals the closure formula": all(
                kt.TABLE_P3[r["m"]] == (r["predicted_untwisted"],
                                        r["predicted_twisted"]) for r in caps),
            "untwisted two-cap modes are exactly the predicted set": all(
                r["pp"]["modes"] == kt.TABLE_P3[r["m"]][0]
                and r["pp"]["dim"] == len(kt.TABLE_P3[r["m"]][0])
                for r in caps),
            "twisted two-cap modes are exactly the predicted set": all(
                r["pm"]["modes"] == kt.TABLE_P3[r["m"]][1]
                and r["pm"]["dim"] == len(kt.TABLE_P3[r["m"]][1])
                for r in caps),
            "only the product of twists matters": all(
                r["mm"] == r["pp"] and r["mp"] == r["pm"] for r in caps),
            "periodic operator reproduces the two-cap count": all(
                k["nullity"] == 2 * len(kt.TABLE_P3[k["m"]][0])
                and k["invariant_nullity"] == len(kt.TABLE_P3[k["m"]][0])
                and k["deck_commutator"] < zero
                and k["operator_symmetry"] < zero for k in kernels),
        },
        "field_freedom": {
            "frozen truncation at m = 7 is uniquely solvable": (
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
        },
        "continuum_free_field": {
            "fixed-geometry refinement matches the continuum survivors": all(
                (r["untwisted"], r["twisted"]) == kt.REFINEMENT_EXPECTED[r["n"]]
                and r["untwisted"] == r["continuum_untwisted"]
                and r["twisted"] == r["continuum_twisted"]
                for r in refinement)
            and [r["n"] for r in refinement] == list(kt.REFINEMENT_SIZES),
            "explicit continuum modes satisfy both caps at tau = 7 pi / 48":
                all(c["exists"] and c["amplitude"] > 0.5
                    and c["cap_0_defect"] < zero
                    and c["cap_tau_defect"] < zero for c in counter)
                and (counter[0]["k"], counter[1]["k"]) == (48, 24),
            "every sampled rational spacing admits the mode k = q": all(
                c["exists"] and c["k"] == c["q"] and c["amplitude"] > 0.5
                and c["cap_0_defect"] < zero and c["cap_tau_defect"] < zero
                for c in rational_scan),
        },
        "attachment_contribution": {
            "orientable closed time selects the same harmonics, both "
            "quadratures": all(
                c["orientable_periodic_modes"] == kt.TABLE_P3[c["m"]][0]
                and c["orientable_antiperiodic_modes"] == kt.TABLE_P3[c["m"]][1]
                and c["orientable_periodic_dim"]
                == 2 * len(kt.TABLE_P3[c["m"]][0])
                and c["orientable_antiperiodic_dim"]
                == 2 * len(kt.TABLE_P3[c["m"]][1])
                and c["klein_untwisted_dim"] == len(kt.TABLE_P3[c["m"]][0])
                and c["klein_twisted_dim"] == len(kt.TABLE_P3[c["m"]][1])
                and c["antiperiodic_deck_commutator"] < zero
                and c["antiperiodic_deck_involution"] < zero
                for c in closure_control),
            "reciprocal time-symmetric response exists without the deck": (
                not response["source_is_deck_invariant"]
                and response["unrestricted_smallest_singular"] > firm
                and response["unrestricted_green_symmetry"] < zero
                and response["time_mirror_defect_without_deck"] < zero
                and response["field_scale"] > 0.1),
            "attachment restricts sources to the deck-invariant half": (
                response["source_space_dim"] == 7 * lat.dim
                and response["deck_invariant_source_dim"]
                == 12 + 3 * lat.dim),
            "attachment halves the solvability conditions at resonance": (
                response["conditions_unrestricted"] == 22
                and response["conditions_on_invariant_sources"] == 11),
        },
    }
    checks = {name: bool(value) for group in groups.values()
              for name, value in group.items()}
    names = {field: list(group) for field, group in groups.items()}
    counterexamples = {}
    if all(groups["continuum_free_field"].values()):
        counterexamples["continuum_free_field"] = COUNTEREXAMPLE_TEXT
    return {
        "public_preregistration": PREREG, "reviewed_commit": REVIEWED,
        "scope": SCOPE, "seed": kt.SEED,
        "lattice_n": lat.n, "period_error": period_error,
        "halving": halving, "orientable_control": control,
        "direct_admissibility": direct, "conjugation": conj,
        "sector_action": sector,
        "two_caps": caps, "closure_control": closure_control,
        "periodic_kernels": kernels,
        "slaved_field": slaved, "fredholm": fredholm,
        "refinement": refinement, "continuum_counterexamples": counter,
        "rational_scan": rational_scan, "response_control": response,
        "continuum_divisor_diagnostic": divisors,
        "frozen_labels_narrowed": dict(kt.FROZEN_LABELS),
        "checks": checks, "check_groups": names,
        "checks_passed": all(checks.values()),
        "verdict": kt.verdict(checks, names, counterexamples),
    }


def render(report):
    v = report["verdict"]
    lines = ["# Klein time attachment of a classical field", "",
             report["scope"], "",
             f"Public freeze: `{report['public_preregistration']}`. "
             f"Corrected after review of `{report['reviewed_commit']}`. "
             f"Lattice `N = {report['lattice_n']}`. Seed `{report['seed']}`.",
             "", "| verdict field | value |", "|---|---|"]
    lines += [f"| {k} | **{v[k]}** |" for k in kt.VERDICT_FIELDS]
    for field, text in v["counterexamples"].items():
        lines += ["", f"**Counterexample ({field}).** {text}"]
    lines += ["", "Frozen labels narrowed after review: "
              + "; ".join(f"`{old}` -> `{kt.CONFIRMED[f]}`"
                          for f, old in report["frozen_labels_narrowed"].items())
              + "."]
    if v["failed_checks"] or v["missing_checks"]:
        lines += ["", "**Checks outside tolerance:** "
                  + ", ".join(f"`{c}`" for c in v["failed_checks"]),
                  "", "**Checks not run:** "
                  + ", ".join(f"`{c}`" for c in v["missing_checks"])]
    lines += ["", "## P1 Phase space", "",
              "| cap | eta | dim + | dim - | max Omega on + | pairing min s.v. |",
              "|---:|---:|---:|---:|---:|---:|"]
    lines += [f"| {r['cap']} | {r['eta']:+d} | {r['dim_plus']} | "
              f"{r['dim_minus']} | {r['omega_on_plus']:.1e} | "
              f"{r['pairing_min_singular']:.4f} |" for r in report["halving"]]
    c = report["orientable_control"]
    lines += ["", f"Orientable control: fixed dimension `{c['dim_fixed']}`, "
              f"restricted Wronskian rank `{c['restricted_rank']}`.",
              "", "## P3 Two caps against orientable closed time", "",
              "| m | tau/pi | harmonics (untwisted / twisted) | "
              "orientable dim | Klein dim |", "|---:|---:|---|---|---|"]
    for r, cc in zip(report["two_caps"], report["closure_control"]):
        lines.append(
            f"| {r['m']} | {r['tau_over_pi']:.4f} | "
            f"`{list(r['pp']['modes'])}` / `{list(r['pm']['modes'])}` | "
            f"{cc['orientable_periodic_dim']} / "
            f"{cc['orientable_antiperiodic_dim']} | "
            f"{cc['klein_untwisted_dim']} / {cc['klein_twisted_dim']} |")
    lines += ["", "## Fixed physical cap separation `tau = 7 pi / 48`", "",
              "| N | m | untwisted survivors | twisted survivors |",
              "|---:|---:|---|---|"]
    lines += [f"| {r['n']} | {r['m']} | `{list(r['untwisted'])}` | "
              f"`{list(r['twisted'])}` |" for r in report["refinement"]]
    s, f, rc = (report["slaved_field"], report["fredholm"],
                report["response_control"])
    lines += ["", "## P4 Closed history, frozen truncation only", "",
              f"`m = {s['m']}`, `N = 24`: smallest singular value "
              f"`{s['smallest_singular_value']:.4f}`, Green symmetry error "
              f"`{s['green_symmetry_error']:.1e}`.", "",
              f"Without the deck constraint, a source that is not deck "
              f"invariant already has a unique response, time-symmetric to "
              f"`{rc['time_mirror_defect_without_deck']:.1e}`. The attachment "
              f"admits `{rc['deck_invariant_source_dim']}` of "
              f"`{rc['source_space_dim']}` source dimensions and reduces the "
              f"resonant solvability conditions from "
              f"`{rc['conditions_unrestricted']}` to "
              f"`{rc['conditions_on_invariant_sources']}`.", "",
              f"Resonant `m = {f['m']}`: nullity `{f['nullity']}`, residual "
              f"`{f['residual_orthogonal_source']:.1e}` for an orthogonal "
              f"source and `{f['residual_generic_source']:.4f}` for a generic "
              f"one.", "",
              "## Continuum divisors (diagnostic, not a required check)", "",
              "| tau/pi | min gap, k <= 1e5 | at k |", "|---|---:|---:|"]
    lines += [f"| {name} | {d['min_gap']:.6f} | {d['at_k']} |"
              for name, d in report["continuum_divisor_diagnostic"].items()]
    lines += ["", "| Required check | Pass |", "|---|---|"]
    lines += [f"| {k} | {val} |" for k, val in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/"
              f"{len(report['checks'])} required checks.", "",
              "These are classical kinematic structures. No quantum rule, "
              "lepton index, counting function or sector pairing follows "
              "from them.", ""]
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
