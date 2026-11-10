"""A localized, responding Klein handle. Freeze ``5d15963``.

Reports separate verdict fields and exits nonzero when any required check
fails. Not Stage 1 of the boundary-pairing protocol.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from geometrodynamics.bulk import klein_handle_response as kh

PREREG = "5d15963"
SCOPE = ("A lattice ring with one handle whose seam may reverse time and "
         "fibre; one waist coordinate per tube solved from the same "
         "Hamiltonian as the field; effective oscillator support; static "
         "prescribed potential. Classical throughout. No probability rule, "
         "spinor matching or statement about the cell BAM occupies.")
ZERO, EXACT = 1e-9, 1e-12
SIZES = {"primary": kh.PRIMARY, "replicate": kh.REPLICATE}


def measure(config: dict) -> dict:
    net = kh.Network(**config)
    rng = np.random.default_rng(kh.SEED)
    out = {
        "config": config, "cover_sites": net.size,
        "H1": kh.localized_report(net, rng),
        "H2": kh.local_bracket_report(net, rng),
        "H3": kh.response_structure_report(net, rng),
        "H4": kh.order_counting_report(net, rng),
        "H5": kh.seam_cell_report(net),
        "H6": kh.potential_parity_report(net, rng),
        "H7": kh.pairing_report(net, rng),
    }
    out["supplementary"] = kh.supplementary_report(
        net, np.random.default_rng(kh.SEED + 7))
    return out


def evaluate(m: dict) -> dict:
    """Required checks for one network size, grouped by verdict field."""
    size = m["cover_sites"]
    h1, h2, h3, h4, h5, h6, h7 = (m[k] for k in
                                  ("H1", "H2", "H3", "H4", "H5", "H6", "H7"))
    cell = {(r["s_t"], r["s_f"]): r for r in h5}
    both = cell[(-1, -1)]
    return {
        "localized_structure": {
            "H1 deck map is an antisymplectic involution": all(
                t["involution_error"] < ZERO and t["antisymplectic_error"] < ZERO
                for t in h1["twists"]),
            "H1 fixed space is a Lagrangian half": all(
                t["fixed_dim"] == size and t["field_dim"] == 2 * size
                and t["omega_on_fixed"] < ZERO for t in h1["twists"]),
            "H1 frozen histories satisfy the seam relation": all(
                t["frozen_seam_defect"] < ZERO for t in h1["twists"]),
            "H1 even and odd spectra are the ordinary and twisted handle": (
                h1["swap_commutator"] < ZERO
                and h1["spectrum_error"]["even"] < ZERO
                and h1["spectrum_error"]["odd"] < ZERO
                and h1["even_odd_spectra_differ"] > 1e-3),
            "H1 time-preserving seam keeps a symplectic subspace": (
                h1["control"]["symplectic_error"] < ZERO
                and h1["control"]["fixed_dim"] == size
                and h1["control"]["restricted_rank"] == size),
            "H3 scheme is reversed by the deck map": (
                h3["reversibility_error"] < ZERO),
            "H3 responding histories satisfy the seam relation at every "
            "amplitude": (
                [r["amplitude"] for r in h3["amplitudes"]]
                == list(kh.AMPLITUDES) + [kh.STRESS_AMPLITUDE]
                and all(r["seam_defect"] < ZERO for r in h3["amplitudes"])
                and h3["moving_mouth_seam_defect"] < ZERO),
            "H3 tangent space is a Lagrangian half": (
                h3["tangent_dim"] * 2 == h3["phase_space_dim"]
                and h3["tangent_rank"] == h3["tangent_dim"]
                and h3["tangent_lagrangian_error"] < ZERO),
            "H3 linearized tangent agrees with finite differences": (
                h3["tangent_fd_error"] < 1e-6),
            "H3 energy closes to 1e-5 for amplitudes up to 1/8": all(
                r["energy_drift"] < 1e-5 for r in h3["amplitudes"]
                if r["amplitude"] <= 0.125),
            "H4 mouth response is second order in amplitude": all(
                abs(s - 2) < 0.15 for s in h4["q_slopes"])
            and len(h4["q_slopes"]) == 3,
            "H4 field correction is third order in amplitude": all(
                abs(s - 3) < 0.3 for s in h4["deviation_slopes"])
            and len(h4["deviation_slopes"]) == 3,
            "H4 the two mouths respond differently and mirror each other":
                all(r["mouth_asymmetry"] > 0 and r["mouth_mirror_error"] < ZERO
                    for r in h4["rows"]),
        },
        "local_bracket": {
            "H2 sheet data are free at the mirror time": (
                h2["sheet_data_rank"] == h2["sheet_data_dim"] == size),
            "H2 sheet Wronskian is symmetric about the mirror time": (
                h2["mirror_symmetry_error"] < ZERO
                and h2["cover_wronskian"] < ZERO),
            "H2 sheet Wronskian is not conserved": (
                h2["relative_variation"] > 1e-6),
        },
        "neutral_sector_pairing": {
            "H7 mechanical response generates no imaginary part": (
                len(h7["mechanical_leak"]) == len(kh.AMPLITUDES) + 1
                and all(v < EXACT for v in h7["mechanical_leak"])),
            "H7 the mouths were responding in those runs": all(
                v > 1e-7 for v in h7["mechanical_mouth_motion"]),
        },
        "charged_sector_pairing": {
            "H5 each deck map is an involution with the predicted "
            "symplectic sign": all(
                r["involution_error"] < ZERO
                and r["symplectic_sign_error"] < ZERO for r in h5)
            and set(cell) == set(kh.CELLS),
            "H5 charge density has parity s_t s_f": all(
                r["charge_parity_error"] < ZERO for r in h5),
            "H5 fixed space is half, Lagrangian iff the seam reverses time":
                all(2 * r["fixed_dim"] == r["phase_space_dim"]
                    and ((r["omega_on_fixed"] < ZERO) if r["s_t"] == -1
                         else (r["omega_rank_on_fixed"] == r["fixed_dim"]))
                    for r in h5),
            "H5 charge vanishes identically in the two mixed cells": all(
                cell[c]["charge_on_fixed"] < EXACT
                for c in ((1, -1), (-1, 1))),
            "H5 twist is a quarter rotation in the antilinear cells only":
                all((r["twist_is_a_rotation_error"] < EXACT)
                    if r["s_f"] == -1
                    else (r["twist_is_a_rotation_error"] > 1.0) for r in h5),
            "H6 matching potential parity admits the seam": (
                h6["matching_seam_defect"] < ZERO),
            "H6 mismatched potential parity breaks it": (
                h6["mismatched_seam_defect"] > 1e-4),
            "H7 charge is nondegenerate with balanced signature": (
                both["charge_positive"] == size
                and both["charge_negative"] == size
                and both["field_fixed_dim"] == 2 * size),
            "H7 charge is the Wronskian pairing of the two sectors": (
                h7["charge_minus_wronskian"] < EXACT
                and h7["re_outside_untwisted"] < ZERO
                and h7["im_outside_twisted"] < ZERO),
            "H7 charge is conserved with potential and moving mouths": (
                h7["charge_drift"] < ZERO),
            "H7 sheet charge obeys Q1(t) + Q1(-t) = Q": (
                h7["sheet_relation_error"] < ZERO),
            "H7 potential couples the sectors at first order": all(
                abs(s - 1) < 0.1 for s in h7["potential_slopes"])
            and len(h7["potential_slopes"]) == 2,
        },
        "positive_response": {
            "H8 charge is indefinite": (
                both["charge_positive"] > 0 and both["charge_negative"] > 0),
            "H8 multiplication by i preserves no admissible field": (
                both["fix_meet_i_fix"] == 0),
        },
    }


def run_probe(progress=lambda s: None):
    measured, per_size = {}, {}
    for name, config in SIZES.items():
        progress(f"{name}: {config}")
        measured[name] = measure(config)
        per_size[name] = evaluate(measured[name])
    groups = {field: list(per_size["primary"][field])
              for field in per_size["primary"]}
    checks = {}
    for field, names in groups.items():
        for check in names:
            checks[check] = all(bool(per_size[s][field][check])
                                for s in SIZES)
    return {
        "public_preregistration": PREREG, "scope": SCOPE, "seed": kh.SEED,
        "measurements": measured,
        "checks_by_size": {s: {c: bool(v) for g in per_size[s].values()
                               for c, v in g.items()} for s in SIZES},
        "checks": checks, "check_groups": groups,
        "checks_passed": all(checks.values()),
        "verdict": kh.verdict(checks, groups),
    }


def render(report):
    v = report["verdict"]
    prim = report["measurements"]["primary"]
    repl = report["measurements"]["replicate"]
    lines = ["# A localized, responding Klein handle", "", report["scope"], "",
             f"Public freeze: `{report['public_preregistration']}`. Seed "
             f"`{report['seed']}`. Cover sites: "
             f"`{prim['cover_sites']}` and `{repl['cover_sites']}`; every "
             "required check must hold at both.", "",
             "| verdict field | value |", "|---|---|"]
    lines += [f"| {k} | **{v[k]}** |" for k in kh.VERDICT_FIELDS]
    if v["failed_checks"] or v["missing_checks"]:
        lines += ["", "**Outside tolerance:** "
                  + ", ".join(f"`{c}`" for c in v["failed_checks"]),
                  "", "**Not run:** "
                  + ", ".join(f"`{c}`" for c in v["missing_checks"])]
    lines += ["", "## Four seam types, charged field (primary size)", "",
              "| s_t | s_f | fixed dim | Omega on fixed | charge on fixed | "
              "signature | Fix meet i Fix |", "|---:|---:|---:|---|---|---|---:|"]
    for r in prim["H5"]:
        omega = ("zero" if r["omega_on_fixed"] < ZERO
                 else f"rank {r['omega_rank_on_fixed']}")
        charge = ("zero" if r["charge_on_fixed"] < EXACT else "nonzero")
        lines.append(
            f"| {r['s_t']:+d} | {r['s_f']:+d} | {r['fixed_dim']} | {omega} | "
            f"{charge} | ({r['charge_positive']}, {r['charge_negative']}) | "
            f"{r['fix_meet_i_fix']} |")
    lines += ["", "## Response (primary size)", "",
              "| amplitude | seam defect | max mouth excursion | energy drift |",
              "|---:|---:|---:|---:|"]
    lines += [f"| {r['amplitude']:.6g} | {r['seam_defect']:.1e} | "
              f"{r['max_q']:.3e} | {r['energy_drift']:.1e} |"
              for r in prim["H3"]["amplitudes"]]
    lines += ["", "Supplementary, not frozen:", ""]
    lines += [f"- amplitude `{r['amplitude']}`: seam defect "
              f"`{r['seam_defect']:.1e}`, mouth excursion `{r['max_q']:.3f}`, "
              f"barrier factor `{r['min_barrier_factor']:.3f}`, energy drift "
              f"`{r['energy_drift']:.1e}`"
              for r in prim["supplementary"]["strong_response"]]
    h4, h7 = prim["H4"], prim["H7"]
    lines += ["", f"Order counting: mouth slopes "
              f"`{[round(s, 4) for s in h4['q_slopes']]}`, field-correction "
              f"slopes `{[round(s, 4) for s in h4['deviation_slopes']]}`.", "",
              f"Sector coupling: imaginary part from mechanical response "
              f"`{max(h7['mechanical_leak']):.1e}`; from the potential, "
              f"slopes `{[round(s, 4) for s in h7['potential_slopes']]}`. "
              "On the ordinary handle the potential gives slopes "
              f"`{[round(s, 4) for s in prim['supplementary']['ordinary_handle_potential_slopes']]}`"
              " (supplementary).", "",
              "| Required check (both sizes) | Pass |", "|---|---|"]
    lines += [f"| {k} | {val} |" for k, val in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/"
              f"{len(report['checks'])} required checks.", "",
              "An indefinite conserved pairing is not a probability, a norm "
              "or a spinor pairing. Which seam type BAM has is not "
              "established here.", ""]
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
