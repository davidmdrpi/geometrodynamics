"""The matching map from the recorded BAM geometry to the seam cells.

Freeze ``3d3ee38``. Reports separate verdict fields and exits nonzero when
any required check fails.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from geometrodynamics.bulk import klein_handle_matching as km
from geometrodynamics.bulk import klein_handle_response as kh

PREREG = "3d3ee38"
SCOPE = ("Part I: the repository's recorded maps on the bulk mouth S^3 acting "
         "on fibre-charged functions. Part II: a lattice multiplet on the "
         "handle cover with a sign around the doubled loop. Part III: the Pin "
         "sign against the trapped-flux sign in round 14's model. Classical "
         "throughout. Which field BAM's matter is, and P against PT, are not "
         "decided here.")
ZERO, EXACT, FIT = 1e-9, 1e-12, 1e-10
SIZES = {"primary": km.PRIMARY, "replicate": km.REPLICATE}
NAMES = ("i", "sigma_x", "sigma_y", "sigma_z")


def measure_sphere() -> dict:
    rng = np.random.default_rng(km.SEED)
    return {"antipode": km.antipode_report(rng),
            "two_spheres": km.two_spheres_report(rng),
            "transport": km.transport_report(rng)}


def measure_lattice(config: dict) -> dict:
    net = kh.Network(**config)
    rng = np.random.default_rng(km.SEED + 1)
    return {"config": config, "sites": net.size,
            "pr129": km.pr129_cell_report(net),
            "consistency": km.consistency_report(net, rng),
            "seam_charges": km.seam_charge_report(net),
            "dynamics": km.transport_seam_dynamics(net, rng),
            "labels": km.label_independence_report(net, rng)}


def evaluate_sphere(m: dict) -> dict:
    a, t, x = m["antipode"], m["two_spheres"], m["transport"]
    minus_i_sigma_y = np.array([[0.0, -1.0], [1.0, 0.0]])
    return {
        "antipodal_map": {
            "M1 antipode is the fibre translation by pi": (
                a["is_fibre_translation_by_pi"] < EXACT
                and a["fibre_commutator"] < EXACT),
            "M1 it fixes the Hopf base": a["hopf_base_moved"] < EXACT,
            "M1 it acts on charge k as the scalar (-1)^k": (
                a["bidegrees"] == 15 and a["scalar_error"] < FIT
                and a["fit_residual"] < FIT),
        },
        "two_spheres": {
            "M2 antipode preserves the spatial brane mouth": (
                t["equator_real_part"] < EXACT
                and t["antipode_leaves_equator"] < EXACT),
            "M2 Hopf map is two-to-one on the spatial brane mouth": (
                t["hopf_two_to_one_error"] < EXACT
                and t["hopf_formula_error"] < EXACT),
            "M2 transport does not preserve the spatial brane mouth": (
                t["sigma_leaves_equator"] > 0.5),
            "M2 transport covers the Hopf-base antipode and reverses the "
            "fibre": (t["sigma_covers_base_antipode"] < EXACT
                      and t["sigma_reverses_fibre"] < EXACT),
            "M2 transport is the recorded rotation with square -I": (
                t["sigma_squared_is_antipode"] < EXACT
                and t["sigma_matches_recorded_matrix"] < EXACT
                and abs(t["sigma_determinant"] - 1) < EXACT),
        },
        "transport_on_fibre_modes": {
            "M3 pullback sends charge k to charge -k": (
                x["bidegrees"] == 15 and x["fit_residual"] < FIT
                and x["smallest_same_charge_residual"] > 0.1),
            "M3 A squares to (-1)^k": x["square_error"] < FIT,
            "M3 on the doublet A is J up to the Pin sign": (
                x["doublet_basis"] == ["z1", "z2"]
                and x["doublet_imaginary_part"] < FIT
                and float(np.abs(np.array(x["doublet_matrix"])
                                 - minus_i_sigma_y).max()) < FIT),
        },
    }


def evaluate_lattice(m: dict) -> dict:
    size = m["sites"]
    con = {(r["components"], r["loop_sign"]): r for r in m["consistency"]}
    seam = {(r["internal"], r["s_t"]): r for r in m["seam_charges"]}
    dyn, lab = m["dynamics"], m["labels"]

    def zero(row, names):
        return all(row[n]["max"] < EXACT for n in names)

    def balanced(row, names):
        return all(row[n]["positive"] == row[n]["negative"] == 2 * size
                   for n in names)

    spin = NAMES[1:]
    return {
        "pr129_cell": {
            "M4 charge is deck odd and vanishes for both twists": all(
                t["charge_on_fixed"] < EXACT and t["charge_parity_error"] < ZERO
                and 2 * t["fixed_dim"] == t["field_dim"]
                and t["omega_on_fixed"] < ZERO for t in m["pr129"]["twists"])
            and len(m["pr129"]["twists"]) == 2,
        },
        "spinor_seam_consistency": {
            "M5 D squares to B^2 W in all four cases": (
                len(con) == 4 and all(
                    r["square_error"] < EXACT and r["mouth_square_error"] < EXACT
                    for r in con.values())),
            "M5 H is invariant under D in all four cases": all(
                r["energy_invariance_error"] < EXACT for r in con.values()),
            "M5 singlet is admissible only on the periodic loop": (
                2 * con[(1, 1)]["fixed_dim"] == con[(1, 1)]["field_dim"]
                and con[(1, -1)]["fixed_dim"] == 0),
            "M5 doublet with J is admissible only on the antiperiodic loop":
                (2 * con[(2, -1)]["fixed_dim"] == con[(2, -1)]["field_dim"]
                 and con[(2, 1)]["fixed_dim"] == 0),
        },
        "transport_seam_charges": {
            "M6 every seam is an involution with the predicted symplectic "
            "sign": len(seam) == 4 and all(
                r["involution_error"] < EXACT
                and r["symplectic_sign_error"] < EXACT
                and all(r[n]["parity_error"] < ZERO for n in NAMES)
                for r in seam.values()),
            "M6 J seam, time reversed: charge survives, spin vanishes": (
                balanced(seam[("J", -1)], ("i",))
                and zero(seam[("J", -1)], spin)
                and seam[("J", -1)]["omega_on_fixed"] < ZERO),
            "M6 J seam, time kept: charge vanishes, spin survives": (
                zero(seam[("J", 1)], ("i",))
                and balanced(seam[("J", 1)], spin)),
            "M6 linear seam, time reversed: both vanish": zero(
                seam[("identity", -1)], NAMES),
            "M6 linear seam, time kept: both survive": balanced(
                seam[("identity", 1)], NAMES),
            "M6 responding histories satisfy the seam relation": (
                len(dyn["rows"]) == 4 and all(
                    r["initial_fixed_error"] < EXACT and r["seam_defect"] < ZERO
                    and r["max_mouth"] > 1e-3 for r in dyn["rows"])),
            "M6 fixed space is a Lagrangian half": (
                2 * dyn["fixed_dim"] == dyn["phase_space_dim"]
                and dyn["omega_on_fixed"] < ZERO),
            "M6 charge and spin are conserved and total spin is zero": all(
                all(v < ZERO for v in r["drift"].values())
                and r["total_spin"] < ZERO and abs(r["charge"]) > 1e-6
                for r in dyn["rows"]),
            "M6 sheet charge obeys Q1(t) + Q1(-t) = Q": all(
                r["sheet_relation_error"] < ZERO
                and r["sheet_charge_variation"] > 1e-6 for r in dyn["rows"]),
        },
        "flux_sign_vs_pin_label": {
            "M7 quarter rotation is a symmetry that leaves the flux alone": (
                lab["quarter_preserves_energy"] < EXACT
                and lab["quarter_preserves_flux"] < EXACT),
            "M7 quarter rotation flips the Pin sign": (
                lab["quarter_flips_pin_sign"] < EXACT
                and lab["quarter_commutes_with_deck"] > 1.0
                and dyn["minus_j_is_quarter_rotated_j"] < EXACT),
            "M7 conjugation is a symmetry that flips the flux and keeps the "
            "Pin sign": (lab["conjugation_preserves_energy"] < EXACT
                         and lab["conjugation_flips_flux"] < EXACT
                         and lab["conjugation_commutes_with_deck"] < EXACT),
            "M7 both map admissible histories to admissible histories": (
                lab["base_seam_defect"] < ZERO
                and lab["rotated_seam_defect_opposite_pin"] < ZERO
                and lab["conjugated_seam_defect"] < ZERO
                and lab["rotated_flux_change"] < ZERO
                and lab["conjugated_flux_sum"] < ZERO
                and lab["flux_scale"] > 1e-3),
            "M7 the rotated history is not admissible for the original Pin "
            "sign": lab["rotated_seam_defect_same_pin"] > 1e-3,
        },
    }


def run_probe(progress=lambda s: None):
    progress("Part I: recorded maps on S^3")
    sphere = measure_sphere()
    groups_values = evaluate_sphere(sphere)
    lattice, per_size = {}, {}
    for name, config in SIZES.items():
        progress(f"Parts II-III, {name}: {config}")
        lattice[name] = measure_lattice(config)
        per_size[name] = evaluate_lattice(lattice[name])
    checks = {c: bool(v) for g in groups_values.values() for c, v in g.items()}
    groups = {f: list(g) for f, g in groups_values.items()}
    for field, names in per_size["primary"].items():
        groups[field] = list(names)
        for check in names:
            checks[check] = all(bool(per_size[s][field][check]) for s in SIZES)
    return {
        "public_preregistration": PREREG, "scope": SCOPE, "seed": km.SEED,
        "sphere": sphere, "lattice": lattice,
        "checks": checks, "check_groups": groups,
        "checks_passed": all(checks.values()),
        "verdict": km.verdict(checks, groups),
    }


def render(report):
    v = report["verdict"]
    sph, prim = report["sphere"], report["lattice"]["primary"]
    lines = ["# The matching map", "", report["scope"], "",
             f"Public freeze: `{report['public_preregistration']}`. Seed "
             f"`{report['seed']}`. Lattice checks must hold at both sizes.",
             "", "| verdict field | value |", "|---|---|"]
    lines += [f"| {k} | **{v[k]}** |" for k in km.VERDICT_FIELDS]
    if v["failed_checks"] or v["missing_checks"]:
        lines += ["", "**Outside tolerance:** "
                  + ", ".join(f"`{c}`" for c in v["failed_checks"]),
                  "", "**Not run:** "
                  + ", ".join(f"`{c}`" for c in v["missing_checks"])]
    t, x = sph["two_spheres"], sph["transport"]
    lines += ["", "## Part I: the two recorded maps", "",
              "| | antipode `-I_4` | transport `sigma` |", "|---|---|---|",
              "| on the fibre | translation by pi | reversal |",
              "| on the Hopf base | identity | antipode |",
              f"| on the spatial brane mouth | antipode | leaves it "
              f"(`{t['sigma_leaves_equator']:.3f}`) |",
              "| on charge k | scalar `(-1)^k` | antilinear, square `(-1)^k` |",
              "", f"Doublet matrix of `A` in the basis `(z1, z2)`: "
              f"`{np.round(x['doublet_matrix'], 12).tolist()}`.", "",
              "## Part II: when the transport seam admits a field "
              "(primary size)", "",
              "| field | B | loop sign W | D squared | fixed dim |",
              "|---|---|---:|---:|---:|"]
    lines += [f"| {'singlet' if r['components'] == 1 else 'doublet'} | "
              f"{r['internal']} | {r['loop_sign']:+d} | "
              f"{r['predicted_square']:+d} | {r['fixed_dim']} |"
              for r in prim["consistency"]]
    lines += ["", "| seam | time | charge signature | spin signature |",
              "|---|---|---|---|"]
    for r in prim["seam_charges"]:
        spin = r["sigma_x"]
        lines.append(
            f"| {r['internal']}, W = {r['loop_sign']:+d} | "
            f"{'reversed' if r['s_t'] == -1 else 'kept'} | "
            f"({r['i']['positive']}, {r['i']['negative']}) | "
            f"({spin['positive']}, {spin['negative']}) |")
    lab = prim["labels"]
    lines += ["", "## Part III: Pin sign against flux sign", "",
              f"Quarter rotation: flips the Pin sign, flux unchanged "
              f"(`{lab['rotated_flux_change']:.1e}`). Conjugation: keeps the "
              f"Pin sign, flux reversed (`{lab['conjugated_flux_sum']:.1e}`). "
              f"A rotated history tested against the original Pin sign "
              f"misses by `{lab['rotated_seam_defect_same_pin']:.3f}`.", "",
              "| Required check | Pass |", "|---|---|"]
    lines += [f"| {k} | {val} |" for k, val in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/"
              f"{len(report['checks'])} required checks.", "",
              "Which cell BAM occupies depends on the field type and on the "
              "time extension; neither is decided here.", ""]
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
