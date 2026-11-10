"""A responding electric field on the Klein handle. Freeze ``4ab8bcd``.

Reports separate verdict fields and exits nonzero when any required check
fails. Not Stage 1 of the boundary-pairing protocol.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from geometrodynamics.bulk import klein_handle_electric as ke
from geometrodynamics.bulk import klein_handle_response as kh

PREREG = "4ab8bcd"
REVIEWED = "562af71"
WITHDRAWN = {
    "response becomes dynamically active only through trapped handle flux":
        "Admissible Klein data with x swap even, y swap odd and Pi = a = E = 0 "
        "have no charge, no field and no flux, and generate sheet charge with "
        "the handle-loop flux zero throughout. Trapped flux is sufficient, "
        "not necessary.",
    "the tangent test was done with Gauss's law imposed":
        "The frozen test perturbed along every deck-fixed direction; those "
        "violate linearized Gauss. It establishes an ambient property only.",
    "the generated imaginary part is a physical orientation":
        "A gauge transformation preserving the deck condition changes it. "
        "The gauge-invariant statements are about link currents, sheet "
        "charge, loop holonomy and the electric field.",
}
SCOPE = ("Lattice scalar electrodynamics on the cover network of a ring with "
         "one handle: charged field, link potentials, electric fields and two "
         "waist coordinates solved from one Hamiltonian. Effective oscillator "
         "support. Classical throughout; no probability rule, norm, "
         "spinor matching or statement about the seam BAM has. Corrected "
         "after review: trapped flux is sufficient but not necessary for "
         "sheet pairing, the frozen tangent test was ambient, and the "
         "imaginary part of the field is gauge dependent.")
ZERO, EXACT = 1e-9, 1e-12
SIZES = {"primary": ke.PRIMARY, "replicate": ke.REPLICATE}
MATCHED = ((1, 1), (-1, -1))


def measure(config: dict) -> dict:
    model = ke.Electric(kh.Network(**config))
    rng = np.random.default_rng(ke.SEED)
    return {
        "config": config, "sites": model.n_sites, "edges": model.n_edges,
        "phase_space_dim": model.dim,
        "deck": ke.deck_report(model, rng),
        "structure": ke.structure_report(model, rng),
        "admissibility": ke.flux_admissibility_report(model),
        "spontaneous": {"klein": ke.spontaneous_report(model, rng, -1, -1),
                        "ordinary": ke.spontaneous_report(model, rng, 1, 1)},
        "flux": {"klein": ke.flux_report(model, rng, -1, -1),
                 "ordinary": ke.flux_report(model, rng, 1, 1)},
        "supplementary": ke.supplementary_report(
            model, np.random.default_rng(ke.SEED + 7)),
        "review": review(model),
    }


def review(model) -> dict:
    """Measurements added after review of ``562af71``; not frozen."""
    rng = np.random.default_rng(ke.SEED + 11)
    return {"matter_phase": ke.matter_phase_report(model, rng),
            "constrained": ke.constrained_tangent_report(model, rng),
            "gauge": ke.gauge_report(model, rng)}


def evaluate(m: dict) -> dict:
    deck, st, adm = m["deck"], m["structure"], m["admissibility"]
    klein, plain = m["flux"]["klein"], m["flux"]["ordinary"]
    cells = {(c["s_t"], c["s_f"]): c for c in st["cells"]}
    rows = [r for c in st["cells"] for r in c["amplitudes"]]
    reversing = [r for c in st["cells"] if c["s_t"] == -1
                 for r in c["amplitudes"]]
    preserving = [r for c in st["cells"] if c["s_t"] == 1
                  for r in c["amplitudes"]]
    flux_cells = {(c["s_t"], c["s_f"]): c for c in adm["cells"]}
    rev = m["review"]
    mp, con, gau = rev["matter_phase"], rev["constrained"], rev["gauge"]
    half_sites = m["sites"] // 2

    def unit(values, centre, band, count):
        return len(values) == count and all(abs(v - centre) < band
                                            for v in values)

    return {
        "structure": {
            "E1 each deck map is an involution with the predicted "
            "symplectic sign": len(deck) == 4 and all(
                r["involution_error"] < ZERO
                and r["symplectic_sign_error"] < ZERO for r in deck),
            "E1 every deck map preserves H and maps the Gauss function "
            "covariantly": all(
                r["energy_invariance_error"] < EXACT
                and r["gauss_covariance_error"] < EXACT for r in deck),
            "E1 fixed space is half, Lagrangian iff the seam reverses time":
                all(2 * r["fixed_dim"] == r["phase_space_dim"]
                    and ((r["omega_on_fixed"] < ZERO) if r["s_t"] == -1
                         else (r["omega_rank_on_fixed"] == r["fixed_dim"]))
                    for r in deck),
            "E1 scheme is reversed by the deck map": (
                st["reversibility_error"] < ZERO),
            "E1 mirror data are deck fixed and on the Gauss surface": (
                set(cells) == set(ke.CELLS) and len(rows) == 16
                and all(r["initial_fixed_error"] < EXACT for r in rows)),
            "E1 responding histories satisfy the seam relation": (
                len(reversing) == 8
                and all(r["symmetry_defect"] < ZERO for r in reversing)),
            "E1 time-preserving fixed set is invariant under the flow": (
                len(preserving) == 8
                and all(r["symmetry_defect"] < ZERO for r in preserving)),
            "E1 tangent space is a Lagrangian half": (
                2 * st["tangent_dim"] == m["phase_space_dim"]
                and st["tangent_rank"] == st["tangent_dim"]
                and st["tangent_lagrangian_error"] < 1e-7),
            "E1 Gauss law holds and total charge vanishes": all(
                r["gauss_max"] < ZERO and r["total_charge"] < EXACT
                for r in rows),
            "E1 energy closes to 1e-5": all(
                r["energy_drift"] < 1e-5 for r in rows),
            "E1 the histories carried charge, field strength and mouth "
            "motion": all(r["local_charge"] > 1e-7
                          and r["max_field_strength"] > 1e-7
                          and r["max_mouth"] > 1e-7 for r in rows),
        },
        "constrained_structure": {
            "R2 frozen tangent directions leave the Gauss surface": (
                con["ambient_linearized_gauss_violation"] > 0.1
                and 2 * con["ambient_fixed_dim"] == con["phase_space_dim"]),
            "R2 the presymplectic kernel on the Gauss surface is the gauge "
            "directions": (
                con["gauss_rank"] == m["sites"]
                and con["surface_dim"] == con["phase_space_dim"] - m["sites"]
                and con["presymplectic_nullity"] == con["gauge_rank"]
                == m["sites"]
                and con["gauge_in_surface_error"] < ZERO
                and con["gauge_in_presymplectic_kernel"] < ZERO),
            "R2 constrained fixed space is isotropic": (
                con["constrained_omega"] < ZERO),
            "R2 it is half the reduced space after removing its gauge "
            "directions": (
                con["gauge_directions_in_fixed"] == half_sites
                and con["reduced_dim"] == con["phase_space_dim"]
                - 2 * m["sites"]
                and 2 * con["reduced_fixed_dim"] == con["reduced_dim"]),
            "R2 constrained tangents stay on the surface and isotropic "
            "along a responding history": (
                con["propagated_gauss_violation"] < 1e-7
                and con["propagated_omega"] < 1e-7
                and con["propagated_rank"] == con["constrained_fixed_dim"]),
        },
        "matter_phase_route": {
            "R1 the data are admissible with no charge, field, momentum or "
            "link potential": all(
                mp["initial"][k] < EXACT for k in (
                    "deck_fixed_error", "charge_density", "field", "momentum",
                    "link_potential")),
            "R1 the gauge-invariant input is a nonzero link current": (
                mp["initial"]["max_link_current"] > 1e-6),
            "R1 sheet charge is generated at second order in amplitude": (
                unit(mp["sheet_charge_slopes"], 2, 0.15, 3)
                and all(r["max_sheet_charge"] > 1e-9 for r in mp["rows"])),
            "R1 handle-loop flux and holonomy stay zero throughout": all(
                r["max_handle_flux"] < EXACT and r["max_loop_holonomy"] < EXACT
                and r["field_swap_odd_error"] < EXACT
                and r["max_field"] > 1e-7 for r in mp["rows"]),
            "R1 the histories are admissible and obey Gauss's law": all(
                r["seam_defect"] < ZERO and r["gauss_max"] < ZERO
                and r["gauss_identity_error"] < ZERO for r in mp["rows"]),
            "R1 charge conjugation reverses the sheet charge": (
                mp["conjugate_sum"] < EXACT),
            "R1 the sheet charge is already there with the coupling off": (
                mp["uncoupled_max_sheet_charge"] > 1e-6
                and mp["uncoupled_max_field"] < EXACT),
            "R1 an independent integrator and a halved step agree": (
                mp["dop853_disagreement"] < 1e-2
                and mp["step_refinement_change"] < 1e-2
                and mp["dop853_gauss_max"] < ZERO),
        },
        "gauge_dependence": {
            "R3 a deck-preserving gauge transformation changes the "
            "imaginary part": (
                gau["gauged_deck_fixed_error"] < EXACT
                and gau["imaginary_before"] < EXACT
                and gau["imaginary_after"] > 1e-3
                and gau["imaginary_history_difference"] > 1e-3),
            "R3 it leaves field, currents, charge, holonomy and energy "
            "unchanged": all(v < ZERO
                             for v in gau["invariant_differences"].values()),
        },
        "flux_admissibility": {
            "E2 cover has three independent cycles": adm["cycle_dim"] == 3,
            "E2 handle loop is source free and swap even": (
                adm["handle_cycle_divergence"] < EXACT
                and adm["handle_cycle_swap_error"] < EXACT),
            "E2 admissible flux dimension is 2 in matched cells, 1 in "
            "mixed": all(flux_cells[c]["admissible_flux_dim"]
                         == (2 if c in MATCHED else 1) for c in ke.CELLS),
            "E2 flux through the tubes is admissible iff the cell is "
            "matched": all(
                (flux_cells[c]["max_on_tube_edges"] > 1e-3) if c in MATCHED
                else (flux_cells[c]["max_on_tube_edges"] < EXACT)
                for c in ke.CELLS),
        },
        "spontaneous_pairing": {
            "E3 no imaginary part or gauge field from real data": all(
                len(m["spontaneous"][k]["rows"]) == 5
                and all(r["imaginary"] < EXACT and r["gauge"] < EXACT
                        for r in m["spontaneous"][k]["rows"])
                for k in ("klein", "ordinary")),
            "E3 the mouths were responding in those runs": all(
                r["max_mouth"] > 1e-7 for k in ("klein", "ordinary")
                for r in m["spontaneous"][k]["rows"]),
            "R3 gauge-invariantly: no current, charge or field from real "
            "data": (gau["quiet_current"] < EXACT
                     and gau["quiet_density"] < EXACT
                     and gau["quiet_field"] < EXACT),
        },
        "flux_pairing": {
            "E4 imaginary part is first order in the trapped flux": unit(
                klein["imaginary_slopes"], 1, 0.1, 2),
            "E4 it is exactly odd in the flux, the real part even": (
                klein["odd_error"] < ZERO and klein["even_error"] < ZERO
                and klein["imaginary_scale"] > 1e-6),
            "E4 the flux histories are admissible": (
                klein["initial_fixed_error"] < EXACT
                and klein["seam_relation_defect"] < ZERO),
            "R3 link current, sheet charge and loop holonomy are first "
            "order in the flux": (
                unit(gau["current_slopes"], 1, 0.1, 2)
                and unit(gau["sheet_charge_slopes"], 1, 0.1, 2)
                and unit(gau["holonomy_slopes"], 1, 0.1, 2)
                and gau["current_scale"] > 1e-7),
            "R3 they are exactly odd in the flux": (
                gau["current_odd_error"] < EXACT
                and gau["sheet_charge_odd_error"] < EXACT
                and gau["holonomy_odd_error"] < EXACT),
        },
        "flux_dynamics": {
            "E5 flux response is second order in field amplitude": unit(
                klein["flux_slopes_in_amplitude"], 2, 0.15, 3),
            "E5 flux response is first order in the flux": unit(
                klein["flux_slopes_in_flux"], 1, 0.1, 2),
        },
        "pairing_carrier": {
            "E6 sheet charge equals the flux through the seam": (
                klein["gauss_identity_error"] < ZERO
                and klein["gauss_max"] < ZERO),
            "E6 total charge vanishes and the sheet charge is odd in time":
                (klein["total_charge"] < EXACT
                 and klein["sheet_charge_antisymmetry_error"] < ZERO),
            "E6 charge separates across the seam": (
                klein["sheet_charge_max"] > 1e-9),
        },
        "seam_specificity": {
            "E7 ordinary handle: flux couples the sectors at first order":
                unit(plain["imaginary_slopes"], 1, 0.1, 2),
            "E7 ordinary handle: exactly odd in the flux": (
                plain["odd_error"] < ZERO and plain["even_error"] < ZERO
                and plain["imaginary_scale"] > 1e-6),
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
    checks = {check: all(bool(per_size[s][field][check]) for s in SIZES)
              for field, names in groups.items() for check in names}
    return {
        "public_preregistration": PREREG, "reviewed_commit": REVIEWED,
        "scope": SCOPE, "seed": ke.SEED,
        "withdrawn_claims": WITHDRAWN,
        "frozen_labels_changed": dict(ke.FROZEN_LABELS),
        "measurements": measured,
        "checks_by_size": {s: {c: bool(v) for g in per_size[s].values()
                               for c, v in g.items()} for s in SIZES},
        "checks": checks, "check_groups": groups,
        "checks_passed": all(checks.values()),
        "verdict": ke.verdict(checks, groups),
    }


def render(report):
    v = report["verdict"]
    prim = report["measurements"]["primary"]
    repl = report["measurements"]["replicate"]
    lines = ["# A responding electric field on the Klein handle", "",
             report["scope"], "",
             f"Public freeze: `{report['public_preregistration']}`. Seed "
             f"`{report['seed']}`. Phase-space dimensions "
             f"`{prim['phase_space_dim']}` and `{repl['phase_space_dim']}`; "
             "every required check must hold at both sizes.", "",
             "| verdict field | value |", "|---|---|"]
    lines += [f"| {k} | **{v[k]}** |" for k in ke.VERDICT_FIELDS]
    if v["failed_checks"] or v["missing_checks"]:
        lines += ["", "**Outside tolerance:** "
                  + ", ".join(f"`{c}`" for c in v["failed_checks"]),
                  "", "**Not run:** "
                  + ", ".join(f"`{c}`" for c in v["missing_checks"])]
    lines += ["", f"Corrected after review of `{report['reviewed_commit']}`. "
              "Withdrawn:", ""]
    lines += [f"- *{claim}.* {why}"
              for claim, why in report["withdrawn_claims"].items()]
    lines += ["", "## Review: sheet pairing without trapped flux", "",
              "| size | amplitude | max sheet charge | max handle flux | "
              "max field |", "|---|---:|---:|---:|---:|"]
    for size, m in (("primary", prim), ("replicate", repl)):
        for r in m["review"]["matter_phase"]["rows"]:
            lines.append(f"| {size} | {r['amplitude']:.6g} | "
                         f"{r['max_sheet_charge']:.3e} | "
                         f"{r['max_handle_flux']:.1e} | {r['max_field']:.3e} |")
    lines += ["", "## Review: the Gauss surface", "",
              "| size | ambient fixed | violation | constrained fixed | "
              "gauge in it | reduced dim | reduced fixed |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for size, m in (("primary", prim), ("replicate", repl)):
        c = m["review"]["constrained"]
        lines.append(
            f"| {size} | {c['ambient_fixed_dim']} | "
            f"{c['ambient_linearized_gauss_violation']:.4f} | "
            f"{c['constrained_fixed_dim']} | {c['gauge_directions_in_fixed']} "
            f"| {c['reduced_dim']} | {c['reduced_fixed_dim']} |")
    g = prim["review"]["gauge"]
    lines += ["", "## Review: gauge dependence (primary size)", "",
              f"Imaginary norm `{g['imaginary_before']:.1e}` before a "
              f"deck-preserving gauge transformation, "
              f"`{g['imaginary_after']:.4f}` after. Largest change in a "
              "gauge-invariant observable: "
              f"`{max(g['invariant_differences'].values()):.1e}`. "
              "Gauge-invariant slopes in the flux: current "
              f"`{[round(x, 4) for x in g['current_slopes']]}`, sheet charge "
              f"`{[round(x, 4) for x in g['sheet_charge_slopes']]}`, holonomy "
              f"`{[round(x, 4) for x in g['holonomy_slopes']]}`.",
              "", "## Admissible trapped flux (primary size)", "",
              "| s_t | s_f | admissible flux dim | max on tube edges |",
              "|---:|---:|---:|---:|"]
    lines += [f"| {c['s_t']:+d} | {c['s_f']:+d} | "
              f"{c['admissible_flux_dim']} | {c['max_on_tube_edges']:.2e} |"
              for c in prim["admissibility"]["cells"]]
    lines += ["", "## Trapped flux, real field", "",
              "| size | seam | Im slopes in flux | odd error | flux slopes "
              "in amplitude | max sheet charge |", "|---|---|---|---:|---|---:|"]
    for size, m in (("primary", prim), ("replicate", repl)):
        for name in ("klein", "ordinary"):
            f = m["flux"][name]
            lines.append(
                f"| {size} | {name} | "
                f"`{[round(s, 4) for s in f['imaginary_slopes']]}` | "
                f"{f['odd_error']:.1e} | "
                f"`{[round(s, 4) for s in f['flux_slopes_in_amplitude']]}` | "
                f"{f['sheet_charge_max']:.2e} |")
    lines += ["", "## Supplementary, not frozen", "",
              "| size | seam | screening at peak | max sheet charge | "
              "smallest per-field max sheet charge |",
              "|---|---|---:|---:|---:|"]
    for size, m in (("primary", prim), ("replicate", repl)):
        sup = m["supplementary"]
        for c in sup["cells"]:
            seam = "klein" if c["s_t"] == -1 else "ordinary"
            lines.append(
                f"| {size} | {seam} | {c['screening_at_peak']} of "
                f"{sup['draws']} | {c['max_sheet_charge']:.2e} | "
                f"{c['min_over_draws_of_max_sheet_charge']:.2e} |")
    lines += ["", "| Required check (both sizes) | Pass |", "|---|---|"]
    lines += [f"| {k} | {val} |" for k, val in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/"
              f"{len(report['checks'])} required checks.", "",
              "The orienting datum is conjugation-odd initial data, trapped "
              "flux or a matter current; neither is derived. An indefinite "
              "pairing is not a norm.",
              ""]
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
