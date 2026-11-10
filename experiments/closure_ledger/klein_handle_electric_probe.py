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
SCOPE = ("Lattice scalar electrodynamics on the cover network of a ring with "
         "one handle: charged field, link potentials, electric fields and two "
         "waist coordinates solved from one Hamiltonian. Effective oscillator "
         "support. Trapped flux enters as initial data. Classical "
         "throughout; no probability rule, norm, spinor matching or "
         "statement about the seam BAM has.")
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
    }


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
        "public_preregistration": PREREG, "scope": SCOPE, "seed": ke.SEED,
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
    lines += ["", "## Admissible trapped flux (primary size)", "",
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
              "Trapped flux is initial data here. An orientation is not a "
              "complex structure, and an indefinite pairing is not a norm.",
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
