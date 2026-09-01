"""Probe: does negative-coupling EGB actually work?

Frozen in `docs/negative_egb_prereg.md` **before this module existed** (commit
`998e328`). It does not: the branch closes on physical grounds.

Run:  python -m experiments.closure_ledger.negative_egb_probe
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from typing import List

import numpy as np

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk.negative_egb import (  # noqa: E402
    measure_the_graviton_degenerates,
    measure_the_throat_matter_is_not_exotic,
    measure_no_coupling_satisfies_both,
    measure_the_bracket_is_continuous_at_the_seam,
    measure_the_critical_exterior_is_empty,
    measure_the_exterior_constrains_alpha_oppositely,
    measure_the_negative_egb_ledger,
)


def run_probe() -> dict:
    checks: List[dict] = []

    opposite = measure_the_exterior_constrains_alpha_oppositely()
    checks.append({
        "id": "E1", "name": "*** the exterior constrains alpha the OTHER way ***",
        "detail": opposite,
        "pass": bool(opposite["ricci_null_is_positive"]
                     and opposite["ratio_is_chi_independent"]
                     and opposite["the_two_bounds_coincide"])})

    scan = measure_no_coupling_satisfies_both()
    checks.append({
        "id": "E2", "name": "*** no open set of couplings survives (a search) ***",
        "detail": scan,
        "pass": bool(scan["the_surviving_set_is_measure_zero"]
                     and scan["survivors_are_at_the_critical_value"])})

    seam = measure_the_bracket_is_continuous_at_the_seam()
    checks.append({
        "id": "E3", "name": "and the two bounds meet by continuity, not by luck",
        "detail": seam,
        "pass": bool(seam["ratio_is_continuous"]
                     and seam["bracket_vanishes_at_the_seam"])})

    empty = measure_the_critical_exterior_is_empty()
    checks.append({
        "id": "E4",
        "name": "the surviving coupling forces a vacuum-form 5D exterior",
        "detail": empty,
        "pass": bool(empty["it_is_exactly_vacuum_energy"]
                     and empty["pressure_is_coupling_independent"])})

    honest = measure_the_throat_matter_is_not_exotic()
    checks.append({
        "id": "E4b",
        "name": "and there the throat matter is NOT exotic (NEC, WEC, DEC)",
        "detail": honest,
        "pass": bool(honest["nec_holds"] and honest["wec_holds"]
                     and honest["dec_holds"])})

    graviton = measure_the_graviton_degenerates()
    checks.append({
        "id": "E4c",
        "name": "*** THE CLOSURE: the graviton kinetic term vanishes there ***",
        "detail": graviton,
        "pass": bool(graviton["law_holds"]
                     and graviton["kinetic_vanishes_at_criticality"]
                     and graviton["superluminal_below_criticality"])})

    ledger = measure_the_negative_egb_ledger()
    checks.append({
        "id": "E5", "name": "the verdict and what is left untested",
        "detail": ledger, "pass": bool(len(ledger["entries"]) >= 5)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    opp, scan = detail("E1"), detail("E2")
    seam, empty = detail("E3"), detail("E4")
    honest, graviton, ledger = detail("E4b"), detail("E4c"), detail("E5")

    L: List[str] = [
        "# Does negative-coupling EGB actually work?", "",
        f"**{summary['passed']}/{summary['total']} checks pass — the branch "
        "closes, on the graviton rather than the matter.**", "",
        "Frozen in `docs/negative_egb_prereg.md` before this module existed.", "",
        "> ## The step the previous round missed", "",
        "> `α_GB` is a **coupling constant in the action**, so the same value "
        "acts in the exterior the throat is glued into. PR #277 analysed the "
        "throat in isolation, and it should not have. The NEC then pins one "
        "coupling — and at exactly that coupling the graviton degenerates.", "",
        f"> **{ledger['verdict']}.**", "",
        "---", "",
        "## E1 — the exterior constrains `α_GB` in the opposite direction", "",
        f"For `ℝ × S⁴_R`: `R_kk = {opp['ricci_null']:.6f} = 3/R²` — **positive**, "
        "because the exterior holds ordinary matter rather than flaring out — "
        f"and `H_kk = {opp['lanczos_null']:.6f} = 12/R⁴`. The ratio "
        f"`H_kk/R_kk = 4μ/f⁴ = {opp['ratio_value']:.6f} = 4/R²` is "
        f"**independent of `χ`** (`{opp['ratio_is_chi_independent']}`), so the "
        "exterior's constraint is a single number.", "",
        "| `α_GB` | `8πG₅T_kk` outside | exterior NEC? |", "|--|--|--|"]
    for row in opp["rows"]:
        L.append(f"| `{row['coupling']:+.6f}` | "
                 f"`{row['exterior_null_stress']:+.6f}` | "
                 f"**{'yes' if row['exterior_nec'] else 'no'}** |")
    L += ["",
          "```",
          f"exterior needs  α_GB ≥ {opp['exterior_needs_alpha_at_least']:.6f}",
          f"throat   needs  α_GB ≤ {opp['throat_needs_alpha_at_most']:.6f}",
          "```", "",
          "> " + opp["why"], "",
          "## E2 — searching for a coupling that works everywhere", "",
          f"A {scan['samples']}-point scan over `{scan['scanned_range']}`. A "
          "surviving *interval* would reopen the branch, so this looks for one "
          "rather than asserting none.", "",
          f"Survivors: **`{scan['couplings_satisfying_both']}`**, a set of width "
          f"`{scan['surviving_width']:.1e}`.", "",
          "> **" + scan["what_it_means"] + "**", "",
          "## E3 — the two bounds meet by continuity", "",
          "| side of the seam | `H_kk/R_kk` |", "|--|--|",
          f"| throat, at `s = S` | `{seam['ratio_inside_at_the_seam']:.9f}` |",
          f"| exterior, at `χ = a` | `{seam['ratio_outside_at_the_seam']:.9f}` |",
          "",
          f"Jump `{seam['ratio_jump']:.1e}`, and at the critical coupling the "
          f"shared bracket is `{seam['bracket_inside']:.1e}` — exactly zero.", "",
          "> " + seam["why"], "",
          "## E4 — what the one surviving coupling costs", "",
          "| `α_GB` | `8πG₅ρ` | `8πG₅p` | `ρ+p` | `w = p/ρ` |",
          "|--|--|--|--|--|"]
    for row in empty["rows"]:
        L.append(f"| `{row['coupling']:+.6f}` | `{row['density']:+.6f}` | "
                 f"`{row['pressure']:+.6f}` | `{row['sum']:+.6f}` | "
                 f"`{row['equation_of_state']:+.4f}` |")
    L += ["",
          "> " + empty["why_the_pressure_cannot_move"], "",
          "> **" + empty["what_this_costs"] + "**", "",
          "## E4b — but the throat matter there is not exotic", "",
          "> **" + honest["what_this_corrects"] + "**", "",
          f"`q` runs `{honest['q_at_the_mouth']:.4f}` (mouth) to "
          f"`{honest['q_at_the_neck']:.4f}` (neck). Minima along the throat: "
          f"`ρ+p_s = {honest['min_nec_radial']:+.1e}`, "
          f"`ρ+p_Ω = {honest['min_nec_angular']:+.4f}`, "
          f"`ρ−|p_s| = {honest['min_dec_radial']:+.1e}`. "
          f"NEC **{honest['nec_holds']}**, WEC **{honest['wec_holds']}**, "
          f"DEC **{honest['dec_holds']}**.", "",
          "## E4c — the closure: the graviton", "",
          "Linearising the full `G_ab + α H_ab` on **this** background — a "
          "product, not a maximally symmetric spacetime, so the textbook "
          "coefficient does not apply:", "",
          "| `α_GB` | `C_kin` (`ω²`) | predicted | `C_spatial` (`κ²`) | `c²` |",
          "|--|--|--|--|--|"]
    for row in graviton["rows"]:
        speed = ("∞" if not np.isfinite(row["speed_squared"])
                 else f"`{row['speed_squared']:.2f}`")
        L.append(f"| `{row['coupling']:+.5f}` | "
                 f"`{row['temporal_coefficient']:+.7f}` | "
                 f"`{row['predicted_kinetic']:+.7f}` | "
                 f"`{row['spatial_coefficient']:+.7f}` | {speed} |")
    L += ["",
          "```", graviton["kinetic_law"],
          "c^2 = 1/(1 + 4 alpha/R^2)", "```", "",
          "> " + graviton["why_it_had_to_be_derived"], "",
          "> **" + graviton["why_this_closes_the_branch"] + "**", "",
          "> **" + graviton["and_it_is_bad_before_criticality"] + "**", "",
          "## E5 — the ledger", "", "| claim | verdict | evidence |",
          "|--|--|--|"]
    for entry in ledger["entries"]:
        L.append(f"| {entry['claim']} | **{entry['verdict']}** | "
                 f"{entry['evidence']} |")
    L += ["", "**What this closes.** " + ledger["what_this_closes"], "",
          "**What remains untested.** " + ledger["what_remains_untested"], "",
          "**The remaining branches.**", ""]
    for key, value in ledger["the_remaining_branches"].items():
        L.append(f"- **{key}** — {value}")
    L += [""]
    return "\n".join(L)


def main() -> int:
    summary = run_probe()
    text = render_markdown(summary)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs",
                          f"{stamp}_negative_egb_probe")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "probe.json"), "w") as handle:
        json.dump(summary, handle, indent=2, default=float)
    with open(os.path.join(outdir, "probe.md"), "w") as handle:
        handle.write(text)
    print(f"\n\nWrote: {os.path.join(outdir, 'probe.json')}")
    print(f"Wrote: {os.path.join(outdir, 'probe.md')}")
    return 0 if summary["passed"] == summary["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
