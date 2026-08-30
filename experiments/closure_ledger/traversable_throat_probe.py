"""Probe: the whole-throat S-matrix of a supported traversable 5D throat.

Run:  python -m experiments.closure_ledger.traversable_throat_probe
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from typing import List

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.tangherlini.traversable_throat import (  # noqa: E402
    measure_the_closure_offsets_disagree,
    measure_the_geometry_is_derived_not_recalled,
    measure_the_null_energy_price,
    measure_the_scattering_is_symmetric_and_unitary,
    measure_the_threshold_law,
    measure_the_traversable_throat_ledger,
    measure_whether_the_loop_can_close,
)


def run_probe() -> dict:
    checks: List[dict] = []

    geometry = measure_the_geometry_is_derived_not_recalled()
    checks.append({
        "id": "G0", "name": "the geometry and its potential, derived",
        "detail": geometry,
        "pass": bool(geometry["throat_value_is_exact"]
                     and geometry["asymptotic_tail_matches_pr_275"]
                     and geometry["the_potential_is_even"]
                     and geometry["the_potential_is_positive"])})

    price = measure_the_null_energy_price()
    checks.append({
        "id": "G1", "name": "*** the GR price of traversability ***",
        "detail": price,
        "pass": bool(price["density_vanishes_identically"]
                     and price["radial_nec_is_violated_everywhere"]
                     and price["the_integral_is_exact"])})

    scattering = measure_the_scattering_is_symmetric_and_unitary()
    checks.append({
        "id": "G2", "name": "the whole-throat S-matrix is unitary and reciprocal",
        "detail": scattering,
        "pass": bool(scattering["unitarity_holds"]
                     and scattering["second_order_in_the_step"])})

    threshold = measure_the_threshold_law()
    checks.append({
        "id": "G3", "name": "the threshold law (NOT T_0(0) = g)",
        "detail": threshold,
        "pass": bool(threshold["the_magnitude_law_holds"]
                     and threshold["the_ratio_is_a_constant_phase"])})

    closure = measure_the_closure_offsets_disagree()
    checks.append({
        "id": "G4",
        "name": "*** phase closure and group closure demand different offsets ***",
        "detail": closure, "pass": bool(closure["the_two_offsets_disagree"])})

    loop = measure_whether_the_loop_can_close()
    checks.append({
        "id": "G5", "name": "*** Lambda = 1 cannot occur at finite frequency ***",
        "detail": loop,
        "pass": bool(loop["transmission_is_below_one_at_finite_frequency"]
                     and loop["the_approach_to_unity_is_exponential"])})

    ledger = measure_the_traversable_throat_ledger()
    checks.append({
        "id": "G6", "name": "the ledger", "detail": ledger,
        "pass": bool(len(ledger["entries"]) >= 8)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    price = detail("G1")
    closure = detail("G4")
    loop = detail("G5")

    L: List[str] = [
        "# The whole-throat S-matrix of a supported traversable 5D throat", "",
        f"**{summary['passed']}/{summary['total']} checks pass.**", "",
        "`transaction/network.py` (PR #216) replaced `handshake.py`'s "
        "phenomenological advanced confirmation with a Morris–Thorne–Yurtsever "
        "mechanism — everything propagates forward in local time, and the "
        "mouths' clock offset carries the return into the exterior past. That "
        "was a real relocation of the postulate. But the geometry was never "
        "supplied: `MouthPort.t`, `r_in`, `r_out` are inputs, and "
        "`closure_offset` *solves* for the `Δ` that makes the loop close.", "",
        "This round supplies the geometry and recomputes the closure without "
        "tuning `Δ` to the answer.", "",
        "> ## Three results", "",
        "> **1.** Traversability costs exactly "
        f"`{price['null_energy_integral_in_units_of_one_over_G5']:.6f}/G₅` "
        "= `−3/(16G₅a)` in the complete null-ray energy integral.", "",
        "> **2.** No single clock offset closes the loop: phase closure and "
        f"group closure disagree by up to "
        f"`{closure['worst_disagreement']:.2f}` over the sampled band.", "",
        "> **3.** `Λ = 1` — PR #216's completed transaction — cannot occur at "
        "any finite frequency, and the deficit closes only exponentially, "
        f"`1 − |T|² ~ exp({loop['log_slope']:.2f} ω)`.", "",
        "---", "",
        "## Why not Tangherlini — an argument, not a computation", "",
        "`supp G_ret(c,s) ⊂ J⁺(s)` and `supp G_adv(c,d) ⊂ J⁻(d)`, so a nonzero "
        "product `G_ret(c,s)G_adv(c,d)` needs `c ∈ J⁺(s) ∩ J⁻(d)`. But then "
        "`s → c → d` is a causal chain, so `d ∈ J⁺(s)`. Contrapositive: if "
        "`d ∉ J⁺(s)` the product vanishes for **every** `c`. **An advanced leg "
        "by itself does not evade ER non-traversability.** PR #275's `T_ℓ` is "
        "exterior→*horizon*, and is exactly that zero for the cross-mouth "
        "channel.", "",
    ]

    geometry = detail("G0")
    L += ["## G0 — the geometry, derived", "",
          f"`{geometry['profile']}`, giving `{geometry['potential']}`.", "",
          "| `ℓ` | `V_ℓ(0)` | exact `[ℓ(ℓ+2)+3/2]/a²` | tail coeff | exact `(ℓ+1)²−¼` |",
          "|--|--|--|--|--|"]
    for row in geometry["rows"]:
        L.append(f"| {row['ell']} | `{row['V_at_throat']:.6f}` | "
                 f"`{row['exact_throat_value']:.6f}` | "
                 f"`{row['asymptotic_coefficient']:.6f}` | "
                 f"`{row['exact_asymptotic_coefficient']:.6f}` |")
    L += ["", "**" + geometry["it_is_one_barrier_not_two"] + "**", "",
          geometry["why_the_same_jost_condition_applies"], "",
          "*" + geometry["the_throat_constant_echoes_pr_275"] + "*", ""]

    L += ["## G1 — the GR price", "",
          "| `s` | `8πG₅ρ` | `8πG₅p_s` | `8πG₅p_Ω` |", "|--|--|--|--|"]
    for s, d, ps, po in zip(price["sample_s"], price["density"],
                            price["radial_pressure"], price["angular_pressure"]):
        L.append(f"| `{s:g}` | `{d:.6f}` | `{ps:+.6f}` | `{po:+.6f}` |")
    L += ["",
          "`ρ = 0` identically, and the radial NEC `ρ + p_s < 0` **everywhere**. "
          "Along a complete radial null geodesic with `k^t̂ = 1`:", "",
          "```",
          f"∫ T_ab k^a k^b dλ = {price['null_energy_integral_in_units_of_one_over_G5']:.6f}/G₅ "
          f"= {price['exact_null_energy_integral']}   (exact)",
          "```", "",
          "The exoticity scales as `1/(G₅a)` — a wider throat is cheaper.", "",
          "> **" + price["what_this_is_NOT"] + "**", ""]

    scattering = detail("G2")
    L += ["## G2 — the whole-throat S-matrix", "",
          "| `ω` | `\\|T\\|` | `\\|R\\|` | `\\|R\\|²+\\|T\\|²−1` |",
          "|--|--|--|--|"]
    for w, t, r, u in zip(scattering["omega"], scattering["transmission_modulus"],
                          scattering["reflection_modulus"],
                          scattering["unitarity_residual"]):
        L.append(f"| `{w:g}` | `{t:.8f}` | `{r:.8f}` | `{u:+.1e}` |")
    L += ["",
          f"Worst unitarity residual `{scattering['worst_unitarity_residual']:.1e}`, "
          "imposed nowhere. " + scattering["reciprocity_is_structural"], ""]

    threshold = detail("G3")
    L += ["## G3 — the threshold law, and the false regression avoided", "",
          "> " + threshold["the_static_conductance_is_not_T_at_zero"], "",
          f"Static: `I₃ = {threshold['static_resistance_I3']:.4f}`, "
          f"`g = {threshold['static_conductance_g']:.6f} = π²a²`. The dynamical "
          "law is instead:", "",
          "```", threshold["law"], "```", "",
          "| `ω` | `\\|T₀\\|/(π/8)(aω)²` | `arg(ratio)/π` |", "|--|--|--|"]
    for w, m, p in zip(threshold["omega"], threshold["magnitude_ratio"],
                       threshold["phase_of_ratio_over_pi"]):
        L.append(f"| `{w:g}` | `{m:.6f}` | `{p:+.4f}` |")
    L += ["", "> " + threshold["the_phase_offset_is_a_convention"], ""]

    L += ["## G4 — no single clock offset closes the loop", "",
          "PR #216 sets `Δ_BA = −(d_A + d_B + τ_th)`, which *solves* for "
          "closure. Once the throat is dispersive there is no unique `τ_th`: "
          "the geometry supplies `δ_ℓ(ω) = arg T_ℓ(ω)`, and phase closure "
          "`Φ = 2πn` and group closure `dΦ/dω = 0` demand **different** "
          "offsets.", "",
          f"With `d_A + d_B = {closure['exterior_legs']:.4f}` (antipodal on the "
          "unit `S³`):", "",
          "| `ω` | `δ_ℓ` | `dδ/dω` (Wigner) | `Δ_phase` | `Δ_group` | gap |",
          "|--|--|--|--|--|--|"]
    for row in closure["rows"]:
        L.append(f"| `{row['omega']:g}` | `{row['delta']:+.4f}` | "
                 f"`{row['wigner_delay']:+.4f}` | `{row['delta_phase']:+.4f}` | "
                 f"`{row['delta_group']:+.4f}` | **`{row['disagreement']:.4f}`** |")
    L += ["",
          "> " + closure["why_that_matters"], "",
          closure["bandwidth_from_dispersion"], "",
          "**" + closure["what_this_replaces"] + "**", ""]

    L += ["## G5 — can the transaction actually complete?", "",
          "| `ω` | `1 − \\|T\\|²` |", "|--|--|"]
    for w, d in zip(loop["omega"][::4], loop["reflection_probability"][::4]):
        L.append(f"| `{w:.2f}` | `{d:.4e}` |")
    L += ["",
          "> **" + loop["so_lambda_cannot_equal_one_exactly"] + "**", "",
          "> " + loop["but_the_deficit_closes_exponentially"], "",
          "| `D_loop` | regime |", "|--|--|"]
    for key, value in loop["causal_classification"].items():
        L.append(f"| `{key}` | {value} |")
    L += [""]

    ledger = detail("G6")
    L += ["## G6 — the ledger", "", "| claim | verdict | evidence |", "|--|--|--|"]
    for e in ledger["entries"]:
        L.append(f"| {e['claim']} | **{e['verdict']}** | {e['evidence']} |")
    L += ["", "**What this round establishes.** "
          + ledger["what_this_round_establishes"], "",
          "**What remains postulated.** " + ledger["what_remains_postulated"], "",
          "**Still open.** " + ledger["still_open"], ""]
    return "\n".join(L)


def main() -> int:
    summary = run_probe()
    text = render_markdown(summary)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs",
                          f"{stamp}_traversable_throat_probe")
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
