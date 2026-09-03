"""Probe: positive coarea count or holonomy-weighted coarea current?

Against `docs/closure_current_prereg.md` (commit `f954e3d`), written before
`geometrodynamics/bulk/closure_current.py`.
Run:  python -m experiments.closure_ledger.closure_current_probe
"""

from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk import closure_current as cc   # noqa: E402


def run_probe() -> dict:
    r1 = cc.pin_loop_reduction()
    c1 = cc.branch_holonomy_is_sign_D()
    r2 = {g: cc.holonomy_weighted_law(g) for g in (0.3, 1.0, 2.0)}
    r3 = cc.sector_prior_control()
    r4 = cc.stationarity_audit()
    r2b = {g: cc.singlet_loop_law(g, n=40001) for g in (0.3, 1.0, 2.0)}
    r3b = cc.oriented_sector_prior_control()
    r4b = cc.integrated_sector_weights(n=40001)
    r4c = cc.oriented_current_audit()
    r5 = cc.pin_label_versus_weight()
    checks = {
        "R1_pin_loop_reduces_to_triangle_with_partner_sign": r1["reduces_to_triangle"],
        "C1_branch_holonomy_is_sign_D": c1["signed_density_is_holonomy_times_coarea"],
        "R2_holonomy_weighted_current_is_born_without_projectors": all(h["E_is_cos"] and h["max_deviation"] < 1e-9 for h in r2.values()),
        "R2b_derived_singlet_loop_gives_minus_cos_directly": all(
            h["E_is_minus_cos"] and h["max_deviation"] < 1e-9 for h in r2b.values()),
        "R3_sector_prior_chosen_marginals_fixed_E_moves": r3["marginals_stay_half"] and r3["E_moves"],
        "R3b_sector_prior_load_bearing_on_the_oriented_branch_too": r3b["quantum_law_only_at_ratio_one"],
        "R4_stationarity_unimplemented_and_proxy_is_disjoint_from_closure": (
            (not r4["implemented"]) and (not r4["proxy_tests_the_repository_condition"])
            and r4["phase_stationarity_proxy"]["closure_and_phase_stationarity_are_disjoint"]),
        "sector_integrals_of_the_oriented_current_are_nonnegative": r4b["all_sector_integrals_nonnegative"],
        "R5_pin_supplies_label_not_weight": r5["selected_by_pin_structure"] == "nothing about the measure",
    }
    v = cc.verdict(r1["reduces_to_triangle"], c1["signed_density_is_holonomy_times_coarea"],
                   checks["R2_holonomy_weighted_current_is_born_without_projectors"],
                   r4["stationarity_decides_the_fork"], False)
    return {"timestamp": datetime.now(timezone.utc).isoformat(),
            "prereg": "docs/closure_current_prereg.md @ f954e3d",
            "R1": r1, "C1": {k: val for k, val in c1.items() if k != "pi_branch_fraction_by_sector"},
            "R2": {g: {"E": h["E"], "max_deviation": h["max_deviation"]} for g, h in r2.items()},
            "R2b": {g: {"E": h["E"], "max_deviation": h["max_deviation"]} for g, h in r2b.items()},
            "R3": r3, "R3b": r3b, "underived_inputs": cc.underived_inputs(),
            "R4": {k: r4[k] for k in ("named_in_module_docstring", "repository_condition",
                                      "implemented", "weight_form", "why")},
            "R4_proxy": r4["phase_stationarity_proxy"],
            "sector_integrals": r4b, "oriented_current_audit": r4c,
            "R5": r5, "checks": checks, "passed": sum(checks.values()), "total": len(checks),
            "verdict": v,
            "dependency_ledger": {
                "loop": "loop( source pair [x, -x via J], geodesic realignment at A, B [chosen], J = L_{-j} [derived] ) -> triangle x -> u -> -v -> x [derived]",
                "closure holonomy": "-sgn D(x, u, -v) [derived from J^2 = -1 and the lift G]",
                "positive coarea": "|D|/(2|u x v|) [derived conditioning; counting measure on sectors: chosen]",
                "oriented current": "e^{i Phi} x coarea [derived label; adopting it as weight: the open step]",
                "sector prior": "counting [chosen]"},
            "open_audit_item": cc.oriented_current_audit(),
            "where_the_gap_is": ("Not the spin structure (derived), not the loop (derived), not Bell "
                                 "(evaded by measurement dependence either way), not the relative "
                                 "outcome law (the oriented current gives it analytically). It is "
                                 "whether observed event frequencies are the positive count of closed "
                                 "histories or their oriented, holonomy-weighted sum. Nothing classical "
                                 "in the repository decides; the rule forbids deciding it by the answer. "
                                 "Corrected after review: this is NOT a single binary choice -- "
                                 "the equal outcome-sector prior moves the correlation on BOTH "
                                 "branches (quantum only at r = 1) and the current-to-frequency "
                                 "readout (linear or quadratic in the integrated current) is a "
                                 "third open item. "
                                 "Sharper form, recorded as an open audit item and not a criterion: "
                                 "if the observable is an integral of a section of the local system "
                                 "the Pin/Hopf data define over the closure locus, the sign is "
                                 "geometrically mandatory; if the physical object is a measure on "
                                 "histories, positivity forces |D|.")}


def render(s: dict) -> str:
    L = [f"# Closure-current fork probe — {s['passed']}/{s['total']}", "",
         f"Pre-registration: `{s['prereg']}`. Verdict: **`{s['verdict']}`**", "",
         "| candidate | E(1) | S_max |", "|--|--|--|"]
    for k, v in s["R5"]["candidates"].items():
        L.append(f"| `{k}` | `{v['E(1)']:+.5f}` | `{v['S_max']:.4f}` |")
    L += ["", "R2b (the derived singlet loop, computed directly): "
          + ", ".join(f"γ={g}: `E = {h['E']:+.6f}`" for g, h in s["R2b"].items())
          + " — `P(s_A,s_B) = (1 − s_A s_B cos γ)/4`, deviations `≤ 1e-9`.",
          "R3b (the sector prior on the oriented branch): "
          + ", ".join(f"r={r['ratio']}: `E_triplet = {r['E_triplet']:+.5f}`, "
                      f"`E_singlet = {r['E_singlet']:+.5f}`" for r in s["R3b"]["rows"])
          + " — the quantum law only at `r = 1`, with marginals `1/2` throughout.",
          f"R1: `q5 = -q0 G^-1` to `{s['R1']['q5_equals_minus_q0_Ginv']:.1e}`; frame holonomy = Ω(x;u,−v) to "
          f"`{s['R1']['frame_holonomy_is_Omega_of_triangle_x_u_minus_v']:.1e}`; lift = cos(Ω/2)+sin(Ω/2)x to "
          f"`{s['R1']['lift_is_cos_half_Omega_plus_sin_half_Omega_x']:.1e}`.",
          f"C1: e^{{iΩ/2}} = sgn D on the closure set to `{s['C1']['exp_i_half_Omega_equals_sign_D']:.1e}`.",
          "R2: holonomy-weighted current gives P_like = (1+cos γ)/4, P_unlike = (1−cos γ)/4: deviations "
          + ", ".join(f"γ={g}: `{v['max_deviation']:.1e}`" for g, v in s["R2"].items()) + ".",
          "R3: sector prior ratio 0.5/1/2 → E(1) = " + ", ".join(f"`{r['E']:.4f}`" for r in s["R3"]["rows"]) + "; marginals 1/2 throughout.",
          f"R4: the repository's condition is *{s['R4']['repository_condition']}*, implemented: "
          f"{s['R4']['implemented']}. The phase-stationarity **proxy** is analytically disjoint from sharp "
          f"closure: `∇Ω|_{{N=0}} = 2(u×v)/D`, minimum norm `{s['R4_proxy']['min_gradient_norm_on_closure_set']:.4f}` "
          f"on the closure set (finite-difference residual `{s['R4_proxy']['finite_difference_residual']:.1e}`); the "
          "only points with `D = 0` are `x = −u, −v`, where the `arg` chart is singular, not stationary.",
          "Oriented-current sector integrals: `∫_Γ D_s dσ = 2π(1+u·v) ≥ 0` (residual "
          f"`{s['sector_integrals']['max_quadrature_residual']:.1e}`) — the cancellation is internal, so the "
          "holonomy weighting is wave interference, not negative event probabilities.", "", "## Checks", ""]
    for k, ok in s["checks"].items():
        L.append(f"* {'PASS' if ok else 'FAIL'} — {k}")
    L += ["", "## Where the gap is", "", s["where_the_gap_is"], "",
          "**Three inputs remain underived, not one binary choice:**", ""]
    for item in s["underived_inputs"]:
        L.append(f"1. {item}")
    L += ["", "## Dependency ledger", ""]
    for k, v in s["dependency_ledger"].items():
        L.append(f"* `{k}` = {v}")
    return "\n".join(L) + "\n"


def main() -> int:
    s = run_probe()
    text = render(s)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs", f"{stamp}_closure_current_probe")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "probe.json"), "w") as h:
        json.dump(s, h, indent=2, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else str(o))
    with open(os.path.join(outdir, "probe.md"), "w") as h:
        h.write(text)
    return 0 if s["passed"] == s["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
