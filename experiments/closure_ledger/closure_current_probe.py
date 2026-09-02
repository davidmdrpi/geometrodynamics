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
    r5 = cc.pin_label_versus_weight()
    checks = {
        "R1_pin_loop_reduces_to_triangle_with_partner_sign": r1["reduces_to_triangle"],
        "C1_branch_holonomy_is_sign_D": c1["signed_density_is_holonomy_times_coarea"],
        "R2_holonomy_weighted_current_is_born_without_projectors": all(h["E_is_cos"] and h["max_deviation"] < 1e-9 for h in r2.values()),
        "R3_sector_prior_chosen_marginals_fixed_E_moves": r3["marginals_stay_half"] and r3["E_moves"],
        "R4_stationarity_unimplemented_and_does_not_decide": (not r4["implemented"]) and r4["stationary_set"]["no_stationary_points"],
        "R5_pin_supplies_label_not_weight": r5["selected_by_pin_structure"] == "nothing about the measure",
    }
    v = cc.verdict(r1["reduces_to_triangle"], c1["signed_density_is_holonomy_times_coarea"],
                   checks["R2_holonomy_weighted_current_is_born_without_projectors"],
                   r4["stationarity_decides_the_fork"], False)
    return {"timestamp": datetime.now(timezone.utc).isoformat(),
            "prereg": "docs/closure_current_prereg.md @ f954e3d",
            "R1": r1, "C1": {k: val for k, val in c1.items() if k != "pi_branch_fraction_by_sector"},
            "R2": {g: {"E": h["E"], "max_deviation": h["max_deviation"]} for g, h in r2.items()},
            "R3": r3, "R4": {k: r4[k] for k in ("named_in_module_docstring", "implemented", "weight_form")},
            "R4_stationary_points": r4["stationary_set"]["n_stationary_points"],
            "R4_sin_distance_to_closure_circle": r4["stationary_set"]["sin_distance_to_closure_circle"],
            "R5": r5, "checks": checks, "passed": sum(checks.values()), "total": len(checks),
            "verdict": v,
            "dependency_ledger": {
                "loop": "loop( source pair [x, -x via J], geodesic realignment at A, B [chosen], J = L_{-j} [derived] ) -> triangle x -> u -> -v -> x [derived]",
                "closure holonomy": "-sgn D(x, u, -v) [derived from J^2 = -1 and the lift G]",
                "positive coarea": "|D|/(2|u x v|) [derived conditioning; counting measure on sectors: chosen]",
                "oriented current": "e^{i Phi} x coarea [derived label; adopting it as weight: the open step]",
                "sector prior": "counting [chosen]"},
            "where_the_gap_is": ("Not the spin structure (derived), not the loop (derived), not Bell "
                                 "(evaded by measurement dependence either way), not the relative "
                                 "outcome law (the oriented current gives it analytically). It is "
                                 "whether observed event frequencies are the positive count of closed "
                                 "histories or their oriented, holonomy-weighted sum. Nothing classical "
                                 "in the repository decides; the rule forbids deciding it by the answer.")}


def render(s: dict) -> str:
    L = [f"# Closure-current fork probe — {s['passed']}/{s['total']}", "",
         f"Pre-registration: `{s['prereg']}`. Verdict: **`{s['verdict']}`**", "",
         "| candidate | E(1) | S_max |", "|--|--|--|"]
    for k, v in s["R5"]["candidates"].items():
        L.append(f"| `{k}` | `{v['E(1)']:+.5f}` | `{v['S_max']:.4f}` |")
    L += ["", f"R1: `q5 = -q0 G^-1` to `{s['R1']['q5_equals_minus_q0_Ginv']:.1e}`; frame holonomy = Ω(x;u,−v) to "
          f"`{s['R1']['frame_holonomy_is_Omega_of_triangle_x_u_minus_v']:.1e}`; lift = cos(Ω/2)+sin(Ω/2)x to "
          f"`{s['R1']['lift_is_cos_half_Omega_plus_sin_half_Omega_x']:.1e}`.",
          f"C1: e^{{iΩ/2}} = sgn D on the closure set to `{s['C1']['exp_i_half_Omega_equals_sign_D']:.1e}`.",
          "R2: holonomy-weighted current gives P_like = (1+cos γ)/4, P_unlike = (1−cos γ)/4: deviations "
          + ", ".join(f"γ={g}: `{v['max_deviation']:.1e}`" for g, v in s["R2"].items()) + ".",
          "R3: sector prior ratio 0.5/1/2 → E(1) = " + ", ".join(f"`{r['E']:.4f}`" for r in s["R3"]["rows"]) + "; marginals 1/2 throughout.",
          f"R4: stationarity named in the docstring, implemented: {s['R4']['implemented']}; stationary points of Ω "
          f"in the source direction: {s['R4_stationary_points']} (Lexell: the level sets are circles through −u and −v; "
          "the closure phase has no critical points).", "", "## Checks", ""]
    for k, ok in s["checks"].items():
        L.append(f"* {'PASS' if ok else 'FAIL'} — {k}")
    L += ["", "## Where the gap is", "", s["where_the_gap_is"], "", "## Dependency ledger", ""]
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
