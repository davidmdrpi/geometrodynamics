"""Probe: does a classical BAM action select the oriented branch?

Against `docs/history_action_prereg.md` (commit `a33a901`), written before
`geometrodynamics/bulk/history_action.py`.
Run:  python -m experiments.closure_ledger.history_action_probe
"""

from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk import history_action as ha   # noqa: E402


def run_probe() -> dict:
    oracle = ha.morse_bott_oracle()
    degen = ha.class_function_degeneracy()
    addit = ha.additive_functionals_have_no_critical_points()
    kappas = {k: ha.saddle_branch_ratio(v) for k, v in
              {"pi/4": math.pi / 4, "3pi/4": 3 * math.pi / 4,
               "1 (hbar=1)": 1.0, "pi/2": math.pi / 2}.items()}
    orbits = ha.sector_orbits()
    fibre = ha.fibre_action_is_weight_blind()
    homog = ha.detector_response_homogeneity()
    quad = ha.quadratic_readout_law()
    radial = ha.radial_action_compatibility()
    gate = ha.source_observable_signalling()
    v = ha.verdicts()

    checks = [
        ("O1-O4  the frozen oracle reproduces", oracle["passes"]),
        ("O4     saddle magnitude IS the positive coarea density",
         oracle["O4_saddle_magnitude_is_coarea_rel"] < 1e-3),
        ("A-ctl  every F(cos theta) shares the critical set",
         degen["all_share_the_critical_set"]),
        ("A1     theta additive, S_H not",
         addit["theta_additive"] and addit["S_H_not_additive"]),
        ("A2     the additive functional is nowhere stationary",
         addit["theta_has_no_critical_point"]),
        ("F1     saddle magnitudes never separate the branches",
         all(abs(r["magnitude"] - 1.0) < 1e-12 for r in kappas.values())),
        ("F2     ratio real iff 4kappa/pi odd",
         kappas["pi/4"]["ratio_is_real"] and kappas["3pi/4"]["ratio_is_real"]
         and not kappas["1 (hbar=1)"]["ratio_is_real"]),
        ("B1     two sector orbits at every angle but pi/2",
         orbits["forced_only_at_right_angle"]),
        ("B2     r not forced at any CHSH angle",
         not orbits["forced_at_any_chsh_angle"]),
        ("B3     fibre symmetries cannot change sector weights",
         fibre["fibre_blind"]),
        ("C1     every BAM coupling is degree 2", homog["all_quadratic"]),
        ("C2     a quadratic readout is superquantum",
         quad["quadratic_exceeds_tsirelson"]),
        ("D1     the radial action is a genuine one-form integral",
         radial["radial_action_is_a_one_form_integral"]
         and radial["both_legs_nonzero"]),
        ("E1     source density is exactly antipodally even",
         gate["density_is_antipodally_even_residual"] < 1e-15),
        ("E2     odd observables blind, even observables signal",
         gate["odd_observables_are_blind"] and gate["even_observables_signal"]),
        ("E3     non-coplanar supports are mutually singular",
         gate["non_coplanar_supports_mutually_singular"]),
    ]
    return {
        "oracle": oracle, "degeneracy": degen, "additivity": addit,
        "kappa": kappas, "orbits": orbits, "fibre": fibre,
        "homogeneity": homog, "quadratic": quad, "radial": radial,
        "gate": gate, "verdicts": v,
        "ledger": ha.dependency_ledger(),
        "checks": checks,
        "passed": sum(1 for _, ok in checks if ok), "total": len(checks),
    }


def render(s: dict) -> str:
    L = ["# Round 7 — does a classical BAM action select the oriented branch?", "",
         f"**{s['passed']}/{s['total']} checks pass.**", "",
         "## Five independent verdicts", ""]
    for k, val in s["verdicts"].items():
        L.append(f"* **{k}** — `{val}`")
    L += ["", "## Checks", ""]
    for name, ok in s["checks"]:
        L.append(f"* {'PASS' if ok else 'FAIL'}  {name}")

    o, a, q, g = s["oracle"], s["additivity"], s["quadratic"], s["gate"]
    L += ["", "## The candidate, and why it is not an action", "",
          f"* `-1/2 Tr G = -cos(Omega/2) = -D/sqrt(D^2+N^2)` to `{o['O1_closed_form']:.1e}`",
          f"* closure set is the critical manifold: `|grad S_H|` <= `{o['O2_grad_on_closure']:.1e}`",
          f"* transverse curvature `sgn(D)|uxv|^2/D^2`, rel `{o['O3_transverse_curvature_rel']:.1e}`;"
          f" index is `sgn D`: {o['O3_index_is_sign_D']}",
          f"* **saddle magnitude = positive coarea density**, rel "
          f"`{o['O4_saddle_magnitude_is_coarea_rel']:.1e}`",
          "",
          f"* `theta` additive to `{a['theta_is_additive_residual']:.1e}`;"
          f" `S_H` additivity defect `{a['S_H_additivity_defect']:.3f}`",
          f"* `min |grad theta|` on closure = `{a['min_|grad theta|_on_closure']:.4f}` > 0",
          ""]
    L += ["| kappa | 4k/pi | real? | selects |", "|---|---|---|---|"]
    for name, r in s["kappa"].items():
        L.append(f"| {name} | {r['4kappa_over_pi']:.4f} | {r['ratio_is_real']} "
                 f"| {r['selects']} |")

    L += ["", "## B — sector orbits", "",
          "| gamma | group order | orbits | r forced |", "|---|---|---|---|"]
    for r in s["orbits"]["rows"]:
        L.append(f"| {r['gamma']:.4f} | {r['group_order']} | {r['n_orbits']} "
                 f"| {r['r_forced']} |")

    L += ["", "## C — measured homogeneity of every BAM coupling", "",
          "| coupling | where | degree |", "|---|---|---|"]
    for r in s["homogeneity"]["rows"]:
        L.append(f"| `{r['coupling']}` | `{r['where']}` | {r['measured_degree']:.6f} |")
    L += ["", f"* linear readout: `S_max = {q['S_max_linear']:.6f}` (Tsirelson)",
          f"* quadratic readout: `S_max = {q['S_max_quadratic']:.6f}` = `8 sqrt2/3`"
          f" at `gamma = {q['argmax_gamma_quadratic']:.4f}` — **superquantum**",
          f"* quadratic marginals stay `1/2` to `{q['quadratic_marginal_deviation']:.1e}`"]

    L += ["", "## E — causality gate", "",
          f"* conditioned density antipodally even to `{g['density_is_antipodally_even_residual']:.1e}`",
          f"* odd observables blind (spread `{g['odd_observable_spread']:.1e}`)",
          f"* even observables signal (spread `{g['even_observable_spread']:.4f}`)",
          f"* non-coplanar total variation `{g['non_coplanar_total_variation']}`"
          " — mutually singular, a one-shot signal",
          "", "## Dependency ledger", ""]
    for row in s["ledger"]:
        L.append(f"* `{row['input']}` — **{row['status']}** ({row['where']})")
    return "\n".join(L) + "\n"


def main() -> int:
    s = run_probe()
    text = render(s)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs", f"{stamp}_history_action_probe")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "probe.json"), "w") as h:
        json.dump(s, h, indent=2, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else str(o))
    with open(os.path.join(outdir, "probe.md"), "w") as h:
        h.write(text)
    return 0 if s["passed"] == s["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
