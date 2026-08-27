"""The local geometry of the quark fit manifold: response Jacobian and SVD.

THE RESULT, FIRST
─────────────────
PR #272 measured one knob at a time, found the quark residuals individually
right and jointly wrong, and named the gap a "missing correlation" -- guessing
it lived between transport and resistance.

THE GUESS WAS WRONG, AND THE OBJECT IS NOT A SCALAR RELATION. It is the
response map

    J_ia = d ln(m_i / m_d) / d ln p_a ,   i in {s, c, b, t}

Three results, in order of what they cost:

1. THE FIT MANIFOLD IS FOUR-DIMENSIONAL AND THE MODEL HAS ELEVEN KNOBS.
   rank J = 4 exactly -- it cannot exceed the number of independently scored
   masses. Eight knobs carry first-order response, so FOUR of those eight
   directions move no observable at all. The masses do not constrain pinhole,
   transport, resistance and N separately; they constrain at most four
   combinations of everything. Every compensation seen since PR #76 is this.

2. THE CORRECTION MOVES AWAY FROM WHAT THE MASSES WANT.
       legacy     cos(Theta) = +0.464  ( 62.4 deg)
       corrected  cos(Theta) = -0.616  (128.0 deg)
   PR #272's three per-knob improvements were MISLEADING in the strict sense:
   every residual moved toward its lock while the displacement they jointly
   produce turned from partly-helpful to actively-harmful. A scalar residual
   cannot see this; only the projection can.

3. THE 0.018% MIN-NORM FIT IS LOCAL FLEXIBILITY, NOT STRUCTURE.
   It is carried 98.3% by the WEAKEST singular direction (sigma_4 = 0.852
   against sigma_1 = 19.2) -- the definition of a compensator. And it does not
   survive leave-one-species-out: fit {s, c, t}, predict b, and the answer is
   -10.4% with the fitted three at 0.06%.

THE NULL STRUCTURE, DERIVED RATHER THAN OBSERVED
────────────────────────────────────────────────
Three directions show zero first derivative and they are DIFFERENT OBJECTS.

  action_base       EXACT INVARIANCE, ALL ORDERS. H(a) = H(0) + a*I, and the
                    min_eigenvalue spectrum-zero subtraction removes it
                    identically -- UPSTREAM of the d-anchor, so it does not
                    reappear in an absolute-scale observable either. Dropped
                    from the first-order parameter space.

  phase             Z2-EVEN, QUADRATIC. Enters only via cos(phase*dk) -- but
                    ONLY because partition_mixing = 0 at the lock. The
                    different-partition element carries e^{i*phase*k}, which is
                    not even. THE Z2 IS A PROPERTY OF THE LOCK, NOT THE MODEL.

  partition_mixing  Z2-EVEN, QUADRATIC, BY UNITARY EQUIVALENCE.
                    H(-p) = D H(p) D† with D = diag(sigma(p)), exact to 0.0,
                    hence isospectral. Stronger than "even function": a
                    discrete gauge symmetry of the spectrum.

Both quadratic directions are handled by the chart q = x^2, under which
d ln m / dq is finite and constant. Their first derivatives in x are not
informative and are never reported as such.

NOT A SLOPPY MODEL
──────────────────
Condition number over the identifiable subspace is 22.6, not 10^3-10^6. The
degeneracy is DIMENSIONAL -- eleven knobs, four observables -- not
ill-conditioning. The fix is more observables, not fewer knobs: the v4
flavor-CP layer already yields CKM angles and J from the same Hamiltonian.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from geometrodynamics.qcd.response_jacobian import (  # noqa: E402
    LOG_KNOBS, SCORED_SPECIES, measure_leave_one_species_out,
    measure_the_geometric_displacement_against_what_the_masses_want,
    measure_the_jacobian_converges, measure_the_jacobian_ledger,
    measure_the_null_structure_is_three_different_objects,
    measure_the_singular_system_and_effective_rank,
    measure_which_directions_repair_the_masses, response_jacobian)


def run_probe() -> dict:
    checks: List[dict] = []

    nul = measure_the_null_structure_is_three_different_objects()
    checks.append({
        "id": "R0",
        "name": "the three apparent nulls are three DIFFERENT objects",
        "detail": nul,
        "pass": bool(nul["action_base"]["H_minus_identity_shift"] < 1e-10
                     and nul["action_base"]["zero_subtracted_spectrum_gap"] < 1e-10
                     and nul["phase"]["H_even_at_the_lock"] == 0.0
                     and nul["phase"]["H_even_with_mixing_switched_on"] > 1e-3
                     and nul["partition_mixing"]["conjugation_gap"] == 0.0
                     and all(v["relative_spread_over_two_decades"] < 1e-3
                             for v in nul["quadratic_chart"].values()))})

    cnv = measure_the_jacobian_converges()
    checks.append({
        "id": "R1", "name": "the Jacobian converges; no adiabatic relabelling noise",
        "detail": cnv, "pass": bool(cnv["all_converged_below_1e-3_relative"])})

    svd = measure_the_singular_system_and_effective_rank()
    checks.append({
        "id": "R2",
        "name": "*** rank 4 against 8 first-order knobs: 4 invisible directions ***",
        "detail": svd,
        "pass": bool(svd["rank"] == 4
                     and svd["rank_is_capped_by_the_observable_count"]
                     and svd["n_invisible_first_order_directions"] == 4
                     and svd["not_a_sloppy_model"])})

    rep = measure_which_directions_repair_the_masses()
    checks.append({
        "id": "R3", "name": "the repair rides the WEAKEST singular direction",
        "detail": rep,
        "pass": bool(rep["the_repair_rides_the_weakest_direction"]
                     and rep["dominant_share"] > 0.9
                     and rep["linear_residual_after_repair"] < 1e-12)})

    prj = measure_the_geometric_displacement_against_what_the_masses_want()
    checks.append({
        "id": "R4",
        "name": "*** the corrected geometry moves AWAY from what the masses want ***",
        "detail": prj,
        "pass": bool(prj["the_correction_flips_the_sign_of_the_overlap"])})

    hld = measure_leave_one_species_out()
    checks.append({
        "id": "R5", "name": "leave-one-out: local flexibility, not structure",
        "detail": hld,
        "pass": bool(hld["it_explodes_under_holdout"]
                     and all(r["fitted_species_max_error_percent_exact"] < 0.1
                             for r in hld["rows"]))})

    led = measure_the_jacobian_ledger()
    checks.append({"id": "R6", "name": "the ledger", "detail": led,
                   "pass": bool(len(led["entries"]) >= 9)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    L: List[str] = [
        "# The local geometry of the quark fit manifold",
        "",
        f"**{summary['passed']}/{summary['total']} checks pass.**",
        "",
        "PR #272 measured one knob at a time and named the gap between "
        "individually-right and jointly-wrong residuals a *missing "
        "correlation*, guessing it lived between `transport` and `resistance`. "
        "**That guess was wrong, and the object is not a scalar relation.** It "
        "is the response map `J_ia = ∂ln(m_i/m_d)/∂ln p_a` and its SVD.",
        "",
    ]

    nul = detail("R0")
    L += ["## R0 — three apparent nulls, three different objects", "",
          "| direction | status | mechanism | check |", "|--|--|--|--|",
          f"| `action_base` | **{nul['action_base']['status']}** | "
          f"`H(a) = H(0) + a·I`, killed by the spectrum-zero subtraction | "
          f"`{nul['action_base']['H_minus_identity_shift']:.1e}` |",
          f"| `phase` | **{nul['phase']['status']}** | enters only as "
          f"`cos(φ·dk)` — *because* `partition_mixing = 0` | "
          f"`{nul['phase']['H_even_at_the_lock']:.1f}` |",
          f"| `partition_mixing` | **{nul['partition_mixing']['status']}** | "
          f"`H(−p) = D H(p) D†`, `D = diag(σ)` — isospectral | "
          f"`{nul['partition_mixing']['conjugation_gap']:.1f}` |", "",
          "`action_base` is removed **upstream of the `d`-anchor** — the "
          "zero-point subtraction does it, and the alternative "
          "`spectrum_zero_mode=\"action_base\"` kills it too. So it does *not* "
          "reappear in an absolute-scale observable of this model. It is "
          "dropped from the first-order parameter space.", "",
          f"The `phase` Z₂ is a property of the **lock**, not the model: "
          f"switch `partition_mixing` on and the evenness breaks "
          f"(`{nul['phase']['H_even_with_mixing_switched_on']:.3f}`).", "",
          "Both quadratic directions are handled by the chart `q = x²`, under "
          "which `∂ln m/∂q` is finite and constant:", "",
          "| knob | `∂ln m/∂q` | spread over two decades in `q` |", "|--|--|--|"]
    for knob, v in nul["quadratic_chart"].items():
        L.append(f"| `{knob}` | `{np.linalg.norm(v['d_ln_m_d_q']):.4f}` | "
                 f"`{v['relative_spread_over_two_decades']:.1e}` |")
    L += ["", f"*{nul['why_they_must_not_share_a_category']}.*", ""]

    cnv = detail("R1")
    L += ["## R1 — the Jacobian converges", "",
          "`extract_physical_spectrum` reassigns species by adiabatic "
          "continuation, so a relabelling inside the difference stencil would "
          "show up as step-dependent noise. None appears — every column norm "
          "is stable across three decades of step size.", "",
          "| knob | column norm | relative spread, `1e-3`…`1e-6` |",
          "|--|--|--|"]
    for row in cnv["rows"]:
        norm = row["column_norm_by_step"]["1e-04"]
        L.append(f"| `{row['knob']}` | `{norm:.4f}` | "
                 f"`{row['relative_spread_over_1e-3_to_1e-6']:.1e}` |")
    L += ["",
          "The largest entry in the whole matrix is "
          "`∂ln m_b/∂ln uplift_asymmetry = −18.96`: **the shell-index axiom "
          "`ε = 1 − 1/k₅² = 24/25` is by far the most load-bearing number in "
          "the quark ladder**, and it is the one PR #272 classified as *exactly "
          "invariant* — correctly, since it reads no radial operator. Being "
          "safe from the operator correction is not the same as being "
          "unimportant, and its elasticity had never been measured.", ""]

    svd = detail("R2")
    L += ["## R2 — effective rank: four observables, eleven knobs", "",
          f"`rank J = {svd['rank']}` against `{svd['n_first_order_knobs']}` "
          f"first-order knobs, so **{svd['n_invisible_first_order_directions']} "
          f"directions move no observable at all** — before any accidental "
          f"degeneracy.", "",
          "| `a` | `σ_a` | `σ_a/σ_1` |", "|--|--|--|"]
    for a, (s, ratio) in enumerate(zip(svd["singular_values"],
                                       svd["sigma_ratios"])):
        L.append(f"| {a+1} | `{s:.4f}` | `{ratio:.4f}` |")
    L += ["",
          f"Condition number `{svd['condition_number']:.1f}`. **This is not a "
          f"sloppy model** in the Sethna sense — a sloppy spectrum spans "
          f"10³–10⁶. The identifiable subspace is well conditioned; the "
          f"degeneracy is *dimensional*: {svd['the_problem_is_dimensional']}.",
          "", "**Right singular vectors** (parameter combinations):", "",
          "| knob | " + " | ".join(f"`v{a+1}`" for a in range(4)) + " |",
          "|--" * 5 + "|"]
    for knob in LOG_KNOBS:
        L.append(f"| `{knob}` | " + " | ".join(
            f"`{svd['right_singular_vectors'][f'v{a+1}'][knob]:+.4f}`"
            for a in range(4)) + " |")
    L += ["", "**Left singular vectors** (species response):", "",
          "| species | " + " | ".join(f"`u{a+1}`" for a in range(4)) + " |",
          "|--" * 5 + "|"]
    for s in SCORED_SPECIES:
        L.append(f"| `{s}` | " + " | ".join(
            f"`{svd['left_singular_vectors'][f'u{a+1}'][s]:+.4f}`"
            for a in range(4)) + " |")
    L += ["",
          "The quadratic directions add no new observable direction — with "
          "rank 4 and four observables the column space is already all of "
          "`ℝ⁴`, so *any* 4-vector projects exactly. The content is **which** "
          "combination each mimics:", "",
          "| knob | dominant equivalent first-order knob |", "|--|--|"]
    for knob, v in svd["quadratic_directions_are_not_new_observable_directions"].items():
        L.append(f"| `{knob}` | `{v['dominant']}` |")
    L.append("")

    rep = detail("R3")
    L += ["## R3 — which singular directions actually repair the masses", "",
          "Reporting `σ_a` alone says nothing. `c_a = u_aᵀr` and `c_a/σ_a` say "
          "everything: a tiny `σ` is harmless if the residual has no "
          "projection on it, and a moderate one can dominate.", "",
          "| `a` | `σ_a` | `c_a = u_aᵀr` | `c_a/σ_a` | share of `\\|δx\\|²` |",
          "|--|--|--|--|--|"]
    for row in rep["rows"]:
        L.append(f"| {row['a']} | `{row['sigma']:.4f}` | `{row['c_a']:+.6f}` | "
                 f"`{row['c_over_sigma']:+.6f}` | "
                 f"`{row['share_of_correction_norm_squared']:.4f}` |")
    L += ["",
          f"**{100*rep['dominant_share']:.1f}% of the repair rides `v"
          f"{rep['dominant_direction']}` — the *weakest* direction.** "
          f"{rep['verdict']}.", "",
          f"The correction itself has `|δ ln p| = {rep['delta_x_min_norm']:.6f}` "
          f"and drives the linear residual to "
          f"`{rep['linear_residual_after_repair']:.1e}` (exact re-run: "
          f"`0.018%`, from the lock's `1.61%`). **That is not a success of the "
          f"physics.** It shows the frozen Hamiltonian has enough local freedom "
          f"to absorb essentially the whole residual — the old `1.61%` floor "
          f"was set by *where the lock sits*, not by the functional form.", ""]

    prj = detail("R4")
    L += ["## R4 — does the geometry move where the data want?", "",
          "The decisive test. Project the operator-induced displacement "
          "`δx_geom` onto the mass-optimal direction `δx_min = J⁺r`.", "",
          "| operator | `\\|δx_geom\\|` | `cos Θ` (parameter) | `cos Θ` "
          "(observable) | invisible fraction |", "|--|--|--|--|--|"]
    for name, label in (("legacy", "legacy"), ("scalar_correct", "corrected")):
        p = prj["per_operator"][name]
        L.append(f"| {label} | `{p['delta_x_geom_norm']:.6f}` | "
                 f"`{p['cos_theta_parameter_space']:+.4f}` | "
                 f"**`{p['cos_theta_observable_space']:+.4f}`** "
                 f"(`{p['angle_observable_space_deg']:.1f}°`) | "
                 f"`{p['fraction_of_displacement_invisible_to_the_masses']:.4f}` |")
    L += ["",
          f"> **{prj['verdict']}**", "",
          f"And {prj['and_most_of_it_is_invisible_anyway']} — so even the part "
          f"that points somewhere mostly points nowhere observable.", ""]

    hld = detail("R5")
    L += ["## R5 — leave-one-species-out", "",
          "Fit the minimum-norm correction on three masses; evaluate the "
          "fourth with **nothing readjusted**.", "",
          "| held out | `\\|δx\\|` | fitted three (exact) | **held-out (exact)** |",
          "|--|--|--|--|"]
    for row in hld["rows"]:
        L.append(f"| `{row['held_out']}` | `{row['correction_norm']:.6f}` | "
                 f"`{row['fitted_species_max_error_percent_exact']:.3f}%` | "
                 f"**`{row['held_out_error_percent_exact']:+.3f}%`** |")
    L += ["", f"{hld['verdict']}. The fitted species land at "
          f"`0.002–0.06%` every time; the withheld one does not come along.", ""]

    led = detail("R6")
    L += ["## R6 — the ledger", "", "| claim | verdict | evidence |",
          "|--|--|--|"]
    for e in led["entries"]:
        L.append(f"| {e['claim']} | **{e['verdict']}** | {e['evidence']} |")
    L += ["", f"**Headline.** {led['headline']}.", "",
          f"**What would settle it.** {led['what_would_settle_it']}.",
          ""]
    return "\n".join(L)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    if isinstance(o, tuple):
        return list(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = (Path(__file__).resolve().parent / "runs"
           / f"{stamp}_quark_response_jacobian_probe")
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2,
                                               default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0 if summary["passed"] == summary["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
