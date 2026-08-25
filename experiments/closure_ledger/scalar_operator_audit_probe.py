"""Correcting the radial scalar operator, and auditing every claim it touches.

THE CORRECTION
──────────────
Until PR #271, `tangherlini.radial.V_tangherlini` carried

    V_legacy = A[ l(l+2)/r^2 + 3 r_h^2/r^4 ] ,      A = 1 - r_h^2/r^2

which is NOT the master potential of a minimally coupled massless scalar. For
ds^2 = -A dt^2 + A^-1 dr^2 + r^2 dOmega_n^2 such a scalar obeys

    (1/r^n) d_r(r^n A R') + (w^2/A - l(l+n-1)/r^2) R = 0 ,

and the unique first-derivative-free Schroedinger form, with dr* = dr/A and
psi = r^(n/2) R, carries

    V_scalar = A[ l(l+n-1)/r^2 + n(n-2)A/(4r^2) + n A'/(2r) ] .

At n = 3 that is A[(l(l+2) + 3/4)/r^2 + (9/4) r_h^2/r^4], so

    V_scalar - V_legacy = 3 A^2 / (4 r^2) .

THIS IS A BUG, NOT A CONVENTION. The old generic name implied the canonical
scalar operator and the implementation was short of it by an l-independent
term. PR #270 discovered, independently validated and ISOLATED it without
touching anything; this round corrects it and audits what moved.

T0  THE CORRECTION, WITH THREE INDEPENDENT CONFIRMATIONS. The gap equals
    3A^2/(4r^2) in closed form (2.4e-15); it carries no l (2.4e-15); and the
    flat limit reproduces the Bessel form (2.2e-16), which is what settles
    WHICH operator is the scalar one -- psi = r^(1/2) J_(l+1)(wr) obeys
    Bessel's equation with V = ((l+1)^2 - 1/4)/r^2 = (l(l+2) + 3/4)/r^2. The
    general-n form was verified symbolically at n = 2..6, so it is standard
    and independently re-derived here rather than cited. `V_scalar_tangherlini`
    agrees BITWISE with `dynamics.master_potential`.

T1  THE EIGENVALUES BARELY MOVE, AND MOVE LESS AS l RISES.
        l     legacy          corrected       shift
        0     1.00065891      1.00198000      +0.1320%
        1     1.05472694      1.05582653      +0.1043%
        2     1.13156946      1.13239953      +0.0734%
        3     1.21908274      1.21966785      +0.0480%
        4     1.30869618      1.30909388      +0.0304%
        5     1.39597349      1.39624205      +0.0192%
    Eigenfunction overlaps exceed 0.999998 everywhere. The monotone fall with
    l is not a coincidence: an eigenvalue AVERAGES the potential against a
    bound state, so an l-independent shift matters least where the centrifugal
    term already dominates. This is the reassuring half of the audit.

T2  *** THE BARRIER SUMS DO NOT, AND THAT IS WHERE THE MEANING MOVES. *** A
    barrier height READS the potential directly instead of averaging it:
        l = 1..5:  22.00824 (-2.19%)  ->  22.33119 (-0.75%)
        l = 0..5:  22.45268 (-0.21%)  ->  22.83642 (+1.50%)
    So the two gamma statements in the tree move in OPPOSITE directions. The
    canonical README claim -- "Pinhole gamma = Sum V_max[1..5] ... -2.2% off
    the locked gamma = 22.5" -- IMPROVES threefold, to -0.75%, with nothing
    tuned. The l = 0..5 statement BREAKS: it overshoots at +1.50%, and the sum
    closest to 22.5 SWAPS from l = 0..5 to l = 1..5.

    The geometric root moves with it:
        R_OUTER at gamma = 22.5, l=0..5:  1.26227 -> 1.24614
        R_OUTER at gamma = 22.5, l=1..5:  1.28737 -> 1.26788

T3  WHAT SURVIVES EXACTLY. Delta V carries no l, so the cross-l perturbation
    operator V_(l+2) - V_l is unchanged to 3.6e-15 -- ALGEBRAICALLY EXACT. Its
    MATRIX ELEMENTS are not, because they are taken between eigenfunctions that
    drift. Structure invariant, numbers shifted; that distinction is the whole
    partition of the audit.

T4  RATIOS ABSORB MOST OF IT. alpha_q(l,0) is a ratio of throat derivatives
    normalised to l = 1, so the common part of the shift divides out and only
    the differential survives. Measured, not assumed -- a ratio is not
    automatically safe.

T5  ACTIONS MOVE MORE THAN EIGENVALUES. Closed-orbit radial actions integrate
    sqrt(w^2 - V) against the potential directly, without the bound state's
    averaging, so they sit between the eigenvalues and the barrier sums. Each
    action is evaluated at ITS OWN operator's ground frequency, so the
    comparison is between two self-consistent orbits.

T6  *** THE LEDGER. *** Every load-bearing radial claim sorted into three
    verdicts -- EXACTLY INVARIANT, NUMERICALLY SHIFTED, INTERPRETATION CHANGED
    -- because the question is not whether the old tests stay green but which
    published claims are algebraically untouched, which keep their meaning with
    different digits, and which no longer say what they said.

    The one INTERPRETATION CHANGED entry that matters: "adding the l = 0 5D
    channel closes the gamma discrepancy" is WITHDRAWN, not replaced. The
    corrected l = 1..5 sum does land nearer 22.5 than the old one did, but that
    is an observation, not a derivation of why 22.5.

WHAT IS DELIBERATELY NOT RE-RUN
───────────────────────────────
The Hopf fibration, the Pin- structure on the exchange RP^2, the odd-k winding
ladder and antipodal parity have NO dependence on the radial scalar operator.
Proximity is not dependence, and re-running them would manufacture the
appearance of an audit without adding one.

    python -m experiments.closure_ledger.scalar_operator_audit_probe
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.tangherlini.operator_audit import (
    LOCKED_GAMMA,
    measure_which_geometry_preserves_the_lepton_ladder,
    measure_the_downstream_ledger,
    measure_the_eigenvalue_shifts,
    measure_the_eigenvector_derived_quantities,
    measure_the_gamma_sums_and_the_r_outer_fixed_point,
    measure_the_two_operators_and_their_exact_gap,
    measure_the_wkb_action_shift,
    measure_what_survives_exactly,
)


def run_probe() -> dict:
    checks: List[dict] = []

    gap = measure_the_two_operators_and_their_exact_gap()
    checks.append({
        "id": "T0", "name": "the correction, with three independent confirmations",
        "detail": gap,
        "pass": bool(gap["gap_matches_the_closed_form"] < 1e-12
                     and gap["the_gap_carries_no_ell"] < 1e-12
                     and gap["flat_limit_matches_bessel"] < 1e-12
                     and gap["agrees_bitwise_with_dynamics_master_potential"])})

    eig = measure_the_eigenvalue_shifts()
    checks.append({
        "id": "T1", "name": "the eigenvalues barely move, and less as l rises",
        "detail": eig,
        "pass": bool(eig["all_shifts_below_a_fifth_of_a_percent"]
                     and eig["sensitivity_falls_with_ell"]
                     and eig["eigenfunctions_barely_move"])})

    gam = measure_the_gamma_sums_and_the_r_outer_fixed_point()
    checks.append({
        "id": "T2", "name": "*** the barrier sums move, and the gamma story swaps ***",
        "detail": gam,
        "pass": bool(gam["the_canonical_readme_claim_improves"]
                     and gam["the_ell_zero_closure_claim_breaks"]
                     and gam["the_closest_channel_set_swaps"])})

    inv = measure_what_survives_exactly()
    checks.append({
        "id": "T3", "name": "the cross-l operator survives exactly; its elements drift",
        "detail": inv,
        "pass": bool(inv["the_cross_ell_operator_is_unchanged"] < 1e-12
                     and inv["structure_invariant_numbers_shifted"])})

    rat = measure_the_eigenvector_derived_quantities()
    checks.append({
        "id": "T4", "name": "ratios absorb most of the shift",
        "detail": rat,
        "pass": bool(rat["the_reference_mode_is_exactly_one"]
                     and rat["ratios_absorb_most_of_the_shift"])})

    act = measure_the_wkb_action_shift()
    checks.append({
        "id": "T5", "name": "actions move more than eigenvalues",
        "detail": act,
        "pass": bool(act["each_action_uses_its_own_ground_frequency"]
                     and act["largest_drift_percent"]
                     > eig["largest_ground_shift_percent"])})

    led = measure_the_downstream_ledger()
    checks.append({
        "id": "T6", "name": "*** the ledger: invariant / shifted / reinterpreted ***",
        "detail": led,
        "pass": bool(led["counts"]["EXACTLY INVARIANT"] >= 2
                     and led["counts"]["NUMERICALLY SHIFTED"] >= 5
                     and led["counts"]["INTERPRETATION CHANGED"] >= 1)})

    lad = measure_which_geometry_preserves_the_lepton_ladder()
    checks.append({
        "id": "T7",
        "name": "*** the one narrow downstream re-derivation: gamma is the selector ***",
        "detail": lad,
        "pass": bool(lad["nothing_was_retuned"]
                     and lad["B_and_C_are_bit_identical"]
                     and lad["gamma_is_the_selector_r_outer_is_downstream"]
                     and lad["fixing_r_outer_breaks_the_ladder"])})

    return {
        "probe": "scalar_operator_audit",
        "question": "the radial potential was short of the minimally coupled "
                    "scalar master operator by 3A^2/(4r^2) -- which published "
                    "claims are algebraically untouched, which keep their "
                    "meaning with different digits, and which no longer say "
                    "what they said?",
        "answer": "the eigenvalues move at the 1e-3 level and less as l rises, "
                  "and the cross-l operator is exactly invariant; but the "
                  "barrier sums move enough that the two gamma statements in "
                  "the tree swap places -- the canonical l = 1..5 claim "
                  "improves threefold to -0.75% with nothing tuned, while the "
                  "claim that the l = 0 channel closes the gap overshoots to "
                  "+1.50% and is withdrawn",
        "headline": {
            "the_gap": gap["the_gap"],
            "omega_1_0": [eig["omega_1_0_legacy"], eig["omega_1_0_correct"]],
            "largest_eigenvalue_shift_percent": eig["largest_ground_shift_percent"],
            "canonical_gamma_residual": [gam["canonical_residual_before"],
                                         gam["canonical_residual_after"]],
            "ell_zero_gamma_residual": [gam["ell_zero_residual_before"],
                                        gam["ell_zero_residual_after"]],
            "closest_channel_set_swaps": gam["the_closest_channel_set_swaps"],
            "cross_ell_operator_unchanged_to": inv["the_cross_ell_operator_is_unchanged"],
            "verdict_counts": led["counts"],
            "B_and_C_bit_identical": lad["B_and_C_are_bit_identical"],
            "d_ln_m_mu_over_d_ln_gamma": lad["d_ln_m_mu_over_d_ln_gamma"],
        },
        "checks": checks,
        "passed": sum(1 for c in checks if c["pass"]),
        "total": len(checks),
    }


def render_markdown(summary: dict) -> str:
    L = ["# Scalar operator audit — correcting it, and pricing the correction", "",
         f"**Question.** {summary['question']}", "",
         f"**Answer.** {summary['answer']}", "",
         f"**{summary['passed']}/{summary['total']} checks pass.**", "",
         "| id | check | result |", "|----|-------|--------|"]
    for c in summary["checks"]:
        L.append(f"| {c['id']} | {c['name']} | {'PASS' if c['pass'] else 'FAIL'} |")

    g = next(c for c in summary["checks"] if c["id"] == "T0")["detail"]
    L += ["", "## T0 — the correction", "",
          "| | operator |", "|--|--|",
          f"| legacy | `{g['legacy']}` |",
          f"| corrected | `{g['scalar_correct']}` |",
          f"| gap | `{g['the_gap']}` |", "",
          f"Three independent confirmations: the gap matches that closed form to "
          f"`{g['gap_matches_the_closed_form']:.1e}`; it carries no `ℓ` to "
          f"`{g['the_gap_carries_no_ell']:.1e}`; and the flat limit reproduces "
          f"the Bessel form to `{g['flat_limit_matches_bessel']:.1e}`, which is "
          f"what settles *which* operator is the scalar one.", "",
          f"`V_scalar_tangherlini` agrees bitwise with "
          f"`dynamics.master_potential`: "
          f"`{g['agrees_bitwise_with_dynamics_master_potential']}`.", "",
          f"**{g['why_it_is_a_bug_not_a_convention'].capitalize()}.**", ""]

    e = next(c for c in summary["checks"] if c["id"] == "T1")["detail"]
    L += ["## T1 — the eigenvalues barely move", "",
          "| `ℓ` | legacy `ω_{ℓ,0}` | corrected | shift | min overlap |",
          "|--|--|--|--|--|"]
    for r in e["rows"]:
        L.append(f"| {r['ell']} | `{r['omega_legacy'][0]:.8f}` | "
                 f"`{r['omega_correct'][0]:.8f}` | "
                 f"`{r['ground_shift_percent']:+.4f}%` | "
                 f"`{r['min_eigenfunction_overlap']:.6f}` |")
    L += ["", f"The monotone fall with `ℓ` is not a coincidence: "
              f"{e['why_so_small']}.", ""]

    gm = next(c for c in summary["checks"] if c["id"] == "T2")["detail"]
    L += ["## T2 — the barrier sums, and where the meaning moves", "",
          "| channels | legacy sum | residual | corrected sum | residual | "
          "`R_OUTER` legacy | corrected |", "|--|--|--|--|--|--|--|"]
    for r in gm["rows"]:
        L.append(f"| `{r['channels']}` | `{r['sum_legacy']:.5f}` | "
                 f"`{r['residual_legacy_percent']:+.2f}%` | "
                 f"`{r['sum_correct']:.5f}` | "
                 f"`{r['residual_correct_percent']:+.2f}%` | "
                 f"`{r['r_outer_legacy']:.5f}` | `{r['r_outer_correct']:.5f}` |")
    L += ["",
          f"The two `γ` statements move in **opposite directions**. The "
          f"canonical README claim improves from "
          f"`{gm['canonical_residual_before']:+.2f}%` to "
          f"`{gm['canonical_residual_after']:+.2f}%` — and "
          f"{gm['and_it_is_not_a_refit']}. The `ℓ = 0…5` claim goes from "
          f"`{gm['ell_zero_residual_before']:+.2f}%` to "
          f"`{gm['ell_zero_residual_after']:+.2f}%`, and the sum closest to "
          f"`22.5` swaps from **{gm['closest_to_gamma_before']}** to "
          f"**{gm['closest_to_gamma_after']}**.", "",
          f"**What has to be reopened:** {gm['what_has_to_be_reopened']}.", ""]

    inv = next(c for c in summary["checks"] if c["id"] == "T3")["detail"]
    L += ["## T3 — what survives exactly", "",
          f"`ΔV` carries no `ℓ`, so `V_{{ℓ+2}} − V_ℓ` is unchanged to "
          f"`{inv['the_cross_ell_operator_is_unchanged']:.1e}` — algebraically "
          f"exact. Its matrix elements are not:", "",
          "| `ℓ` | element legacy | corrected | drift |", "|--|--|--|--|"]
    for r in inv["matrix_elements"]:
        L.append(f"| {r['ell']} | `{r['element_legacy']:.6e}` | "
                 f"`{r['element_correct']:.6e}` | `{r['drift_percent']:+.3f}%` |")
    L += ["", f"**{inv['the_partition'].capitalize()}.**", ""]

    rat = next(c for c in summary["checks"] if c["id"] == "T4")["detail"]
    L += ["## T4 — ratios absorb most of the shift", "",
          "| `ℓ` | `α_q` legacy | corrected | drift |", "|--|--|--|--|"]
    for r in rat["rows"]:
        L.append(f"| {r['ell']} | `{r['alpha_q_legacy']:+.6f}` | "
                 f"`{r['alpha_q_correct']:+.6f}` | `{r['drift_percent']:+.3f}%` |")
    L += ["", f"*{rat['caveat'].capitalize()}.*", ""]

    act = next(c for c in summary["checks"] if c["id"] == "T5")["detail"]
    L += ["## T5 — actions move more than eigenvalues", "",
          "| `ℓ` | action legacy | corrected | drift |", "|--|--|--|--|"]
    for r in act["rows"]:
        L.append(f"| {r['ell']} | `{r['action_legacy']:.6f}` | "
                 f"`{r['action_correct']:.6f}` | `{r['drift_percent']:+.3f}%` |")
    L += ["", f"{act['why_it_moves_more_than_omega'].capitalize()}.", ""]

    led = next(c for c in summary["checks"] if c["id"] == "T6")["detail"]
    L += ["## T6 — the ledger", "",
          f"{led['the_question_this_answers'].capitalize()}.", "",
          "| claim | verdict | evidence |", "|--|--|--|"]
    for x in led["entries"]:
        L.append(f"| {x['claim']} | **{x['verdict']}** | {x['evidence']} |")
    L += ["", "| verdict | count |", "|--|--|"]
    for k, v in led["counts"].items():
        L.append(f"| {k} | {v} |")
    L += ["", f"**Not re-run, and why.** {led['not_re_run_and_why'].capitalize()}.",
          "", f"**Still open.** {led['what_is_still_open'].capitalize()}.", ""]
    lad = next(c for c in summary["checks"] if c["id"] == "T7")["detail"]
    L += ["## T7 — the one narrow downstream re-derivation", "",
          f"{lad['the_three_geometries']}, each passed through the **locked** "
          f"lepton Hamiltonian with **nothing retuned**.", "",
          "| case | `R_OUTER` | `γ` | `m_μ` err | `m_τ` err |",
          "|--|--|--|--|--|"]
    for r in lad["rows"]:
        ro = "—" if r["r_outer"] is None else f"`{r['r_outer']:.5f}`"
        L.append(f"| {r['case']} | {ro} | `{r['gamma']:.5f}` | "
                 f"`{r['mu_error_percent']:+.2f}%` | "
                 f"`{r['tau_error_percent']:+.2f}%` |")
    L += ["",
          f"**B and C are bit-identical.** {lad['why'].capitalize()}. So the "
          f"channel-set choice leaves no trace in any observable once `γ` is "
          f"enforced — and the comparison cannot decide it.", "",
          f"> **`γ = 22.5` is the selector; `R_OUTER` is downstream of it.**", "",
          f"Fixing `R_OUTER` and letting `γ` float is what breaks the ladder: "
          f"`{lad['corrected_fixed_R_mu_errors'][0]:+.2f}%` and "
          f"`{lad['corrected_fixed_R_mu_errors'][1]:+.2f}%`, against the legacy "
          f"geometry's `{lad['legacy_fixed_R_mu_error']:+.2f}%`. So the "
          f"correction **weakens** the geometry-supplies-`γ` story even while "
          f"improving the `1..5` residual in isolation.", "",
          f"The reason is sensitivity: `d ln m_μ / d ln γ = "
          f"{lad['d_ln_m_mu_over_d_ln_gamma']:.2f}`, so a sub-percent geometric "
          f"residual is **not** a small residual in this chain.", "",
          f"**What this does not settle:** {lad['what_this_does_not_settle']}.", ""]
    return "\n".join(L)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = (Path(__file__).resolve().parent / "runs"
           / f"{stamp}_scalar_operator_audit_probe")
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2,
                                               default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0 if summary["passed"] == summary["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
