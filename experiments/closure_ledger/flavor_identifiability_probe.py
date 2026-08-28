"""Is the v4 CKM realization a prediction, or a locally flexible fit?

THE ANSWER
──────────
rank K = 4. THE CKM REALIZATION IS A FIT, NOT A PREDICTION.

The decisive object is not another pseudoinverse. Build two response maps over
one common parameter chart x,

    J_M = dy_M/dx ,   y_M = ln(m_{s,c,b,t} / m_d)
    J_F = dy_F/dx ,   y_F = (theta12, theta23, theta13, delta)

take the mass-preserving tangent space N_M = ker J_M, and form

    K = J_F N_M .

rank K counts the physically independent CKM directions reachable WITHOUT
disturbing the masses. If it equals 4 -- the full dimension of the physical
flavor space of a unitary 3x3 matrix -- the mass-preserving freedom already
spans everything the CKM can be, and fitting it predicts nothing.

It equals 4, robustly: stable over finite-difference steps 1e-5..1e-7 and
null-space cutoffs 1e-6..1e-10 (nine cells), singular values spread over only
379x, kernel verified to 5e-15. Confirmed by direct construction: an arbitrary
target dy_F is realized to ~1e-14 with the masses held to ~1e-14.

There is therefore NO left-null vector, hence NO first-order invariant relation
w^T dy_F = 0 to compare against experiment. The round was built to extract one
if it existed.

THE phi_h A/B TEST
──────────────────
    phi_h RELEASED   rank K = 4   sv = [2.611, 0.468, 0.00827, 0.00689]
    phi_h FIXED      rank K = 4   sv = [0.542, 0.366, 0.00793, 0.00357]

Fixing the derived Hopf holonomy does NOT lower the rank. The other fitted
matrix elements absorb arbitrary CKM data on their own, so the derived phase
produces no flavor prediction by itself. It is the most EFFICIENT CP handle --
releasing it multiplies the leading singular value by 4.8 -- but efficiency is
not identifiability.

This confirms and sharpens PR #173, which found that ADDING phi_h left the
observable rank unchanged. The stronger statement: with phi_h held FIXED, the
mass-preserving freedom still covers all four physical flavor coordinates.

THE CENSUS, MEASURED RATHER THAN COUNTED
────────────────────────────────────────
"New symbol count" and "calibration degree-of-freedom count" are different
numbers. Measuring the second as rank J_F restricted to each group:

    group                          symbols   measured flavor rank
    v4 targeted etas                  3              3
    eta_k3k5_minus retune             1              1
    diag_shift_plus                   3              2
    diag_shift_minus                  3              2
    phi_h                             1              1
    v4 additions, phi_h fixed         9              4

Two structural facts fall out. The TRACE direction of each diagonal-shift
triple is an EXACT CKM GAUGE (|J_F.1| ~ 2e-10 against |J_M.1| = 12.5), so a
uniform shift within a block moves masses and cannot have been selected on
flavor data. And both realized triples are TRACELESS to ~1e-10, with
diag_shift_plus collapsing further to (+d, -d, 0).

THE "+3 PARAMETERS FOR +5 INDEPENDENT OBSERVABLES" CLAIM
────────────────────────────────────────────────────────
It cannot hold, for a reason needing no computation: a unitary 3x3 CKM has
exactly FOUR physical parameters. The nine quoted flavor-CP observables are all
functions of those four, so at most four are independent. "+5 independent
observables" exceeds the ceiling. Against that ceiling the measured calibration
dimension is 4, so the net predictive surplus is at most ZERO, not +2.

WHAT THIS DOES NOT SAY
──────────────────────
The v4 numbers are not wrong and the <=1% agreement across nine observables is
real. What the rank shows is that the agreement is not evidence FOR the
Hamiltonian, because the same Hamiltonian could have reproduced any LOCALLY
NEIGHBOURING CKM equally well at the same masses.

RESTRICTED TO THE ACTUAL v4 CALIBRATION, TOO. The headline rank lets every
coordinate move, including the frozen v3 lock. The decisive surplus question is
K_v4 = J_F[:,G] . ker(J_M[:,G]) with the v3 knobs held fixed -- and over the 10
quantities selected in the re-lock (three targeted etas, the eta35 RETUNE, six
diagonal shifts) with phi_h also fixed, rank K_v4 = 4. Including the retune
improves the conditioning 8.0e-6 -> 3.1e-3.

AND IT IS NONLINEARLY TRUE. The tangent-space construction is algebra on the
same matrices; the independent check scales a direction by eps and re-runs the
Hamiltonian, giving clean O(eps^2) (ratios 4.00).

Scope: a LOCAL, first-order statement at the v4 lock. A rank says nothing about
how far the mass-preserving surface extends before nonlinearity or positivity
bites, and the parameter excursions required are not small (|dx| ~ 50-80 for a
unit dy_F) -- so the claim licenses INFINITESIMAL CKM directions, not arbitrary
finite CKM matrices.
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

from geometrodynamics.qcd.flavor_identifiability import (  # noqa: E402
    FLAVOR_OBSERVABLES, PHYSICAL_FLAVOR_DIMENSION,
    measure_both_jacobians_converge, measure_the_calibration_dof_census,
    measure_the_counting_claim, measure_the_flavor_ledger,
    measure_the_mass_preserving_flavor_rank,
    measure_the_observables_are_four_and_rephasing_invariant,
    measure_the_phi_h_ab_test,
    measure_the_reachability_is_nonlinearly_true,
    measure_the_restricted_v4_calibration_rank)


def run_probe() -> dict:
    checks: List[dict] = []

    obs = measure_the_observables_are_four_and_rephasing_invariant()
    checks.append({
        "id": "F0",
        "name": "four independent, rephasing-invariant flavor coordinates",
        "detail": obs,
        "pass": bool(len(obs["observables"]) == PHYSICAL_FLAVOR_DIMENSION
                     and obs["ckm_unitarity_defect"] < 1e-12
                     and obs["max_drift_under_random_rephasing"] < 1e-12)})

    cnv = measure_both_jacobians_converge()
    checks.append({
        "id": "F1", "name": "both Jacobians converge and their ranks are stable",
        "detail": cnv,
        "pass": bool(cnv["both_converged"] and cnv["ranks_stable"])})

    cen = measure_the_calibration_dof_census()
    checks.append({
        "id": "F2",
        "name": "calibration DOF measured, not counted; the trace is a CKM gauge",
        "detail": cen,
        "pass": bool(cen["the_trace_direction_is_flavor_blind"]
                     and cen["both_realised_triples_are_traceless"])})

    rnk = measure_the_mass_preserving_flavor_rank()
    checks.append({
        "id": "F3",
        "name": "*** rank K = 4: the CKM is a fit, not a prediction ***",
        "detail": rnk,
        "pass": bool(rnk["spans_the_whole_physical_flavor_space"]
                     and rnk["rank_stable_across_the_grid"]
                     and rnk["kernel_is_a_kernel"] < 1e-10
                     and rnk["left_null_vectors"] == 0
                     and all(t["flavor_miss"] < 1e-10
                             and t["mass_disturbance"] < 1e-10
                             for t in rnk["tangent_space_construction"]))})

    restricted = measure_the_restricted_v4_calibration_rank()
    checks.append({
        "id": "F3b",
        "name": "*** rank K_v4 = 4 with the FROZEN v3 lock and phi_h fixed ***",
        "detail": restricted,
        "pass": bool(restricted["every_group_reaches_the_full_flavor_space"]
                     and restricted["rank_K_v4"] == PHYSICAL_FLAVOR_DIMENSION
                     and restricted["including_the_retune_improves_conditioning"])})

    nonlinear = measure_the_reachability_is_nonlinearly_true()
    checks.append({
        "id": "F3c",
        "name": "the reachability is nonlinearly true: O(eps^2) on a real re-run",
        "detail": nonlinear,
        "pass": bool(nonlinear["both_are_second_order"])})

    ab = measure_the_phi_h_ab_test()
    checks.append({
        "id": "F4",
        "name": "*** fixing the derived phi_h does NOT lower the flavor rank ***",
        "detail": ab,
        "pass": bool(not ab["fixing_phi_h_lowers_the_rank"]
                     and ab["per_case"]["phi_h_fixed"]["rank_K"]
                     == PHYSICAL_FLAVOR_DIMENSION)})

    cnt = measure_the_counting_claim()
    checks.append({
        "id": "F5", "name": "the '+3 for +5, net +2' counting claim",
        "detail": cnt,
        "pass": bool(cnt["the_observable_claim_exceeds_the_ceiling"]
                     and cnt["measured_net_surplus_phi_h_fixed"] <= 0)})

    led = measure_the_flavor_ledger()
    checks.append({"id": "F6", "name": "the ledger", "detail": led,
                   "pass": bool(len(led["entries"]) >= 6)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    L: List[str] = [
        "# Is the v4 CKM realization a prediction, or a locally flexible fit?",
        "", f"**{summary['passed']}/{summary['total']} checks pass.**", "",
        "The decisive object is not another pseudoinverse. Build two response "
        "maps over one common parameter chart `x` — `J_M = ∂y_M/∂x` for the "
        "mass ratios and `J_F = ∂y_F/∂x` for `(θ₁₂, θ₂₃, θ₁₃, δ)` — take the "
        "mass-preserving tangent space `N_M = ker J_M`, and form **`K = J_F "
        "N_M`**. Its rank counts the physically independent CKM directions "
        "reachable *without disturbing the masses*.", "",
    ]

    obs = detail("F0")
    L += ["## F0 — four coordinates, not ten redundant ones", "",
          "| observable | rad | deg |", "|--|--|--|"]
    for name, r, d in zip(obs["observables"], obs["values_rad"],
                          obs["values_deg"]):
        L.append(f"| `{name}` | `{r:.6f}` | `{d:.3f}°` |")
    L += ["",
          f"Jarlskog `J = {obs['jarlskog']:.6e}`, CKM unitary to "
          f"`{obs['ckm_unitarity_defect']:.1e}`, and the extraction is "
          f"invariant under random rephasing of both eigenbases to "
          f"`{obs['max_drift_under_random_rephasing']:.1e}` — it reads only "
          f"`|V_ij|` and `J`, never LAPACK's phase convention.", "",
          f"> {obs['why_this_matters']}.", "",
          f"*{obs['validation_outputs_not_rows']}.*", ""]

    cen = detail("F2")
    L += ["## F2 — the calibration-DOF census, measured", "",
          "The measured calibration dimension of a knob group is `rank J_F` "
          "restricted to it — how many independent physical flavor directions "
          "it can actually move. **Symbol count and DOF count are different "
          "numbers.**", "",
          "| group | symbols | measured flavor rank |", "|--|--|--|"]
    for row in cen["rows"]:
        L.append(f"| {row['group']} | {row['symbols']} | "
                 f"**{row['measured_flavor_rank']}** |")
    L += ["", "**The trace direction of each diagonal-shift triple is an exact "
          "CKM gauge:**", "",
          "| triple | `\\|J_F·1\\|` | `\\|J_M·1\\|` | realised value | traceless? |",
          "|--|--|--|--|--|"]
    for name, v in cen["diagonal_shift_traces"].items():
        L.append(f"| `{name}` | `{v['flavor_response_norm']:.1e}` | "
                 f"`{v['mass_response_norm']:.2f}` | "
                 f"`{tuple(round(x, 10) for x in v['realised_value'])}` | "
                 f"{'**yes**' if v['realised_is_traceless'] else 'no'} |")
    L += ["",
          f"So a uniform shift within a block moves masses and **cannot have "
          f"been selected on flavor data**. {cen['symbol_count_is_not_dof_count']}. "
          f"And {cen['diag_shift_plus_collapses_further']}.", ""]

    rnk = detail("F3")
    L += ["## F3 — the decisive object: `rank K`", "",
          f"**`rank K = {rnk['rank_K']}`**, against a physical flavor space of "
          f"dimension {PHYSICAL_FLAVOR_DIMENSION}. `ker J_M` is "
          f"`{rnk['kernel_dimension']}`-dimensional and verified to "
          f"`{rnk['kernel_is_a_kernel']:.1e}`.", "",
          f"Singular values of `K`: "
          f"`{[round(s, 6) for s in rnk['singular_values']]}` — spread "
          f"`{rnk['singular_value_spread']:.0f}×`, no near-degeneracy.", "",
          "Rank is invariant under invertible reparameterisation of `x`, which "
          "is exactly why it is the right object here: no metric, no "
          "pseudoinverse, no chart is being smuggled in. And it is stable "
          "across the whole grid:", "",
          "| step | rcond | dim ker | rank K |", "|--|--|--|--|"]
    for row in rnk["grid"]:
        L.append(f"| `{row['step']:.0e}` | `{row['rcond']:.0e}` | "
                 f"{row['kernel_dimension']} | **{row['rank_K']}** |")
    L += ["", "**Tangent-space construction** — hit an arbitrary target `δy_F` "
          "using only mass-preserving directions. This is algebra on the "
          "matrices above, *not* independent evidence; the nonlinear check is "
          "F3c:", "",
          "| trial | flavor miss | mass disturbance | `\\|δx\\|` |", "|--|--|--|--|"]
    for i, t in enumerate(rnk["tangent_space_construction"]):
        L.append(f"| {i} | `{t['flavor_miss']:.1e}` | "
                 f"`{t['mass_disturbance']:.1e}` | "
                 f"`{t['parameter_excursion']:.1f}` |")
    L += ["",
          f"> **{rnk['verdict']}.**", "",
          f"{rnk['so_there_is_no_invariant_relation']}.", ""]

    restricted = detail("F3b")
    L += ["## F3b — the restricted question, which is the decisive one", "",
          "The headline `rank K` lets *every* coordinate move, including the "
          "eight v3 knobs. But v4 was built **additively over a frozen v3 "
          "lock**, so the surplus question is whether the v4 calibration "
          "freedoms *alone* span the CKM at fixed masses:", "",
          "```", "K_v4 = J_F[:,G] · ker(J_M[:,G])    (v3 knobs held fixed)",
          "```", "",
          "| group | coordinates | `dim ker J_M` | **rank K_v4** | smallest σ |",
          "|--|--|--|--|--|"]
    for row in restricted["rows"]:
        L.append(f"| {row['group']} | {row['coordinates']} | "
                 f"{row['kernel_dimension']} | **{row['rank_K_v4']}** | "
                 f"`{row['smallest_singular_value']:.2e}` |")
    L += ["",
          f"> **{restricted['verdict']}.**", "",
          f"The group includes the `eta_k3k5_minus` **retune** (`5.0 → 5.586`) "
          f"— the round's own rule was to count every quantity selected using "
          f"flavor data regardless of whether the symbol already existed, and "
          f"the first draft's nine-symbol group broke it. Including it also "
          f"**improves the conditioning**: the smallest singular value rises "
          f"`{restricted['conditioning_improvement'][0]:.1e}` → "
          f"`{restricted['conditioning_improvement'][1]:.1e}`.", "",
          f"*{restricted['what_the_first_draft_reported_instead']}.*",
          ""]

    nl = detail("F3c")
    L += ["## F3c — and it is nonlinearly true", "",
          "The tangent-space construction above is algebra on the same "
          "first-order matrices whose rank was just computed — **not** "
          "independent evidence, and the first draft presented it as though it "
          "were. The real check scales one of those directions by `ε`, re-runs "
          "the actual Hamiltonian, and asks whether both errors fall as "
          "`O(ε²)`.", "",
          "| `ε` | mass error | flavor miss |", "|--|--|--|"]
    for row in nl["rows"]:
        L.append(f"| `{row['epsilon']:.2e}` | `{row['mass_error']:.3e}` | "
                 f"`{row['flavor_miss']:.3e}` |")
    L += ["",
          f"Halving `ε` quarters both: mass ratios "
          f"`{[round(r, 3) for r in nl['mass_error_ratios']]}`, miss ratios "
          f"`{[round(r, 3) for r in nl['flavor_miss_ratios']]}` — clean second "
          f"order. So the reachability is a property of the model, not only of "
          f"its Jacobian.", "",
          f"**What it licenses.** {nl['what_it_licenses']}.", ""]

    ab = detail("F4")
    L += ["## F4 — the `φ_h` A/B test", "",
          "The library treats `φ_h = π/k₅` as *derived* and as the source of "
          "CP structure. If holding it fixed dropped the flavor rank from 4 to "
          "3, and the missing direction were the observed CP relation, that "
          "would be evidence that topology removes one calibration freedom.",
          "",
          "| case | coordinates | dim ker `J_M` | **rank K** | singular values |",
          "|--|--|--|--|--|"]
    for label, case in ab["per_case"].items():
        L.append(f"| `{label}` | {case['n_coordinates']} | "
                 f"{case['kernel_dimension']} | **{case['rank_K']}** | "
                 f"`{[round(s, 6) for s in case['singular_values']]}` |")
    d = ab["phi_h_response_direction"]
    L += ["", f"> **{ab['verdict']}.**", "",
          f"**Is it actually a CP handle?** Its `J_F` column is "
          f"`{d['delta_share']:.5f}` along `δ` — "
          + ", ".join(f"`{k} {v:+.4f}`" for k, v in d["by_observable"].items())
          + f". That share **is** chart-independent: "
          f"{ab['why_the_direction_claim_survives_the_chart']}.", "",
          f"**But the `{ab['leading_singular_value_ratio']:.1f}×` is not.** "
          f"{ab['the_ratio_is_chart_dependent']}.", "",
          f"*{ab['relation_to_pr_173']}.*", ""]

    cnt = detail("F5")
    L += ["## F5 — the counting claim", "",
          f"| | claimed | measured |", "|--|--|--|",
          f"| new parameters | {cnt['claimed_new_parameters']} | "
          f"{cnt['symbols_added_phi_h_fixed']} symbols, "
          f"**{cnt['measured_calibration_dimension_phi_h_fixed']}** "
          f"independent flavor directions |",
          f"| new independent observables | "
          f"{cnt['claimed_new_independent_observables']} | "
          f"**≤ {cnt['physical_flavor_dimension']}** (the ceiling) |",
          f"| net predictive surplus | +{cnt['claimed_net_surplus']} | "
          f"**{cnt['measured_net_surplus_phi_h_fixed']}** |", "",
          f"> **{cnt['verdict']}.**", ""]

    led = detail("F6")
    L += ["## F6 — the ledger", "", "| claim | verdict | evidence |",
          "|--|--|--|"]
    for e in led["entries"]:
        L.append(f"| {e['claim']} | **{e['verdict']}** | {e['evidence']} |")
    L += ["", f"**Headline.** {led['headline'].capitalize()}.", "",
          f"**What this does not say.** {led['what_this_does_not_say']}", "",
          f"**Scope.** {led['scope']}", "",
          f"**The recommendation.** {led['the_recommendation']}", ""]
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
           / f"{stamp}_flavor_identifiability_probe")
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2,
                                               default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0 if summary["passed"] == summary["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
