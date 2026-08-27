"""The local geometry of the quark fit manifold: response Jacobian and SVD.

THE DURABLE RESULT, FIRST
─────────────────────────
PR #272's residual round measured one knob at a time, found the quark residuals
individually right and jointly wrong, and guessed the gap was a scalar relation
between transport and resistance. That guess was wrong, and so was the object.
The right one is the response map

    J_ia = d ln(m_i / m_d) / d ln p_a ,   i in {s, c, b, t}

Three statements survive review. None depends on a choice of metric in
parameter space, which is why they are the ones to build on:

1. FOUR SCORED MASSES CANNOT IDENTIFY THE CURRENT PARAMETERISATION.
   rank J = 4 exactly -- capped by the observable count -- while eight knobs
   carry first-order response and three more are gauge or quadratic. J^T J has
   four exact zero eigenvalues. pinhole, transport, resistance and N are not
   separately constrained by the mass ladder; at most four combinations are.

2. THE POSITIVE-KNOB JACOBIAN IS NUMERICALLY FULL ROW RANK, converged to five
   digits across three decades of step size, with no adiabatic-relabelling
   noise. It is a real derivative, not a difference artefact.

3. PER-KNOB PROXIMITY TO A FITTED VALUE DOES NOT DETERMINE THE EFFECT ON THE
   SPECTRUM. The three radial residuals do not compose linearly into the
   ladder, and their individual agreement carries no information about the sign
   or size of their joint effect.

WHAT THIS ROUND HAD TO WITHDRAW
───────────────────────────────
Four things in the first draft were wrong or overclaimed. Recorded because the
programme's rule is that wrong first drafts get written down.

  THE phase Z2 IS MODEL-WIDE, NOT LOCK-ONLY. H(-phi, p) = conj(H(phi, p)) for
  ARBITRARY p, and H Hermitian means conj(H) = H^T is isospectral. The spectrum
  is even in phi at every mixing (0.0 on eigenvalues AND on anchored masses at
  p = 0.05, 0.3). The first draft mistook "matrix not equal" for "spectrum not
  even". The two Z2s are of different kinds: phase is ANTIUNITARY on the
  spectrum, partition_mixing is UNITARY-CONJUGATION.

  THE HEADLINE PROJECTION USED THE WRONG DISPLACEMENT. g_legacy and g_corrected
  are candidate displacements FROM THE LOCK. The move the correction makes is
  Delta g = g_corrected - g_legacy, and cos(J.Delta g, r - J.g_legacy) = +0.873
  -- TOWARD the residual the legacy triple leaves, not away. The first draft's
  cos = -0.616 is true about g_corrected and an invalid basis for "the
  correction moves away from the data". WITHDRAWN.

  AND THE TWO OBJECTIVES DISAGREE, which is substance, not technicality:
      L2 log-residual   0.0548 -> 0.0433   IMPROVES
      max rel error      3.19% ->  3.80%   WORSENS
  Both reported; neither privileged.

  THE MIN-NORM GEOMETRY IS METRIC-DEPENDENT. delta_x_min, the 98.35% share, the
  right singular vectors and the "invisible fraction" are Euclidean in the
  chosen eight log coordinates and change under column rescaling. Scoped
  throughout; the physical conclusions ride the observable-space quantities.

  LEAVE-ONE-SPECIES-OUT WAS NOT A PREDICTIVE HOLDOUT. ker(J_keep) is
  5-dimensional and the held-out row is reachable from it in every case: a
  kernel shift fits the withheld mass to ~1e-15 at 1.02-4.77x the norm. The
  -10.4% measures the MINIMUM-LOG-NORM REGULARISER, not the Hamiltonian.
  Reported now as a regulariser diagnostic.

THE NULL STRUCTURE, DERIVED RATHER THAN OBSERVED
────────────────────────────────────────────────
  action_base       EXACT INVARIANCE, ALL ORDERS. H(a) = H(0) + a*I, removed by
                    the min_eigenvalue spectrum-zero subtraction -- UPSTREAM of
                    the d-anchor, so it does not reappear in an absolute-scale
                    observable either. Dropped from the parameter space.

  phase             ANTIUNITARY Z2 of the spectrum, quadratic at zero.
  partition_mixing  UNITARY-CONJUGATION Z2, H(-p) = D H(p) D†, quadratic.

Both handled by the chart q = x^2, licensed by an exact 10x -> 100x scaling.

ON "SLOPPINESS"
───────────────
The four NONZERO singular values span only 22.6x, so there is no long hierarchy
among identifiable directions -- but the rank deficiency is exact. The dominant
pathology is STRUCTURAL NON-IDENTIFIABILITY, not ill-conditioning. And which
knobs drift is itself structure: identifiable share runs from 1.0000
(uplift_asymmetry) to 0.0003 (eta_k3k5_minus), so "any knob would drift" is
false.
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
    LOG_KNOBS, SCORED_SPECIES,
    measure_the_minimum_norm_regulariser_does_not_predict_a_held_out_mass,
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
                     # the SPECTRUM is even in phase at EVERY mixing ...
                     and nul["phase"]["max_spectrum_and_anchored_mass_gap"] == 0.0
                     and nul["phase"]["max_conjugation_gap_over_all_mixings"] == 0.0
                     # ... even though the MATRIX is not
                     and nul["phase"]["max_matrix_gap_over_all_mixings"] > 1e-3
                     and nul["phase"]["the_symmetry_is_model_wide_not_lock_only"]
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
        "name": "*** rank 4 against 8 first-order knobs: 4 null directions ***",
        "detail": svd,
        "pass": bool(svd["rank"] == 4
                     and svd["rank_is_capped_by_the_observable_count"]
                     and svd["n_first_order_null_directions"] == 4
                     and svd["gram_exact_zeros"] == 4
                     and svd["no_long_hierarchy_among_nonzero_singular_values"])})

    rep = measure_which_directions_repair_the_masses()
    checks.append({
        "id": "R3", "name": "the repair rides the WEAKEST singular direction",
        "detail": rep,
        "pass": bool(rep["the_repair_rides_the_weakest_direction"]
                     and rep["dominant_share"] > 0.9
                     and rep["linear_residual_after_repair"] < 1e-12
                     # the exact nonlinear consequence, computed not quoted
                     and rep["exact_max_rel_err_percent_after_repair"] < 0.05
                     and rep["largest_single_knob_move_percent"] < 2.0)})

    prj = measure_the_geometric_displacement_against_what_the_masses_want()
    checks.append({
        "id": "R4",
        "name": "*** the operator-induced move points TOWARD the residual ***",
        "detail": prj,
        "pass": bool(
            prj["the_operator_induced_move"]["cos_move_vs_residual_after_legacy"]
            > 0.8
            and prj["the_two_objectives_disagree"]["l2_improves"]
            and prj["the_two_objectives_disagree"]["max_rel_err_worsens"])})

    hld = measure_the_minimum_norm_regulariser_does_not_predict_a_held_out_mass()
    checks.append({
        "id": "R5",
        "name": "the holdout measures the REGULARISER, not the Hamiltonian",
        "detail": hld,
        "pass": bool(hld["every_held_out_species_is_reachable"]
                     and hld["a_kernel_shift_fits_the_held_out_mass_exactly"])})

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
          f"| `phase` | **{nul['phase']['status']}** | "
          f"`H(−φ, p) = H(φ, p)*` for *arbitrary* `p`; `H` Hermitian ⟹ "
          f"isospectral | "
          f"`{nul['phase']['max_conjugation_gap_over_all_mixings']:.1f}` |",
          f"| `partition_mixing` | **{nul['partition_mixing']['status']}** | "
          f"`H(−p) = D H(p) D†`, `D = diag(σ)` | "
          f"`{nul['partition_mixing']['conjugation_gap']:.1f}` |", "",
          "`action_base` is removed **upstream of the `d`-anchor** — the "
          "zero-point subtraction does it, and `spectrum_zero_mode="
          "\"action_base\"` kills it too. So it does *not* reappear in an "
          "absolute-scale observable of this model, and it is dropped from the "
          "first-order parameter space.", "",
          "**Corrected from the first draft.** That draft called the `phase` "
          "`Z₂` *a property of the lock*, because the **matrix** stops "
          "satisfying `H(−φ) = H(φ)` once `partition_mixing ≠ 0` "
          f"(`{nul['phase']['max_matrix_gap_over_all_mixings']:.3f}` at "
          f"`φ = 0.37, p = 0.3`). But the **spectrum** is even at every "
          f"mixing — `{nul['phase']['max_spectrum_and_anchored_mass_gap']:.1f}` "
          "on the eigenvalues *and* on the anchored masses. Matrix inequality "
          "is not spectral asymmetry.", "",
          f"> {nul['the_two_Z2s_are_different_kinds']}.", "",
          "Both are handled by the chart `q = x²`, under which `∂ln m/∂q` is "
          "finite and constant:", "",
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
          f"first-order knobs, so **{svd['n_first_order_null_directions']} "
          f"directions in these coordinates move no observable at all**; "
          f"`JᵀJ` carries `{svd['gram_exact_zeros']}` exact zeros. The "
          f"identifiable tangent space is 4-dimensional — the *model* has more "
          f"local non-identifiability still, since the two quadratic "
          f"coordinates are held fixed here.", "",
          "| `a` | `σ_a` | `σ_a/σ_1` |", "|--|--|--|"]
    for a, (s, ratio) in enumerate(zip(svd["singular_values"],
                                       svd["sigma_ratios"])):
        L.append(f"| {a+1} | `{s:.4f}` | `{ratio:.4f}` |")
    L += ["",
          f"The four **nonzero** singular values span only "
          f"`{svd['nonzero_singular_value_ratio']:.1f}×` — no long hierarchy "
          f"among the identifiable directions. That is *not* the same as "
          f"\"not a sloppy model\": "
          f"{svd['the_dominant_pathology_is_structural_non_identifiability']}."
          f"\n\n*Scope: {svd['metric_scope']}. Rescaling a column changes "
          f"every right-space quantity below.*",
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
    for knob, v in svd["quadratic_directions_add_no_observable_direction"].items():
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
          f"*(Scope: {rep['metric_caveat']}.)*", "",
          f"The correction has `|δ ln p| = {rep['delta_x_min_norm']:.6f}`, its "
          f"largest single knob moving "
          f"`{rep['largest_single_knob_move_percent']:.2f}%`, and an **exact "
          f"nonlinear re-run** — `exp(δ)` applied to the lock, spectrum "
          f"re-solved — gives "
          f"`{rep['exact_max_rel_err_percent_after_repair']:.4f}%` from the "
          f"lock's `{rep['locked_max_rel_err_percent']:.2f}%`.", "",
          f"> {rep['verdict']}.", ""]

    prj = detail("R4")
    move = prj["the_operator_induced_move"]
    obj = prj["the_two_objectives_disagree"]
    L += ["## R4 — does the operator correction move toward the data?", "",
          "Two different questions live here, and the first draft ran them "
          "together.", "",
          "**(a) Which candidate displacement from the lock lands nearer the "
          "data?**", "",
          "| operator | `\\|δx_geom\\|` | `cos(J·g, r)` | L2 log-residual | "
          "max rel err | null-space fraction |", "|--|--|--|--|--|--|"]
    for name, label in (("legacy", "legacy"), ("scalar_correct", "corrected")):
        o = prj["per_operator"][name]
        L.append(f"| {label} | `{o['delta_x_geom_norm']:.6f}` | "
                 f"`{o['cos_to_r_observable_space']:+.4f}` | "
                 f"`{o['l2_log_residual']:.4f}` | "
                 f"`{o['max_rel_err_percent_linear']:.2f}%` | "
                 f"`{o['fraction_of_displacement_in_the_null_space']:.4f}` |")
    L += ["",
          "**(b) Which way does the *correction itself* move?** That is "
          "`Δg = g_corrected − g_legacy`, compared against the residual the "
          "legacy triple leaves.", "",
          f"| `\\|Δg\\|` | `\\|r − J·g_legacy\\|` | `cos(J·Δg, r − J·g_legacy)` |",
          "|--|--|--|",
          f"| `{move['delta_g_norm']:.6f}` | "
          f"`{move['residual_left_by_the_legacy_triple_norm']:.6f}` | "
          f"**`{move['cos_move_vs_residual_after_legacy']:+.4f}`** "
          f"(`{move['angle_deg']:.1f}°`) |", "",
          f"> **{prj['withdrawn']}.**", "",
          f"And the two objectives disagree, which is substance rather than "
          f"bookkeeping: the L2 log-residual **improves** "
          f"(`{obj['l2_log_residual']['legacy']:.4f}` → "
          f"`{obj['l2_log_residual']['scalar_correct']:.4f}`) while the "
          f"repository's max-relative-error score **worsens** "
          f"(`{obj['max_rel_err_percent_linear']['legacy']:.2f}%` → "
          f"`{obj['max_rel_err_percent_linear']['scalar_correct']:.2f}%`). "
          f"Both are true; the direction of \"improvement\" is "
          f"objective-dependent and no metric-free claim is available.", "",
          f"**What survives.** {prj['what_survives']}.", ""]

    hld = detail("R5")
    L += ["## R5 — the holdout measures the regulariser, not the Hamiltonian",
          "",
          "For each holdout `J_keep` is `3×8` with a five-dimensional kernel, "
          "and because the full `J` has rank 4 the held-out row is **not** in "
          "the span of the other three. So a `z ∈ ker(J_keep)` moving the "
          "withheld species always exists.", "",
          "| held out | `dim ker` | `\\|j_h·Z\\|` | min-norm miss | after a "
          "kernel shift | norm cost |", "|--|--|--|--|--|--|"]
    for row in hld["rows"]:
        rep2 = row["after_a_kernel_shift"]
        L.append(f"| `{row['held_out']}` | `{row['kernel_dimension']}` | "
                 f"`{row['held_out_row_reachable_from_kernel']:.4f}` | "
                 f"`{row['min_norm_held_out_error_percent']:+.3f}%` | "
                 f"**`{rep2['held_out_error_percent']:+.1e}%`** | "
                 f"`{rep2['norm_ratio_to_min_norm']:.2f}×` |")
    L += ["",
          f"> **{hld['verdict']}.**", "",
          "The first draft read the `−10.4%` as a failed prediction. It is "
          "not one: `δ + λz` fits the withheld mass to machine precision while "
          "holding the other three exact, at between "
          f"`{hld['at_norm_cost_between'][0]:.2f}×` and "
          f"`{hld['at_norm_cost_between'][1]:.2f}×` the minimum norm.", ""]

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
