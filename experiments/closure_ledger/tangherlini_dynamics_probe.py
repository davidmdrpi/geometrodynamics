"""Full classical 4+1 spherically symmetric Einstein-scalar dynamics.

THE QUESTION
────────────
Every gravity result in this repository so far is stationary, weak-field or
linearized, and THESIS.md has said in the same words across five rounds that
"the strong-field endpoint (a horizon / a resolved throat) is left for full
numerical relativity". Nothing in the tree evolves the Einstein equations in
time. This probe is the first that does, at the highest-symmetry 4+1 problem:
D = 5, spherical symmetry (S^3 angular sector), one minimally coupled massless
scalar, in horizon-penetrating coordinates.

Vacuum is not an option -- Birkhoff in D dimensions makes Tangherlini the
unique spherically symmetric vacuum, so there is nothing to evolve. The scalar
is the dynamical content.

THE GAUGE
─────────
Ingoing Eddington-Finkelstein with areal radius,

    ds^2 = -A(v,r) e^{2 delta(v,r)} dv^2 + 2 e^{delta(v,r)} dv dr + r^2 dOmega_n^2

with n = D - 2 = 3. Ingoing null cones v = const cross the future horizon, so
the chart is horizon-penetrating by construction -- no excision and no
singularity-avoiding lapse is needed to reach it.

T0  THE FIELD EQUATIONS ARE DERIVED, NOT RECALLED. The metric, connection,
    Ricci tensor and Einstein tensor are built in sympy for GENERAL n, the
    massless-scalar stress tensor is added, and the independent components come
    out. Three self-checks, all passing at BOTH n = 3 and n = 2 (the known
    D = 4 system, which is what validates the general-n derivation):
      rr    ->  d_r delta = (kappa/n) r (d_r phi)^2
      vr    ->  (r^{n-1} e^delta A)' = (n-1) r^{n-2} e^delta
      vac   ->  A = 1 - 2m/r^{n-1} with m constant solves vr (Birkhoff)
    The vr result is the surprise: it is not an ODE to be integrated alongside
    A, it is an exact QUADRATURE. Each slice costs three cumulative integrals
    and no implicit solve.

T1  THE HIERARCHY REPRODUCES AN EXACT SOLUTION. phi = cos(w(v-r)) J_1(wr)/r
    solves the flat D = 5 wave equation exactly in these coordinates, so d_v phi
    is known in closed form. The metric comes back flat to 1e-15 and the psi
    quadrature converges at SECOND ORDER (rate 2.010 -> 2.003). That is what
    pins the scheme's order; nothing else in the round does.

T2  TANGHERLINI IS AN EXACT FIXED POINT. With phi = 0 and an interior mass the
    A quadrature is int 2s ds = r^2, so A = 1 - r_h^2/r^2 must come back at
    MACHINE PRECISION rather than to some tolerance that hides a wrong power of
    r. Metric error 1.6e-15, delta identically 0, psi identically 0.

T3  *** THE UNUSED EINSTEIN EQUATION CONVERGES AT SECOND ORDER. *** The
    hierarchy SOLVES rr and vr on every slice, so their residuals are
    identically zero and testing them would be circular. vv is the one
    independent component left over, it contains d_v A, and the code never
    forms d_v A for any other purpose -- so its residual tests the evolution
    rather than restating it. Measured rate 1.989 -> 1.997 -> 1.999 over a
    four-fold refinement. This is the characteristic-scheme analogue of a
    Hamiltonian/momentum constraint test, and the analogy is stated rather than
    assumed.

    Two things this found. A one-sided d_v A capped the measured rate at 1.05;
    it is centred now. And an imposed outer boundary condition -- freezing phi
    at r = R -- left an O(1) vv residual there that did not converge AT ALL:
    the characteristic closure spends its three constants on the gauge and on
    central regularity, so psi(v, R) is already determined and nothing may be
    imposed on it.

T4  *** A REGULAR CENTRE FORBIDS A TRAPPED SURFACE. *** The vr quadrature with
    a regular centre reads

        r^{n-1} e^{delta(r)} A(r) = (n-1) int_0^r s^{n-2} e^{delta(s)} ds ,

    a positive integrand over a positive interval. So A > 0 strictly, for
    r > 0, IDENTICALLY -- and no trapped surface can sit on a regular-centred
    ingoing null slice. Four profile families driven to amplitudes where min A
    falls to 5.6e-03 confirm the code obeys the proof; A never crosses.

    SO HORIZON FORMATION IS NOT OBSERVABLE IN THIS GAUGE with a regular centre,
    and the criterion has to be posed as the loss of central regularity rather
    than as A changing sign. This is a statement about the chart, not about the
    physics: collapse still happens, and the trapped region is reached once the
    centre stops being regular, at which point the slice carries a nonzero
    interior mass constant and this quadrature no longer applies. It is the
    reason production characteristic-collapse codes use OUTGOING null cones or
    excise the centre.

T5  A DISCREPANCY FOUND IN PASSING, REPORTED AND NOT ACTED ON. The
    Schroedinger-form potential for a minimally coupled massless scalar with
    psi = r^{n/2} phi is, at n = 3,
        V = A[(l(l+2) + 3/4)/r^2 + (9/4) r_h^2/r^4]
    and tangherlini.radial.V_tangherlini carries A[l(l+2)/r^2 + 3 r_h^2/r^4].
    The difference is exactly 3A^2/(4r^2), reproduced here to 5.4e-16.

    The flat limit settles which is which with no appeal to anything: at
    r_h -> 0 the regular solution is phi = J_{l+1}(wr)/r, so psi =
    r^{1/2} J_{l+1}(wr), and Bessel's equation gives V -> ((l+1)^2 - 1/4)/r^2 =
    (l(l+2) + 3/4)/r^2 -- the derived form, matched here to 4.3e-16.

    NOTHING IS CHANGED. V_tangherlini is consumed by roughly fifty probes and
    by several derived constants; replacing it is a decision about the
    repository's published numbers and not a side effect of a dynamics round.

T6  WHAT THIS ROUND DID NOT EARN. Two horizon-penetrating time-domain
    constructions were built for a test scalar on a fixed Tangherlini
    background -- a Kerr-Schild slicing of the same chart, and a tortoise
    (t, r*) evolution with the derived potential. Both are stable, both
    CONVERGE, and they do not agree: real parts within 0.3% at l = 1, damping
    rates apart by 37%. A frequency-domain shooting solve did not converge on
    the damping rates either.

    SO NO QUASINORMAL FREQUENCY IS REPORTED, and the retarded outer-to-inner
    transfer function is not built, because it is a ratio of the same two
    signals and would inherit the same unresolved error. Two converged numbers
    that disagree mean a systematic error in at least one construction, and the
    arc has a standing entry for exactly this: A CONVERGED NUMBER IS NOT A
    CORRECT NUMBER.

WHAT IS PUT IN, AND WHAT IS NOT CLAIMED
───────────────────────────────────────
Classical, spherically symmetric, one massless scalar, second-order accurate
and stated as such. No matter model from the rest of the arc appears; no
charge, no winding, no S^3 harmonics above the monopole in the nonlinear
sector. Horizon PERSISTENCE is shown only on a seeded background, where it is
exact; a dynamically formed horizon is not evolved, for the reason T4 gives.
The perturbation spectrum and the transfer function are NOT delivered.

    python -m experiments.closure_ledger.tangherlini_dynamics_probe
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.tangherlini.dynamics import (
    derive_the_field_equations,
    measure_a_regular_centre_forbids_a_trapped_surface,
    measure_the_hierarchy_reproduces_the_exact_flat_mode,
    measure_the_master_potential_disagrees_with_the_repo,
    measure_the_spectrum_is_not_yet_cross_validated,
    measure_the_unused_einstein_equation_converges,
    measure_tangherlini_is_a_fixed_point,
)


def run_probe() -> dict:
    checks: List[dict] = []

    d5 = derive_the_field_equations(3)
    d4 = derive_the_field_equations(2)
    derivation = {
        "spacetime_dim": 5,
        "delta_law": d5["delta_law"],
        "A_law": d5["A_law"],
        "wave_law": d5["wave_law"],
        "why_vv_is_the_monitor": d5["why_vv_is_the_monitor"],
        "checks_at_n_3": {
            "rr_is_the_delta_quadrature": d5["the_rr_equation_is_the_delta_quadrature"],
            "vr_is_the_A_quadrature": d5["the_vr_equation_is_the_A_quadrature"],
            "birkhoff": d5["tangherlini_solves_vr_with_constant_mass"]},
        "checks_at_n_2_the_known_D4_system": {
            "rr_is_the_delta_quadrature": d4["the_rr_equation_is_the_delta_quadrature"],
            "vr_is_the_A_quadrature": d4["the_vr_equation_is_the_A_quadrature"],
            "birkhoff": d4["tangherlini_solves_vr_with_constant_mass"]},
    }
    checks.append({
        "id": "T0", "name": "the field equations are derived, and reduce to D = 4",
        "detail": derivation,
        "pass": bool(all(derivation["checks_at_n_3"].values())
                     and all(derivation["checks_at_n_2_the_known_D4_system"].values()))})

    flat = measure_the_hierarchy_reproduces_the_exact_flat_mode()
    checks.append({
        "id": "T1", "name": "the hierarchy reproduces an exact flat-space mode",
        "detail": flat,
        "pass": bool(flat["the_metric_is_exactly_flat"]
                     and flat["converges_at_second_order"])})

    tang = measure_tangherlini_is_a_fixed_point()
    checks.append({
        "id": "T2", "name": "Tangherlini is an exact fixed point",
        "detail": tang,
        "pass": bool(tang["the_metric_is_exact"] and tang["the_slice_is_static"]
                     and tang["delta_is_identically_zero"])})

    conv = measure_the_unused_einstein_equation_converges()
    checks.append({
        "id": "T3", "name": "*** the unused Einstein equation converges at 2nd order ***",
        "detail": conv,
        "pass": bool(conv["converges_at_second_order"])})

    trap = measure_a_regular_centre_forbids_a_trapped_surface()
    checks.append({
        "id": "T4", "name": "*** a regular centre forbids a trapped surface ***",
        "detail": trap,
        "pass": bool(trap["no_trapped_surface_anywhere"]
                     and trap["A_is_positive_by_the_quadrature"])})

    pot = measure_the_master_potential_disagrees_with_the_repo()
    checks.append({
        "id": "T5", "name": "the master potential disagrees with the repository",
        "detail": pot,
        "pass": bool(pot["gap_matches_the_closed_form"] < 1e-12
                     and pot["flat_limit_matches_bessel"] < 1e-12
                     and pot["nothing_was_changed"])})

    spec = measure_the_spectrum_is_not_yet_cross_validated()
    checks.append({
        "id": "T6", "name": "the spectrum is NOT cross-validated, and is not reported",
        "detail": spec,
        "pass": bool(spec["no_frequency_is_reported"] and spec["both_are_converged"]
                     and len(spec["not_delivered"]) == 2)})

    return {
        "probe": "tangherlini_dynamics",
        "question": "what does the highest-symmetry 4+1 Einstein-scalar system "
                    "do when it is actually evolved, in a horizon-penetrating "
                    "chart -- and does the evolution satisfy the Einstein "
                    "equation it never solves?",
        "answer": "the vv equation, which the hierarchy never uses, converges at "
                  "second order (1.989 -> 1.999); and a regular centre makes "
                  "A > 0 identically, so no trapped surface can sit on a "
                  "regular-centred ingoing null slice -- horizon formation is "
                  "not observable in this chart, which is a statement about the "
                  "chart and not about the physics",
        "headline": {
            "constraint_rate": conv["final_rate"],
            "tangherlini_error": tang["rows"][0]["metric_error"],
            "flat_mode_rate": flat["final_rate"],
            "smallest_A_reached": trap["smallest_A_reached"],
            "a_trapped_surface_ever_appeared": not trap["no_trapped_surface_anywhere"],
            "spectrum_delivered": False,
            "transfer_function_delivered": False,
        },
        "checks": checks,
        "passed": sum(1 for c in checks if c["pass"]),
        "total": len(checks),
    }


def render_markdown(summary: dict) -> str:
    L = ["# Tangherlini dynamics probe — the first evolved Einstein equations", "",
         f"**Question.** {summary['question']}", "",
         f"**Answer.** {summary['answer']}", "",
         f"**{summary['passed']}/{summary['total']} checks pass.**", "",
         "| id | check | result |", "|----|-------|--------|"]
    for c in summary["checks"]:
        L.append(f"| {c['id']} | {c['name']} | {'PASS' if c['pass'] else 'FAIL'} |")

    d = next(c for c in summary["checks"] if c["id"] == "T0")["detail"]
    L += ["", "## T0 — the system, derived for general `n`", "",
          f"* `{d['delta_law']}`", f"* `{d['A_law']}`", f"* `{d['wave_law']}`", "",
          "| check | `n = 3` (D = 5) | `n = 2` (the known D = 4 system) |",
          "|--|--|--|"]
    for k in d["checks_at_n_3"]:
        L.append(f"| {k.replace('_', ' ')} | "
                 f"{'yes' if d['checks_at_n_3'][k] else 'NO'} | "
                 f"{'yes' if d['checks_at_n_2_the_known_D4_system'][k] else 'NO'} |")
    L += ["", f"The `vv` equation is the monitor because {d['why_vv_is_the_monitor']}.",
          ""]

    f = next(c for c in summary["checks"] if c["id"] == "T1")["detail"]
    L += ["## T1 — the exact flat-space mode", "",
          f"`{f['the_exact_solution']}` solves the flat `D = 5` wave equation in "
          f"these coordinates, so `d_v phi` is known in closed form.", "",
          "| points | flat-metric error | `psi` relative error | rate |",
          "|--|--|--|--|"]
    for r in f["rows"]:
        rt = "—" if r["rate"] is None else f"`{r['rate']:.3f}`"
        L.append(f"| {r['points']} | `{r['flat_metric_error']:.1e}` | "
                 f"`{r['psi_relative_error']:.3e}` | {rt} |")
    L += ["", f"*{f['why_only_second_order'].capitalize()}.*", ""]

    t = next(c for c in summary["checks"] if c["id"] == "T2")["detail"]
    L += ["## T2 — Tangherlini is an exact fixed point", "",
          "| points | metric error | max abs `delta` | max abs `psi` |",
          "|--|--|--|--|"]
    for r in t["rows"]:
        L.append(f"| {r['points']} | `{r['metric_error']:.3e}` | "
                 f"`{r['max_abs_delta']:.1e}` | `{r['max_abs_psi']:.1e}` |")
    L += ["", f"{t['what_it_pins'].capitalize()}.", ""]

    c3 = next(c for c in summary["checks"] if c["id"] == "T3")["detail"]
    L += ["## T3 — the Einstein equation the code never solves", "",
          f"The monitored equation is **{c3['the_equation']}**.", "",
          "| points | spacing | max abs `vv` residual | at radius | rate |",
          "|--|--|--|--|--|"]
    for r in c3["rows"]:
        rt = "—" if r["rate"] is None else f"**`{r['rate']:.3f}`**"
        L.append(f"| {r['points']} | `{r['spacing']:.4f}` | "
                 f"`{r['max_abs_vv_residual']:.4e}` | `{r['at_radius']:.2f}` | {rt} |")
    L += ["", f"**What this is not.** {c3['what_it_is_not'].capitalize()}.", "",
          f"**What an imposed outer condition did.** "
          f"{c3['what_an_imposed_outer_condition_did'].capitalize()}.", ""]

    t4 = next(c for c in summary["checks"] if c["id"] == "T4")["detail"]
    L += ["## T4 — a regular centre forbids a trapped surface", "",
          f"    {t4['the_identity']}", "",
          "a positive integrand over a positive interval, so `A > 0` strictly "
          "for `r > 0`, identically. The scan below is not the proof — it is the "
          "check that the code obeys the proof.", "",
          "| profile | amplitude | min `A` | at radius | at `v` | trapped? |",
          "|--|--|--|--|--|--|"]
    for r in t4["rows"]:
        L.append(f"| {r['profile']} | {r['amplitude']:.1f} | "
                 f"`{r['min_A']:.4e}` | `{r['at_radius']:.3f}` | "
                 f"`{r['at_v']:.2f}` | {'YES' if r['trapped'] else 'no'} |")
    L += ["", f"**The consequence.** {t4['the_consequence'].capitalize()}.", "",
          f"*Why it is a chart statement:* {t4['why_it_is_a_chart_statement']}. "
          f"{t4['what_production_codes_do'].capitalize()}.", ""]

    t5 = next(c for c in summary["checks"] if c["id"] == "T5")["detail"]
    L += ["## T5 — a discrepancy found in passing", "",
          f"| | potential |", "|--|--|",
          f"| derived here | `{t5['derived']}` |",
          f"| `tangherlini.radial.V_tangherlini` | `{t5['in_the_repository']}` |",
          f"| difference | `{t5['the_difference']}` |", "",
          f"The difference matches that closed form to "
          f"`{t5['gap_matches_the_closed_form']:.1e}`, and the flat limit matches "
          f"Bessel to `{t5['flat_limit_matches_bessel']:.1e}` — "
          f"{t5['the_flat_limit_is_the_proof']}.", "",
          f"**Nothing was changed.** {t5['why_not'].capitalize()}.", "",
          f"*Caveat:* {t5['caveat']}.", ""]

    t6 = next(c for c in summary["checks"] if c["id"] == "T6")["detail"]
    L += ["## T6 — what this round did not earn", "",
          "| asked for | delivered? |", "|--|--|"]
    for item in t6["what_was_asked"]:
        got = any(item.split(" (")[0] in d_ for d_ in t6["delivered"])
        L.append(f"| {item} | {'yes' if got else '**no**'} |")
    L += ["",
          f"Two horizon-penetrating time-domain constructions, both stable and "
          f"both converged, disagree: Kerr–Schild `{t6['kerr_schild_ell_1']}` "
          f"against tortoise `{t6['tortoise_ell_1']}` at `l = 1`. Real parts "
          f"within {t6['real_parts_agree_to']}; damping rates apart by "
          f"**{t6['damping_rates_differ_by']}**.", "",
          "So **no quasinormal frequency is reported**, "
          "and the transfer function is not built — it is a ratio of the same "
          f"two signals and would inherit the same unresolved error.", "",
          f"*First thing to chase:* {t6['the_first_thing_to_chase']}.", "",
          f"**{t6['the_standing_lesson'].capitalize()}.**", ""]
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
           / f"{stamp}_tangherlini_dynamics_probe")
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2,
                                               default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0 if summary["passed"] == summary["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
