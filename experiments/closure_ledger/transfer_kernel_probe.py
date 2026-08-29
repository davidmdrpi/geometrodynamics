"""Probe: the retarded outer→inner transfer kernel on the D = 5 background.

Run:  python -m experiments.closure_ledger.transfer_kernel_probe
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from typing import List

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.tangherlini.transfer_kernel import (  # noqa: E402
    measure_the_exact_background_anchors,
    measure_the_kernel_against_the_time_domain_evolution,
    measure_the_kernel_is_causal,
    measure_the_kernel_reproduces_the_published_ringdown,
    measure_the_scattering_is_well_conditioned,
    measure_the_transfer_is_not_rigid,
    measure_the_transfer_kernel_ledger,
    measure_what_the_causality_gate_caught,
)


def run_probe() -> dict:
    checks: List[dict] = []

    anchors = measure_the_exact_background_anchors()
    checks.append({
        "id": "K0", "name": "three exact D=5 anchors, derived not recalled",
        "detail": anchors,
        "pass": bool(anchors["ell_1_peak_is_exactly_100_over_81"]
                     and anchors["ell_1_peak_radius_squared_is_exactly_9_over_5"]
                     and anchors["c_of_ell_1_is_exactly_nine_quarters"])})

    conditioning = measure_the_scattering_is_well_conditioned()
    checks.append({
        "id": "K1",
        "name": "the real-frequency scattering is well conditioned (unitarity)",
        "detail": conditioning,
        "pass": bool(conditioning["unitarity_holds"]
                     and conditioning["second_order_in_the_step"])})

    causal = measure_the_kernel_is_causal()
    checks.append({
        "id": "K2", "name": "GATE 1 - the kernel is causal",
        "detail": causal, "pass": bool(causal["the_kernel_is_causal"])})

    ringdown = measure_the_kernel_reproduces_the_published_ringdown()
    checks.append({
        "id": "K3",
        "name": "GATE 3 - the kernel carries the published ringdown",
        "detail": ringdown,
        "pass": bool(ringdown["the_ringdown_is_the_published_one"])})

    cross = measure_the_kernel_against_the_time_domain_evolution()
    checks.append({
        "id": "K4",
        "name": "independent check - kernel convolution vs time-domain evolution",
        "detail": cross, "pass": bool(cross["the_two_methods_agree"])})

    rigid = measure_the_transfer_is_not_rigid()
    checks.append({
        "id": "K5", "name": "*** the transfer is NOT rigid ***",
        "detail": rigid,
        "pass": bool(rigid["the_sum_rule_holds"]
                     and rigid["the_kernel_is_not_rigid"])})

    caught = measure_what_the_causality_gate_caught()
    checks.append({
        "id": "K6", "name": "what the causality gate caught",
        "detail": caught, "pass": True})

    ledger = measure_the_transfer_kernel_ledger()
    checks.append({
        "id": "K7", "name": "the ledger", "detail": ledger,
        "pass": bool(len(ledger["entries"]) >= 8)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    rigid = detail("K5")
    L: List[str] = [
        "# The retarded outer→inner transfer kernel, `D = 5` Tangherlini", "",
        f"**{summary['passed']}/{summary['total']} checks pass.**", "",
        "PR #270 deferred this object because it is a ratio of two signals it "
        "could not trust. PR #274 settled which signal was right, against a "
        "published spectrum. This round builds the kernel.", "",
        "A wave sent in from the far region reaches the horizon filtered. That "
        "filter is `T_ℓ(ω)`; in time it is a convolution kernel with "
        "`ψ_trans(v) = (K_ℓ ⋆ ψ_inc)(v)`, `v = t + r*`. Since `T → 1` at high "
        "frequency, `K_ℓ(t) = δ(t) + K_ℓ^reg(t)` — instantaneous part plus "
        "memory.", "",
        "> ## The transfer is not rigid, and not marginally.", "",
        f"| quantity | value | |", "|--|--|--|",
        f"| `∫K_reg dt` | `{rigid['sum_rule_measured']:.6f}` | exact value `−1`, "
        f"a **sum rule** from `T(0) = 0` |",
        f"| `∫\\|K_reg\\| dt` | `{rigid['memory_absolute_mass']:.4f}` | against the "
        f"`δ`'s weight of `1` |",
        f"| `T(ω→0)` | `{rigid['transmission_at_dc']:.3e}` | the barrier blocks "
        f"DC completely |", "",
        "A rigid exchange kernel is `δ(t)`, possibly delayed: whatever enters "
        "leaves undistorted, and a static signal passes perfectly. The real "
        "geometry **blocks a static signal completely**, and does so *entirely* "
        "through the memory term, which exactly cancels the instantaneous one "
        "at DC. In absolute mass the memory is about twice the delta. It is "
        "not a correction to rigid exchange — it is the same size as the thing "
        "it would correct.", "",
    ]

    anchors = detail("K0")
    L += ["## K0 — three exact anchors", "",
          "**1. The barrier peak at `ℓ = 1` is rational.**", "",
          "| `ℓ` | `r²` at peak | `r*` | `V_max` |", "|--|--|--|--|"]
    for row in anchors["barrier_peaks"]:
        L.append(f"| {row['ell']} | `{row['r_squared']:.6f}` | "
                 f"`{row['r_star']:.4f}` | `{row['V_max']:.9f}` |")
    L += ["",
          "At `ℓ = 1`: **`r² = 9/5` and `V_max = 100/81`, exactly.** " +
          anchors["only_ell_1_is_rational"], "",
          "Note `r² → 2`, the photon sphere PR #274 pinned exactly.", "",
          "**2. The integrated potential is exact:** `∫V dr* = ℓ(ℓ+2) + 3/2`. "
          + anchors["the_integral_is_elementary"], "",
          "| `ℓ` | exact | quadrature (truncated) | predicted missing tail |",
          "|--|--|--|--|"]
    for row in anchors["integrated_potential"]:
        L.append(f"| {row['ell']} | `{row['exact']:.3f}` | "
                 f"`{row['quadrature_on_the_truncated_domain']:.4f}` | "
                 f"`{row['missing_tail_beyond_the_outer_edge']:.4f}` |")
    L += ["",
          "The deficit matches the predicted tail in every row, so the domain "
          "truncation is accounted for rather than assumed small.", "",
          "**3. Hence the high-frequency phase is closed-form:** "
          "`T_ℓ(ω) → exp(−i c_ℓ/ω)` with `c_ℓ = (ℓ(ℓ+2)+3/2)/2`, so "
          "**`c₁ = 9/4` exactly**. This is what makes the kernel computable at "
          "all: `T − 1 ~ −i c/ω` decays too slowly to transform numerically, "
          "and knowing `c` lets it be removed analytically instead of windowed "
          "away.", ""]

    cond = detail("K1")
    L += ["## K1 — GATE 2: flux conservation, and why this is not the failed shoot",
          "", "| `ω` | `\\|T\\|` | `\\|R\\|` | `\\|R\\|²+\\|T\\|²−1` |",
          "|--|--|--|--|"]
    for w, t, r, u in zip(cond["omega"], cond["transmission_modulus"],
                          cond["reflection_modulus"], cond["unitarity_residual"]):
        L.append(f"| `{w:g}` | `{t:.8f}` | `{r:.8f}` | `{u:+.1e}` |")
    L += ["",
          f"Worst residual **`{cond['worst_unitarity_residual']:.1e}`**, and "
          "unitarity is imposed nowhere — it is a consequence of the "
          "computation, so it measures it.", "",
          "> " + cond["why_this_is_not_the_failed_shoot"], "",
          "Second order in the spatial step: `|T|` at `ω = 1` gives "
          + ", ".join(f"`{r['transmission_modulus']:.8f}`"
                      for r in cond["step_refinement"])
          + ", successive differences "
          + ", ".join(f"`{d:.1e}`" for d in cond["successive_differences"])
          + ".", ""]

    causal = detail("K2")
    L += ["## K2 — GATE 1: causal support", "",
          "| `t` | `K_reg(t)` |", "|--|--|"]
    for t, k in zip(causal["times_before_zero"], causal["kernel_before_zero"]):
        L.append(f"| `{t:g}` | `{k:+.3e}` |")
    L += ["",
          f"Worst acausal value **`{causal['worst_acausal_value']:.2e}`**, and "
          f"**`{causal['noise_floor_far_from_the_front']:.1e}`** away from the "
          "front.", "",
          "> " + causal["the_exact_zero_is_a_free_error_bar"], "",
          "| `t` | `K_reg(t)` |", "|--|--|"]
    for t, k in zip(causal["times_after_zero"], causal["kernel_after_zero"]):
        L.append(f"| `{t:g}` | `{k:+.3e}` |")
    L += [""]

    ring = detail("K3")
    L += ["## K3 — GATE 3: the kernel carries the published ringdown", "",
          f"Reference (external): `{ring['published'][0]:.8f}"
          f"{ring['published'][1]:+.8f}i`. Source: {ring['source']}.", "",
          "| `dt` | window | fitted `ω` | real err | damping err |",
          "|--|--|--|--|--|"]
    for row in ring["rows"]:
        L.append(f"| `{row['sample_spacing']}` | `{tuple(row['window'])}` | "
                 f"`{row['omega'][0]:.6f}{row['omega'][1]:+.6f}i` | "
                 f"`{100*row['real_relative_error']:.3f}%` | "
                 f"`{100*row['damping_relative_error']:.3f}%` |")
    L += ["",
          f"Band: real part `{100*ring['real_part_band']['min']:.3f}%`–"
          f"`{100*ring['real_part_band']['max']:.3f}%`, damping "
          f"`{100*ring['damping_band']['min']:.3f}%`–"
          f"`{100*ring['damping_band']['max']:.3f}%`.", "",
          "> " + ring["the_band_is_reported_not_the_best_row"], "",
          "> " + ring["scored_against_an_external_value"], ""]

    cross = detail("K4")
    L += ["## K4 — an independent method reproduces the kernel", "",
          "Deep inside, the transmitted wave as a function of `v = t + r*` is "
          "exactly `K ⋆ g`. PR #274's characteristic evolution shares no code "
          "with the transfer matrix, so this is real cross-validation — and it "
          "exposed a subtlety about where the pulse may be launched.", "",
          "| pulse centre | launch `r*` | `V` at launch | max diff | rms diff |",
          "|--|--|--|--|--|"]
    for row in cross["rows"]:
        L.append(f"| `{row['pulse_centre']:g}` | `{row['launch_r_star']:g}` | "
                 f"`{row['potential_at_launch']:.2e}` | "
                 f"`{100*row['peak_relative_max_difference']:.2f}%` | "
                 f"`{100*row['peak_relative_rms_difference']:.2f}%` |")
    L += ["", "> " + cross["what_this_exposed"], ""]

    caught = detail("K6")
    L += ["## K6 — what the causality gate caught", ""]
    for key in ("missing_dc_cell", "gibbs_ringing_from_the_slow_tail"):
        item = caught[key]
        L += [f"**{key.replace('_', ' ')}.** Symptom: {item['symptom']}. "
              f"Cause: {item['cause']}. Fix: {item['fix']}. "
              f"*Would have been read as: {item['would_have_been_read_as']}.*", ""]
    L += ["> " + caught["the_general_point"], ""]

    ledger = detail("K7")
    L += ["## K7 — the ledger", "", "| claim | verdict | evidence |",
          "|--|--|--|"]
    for e in ledger["entries"]:
        L.append(f"| {e['claim']} | **{e['verdict']}** | {e['evidence']} |")
    L += ["", f"**The lesson this round adds.** {ledger['the_lesson_this_round_adds']}",
          "", f"**Scope of the headline.** {rigid['scope_of_that_claim']}",
          "", f"**The next object.** {ledger['the_next_object']}",
          "", f"**Still blocked.** {ledger['still_blocked']}", ""]
    return "\n".join(L)


def main() -> int:
    summary = run_probe()
    text = render_markdown(summary)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs",
                          f"{stamp}_transfer_kernel_probe")
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
