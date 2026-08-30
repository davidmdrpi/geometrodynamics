"""Probe: PR #216's loop driven by a derived geometry.

Run:  python -m experiments.closure_ledger.derived_network_probe
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from typing import List

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.transaction.derived_network import (  # noqa: E402
    measure_no_perfect_transmission_resonance,
    measure_the_closure_residual_has_an_analytic_ultraviolet_law,
    measure_the_closure_root_is_numerically_converged,
    measure_the_derived_network_ledger,
    measure_the_loop_is_driven_by_the_geometry,
    measure_the_ultraviolet_slope_matches_born,
    measure_what_survives_a_rephasing,
    measure_where_the_loop_closes,
)


def run_probe() -> dict:
    checks: List[dict] = []

    driven = measure_the_loop_is_driven_by_the_geometry()
    legacy = driven["legacy_apis_see_the_derived_transfer"]
    checks.append({
        "id": "N0", "name": "every existing network API is driven by the geometry",
        "detail": driven,
        "pass": bool(driven["lambda_modulus_equals_transmission"]
                     and driven["the_backend_is_the_derived_one"]
                     and driven["the_batched_path_equals_the_scalar_network_path"]
                     and max(legacy.values()) < 1e-12)})

    closure = measure_where_the_loop_closes()
    checks.append({
        "id": "N1",
        "name": "*** where the loop closes, on the DECLARED topology ***",
        "detail": closure,
        "pass": bool(closure["scalar_channel_closes"]
                     and closure["the_two_channels_disagree"])})

    gauge = measure_what_survives_a_rephasing()
    checks.append({
        "id": "N1c",
        "name": "*** and how much of that verdict is gauge, not geometry ***",
        "detail": gauge,
        "pass": bool(gauge["the_derivative_identity_holds"]
                     and gauge["the_swing_is_less_than_one_branch"])})

    law = measure_the_closure_residual_has_an_analytic_ultraviolet_law()
    checks.append({
        "id": "N1b", "name": "the closure function's UV tail law is analytic",
        "detail": law,
        "pass": bool(law["the_limit_is_reached"] and law["the_approach_is_monotone"])})

    resonance = measure_no_perfect_transmission_resonance()
    checks.append({
        "id": "N2", "name": "direct search for R(w) = 0 - a finding, not a theorem",
        "detail": resonance,
        "pass": bool(resonance["no_perfect_transmission_point_found"])})

    ultraviolet = measure_the_ultraviolet_slope_matches_born()
    checks.append({
        "id": "N3", "name": "the UV falloff matches the no-fit Born prediction",
        "detail": ultraviolet,
        "pass": bool(ultraviolet["the_slope_approaches_the_born_value"])})

    converged = measure_the_closure_root_is_numerically_converged()
    checks.append({
        "id": "N5", "name": "the PHASE-sensitive root has its own convergence study",
        "detail": converged,
        "pass": bool(converged["the_root_is_stable"])})

    ledger = measure_the_derived_network_ledger()
    checks.append({
        "id": "N4", "name": "the ledger", "detail": ledger,
        "pass": bool(len(ledger["entries"]) >= 10)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    driven, closure = detail("N0"), detail("N1")
    gauge, law, resonance = detail("N1c"), detail("N1b"), detail("N2")
    ultraviolet, converged, ledger = detail("N3"), detail("N5"), detail("N4")
    scalar = next(r for r in closure["rows"] if r["channel"] == "scalar")

    L: List[str] = [
        "# PR #216's loop, driven by a derived geometry", "",
        f"**{summary['passed']}/{summary['total']} checks pass.**", "",
        "`traversable_throat` computes `T_ℓ(ω)` from a supported traversable 5D "
        "metric. `network.py` carries PR #216's loop eigenvalue. This wires "
        "them — through the APIs that already existed — so the closure "
        "questions are statements about **the BAM module itself**.", "",
        "```",
        "Λ_ℓ(ω, Δ) = η_topo · T_ℓ(ω) · e^{iω(d_A + d_B + Δ)}",
        "```", "",
        "`η_topo` is `NetworkThroat.topological_factor`, and its orientations "
        "are **read off `embedding.topology.make_singlet_pair()`** rather than "
        "chosen. That is not cosmetic: `ConjugatePair` asserts the two mouths "
        "of one throat carry *opposite* orientations, so the scalar channel "
        "has `η_topo = −1`.", "",
        "> ## The result, which reverses an earlier draft of this round", "",
        "> **One clock offset does serve both the carrier and the packet — at "
        f"`ω = {converged['quoted_root']}`.** An earlier draft searched with "
        "`η_topo = +1`, a *chosen* sign, and reported no root at all; the "
        "declared topology shifts `Ψ = θ − ωθ′` by `π` and a root appears.", "",
        "> **But the verdict is gauge dependent, and this says so.** `Ψ` sweeps "
        f"only `{gauge['total_variation']:.4f}` across the band, less than "
        f"`2π = {gauge['two_pi']:.4f}`, so a constant rephasing of the Jost "
        "basis can create or remove the root. Neither *closes* nor *never "
        "closes* is a property of the geometry alone. What is invariant is "
        "`dΨ/dω = −ωθ″` and the total variation.", "",
        "---", "",
        "## N0 — every existing API is the one being driven", "",
        "| `ω` | `\\|T\\|` | `\\|Λ\\|` | difference |", "|--|--|--|--|"]
    for row in driven["rows"]:
        L.append(f"| `{row['omega']:g}` | `{row['transmission_modulus']:.10f}` | "
                 f"`{row['lambda_modulus']:.10f}` | `{row['difference']:.1e}` |")
    legacy = driven["legacy_apis_see_the_derived_transfer"]
    L += ["",
          "`|η_topo| = 1` ⟹ `|Λ| = |T|`. And the dispatch is in `t_AB`, so the "
          "pre-existing entry points agree with the derived transfer exactly:",
          "",
          "| API | disagreement with the derived `T` |", "|--|--|"]
    for name, value in legacy.items():
        L.append(f"| `{name}` | `{value:.1e}` |")
    L += ["",
          "> " + driven["the_dispatch_is_in_t_AB_not_a_parallel_function"], "",
          "> " + driven["the_double_count_is_structurally_impossible"], ""]

    L += ["## N1 — where the loop closes", "",
          "Eliminating `Δ` between phase closure `ω(D+Δ) + θ = 2πn` and "
          "group-delay closure `D + Δ + θ′ = 0` gives", "",
          "```", "Ψ_ℓ(ω) = θ_ℓ − ω θ_ℓ′ = 2πn ,   θ = arg(η_topo T_ℓ)", "```", "",
          "which searches over `n` automatically. (`dΦ/dω = 0` is group-delay "
          "closure *at the carrier* — the necessary first-order condition for "
          "a finite-band packet, not exact packet closure, which would also "
          "constrain the amplitude and every higher derivative.)", "",
          "| channel | `η_topo` | roots | closest approach |", "|--|--|--|--|"]
    for row in closure["rows"]:
        roots = ", ".join(f"`{r:.4f}`" for r in row["roots"]) or "none"
        L.append(f"| {row['channel']} | `{row['topological_factor']:+.0f}` | "
                 f"{roots} | `{row['smallest_distance']:.2e}` |")
    L += ["",
          "> **" + closure["what_this_reverses"] + "**", "",
          "> " + closure["why_two_channels"], "",
          "> " + closure["how_it_is_searched"], ""]

    L += ["## N1c — what a rephasing can and cannot move", "",
          f"`Ψ` runs over `[{gauge['psi_min']:.4f}, {gauge['psi_max']:.4f}]`, a "
          f"total variation of `{gauge['total_variation']:.4f}` against "
          f"`2π = {gauge['two_pi']:.4f}`.", "",
          "> " + gauge["what_is_not_invariant"], "",
          "> " + gauge["what_is_invariant"], "",
          "`dΨ/dω = −ωθ″` verified against its own refinement:", "",
          "| step | relative error | ratio |", "|--|--|--|"]
    for row in gauge["derivative_convergence"]:
        ratio = ("—" if row["ratio_to_previous"] is None
                 else f"`{row['ratio_to_previous']:.2f}`")
        L.append(f"| `{row['step']:.0e}` | `{row['relative_error']:.2e}` | {ratio} |")
    L += ["",
          "> " + gauge["how_much_of_the_constant_is_now_derived"], ""]

    L += ["## N1b — the UV tail law is analytic, not fitted", "",
          f"`{law['closed_form']}`, so `ω[Ψ − arg η_topo] → −∫V_ℓ ds`. The "
          "topological constant must come out first — otherwise `ωΨ` diverges "
          "linearly instead of tending to a limit:", "",
          "| `ω` | `Ψ` | `Ψ − arg η_topo` | `ω(Ψ − arg η_topo)` |",
          "|--|--|--|--|"]
    for w, p, c, sc in zip(law["omega"], law["psi"], law["residual"],
                           law["omega_times_residual"]):
        L.append(f"| `{w:g}` | `{p:+.6f}` | `{c:+.6f}` | `{sc:+.6f}` |")
    L += ["",
          f"Predicted `{law['predicted_limit_of_omega_times_C']:.6f}`, reached to "
          f"`{100*law['relative_error_at_the_top']:.2f}%` by `ω = 20`, "
          "monotonically.", "",
          "> **" + law["why_the_constant_must_be_subtracted"] + "**", "",
          "> " + law["what_this_settles"], "",
          "> **" + law["what_this_does_NOT_settle"] + "**", ""]

    L += ["## N2 — searching for perfect transmission, not assuming it away", "",
          "> " + resonance["this_is_a_finding_not_a_theorem"], "",
          f"Interior minima of `|R|` found: **{resonance['interior_minima_found']}**. "
          f"Smallest `|R|` = `{resonance['smallest_reflection_modulus']:.3e}` at "
          f"`ω = {resonance['smallest_reflection_at_omega']:.2f}`, which is the "
          f"top of the scanned range ({resonance['it_is_at_the_top_of_the_range']}) "
          "— i.e. the minimum is the UV limit, not an interior resonance.", "",
          "> " + resonance["consequence_for_the_loop"], ""]

    L += ["## N3 — the UV falloff against a no-fit oracle", "",
          f"`{ultraviolet['fourier_transform']}`, so first-Born reflection gives "
          f"`1 − |T|² ~ e^{{−4aω}}`: slope `{ultraviolet['predicted_slope']:.1f}` "
          "with nothing fitted.", "",
          "| `ω` | local slope |", "|--|--|"]
    for w, sl in list(zip(ultraviolet["omega"], ultraviolet["local_slope"]))[::3]:
        L.append(f"| `{w:.2f}` | `{sl:+.4f}` |")
    L += ["", "> " + ultraviolet["why_this_is_a_no_fit_oracle"], ""]

    L += ["## N5 — the phase-sensitive root, refined", "",
          "> " + converged["why_unitarity_is_not_enough"], "",
          "| variant | `edge` | `steps` | fd step | root | shift |",
          "|--|--|--|--|--|--|"]
    for row in converged["rows"]:
        L.append(f"| {row['variant']} | `{row['edge']:g}` | `{row['steps']}` | "
                 f"`{row['fd_step']:.0e}` | `{row['root']:.9f}` | "
                 f"`{row['shift_from_baseline']:.1e}` |")
    L += ["", "> " + converged["what_the_spread_licenses"], ""]

    L += ["## N4 — the ledger", "", "| claim | verdict | evidence |", "|--|--|--|"]
    for e in ledger["entries"]:
        L.append(f"| {e['claim']} | **{e['verdict']}** | {e['evidence']} |")
    L += ["", "**What the wiring changes.** " + ledger["what_the_wiring_changes"],
          "", "**Scope.** " + ledger["the_geometry_is_still_an_oracle_not_a_glued_solution"],
          "", "**Still open.** " + ledger["still_open"], ""]
    return "\n".join(L)


def main() -> int:
    summary = run_probe()
    text = render_markdown(summary)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs",
                          f"{stamp}_derived_network_probe")
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
