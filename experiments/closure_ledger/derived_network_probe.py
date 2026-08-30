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
    measure_the_closure_residual_has_no_root,
    measure_the_derived_network_ledger,
    measure_the_loop_is_driven_by_the_geometry,
    measure_the_ultraviolet_slope_matches_born,
)


def run_probe() -> dict:
    checks: List[dict] = []

    driven = measure_the_loop_is_driven_by_the_geometry()
    checks.append({
        "id": "N0", "name": "the network loop is driven by the geometry",
        "detail": driven,
        "pass": bool(driven["lambda_modulus_equals_transmission"]
                     and driven["the_backend_is_the_derived_one"]
                     and driven["the_batched_path_equals_the_scalar_network_path"])})

    closure = measure_the_closure_residual_has_no_root()
    checks.append({
        "id": "N1",
        "name": "*** no finite frequency closes carrier and packet together ***",
        "detail": closure,
        "pass": bool(closure["there_is_no_simultaneous_closure"])})

    law = measure_the_closure_residual_has_an_analytic_ultraviolet_law()
    checks.append({
        "id": "N1b", "name": "the closure residual's UV law is analytic",
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

    ledger = measure_the_derived_network_ledger()
    checks.append({
        "id": "N4", "name": "the ledger", "detail": ledger,
        "pass": bool(len(ledger["entries"]) >= 7)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    driven, closure = detail("N0"), detail("N1")
    law, resonance = detail("N1b"), detail("N2")
    ultraviolet, ledger = detail("N3"), detail("N4")

    L: List[str] = [
        "# PR #216's loop, driven by a derived geometry", "",
        f"**{summary['passed']}/{summary['total']} checks pass.**", "",
        "`traversable_throat` computes `T_ℓ(ω)` from a supported traversable 5D "
        "metric. `network.py` carries PR #216's loop eigenvalue. This wires "
        "them, so the closure questions are statements about **the BAM module "
        "itself** rather than a reconstruction beside it.", "",
        "```",
        "Λ_ℓ(ω, Δ) = η_topo · T_ℓ(ω) · e^{iω(d_A + d_B + Δ)}",
        "```", "",
        "`η_topo` is `NetworkThroat.topological_factor` — the module's own deck "
        "orientations and mouth phases. **No separate `tau_th` phase**: a "
        "whole-throat `T` already carries the transit in `arg T`, and adding "
        "one would double-count the Wigner delay.", "",
        "> ## The result", "",
        "> **No finite frequency lets one clock offset serve both a "
        "monochromatic carrier and a localised packet.** The branch-free "
        "residual `C = Arg exp[i(θ − ωθ′)]` has **no root** over "
        f"`{closure['range']}` at {closure['samples']} points, and its decay is "
        f"analytic: `ωC → −∫V_ℓ ds = {law['predicted_limit_of_omega_times_C']:.6f}` "
        "`= −9π/8`. So `C` vanishes as `1/ω` — simultaneous closure is a UV "
        "limit, the same limit in which `|T| → 1`.", "",
        "---", "",
        "## N0 — the loop is the module's own", "",
        "| `ω` | `\\|T\\|` | `\\|Λ\\|` | difference |", "|--|--|--|--|"]
    for row in driven["rows"]:
        L.append(f"| `{row['omega']:g}` | `{row['transmission_modulus']:.10f}` | "
                 f"`{row['lambda_modulus']:.10f}` | `{row['difference']:.1e}` |")
    L += ["",
          f"`|η_topo| = 1` ⟹ `|Λ| = |T|`, confirmed to `1e-12`. The batched scan "
          "path is asserted equal to the scalar `network.py` path, so the "
          "continuous searches below the same object the module exposes.", "",
          "> " + driven["no_extra_transit_phase_is_applied"], ""]

    L += ["## N1 — the branch-free closure residual has no root", "",
          "Eliminating `Δ` between phase closure `ω(D+Δ) + θ = 2πn` and group "
          "closure `D + Δ + θ′ = 0` gives `θ − ωθ′ = 2πn`, so", "",
          "```", "C_ℓ(ω) = Arg exp[i(θ_ℓ − ω θ_ℓ′)]", "```", "",
          "vanishes exactly when one offset serves both. It searches over `n` "
          "automatically.", "",
          "| `ω` | `C_ℓ(ω)` |", "|--|--|"]
    for row in closure["residual_at_probes"]:
        L.append(f"| `{row['omega']:.2f}` | `{row['residual']:+.5f}` |")
    L += ["",
          f"Roots found: **{closure['roots'] or 'none'}**. Smallest `|C|` = "
          f"`{closure['smallest_absolute_residual']:.5f}`.", "",
          "> " + closure["why_this_is_the_invariant_statement"], "",
          "> " + closure["what_it_still_depends_on"], ""]

    L += ["## N1b — and its decay is analytic, not fitted", "",
          f"`{law['closed_form']}`, so `ωC → −∫V_ℓ ds`:", "",
          "| `ω` | `C` | `ωC` |", "|--|--|--|"]
    for w, c, sc in zip(law["omega"], law["residual"], law["omega_times_residual"]):
        L.append(f"| `{w:g}` | `{c:+.6f}` | `{sc:+.6f}` |")
    L += ["",
          f"Predicted `{law['predicted_limit_of_omega_times_C']:.6f}`, reached to "
          f"`{100*law['relative_error_at_the_top']:.2f}%` by `ω = 20`, "
          "monotonically.", "",
          "> **" + law["what_this_settles"] + "**", ""]

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
