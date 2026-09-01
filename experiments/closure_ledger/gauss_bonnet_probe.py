"""Probe: does Gauss–Bonnet reopen the finite-mouth throat?

Every prediction was frozen in `docs/gauss_bonnet_prereg.md` **before this
module existed** (commit `d47df40`). The branch is **narrowed, not closed**:
G6 was withdrawn after review, and its replacement is a stronger constraint.

Run:  python -m experiments.closure_ledger.gauss_bonnet_probe
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from typing import List

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk.gauss_bonnet import (  # noqa: E402
    measure_gauss_bonnet_reinforces_einstein,
    measure_the_negative_coupling_requirement,
    measure_the_gauss_bonnet_ledger,
    measure_the_lanczos_tensor_is_correct,
    measure_the_required_coupling,
)


def run_probe() -> dict:
    checks: List[dict] = []

    validation = measure_the_lanczos_tensor_is_correct()
    checks.append({
        "id": "B0", "name": "the Lanczos tensor vanishes in D = 4, as it must",
        "detail": validation,
        "pass": bool(validation["topological_in_four_dimensions"])})

    sign = measure_gauss_bonnet_reinforces_einstein()
    checks.append({
        "id": "B1",
        "name": "*** Gauss-Bonnet has the SAME sign as Einstein at the neck ***",
        "detail": sign,
        "pass": bool(sign["reinforces"] and sign["ratio_is_everywhere_positive"]
                     and sign["worst_relative_error"] < 1e-2)})

    coupling = measure_the_required_coupling()
    checks.append({
        "id": "B2", "name": "*** the matter NEC would need a NEGATIVE coupling ***",
        "detail": coupling,
        "pass": bool(coupling["no_positive_coupling_works"]
                     and coupling["heterotic_makes_it_worse"])})

    expansion = measure_the_negative_coupling_requirement()
    checks.append({
        "id": "B3",
        "name": "*** a neck-only cancellation is not a wormhole ***",
        "detail": expansion,
        "pass": bool(not expansion["neck_only_satisfies_nec_globally"]
                     and expansion["global_satisfies_nec"]
                     and expansion["tower_terminates_at_gauss_bonnet"])})

    ledger = measure_the_gauss_bonnet_ledger()
    checks.append({
        "id": "B4", "name": "the verdict and what remains",
        "detail": ledger, "pass": bool(len(ledger["entries"]) >= 5)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    val, sign = detail("B0"), detail("B1")
    coup, exp, ledger = detail("B2"), detail("B3"), detail("B4")

    L: List[str] = [
        "# Does Gauss–Bonnet reopen the throat?", "",
        f"**{summary['passed']}/{summary['total']} checks pass — the branch is "
        "narrowed, not closed.**", "",
        "Frozen in `docs/gauss_bonnet_prereg.md` before this module existed.", "",
        "> ## The verdict", "",
        f"> **{ledger['verdict']}.**", "",
        "> Raychaudhuri does not care about the gravitational action, so "
        "flare-out still forces `R_kk < 0`; only *which tensor supplies it* "
        "changes. The branch's hope was that `α_GB H_kk` supplies the negative "
        "part geometrically. Instead `H_kk` has the **same sign** as `R_kk` at "
        "every neck, so the best-motivated coupling makes the violation worse.",
        "",
        "> **An earlier draft of this probe said the branch was closed. It is "
        "not.** A sufficiently negative constant coupling *does* satisfy the "
        "matter NEC — see B3, which replaces a withdrawn claim that the "
        "derivative expansion had broken down. What is closed is the cheap "
        "version.", "",
        "---", "",
        "## B0 — validate the implementation where the answer is known", "",
        "In `D = 4` the Gauss–Bonnet invariant is topological, so `H_ab` must "
        "vanish identically:", "",
        "| metric | `R_kk` | `H_kk` |", "|--|--|--|"]
    for row in val["rows"]:
        L.append(f"| {row['metric']} | `{row['ricci_kk']:+.2e}` | "
                 f"`{row['lanczos_kk']:+.2e}` |")
    L += ["", "> " + val["why"], "",
          "## B1 — the decisive sign", "",
          "```",
          "R_kk = −3(N f″ − N′f′)/(N f)",
          "H_kk = 12(f′² − 1)(N f″ − N′f′)/(N f³)",
          "```", "",
          "The lapse drops out at the neck for the same reason it did in P4: "
          "`N′` multiplies `f′`, and `f′(0) = 0` is what *makes* it a neck. So "
          "for any `N`:", "",
          "```",
          f"R_kk = {sign['ricci_at_neck']:.6f} = −3f″/f₀",
          f"H_kk = {sign['lanczos_at_neck']:.6f} = −12f″/f₀³",
          f"H_kk/R_kk = {sign['ratio_at_neck']:.6f} = 4/f₀²  > 0",
          "```", "",
          "Against an independent Riemann tensor built by numerical "
          "differentiation, sharing no algebra with those closed forms:", "",
          "| `s` | `R_kk` closed | `R_kk` numeric | `H_kk` closed | `H_kk` numeric |",
          "|--|--|--|--|--|"]
    for row in sign["independent_checks"]:
        L.append(f"| `{row['s']:.3f}` | `{row['ricci_closed']:+.4f}` | "
                 f"`{row['ricci_numeric']:+.4f}` | "
                 f"`{row['lanczos_closed']:+.2f}` | "
                 f"`{row['lanczos_numeric']:+.2f}` |")
    L += ["",
          f"Worst relative error `{sign['worst_relative_error']:.1e}` "
          "(finite-difference limited).", "",
          "> **" + sign["why_this_decides_the_branch"] + "**", "",
          "> " + sign["the_ratio_is_the_misner_sharp_parameter"], "",
          "## B2 — the coupling the matter NEC would demand", "",
          "| case | `α_GB` | `8πG₅T_kk` at the neck | NEC? |", "|--|--|--|--|"]
    for row in coup["rows"]:
        L.append(f"| {row['case']} | `{row['coupling']:+.9f}` | "
                 f"`{row['matter_null_stress']:+.4f}` | "
                 f"**{'yes' if row['nec_satisfied'] else 'no'}** |")
    L += ["",
          f"`{coup['threshold_formula']}`, i.e. "
          f"`α_GB ≤ {coup['threshold']:.9f}` here.", "",
          "> **" + coup["why_the_sign_matters"] + "**", "",
          "## B3 — a neck-only cancellation is not a wormhole", "",
          "**Withdrawn first.** An earlier draft read `α_GB H_kk/R_kk = −1` as "
          "the derivative expansion breaking down. In `D = 5` the Lovelock "
          "tower already terminates at Gauss–Bonnet:", "",
          "| Lovelock `k` | status in `D = 5` |", "|--|--|"]
    for row in exp["lovelock_in_five_dimensions"]:
        L.append(f"| {row['order']} | {row['status']} |")
    L += ["", "> **" + exp["what_was_withdrawn"] + "**", "",
          "**What replaces it.** The NEC must hold *along* the throat, and the "
          "neck is its easiest point:", "",
          "| requirement | `α_GB` | `min T_kk` over throat | NEC everywhere? |",
          "|--|--|--|--|",
          f"| neck only, `−b²/4` | `{exp['neck_threshold']:+.9f}` | "
          f"`{exp['neck_only_min_over_throat']:+.2f}` | "
          f"**{'yes' if exp['neck_only_satisfies_nec_globally'] else 'no'}** |",
          f"| whole throat, `−R²/4` | `{exp['global_threshold']:+.9f}` | "
          f"`{exp['global_min_over_throat']:+.6f}` | "
          f"**{'yes' if exp['global_satisfies_nec'] else 'no'}** |", "",
          f"The mouth sets `f_m⁴/(4b²) = R²/4` **exactly**, independent of `a`, "
          f"so the global requirement is `{exp['ratio_global_to_neck']:.0f}×` "
          "the neck-only one here — and the Gauss–Bonnet length must be "
          f"`√|α| ≥ {exp['length_over_bulk_radius']:.2f} R`, half the radius of "
          "the closed universe rather than a short-distance scale.", "",
          "> " + exp["what_replaces_it"], "",
          "> **" + exp["the_honest_verdict"] + "**", "",
          "## B4 — the ledger", "", "| claim | verdict | evidence |",
          "|--|--|--|"]
    for entry in ledger["entries"]:
        L.append(f"| {entry['claim']} | **{entry['verdict']}** | "
                 f"{entry['evidence']} |")
    L += ["", "**What this narrows.** " + ledger["what_this_narrows"], "",
          "**What remains.**", ""]
    for key, value in ledger["what_remains"].items():
        L.append(f"- **{key}** — {value}")
    L += ["", "**Not refuted here.** " + ledger["not_refuted_here"], ""]
    return "\n".join(L)


def main() -> int:
    summary = run_probe()
    text = render_markdown(summary)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs",
                          f"{stamp}_gauss_bonnet_probe")
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
