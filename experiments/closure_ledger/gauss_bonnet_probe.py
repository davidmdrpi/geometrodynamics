"""Probe: does Gauss–Bonnet reopen the finite-mouth throat?

Every prediction was frozen in `docs/gauss_bonnet_prereg.md` **before this
module existed** (commit `d47df40`). The answer is no, and it fails in the
direction opposite to the one the branch was invoked for.

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
    measure_the_expansion_breaks_down,
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

    expansion = measure_the_expansion_breaks_down()
    checks.append({
        "id": "B3", "name": "and at that size the truncation is unjustified",
        "detail": expansion,
        "pass": bool(expansion["it_equals_the_leading_term"])})

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
        f"**{summary['passed']}/{summary['total']} checks pass — and the "
        "branch closes.**", "",
        "Frozen in `docs/gauss_bonnet_prereg.md` before this module existed.", "",
        "> ## The verdict", "",
        f"> **{ledger['verdict']}.**", "",
        "> Raychaudhuri does not care about the gravitational action, so "
        "flare-out still forces `R_kk < 0`; only *which tensor supplies it* "
        "changes. The branch's hope was that `α_GB H_kk` supplies the negative "
        "part geometrically. Instead `H_kk` has the **same sign** as `R_kk` at "
        "every neck, so the best-motivated coupling makes the violation worse.",
        "",
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
          "## B3 — and even then, outside the theory's own regime", "",
          "| `α_GB` | `α_GB H_kk / R_kk` |", "|--|--|"]
    for row in exp["rows"]:
        L.append(f"| `{row['coupling']:+.9f}` | `{row['relative_size']:+.4f}` |")
    L += ["",
          f"The required Gauss–Bonnet length is "
          f"`√|α| = {exp['gauss_bonnet_length']:.6f}`, exactly "
          f"`{exp['length_ratio']:.4f} × b` — tied to the throat radius rather "
          "than to a separate short scale.", "",
          "> " + exp["why"], "",
          "## B4 — the ledger", "", "| claim | verdict | evidence |",
          "|--|--|--|"]
    for entry in ledger["entries"]:
        L.append(f"| {entry['claim']} | **{entry['verdict']}** | "
                 f"{entry['evidence']} |")
    L += ["", "**What this closes.** " + ledger["what_this_closes"], "",
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
