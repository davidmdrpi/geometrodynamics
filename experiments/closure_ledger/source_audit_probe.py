"""Probe: can any classical field already in BAM keep the finite mouth open?

Every prediction was frozen in `docs/source_audit_prereg.md` **before this
module existed** (commit `765dbaa`). The question is binary, and the answer this
round produces is a **negative result**.

Run:  python -m experiments.closure_ledger.source_audit_probe
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from typing import List

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk.source_audit import (  # noqa: E402
    measure_every_candidate_stress,
    measure_the_dynamic_escape_fails,
    measure_the_flare_out_requirement,
    measure_the_nonminimal_loophole,
    measure_the_source_audit_ledger,
)


def run_probe() -> dict:
    checks: List[dict] = []

    requirement = measure_the_flare_out_requirement()
    checks.append({
        "id": "S0",
        "name": "the neck price is a flare-out theorem, not an artefact of N=1",
        "detail": requirement,
        "pass": bool(requirement["identity_holds"]
                     and requirement["flare_out_is_positive"]
                     and requirement["matches_the_component_result"])})

    candidates = measure_every_candidate_stress()
    checks.append({
        "id": "S1", "name": "*** every existing BAM field gives T_kk >= 0 ***",
        "detail": candidates, "pass": bool(candidates["all_non_negative"])})

    loophole = measure_the_nonminimal_loophole()
    checks.append({
        "id": "S2",
        "name": "*** the nonminimal loophole closes in every dimension ***",
        "detail": loophole, "pass": bool(loophole["conformal_never_flips"])})

    dynamic = measure_the_dynamic_escape_fails()
    checks.append({
        "id": "S3", "name": "*** nonzero K_ij cannot rescue it either ***",
        "detail": dynamic,
        "pass": bool(not dynamic["any_turning_point"]
                     and dynamic["all_focused_within_the_analytic_bound"])})

    ledger = measure_the_source_audit_ledger()
    checks.append({
        "id": "S4", "name": "the verdict and the branches it leaves",
        "detail": ledger, "pass": bool(len(ledger["entries"]) >= 7)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    req, cand = detail("S0"), detail("S1")
    loop, dyn, ledger = detail("S2"), detail("S3"), detail("S4")

    L: List[str] = [
        "# Can any classical field already in BAM keep the throat open?", "",
        f"**{summary['passed']}/{summary['total']} checks pass — and the "
        "result is negative.**", "",
        "Every prediction was frozen in `docs/source_audit_prereg.md` before "
        "this module existed.", "",
        "> ## The verdict", "",
        f"> **{ledger['verdict']}**", "",
        f"> The neck demands `8πG₅T_kk = {req['required_null_stress']:.6f}` at "
        "`R = 1, a = 0.3`. Nine candidates give exactly zero or a manifest "
        "square; the tenth — a nonminimally coupled scalar, the only one whose "
        "sign is not fixed a priori — closes at the defect core in **every** "
        "dimension.", "",
        "---", "",
        "## S0 — the requirement is a flare-out theorem", "",
        "The radial null congruence has a three-dimensional screen in `D = 5`, "
        "so `θ = 3f′/f` and", "",
        "```", "dθ/dλ = −θ²/3 − σ² − R_ab k^a k^b", "```", "",
        f"At the neck `θ = {req['theta_at_neck']:.1e}` and `σ = 0` by spherical "
        f"symmetry, while `dθ/dλ = {req['dtheta_at_neck']:.6f}` — the "
        f"flare-out — against the expected `3/b² = {req['expected_dtheta']:.6f}`. "
        f"The Raychaudhuri identity holds to `{req['raychaudhuri_residual']:.1e}`, "
        "so", "",
        f"```", f"8πG₅T_kk = R_kk = {req['ricci_kk_at_neck']:.6f} = −3/b²",
        "```", "",
        "> " + req["why_this_is_stronger"], "",
        "> **" + req["attribution"] + "**", "",
        "## S1 — every candidate, from its actual stress tensor", "",
        "| | candidate | min `T_kk` | sign | reason |", "|--|--|--|--|--|"]
    for row in cand["rows"]:
        L.append(f"| `{row['id']}` | {row['candidate']} | "
                 f"`{row['null_stress_min']:+.2e}` | **{row['sign']}** | "
                 f"{row['reason']} |")
    L += ["",
          f"Candidates with negative null stress: "
          f"**{cand['candidates_with_negative_null_stress'] or 'none'}**.", "",
          "> **" + cand["the_order_field_cannot_pay_its_own_bill"] + "**", "",
          "## S2 — the nonminimal loophole", "",
          "At a node `q = 0` the prefactor `1 − 8πG₅ξq²` is exactly `1` and "
          "`d²(q²)/dλ² = 2(dq/dλ)²`, so the sign is `sign(1 − 2ξ)`.", "",
          "| coupling | `ξ` | `1 − 2ξ` | null stress |", "|--|--|--|--|"]
    for row in loop["rows"]:
        L.append(f"| {row['coupling']} | `{row['xi']:.4f}` | "
                 f"`{row['one_minus_two_xi']:+.4f}` | {row['null_stress_sign']} |")
    L += ["", f"And the closed form `{loop['closed_form']}`:", "",
          "| `D` | `ξ_c` | `1 − 2ξ_c` |", "|--|--|--|"]
    for row in loop["dimensions"]:
        L.append(f"| {row['D']} | `{row['xi_conformal']:.6f}` | "
                 f"`{row['one_minus_two_xi']:.6f}` |")
    L += ["", "> " + loop["why_the_node_matters"], "",
          "> **" + loop["the_loophole_closes"] + "**", "",
          "## S3 — the dynamic escape", "",
          "If `T_kk ≥ 0` every term on the right of Raychaudhuri is "
          "non-positive, so `θ` is non-increasing and a ray that enters "
          "converging cannot flare back out. Each ray is integrated past its "
          f"own analytic caustic bound, `{dyn['focusing_bound_formula']}`:", "",
          "| profile | `θ₀` | bound | measured caustic | turned `−→+`? |",
          "|--|--|--|--|--|"]
    for row in dyn["rows"]:
        caustic = ("none" if row["measured_caustic"] is None
                   else f"`{row['measured_caustic']:.2f}`")
        L.append(f"| {row['profile']} | `{row['theta_initial']:+.2f}` | "
                 f"`{row['focusing_bound']:.2f}` | {caustic} | "
                 f"**{'YES' if row['turned_positive'] else 'no'}** |")
    L += ["",
          f"Turning points found: **{dyn['turning_cases'] or 'none'}**.", "",
          "> " + dyn["why_the_window_matters"], "",
          "> **" + dyn["the_stronger_statement"] + "**", "",
          "## S4 — the ledger", "", "| claim | verdict | evidence |",
          "|--|--|--|"]
    for entry in ledger["entries"]:
        L.append(f"| {entry['claim']} | **{entry['verdict']}** | "
                 f"{entry['evidence']} |")
    L += ["", "**What this forces.** " + ledger["what_this_forces"], "",
          "**The remaining branches.**", ""]
    for key, value in ledger["the_remaining_branches"].items():
        L.append(f"- **{key}** — {value}")
    L += ["", "**What the audit does not do.** "
          + ledger["what_the_audit_does_not_do"], ""]
    return "\n".join(L)


def main() -> int:
    summary = run_probe()
    text = render_markdown(summary)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs",
                          f"{stamp}_source_audit_probe")
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
