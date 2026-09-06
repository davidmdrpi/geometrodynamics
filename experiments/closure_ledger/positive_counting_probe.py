"""Probe: can positive classical history counting reach the quantum correlation?

Against `docs/positive_counting_prereg.md` (commit `a987342`), written before
`geometrodynamics/bulk/positive_counting.py`.
Run:  python -m experiments.closure_ledger.positive_counting_probe
"""

from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk import positive_counting as pc   # noqa: E402


def run_probe() -> dict:
    red = pc.reduction_residual()
    marg = pc.marginals_are_automatic()
    q1 = pc.threshold_phi_reaches_four()
    minimal = pc.solve_minimal_degree()
    glob = pc.no_global_polynomial()
    sym = pc.symmetry_conditions(6)
    frozen = pc.downstream_numbers_unchanged()
    v = pc.verdict()

    angles = (0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
    quantum_rows = [{"gamma": g,
                     "E": pc.correlation_from_phi(pc.quantum_phi, g),
                     "target": -math.cos(g),
                     "err": abs(pc.correlation_from_phi(pc.quantum_phi, g)
                                + math.cos(g))} for g in angles]
    known = {"|D| (round 5)": pc.chsh_from_phi(np.abs),
             "D (round 6, signed)": pc.chsh_from_phi(lambda d: d),
             "D^2(1-D/5) (round 9)": pc.chsh_from_phi(pc.quantum_phi)}

    checks = [
        ("O1  reduction D = t + sqrt(2t) cos psi holds as a law",
         red["reduction_holds"] and red["t_pair_sums_to_two"]),
        ("O5  marginals are 1/2 for every member of the class",
         marg["marginals_forced_half"]),
        ("Q1  a nonnegative threshold weight reaches CHSH = 4",
         q1["reaches_algebraic_maximum"]),
        ("Q2  an explicit nonnegative Phi reproduces -cos gamma",
         max(r["err"] for r in quantum_rows) < 1e-12),
        ("Q2  and gives exactly Tsirelson",
         abs(known["D^2(1-D/5) (round 9)"] - 2*math.sqrt(2)) < 1e-12),
        ("S1  degree 3 is the minimal nonnegative solution",
         minimal["combined_is_even"] and minimal["phi_nonnegative_on_range"]
         and minimal["degree_2_forces_a2_zero"]),
        ("S2  no polynomial solution is globally nonnegative",
         glob["every_even_top_degree_is_forced_zero"]),
        ("S3  G_n has degree n-1 with leading coefficient 1",
         sym["leading_coefficients_are_one"] and sym["G1_is_constant"]),
        ("R5  nothing in rounds 5-8 or #283 moved",
         frozen["nothing_downstream_moved"]),
    ]
    return {"reduction": red, "marginals": marg, "q1": q1, "minimal": minimal,
            "global": glob, "symmetry_n_conditions": sym["n_conditions"],
            "quantum_rows": quantum_rows, "known_chsh": known,
            "frozen": frozen, "verdict": v, "ledger": pc.dependency_ledger(),
            "checks": checks,
            "passed": sum(1 for _, ok in checks if ok), "total": len(checks)}


def render(s: dict) -> str:
    L = ["# Round 9 — can positive counting reach the quantum correlation?", "",
         f"**{s['passed']}/{s['total']} checks pass.**", "", "## Verdict", ""]
    for k, val in s["verdict"].items():
        L.append(f"* **{k}** — `{val}`")
    L += ["", "## Checks", ""]
    for name, ok in s["checks"]:
        L.append(f"* {'PASS' if ok else 'FAIL'}  {name}")

    L += ["", "## The class spans everything", "",
          "| weight `Φ(D)` | nonnegative? | standard-angle CHSH |",
          "|---|---|---:|"]
    signs = {"|D| (round 5)": "yes", "D (round 6, signed)": "**no** (comparison)",
             "D^2(1-D/5) (round 9)": "yes"}
    for k, val in s["known_chsh"].items():
        L.append(f"| `{k}` | {signs[k]} | {val:.10f} |")
    L += [f"| threshold bump at `1+√2` | yes | "
          f"{s['q1']['rows'][0]['CHSH']:.10f} |", "",
          "Every **nonnegative** row has marginals exactly `1/2`; the signed row "
          "is listed only to show the quantum value was previously reached only "
          "that way.", "",
          "## The explicit quantum witness `Φ(D) = D²(1 − D/5)`", "",
          "| `γ` | `E` | `−cos γ` | error |", "|---|---:|---:|---:|"]
    for r in s["quantum_rows"]:
        L.append(f"| {r['gamma']} | {r['E']:+.12f} | {r['target']:+.12f} "
                 f"| {r['err']:.1e} |")
    m = s["minimal"]
    L += ["", f"* minimal degree is 3: `G_2` odd part `{m['G2_odd_u']}`, `G_3` odd "
          f"part `{m['G3_odd_u']}`, condition `{m['condition']}`",
          f"* combined `G(u) = {m['combined_G_in_u'][0]} + "
          f"({m['combined_G_in_u'][2]}) u²` — even about `t = 1`",
          f"* `min Φ` on `[-1/2, 4]` = `{m['min_phi_on_physical_range']:.2e}`",
          f"* {s['global']['conclusion']}", "", "## Dependency ledger", ""]
    for row in s["ledger"]:
        L.append(f"* `{row['input']}` — **{row['status']}** ({row['where']})")
    return "\n".join(L) + "\n"


def main() -> int:
    s = run_probe()
    text = render(s)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs",
                          f"{stamp}_positive_counting_probe")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "probe.json"), "w") as h:
        json.dump(s, h, indent=2,
                  default=lambda o: o.tolist() if isinstance(o, np.ndarray) else str(o))
    with open(os.path.join(outdir, "probe.md"), "w") as h:
        h.write(text)
    return 0 if s["passed"] == s["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
