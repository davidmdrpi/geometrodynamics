"""Probe: is the `|D|` density forced by the closure set?

Against `docs/conditioning_variable_prereg.md` (commit `39be3e3`), written
before `geometrodynamics/bulk/conditioning_variable.py`.
Run:  python -m experiments.closure_ledger.conditioning_variable_probe
"""

from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk import conditioning_variable as cv   # noqa: E402


def run_probe() -> dict:
    sector = cv.n_window_is_sector_blind()
    grads = cv.gradient_magnitudes_on_closure()
    families = cv.two_conditioning_families(n=2000000)
    axiom = cv.phase_axiom_citation()
    frozen = cv.downstream_numbers_unchanged()
    velocity = cv.velocity_field_non_uniqueness()
    kappa = cv.kappa_is_route_local()
    v = cv.verdict()

    checks = [
        ("O1   |N| is identical in all four outcome sectors",
         sector["counts_all_equal"] and sector["cross_product_is_pm_q_residual"] < 1e-12),
        ("O2   |grad N| is constant on Gamma -> uniform arclength",
         grads["grad_N_is_constant_on_closure"]),
        ("O3   |grad theta| = |q|/|D| -> the |D| density",
         grads["theta_density_varies"]),
        ("O4   the N-window gives E = 0 exactly", sector["E_is_zero"]),
        ("O5   the phase window tracks the repository law",
         families["phase_window_tracks_the_repository_law"]),
        ("BK   same support, two limits (Borel-Kolmogorov)",
         families["same_support_different_limits"]
         and families["N_window_is_zero_at_every_eps"]),
        ("R1   the justification is in the repository, cited from source",
         axiom["axiom_is_stated_on_phase"]),
        ("R3   NO round 5-7 number moved", frozen["nothing_downstream_moved"]),
        ("C1   div K = 0 leaves continuity unchanged",
         velocity["K_is_divergence_free"] and velocity["continuity_unchanged"]),
        ("C2   int K = 0, so a mean-velocity check cannot exclude it",
         velocity["mean_velocity_check_cannot_exclude_it"]),
    ]
    return {"sector": sector, "gradients": grads, "families": families,
            "axiom": axiom, "frozen": frozen, "velocity": velocity,
            "kappa": kappa, "verdict": v, "ledger": cv.dependency_ledger(),
            "checks": checks,
            "passed": sum(1 for _, ok in checks if ok), "total": len(checks)}


def render(s: dict) -> str:
    L = ["# Round 8 — is the `|D|` density forced by the closure set?", "",
         f"**{s['passed']}/{s['total']} checks pass.**", "",
         "## Verdict", ""]
    for k, val in s["verdict"].items():
        L.append(f"* **{k}** — `{val}`")
    L += ["", "## Checks", ""]
    for name, ok in s["checks"]:
        L.append(f"* {'PASS' if ok else 'FAIL'}  {name}")

    L += ["", "## The two conditioning families", "",
          "| eps | E, `\\|N\\| < eps` window | E, phase window |", "|---|---:|---:|"]
    for r in s["families"]["rows"]:
        L.append(f"| {r['eps']} | {r['E_N_window']:+.6f} | {r['E_phase_window']:+.6f} |")
    L += ["", f"Limit of the phase window (closed form): "
          f"`{s['families']['phase_window_limit_closed_form']:.10f}`. The `N`-window "
          "limit is `0` at every `eps`, because `|N|` does not depend on the "
          "outcome signs.", "",
          f"Counts at `eps = 0.01`: `{s['sector']['counts']}` — equal, not "
          "approximately equal.", ""]

    L += ["## The justification, read out of the source", "",
          f"* `{s['axiom']['file']}` line {s['axiom']['matches'][0][0]}: "
          f"`{s['axiom']['text']}`",
          f"* {s['axiom']['note']}", ""]

    L += ["## Rule 3 — nothing downstream moved", "",
          "| quantity | expected | got | delta |", "|---|---:|---:|---:|"]
    for c in s["frozen"]["checks"]:
        L.append(f"| {c['quantity']} | {c['expected']:+.10f} | {c['got']:+.10f} "
                 f"| {c['delta']:.1e} |")

    vel = s["velocity"]
    L += ["", "## C — the velocity field is not unique", "",
          f"* `max|div K| / max|div J|` = `{vel['max|div K| / max|div J|']:.1e}`",
          f"* `max|div(J+K) - div J| / max|div J|` = "
          f"`{vel['max|div(J+K) - div J| / max|div J|']:.1e}`",
          f"* `int K d^3x` relative to `int|J|` = "
          f"`{vel['integrated_K_relative']:.1e}` — a mean-velocity check "
          "cannot exclude it either",
          "", "## Corrected dependency ledger", ""]
    for row in s["ledger"]:
        L.append(f"* `{row['input']}` — **{row['status']}** ({row['where']})")
    return "\n".join(L) + "\n"


def main() -> int:
    s = run_probe()
    text = render(s)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs",
                          f"{stamp}_conditioning_variable_probe")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "probe.json"), "w") as h:
        json.dump(s, h, indent=2,
                  default=lambda o: o.tolist() if isinstance(o, np.ndarray) else str(o))
    with open(os.path.join(outdir, "probe.md"), "w") as h:
        h.write(text)
    return 0 if s["passed"] == s["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
