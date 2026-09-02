"""Probe: derive or refute the Born rule from classical BAM measure.

Every number is compared with `docs/classical_born_prereg.md` (commit
`7ff6e41`), written before `geometrodynamics/bulk/mouth_measurement.py`
existed. Run:  python -m experiments.closure_ledger.classical_born_probe
"""

from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk import mouth_measurement as mm   # noqa: E402

TH = np.linspace(0.0, math.pi, 721)
PRE = {"linear_best_miss": 0.109, "two_harmonic_min_miss": 0.50,
       "archimedes": {0.5: 0.25, 0.9: 0.050, 1.0: 0.0, 1.1: 0.046, 2.0: 0.25},
       "mc": {0.3: 0.9777, 1.0: 0.7702, 2.0: 0.2919}}


def run_probe() -> dict:
    b = mm.born(TH)
    h1 = mm.classification_theorem()
    c1 = mm.linear_family_best_fit()
    c2 = {"max_miss": float(np.max(np.abs(mm.induced_probability(mm.intensity_detector(), TH, n=500) - b)))}
    c3 = mm.two_harmonic_natural_weightings()
    c4 = mm.repository_winding_detector()
    arch = {k: float(np.max(np.abs(mm.archimedes_probability(k, TH) - b))) for k in (0.5, 0.9, 1.0, 1.1, 2.0)}
    mc = {t: mm.archimedes_monte_carlo(1.0, t) for t in (0.3, 1.0, 2.0)}
    uni = mm.archimedes_uniformity()
    basin = float(np.max(np.abs(mm.induced_probability(mm.tuned_born_basin(), TH[1:-1]) - b[1:-1])))
    meas = mm.measure_control()
    rev = mm.reversal_control({"born": b, "archimedes": mm.archimedes_probability(1.0, TH),
                               "linear_best": mm.linear_threshold_closed_form(c1["best_alpha_over_rho"], TH)})
    checks = {
        "H1_symmetry_plus_haar_do_not_select_born": h1["all_satisfy_reversal"] and h1["every_f_is_realisable"],
        "C1_linear_family_misses_by_pre_registered_amount": abs(c1["best_max_miss"] - PRE["linear_best_miss"]) < 2e-3 and not c1["reaches_born"],
        "C2_intensity_is_the_step": abs(c2["max_miss"] - 0.5) < 1e-2,
        "C3_two_harmonic_misses_by_at_least_half": min(r["max_miss"] for r in c3) >= PRE["two_harmonic_min_miss"] - 1e-9,
        "C4_repository_detector_is_the_step": c4["induced_f_is_step"],
        "C5_archimedes_born_iff_kappa_one": arch[1.0] < 1e-12 and all(abs(arch[k] - PRE["archimedes"][k]) < 2e-3 for k in arch),
        "C5_monte_carlo": all(abs(mc[t] - PRE["mc"][t]) < 3e-3 for t in mc) and uni["uniform"],
        "basin_control_tuned_reproduces_born": basin < 2e-4,
        "measure_control_fires": meas["measure_matters"],
        "reversal_control": all(v < 1e-12 for v in rev.values()),
    }
    v = mm.verdict(c1["best_max_miss"], [r["max_miss"] for r in c3], arch[1.0], [arch[k] for k in (0.5, 0.9, 1.1, 2.0)])
    return {"timestamp": datetime.now(timezone.utc).isoformat(),
            "prereg": "docs/classical_born_prereg.md @ 7ff6e41",
            "classification": h1, "C1": c1, "C2": c2, "C3": c3, "C4": c4,
            "C5": {"max_miss_by_kappa": arch, "monte_carlo_kappa1": mc, "uniformity": uni},
            "controls": {"tuned_basin_miss": basin, "measure": meas, "reversal": rev},
            "checks": checks, "passed": sum(checks.values()), "total": len(checks),
            "verdict": v,
            "typology": ("C: deterministic hidden outcome. The only classical route to Born found "
                         "(C5) is Bell 1964 / Kochen-Specker 1967: outcomes fixed by (x, y), "
                         "probabilities from ignorance of y. Outcome D (setting-dependent ensemble "
                         "from a derived global boundary problem) was not found and nothing here "
                         "supplies it."),
            "dependency_ledger": {
                "f under fibre Haar": "f( rotational covariance [derived], Haar on S^1 [derived], "
                                      "basin shape [from the coupling: derived for C1-C4; tuned for the control] )",
                "C5 Born": "( Haar on S^2 [chosen: the detector mouth is unprepared], kappa = 1 [chosen], "
                           "D = sign(a.(x + kappa y)) [chosen] )"},
            "what_is_not_claimed": ("Nothing about composition or CHSH. Nothing about the field "
                                    "sector beyond the spin-frame degree of freedom.")}


def render(s: dict) -> str:
    L = [f"# Classical Born-rule probe — {s['passed']}/{s['total']}", "",
         f"Pre-registration: `{s['prereg']}`. Verdict: **`{s['verdict']}`**", "",
         "| candidate | coupling | induced f(θ) under fibre Haar | max miss from Born |", "|--|--|--|--|",
         f"| C1 | sign of a linear functional of the frame | arccos family with plateaus | `{s['C1']['best_max_miss']:.3f}` (best, α/ρ = {s['C1']['best_alpha_over_rho']:.2f}) |",
         f"| C2 | classical Malus intensities | step | `{s['C2']['max_miss']:.3f}` |",
         f"| C3 | spinor + frame harmonics | irregular | `{min(r['max_miss'] for r in s['C3']):.3f}` (best of four) |",
         f"| C4 | repository winding Stern–Gerlach | step | `{s['C4']['max_miss_from_born']:.3f}` |",
         "| C5 | sign(a·(x + κ y)), y Haar on S² | clip((1 + cos θ/κ)/2) | " +
         ", ".join(f"κ={k}: `{v:.3f}`" for k, v in s['C5']['max_miss_by_kappa'].items()) + " |", "",
         f"Monte Carlo at κ = 1: " + ", ".join(f"θ={t}: `{p:.4f}`" for t, p in s['C5']['monte_carlo_kappa1'].items()), "",
         "## Checks", ""]
    for k, ok in s["checks"].items():
        L.append(f"* {'PASS' if ok else 'FAIL'} — {k}")
    L += ["", "## Typology", "", s["typology"], "", "## Dependency ledger", ""]
    for k, v in s["dependency_ledger"].items():
        L.append(f"* `{k}` = {v}")
    L += ["", s["what_is_not_claimed"], ""]
    return "\n".join(L)


def main() -> int:
    s = run_probe()
    text = render(s)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs", f"{stamp}_classical_born_probe")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "probe.json"), "w") as h:
        json.dump(s, h, indent=2, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else str(o))
    with open(os.path.join(outdir, "probe.md"), "w") as h:
        h.write(text)
    return 0 if s["passed"] == s["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
