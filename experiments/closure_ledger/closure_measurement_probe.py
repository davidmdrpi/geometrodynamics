"""Probe: sharp phase closure and measurement dependence at the source.

Compared with `docs/closure_measurement_dependence_prereg.md` (commit
`1b0144e`), written before `geometrodynamics/bulk/closure_measurement.py`.
Run:  python -m experiments.closure_ledger.closure_measurement_probe
"""

from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk import closure_measurement as cm   # noqa: E402

PRE = {0.3: 0.82095, 0.25 * math.pi: 0.53557, 1.0: 0.39850, 0.5 * math.pi: 0.0, 2.0: -0.30350}


def run_probe() -> dict:
    E = {g: cm.correlation(g) for g in PRE}
    closed = {g: float(cm.closed_form_correlation(g)) for g in PRE}
    marg = {g: cm.marginals(g) for g in (0.2, 1.0, 2.5)}
    marg_dev = max(abs(v - 0.5) for m in marg.values() for v in m.values())
    S_std = cm.chsh()
    S_max = cm.chsh_max(n_grid=41)
    barg = cm.bargmann_identity()
    win = {eps: cm.window_monte_carlo(1.0, eps) for eps in (0.4, 0.2, 0.1, 0.05)}
    md = {"non_coplanar": cm.measurement_dependence((0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (1.0, 0.0)),
          "coplanar": cm.measurement_dependence((0.0, 0.0), (1.0, 0.0), (0.5, 0.0), (1.0, 0.0)),
          "coplanar_standard_pairs": cm.measurement_dependence((0.0, 0.0), (0.25 * math.pi, 0.0),
                                                                (0.5 * math.pi, 0.0), (-0.25 * math.pi, 0.0))}
    ctrl = {"two_leg": cm.two_leg_loop_control(), "local": cm.local_detector_control(),
            "width": cm.gaussian_width_control(), "sign": cm.correlation(1.0),
            "pos_variant": {"E(1)": cm.correlation(1.0, "pos"), "S": cm.chsh(variant="pos")}}
    dev_cos = max(abs(E[g] - math.cos(g)) for g in PRE)
    checks = {
        "P1_measurement_dependence_non_coplanar_maximal": md["non_coplanar"]["total_variation"] == 1.0,
        "P1_measurement_dependence_coplanar_in_the_density": md["coplanar"]["total_variation"] > 1e-2,
        "P2_no_signalling": marg_dev < 1e-9,
        "P3_quadrature_equals_closed_form": max(abs(E[g] - closed[g]) for g in PRE) < 1e-4,
        "P3_pre_registered_values": max(abs(closed[g] - PRE[g]) for g in PRE) < 1e-4,
        "P3_not_cos": dev_cos > 0.1,
        "P4_bell_violation_standard_angles": abs(S_std - 2.1423) < 2e-3,
        "P4_below_tsirelson": S_max["below_tsirelson"],
        "P5_signed_is_bargmann_and_cos": barg["max_residual"] < 1e-12 and barg["signed_variant_is_cos"],
        "P6_strict_zero_variant": abs(ctrl["pos_variant"]["E(1)"] - 0.46495) < 1e-4,
        "P7_window_converges": abs(win[0.05] - closed[1.0]) < 5e-3 and abs(win[0.4] - closed[1.0]) > abs(win[0.1] - closed[1.0]),
        "control_two_leg_loop_no_conditioning": ctrl["two_leg"]["closure_automatic"],
        "control_local_detectors_respect_bell": ctrl["local"]["local_bound_respected"],
        "control_gaussian_width_matters": ctrl["width"]["depends_on_sigma"],
    }
    v = cm.verdict(dev_cos, S_max["S_max"], marg_dev, md["coplanar"]["total_variation"])
    return {"timestamp": datetime.now(timezone.utc).isoformat(),
            "prereg": "docs/closure_measurement_dependence_prereg.md @ 1b0144e",
            "E": E, "closed_form": closed, "marginals": marg, "S_standard": S_std, "S_max": S_max,
            "bargmann": barg, "window": win, "measurement_dependence": md, "controls": ctrl,
            "checks": checks, "passed": sum(checks.values()), "total": len(checks), "verdict": v,
            "dependency_ledger": {
                "rho(x|a,b)": "rho( Haar on S^2 [invariant prior; pair direction physical: chosen], "
                              "closure axiom Omega = 0 mod 2pi [repository axiom], geodesic-realignment "
                              "detection model [chosen], coarea conditioning [derived; window limit] )",
                "E(gamma)": "E( rho(x|a,b) [above], outcome signs as history boundary data [chosen: D-type] )",
                "Born": "signed closure measure [identity: Re Bargmann] -- not a probability; imported if used"},
            "typology": ("D: the admissible source ensemble depends on both future settings through the "
                         "closure constraint (support = the great circle through a and b), with exact "
                         "no-signalling. It is not a deterministic detector: the outcome is which closed "
                         "history is realised. Bell's theorem is evaded by measurement dependence, S = "
                         f"{S_std:.4f} at the standard angles, but the correlation is not the quantum one; "
                         "the quantum one is the signed closure measure.")}


def render(s: dict) -> str:
    L = [f"# Closure measurement-dependence probe — {s['passed']}/{s['total']}", "",
         f"Pre-registration: `{s['prereg']}`. Verdict: **`{s['verdict']}`**", "",
         "| γ | E (quadrature) | E (closed form) | cos γ | P(A=+) |", "|--|--|--|--|--|"]
    for g in s["E"]:
        m = s["marginals"].get(g, {"P(A=+)": float("nan")})
        L.append(f"| `{g:.4f}` | `{s['E'][g]:+.5f}` | `{s['closed_form'][g]:+.5f}` | `{math.cos(g):+.5f}` | "
                 f"`{s['marginals'][1.0]['P(A=+)']:.4f}` |")
    L += ["", f"CHSH at (0, π/2, π/4, −π/4): `{s['S_standard']:.4f}`; maximum over settings: "
          f"`{s['S_max']['S_max']:.4f}` at `{tuple(round(v, 3) for v in s['S_max']['settings'])}` "
          f"(below 2√2 = 2.8284: {s['S_max']['below_tsirelson']}).",
          f"Signed variant: `E = cos γ` (deviation `{s['bargmann']['signed_variant_max_deviation_from_cos']:.1e}`), "
          f"`D/4 = Re Tr(P_x P_u P_v)` to `{s['bargmann']['max_residual']:.1e}`; negative-weight fraction of the "
          "closure circle for like outcomes: " + ", ".join(f"γ={g}: `{f:.4f}`" for g, f in
                                                             s['bargmann']['negative_weight_fraction_of_circle_like_outcomes'].items()) + ".",
          "Window Monte Carlo at γ = 1: " + ", ".join(f"ε={e}: `{v:.4f}`" for e, v in s["window"].items())
          + f" → coarea `{s['closed_form'][1.0]:.4f}`.",
          f"Measurement dependence: non-coplanar settings TV `{s['measurement_dependence']['non_coplanar']['total_variation']:.3f}`; "
          f"coplanar settings share the circle, TV of the densities `{s['measurement_dependence']['coplanar']['total_variation']:.4f}` "
          f"(standard pairs `{s['measurement_dependence']['coplanar_standard_pairs']['total_variation']:.4f}`).",
          f"Local-detector control: CHSH `{s['controls']['local']['CHSH_standard']:.4f}`; Gaussian-width control: "
          + ", ".join(f"σ={r['sigma']}: `{r['E']:.4f}`" for r in s['controls']['width']['rows']) + ".", "",
          "## Checks", ""]
    for k, ok in s["checks"].items():
        L.append(f"* {'PASS' if ok else 'FAIL'} — {k}")
    L += ["", "## Typology", "", s["typology"], "", "## Dependency ledger", ""]
    for k, v in s["dependency_ledger"].items():
        L.append(f"* `{k}` = {v}")
    return "\n".join(L) + "\n"


def main() -> int:
    s = run_probe()
    text = render(s)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs", f"{stamp}_closure_measurement_probe")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "probe.json"), "w") as h:
        json.dump(s, h, indent=2, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else str(o))
    with open(os.path.join(outdir, "probe.md"), "w") as h:
        h.write(text)
    return 0 if s["passed"] == s["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
