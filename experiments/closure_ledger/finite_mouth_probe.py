"""Probe: the finite-mouth scalar-flat handle, against six frozen predictions.

Every check below compares a computation against a number written down in
`docs/finite_mouth_prereg.md` **before this module existed** (commit `ca07204`).
That is the difference between this probe and the 45 all-passing runs the audit
found: these checks can fail, and P1 and P4 additionally run active
falsification attempts rather than confirmations.

Run:  python -m experiments.closure_ledger.finite_mouth_probe
"""

from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone
from typing import List

import numpy as np

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk import finite_mouth as fm  # noqa: E402

CASES = ((1.0, 0.30), (1.0, 0.80), (2.5, 0.30), (1.0, 0.05))


def check_p1() -> dict:
    """P1 — the geometry is parameter-free, and perturbing it breaks matching."""
    rows, worst = [], 0.0
    for R, a in CASES:
        g = fm.geometry(R, a)
        rows.append({
            "R": R, "a": a,
            "b": g.neck_radius, "b_exact": R * math.sin(a) ** 2,
            "S": g.half_length, "S_exact": R * math.sin(a) * math.cos(a),
            "L": g.proper_length, "L_exact": R * math.sin(2.0 * a),
            "f_m": g.mouth_radius,
            "match_radius": abs(float(fm.throat_radius(g.half_length, g.neck_radius))
                                - R * math.sin(a)),
            "match_slope": abs(float(fm.throat_radius_derivative(
                g.half_length, g.neck_radius)) - math.cos(a)),
        })
        worst = max(worst, rows[-1]["match_radius"], rows[-1]["match_slope"],
                    abs(rows[-1]["L"] - rows[-1]["L_exact"]))
    # active falsifier: is there any OTHER (b, S) that also matches?
    R, a = 1.0, 0.30
    g = fm.geometry(R, a)
    perturbed = []
    for db in (0.05, -0.05):
        b2 = g.neck_radius * (1.0 + db)
        # keep the areal radius correct, then the slope must fail
        s2 = math.sqrt(max(g.mouth_radius ** 2 - b2 ** 2, 0.0))
        perturbed.append({
            "delta_b_fraction": db,
            "slope_error": abs(s2 / g.mouth_radius - math.cos(a))})
    return {
        "rows": rows, "worst_residual": worst,
        "holds": bool(worst < 1e-12),
        "perturbations_break_matching": perturbed,
        "no_free_parameter_survives": bool(
            all(p["slope_error"] > 1e-3 for p in perturbed)),
        "why": ("Darmois matching is two conditions (areal radius and its "
                "normal derivative) on two constants (b, S), so it has a unique "
                "solution. Holding the radius right and moving b by 5% leaves a "
                "slope error of order 1e-2: there is no second matching pair, "
                "and therefore no tube area, neck radius or length left to "
                "choose. PR #263-#265 spent three rounds discovering that a "
                "chosen area was carrying the answer."),
    }


def check_p2() -> dict:
    """P2 — the Misner-Sharp mass parameter does not jump at the seam."""
    rows, worst = [], 0.0
    for R, a in CASES:
        g = fm.geometry(R, a)
        s = np.linspace(-g.half_length, g.half_length, 401)
        inside = fm.misner_sharp(s, g.neck_radius)
        outside = float(fm.misner_sharp_exterior(a, R))
        spread = float(np.max(inside) - np.min(inside))
        jump = abs(float(inside[-1]) - outside)
        rows.append({"R": R, "a": a, "mu_inside": float(inside[0]),
                     "mu_constant_to": spread, "mu_exterior_at_a": outside,
                     "jump": jump})
        worst = max(worst, spread, jump)
    return {"rows": rows, "worst": worst, "holds": bool(worst < 1e-12),
            "why": ("mu = f^2(1-f'^2) equals b^2 everywhere inside and "
                    "R^2 sin^4 chi outside; they agree at chi = a because "
                    "b = R sin^2 a. The seam is smooth exactly when the "
                    "quasi-local mass parameter is continuous -- the 5D lift "
                    "of PR #265's Hawking-mass matching.")}


def check_p3() -> dict:
    """P3 — no Israel surface layer, and the normal pressure agrees."""
    rows, worst = [], 0.0
    for R, a in CASES:
        j = fm.junction_jumps(R, a)
        rows.append({"R": R, "a": a,
                     "metric_jump": j["induced_metric_jump"],
                     "curvature_jump": j["extrinsic_curvature_jump"],
                     "p_normal_in": j["normal_pressure_inside"],
                     "p_normal_out": j["normal_pressure_outside"],
                     "p_jump": j["normal_pressure_jump"]})
        worst = max(worst, j["induced_metric_jump"],
                    j["extrinsic_curvature_jump"], j["normal_pressure_jump"])
    return {"rows": rows, "worst": worst, "holds": bool(worst < 1e-12),
            "c1_not_c2": fm.junction_jumps()["second_derivative_jumps"],
            "why": ("[h_ab] = [K_ab] = 0 gives S_ab = 0. The normal pressure "
                    "agreeing at -3/(8 pi G_5 R^2) on both sides is the "
                    "Gauss-Codazzi constraint that must hold when no shell is "
                    "present, and it is computed from the two sides separately.")}


def check_p4() -> dict:
    """P4 — every smooth traversable lapse pays the same neck NEC price.

    This is a **falsification attempt**, not a confirmation: it searches a
    deliberately hostile family of lapses -- asymmetric, rapidly oscillating,
    exponentially growing, nearly vanishing -- for any that evades the price.
    """
    g = fm.geometry()
    predicted = fm.null_energy_at_neck(g.neck_radius)
    families = {
        "N = 1 (ultrastatic)": lambda s, b: np.ones_like(s),
        "N = 1 + 0.7 s (asymmetric)": lambda s, b: 1.0 + 0.7 * s,
        "N = 1 + 3 s^2": lambda s, b: 1.0 + 3.0 * s ** 2,
        "N = 1 - 2 s^2 + 5 s^3": lambda s, b: 1.0 - 2.0 * s ** 2 + 5.0 * s ** 3,
        "N = exp(4 s)": lambda s, b: np.exp(4.0 * s),
        "N = 2 + cos(9 s)": lambda s, b: 2.0 + np.cos(9.0 * s),
        "N = 0.05 + 8 s^2 (nearly null)": lambda s, b: 0.05 + 8.0 * s ** 2,
    }
    rows, evaded = [], []
    for name, lapse in families.items():
        value = float(fm.stress_tensor(np.array([0.0]), g.neck_radius,
                                       lapse=lapse)["radial_nec"][0])
        rows.append({"lapse": name, "nec_at_neck": value,
                     "deviation": abs(value - predicted)})
        if value >= 0.0:
            evaded.append(name)
    # the vacuum escape: N(0) = 0
    vac = float(fm.lapse_vacuum(np.array([0.0]), g.neck_radius)[0])
    return {
        "predicted": predicted, "rows": rows,
        "worst_deviation": max(r["deviation"] for r in rows),
        "no_lapse_evades_it": bool(not evaded),
        "holds": bool(not evaded
                      and max(r["deviation"] for r in rows) < 1e-9),
        "vacuum_lapse_at_neck": vac,
        "why": ("The lapse enters p_s only through 3 f'N'/(fN), and f'(0) = 0 "
                "is what MAKES s = 0 a neck. So the term vanishes there "
                "whatever N' does. This is stronger than the proposal stated: "
                "it needs no reflection symmetry, which is why an asymmetric "
                "and an oscillating lapse give the identical value."),
        "the_escape": ("N(0) = 0 -- the Tangherlini horizon branch, which is "
                       "vacuum and non-traversable. Smooth AND traversable "
                       "implies radial NEC violation at the neck."),
    }


def check_p5() -> dict:
    """P5 — the closed-form admittance, against an independent BVP solve."""
    rows, worst_rel, worst_row = [], 0.0, 0.0
    for ell in (0, 1, 2, 3, 5):
        closed = fm.static_admittance(ell)
        numeric = fm.solve_admittance(ell, steps=4000)
        rel = float(np.max(np.abs(numeric - closed)) / np.max(np.abs(closed)))
        rows.append({"ell": ell, "closed_diag": float(closed[0, 0]),
                     "closed_off": float(closed[0, 1]),
                     "relative_error": rel})
        worst_rel = max(worst_rel, rel)
    # second-order convergence of the independent solver
    convergence = []
    prev = None
    closed = fm.static_admittance(2)
    for steps in (1000, 2000, 4000):
        err = float(np.max(np.abs(fm.solve_admittance(2, steps=steps) - closed))
                    / np.max(np.abs(closed)))
        convergence.append({"steps": steps, "relative_error": err,
                            "ratio": None if prev is None else prev / err})
        prev = err
    # monopole: row sums must vanish exactly, and G must match two ways
    g_closed = fm.monopole_conductance()
    y0 = fm.static_admittance(0)
    worst_row = float(np.max(np.abs(y0.sum(axis=1))))
    resistance = fm.static_resistance()
    return {
        "rows": rows, "worst_relative_error": worst_rel,
        "convergence": convergence,
        "monopole_conductance": g_closed,
        "monopole_from_resistance": 2.0 * math.pi ** 2 / resistance,
        "conductance_agreement": abs(g_closed - 2.0 * math.pi ** 2 / resistance),
        "static_resistance": resistance,
        "worst_monopole_row_sum": worst_row,
        "row_sums_vanish": bool(worst_row < 1e-12),
        "holds": bool(worst_rel < 1e-5 and worst_row < 1e-12
                      and abs(g_closed - 2.0 * math.pi ** 2 / resistance) < 1e-12
                      and convergence[-1]["ratio"] > 3.5),
        "why": ("The BVP solve never uses the sinh/cosh reduction, so this is a "
                "regression and not an identity. A shooting basis was tried "
                "first and rejected: its two solutions span e^{+-kx} over a "
                "rapidity 2kX ~ 23 at l = 5, so the far boundary condition "
                "costs ten digits and the error grew under refinement instead "
                "of shrinking. The conservative tridiagonal form converges at "
                "second order at every l tested."),
        "no_static_shunt": ("Row sums vanish exactly: a constant is in the "
                            "kernel, so there is no static monopole shunt "
                            "through the handle."),
    }


def check_p6() -> dict:
    """P6 — the vacuum control comes from the lapse alone."""
    g = fm.geometry()
    s = np.linspace(-g.half_length, g.half_length, 401)
    # identical spatial profile for both branches
    profile = fm.throat_radius(s, g.neck_radius)
    ricci = float(np.max(np.abs(fm.spatial_ricci_scalar(s, g.neck_radius))))
    n_vac = fm.lapse_vacuum(s, g.neck_radius)
    n_trav = fm.lapse_ultrastatic(s, g.neck_radius)
    # Tangherlini identification: F = 1 - b^2/r^2 with r^2 = s^2 + b^2
    r_sq = s * s + g.neck_radius ** 2
    tangherlini = 1.0 - g.neck_radius ** 2 / r_sq
    lapse_matches = float(np.max(np.abs(n_vac ** 2 - tangherlini)))
    # vacuum branch must have vanishing radial NEC away from the horizon
    interior = np.abs(s) > 0.15 * g.half_length
    vac_stress = fm.stress_tensor(s[interior], g.neck_radius,
                                  lapse=fm.lapse_vacuum)
    return {
        "spatial_ricci_max": ricci,
        "profile_is_shared": True,
        "lapse_vacuum_at_neck": float(n_vac[len(s) // 2]),
        "lapse_travers_at_neck": float(n_trav[len(s) // 2]),
        "tangherlini_lapse_residual": lapse_matches,
        "vacuum_radial_nec_max": float(np.max(np.abs(vac_stress["radial_nec"]))),
        "vacuum_is_stress_free": bool(
            np.max(np.abs(vac_stress["radial_nec"])) < 1e-6),
        "holds": bool(ricci < 1e-10 and lapse_matches < 1e-12
                      and abs(float(n_vac[len(s) // 2])) < 1e-15
                      and np.max(np.abs(vac_stress["radial_nec"])) < 1e-6),
        "why": ("N_vac^2 = |s|^2/(s^2+b^2) is exactly the Tangherlini F(r) = "
                "1 - b^2/r^2 with r^2 = s^2+b^2, on the SAME spatial profile "
                "that the ultrastatic branch uses. The repository's vacuum "
                "throat and transaction throat are one spatial geometry with "
                "two lapses; the entire physical fork is the number N(0)."),
    }


def run_probe() -> dict:
    checks: List[dict] = []
    for cid, name, fn in (
            ("P1", "*** the geometry is parameter-free (falsifier run) ***", check_p1),
            ("P2", "Misner-Sharp continuity across the seam", check_p2),
            ("P3", "no Israel shell, and the normal pressure agrees", check_p3),
            ("P4", "*** every smooth traversable lapse pays -3/b^2 (falsifier run) ***",
             check_p4),
            ("P5", "*** the closed-form admittance vs an independent BVP solve ***",
             check_p5),
            ("P6", "the vacuum control is the same metric with N(0) = 0", check_p6)):
        detail = fn()
        extra = detail.get("no_free_parameter_survives", True) and \
            detail.get("no_lapse_evades_it", True)
        checks.append({"id": cid, "name": name, "detail": detail,
                       "pass": bool(detail["holds"] and extra)})
    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    p1, p2, p3 = detail("P1"), detail("P2"), detail("P3")
    p4, p5, p6 = detail("P4"), detail("P5"), detail("P6")
    g = fm.geometry()

    L: List[str] = [
        "# The finite-mouth scalar-flat handle", "",
        f"**{summary['passed']}/{summary['total']} pre-registered predictions "
        "reproduced.**", "",
        "Every number below was frozen in `docs/finite_mouth_prereg.md` "
        "**before this module existed**. P1 and P4 are falsification attempts "
        "rather than confirmations.", "",
        "> ## The construction", "",
        "> One assumption: the closed `S³` universe is the totally geodesic "
        "equator of a round `S⁴_R` spatial bulk. It adds no scale and no shape "
        "parameter. Demanding `⁴R = 0` in the throat forces "
        "`f = √(s²+b²)`, and Darmois matching to the exterior then fixes "
        "**both** constants:", "",
        "> ```",
        "> b = R sin²a ,   S = R sin a cos a ,   L = R sin 2a",
        "> ```", "",
        "> There is no tube area, neck radius or throat length left to choose "
        "— which is exactly the freedom that carried the answer through "
        "PR #263–#265.", "",
        "---", "",
        "## P1 — the geometry is parameter-free", "",
        "| `R` | `a` | `b` | `S` | `L` | matching residual |",
        "|--|--|--|--|--|--|"]
    for r in p1["rows"]:
        L.append(f"| `{r['R']:g}` | `{r['a']:g}` | `{r['b']:.9f}` | "
                 f"`{r['S']:.9f}` | `{r['L']:.9f}` | "
                 f"`{max(r['match_radius'], r['match_slope']):.1e}` |")
    L += ["",
          "**Falsifier.** Holding the areal radius correct and moving `b` by "
          "±5% leaves a slope error of "
          f"`{p1['perturbations_break_matching'][0]['slope_error']:.2e}`: there "
          "is no second matching pair.", "",
          "> " + p1["why"], ""]

    L += ["## P2 — the quasi-local mass does not jump", "",
          "| `R` | `a` | `μ` inside | `μ_ext(a)` | jump |", "|--|--|--|--|--|"]
    for r in p2["rows"]:
        L.append(f"| `{r['R']:g}` | `{r['a']:g}` | `{r['mu_inside']:.9f}` | "
                 f"`{r['mu_exterior_at_a']:.9f}` | `{r['jump']:.1e}` |")
    L += ["", "> " + p2["why"], ""]

    L += ["## P3 — no Israel shell", "",
          "| `R` | `a` | `[h]` | `[K]` | `p_n` inside | `p_n` outside |",
          "|--|--|--|--|--|--|"]
    for r in p3["rows"]:
        L.append(f"| `{r['R']:g}` | `{r['a']:g}` | `{r['metric_jump']:.1e}` | "
                 f"`{r['curvature_jump']:.1e}` | `{r['p_normal_in']:.6f}` | "
                 f"`{r['p_normal_out']:.6f}` |")
    L += ["", "> " + p3["why"], "", "> " + p3["c1_not_c2"], ""]

    L += ["## P4 — the neck NEC price, attacked", "",
          f"Predicted `8πG₅(ρ+p_s)|₀ = −3/b² = {p4['predicted']:.9f}` for "
          "**every** smooth lapse with `N(0) > 0`. Seven hostile lapses:", "",
          "| lapse | `(ρ+p_s)` at the neck | deviation |", "|--|--|--|"]
    for r in p4["rows"]:
        L.append(f"| `{r['lapse']}` | `{r['nec_at_neck']:.9f}` | "
                 f"`{r['deviation']:.1e}` |")
    L += ["",
          f"**None evades it.** Worst deviation `{p4['worst_deviation']:.1e}`.", "",
          "> " + p4["why"], "", "> **" + p4["the_escape"] + "**", ""]

    L += ["## P5 — the admittance, against an independent solver", "",
          "| `ℓ` | closed-form diagonal | off-diagonal | relative error |",
          "|--|--|--|--|"]
    for r in p5["rows"]:
        L.append(f"| {r['ell']} | `{r['closed_diag']:.9f}` | "
                 f"`{r['closed_off']:.9f}` | `{r['relative_error']:.1e}` |")
    L += ["", "Second-order convergence of the independent solve at `ℓ = 2`:", "",
          "| steps | relative error | ratio |", "|--|--|--|"]
    for c in p5["convergence"]:
        ratio = "—" if c["ratio"] is None else f"`{c['ratio']:.2f}`"
        L.append(f"| `{c['steps']}` | `{c['relative_error']:.2e}` | {ratio} |")
    L += ["",
          f"Monopole: `G = {p5['monopole_conductance']:.9f}` = `π²R²sin⁴a/cos a`, "
          f"and `2π²/I₃` agrees to `{p5['conductance_agreement']:.1e}` with "
          f"`I₃ = {p5['static_resistance']:.9f}`.", "",
          f"Row sums vanish to `{p5['worst_monopole_row_sum']:.1e}`. "
          + p5["no_static_shunt"], "",
          "> " + p5["why"], ""]

    L += ["## P6 — one spatial geometry, two lapses", "",
          "| branch | `N(0)` | stress | causal character |", "|--|--|--|--|",
          f"| Tangherlini | `{p6['lapse_vacuum_at_neck']:.1e}` | vacuum "
          f"(radial NEC `{p6['vacuum_radial_nec_max']:.1e}`) | horizon, "
          "non-traversable |",
          f"| ultrastatic | `{p6['lapse_travers_at_neck']:.1f}` | anisotropic, "
          "NEC-violating | traversable |", "",
          f"`⁴R` vanishes on the shared profile to "
          f"`{p6['spatial_ricci_max']:.1e}`, and `N_vac²` reproduces the "
          f"Tangherlini `F(r)` to `{p6['tangherlini_lapse_residual']:.1e}`.", "",
          "> " + p6["why"], "",
          "---", "",
          "## What this does not settle", "",
          "The ANEC cost of the finite ultrastatic throat is "
          f"`{fm.null_energy_integral(1.0, 0.30):.6f}/(8πG₅)` at `R = 1, "
          "a = 0.3`, diverging as `−3π/(2Ra²)` for small mouths — so the "
          "point-mouth limit is singularly expensive. Whether any classical "
          "BAM degree of freedom can supply that stress, and therefore whether "
          "`N(0) > 0` is available at all, is untouched here. If none can, the "
          "geometry collapses onto the Tangherlini branch and the MTY "
          "transaction mechanism is unavailable.", "",
          "The discrete BAM identification is deliberately absent: "
          "`Φ_spatial`, `(−1)^ℓ`, `η_orientation`, `η_wrap` and `U_spin` "
          "remain five separate objects, none folded into the metric.", ""]
    return "\n".join(L)


def main() -> int:
    summary = run_probe()
    text = render_markdown(summary)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs",
                          f"{stamp}_finite_mouth_probe")
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
