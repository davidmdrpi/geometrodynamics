"""
The restored geometry: one continuous S², warped by the wave it carries

> Framing: a *display* of a solved classical field as a radial displacement of
> a fixed closed surface.  Nothing here is backreaction — the wave does not
> feel the geometry it is deforming.

WHAT THIS IS FOR
────────────────
The archive rendered the BAM picture as ONE CONTINUOUS SURFACE whose radius
carried the field, nested between two fixed shells like Russian dolls.  The
recent work replaced that with maps, strips and meridional sections: correct,
measurable, and no longer the thing you could look at.  This probe pins the
restored object so the picture cannot quietly drift again.

    r(θ, φ, t) = R_mid + Δ · tanh( g · u(θ, φ, t) / u_ref )

    R_inner = 0.74  ·······  the inner doll
    R_mid   = 1.00  ·······  the surface the wave lives on
    R_outer = 1.26  ·······  the outer doll

The radii are deliberately the vacuole of `ring_caustic_defect_probe`, so the
shell that probe's ring caustic lands on is the inner doll drawn here.

WHAT IS CHECKED
───────────────
T2  ONE CLOSED SURFACE.  Poles single-valued, φ seam matched to machine
    precision, and the mesh a genuine topological sphere — nothing cut out.
    This is what makes "a pulse fills its own void" a statement about a closed
    manifold rather than about a patch.

T3  IT NEVER TOUCHES A DOLL.  Over a full return time the warped radius stays
    strictly inside (R_inner, R_outer), and the margins are reported.

T4  THE DISPLAY IS HONEST.  tanh is strictly increasing, so the displacement
    preserves the sign of the field everywhere and the ORDERING of every pair
    of amplitudes.  It does not preserve ratios, and the probe says so rather
    than pretending otherwise.  Raising the gain cannot flip a sign.

T5  THE DEEPEST WARP IS AT THE ANTIPODE.  Measured, not imposed: the largest
    displacement over the run occurs at geodesic distance π from the source,
    at t slightly under π.  The shortfall is the pulse's own width — it
    shrinks as the pulse narrows, which the probe demonstrates by scanning.

T6  THE FOCUS RINGS.  The arrival is not a single mound: the surface is driven
    outward, then inverts and is pulled inward past the mid radius.  Both
    excursions are reported, and it is the INWARD one — toward the inner doll,
    toward a throat — that the geometry ends up in.

SCOPE
─────
A display, on a fixed background.  The wave is linear and does not act on the
surface; the radial displacement is a rendering of u, chosen bounded so the
picture stays inside the vacuole.  Nothing here shows a throat forming, and
nothing here is dynamical geometry.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.warped_sphere import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    NestedShells,
    WarpedSphere,
    measure_containment,
    measure_focus,
)

N_THETA, N_PHI = 61, 91
PULSE_WIDTH = 0.18
GAIN = 1.6


def _surface(pulse_width: float = PULSE_WIDTH, gain: float = GAIN) -> WarpedSphere:
    return WarpedSphere(n_theta=N_THETA, n_phi=N_PHI,
                        pulse_width=pulse_width, gain=gain)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    sh = NestedShells()
    return {
        "name": "T1_goal",
        "claim": ("restore the archive's geometry — one continuous S² whose "
                  "radius carries a SOLVED wave, nested between two shells"),
        "r_inner": sh.r_inner,
        "r_mid": sh.r_mid,
        "r_outer": sh.r_outer,
        "half_gap": sh.delta,
        "same_vacuole_as_ring_caustic_probe": True,
        "pass": True,
    }


def t2_one_closed_surface() -> dict:
    s = _surface()
    s.advance_to(1.0)
    X, Y, Z = s.mesh()
    seam = float(np.max(np.abs(np.stack([X[:, 0] - X[:, -1],
                                         Y[:, 0] - Y[:, -1],
                                         Z[:, 0] - Z[:, -1]]))))
    pole_spread = max(float(np.ptp(Z[0, :])), float(np.ptp(Z[-1, :])))
    finite = bool(np.all(np.isfinite(X)) and np.all(np.isfinite(Y))
                  and np.all(np.isfinite(Z)))
    return {
        "name": "T2_one_closed_surface",
        "seam_mismatch": seam,
        "pole_spread": pole_spread,
        "carries_poles": bool(s.theta[0] == 0.0 and s.theta[-1] == math.pi),
        "all_finite": finite,
        "is_closed": bool(s.is_closed()),
        "pass": bool(s.is_closed() and finite and seam < 1e-12
                     and pole_spread < 1e-12),
    }


def t3_never_touches_a_doll() -> dict:
    s = _surface()
    r = measure_containment(s, t_end=RETURN_TIME, frames=160)
    return {
        "name": "T3_never_touches_a_doll",
        **{k: r[k] for k in ("r_min", "r_max", "r_inner", "r_outer",
                             "contained", "closest_approach_inner",
                             "closest_approach_outer")},
        "pass": bool(r["contained"]),
    }


def t4_display_is_order_preserving() -> dict:
    s = _surface()
    s.advance_to(1.3)
    u = s.field().ravel()
    d = s.displacement().ravel()
    order = np.argsort(u)
    worst_inversion = float(np.min(np.diff(d[order]))) if u.size > 1 else 0.0
    sign_ok = bool(np.all(np.sign(d) == np.sign(u)))

    lo, hi = _surface(gain=0.5), _surface(gain=3.0)
    for x in (lo, hi):
        x.advance_to(0.9)
    same_sign = bool(np.all(np.sign(lo.displacement())
                            == np.sign(hi.displacement())))
    # ratios are NOT preserved — state the distortion rather than hide it
    nz = np.abs(u) > 1e-9
    ratio_spread = float(np.ptp((np.abs(d[nz]) / np.abs(u[nz]))))
    return {
        "name": "T4_display_is_order_preserving",
        "sign_preserved": sign_ok,
        "worst_order_inversion": worst_inversion,
        "gain_cannot_flip_a_sign": same_sign,
        "ratio_distortion_spread": ratio_spread,
        "note": ("tanh is strictly increasing: signs and orderings survive, "
                 "ratios do not — the spread above is that distortion"),
        "pass": bool(sign_ok and same_sign and worst_inversion >= -1e-15),
    }


def t5_deepest_warp_is_at_the_antipode() -> dict:
    r = measure_focus(_surface(), t_end=1.15 * ANTIPODAL_TIME, frames=200)
    scan = []
    for w in (0.24, 0.18, 0.12, 0.08):
        f = measure_focus(_surface(pulse_width=w), t_end=1.15 * ANTIPODAL_TIME,
                          frames=200)
        scan.append({"pulse_width": w, "peak_time": f["peak_time"],
                     "time_error": f["time_error"]})
    monotone = all(scan[i]["time_error"] >= scan[i + 1]["time_error"] - 1e-12
                   for i in range(len(scan) - 1))
    return {
        "name": "T5_deepest_warp_is_at_the_antipode",
        "peak_distance": r["peak_distance"],
        "distance_error": r["distance_error"],
        "peak_time": r["peak_time"],
        "antipodal_time": ANTIPODAL_TIME,
        "time_error": r["time_error"],
        "pulse_width_scan": scan,
        "bias_shrinks_with_pulse_width": monotone,
        "pass": bool(r["distance_error"] < 1e-9 and r["time_error"] < 0.25
                     and monotone),
    }


def t6_the_focus_rings() -> dict:
    s = _surface()
    out_peak = {"fraction": -math.inf, "t": 0.0}
    in_peak = {"fraction": -math.inf, "t": 0.0}
    frames = 260
    for i in range(frames):
        t = (i + 1) * 1.15 * ANTIPODAL_TIME / frames
        s.advance_to(t)
        ex = s.excursion()
        if ex["outward_fraction"] > out_peak["fraction"] and t > 0.5 * ANTIPODAL_TIME:
            out_peak = {"fraction": ex["outward_fraction"], "t": s.t}
        if ex["inward_fraction"] > in_peak["fraction"] and t > 0.5 * ANTIPODAL_TIME:
            in_peak = {"fraction": ex["inward_fraction"], "t": s.t}
    return {
        "name": "T6_the_focus_rings",
        "outward_peak_fraction": out_peak["fraction"],
        "outward_peak_time": out_peak["t"],
        "inward_peak_fraction": in_peak["fraction"],
        "inward_peak_time": in_peak["t"],
        "inverts_after_the_mound": bool(in_peak["t"] > out_peak["t"]),
        "pass": bool(out_peak["fraction"] > 0.4 and in_peak["fraction"] > 0.4
                     and in_peak["t"] > out_peak["t"]),
    }


def t7_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T7_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_one_closed_surface(), t3_never_touches_a_doll(),
             t4_display_is_order_preserving(),
             t5_deepest_warp_is_at_the_antipode(), t6_the_focus_rings()]
    tests.append(t7_assessment(tests))
    t2, t3, t4, t5, t6 = tests[1], tests[2], tests[3], tests[4], tests[5]
    sh = NestedShells()

    if all(t["pass"] for t in tests):
        verdict_class = "THE_SURFACE_ITSELF_CARRIES_THE_WAVE"
        verdict = (
            "THE SURFACE ITSELF CARRIES THE WAVE. The archive's geometry is "
            "back, and this time the warp is solved rather than drawn.\n\n"
            "ONE CLOSED SURFACE. The mesh carries its own poles and its φ seam "
            f"matches to {t2['seam_mismatch']:.1e}, so it is a single "
            "manifold with nothing cut out of it. That is what makes 'a pulse "
            "sweeps every point once and fills its own void' a statement "
            "about a closed surface rather than about a patch — and it is "
            "exactly why a ring, not a pulse, is what a throat needs.\n\n"
            "NESTED, NEVER TOUCHING. Over a full return the radius stays "
            f"within [{t3['r_min']:.4f}, {t3['r_max']:.4f}], clearing the "
            f"inner doll by {t3['closest_approach_inner']:.4f} and the outer "
            f"by {t3['closest_approach_outer']:.4f}. The bound is structural: "
            "tanh cannot leave the gap.\n\n"
            "AN HONEST DISPLAY. The displacement preserves the sign of the "
            "field everywhere and the ordering of every pair of amplitudes "
            f"(worst inversion {t4['worst_order_inversion']:.1e}), and no "
            "gain can flip a sign. It does NOT preserve ratios — the "
            f"distortion spans {t4['ratio_distortion_spread']:.3f} across the "
            "surface — which is the price of keeping the picture inside the "
            "vacuole, and is stated rather than hidden.\n\n"
            "THE FOCUS IS MEASURED, NOT KEYED TO THE CLOCK. The deepest "
            "deformation over the run sits at geodesic distance "
            f"{t5['peak_distance']:.6f} from the source — the antipode, to "
            f"{t5['distance_error']:.1e} — at t = {t5['peak_time']:.4f} "
            f"against π = {ANTIPODAL_TIME:.4f}. The shortfall of "
            f"{t5['time_error']:.4f} is the pulse's own width, not solver "
            "error: narrowing the pulse from 0.24 to 0.08 shrinks it "
            "monotonically. This is the difference from the archive, which "
            "grew its mound on a growth function tied to the frame number.\n\n"
            "AND IT RINGS. The arrival is not one mound. The surface is "
            f"driven {100 * t6['outward_peak_fraction']:.1f}% of the way to "
            f"the outer doll at t = {t6['outward_peak_time']:.3f}, then "
            f"inverts and is pulled {100 * t6['inward_peak_fraction']:.1f}% "
            f"of the way to the inner doll at t = {t6['inward_peak_time']:.3f}. "
            "The focus ends up pulling the geometry INWARD — toward the inner "
            "shell, which is the one the ring caustic lands on.\n\n"
            "SCOPE. This is a display of a solved field as a radial "
            "displacement of a FIXED surface, not backreaction: the wave does "
            "not feel what it is deforming, so no throat forms here and "
            "nothing about nucleation follows. What it restores is the object "
            "the intuition was built on — and it now moves because a solver "
            "says so."
        )
    else:
        verdict_class = "WARPED_SPHERE_INCONCLUSIVE"
        verdict = ("INCONCLUSIVE. A check failed; review closure of the mesh, "
                   "containment between the shells, monotonicity of the "
                   "display, the antipodal focus, or its inversion.")

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "one continuous S² whose radius carries a solved wave, nested "
            "between the two shells of the ring-caustic vacuole"
        ),
        "the_surface": "closed — its own poles, its seam matched, nothing cut out",
        "the_dolls": "never touched; the display is bounded by construction",
        "the_focus": "deepest warp at the antipode, at t just under π by the pulse width",
        "the_inversion": "mound first, then pulled inward past R_mid",
        "geometry": {
            "r_inner": sh.r_inner, "r_mid": sh.r_mid, "r_outer": sh.r_outer,
            "half_gap": sh.delta, "gain": GAIN, "pulse_width": PULSE_WIDTH,
            "mesh": [N_THETA, N_PHI],
            "antipodal_time": ANTIPODAL_TIME, "return_time": RETURN_TIME,
        },
        "tests": tests,
        "n_passed": tests[-1]["n_passed"],
        "n_total": tests[-1]["n_total"],
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    g = s["geometry"]
    t = {x["name"]: x for x in s["tests"]}
    out: List[str] = []
    out.append("# The restored geometry: one continuous S², warped by the wave "
               "it carries\n")
    out.append(f"_{s['timestamp_utc']}_\n")
    out.append("> A **display** of a solved classical field as a radial "
               "displacement of a fixed closed surface. Not backreaction.\n")

    out.append("## The dolls\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| `R_inner` | {g['r_inner']:.4f} |")
    out.append(f"| `R_mid` | {g['r_mid']:.4f} |")
    out.append(f"| `R_outer` | {g['r_outer']:.4f} |")
    out.append(f"| half-gap `Δ` | {g['half_gap']:.4f} |")
    out.append(f"| display gain `g` | {g['gain']:.2f} |")
    out.append(f"| pulse width | {g['pulse_width']:.2f} |")
    out.append(f"| mesh | {g['mesh'][0]} × {g['mesh'][1]} |\n")

    c = t["T2_one_closed_surface"]
    out.append("## One closed surface\n")
    out.append("| check | value |")
    out.append("|---|---:|")
    out.append(f"| seam mismatch | `{c['seam_mismatch']:.1e}` |")
    out.append(f"| pole spread | `{c['pole_spread']:.1e}` |")
    out.append(f"| carries the poles | {c['carries_poles']} |\n")

    c = t["T3_never_touches_a_doll"]
    out.append("## Nested, never touching\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| `r_min` over the run | {c['r_min']:.4f} |")
    out.append(f"| `r_max` over the run | {c['r_max']:.4f} |")
    out.append(f"| clearance, inner | {c['closest_approach_inner']:.4f} |")
    out.append(f"| clearance, outer | {c['closest_approach_outer']:.4f} |\n")

    c = t["T5_deepest_warp_is_at_the_antipode"]
    out.append("## The focus is measured\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| deepest warp at distance | {c['peak_distance']:.6f} |")
    out.append(f"| distance error from π | `{c['distance_error']:.1e}` |")
    out.append(f"| at time | {c['peak_time']:.4f} |")
    out.append(f"| against π | {c['antipodal_time']:.4f} |")
    out.append(f"| shortfall | {c['time_error']:.4f} |\n")
    out.append("The shortfall is the pulse's own width, not solver error:\n")
    out.append("| pulse width | peak time | shortfall |")
    out.append("|---:|---:|---:|")
    for r in c["pulse_width_scan"]:
        out.append(f"| {r['pulse_width']:.2f} | {r['peak_time']:.4f} | "
                   f"{r['time_error']:.4f} |")
    out.append("")

    c = t["T6_the_focus_rings"]
    out.append("## And it rings\n")
    out.append("| | fraction of the gap | at time |")
    out.append("|---|---:|---:|")
    out.append(f"| driven out toward `R_outer` | "
               f"{100 * c['outward_peak_fraction']:.1f}% | "
               f"{c['outward_peak_time']:.3f} |")
    out.append(f"| pulled in toward `R_inner` | "
               f"{100 * c['inward_peak_fraction']:.1f}% | "
               f"{c['inward_peak_time']:.3f} |\n")

    out.append("## Verdict\n")
    out.append(f"**{s['n_passed']}/{s['n_total']} checks passed.**\n")
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
    return "\n".join(out)


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_warped_sphere_geometry_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
