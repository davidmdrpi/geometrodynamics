"""
Draw the wave as vectors normal to the surface, and they intersect

> Framing: still a REPRESENTATION of a scalar field on a fixed round S².  What
> changes is which object is drawn — the vectors, not their tips.

THE ABSTRACTION THAT WAS MISSING
────────────────────────────────
Every earlier round drew the field as the tip of a radial displacement, the
graph r = f(σ).  A graph is embedded by construction, which is exactly why it
could not wind (circle_slice), could not self-intersect at any gap or gain
(slice_folding), and needed an invented tangential freedom before anything
would fold.

Draw the VECTORS instead of their tips and the obstruction is gone, for a
classical reason: neighbouring normals to a curve meet at its centre of
curvature.  The normal family has an envelope — the evolute — and a normal of
length L crosses its neighbours as soon as L exceeds the local radius of
curvature ρ = 1/κ.  Nothing is added by hand.

WHAT IS CHECKED
───────────────
T2  THE SAME WAVE, DRAWN TWO WAYS.  At one instant the graph of the tips gives
    0 self-intersections and the normal field gives hundreds.  Nothing about
    the field changed; only which object was drawn.

T3  THE THRESHOLD IS THE RADIUS OF CURVATURE.  ρ_min falls 0.0755 → 0.0338 as
    the ring converges — the wave sharpens its own surface — and the first
    drawn crossing falls with it, 0.452 → 0.182.  The envelope bounds the
    drawn crossing from below at every time, as it must: a finite sampling
    stride can only ever lag the continuous envelope.

T4  THE RESET IS A SECOND, SEPARABLE MECHANISM.  A normal leaving through
    R_outer re-enters at R_inner at the angle where it left, and that stub
    shoots outward from deep inside the annulus across vectors it could never
    otherwise reach.  At the focus with L = 0.35: 304 crossings among the
    normals, 402 with the reset.

T5  AND THE GAP MATTERS AGAIN.  slice_folding established that the fold
    threshold was independent of the gap, so shrinking the vacuole could never
    buy an intersection.  For normals the vector length IS what spans the gap,
    so the two are one knob — and at the tightest gap the reset dominates
    completely: 0 crossings from the normals alone, 522 once they re-enter.

SCOPE
─────
The vector length L is a display choice like every gain in this series; the
directions and the curvature are the surface's own.  What is derived is that
the crossing threshold is the radius of curvature, that the focusing drives it
down, and that the reset is an independent mechanism on top.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.normal_field import (
    measure_normals_intersect_where_a_graph_cannot,
    measure_the_gap_matters_again,
    measure_the_reset_adds_intersections,
    measure_the_threshold_is_the_radius_of_curvature,
)

FOCUS_TIME = 3.0
GAIN = 0.30


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "claim": ("draw the field as vectors normal to the surface rather than "
                  "as the graph of their tips, and ask again whether anything "
                  "intersects"),
        "why_it_can": "neighbouring normals meet at the centre of curvature",
        "threshold": "L > ρ = 1/κ",
        "pass": True,
    }


def t2_the_same_wave_drawn_two_ways() -> dict:
    r = measure_normals_intersect_where_a_graph_cannot(t=FOCUS_TIME, gain=GAIN)
    return {"name": "T2_the_same_wave_drawn_two_ways", **r,
            "pass": bool(r["the_graph_never_crosses"] and r["the_normals_do"])}


def t3_the_threshold_is_the_radius_of_curvature() -> dict:
    r = measure_the_threshold_is_the_radius_of_curvature(gain=GAIN)
    return {"name": "T3_the_threshold_is_the_radius_of_curvature", **r,
            "pass": bool(r["the_focus_sharpens_the_surface"]
                         and r["the_envelope_bounds_the_drawn_crossing"])}


def t4_the_reset_is_a_second_mechanism() -> dict:
    r = measure_the_reset_adds_intersections(t=FOCUS_TIME, gain=GAIN)
    return {"name": "T4_the_reset_is_a_second_mechanism", **r,
            "pass": bool(r["the_reset_adds_crossings"]
                         and r["it_grows_with_the_stub"])}


def t5_the_gap_matters_again() -> dict:
    r = measure_the_gap_matters_again(t=FOCUS_TIME, gain=GAIN)
    return {"name": "T5_the_gap_matters_again", **r,
            "pass": bool(r["the_gap_changes_the_count"]
                         and r["every_gap_crosses"])}


def t6_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T6_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_the_same_wave_drawn_two_ways(),
             t3_the_threshold_is_the_radius_of_curvature(),
             t4_the_reset_is_a_second_mechanism(),
             t5_the_gap_matters_again()]
    tests.append(t6_assessment(tests))
    t2, t3, t4, t5 = tests[1], tests[2], tests[3], tests[4]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_NORMALS_INTERSECT"
        verdict = (
            "THE NORMALS INTERSECT. Four rounds of this series established that "
            "a height field cannot wind and cannot cross itself. All four were "
            "about the same object — the graph of the displacement's tips — and "
            "drawing the displacement itself removes the obstruction "
            "entirely.\n\n"
            "THE SAME WAVE, DRAWN TWO WAYS. At one instant the graph of the "
            f"tips gives {t2['graph_self_intersections']} self-intersections "
            f"and the normal field gives {t2['most_normal_crossings']}. Nothing "
            "about the field changed and nothing was added by hand; only which "
            "object was drawn. The reason is classical: neighbouring normals to "
            "a curve meet at its centre of curvature, so the normal family has "
            "an envelope and a normal longer than the local radius of curvature "
            "has already crossed its neighbours.\n\n"
            "THE THRESHOLD IS THE RADIUS OF CURVATURE, AND THE WAVE DRIVES IT "
            f"DOWN. ρ_min falls from {t3['rho_at_the_start']:.4f} in mid-flight "
            f"to {t3['rho_at_the_focus']:.4f} at the focus, a factor of "
            f"{t3['sharpening_factor']:.2f} — the converging ring sharpens its "
            "own surface — and the first drawn crossing falls with it, "
            f"{t3['rows'][0]['first_drawn_crossing']:.3f} to "
            f"{t3['rows'][-1]['first_drawn_crossing']:.3f}. The envelope bounds "
            "the drawn crossing from below at every time, as it must: a finite "
            "sampling stride can only lag the continuous envelope, never "
            "precede it.\n\n"
            "AND THE RESET IS A SECOND MECHANISM ON TOP. A normal long enough "
            "to leave through R_outer re-enters at R_inner at the angle where "
            "it left, and that stub shoots outward from deep inside the annulus "
            "across vectors it could never otherwise have reached. At the focus "
            f"the count goes from {t4['rows'][2]['without_reset']} to "
            f"{t4['rows'][2]['with_reset']} at L = 0.35, and the addition grows "
            f"with the stub, reaching {t4['most_added']} at the longest length "
            "tested. The two mechanisms are separable and both are real.\n\n"
            "AND THE GAP MATTERS AGAIN. This is what the height representation "
            "had severed: slice_folding showed the fold threshold did not know "
            "the gap existed, so shrinking the vacuole could never buy an "
            "intersection. For normals the vector length IS what spans the gap, "
            "so they are one knob — and at the tightest gap the reset dominates "
            f"completely: {t5['rows'][-1]['crossings']} crossings from the "
            f"normals alone against {t5['rows'][-1]['with_reset']} once they "
            "re-enter. Reducing the distance between the shells now produces "
            "intersections rather than being unable to.\n\n"
            "SCOPE. The vector length is a display choice like every gain in "
            "this series; the directions and the curvature are the surface's "
            "own. What is derived is that the crossing threshold is the radius "
            "of curvature, that the focusing drives it down, and that the reset "
            "adds an independent mechanism."
        )
    else:
        verdict_class = "NORMAL_FIELD_INCONCLUSIVE"
        verdict = ("INCONCLUSIVE. A check failed; review the graph comparison, "
                   "the curvature threshold, the reset, or the gap sweep.")

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "the slice's normal vectors drawn as vectors, and where they cross"
        ),
        "the_object": "the displacement itself, not the graph of its tips",
        "the_threshold": "the radius of curvature of the deformed surface",
        "the_second_mechanism": "the stub that re-enters at R_inner",
        "geometry": {"focus_time": FOCUS_TIME, "gain": GAIN},
        "tests": tests,
        "n_passed": tests[-1]["n_passed"],
        "n_total": tests[-1]["n_total"],
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    t = {x["name"]: x for x in s["tests"]}
    out: List[str] = []
    out.append("# Draw the wave as vectors, and they intersect\n")
    out.append(f"_{s['timestamp_utc']}_\n")

    c = t["T2_the_same_wave_drawn_two_ways"]
    out.append("## The same wave, drawn two ways\n")
    out.append(f"Graph of the tips: **{c['graph_self_intersections']}** "
               "self-intersections.\n")
    out.append("| vector length | normal crossings |")
    out.append("|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['length']:.2f} | {r['normal_crossings']} |")
    out.append("")

    c = t["T3_the_threshold_is_the_radius_of_curvature"]
    out.append("## The threshold is the radius of curvature\n")
    out.append("| t | `ρ_min` | first drawn crossing |")
    out.append("|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['t']:.1f} | {r['rho_min']:.4f} | "
                   f"{r['first_drawn_crossing']:.4f} |")
    out.append(f"\nThe focus sharpens the surface by "
               f"{c['sharpening_factor']:.2f}×.\n")

    c = t["T4_the_reset_is_a_second_mechanism"]
    out.append("## The reset adds more\n")
    out.append("| L | wrapped | without reset | with reset | added |")
    out.append("|---:|---:|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['length']:.2f} | {r['wrapped']} | "
                   f"{r['without_reset']} | {r['with_reset']} | "
                   f"**{r['added']}** |")
    out.append("")

    c = t["T5_the_gap_matters_again"]
    out.append("## And the gap matters again\n")
    out.append("| δ | L = δ | normals alone | with reset |")
    out.append("|---:|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['delta']:.2f} | {r['length']:.2f} | "
                   f"{r['crossings']} | **{r['with_reset']}** |")
    out.append("")

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
    out = here / "runs" / f"{ts}_normal_field_probe"
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
