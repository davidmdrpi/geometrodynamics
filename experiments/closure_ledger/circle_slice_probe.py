"""
A circle slice of the scalar wave, and a bulk it can wrap through

> Framing: a REPRESENTATION of the scalar field on a fixed round S², cut by one
> great circle and drawn inside a glued annulus — not a derived boundary
> condition, and not backreaction.

THE SLICE
─────────
Cut the sphere with the great circle through the source and its antipode.
Parametrised by σ ∈ [−π, π), the geodesic distance from the source along that
circle is d = |σ|, so ONE circle carries BOTH halves of the wave: two lobes
running in opposite directions, meeting head-on at σ = ±π.

Nothing is re-solved.  The field is the 2-D solve sampled at d(σ), so the slice
inherits the real caustic rather than the 2× superposition that a 1-D wave on a
circle would give.  T2 checks that identification through the full (θ, φ)
machinery rather than assuming it.

THE BULK, AND ITS CROSSING RULE
───────────────────────────────
The slice lives in the vacuole: the annulus between R_inner and R_outer.  The
crossing rule glues them — a radius past R_outer re-enters at R_inner — so the
wave that reaches up into the bulk comes back INSIDE the circle.  The radial
coordinate becomes periodic with period gap = R_outer − R_inner, and the space
the curve lives on is a torus S¹_σ × S¹_ρ.

WHAT IS CHECKED
───────────────
T2  THE SLICE IS THE SPHERE.  The field at σ equals the sphere's field at
    geodesic distance |σ| reached through (θ, φ), and both arms are identical.

T3  THE LOBES MEET AT THE ANTIPODE, on time, symmetrically.

T4  THE WRAP THRESHOLD IS THE HALF-GAP OVER THE PEAK.  ε_crit = (R_outer −
    R_mid)/max|u|, predicted and then found by bisection.

T5  THE CURVE CLOSES.  σ = −π and σ = +π are the same point, wrapped or not —
    which is what makes the winding number a well-posed integer at all.

T6  THE SEAM IS CROSSED IN PAIRS.  The load-bearing NEGATIVE result: amplitude
    buys unsigned crossings (0 → 2 → 4 → 6), and the signed total is exactly
    zero at every gain and every time.  The curve is a graph r = f(σ) with f
    single-valued, so its radial winding number is identically zero.  A height
    field cannot wind.  The crossing ledger and an independently-computed
    degree are cross-checked against each other.

T7  DIFFERENT WAVES DIFFER IN ARC, NOT IN WINDING.  Driven at one common gain,
    pulses from 0.36 to 0.08 cross the seam the same number of times but put
    0.155 → 0.033 of the circle on the far sheet — and that arc is 2.61 × the
    pulse width for all of them.  None of them winds.

SCOPE
─────
The crossing rule is a representation choice, not a derived boundary condition:
nothing here makes the wave dynamically aware of the seam.  The zero winding
number is a fact about GRAPHS, so it constrains this construction and any other
that draws the bulk excursion as a height — which is the useful part.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.circle_slice import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    CircleSlice,
    measure_different_waves_wrap_differently,
    measure_the_curve_closes,
    measure_the_lobes_meet_at_the_antipode,
    measure_the_seam_is_crossed_in_pairs,
    measure_the_slice_is_the_sphere,
    measure_the_wrap_threshold,
)

PULSE_WIDTH = 0.18
N_SIGMA = 1441
N_RADIAL = 900
COMMON_GAIN = 0.80
WIDTHS = (0.36, 0.24, 0.14, 0.08)


def _slice() -> CircleSlice:
    return CircleSlice(pulse_width=PULSE_WIDTH, n_sigma=N_SIGMA,
                       n_radial=N_RADIAL)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "claim": ("cut the scalar wave with one great circle, and give it a "
                  "bulk whose outer boundary is glued to its inner one"),
        "slice": "d = |σ| along the circle through source and antipode",
        "crossing_rule": "a radius past R_outer re-enters at R_inner",
        "resulting_space": "a torus S¹_σ × S¹_ρ",
        "pass": True,
    }


def t2_the_slice_is_the_sphere() -> dict:
    r = measure_the_slice_is_the_sphere(_slice())
    return {
        "name": "T2_the_slice_is_the_sphere",
        **r,
        "pass": bool(r["is_a_cut_of_the_sphere"] and r["both_arms_are_the_same"]),
    }


def t3_the_lobes_meet_at_the_antipode() -> dict:
    r = measure_the_lobes_meet_at_the_antipode(_slice())
    return {
        "name": "T3_the_lobes_meet_at_the_antipode",
        **r,
        "pass": bool(r["lobes_are_symmetric"] and r["arrives_on_time"]),
    }


def t4_the_wrap_threshold() -> dict:
    r = measure_the_wrap_threshold(_slice())
    return {
        "name": "T4_the_wrap_threshold",
        **r,
        "pass": bool(r["threshold_is_the_half_gap_over_the_peak"]),
    }


def t5_the_curve_closes() -> dict:
    r = measure_the_curve_closes(_slice())
    return {"name": "T5_the_curve_closes", **r,
            "pass": bool(r["closes_exactly"])}


def t6_the_seam_is_crossed_in_pairs() -> dict:
    r = measure_the_seam_is_crossed_in_pairs(_slice())
    return {
        "name": "T6_the_seam_is_crossed_in_pairs",
        **r,
        "pass": bool(r["amplitude_buys_crossings"]
                     and r["the_signed_count_is_always_zero"]
                     and r["a_height_field_cannot_wind"]),
    }


def t7_different_waves_differ_in_arc() -> dict:
    r = measure_different_waves_wrap_differently(widths=WIDTHS,
                                                 gain=COMMON_GAIN)
    return {
        "name": "T7_different_waves_differ_in_arc",
        **r,
        "pass": bool(r["they_differ_in_arc_not_in_count"]
                     and r["the_far_sheet_arc_is_the_pulse"]
                     and r["none_of_them_winds"]),
    }


def t8_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T8_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_the_slice_is_the_sphere(),
             t3_the_lobes_meet_at_the_antipode(), t4_the_wrap_threshold(),
             t5_the_curve_closes(), t6_the_seam_is_crossed_in_pairs(),
             t7_different_waves_differ_in_arc()]
    tests.append(t8_assessment(tests))
    t2, t3, t4, t6, t7 = tests[1], tests[2], tests[3], tests[5], tests[6]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_SEAM_IS_CROSSED_IN_PAIRS"
        verdict = (
            "THE SEAM IS CROSSED IN PAIRS. Gluing the vacuole's outer boundary "
            "to its inner one lets the wave that reaches into the bulk come "
            "back inside the circle, and makes the space the curve lives on a "
            "torus — but it does not give the curve anywhere to wind.\n\n"
            "THE SLICE IS THE SPHERE, NOT A NEW PROBLEM. The field on the "
            "circle matches the sphere's field at geodesic distance |σ|, "
            f"reached through the full (θ, φ) route, to "
            f"{t2['worst_mismatch_against_the_sphere']:.1e}, and the two arms "
            f"agree to {t2['worst_mirror_asymmetry']:.1e}. So the slice "
            "inherits the real caustic rather than the 2× superposition a 1-D "
            "wave on a circle would give. The lobes leave in opposite "
            f"directions and arrive together at t = {t3['arrival_time']:.3f}, "
            f"symmetric to {abs(t3['final_left_lobe'] + t3['final_right_lobe']):.1e}.\n\n"
            "THE WRAP THRESHOLD IS EXACTLY THE HALF-GAP OVER THE PEAK. "
            f"ε_crit = {t4['predicted_threshold']:.6f} predicted, "
            f"{t4['measured_threshold']:.6f} found by bisection — a relative "
            f"error of {t4['relative_error']:.1e}. Below it the slice stays in "
            "the annulus; above it, it wraps.\n\n"
            "AND THEN NOTHING ACCUMULATES. Driving the gain up buys unsigned "
            f"crossings — {t6['most_unsigned_crossings']} of them at the "
            "hardest drive tested — while the SIGNED total stays at "
            f"{t6['worst_signed_total']} and the winding number at "
            f"{t6['worst_winding_number']}, at every gain and every time, with "
            "the crossing ledger and an independently computed degree agreeing. "
            "The reason is not subtle and is worth stating plainly: the curve "
            "is a GRAPH r = f(σ) over the circle with f single-valued, so every "
            "outward crossing of the seam is paid for by an inward one. A "
            "height field cannot wind.\n\n"
            "DIFFERENT WAVES DIFFER IN ARC, NOT IN TOPOLOGY. At one common "
            f"gain, pulses from {WIDTHS[0]} down to {WIDTHS[-1]} cross the seam "
            "the same number of times but put "
            f"{t7['arc_range'][1]:.3f} down to {t7['arc_range'][0]:.3f} of the "
            f"circle on the far sheet — a {t7['arc_ratio']:.1f}× spread. "
            "Divided by the pulse width that is "
            f"{t7['mean_arc_over_pulse_width']:.2f} for all of them, spread "
            f"{t7['arc_over_pulse_width_spread']:.2f}: the far-sheet arc simply "
            "IS the pulse. The winding column is zero all the way down.\n\n"
            "WHAT THIS RULES OUT. A crossing rule of this kind cannot "
            "manufacture topological charge from the amplitude of a scalar "
            "height, however hard it is driven and however many times the seam "
            "is crossed. A stable topological object would have to come from "
            "somewhere else — from a curve that is free to stop being a graph, "
            "which is a statement about what the next construction needs "
            "rather than a defect in this one.\n\n"
            "SCOPE. The crossing rule is a representation choice, not a derived "
            "boundary condition: nothing here makes the wave dynamically aware "
            "of the seam, and the field is a linear scalar on a fixed round "
            "background."
        )
    else:
        verdict_class = "CIRCLE_SLICE_INCONCLUSIVE"
        verdict = ("INCONCLUSIVE. A check failed; review the slice identity, "
                   "the lobe arrival, the wrap threshold, the closure, the "
                   "crossing ledger, or the arc comparison.")

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "a great-circle slice of the scalar wave, drawn in an annulus whose "
            "boundaries are glued"
        ),
        "the_slice": "d = |σ|, both halves of the wave on one circle",
        "the_rule": "R_outer is glued to R_inner; the radial direction is S¹",
        "the_result": "signed crossings and winding number are identically zero",
        "the_reason": "the curve is a graph r = f(σ) with f single-valued",
        "geometry": {
            "pulse_width": PULSE_WIDTH, "n_sigma": N_SIGMA,
            "n_radial": N_RADIAL, "common_gain": COMMON_GAIN,
            "widths": list(WIDTHS),
            "antipodal_time": ANTIPODAL_TIME, "return_time": RETURN_TIME,
        },
        "tests": tests,
        "n_passed": tests[-1]["n_passed"],
        "n_total": tests[-1]["n_total"],
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    t = {x["name"]: x for x in s["tests"]}
    out: List[str] = []
    out.append("# A circle slice, and a bulk the wave can wrap through\n")
    out.append(f"_{s['timestamp_utc']}_\n")
    out.append("> A **representation** of the scalar field on a fixed round "
               "`S²`, cut by one great circle — not a derived boundary "
               "condition.\n")

    c = t["T2_the_slice_is_the_sphere"]
    out.append("## The slice is the sphere\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| against the (θ, φ) route | `{c['worst_mismatch_against_the_sphere']:.1e}` |")
    out.append(f"| mirror asymmetry `σ → −σ` | `{c['worst_mirror_asymmetry']:.1e}` |")
    d = t["T3_the_lobes_meet_at_the_antipode"]
    out.append(f"| lobe arrival time | {d['arrival_time']:.4f} |")
    out.append(f"| lobe positions | {d['final_left_lobe']:+.4f} / "
               f"{d['final_right_lobe']:+.4f} |\n")

    c = t["T4_the_wrap_threshold"]
    out.append("## The wrap threshold\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| run peak `max\\|u\\|` | {c['run_peak']:.6f} |")
    out.append(f"| predicted `(R_outer − R_mid)/peak` | {c['predicted_threshold']:.6f} |")
    out.append(f"| found by bisection | {c['measured_threshold']:.6f} |")
    out.append(f"| relative error | `{c['relative_error']:.1e}` |\n")

    c = t["T6_the_seam_is_crossed_in_pairs"]
    out.append("## Amplitude buys crossings, never charge\n")
    out.append("| gain / threshold | unsigned | signed | winding | sheets |")
    out.append("|---:|---:|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['gain_over_threshold']:.1f} | {r['unsigned']} | "
                   f"{r['signed']:+d} | {r['winding']:+d} | {r['sheets']} |")
    out.append("\nThe curve is a graph `r = f(σ)` with `f` single-valued, so "
               "every outward crossing is paid for by an inward one. **A height "
               "field cannot wind.**\n")

    c = t["T7_different_waves_differ_in_arc"]
    out.append("## Different waves, one bulk, one gain\n")
    out.append("| pulse width | crossings | arc on the far sheet | arc / width | winding |")
    out.append("|---:|---:|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['pulse_width']:.2f} | {r['max_crossings']} | "
                   f"{r['arc_on_the_far_sheet']:.4f} | "
                   f"{r['arc_over_pulse_width']:.3f} | {r['winding']:+d} |")
    out.append(f"\nA {c['arc_ratio']:.1f}× spread in arc, "
               f"{c['mean_arc_over_pulse_width']:.2f} per unit pulse width for "
               "all of them — the far-sheet arc **is** the pulse.\n")

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
    out = here / "runs" / f"{ts}_circle_slice_probe"
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
