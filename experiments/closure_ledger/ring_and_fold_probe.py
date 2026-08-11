"""
Can the converging ring reach across the gap — and can it ever intersect?

> Framing: a REPRESENTATION of the scalar field on a fixed round S², drawn in a
> vacuole.  The tangential mixing λ is a modelling choice and is reported as one.

THE QUESTION
────────────
The scalar wave is not only its focal pulse.  A ring leaves the source, thins to
a minimum at the equator, and then GROWS as it converges on the antipode.  So it
is fair to ask whether that ring — rather than the pulse — is what first reaches
across the vacuole, and whether shrinking the gap or raising the energy can make
it intersect.

Those are two questions with different answers, and the useful result is that
they are controlled by different knobs.

WHAT IS CHECKED
───────────────
T2  THE RING IS REAL AND IT GROWS.  0.156 at the equator rising to 0.933 at the
    focus — a factor of 5.97 — following the 1/√(sin d) law for a closing ring
    to a mean ratio of 1.034.

T3  IT REACHES ACROSS, EASILY.  The threshold is exactly δ/max|u|, so both knobs
    buy it.  Shrinking the gap also buys LEAD: at δ = 0.09, ε = 0.80 the ring
    spans from d = 1.83 — just past the equator — for the whole converging leg,
    a lead of 1.41 on the focal pulse.  Not an instant at the focus, a state.

T4  AND IT NEVER INTERSECTS.  Swept over gap and gain with a real segment
    intersection test: zero, at up to 346 seam crossings.  A curve r = f(σ) with
    f single-valued is EMBEDDED by construction.  This is the winding
    obstruction seen from the side, and neither knob touches it.

T5  THE DETECTOR WORKS.  A control, so T4 is not a broken counter: a circle
    gives 0, a limaçon with an inner loop, a lemniscate and a folded loop all
    give more.

T6  WHAT DOES FOLD.  Give the material a tangential displacement,
    σ = σ₀ + λ ε ∂_σu, and the map folds where 1 + λ ε ∂²_σu < 0.  Predicted
    threshold λε = 1/max(−∂²_σu), found by bisection to a relative 1e-12.  It
    folds first ON THE CONVERGING RING, 0.016 from the antipode.

T7  FOLDING IS NECESSARY, NOT SUFFICIENT.  A crossing always comes with a fold
    (crossing-without-fold is zero at every drive — the graph theorem again);
    a fold need not cross, because the two branches can stay radially apart.

T8  THE TWO KNOBS ARE ORTHOGONAL.  The fold threshold is INDEPENDENT of the gap
    — spread 0.0 across δ = 0.26, 0.12, 0.05 — while the spanning threshold
    scales with δ directly.  Shrinking the vacuole can never buy an
    intersection.  What the fold threshold does scale with is the pulse:
    λε ≈ 0.385 w², spread 3.7%, and it converges under refinement.

SCOPE
─────
λ is not derived from the scalar equation; λ = 0 is exactly the height field of
the earlier probes.  What is derived is the threshold given λ, its independence
from the gap, and its w² scaling.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.slice_folding import (
    measure_a_graph_never_intersects,
    measure_folding_ignores_the_gap,
    measure_folding_is_necessary_for_crossing,
    measure_the_fold_threshold,
    measure_the_fold_threshold_converges,
    measure_the_fold_threshold_scales_with_the_pulse,
    measure_the_intersection_test_can_see_one,
    measure_the_ring_grows_as_it_converges,
    measure_when_the_ring_spans_the_gap,
)

PULSE_WIDTH = 0.18
DELTAS = (0.26, 0.16, 0.09, 0.05, 0.03)
GAINS = (0.10, 0.20, 0.40, 0.80)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "claim": ("determine whether the converging ring, not the focal pulse, "
                  "can reach across the vacuole — and whether it can intersect"),
        "reaching": "ε·max|u| ≥ δ",
        "folding": "1 + λ ε ∂²_σu < 0",
        "pass": True,
    }


def t2_the_ring_grows() -> dict:
    r = measure_the_ring_grows_as_it_converges(pulse_width=PULSE_WIDTH)
    return {"name": "T2_the_ring_grows", **r,
            "pass": bool(r["the_ring_thins_then_grows"]
                         and r["follows_one_over_root_sin"])}


def t3_it_reaches_across() -> dict:
    r = measure_when_the_ring_spans_the_gap(deltas=DELTAS, gains=GAINS,
                                            pulse_width=PULSE_WIDTH)
    return {"name": "T3_it_reaches_across", **r,
            "pass": bool(r["threshold_is_delta_over_peak"]
                         and r["shrinking_the_gap_buys_dwell"]
                         and r["shrinking_the_gap_buys_lead"]
                         and r["the_converging_ring_spans_well_before_the_pulse"])}


def t4_it_never_intersects() -> dict:
    r = measure_a_graph_never_intersects(deltas=(0.26, 0.09, 0.03),
                                         gains=(0.5, 2.0, 10.0), frames=50,
                                         n_sigma=801)
    return {"name": "T4_it_never_intersects", **r,
            "pass": bool(r["a_graph_is_embedded"]
                         and r["and_it_was_genuinely_driven"])}


def t5_the_detector_works() -> dict:
    r = measure_the_intersection_test_can_see_one()
    return {"name": "T5_the_detector_works", **r,
            "pass": bool(r["the_detector_works"])}


def t6_what_does_fold() -> dict:
    r = measure_the_fold_threshold(frames=120)
    c = measure_the_fold_threshold_converges(pulse_width=PULSE_WIDTH)
    return {"name": "T6_what_does_fold", **r,
            "convergence": c["rows"], "last_step_drift": c["last_step_drift"],
            "pass": bool(r["the_threshold_is_front_curvature"]
                         and r["past_it_the_curve_self_intersects"]
                         and r["below_it_the_curve_does_not"]
                         and r["it_folds_at_the_converging_ring"]
                         and c["it_converges"])}


def t7_folding_is_necessary_not_sufficient() -> dict:
    r = measure_folding_is_necessary_for_crossing()
    return {"name": "T7_folding_is_necessary_not_sufficient", **r,
            "pass": bool(r["a_crossing_always_folds"]
                         and r["nothing_crosses_below_threshold"])}


def t8_the_knobs_are_orthogonal() -> dict:
    g = measure_folding_ignores_the_gap(frames=80)
    w = measure_the_fold_threshold_scales_with_the_pulse(frames=100)
    return {"name": "T8_the_knobs_are_orthogonal", **g,
            "width_rows": w["rows"], "width_constant": w["mean_constant"],
            "width_spread": w["spread"],
            "pass": bool(g["fold_threshold_is_gap_independent"]
                         and g["span_threshold_scales_with_the_gap"]
                         and w["scales_as_the_width_squared"]
                         and w["narrower_folds_sooner"])}


def t9_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T9_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_the_ring_grows(), t3_it_reaches_across(),
             t4_it_never_intersects(), t5_the_detector_works(),
             t6_what_does_fold(), t7_folding_is_necessary_not_sufficient(),
             t8_the_knobs_are_orthogonal()]
    tests.append(t9_assessment(tests))
    t2, t3, t4, t6, t7, t8 = (tests[1], tests[2], tests[3], tests[5],
                              tests[6], tests[7])

    if all(t["pass"] for t in tests):
        verdict_class = "THE_RING_REACHES_BUT_ONLY_A_FOLD_CROSSES"
        verdict = (
            "THE RING REACHES ACROSS, BUT ONLY A FOLD CROSSES. The converging "
            "ring really does get there first, and the gap really is a knob — "
            "but it is not the knob that buys an intersection, and no setting "
            "of it ever will be.\n\n"
            "THE RING IS REAL AND IT GROWS. It thins to "
            f"{t2['equator_height']:.3f} at the equator and rises to "
            f"{t2['peak_height']:.3f} at the focus, a factor of "
            f"{t2['growth_from_the_equator']:.2f}, following the 1/√(sin d) law "
            f"for a closing ring to a mean ratio of {t2['law_ratio_mean']:.3f} "
            f"(spread {t2['law_ratio_spread']:.3f}). All of that is before the "
            "focal pulse.\n\n"
            "AND IT REACHES ACROSS EASILY. The threshold is exactly δ/max|u|, "
            "so raising the energy and shrinking the gap both buy it. Shrinking "
            "the gap buys something extra — LEAD. At the tightest setting "
            "tested the ring spans from just past the equator for the whole "
            f"converging leg, a lead of {t3['longest_lead_before_the_focus']:.2f} "
            f"on the pulse and a dwell of {t3['max_dwell_fraction']:.2f} of the "
            "run. The reaching stops being an instant at the focus and becomes "
            "a sustained state.\n\n"
            "AND IT STILL NEVER INTERSECTS. Swept over gap and gain with a real "
            "segment-intersection test — validated against a limaçon, a "
            "lemniscate and a folded loop so it is not a broken counter — the "
            f"answer is {t4['worst_self_intersections']}, at up to "
            f"{t4['most_seam_crossings']} seam crossings. A curve r = f(σ) with "
            "f single-valued puts exactly one radius at each direction, so it "
            "is embedded by construction and two of its points cannot occupy "
            "the same place. This is the winding obstruction seen from the "
            "side, and neither knob touches it.\n\n"
            "WHAT DOES FOLD IS TANGENTIAL FREEDOM. Let each material point move "
            "sideways as well as outward and the map σ₀ ↦ σ can fold. The "
            f"threshold λε = {t6['predicted_gain']:.5f} is predicted from the "
            "front's curvature alone and found by bisection to a relative "
            f"{t6['relative_error']:.1e}; it converges under joint refinement of "
            f"the angular and radial grids (last step {t6['last_step_drift']:.1e}). "
            "Past it the curve self-intersects, below it it does not — and it "
            f"folds first at {t6['distance_to_the_antipode']:.4f} from the "
            "antipode, ON the converging ring, at the moment of tightest focus. "
            "The intuition was right about where to look.\n\n"
            "FOLDING IS NECESSARY, NOT SUFFICIENT. A crossing always comes with "
            f"a fold — crossing-without-fold is {t7['crossing_without_fold']} at "
            "every drive tested, which is the embedded-graph theorem again — "
            f"while fold-without-crossing happens {t7['fold_without_crossing']} "
            "times, because the two branches of a folded map can stay radially "
            "apart and never meet in the plane.\n\n"
            "AND THE TWO KNOBS ARE ORTHOGONAL. The fold threshold does not know "
            f"the gap exists: spread {t8['fold_threshold_spread']:.1e} across a "
            "fivefold range of δ, while the spanning threshold scales with δ "
            "directly. What the fold threshold does scale with is the pulse — "
            f"λε ≈ {t8['width_constant']:.3f} w², spread "
            f"{t8['width_spread']:.3f} — so narrow fronts fold sooner because "
            "folding is about how sharply the front turns, not how tall it is. "
            "Reducing the distance between the shells changes when the wave "
            "arrives at them; it changes nothing about whether it can cross "
            "itself.\n\n"
            "SCOPE. λ is a modelling choice, not derived from the scalar "
            "equation, and λ = 0 is exactly the height field of the earlier "
            "probes. What is derived is the threshold given λ, its independence "
            "from the gap, and its w² scaling."
        )
    else:
        verdict_class = "RING_AND_FOLD_INCONCLUSIVE"
        verdict = ("INCONCLUSIVE. A check failed; review the ring growth, the "
                   "spanning sweep, the intersection sweep, the detector "
                   "control, the fold threshold, or the orthogonality of the "
                   "two knobs.")

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "the converging ring against the vacuole gap, and the tangential "
            "freedom that finally folds the curve"
        ),
        "reaching": "set by the gap",
        "folding": "set by the front's curvature, independent of the gap",
        "geometry": {"pulse_width": PULSE_WIDTH, "deltas": list(DELTAS),
                     "gains": list(GAINS)},
        "tests": tests,
        "n_passed": tests[-1]["n_passed"],
        "n_total": tests[-1]["n_total"],
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    t = {x["name"]: x for x in s["tests"]}
    out: List[str] = []
    out.append("# The ring reaches across, but only a fold crosses\n")
    out.append(f"_{s['timestamp_utc']}_\n")

    c = t["T2_the_ring_grows"]
    out.append("## The ring grows as it converges\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| equator height (d = {c['equator_distance']:.3f}) | {c['equator_height']:.4f} |")
    out.append(f"| peak height | {c['peak_height']:.4f} |")
    out.append(f"| growth | {c['growth_from_the_equator']:.2f}× |")
    out.append(f"| against `1/√(sin d)` | {c['law_ratio_mean']:.3f} ± {c['law_ratio_spread']:.3f} |\n")

    c = t["T3_it_reaches_across"]
    out.append("## Reaching across — and the lead it gets on the pulse\n")
    out.append("| δ | ε | dwell | spans from d | lead |")
    out.append("|---:|---:|---:|---:|---:|")
    for r in c["rows"]:
        if r["converging_span_distance"] is None:
            continue
        out.append(f"| {r['delta']:.2f} | {r['gain']:.2f} | "
                   f"{r['dwell_fraction']:.3f} | "
                   f"{r['converging_span_distance']:.3f} | "
                   f"{r['lead_before_the_focus']:.3f} |")
    out.append("")

    c = t["T4_it_never_intersects"]
    out.append("## ...and never intersecting\n")
    out.append("| δ | ε | seam crossings | self-intersections |")
    out.append("|---:|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['delta']:.2f} | {r['gain']:.1f} | "
                   f"{r['seam_crossings']} | **{r['self_intersections']}** |")
    out.append("")

    c = t["T6_what_does_fold"]
    out.append("## The fold threshold\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| predicted `λε` | {c['predicted_gain']:.6f} |")
    out.append(f"| found by bisection | {c['measured_gain']:.6f} |")
    out.append(f"| relative error | `{c['relative_error']:.1e}` |")
    out.append(f"| folds at, from the antipode | {c['distance_to_the_antipode']:.4f} |")
    out.append(f"| convergence drift | `{c['last_step_drift']:.1e}` |\n")

    c = t["T8_the_knobs_are_orthogonal"]
    out.append("## The two knobs are orthogonal\n")
    out.append("| δ | span threshold | fold threshold |")
    out.append("|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['delta']:.2f} | {r['span_threshold_gain']:.3f} | "
                   f"{r['fold_threshold_product']:.6f} |")
    out.append(f"\nFold-threshold spread `{c['fold_threshold_spread']:.1e}`, "
               f"and `λε ≈ {c['width_constant']:.3f} w²` across the pulse "
               f"widths (spread {c['width_spread']:.3f}).\n")

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
    out = here / "runs" / f"{ts}_ring_and_fold_probe"
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
