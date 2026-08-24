#!/usr/bin/env python3
"""
Hyperspherical probe — what higher dimension does to the bulk picture
=====================================================================

THE PREMISE. Three-dimensional intuition models an extra dimension as "the same
sphere, with another direction available". That is wrong in ways that matter
here, and the point of this probe is to measure them rather than assert them.

The one that matters most for the finite bearing is T5: a packet can keep
EXACTLY the same angular footprint while its physical overlap shrinks as f0^n.
PR #268's review had already forced the overlap SIZE law to be scoped by
dimension; this follows that scoping to its conclusion.

T0  A CORRECTION FIRST, BECAUSE IT IS THE USUAL ONE. The familiar "spheres peak
    at 5D" is about the VOLUME OF THE UNIT BALL,
        V_d = pi^(d/2)/Gamma(d/2+1) ,  peaking at d = 5 ,
    not the surface measure of the unit sphere,
        A_(d-1) = 2 pi^(d/2)/Gamma(d/2) ,  peaking at d = 7, i.e. S^6 in R^7.
    AND NEITHER IS A FACT ABOUT DIMENSION ALONE. V_d and A_(d-1) carry
    different units at different d, so comparing across d picks a length scale
    and the peak follows it: at R = 0.5 the ball peaks at d = 1, at R = 2 at
    d = 24. "The peak at 5" is a statement about the UNIT ball, nothing more.

T1  *** ALMOST ALL OF A SPHERE IS AT THE EQUATOR OF ANY CHOSEN POINT. *** The
    shell at geodesic angle chi from a pole has measure ~ sin^(n-1)(chi) dchi.
    Writing chi = pi/2 + delta gives sin^(n-1) ~ exp(-(n-1) delta^2/2), so the
    occupied band has width ~1/sqrt(n). MEASURED: std(chi) x sqrt(n) runs
    0.9669 -> 0.9993 -> 1.000000 over n = 2 ... 1000. Almost every point is
    nearly 90 degrees from any chosen point -- concentrated, not scattered.

T2  *** SO THE ANTIPODE IS AN EXTREMELY NON-GENERIC RELATION. *** For random
    unit vectors x.y has mean 0 and width 1/sqrt(n), so alpha -> pi/2. Random
    points do NOT tend toward alpha = pi. The fraction with alpha > 0.99 pi is
    3.2e-04 on S^2 and indistinguishable from zero by n = 10.
    THAT SHARPENS WHAT THE ANTIPODAL IDENTIFICATION IS CLAIMING: selecting
    x <-> -x is not "pairing a point with the far one". It picks a
    vanishing-measure relation out of an enormous nearly-orthogonal majority,
    and gets MORE non-generic as dimension rises. (It does not make the
    identification correct -- only not bland.)

T3  THE WAVEFRONT IS A POINT, THEN AN EQUATORIAL SHELL, THEN A POINT. A wave
    launched at a pole of S^n has transverse measure |S^(n-1)| (R sin chi)^(n-1)
    -- zero at the pole, maximal at the equator, zero again at the antipode.
    Higher n makes the middle dominate more (10.0% of the measure within 0.1
    rad of the equator on S^2, 20.2% on S^7) while the poles get relatively
    tinier. So the dimensional-descent picture STRENGTHENS with dimension, and
    with it the contrast between where the propagation measure lives and where
    the congruence reconverges.

T4  THE EMBEDDING CENTRE IS A DIRECTION SPACE COLLAPSED TO A POINT. An
    intrinsic S^n has no centre; embedded in R^(n+1) it has exactly one, and
    that point is the unique fixed point of SO(n+1). Every direction gives a
    radial line through it. REGULARISING f = 0 -> f0 > 0 RESTORES THAT
    DIRECTION SPACE AT FINITE SIZE, with measure f0^n |S^n|. Offered as an
    INTERPRETATION of PR #268, which stands on its own geometry; nothing in
    regularized_center depends on it.

T5  *** THE COLLAPSE IS f^n, NOT f. *** For ds^2 + f^2 dOmega_n^2 the
    transverse measure of an angular patch is f^n dOmega_n, so squeezing from
    F to f0 costs (f0/F)^n:
        f0/F = 1e-03  ->  1e-03 (S^1), 1e-06 (S^2), 1e-09 (S^3), 1e-12 (S^4).
    The 2-D drawing's l ~ f understates the S^3 case by a factor of a MILLION.
    THE ANGULAR OVERLAP CAN STAY FINITE WHILE THE PHYSICAL OVERLAP COLLAPSES AS
    f0^n. That changes what the finite-bearing picture is a picture OF: much
    less like two ribbons squeezing together, much more like a vast angular
    configuration space packed into an extremely small proper region. The
    yes/no overlap criterion is untouched -- it was always angular.

T6  *** ORIENTABILITY FLIPS WITH PARITY -- AND THIS REPO USES TWO QUOTIENTS
    THAT ARE ALWAYS OPPOSITE. *** RP^n is orientable iff n is ODD, because the
    antipodal map is the restriction of -I on R^(n+1) and det(-I) = (-1)^(n+1).
    A parity effect, not a size effect. The consequence here is specific: the
    repo carries TWO antipodal quotients at CONSECUTIVE n --
        spatial slice        S^d/+-  = RP^d      (d = 3 -> RP^3, ORIENTABLE)
        two-body exchange  (R^d\\0)/+- ~ RP^(d-1)  (d = 3 -> RP^2, NON-orientable)
    and it is the RP^2 whose Pin- structure makes the throat a spinor. Being
    one apart they are ALWAYS of opposite parity, so raising the spatial
    dimension by one SWAPS WHICH IS WHICH: at d = 4 the exchange space is the
    orientable RP^3 and the spin-statistics mechanism would have to be
    re-derived rather than carried over. d = 2 fails differently: RP^1 ~ S^1
    has pi_1 = Z, the braid group, not Z_2.
    Stated as a consequence IF the dimension moved. Not an argument that it
    should.

T7  AND S^3 IS NOT A GENERIC SPHERE, SO NOTHING ABOVE EXTRAPOLATES SMOOTHLY.
    S^3 = SU(2), the unit quaternions; it is PARALLELIZABLE, and the frame
    (q i, q j, q k) is verified orthonormal and tangent to 6.7e-16 at every
    sampled point; and it carries the Hopf fibration S^1 -> S^3 -> S^2. Only
    S^1, S^3 and S^7 are parallelizable -- S^2 admits no nowhere-zero tangent
    field at all. Dimensions do not simply add room: some are algebraically
    special, and the one this repository is built on is one of them.

WHAT IS PUT IN, AND WHAT IS NOT CLAIMED
───────────────────────────────────────
Every measurement here is about MEASURE, ORIENTATION AND FRAMES ON ROUND
SPHERES. No field equation is solved, nothing evolves, and no throat profile
appears except through the f^n weighting that T5 shares with #268. Nothing here
argues the model should move to a different dimension -- T6 says what would
have to be redone if it did. And T4 is explicitly an interpretation, not a
result.

    python -m experiments.closure_ledger.hyperspherical_probe
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.hyperspherical import (
    measure_s3_is_not_a_generic_sphere,
    measure_the_antipode_is_a_vanishing_measure_relation,
    measure_the_bearing_is_the_blown_up_direction_space,
    measure_the_front_is_a_point_then_a_shell_then_a_point,
    measure_the_patch_collapses_as_f_to_the_n,
    measure_the_peak_is_the_ball_not_the_sphere_and_needs_a_unit,
    measure_the_quotient_flips_orientability_with_parity,
    measure_the_shell_measure_concentrates_at_the_equator,
)


def run_probe() -> dict:
    checks: List[dict] = []

    peak = measure_the_peak_is_the_ball_not_the_sphere_and_needs_a_unit()
    checks.append({
        "id": "T0", "name": "the peak is the ball's, and it needs a unit",
        "detail": peak,
        "pass": bool(peak["the_ball_peaks_at_five"]
                     and peak["the_sphere_peaks_two_later"]
                     and peak["the_peak_moves_with_the_unit"])})

    shell = measure_the_shell_measure_concentrates_at_the_equator()
    checks.append({
        "id": "T1", "name": "*** almost all of a sphere is at the equator ***",
        "detail": shell,
        "pass": bool(shell["the_mean_is_always_the_equator"]
                     and shell["std_times_sqrt_n_tends_to_one"]
                     and shell["the_band_narrows_as_one_over_sqrt_n"])})

    anti = measure_the_antipode_is_a_vanishing_measure_relation()
    checks.append({
        "id": "T2", "name": "*** the antipode is a vanishing-measure relation ***",
        "detail": anti,
        "pass": bool(anti["std_scales_as_one_over_sqrt_n"]
                     and anti["the_mean_angle_is_a_right_angle"]
                     and anti["near_antipodal_fraction_collapses"])})

    front = measure_the_front_is_a_point_then_a_shell_then_a_point()
    checks.append({
        "id": "T3", "name": "point -> equatorial shell -> point",
        "detail": front,
        "pass": bool(front["the_peak_is_always_the_equator"]
                     and front["the_poles_vanish"]
                     and front["the_middle_dominates_more_with_dimension"])})

    centre = measure_the_bearing_is_the_blown_up_direction_space()
    checks.append({
        "id": "T4", "name": "the bearing as the blown-up direction space",
        "detail": centre,
        "pass": bool(centre["every_bearing_measure_is_positive"]
                     and centre["and_shrinks_faster_in_higher_dimension"])})

    patch = measure_the_patch_collapses_as_f_to_the_n()
    checks.append({
        "id": "T5", "name": "*** the collapse is f^n, not f ***",
        "detail": patch,
        "pass": bool(patch["the_angular_footprint_is_untouched"]
                     and patch["the_drawn_2d_length_understates_it"]
                     and patch["the_collapse_steepens_with_dimension"])})

    parity = measure_the_quotient_flips_orientability_with_parity()
    checks.append({
        "id": "T6", "name": "*** orientability flips with parity; two quotients ***",
        "detail": parity,
        "pass": bool(parity["orientable_iff_odd"]
                     and parity["the_two_quotients_are_always_opposite"]
                     and parity["raising_the_dimension_swaps_them"])})

    special = measure_s3_is_not_a_generic_sphere()
    checks.append({
        "id": "T7", "name": "S^3 is not a generic sphere",
        "detail": special,
        "pass": bool(special["the_frame_is_global_and_nowhere_zero"]
                     and special["the_fibre_is_a_circle_not_a_point"])})

    return {
        "probe": "hyperspherical",
        "question": "what does higher dimension actually do to the bulk and "
                    "intersection picture -- as opposed to what 3-D intuition "
                    "suggests it does?",
        "answer": "a packet can keep exactly the same angular footprint while "
                  "its physical overlap collapses as f0^n; almost every point "
                  "is at the equator of any other, so the antipodal relation "
                  "is vanishingly non-generic; and orientability of the "
                  "antipodal quotient flips with dimension parity",
        "headline": {
            "patch_collapse_at_squeeze_1e_minus_3": {
                f"S^{n}": (1e-3) ** n for n in (1, 2, 3, 4)},
            "understatement_of_the_2d_drawing_at_n_3":
                patch["understatement_factor_at_n_three"],
            "band_width_law": "Delta chi ~ 1/sqrt(n)",
            "orientable_quotients": "RP^n orientable iff n odd",
        },
        "checks": checks,
        "passed": sum(1 for c in checks if c["pass"]),
        "total": len(checks),
    }


def render_markdown(summary: dict) -> str:
    h = summary["headline"]
    lines = [
        "# Hyperspherical probe — what higher dimension does to the picture",
        "",
        f"**Question.** {summary['question']}",
        "",
        f"**Answer.** {summary['answer']}",
        "",
        f"**{summary['passed']}/{summary['total']} checks pass.**",
        "",
        "| id | check | result |",
        "|----|-------|--------|",
    ]
    for c in summary["checks"]:
        lines.append(f"| {c['id']} | {c['name']} | "
                     f"{'PASS' if c['pass'] else 'FAIL'} |")

    peak = next(c for c in summary["checks"] if c["id"] == "T0")["detail"]
    lines += ["", "## T0 — the peak is the ball's, and it needs a unit", "",
              f"The unit **ball** peaks at `d = {peak['ball_peaks_at']}`; the "
              f"unit **sphere**'s surface peaks at "
              f"`d = {peak['sphere_peaks_at']}` — that is "
              f"`{peak['the_peaking_sphere_is']}`.", "",
              "| `d` | `V_d` (ball) | `A_{d−1}` (sphere) | sphere |",
              "|--|--|--|--|"]
    for r in peak["rows"]:
        lines.append(f"| {r['dim']} | `{r['ball_volume']:.5f}` | "
                     f"`{r['sphere_area']:.5f}` | `{r['sphere']}` |")
    lines += ["", "And the peak is not a fact about dimension alone:", "",
              "| `R` | ball peaks at | sphere peaks at |", "|--|--|--|"]
    for s in peak["radius_scan"]:
        lines.append(f"| `{s['radius']}` | `d = {s['ball_peak_dim']}` | "
                     f"`d = {s['sphere_peak_dim']}` |")
    lines += ["", f"**Why.** {peak['why']} — so the familiar claim is "
                  f"{peak['so_the_familiar_claim_is']}.", ""]

    shell = next(c for c in summary["checks"] if c["id"] == "T1")["detail"]
    lines += ["## T1 — the band narrows as `1/√n`", "",
              "| `n` | `std(χ)` | `std·√n` | mass within `1/√n` of the equator |",
              "|--|--|--|--|"]
    for r in shell["rows"]:
        lines.append(f"| {r['sphere_dim']} | `{r['std_chi']:.6f}` | "
                     f"`{r['std_times_sqrt_n']:.6f}` | "
                     f"`{r['mass_within_one_over_sqrt_n']:.4f}` |")
    lines += ["", f"**{shell['what_it_means'].capitalize()}.**", ""]

    anti = next(c for c in summary["checks"] if c["id"] == "T2")["detail"]
    lines += ["## T2 — so the antipode is non-generic", "",
              "| ambient `n` | `std(x·y)·√n` | mean `α/π` | fraction with `α > 0.99π` |",
              "|--|--|--|--|"]
    for r in anti["rows"]:
        lines.append(f"| {r['ambient_dim']} | "
                     f"`{r['std_dot_times_sqrt_n']:.6f}` | "
                     f"`{r['mean_angle_over_pi']:.4f}` | "
                     f"`{r['fraction_near_antipodal']:.2e}` |")
    lines += ["", f"{anti['so_the_antipodal_identification_is_non_generic'].capitalize()}.",
              "", f"**What is not claimed:** {anti['what_is_not_claimed']}.", ""]

    patch = next(c for c in summary["checks"] if c["id"] == "T5")["detail"]
    lines += ["## T5 — the collapse is `fⁿ`, not `f`", "",
              f"`{patch['law']}`", "",
              "| squeeze `f₀/F` | `S¹` | `S²` | `S³` | `S⁴` |", "|--|--|--|--|--|"]
    for r in patch["rows"]:
        lines.append(f"| `{r['squeeze']:.0e}` | `{r['n_1']:.1e}` | "
                     f"`{r['n_2']:.1e}` | `{r['n_3']:.1e}` | `{r['n_4']:.1e}` |")
    lines += ["",
              f"The 2-D drawing understates the `S³` case by "
              f"`{h['understatement_of_the_2d_drawing_at_n_3']:.1e}`.", "",
              f"**So the picture is of** {patch['what_the_picture_is_of']} — "
              f"and {patch['the_yes_no_criterion_is_unaffected']}.", ""]

    parity = next(c for c in summary["checks"] if c["id"] == "T6")["detail"]
    lines += ["## T6 — parity, and two quotients that are always opposite", "",
              "| spatial `d` | spatial `ℝP^d` | exchange `ℝP^{d−1}` | `π₁`(exchange) | opposite? |",
              "|--|--|--|--|--|"]
    for r in parity["rows"]:
        lines.append(
            f"| {r['spatial_dim']} | `{r['spatial_quotient']}` "
            f"{'orientable' if r['spatial_orientable'] else '**non-orientable**'} | "
            f"`{r['exchange_quotient']}` "
            f"{'orientable' if r['exchange_orientable'] else '**non-orientable**'} | "
            f"`{r['exchange_pi1']}` | {'yes' if r['opposite_parity'] else 'no'} |")
    at = parity["at_the_repos_dimension"]
    lines += ["", f"At this repository's dimension: spatial **{at['spatial']}**, "
                  f"exchange **{at['exchange']}**.", "",
              f"**What is not claimed:** {parity['what_is_not_claimed']}.", ""]

    centre = next(c for c in summary["checks"] if c["id"] == "T4")["detail"]
    lines += ["## T4 — the bearing as the blown-up direction space", "",
              f"{centre['what_the_origin_loses'].capitalize()}; the bearing "
              f"{centre['what_the_bearing_restores']}.", "",
              "| `n` | `f₀` | `|S^n|` | bearing measure `f₀ⁿ|S^n|` |",
              "|--|--|--|--|"]
    for r in centre["rows"]:
        lines.append(f"| {r['angular_dim']} | `{r['neck']:.0e}` | "
                     f"`{r['unit_direction_space_measure']:.4f}` | "
                     f"`{r['bearing_measure']:.3e}` |")
    lines += ["", f"*{centre['this_is_an_interpretation_not_a_result'].capitalize()}.*",
              ""]
    return "\n".join(lines)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = (Path(__file__).resolve().parent / "runs"
           / f"{stamp}_hyperspherical_probe")
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0 if summary["passed"] == summary["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
