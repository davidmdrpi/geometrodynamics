"""Two waves on the circle slice: do the antipodal pulses connect inner to outer?

> Scope: exactly v46's construction, with a second wave added. A LINEAR scalar
> on a FIXED round background; the crossing rule that glues R_outer to R_inner is
> a REPRESENTATION choice, not a derived boundary condition; the two waves do not
> interact -- they are drawn on the same torus and the question is only whether
> their images meet. The gain is a DISPLAY amplitude.

THE QUESTION
────────────
v46 put ONE scalar wave on the great circle through a source and its antipode,
drew it as a radial height r = R_mid + eps u(|sigma|) in the vacuole, and glued
R_outer to R_inner so the radial direction is a circle. Its result was negative
and sharp:

    the curve is a GRAPH r = f(sigma), so its radial winding number is
    identically zero. Every outward crossing of the seam is paid for by an
    inward one. A HEIGHT FIELD CANNOT WIND, and one wave running to its
    antipode never meets itself.

The arc since then has needed TWO waves -- two mouths, and an interference term
that is bilinear and vanishes unless both are there. So:

    one wave pulsing OUTWARD and one pulsing INWARD, both refocusing at the
    antipode -- DO THEY CONNECT, INNER TO OUTER, AND WHERE?

THE ANSWER: YES. AT THE ANTIPODE, ON THE SEAM, AT THE REFOCUS.

And the threshold is not a new number. A single wave crosses the seam when
eps u = gap/2. The pair spans |delta| = 2 eps |u| of the radial circle and
touches through the seam when THAT reaches gap -- the same inequality. v46's
"the wave comes back inside the circle" and this round's "the two pulses connect
inner to outer" are ONE EVENT DESCRIBED TWICE.

WHAT IS CHECKED
───────────────
T1  *** THE IDENTITY. *** The single-wave wrap gain and the pair-contact gain
    are computed independently from the run's own peak and compared. They agree
    to ZERO -- not to a tolerance, because it is an identity: both are
    eps u = gap/2.

T2  *** WHERE. *** Driven to threshold at t = pi, one curve sits at exactly
    R_outer and the other at exactly R_inner, at sigma = +-pi. Glued, those are
    the same point. AND THE ROLES ARE THE OTHER WAY ROUND FROM THE GUESS: the
    antipodal refocus is a RAREFACTION, so it is the INWARD-driven wave that
    bulges out to R_outer and the outward-driven one that dips to R_inner.

T3  *** WHAT TWO WAVES CAN DO THAT ONE CANNOT. *** Below threshold nothing
    connects. At threshold there is a single tangency. Above it that point opens
    into ONE ARC, bounded by two genuine crossings, on which the band between
    the two curves covers the WHOLE radial circle -- no radius at those sigma is
    outside the pair. A single wave past its own wrap threshold does nothing of
    the kind: it crosses, reappears inside, and still leaves every radius
    outside itself, because it is a graph. Two graphs bound a band, and a band
    can be radially surjective. The v46 winding result is re-checked here at
    four gains and still holds.

T4  WHERE ELSE THEY COULD HAVE MET, AND WHY THEY DO NOT. With the sources at
    opposite ends the two travelling pulses cross at the quarter points at
    t = pi/2. That is the obvious place for a connection and it is the WORST
    one: the pulses partially cancel, so |u_A + u_B| is smaller there than
    either alone, and the threshold is 7.4x the co-located best and 9.0x the
    antipodal one. (A first draft guessed 4x from a partial scan.) The cheapest
    connection is always at a REFOCUS, in both configurations -- and a
    co-located pair gets there at about HALF the gain, because at a refocus both
    of a co-located pair are at peak while only one of an antipodal pair is.

WHAT IS STILL PUT IN
────────────────────
Everything v46 put in. The crossing rule is a representation choice; nothing
makes either wave dynamically aware of the seam or of the other wave; the field
is linear on a fixed background, so this measures whether two DRAWN curves meet,
not whether two physical waves interact. The gain is a display amplitude, so
every threshold quoted is a statement about the picture.

    python -m experiments.closure_ledger.two_wave_slice_probe
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.two_wave_slice import (
    measure_like_signs_are_one_case_not_two,
    measure_meeting_mid_flight_is_harder,
    measure_only_the_opposed_pair_connects_on_the_bisector,
    measure_the_bisector_is_degenerate_for_like_signs,
    measure_the_contact_is_an_arc_the_band_covers,
    measure_the_offset_slides_the_connection,
    measure_the_pair_touches_exactly_where_one_wave_wraps,
    measure_where_the_two_pulses_connect,
)


def run_probe() -> dict:
    checks: List[dict] = []

    ident = measure_the_pair_touches_exactly_where_one_wave_wraps()
    checks.append({"id": "T1", "name": "*** one wave wrapping IS two touching ***",
                   "detail": ident,
                   "pass": bool(ident["they_are_the_same_threshold"])})

    where = measure_where_the_two_pulses_connect()
    checks.append({"id": "T2", "name": "*** where: antipode, seam, refocus ***",
                   "detail": where,
                   "pass": bool(where["connected"]
                                and where["the_contact_is_on_the_seam"]
                                and where["the_contact_is_at_the_antipode"])})

    arc = measure_the_contact_is_an_arc_the_band_covers()
    checks.append({"id": "T3",
                   "name": "*** the band covers what one graph cannot ***",
                   "detail": arc,
                   "pass": bool(arc["nothing_connects_below_threshold"]
                                and arc["it_touches_at_threshold"]
                                and arc["one_arc_above"])})

    mid = measure_meeting_mid_flight_is_harder()
    checks.append({"id": "T4", "name": "meeting mid-flight is harder",
                   "detail": mid,
                   "pass": bool(mid["mid_flight_is_harder"]
                                and mid["both_are_cheapest_at_a_refocus"])})

    signs = measure_like_signs_are_one_case_not_two()
    checks.append({"id": "T5", "name": "two cases, not four (a limitation)",
                   "detail": signs,
                   "pass": bool(signs["there_are_two_cases_not_four"]
                                and signs["the_drawn_residue_is_one_ulp_of_r_mid"])})

    bis = measure_the_bisector_is_degenerate_for_like_signs()
    checks.append({"id": "T6",
                   "name": "*** the bisector is degenerate for like signs ***",
                   "detail": bis,
                   "pass": bool(
                       bis["the_like_signed_pair_never_separates_on_a_bisector"]
                       and bis["the_opposed_pair_always_does"]
                       and bis["the_far_bisector_is_the_cheaper_one"])})

    excl = measure_only_the_opposed_pair_connects_on_the_bisector()
    checks.append({"id": "T7",
                   "name": "*** an off-antipodal connection only one pair has ***",
                   "detail": excl,
                   "pass": bool(excl["every_offset_opens_an_arc"]
                                and excl["the_like_signed_pair_reaches_none_of_it"]
                                and excl["the_arc_is_centred_on_the_bisector"]
                                and excl["it_is_off_both_axes"])})

    slide = measure_the_offset_slides_the_connection()
    checks.append({"id": "T8", "name": "the slider: it moves, and it costs",
                   "detail": slide,
                   "pass": bool(slide["the_timing_is_the_pulse_crossing"]
                                and slide["it_rises_monotonically_except_at_the_endpoint"]
                                and slide["the_cheapest_connection_sits_on_one_of_the_four_axes"]
                                and slide["exclusive_is_not_cheap"])})

    return {
        "probe": "two_wave_slice",
        "question": "one wave pulsing outward and one pulsing inward, both "
                    "refocusing at the antipode -- do they connect, inner to "
                    "outer, and where?",
        "answer": "yes: at the antipode, on the seam, at the refocus -- and at "
                  "exactly the amplitude where one wave would have wrapped",
        "contact_gain": where["contact_gain"],
        "single_wave_wrap_gain": ident["single_wave_wrap_gain"],
        "sigma_over_pi": where["sigma_over_pi"],
        "off_antipodal_question": "does an inner-outer pair possess "
                                  "off-antipodal intersections that neither "
                                  "inner-inner nor outer-outer pairs possess?",
        "off_antipodal_answer": excl["answer"],
        "checks": checks,
        "passed": sum(1 for c in checks if c["pass"]),
        "total": len(checks),
    }


def render_markdown(summary: dict) -> str:
    where = next(c for c in summary["checks"] if c["id"] == "T2")["detail"]
    ident = next(c for c in summary["checks"] if c["id"] == "T1")["detail"]
    arc = next(c for c in summary["checks"] if c["id"] == "T3")["detail"]
    mid = next(c for c in summary["checks"] if c["id"] == "T4")["detail"]
    signs = next(c for c in summary["checks"] if c["id"] == "T5")["detail"]
    bis = next(c for c in summary["checks"] if c["id"] == "T6")["detail"]
    excl = next(c for c in summary["checks"] if c["id"] == "T7")["detail"]
    slide = next(c for c in summary["checks"] if c["id"] == "T8")["detail"]
    lines = [
        "# Two-wave slice probe — do the antipodal pulses connect?",
        "",
        f"**Question.** {summary['question']}",
        "",
        f"**Answer.** {summary['answer']}",
        "",
        "| | |",
        "|--|--|",
        f"| where | `σ = {where['sigma_over_pi']:+.0f}π` — the antipode |",
        f"| at what radius | outward-driven → `{where['radius_of_the_outward_wave']:.3f}` "
        f"= `R_inner`, inward-driven → `{where['radius_of_the_inward_wave']:.3f}` "
        f"= `R_outer` |",
        f"| distance to the seam | `{where['distance_to_the_seam']:.1e}` |",
        f"| contact gain there | `{where['contact_gain']:.6f}` |",
        f"| single-wave wrap gain | `{ident['single_wave_wrap_gain']:.6f}` "
        f"(at the run peak, `t = {ident['at_time']:.3f}`) |",
        f"| the two thresholds differ by | `{ident['relative_difference']:.1e}` |",
        "",
        f"**{summary['passed']}/{summary['total']} checks pass.**",
        "",
        "| id | check | result |",
        "|----|-------|--------|",
    ]
    for c in summary["checks"]:
        lines.append(f"| {c['id']} | {c['name']} | "
                     f"{'PASS' if c['pass'] else 'FAIL'} |")
    lines += ["", "## Below, at, and above threshold", "",
              "| gain / threshold | covered fraction | connected | arcs |",
              "|--|--|--|--|"]
    for r in arc["rows"]:
        lines.append(f"| `{r['gain_over_threshold']:.3f}` | "
                     f"`{r['covered_fraction']:.4f}` | "
                     f"{'yes' if r['connected'] else 'no'} | {r['arcs']} |")
    lines += [
        "",
        f"> {arc['what_one_wave_cannot_do']}",
        "",
        "## Where else they could have met",
        "",
        "| configuration | cheapest gain | at `t/π` | mid-flight penalty |",
        "|--|--|--|--|",
    ]
    for key in ("co_located", "antipodal_sources"):
        d = mid[key]
        lines.append(f"| {key.replace('_', ' ')} | "
                     f"`{d['cheapest_contact_gain']:.4f}` | "
                     f"`{d['at_time_over_pi']:.2f}` | "
                     f"`{d['mid_flight_penalty']:.1f}×` |")
    lines += ["", "The cheapest connection is always at a **refocus**. A "
                  "co-located pair reaches it at "
                  f"`{mid['antipodal_over_co_located']:.2f}×` less gain than an "
                  "antipodal one.", ""]

    lines += [
        "## Off the degenerate axis",
        "",
        f"**Question.** {summary['off_antipodal_question']}",
        "",
        f"**Answer.** {summary['off_antipodal_answer']}",
        "",
        "### There are two configurations, not four",
        "",
        "| | |",
        "|--|--|",
        f"| `(out,out)` vs `(in,in)`, as fields | `{signs['worst_as_fields']:.1e}` |",
        f"| the same, through the drawn radii | "
        f"`{signs['drawn_residue_in_ulps']:.1f}` ulp of `R_mid` |",
        "",
        f"> {signs['the_limitation']}",
        "",
        "### The bisector",
        "",
        "| offset `α/π` | bisector `σ/π` | like-signed `|δ|` | opposed threshold |"
        " far bisector | its threshold |",
        "|--|--|--|--|--|--|",
    ]
    for r in bis["rows"]:
        lines.append(
            f"| `{r['offset_over_pi']:.2f}` | `{r['near_bisector_over_pi']:+.4f}` | "
            f"`{r['near_like_signed_separation']:.1e}` | "
            f"`{r['near_opposed_threshold']:.4f}` | "
            f"`{r['far_bisector_over_pi']:+.4f}` | "
            f"`{r['far_opposed_threshold']:.4f}` |")
    lines += [
        "",
        f"> {bis['why']}",
        "",
        "### The exclusive arc",
        "",
        f"Driven to `{excl['drive_over_threshold']:.2f}×` the opposed pair's own "
        "bisector threshold:",
        "",
        "| offset `α/π` | bisector `σ/π` | arc (rad) | centre − bisector |"
        " like-signed samples on it | to nearest source | to nearest antipode |",
        "|--|--|--|--|--|--|--|",
    ]
    for r in excl["rows"]:
        lines.append(
            f"| `{r['offset_over_pi']:.2f}` | `{r['bisector_over_pi']:+.4f}` | "
            f"`{r['opposed_arc']:.4f}` | "
            f"`{r['centre_minus_bisector']:.1e}` | "
            f"**{r['like_signed_samples_on_that_arc']}** | "
            f"`{r['distance_to_the_nearest_source']:.3f}π` | "
            f"`{r['distance_to_the_nearest_antipode']:.3f}π` |")
    lines += [
        "",
        "### What the slider does",
        "",
        "| | |",
        "|--|--|",
        f"| threshold at `α = 0` | `{slide['threshold_at_zero_offset']:.4f}` |",
        f"| threshold at `α = π` | `{slide['threshold_at_pi']:.4f}` |",
        f"| over the sweep | `{slide['threshold_range']:.2f}×` |",
        f"| timing vs `t = α/2` | `{slide['worst_timing_error_over_pi']:.4f}π` |",
        f"| price of exclusivity | "
        f"`{slide['price_of_exclusivity_range'][0]:.2f}–"
        f"{slide['price_of_exclusivity_range'][1]:.2f}×` |",
        f"| cheapest point pins to an axis from | "
        f"`α = {slide['it_pins_to_an_axis_from_offset_over_pi']:.3f}π` |",
        "",
    ]
    if slide["non_monotone_steps"]:
        d = slide["non_monotone_steps"][0]
        lines += [
            f"It rises monotonically but for a turn-over of "
            f"`{100 * d['relative']:.3f}%` into the symmetric endpoint "
            f"`α = {d['offset_over_pi']:.2f}π`. A first pass put that at "
            "`0.08%` and confirmed it by refining the *time* grid — the wrong "
            "axis. The bisector is evaluated off-grid here.",
            "",
        ]
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
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_two_wave_slice_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
