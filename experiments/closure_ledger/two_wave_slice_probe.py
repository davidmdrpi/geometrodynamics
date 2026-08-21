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
    measure_meeting_mid_flight_is_harder,
    measure_the_contact_is_an_arc_the_band_covers,
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
        "checks": checks,
        "passed": sum(1 for c in checks if c["pass"]),
        "total": len(checks),
    }


def render_markdown(summary: dict) -> str:
    where = next(c for c in summary["checks"] if c["id"] == "T2")["detail"]
    ident = next(c for c in summary["checks"] if c["id"] == "T1")["detail"]
    arc = next(c for c in summary["checks"] if c["id"] == "T3")["detail"]
    mid = next(c for c in summary["checks"] if c["id"] == "T4")["detail"]
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
