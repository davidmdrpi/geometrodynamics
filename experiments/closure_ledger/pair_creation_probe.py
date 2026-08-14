"""
Pair creation belongs to the collision topology, not to focusing

> Scope: kinematics on a fixed round sphere, c = 1, no metric evolution and no
> backreaction.  The Breit-Wheeler threshold and cross-section are IMPORTED QED,
> not derived — the cross-section is the textbook closed form and is checked
> against its known peak.  Treating a wavefront's null rays as photons is a
> CORRESPONDENCE, stated rather than justified.  No rate is computed: a rate
> needs a photon number density, which a classical amplitude does not supply.
> Calling the crossing chord a throat is this program's INTERPRETATION, and
> shells.junction (PR #249) priced that throat — the bill is inherited, unpaid.

WHAT THIS ROUND CORRECTS
────────────────────────
Every earlier wave round in this arc drew ONE wavefront converging on its
antipode and treated the caustic as the interesting event.  A caustic is where
the AMPLITUDE gets large; pair creation is a threshold on an INVARIANT, and the
two are not the same quantity:

    gamma gamma -> e+ e-      s = 2 E1 E2 (1 - cos theta) >= (2m)^2

E is what focusing raises.  Theta is what a collision supplies, and focusing
does not create one.

WHAT IS CHECKED
───────────────
T2  FOCUSING IS NEITHER SUFFICIENT NOR NECESSARY.  Collinear momenta have s = 0
    identically — amplify by 1e12 and it is still zero, so no caustic reaches
    threshold by being bright.  And two beams crossed head-on with NO focusing
    anywhere clear threshold as soon as E >= m.  Both halves, because the
    convenient half alone would be a slogan.  Also recorded: a spherically
    converging front does contain diametrically opposed rays, so its
    self-invariant is not zero — the distinction is INDEPENDENCE of the sources,
    not brightness, and saying otherwise would be the easy overclaim.

T3  THE INVARIANT IS A TRIANGLE IDENTITY.  1 - cos theta = (1 - cos delta)/sin^2 t,
    checked against embedded tangent vectors that never use the law of cosines:
    the crossing point is solved as a linear system and the momenta are
    great-circle tangents in the embedding.  Agrees to 1e-14 on S2 AND on S3 —
    the dimension does not matter, because a geodesic triangle lies in a great
    2-sphere whatever it is embedded in.

T4  THE COLLISION IS HEAD-ON TWICE.  s(t) = 4 E1E2 sin^2(delta/2)/sin^2 t is
    U-SHAPED: maximal (4 E1E2, head-on) at BOTH ends of the crossing window and
    minimal at the equator, where the opening angle is exactly delta.  The
    moment the wavefronts are largest is the moment the invariant is smallest.

T5  SO THE THRESHOLD OPENS TWO WINDOWS, NEVER ONE.  Closed form against a
    40,000-point scan.  Below E = m there is no window at all; at E = m exactly
    the window has zero width (threshold touched only at the two isolated
    head-on instants); the two merge only above E = m/sin(delta/2).

T6  AND ONLY THE FAR WINDOW IS A COLLISION OF INDEPENDENT WAVES.  The near
    window sits on top of the sources — the fronts have not separated yet.  The
    far one is reached after each front has crossed a half-circumference.  That
    is why the second interaction has to be ANTIPODAL, and it is measured as a
    ratio of path lengths rather than asserted.

T7  THE PROJECTED ANGLE IS NOT THE OPENING ANGLE.  The momenta are exact —
    perpendicular to their own front to 1e-15, angle matching the closed form to
    2e-13 — but a figure shows their PROJECTION, and measured off the picture
    the opening angle is wrong by up to 67 degrees, and by up to 56 degrees
    DIFFERENTLY at the two crossing points, whose true opening angle is
    identical.  So the renderer draws the angle in the plane the two momenta
    span, and the arrows on the sphere are labelled as decoration.

T8  THE CROSS-SECTION IS IMPORTED AND CHECKED.  The textbook Breit-Wheeler
    sigma(s), verified against its known peak: beta = 0.701, sqrt(s) = 1.40 (2m),
    sigma = 0.256 sigma_T.  It vanishes at threshold and FALLS at large s — the
    most violent part of the crossing is not the most productive.

T9  ONE UP, ONE DOWN.  The orientation labels cancel and the crossing locus
    supplies exactly the two points needed to carry them — two on S2, a circle
    on S3.  Nothing here derives charge from geometry: the label is carried and
    checked, not produced.

THE LESSON
──────────
The same trap as the previous round, one level up.  Last time a caption said
"collapse" where the geometry said "arrive".  This time the whole arc had been
drawing a caustic and calling it a creation event — and the correction is not
that the drawing was inaccurate but that the quantity was wrong.  Amplitude is
not an invariant.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.pair_creation import (
    measure_focusing_alone_creates_no_invariant_mass,
    measure_only_the_far_window_is_independent,
    measure_the_collision_is_head_on_twice,
    measure_the_cross_section_is_imported_and_checked,
    measure_the_invariant_is_a_triangle_identity,
    measure_the_pair_conserves_orientation,
    measure_the_projected_angle_is_not_the_opening_angle,
    measure_the_threshold_opens_two_windows,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("earlier rounds drew a single wavefront refocusing at its "
                     "antipode and treated the caustic as a particle-creation "
                     "event; is it one, and if not, what does the geometry "
                     "actually require?"),
        "scope": ("kinematics on a fixed round sphere. The Breit-Wheeler "
                  "threshold and cross-section are imported QED; "
                  "rays-as-photons is a correspondence; no rate is computed; "
                  "the crossing chord as a throat is interpretation, and PR "
                  "#249 priced that throat without this round paying it."),
        "pass": True,
    }


def t2_focusing_is_neither_sufficient_nor_necessary() -> dict:
    r = measure_focusing_alone_creates_no_invariant_mass()
    return {"name": "T2_focusing_is_neither_sufficient_nor_necessary", **r,
            "pass": bool(r["focusing_is_not_sufficient"]
                         and r["focusing_is_not_necessary"]
                         and r["threshold_is_at_energy_equal_mass"])}


def t3_the_invariant_is_a_triangle_identity() -> dict:
    """Closed form against embedded tangent vectors, in two dimensions."""
    r = measure_the_invariant_is_a_triangle_identity()
    return {"name": "T3_the_invariant_is_a_triangle_identity", **r,
            "pass": bool(r["the_closed_form_is_confirmed"]
                         and r["and_it_is_dimension_independent"])}


def t4_the_collision_is_head_on_twice() -> dict:
    r = measure_the_collision_is_head_on_twice()
    return {"name": "T4_the_collision_is_head_on_twice", **r,
            "pass": bool(r["both_ends_are_head_on"]
                         and r["the_equator_angle_is_exactly_the_separation"]
                         and r["the_minimum_is_at_the_equator"])}


def t5_the_threshold_opens_two_windows() -> dict:
    r = measure_the_threshold_opens_two_windows()
    return {"name": "T5_the_threshold_opens_two_windows", **r,
            "pass": bool(r["the_scan_agrees_with_the_closed_form"]
                         and r["below_E_equals_m_there_is_no_window"]
                         and r["there_are_exactly_two_windows_in_between"]
                         and r["and_they_merge_above"])}


def t6_only_the_far_window_is_independent() -> dict:
    """Why the second interaction has to be antipodal, as a path length."""
    r = measure_only_the_far_window_is_independent()
    return {"name": "T6_only_the_far_window_is_independent", **r,
            "pass": bool(r["near_collision_is_within_the_source_region"]
                         and r["far_collision_is_past_a_quarter_turn"]
                         and r["both_windows_are_head_on_at_their_outer_edge"])}


def t7_the_projected_angle_is_not_the_opening_angle() -> dict:
    """The arrows in the figure cannot carry the claim; this says by how much."""
    r = measure_the_projected_angle_is_not_the_opening_angle()
    return {"name": "T7_the_projected_angle_is_not_the_opening_angle", **r,
            "pass": bool(r["momenta_are_perpendicular_to_the_sphere"]
                         and r["the_projection_misreads_the_angle"]
                         and r["and_it_misreads_the_two_crossings_differently"])}


def t8_the_cross_section_is_imported_and_checked() -> dict:
    r = measure_the_cross_section_is_imported_and_checked()
    return {"name": "T8_the_cross_section_is_imported_and_checked", **r,
            "pass": bool(r["matches_the_textbook_peak"]
                         and r["zero_at_threshold"]
                         and r["falls_at_large_s"])}


def t9_the_pair_conserves_orientation() -> dict:
    r = measure_the_pair_conserves_orientation()
    return {"name": "T9_the_pair_conserves_orientation", **r,
            "pass": bool(r["the_labels_cancel"]
                         and r["above_threshold_there"]
                         and r["crossing_locus_size_on_S2"] == 2)}


def t10_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T10_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_focusing_is_neither_sufficient_nor_necessary(),
             t3_the_invariant_is_a_triangle_identity(),
             t4_the_collision_is_head_on_twice(),
             t5_the_threshold_opens_two_windows(),
             t6_only_the_far_window_is_independent(),
             t7_the_projected_angle_is_not_the_opening_angle(),
             t8_the_cross_section_is_imported_and_checked(),
             t9_the_pair_conserves_orientation()]
    tests.append(t10_assessment(tests))
    t2, t3, t5, t6 = tests[1], tests[2], tests[4], tests[5]
    t7, t8 = tests[6], tests[7]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_SECOND_INTERACTION_MUST_BE_ANTIPODAL"
        verdict = (
            "PAIR CREATION IS A COLLISION, NOT A FOCUS — AND THE GEOMETRY THEN "
            "FORCES A SECOND INTERACTION AT THE ANTIPODE. A caustic is where "
            "the AMPLITUDE gets large; Breit-Wheeler is a threshold on an "
            "INVARIANT, s = 2 E1E2 (1 - cos theta) >= (2m)^2, and those are "
            "different quantities. E is what focusing raises; theta is what a "
            "collision supplies. Focusing is therefore NEITHER SUFFICIENT — "
            "collinear momenta have s = 0 identically, still zero after "
            f"amplifying by {t2.get('largest_amplification_tried', 0):.0e} — "
            "NOR NECESSARY, since two beams crossed head-on with no focusing "
            "anywhere clear threshold as soon as E >= m. The honest "
            "complication is recorded rather than buried: a spherically "
            "converging front DOES contain diametrically opposed rays, so its "
            "self-invariant is not zero, and the real distinction is "
            "INDEPENDENCE of the sources rather than brightness. Put two "
            "sources a geodesic distance delta apart and the crossing obeys "
            "1 - cos theta = (1 - cos delta)/sin^2 t, an identity of geodesic "
            "triangles verified against embedded tangent vectors that never "
            "use the law of cosines — to "
            f"{t3.get('worst_over_all_dimensions', 0):.1e}, on S2 AND on S3, "
            "because a geodesic triangle lies in a great 2-sphere whatever it "
            "is embedded in. So s(t) = 4 E1E2 sin^2(delta/2)/sin^2 t, which is "
            "U-SHAPED: head-on at BOTH ends of the crossing window and minimal "
            "at the equator, where the opening angle is exactly delta. The "
            "moment the wavefronts are largest is the moment the invariant is "
            "smallest, which is the opposite of the intuition the "
            "single-caustic pictures encouraged. A threshold therefore cuts "
            "TWO DISJOINT WINDOWS AND NEVER ONE, confirmed against a "
            "40,000-point scan; below E = m there is none, at E = m the window "
            "has zero width, and they merge only above "
            f"E = m/sin(delta/2) = {t5.get('merge_energy_over_mass', 0):.3f} m. "
            "AND ONLY THE FAR WINDOW IS A COLLISION OF INDEPENDENT WAVES: the "
            "near one sits on top of the sources, the fronts having travelled "
            f"{t6.get('path_before_near_collision', 0):.3f} against a "
            f"separation of {t6.get('separation_of_the_sources', 0):.2f}, "
            "while the far one is reached only after each front has crossed a "
            "half-circumference — a factor of "
            f"{t6.get('ratio_of_path_lengths', 0):.1f} in path length. That is "
            "why the second interaction has to be antipodal, and it is derived "
            "rather than staged. ONE FURTHER TRAP, CAUGHT: the momenta are "
            "exact, perpendicular to their own wavefront to 1e-15 and matching "
            "the closed form to 2e-13, but a FIGURE shows their projection, "
            "and projection does not preserve angles. Measured off the picture "
            "the opening angle is wrong by up to "
            f"{t7.get('worst_projected_error_deg', 0):.1f} degrees and differs "
            f"by up to {t7.get('worst_disagreement_between_the_two_crossings_deg', 0):.1f} "
            "degrees between the two crossing points, whose true opening angle "
            "is IDENTICAL — so the renderer draws the angle in the plane the "
            "two momenta span and labels the arrows on the sphere as "
            "decoration. WHAT IS IMPORTED, PLAINLY: the Breit-Wheeler "
            "threshold and cross-section are QED. The cross-section is the "
            "textbook closed form and is checked against its known peak — "
            f"beta = {t8.get('peak_beta', 0):.3f}, sqrt(s) = "
            f"{t8.get('peak_sqrt_s_over_2m', 0):.2f} (2m), sigma = "
            f"{t8.get('peak_sigma_over_thomson', 0):.3f} sigma_T — and it "
            "FALLS at large s, so the most violent part of the crossing is not "
            "the most productive. Rays-as-photons is a correspondence, no rate "
            "is computed, the orientation labels are carried rather than "
            "derived, and calling the crossing chord a throat through the bulk "
            "is this program's reading, whose exotic-matter bill was priced in "
            "PR #249 and is inherited here unpaid.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "pair_creation", "tests": tests,
            "verdict_class": verdict_class, "verdict": verdict,
            "generated_utc": datetime.now(timezone.utc).isoformat()}


# ════════════════════════════════════════════════════════════════════════════
def render_markdown(s: dict) -> str:
    lines = [f"# probe — {s['probe']}", "", f"_{s['generated_utc']}_", ""]
    for t in s["tests"]:
        lines.append(f"## {t['name']} — {'PASS' if t['pass'] else 'FAIL'}")
        lines.append("")
        for k, v in t.items():
            if k in ("name", "pass"):
                continue
            if isinstance(v, list):
                lines.append(f"- **{k}**:")
                for row in v:
                    if isinstance(row, dict):
                        lines.append("    - " + ", ".join(
                            f"{a}={_fmt(b)}" for a, b in row.items()))
                    else:
                        lines.append(f"    - {_fmt(row)}")
            else:
                lines.append(f"- **{k}**: {_fmt(v)}")
        lines.append("")
    lines += [f"## verdict — {s['verdict_class']}", "", s["verdict"], ""]
    return "\n".join(lines)


def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, dict):
        return ", ".join(f"{a}={_fmt(b)}" for a, b in v.items())
    return str(v)


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_pair_creation_probe"
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
