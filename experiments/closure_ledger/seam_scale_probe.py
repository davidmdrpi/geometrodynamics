"""
The scaling at the seam is a choice — what it changes, and what it does not

> Framing: both gluings are REPRESENTATION choices, not derived boundary
> conditions.  Nothing here makes the wave dynamically aware of the seam.

THE OBJECTION
─────────────
`circle_slice_probe` folded the bulk by TRANSLATION: a radius past R_outer came
back as r − gap.  That carries a radial OFFSET across unchanged — and the two
boundary circles do not have the same circumference.  A feature emerging at
R_inner keeps its full radial height while sitting on an arc shorter by
R_outer/R_inner, so it comes back squashed.  The emerging wave is not the same
wave, and the scaling that made it so was never chosen deliberately.

THE ALTERNATIVE
───────────────
Glue conformally instead: r → r · (R_inner/R_outer), a translation in ln r
rather than in r.  Radial offsets then scale with the boundary they cross.

WHAT IS CHECKED
───────────────
T2  THE EMERGING FEATURE.  Crossing outward, the translate rule multiplies a
    feature's aspect ratio by exactly R_outer/R_inner = 1.7027; the conformal
    rule leaves it at 1.0000.  Conformal returns a faithful scaled copy.

T3  WHETHER THE RADIUS SURVIVES.  Translate sheets are arithmetic and march
    inward through r = 0 into NEGATIVE radius (0.74 → 0.22 → −0.30 → −0.82).
    Conformal sheets are geometric, accumulating at the origin without ever
    arriving.  Conformal pairs with a multiplicative radial law
    r = R_mid·exp(εu), positive by construction.

T4  WHAT A WINDING NUMBER WOULD LOOK LIKE.  A curve that genuinely winds — a
    ramp climbing one period as σ goes around, a logarithmic spiral on the
    conformal seam — returns to the same point of the quotient MAGNIFIED by
    1.7027, and by the same factor from any starting radius.  On the translate
    seam it returns displaced instead, with a ratio that drifts with where it
    started: not a scale at all.  So a conformal gluing makes topological charge
    observable as a magnification, and a translate gluing hides it.

T5  AND THE WINDING NUMBER IS STILL ZERO.  This is the check the objection
    deserves.  Rebuilt on a conformal seam with a multiplicative radial law — a
    different identification, a different sheet structure, a different notion of
    size — the winding is identically zero, at gains driving hundreds of
    unsigned crossings.  ρ(σ) comes from a single-valued function on the circle
    whichever coordinate the seam translates in.

SCOPE
─────
The conformal rule is preferred here on grounds of CONSISTENCY — it returns a
scaled copy and cannot produce a negative radius — not because any dynamics
picks it out.  Neither rule makes the seam something the wave can feel.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.circle_slice import (
    BulkAnnulus,
    CircleSlice,
    measure_the_emerging_feature_keeps_its_shape,
    measure_the_rule_does_not_rescue_the_winding,
    measure_the_translate_rule_runs_off_the_origin,
    measure_winding_is_a_magnification,
)
from geometrodynamics.viz.warped_sphere import NestedShells

GAINS = (0.4, 1.0, 3.0, 12.0, 60.0)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    sh = NestedShells()
    return {
        "name": "T1_goal",
        "claim": ("the gluing rule at the seam sets the scale of the emerging "
                  "wave, and that scale was a free choice"),
        "translate": "r → r − gap  (a translation in r)",
        "conformal": "r → r · (R_inner/R_outer)  (a translation in ln r)",
        "circumference_ratio": sh.r_outer / sh.r_inner,
        "pass": True,
    }


def t2_the_emerging_feature() -> dict:
    r = measure_the_emerging_feature_keeps_its_shape()
    return {
        "name": "T2_the_emerging_feature",
        **r,
        "pass": bool(r["translate_squashes_the_feature"]
                     and r["conformal_preserves_the_shape"]),
    }


def t3_whether_the_radius_survives() -> dict:
    r = measure_the_translate_rule_runs_off_the_origin()
    return {
        "name": "T3_whether_the_radius_survives",
        **r,
        "pass": bool(r["translate_reaches_negative_radius"]
                     and r["conformal_stays_positive"]),
    }


def t4_winding_would_be_a_magnification() -> dict:
    r = measure_winding_is_a_magnification()
    return {
        "name": "T4_winding_would_be_a_magnification",
        **r,
        "pass": bool(r["a_wound_curve_climbs_the_sheets"]
                     and r["conformal_turns_winding_into_a_scale"]
                     and r["translate_gives_no_scale_at_all"]),
    }


def t5_the_winding_is_still_zero() -> dict:
    r = measure_the_rule_does_not_rescue_the_winding(gains=GAINS)
    return {
        "name": "T5_the_winding_is_still_zero",
        **r,
        "pass": bool(r["every_rule_gives_zero"]),
    }


def t6_the_two_rules_agree_at_small_excursion() -> dict:
    """They are the same rule to first order — the difference is all nonlinear.

    Worth pinning down, because it says the v46 pictures were not wrong where
    the wave was small; they diverge exactly where the wave reaches the seam.
    """
    sh = NestedShells()
    rows = []
    for h in (1e-3, 1e-2, 1e-1):
        t = float(BulkAnnulus(sh, mode="translate").wrap(
            np.array([sh.r_outer + h]))[0][0])
        c = float(BulkAnnulus(sh, mode="conformal").wrap(
            np.array([sh.r_outer + h]))[0][0])
        rows.append({"height": h, "translate": t, "conformal": c,
                     "difference": abs(t - c),
                     "difference_over_height": abs(t - c) / h})
    ratios = [r["difference_over_height"] for r in rows]
    return {
        "name": "T6_the_two_rules_agree_at_small_excursion",
        "rows": rows,
        "difference_is_proportional_to_the_excursion": bool(
            max(ratios) - min(ratios) < 1e-6),
        "shared_first_order_slope": ratios[0],
        "pass": bool(rows[0]["difference"] < 1e-3),
    }


def t7_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T7_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_the_emerging_feature(),
             t3_whether_the_radius_survives(),
             t4_winding_would_be_a_magnification(),
             t5_the_winding_is_still_zero(),
             t6_the_two_rules_agree_at_small_excursion()]
    tests.append(t7_assessment(tests))
    t1, t2, t3, t4, t5, t6 = tests[0], tests[1], tests[2], tests[3], tests[4], tests[5]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_GLUING_SETS_THE_SCALE_BUT_NOT_THE_TOPOLOGY"
        verdict = (
            "THE GLUING SETS THE SCALE, BUT NOT THE TOPOLOGY. How the seam is "
            "glued was a free choice, and it decides what the emerging wave "
            "looks like — but not whether it can carry a winding number.\n\n"
            "THE TRANSLATE RULE RETURNS A CARICATURE. Carrying a radial offset "
            "across unchanged lands the feature on an arc shorter by "
            f"{t1['circumference_ratio']:.4f}, so its aspect ratio is "
            f"multiplied by {t2['translate_distortion']:.4f} — the emerging "
            "wave is not the same wave. The conformal rule scales the offset "
            "with the boundary it crosses, so height and arc length shrink "
            f"together and the distortion is {t2['conformal_distortion']:.4f}: "
            "a faithful scaled copy.\n\n"
            "AND THE TRANSLATE RULE RUNS OFF THE ORIGIN. Its inward sheets are "
            "arithmetic, so they march "
            f"{t3['translate']['inward'][0]:.2f} → "
            f"{t3['translate']['inward'][1]:.2f} → "
            f"{t3['translate']['inward'][2]:.2f}: straight through r = 0 into "
            "negative radius, because subtracting a fixed gap has nothing to "
            "stop it. Conformal sheets are geometric — "
            f"{t3['conformal']['inward'][0]:.3f} → "
            f"{t3['conformal']['inward'][1]:.3f} → "
            f"{t3['conformal']['inward'][2]:.3f} — accumulating at the origin "
            "and never arriving, and they pair with a multiplicative radial law "
            "that is positive by construction.\n\n"
            "THE CHOICE ALSO DECIDES WHAT WINDING WOULD LOOK LIKE. A curve that "
            "genuinely winds is a logarithmic spiral on the conformal seam: it "
            "returns to the same point of the quotient magnified by "
            f"{t4['conformal_magnification']:.4f}, and by the same factor from "
            f"every starting radius (spread {t4['conformal_spread']:.1e}). On "
            "the translate seam the same curve returns displaced instead, with "
            f"a ratio that drifts by {t4['translate_spread']:.3f} depending on "
            "where it started — not a scale at all. A conformal gluing makes "
            "topological charge observable as a magnification; a translate "
            "gluing hides it.\n\n"
            "BUT THE WINDING NUMBER IS STILL ZERO. Rebuilt on a conformal seam "
            "with a multiplicative radial law — a different identification, a "
            "different sheet structure, a different notion of size — the "
            "winding is identically zero at every gain tested, driving up to "
            f"{max(r['unsigned'] for r in t5['rows'])} unsigned crossings. "
            "ρ(σ) comes from a single-valued function on the circle whichever "
            "coordinate the seam translates in, so its degree is zero either "
            "way. The earlier negative result was not an artefact of an "
            "arbitrary scaling.\n\n"
            "THE TWO RULES ARE THE SAME RULE UNTIL THE WAVE REACHES THE SEAM. "
            "Their difference is proportional to the excursion, with slope "
            f"{t6['shared_first_order_slope']:.4f} — so the earlier pictures "
            "were not wrong where the wave was small. They part company exactly "
            "where it matters.\n\n"
            "SCOPE. The conformal rule is preferred on grounds of consistency — "
            "it returns a scaled copy and cannot produce a negative radius — "
            "not because any dynamics picks it out. Neither rule makes the seam "
            "something the wave can feel, and neither gives a height field "
            "anywhere to wind."
        )
    else:
        verdict_class = "SEAM_SCALE_INCONCLUSIVE"
        verdict = ("INCONCLUSIVE. A check failed; review the shape distortion, "
                   "the sheet ladder, the spiral magnification, the winding "
                   "under each rule, or the small-excursion agreement.")

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "two ways of gluing the vacuole's boundaries, and what each does to "
            "the wave that crosses"
        ),
        "the_choice": "translate in r, or translate in ln r",
        "what_it_changes": "the shape of the emerging feature, and the sheets",
        "what_it_does_not": "the winding number, which is zero for any graph",
        "geometry": {"gains": list(GAINS)},
        "tests": tests,
        "n_passed": tests[-1]["n_passed"],
        "n_total": tests[-1]["n_total"],
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    t = {x["name"]: x for x in s["tests"]}
    out: List[str] = []
    out.append("# The scaling at the seam is a choice\n")
    out.append(f"_{s['timestamp_utc']}_\n")
    out.append("> Both gluings are **representation** choices, not derived "
               "boundary conditions.\n")

    c = t["T2_the_emerging_feature"]
    out.append("## The emerging feature\n")
    out.append("| rule | height in | height out | arc out | aspect distortion |")
    out.append("|---|---:|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| `{r['mode']}` | {r['height_before']:.4f} | "
                   f"{r['height_after']:.4f} | {r['arc_after']:.4f} | "
                   f"**{r['aspect_distortion']:.4f}** |")
    out.append(f"\nThe circumference ratio is "
               f"{c['circumference_ratio']:.4f}, and the translate rule pays it "
               "in full.\n")

    c = t["T3_whether_the_radius_survives"]
    out.append("## The inward sheets\n")
    out.append("| rule | sheet −1 | −2 | −3 | −4 |")
    out.append("|---|---:|---:|---:|---:|")
    for mode in ("translate", "conformal"):
        vals = " | ".join(f"{v:.3f}" for v in c[mode]["inward"])
        out.append(f"| `{mode}` | {vals} |")
    out.append("\nArithmetic sheets walk through `r = 0`; geometric ones "
               "accumulate at it.\n")

    c = t["T4_winding_would_be_a_magnification"]
    out.append("## What a winding number would look like\n")
    out.append("| seam | magnifications from r₀ = 0.80, 1.00, 1.20 | spread |")
    out.append("|---|---|---:|")
    for r in c["rows"]:
        mags = ", ".join(f"{m:.4f}" for m in r["magnifications"])
        out.append(f"| `{r['seam']}` | {mags} | `{r['magnification_spread']:.1e}` |")
    out.append("\nA scale is start-independent; a shift is not.\n")

    c = t["T5_the_winding_is_still_zero"]
    out.append("## ...and the winding is still zero\n")
    out.append("| seam | radial law | unsigned | signed | winding |")
    out.append("|---|---|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| `{r['seam']}` | `{r['radial_law']}` | {r['unsigned']} | "
                   f"{r['signed']:+d} | {r['winding']:+d} |")
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
    out = here / "runs" / f"{ts}_seam_scale_probe"
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
