"""
What a continuous trace-free deformation does at its own antipodal refocus

> Framing: a faithful REPRESENTATION of a spin-2 field in the ℝ³ embedding of a
> fixed S² — not backreaction, and not linearised GR on a spacetime.

THE QUESTION
────────────
`embedded_wave_probe` established the projection: one potential W carries both
the shape ρ = −½ΔW and the shear h_ab = [2∇₍ₐ∇_b₎W]^TF, so the drawn surface
HAS the solved metric perturbation as its own geometry.  Every point of the
sphere then runs its own principal-strain history, and on a compact S² those
histories all reconverge at the antipode at t = π.

So: what happens to a continuous trace-free deformation where they refocus?

WHAT IS CHECKED
───────────────
T2  THE FOCUS IS A NODE.  h = sin²d·q vanishes at the poles for every q, so
    ḣ does too, and the effective density ∝ ḣ_ab ḣ^ab cannot pile onto the
    focal POINT.  Measured on the antipode against the peak.

T3  IT PILES INTO A RING instead, whose radius tracks the pulse width.  The
    ratio is measured across five widths rather than fitted once.

T4  THE AMPLIFICATION IS NOT A SPIN-2 EFFECT.  It is tempting to read the
    finite ~2× peak as the spin protecting itself from a singularity.  A
    SCALAR pulse refocused the same way amplifies by the same O(1) factor,
    and neither runs away as the pulse narrows.  Launch and focus are the
    same situation on a sphere.  What belongs to the spin is T2 and T3.

T5  A MATERIAL PATCH REPORTS THE EIGENVECTOR — in the limit of a point.  Near
    the focus the shear turns over on the scale of the pulse, so a patch wide
    enough to straddle the focal ring averages a sign change; shrinking it
    recovers the local stretch axis.  The convergence is the check that the
    disagreement is the patch's SIZE and not the construction.

T6  THE PATCH CHANGES SHAPE WITHOUT CHANGING SIZE, at a gain where the
    first-order statement is the whole statement.

T7  AND THE AREA LAW FAILS FIRST AT THE FOCUS.  Push to the display gain and
    the second-order residual — which carries the LOCAL GRADIENT of the field
    — grows by an order of magnitude more on the focal ring than at
    mid-latitude, because there the gradient scale is the pulse width rather
    than the wavelength.  The fitted exponent confirms it is O(ε²), i.e. the
    linear description running out rather than a numerical defect.  This is
    the honest boundary of what h_ab alone can say about a refocusing wave.

SCOPE
─────
No singularity forms here and none can: this is a linear field on a fixed
round background, with no backreaction and no bulk crossing rule.  T7 locates
where such a rule would have to act; it does not supply one.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.embedded_wave import (
    EmbeddedTidalSurface,
    measure_patch_area_invariance,
    measure_patch_shape_history,
    measure_patch_size_convergence,
    measure_where_the_area_law_fails,
)
from geometrodynamics.viz.spin2_tidal import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    Spin2WaveSim,
    TidalField,
    measure_amplification_is_not_protection,
    measure_focal_energy,
)

N_RADIAL = 1200
PULSE_WIDTH = 0.18
SMALL_GAIN = 1e-2
WIDTHS = (0.24, 0.18, 0.12, 0.09, 0.06)


def _surface() -> EmbeddedTidalSurface:
    return EmbeddedTidalSurface(
        sim=Spin2WaveSim(n=N_RADIAL, pulse_width=PULSE_WIDTH),
        n_theta=121, n_phi=181)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "claim": ("determine what a continuous trace-free deformation does "
                  "where all of its principal-strain histories refocus"),
        "energy_density": "ρ_E ∝ ḣ_ab ḣ^ab = 2(ḣ₊² + ḣ_ˣ²)",
        "why_a_node": "h = sin²d·q vanishes at both poles for every q",
        "the_patch": "a disc of labelled particles, carried by the displacement",
        "pass": True,
    }


def t2_the_focus_is_a_node() -> dict:
    r = measure_focal_energy(TidalField(sim=Spin2WaveSim(
        n=N_RADIAL, pulse_width=PULSE_WIDTH)))
    return {
        "name": "T2_the_focus_is_a_node",
        **{k: r[k] for k in ("launch_peak_density", "focal_peak_density",
                             "amplification", "focal_time", "focal_distance",
                             "ring_radius", "density_on_the_antipode",
                             "antipode_over_peak", "invariant_drift",
                             "kinetic_swing")},
        "pass": bool(r["antipode_over_peak"] < 1e-4
                     and abs(r["invariant_drift"]) < 1e-8),
    }


def t3_it_piles_into_a_ring() -> dict:
    rows = []
    for w in WIDTHS:
        r = measure_focal_energy(TidalField(sim=Spin2WaveSim(
            n=N_RADIAL, pulse_width=float(w))))
        rows.append({"pulse_width": float(w),
                     "ring_radius": r["ring_radius"],
                     "ratio": r["ring_radius"] / float(w),
                     "amplification": r["amplification"]})
    ratios = [x["ratio"] for x in rows]
    spread = max(ratios) - min(ratios)
    return {
        "name": "T3_it_piles_into_a_ring",
        "rows": rows,
        "mean_ratio": sum(ratios) / len(ratios),
        "ratio_spread": spread,
        "tracks_the_pulse_width": bool(spread < 0.12),
        "pass": bool(spread < 0.12 and min(ratios) > 0.7),
    }


def t4_amplification_is_not_protection() -> dict:
    r = measure_amplification_is_not_protection(widths=WIDTHS)
    return {
        "name": "T4_amplification_is_not_protection",
        **r,
        "pass": bool(r["both_amplify_by_order_one"]
                     and r["neither_runs_away_as_the_pulse_narrows"]
                     and r["amplification_is_not_a_spin_2_effect"]
                     and r["but_the_focal_node_is"]),
    }


def t5_the_patch_reports_the_eigenvector() -> dict:
    focal = math.pi - 0.94 * PULSE_WIDTH
    conv = measure_patch_size_convergence(centre_distance=focal)
    smooth = measure_patch_shape_history(_surface(), centre_distance=1.20)
    return {
        "name": "T5_the_patch_reports_the_eigenvector",
        "focal_distance": focal,
        "rows": conv["rows"],
        "smallest_patch_alignment": conv["smallest_patch_alignment"],
        "worst_alignment": conv["worst_alignment"],
        "improves_as_the_patch_shrinks": conv["improves_as_the_patch_shrinks"],
        "smooth_alignment": smooth["long_axis_alignment"],
        "smooth_max_aspect_ratio": smooth["max_aspect_ratio"],
        "pass": bool(conv["converges_to_the_eigenvector"]
                     and conv["improves_as_the_patch_shrinks"]
                     and smooth["aligns_with_the_stretch_axis"]),
    }


def t6_shape_without_size() -> dict:
    r = measure_patch_area_invariance(gain=SMALL_GAIN)
    return {
        "name": "T6_shape_without_size",
        **{k: r[k] for k in ("gain", "reference_area", "relative_area_swing",
                             "max_aspect_ratio",
                             "max_aspect_ratio_at_display_gain",
                             "display_gain")},
        "pass": bool(r["area_is_invariant"]),
    }


def t7_the_area_law_fails_first_at_the_focus() -> dict:
    r = measure_where_the_area_law_fails(_surface())
    return {
        "name": "T7_the_area_law_fails_first_at_the_focus",
        **{k: r[k] for k in (
            "smooth_distance", "focal_distance", "display_gain",
            "smooth_worst_area_change", "smooth_max_aspect_ratio",
            "focal_worst_area_change", "focal_max_aspect_ratio",
            "gains", "focal_area_residuals", "residual_exponent",
            "focal_over_smooth_area_change", "focal_over_smooth_distortion")},
        "pass": bool(r["residual_is_second_order"]
                     and r["the_area_law_fails_first_at_the_focus"]
                     and r["the_focus_distorts_harder"]),
    }


def t8_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T8_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_the_focus_is_a_node(), t3_it_piles_into_a_ring(),
             t4_amplification_is_not_protection(),
             t5_the_patch_reports_the_eigenvector(), t6_shape_without_size(),
             t7_the_area_law_fails_first_at_the_focus()]
    tests.append(t8_assessment(tests))
    t2, t3, t4, t5, t6, t7 = tests[1], tests[2], tests[3], tests[4], tests[5], tests[6]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_STRAINS_REFOCUS_ON_A_RING"
        verdict = (
            "THE STRAINS REFOCUS ON A RING, NEVER ON THE POINT. Every "
            "principal-strain history on the sphere reconverges at the "
            "antipode, and what they build there is an annulus with a hole in "
            "the middle.\n\n"
            "THE FOCUS IS A NODE. Because h = sin²d·q vanishes at both poles "
            "for every q, so does ḣ, and the effective density ∝ ḣ_ab ḣ^ab "
            "measured ON the antipode is "
            f"{t2['antipode_over_peak']:.1e} of the peak. There is nowhere for "
            "a spin-2 field to sit at its own focus — the same fact that "
            "forbids a spin-2 point SOURCE, seen at the other end of the "
            "trip.\n\n"
            "IT PILES INTO A RING WHOSE RADIUS IS THE PULSE. Across pulse "
            f"widths {WIDTHS[0]} down to {WIDTHS[-1]} the focal ring radius "
            f"stays at {t3['mean_ratio']:.3f} × the width, spread "
            f"{t3['ratio_spread']:.3f}. The hole is not a numerical floor; it "
            "is the wave's own length scale.\n\n"
            "AND THE AMPLIFICATION IS NOT A SPIN-2 EFFECT. The peak density "
            f"grows by {t4['tensor_amplification_range'][0]:.2f}–"
            f"{t4['tensor_amplification_range'][1]:.2f}× — finite, and "
            "tempting to read as the spin protecting itself from a "
            "singularity. It is not. A scalar pulse refocused the same way "
            f"amplifies by {t4['scalar_amplification_range'][0]:.2f}–"
            f"{t4['scalar_amplification_range'][1]:.2f}×, and neither runs "
            "away as the pulse is narrowed, because launch and focus are "
            "geometrically the same situation on a sphere. What belongs to "
            "the spin is the node and the ring, not the factor.\n\n"
            "THE PATCH IS THE PICTURE. A disc of labelled particles carried "
            "by the displacement reports the local stretch axis to "
            f"{t5['smooth_alignment']:.3f} where the field is smooth, and to "
            f"{t5['smallest_patch_alignment']:.4f} on the focal ring once it "
            "is small enough not to straddle the sign change there — the "
            "convergence being the check that the disagreement is the patch's "
            "size rather than the construction. At a gain where the "
            "first-order statement is the whole statement it changes shape "
            f"without changing size, holding its area to "
            f"{t6['relative_area_swing']:.1e}.\n\n"
            "THE AREA LAW FAILS FIRST, AND HARDEST, AT THE FOCUS. Pushed to "
            f"the display gain ε = {t7['display_gain']:.2f}, the same patch "
            f"moves its area by {100 * t7['smooth_worst_area_change']:.1f}% at "
            f"mid-latitude and {100 * t7['focal_worst_area_change']:.1f}% on "
            f"the focal ring — {t7['focal_over_smooth_area_change']:.0f}× more "
            "— while distorting "
            f"{t7['focal_over_smooth_distortion']:.1f}× harder. The residual "
            f"scales as ε^{t7['residual_exponent']:.2f}, so this is the "
            "second-order term carrying the LOCAL GRADIENT of the field: away "
            "from the focus that scale is the wavelength, on the focal ring it "
            "is the pulse width. The linear description runs out exactly where "
            "the wave reconverges.\n\n"
            "SCOPE. No singularity forms here and none can: a linear field on "
            "a fixed round background, with no backreaction and no bulk "
            "crossing rule. What this establishes is the geometry such a rule "
            "would have to act on — an annulus of finite radius with a node at "
            "its centre — and the amplitude at which the linear description "
            "stops being trustworthy. It does not supply the rule."
        )
    else:
        verdict_class = "FOCAL_REFOCUS_INCONCLUSIVE"
        verdict = ("INCONCLUSIVE. A check failed; review the focal node, the "
                   "ring scaling, the scalar comparison, the patch "
                   "convergence, the area invariance, or the second-order "
                   "residual.")

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "the antipodal refocus of a continuous trace-free deformation of "
            "the embedded sphere"
        ),
        "the_node": "h = sin²d·q vanishes at both poles, so ḣ does too",
        "the_ring": "focal radius ≈ 0.94 × pulse width",
        "the_caveat": "the O(1) amplification is shared with the scalar",
        "the_boundary": "the O(ε²) area residual is largest on the focal ring",
        "geometry": {
            "n_radial": N_RADIAL, "pulse_width": PULSE_WIDTH,
            "small_gain": SMALL_GAIN, "widths": list(WIDTHS),
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
    out.append("# The antipodal refocus of a trace-free deformation\n")
    out.append(f"_{s['timestamp_utc']}_\n")
    out.append("> A faithful **representation** of a spin-2 field in the ℝ³ "
               "embedding of a fixed `S²` — not backreaction.\n")

    c = t["T2_the_focus_is_a_node"]
    out.append("## The focus is a node\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| peak density at the focus | {c['focal_peak_density']:.6f} |")
    out.append(f"| amplification | {c['amplification']:.4f}× |")
    out.append(f"| focal time | {c['focal_time']:.4f} |")
    out.append(f"| ring radius | {c['ring_radius']:.4f} |")
    out.append(f"| density **on** the antipode | "
               f"`{c['antipode_over_peak']:.1e}` of peak |")
    out.append(f"| conserved invariant drift | `{c['invariant_drift']:.1e}` |\n")
    out.append("`∫ρ_E dA` is the *kinetic* half and oscillates against the "
               f"gradient half (swing `{c['kinetic_swing']:.2f}`), so the "
               "solver's invariant is the conservation check, not that.\n")

    c = t["T3_it_piles_into_a_ring"]
    out.append("## The ring radius is the pulse width\n")
    out.append("| pulse width | ring radius | ratio | amplification |")
    out.append("|---:|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['pulse_width']:.2f} | {r['ring_radius']:.4f} | "
                   f"{r['ratio']:.3f} | {r['amplification']:.3f}× |")
    out.append(f"\nMean ratio {c['mean_ratio']:.3f}, spread "
               f"{c['ratio_spread']:.3f}.\n")

    c = t["T4_amplification_is_not_protection"]
    out.append("## ...and the amplification is *not* a spin-2 effect\n")
    out.append("| pulse width | tensor | scalar |")
    out.append("|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['pulse_width']:.2f} | "
                   f"{r['tensor_amplification']:.3f}× | "
                   f"{r['scalar_amplification']:.3f}× |")
    out.append("\nBoth are `O(1)` and neither runs away as the pulse narrows. "
               "What belongs to the spin is the **node** and the **ring**.\n")

    c = t["T5_the_patch_reports_the_eigenvector"]
    out.append("## The patch reports the eigenvector, in the limit of a point\n")
    out.append("| patch radius | aspect ratio | alignment |")
    out.append("|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['radius']:.2f} | {r['aspect_ratio']:.4f} | "
                   f"{r['alignment']:.4f} |")
    out.append(f"\nAt mid-latitude the alignment is "
               f"{c['smooth_alignment']:.3f} at any of these sizes; only near "
               "the focus does the patch straddle a sign change.\n")

    c = t["T7_the_area_law_fails_first_at_the_focus"]
    out.append("## Where the area law runs out\n")
    out.append("| quantity | mid-latitude | focal ring |")
    out.append("|---|---:|---:|")
    out.append(f"| worst area change | "
               f"{100 * c['smooth_worst_area_change']:.2f}% | "
               f"{100 * c['focal_worst_area_change']:.2f}% |")
    out.append(f"| max aspect ratio | {c['smooth_max_aspect_ratio']:.4f} | "
               f"{c['focal_max_aspect_ratio']:.4f} |")
    out.append(f"\nFitted residual exponent `ε^{c['residual_exponent']:.2f}` — "
               "second order, at display gain "
               f"`{c['display_gain']:.2f}`.\n")

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
    out = here / "runs" / f"{ts}_focal_refocus_probe"
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
