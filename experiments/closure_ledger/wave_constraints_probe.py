"""
What a drawn wave has to obey, and what a height field cannot

> Framing: this is an audit of the REPRESENTATION used by the circle-slice work,
> against constraints the physics it stands in for actually imposes.

THREE COMPLAINTS, AND THEY DO NOT ALL LAND THE SAME WAY
───────────────────────────────────────────────────────
1. "The antipode moves before the wavefront arrives."
2. "Height should increase as the ring compresses; it is not shown."
3. "So a deforming surface is not the right representation."

The first is aimed one step away from a real defect.  The second is a fair
complaint about the picture and a wrong expectation about the geometry.  The
third is right, and the sharpest reason for it is measurable.

WHAT IS CHECKED
───────────────
T2  THE FRONT IS CAUSAL.  The antipode sits at 5e-133 at launch, 6e-40 at
    t = 1.4 and 7e-17 at t = 2.0, reaching O(1) only at t ≈ π.  Nothing outruns
    the front.  The early motion is not in the solve.

T3  BUT A CONSTANT NEVER LEAVES.  A Gaussian of width w carries a monopole
    w²/4 (0.008056 measured against 0.008100 predicted), and ℓ = 0 has ω = 0,
    so it never propagates, radiates or decays.  It is not an early response —
    ahead of the front the higher modes cancel it exactly — it is a PERMANENT
    displacement.  Time-averaged displacement is 0.0081 at every point, and the
    quietest instant of a whole run still leaves max|u| = 0.094.  The surface
    never returns home.

T4  AND THE MODE IS FORBIDDEN ANYWAY.  Electromagnetism has no monopole
    radiation; gravity has none at ℓ = 0 or ℓ = 1.  The offending mode is
    exactly the one real radiation cannot have.  The spin-2 field of
    spin2_tidal is already clean: h = sin²d·q kills ℓ = 0 and ℓ = 1.

T5  THE FIX HAS TO STAY INSIDE THE PULSE.  Subtracting the mean moves the
    problem rather than removing it — it leaves −w²/4 at the antipode.
    Cancelling with a wider Gaussian fixes the mean and costs four orders of
    far-side quiet, because the corrector's tail is fatter than the pulse.
    Built from compactly supported bumps both conditions hold at once, and the
    far side is quieter than the original: exactly zero rather than merely
    small.

T6  THE RING DOES CONCENTRATE.  A²·(circumference) is conserved along the
    converging front to 8.4%, which is A ∝ 1/√(sin d).  The physics is right.

T7  BUT YOU CANNOT SEE IT, AND IT CANNOT WIN.  1/√(sin d) is flat across the
    middle: 72% of the trip stays within 50% of the equatorial minimum, and the
    growth is all in the last few percent.  Worse for the expectation, on a
    compact surface the launch is itself a focus — the focal height is 0.93× the
    LAUNCH height, so the focus never beats the source.  Unbounded growth
    belongs to an open geometry.

T8  HEIGHT IS NOT ENERGY.  The energy density is u̇² + |∇u|², so the constant
    offset displaces every point of the surface and carries exactly zero
    energy.  The most global feature of the drawn shape is invisible to the
    physics — which is the measurable core of complaint 3.

SCOPE
─────
This audits one representation; it does not propose its replacement.  What it
establishes is which of the objections are about the solve (none), which are
about the initial data (the monopole, fixable), and which are about drawing a
field as a displacement at all (height-is-not-energy, not fixable).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.wave_constraints import (
    measure_a_localised_correction_is_the_only_one_that_works,
    measure_launch_and_focus_are_the_same,
    measure_the_clean_source_leaves_the_antipode_alone,
    measure_the_constant_monopole,
    measure_the_front_is_causal,
    measure_the_growth_is_invisible_for_most_of_the_trip,
    measure_the_monopole_carries_no_energy,
    measure_the_offset_never_leaves,
    measure_the_ring_bookkeeping,
)

PULSE_WIDTH = 0.18


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "claim": ("audit the height-field representation against the "
                  "constraints the physics it stands in for imposes"),
        "complaints": ["the antipode moves early",
                       "the ring's growth is not shown",
                       "a deforming surface is the wrong representation"],
        "pass": True,
    }


def t2_the_front_is_causal() -> dict:
    r = measure_the_front_is_causal(pulse_width=PULSE_WIDTH)
    return {"name": "T2_the_front_is_causal", **r,
            "pass": bool(r["nothing_outruns_the_front"]
                         and r["the_antipode_stays_dark_until_the_front"])}


def t3_a_constant_never_leaves() -> dict:
    m = measure_the_constant_monopole(pulse_width=PULSE_WIDTH)
    o = measure_the_offset_never_leaves(pulse_width=PULSE_WIDTH, frames=400)
    return {"name": "T3_a_constant_never_leaves", **m,
            "time_averaged_displacement": o["time_averaged_displacement"],
            "far_side_average": o["far_side_average"],
            "quietest_moment": o["quietest_moment"],
            "pass": bool(m["the_monopole_is_w_squared_over_four"]
                         and m["it_never_changes"]
                         and o["every_point_is_offset"]
                         and o["the_surface_never_returns_home"])}


def t4_the_mode_is_forbidden() -> dict:
    """The checkable form: does the field admit a zero-frequency mode at all?

    "A spin-2 field has no ``ℓ = 0`` or ``ℓ = 1``" is a statement about the
    *spin-weighted* decomposition, not about ``∫h dA`` — ``h`` is a tensor
    component in a frame, not a scalar, so its plain sphere-average has no
    reason to vanish and measuring that would be measuring the wrong thing.

    What matters here is the property that lets the scalar keep a permanent
    offset: ``ℓ = 0`` has ``ω = 0``, so a DC component is frozen in.  The spin-2
    operator's spectrum starts at ``ℓ = 2``, so it has no zero-frequency mode and
    can carry no DC at all.  That is measured directly, as the time-averaged
    value at fixed points over several returns.
    """
    from geometrodynamics.viz.spin2_tidal import Spin2WaveSim
    from geometrodynamics.viz.throat_wavefront import BareSphereSim

    t_end, frames = 8.0 * math.pi, 1200
    probes = np.array([0.6, 1.4, 2.2, 2.9])

    tensor = Spin2WaveSim(n=1200, pulse_width=PULSE_WIDTH)
    tensor.reset()
    acc_t = np.zeros(len(probes))
    peak = 0.0
    for i in range(frames):
        tensor.advance_to((i + 1) * t_end / frames)
        acc_t += np.interp(probes, tensor.d, tensor.h)
        peak = max(peak, float(np.max(np.abs(tensor.h))))

    scalar = BareSphereSim(n_theta=8, n_phi=8, pulse_width=PULSE_WIDTH,
                           n_radial=1800)
    scalar.reset()
    acc_s = np.zeros(len(probes))
    for i in range(frames):
        scalar.advance_to((i + 1) * t_end / frames)
        acc_s += scalar.field_at_distance(probes)

    dc_t = np.abs(acc_t / frames)
    dc_s = np.abs(acc_s / frames)
    return {
        "name": "T4_the_mode_is_forbidden",
        "em_forbids": "monopole radiation",
        "gravity_forbids": "monopole and dipole radiation",
        "probe_distances": probes.tolist(),
        "scalar_dc": dc_s.tolist(),
        "spin2_dc": dc_t.tolist(),
        "spin2_peak": peak,
        "scalar_worst_dc": float(np.max(dc_s)),
        "spin2_worst_dc": float(np.max(dc_t)),
        "ratio": float(np.max(dc_s) / max(np.max(dc_t), 1e-30)),
        "the_scalar_holds_a_dc_offset": bool(np.max(dc_s) > 1e-3),
        "the_spin_two_field_cannot": bool(np.max(dc_t) < 1e-3 * peak),
        "pass": bool(np.max(dc_s) > 1e-3 and np.max(dc_t) < 1e-3 * peak),
    }


def t5_the_fix_stays_inside_the_pulse() -> dict:
    a = measure_a_localised_correction_is_the_only_one_that_works(
        pulse_width=PULSE_WIDTH)
    b = measure_the_clean_source_leaves_the_antipode_alone(
        pulse_width=PULSE_WIDTH)
    return {"name": "T5_the_fix_stays_inside_the_pulse", **a,
            "launch_rows": b["rows"],
            "compact_monopole": b["compact_monopole"],
            "compact_far_side": b["compact_far_side"],
            "wide_gaussian_far_side": b["wide_gaussian_far_side"],
            "pass": bool(a["subtracting_the_mean_just_moves_it"]
                         and a["the_localised_source_is_clean"]
                         and b["a_wider_corrector_costs_far_side_quiet"]
                         and b["the_compact_source_wins_on_both"])}


def t6_the_ring_concentrates() -> dict:
    r = measure_the_ring_bookkeeping(pulse_width=PULSE_WIDTH)
    return {"name": "T6_the_ring_concentrates", **r,
            "pass": bool(r["energy_on_the_ring_is_conserved"])}


def t7_but_you_cannot_see_it() -> dict:
    g = measure_the_growth_is_invisible_for_most_of_the_trip(
        pulse_width=PULSE_WIDTH)
    l = measure_launch_and_focus_are_the_same()
    return {"name": "T7_but_you_cannot_see_it", **g,
            "launch_rows": l["rows"], "mean_focal_ratio": l["mean_ratio"],
            "pass": bool(g["most_of_the_trip_shows_nothing"]
                         and l["the_focus_does_not_beat_the_launch"])}


def t8_height_is_not_energy() -> dict:
    r = measure_the_monopole_carries_no_energy(pulse_width=PULSE_WIDTH)
    return {"name": "T8_height_is_not_energy", **r,
            "pass": bool(r["it_carries_no_energy"])}


def t9_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T9_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_the_front_is_causal(), t3_a_constant_never_leaves(),
             t4_the_mode_is_forbidden(), t5_the_fix_stays_inside_the_pulse(),
             t6_the_ring_concentrates(), t7_but_you_cannot_see_it(),
             t8_height_is_not_energy()]
    tests.append(t9_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8 = tests[1:8]

    if all(t["pass"] for t in tests):
        verdict_class = "NOTHING_MOVES_EARLY_AND_NOTHING_EVER_SETTLES"
        verdict = (
            "NOTHING MOVES EARLY, AND NOTHING EVER SETTLES. The causality "
            "objection is aimed one step away from a real defect, and the real "
            "defect is worse than the one that was suspected.\n\n"
            "THE FRONT IS CAUSAL. The amplitude at the antipode is "
            f"{t2['antipode_at_launch']:.0e} at launch and "
            f"{t2['antipode_at_two']:.0e} at t = 2.0, reaching "
            f"{t2['antipode_at_pi']:.3f} only at t ≈ π; signal ahead of the "
            f"front sits at {t2['worst_ahead_of_the_front']:.0e}, the scheme's "
            "own noise floor. Nothing outruns the front, and the early motion "
            "is not in the solve.\n\n"
            "BUT A CONSTANT NEVER LEAVES. A Gaussian of width w carries a "
            f"monopole of {t3['monopole_at_launch']:.6f} against a predicted "
            f"w²/4 = {t3['predicted_w_squared_over_four']:.6f}, and ℓ = 0 has "
            "ω = 0, so it never propagates, radiates or decays. It is NOT an "
            "early response — ahead of the front the higher modes cancel it "
            "exactly — it is a permanent one. The time-averaged displacement "
            f"is {t3['far_side_average']:.6f} at the far side, and the "
            "quietest instant of a whole run still leaves max|u| = "
            f"{t3['quietest_moment']:.3f}. The surface never returns home.\n\n"
            "AND THE MODE IS FORBIDDEN ANYWAY. Electromagnetism has no monopole "
            "radiation and gravity has none at ℓ = 0 or ℓ = 1, so the offending "
            "mode is precisely the one real radiation cannot have. That is not "
            "a blemish on the analogy, it is outside it. The checkable form is "
            "the spectrum: the spin-2 operator starts at ℓ = 2, so it admits no "
            "zero-frequency mode and can hold no DC at all. Measured as the "
            "time-averaged value over four returns, "
            f"{t4['spin2_worst_dc']:.1e} against the scalar's "
            f"{t4['scalar_worst_dc']:.1e} — a factor of {t4['ratio']:.0f}.\n\n"
            "AND THE FIX HAS TO STAY INSIDE THE PULSE. Subtracting the mean "
            "moves the problem rather than removing it, leaving "
            f"{t5['mean_subtracted_antipode']:.5f} at the antipode — exactly "
            "the offset it was meant to remove. Cancelling with a wider "
            "Gaussian fixes the mean and costs four orders of far-side quiet, "
            f"{t5['wide_gaussian_far_side']:.0e} against "
            f"{t5['launch_rows'][0]['antipode_before_arrival']:.0e}, because "
            "the corrector's tail is fatter than the pulse it corrects. Built "
            "from compactly supported bumps both conditions hold at once: "
            f"monopole {t5['compact_monopole']:.1e} and a far side of "
            f"{t5['compact_far_side']:.0e} — quieter than the original, "
            "because it is exactly zero rather than merely small.\n\n"
            "THE RING DOES CONCENTRATE, CORRECTLY. A²·(circumference) is "
            f"conserved along the converging front to {t6['spread']:.3f}, which "
            "is A ∝ 1/√(sin d). The solve is right.\n\n"
            "YOU JUST CANNOT SEE IT, AND IT CANNOT WIN. 1/√(sin d) is flat "
            f"across the middle: {100 * t7['flat_fraction_of_the_trip']:.0f}% of "
            "the trip stays within 50% of the equatorial minimum and the growth "
            f"is all after t = {t7['growth_happens_after_t']:.2f}. Worse for "
            "the expectation, on a compact surface the launch is itself a "
            f"focus — the focal height is {t7['mean_focal_ratio']:.2f}× the "
            "LAUNCH height, so the focus never beats the source. The 5.97× "
            "growth is relative to the equator. Unbounded growth belongs to an "
            "open geometry, not this one.\n\n"
            "AND HEIGHT IS NOT ENERGY. The energy density is u̇² + |∇u|², so "
            "the constant offset displaces every point of the surface and "
            f"carries exactly {t8['energy_of_the_offset']:.1f} energy against "
            f"{t8['energy_of_the_pulse']:.3f} for the pulse. The most global "
            "feature of the drawn shape is invisible to the physics. That is "
            "the measurable core of the objection to drawing a field as a "
            "displacement: the eye is on height, the physics is in the "
            "gradient.\n\n"
            "SCOPE. This audits one representation; it does not propose its "
            "replacement. What it establishes is which objections are about the "
            "solve — none — which are about the initial data, and which are "
            "about drawing a field as a displacement at all."
        )
    else:
        verdict_class = "WAVE_CONSTRAINTS_INCONCLUSIVE"
        verdict = ("INCONCLUSIVE. A check failed; review the causal ladder, "
                   "the monopole, the launch comparison, the ring bookkeeping, "
                   "or the energy of the offset.")

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "an audit of the height-field representation against causality, "
            "mode content, ring bookkeeping and energy"
        ),
        "the_front": "causal",
        "the_defect": "a permanent ℓ = 0 offset the closed geometry cannot shed",
        "the_fix": "a compactly supported, monopole-free source",
        "the_unfixable": "height is not energy",
        "geometry": {"pulse_width": PULSE_WIDTH},
        "tests": tests,
        "n_passed": tests[-1]["n_passed"],
        "n_total": tests[-1]["n_total"],
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    t = {x["name"]: x for x in s["tests"]}
    out: List[str] = []
    out.append("# What a drawn wave has to obey\n")
    out.append(f"_{s['timestamp_utc']}_\n")

    c = t["T2_the_front_is_causal"]
    out.append("## The front is causal\n")
    out.append("| t | amplitude at the antipode | ahead of the front |")
    out.append("|---:|---:|---:|")
    for r in c["rows"]:
        ahead = ("—" if r["ahead_of_the_front"] != r["ahead_of_the_front"]
                 else f"{r['ahead_of_the_front']:.1e}")
        out.append(f"| {r['t']:.2f} | `{r['antipode']:.1e}` | `{ahead}` |")
    out.append("\nNothing outruns the front.\n")

    c = t["T3_a_constant_never_leaves"]
    out.append("## But a constant never leaves\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| monopole at launch | {c['monopole_at_launch']:.6f} |")
    out.append(f"| predicted `w²/4` | {c['predicted_w_squared_over_four']:.6f} |")
    out.append(f"| drift over a full return | `{c['monopole_drift_over_a_full_return']:.1e}` |")
    out.append(f"| time-averaged displacement, far side | {c['far_side_average']:.6f} |")
    out.append(f"| quietest instant, `max\\|u\\|` | {c['quietest_moment']:.4f} |\n")

    c = t["T5_the_fix_stays_inside_the_pulse"]
    out.append("## Only a localised, compact correction works\n")
    out.append("| source | monopole | far side before arrival |")
    out.append("|---|---:|---:|")
    for r in c["launch_rows"]:
        out.append(f"| {r['source']} | `{r['monopole']:+.1e}` | "
                   f"`{r['antipode_before_arrival']:.1e}` |")
    out.append("")

    c = t["T7_but_you_cannot_see_it"]
    out.append("## The ring concentrates where you cannot see it\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| growth from the equator | {c['growth_from_the_equator']:.2f}× |")
    out.append(f"| flat fraction of the trip | {c['flat_fraction_of_the_trip']:.2f} |")
    out.append(f"| growth happens after | t = {c['growth_happens_after_t']:.2f} |")
    out.append(f"| focal height / **launch** height | {c['mean_focal_ratio']:.3f} |\n")

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
    out = here / "runs" / f"{ts}_wave_constraints_probe"
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
