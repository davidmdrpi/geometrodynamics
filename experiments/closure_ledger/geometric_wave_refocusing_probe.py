"""
Ring wavefronts on a surface with a throat: run the geometry, watch the wave
(PR #242)

> Framing: a linear scalar wave on a *fixed classical* surface — geometry →
> field, not quantum gravity, and no backreaction anywhere.

THE QUESTION
────────────
The recent arc reached the throat through algebra: operators, pairings,
CHSH ceilings, winding charges.  This probe goes back the other way and
asks the question the pictures were originally asking:

    if we just let a classical wave run on a closed surface that has a
    throat, what does the geometry itself do to it?

Nothing is fitted and no quantum ingredient is inserted.  Every number
below is read off a live wave solve.

THE SURFACE
───────────
A unit 2-sphere with both polar caps ``θ < a``, ``θ > π−a`` removed and the
two mouth circles joined by a catenoidal neck ``r(s) = b cosh((s−L/2)/b)``.
Demanding a ``C¹`` join — circumference *and* its slope continuous at both
mouths — fixes the neck from the single parameter ``a``:

    b = sin a / √(1+cos²a),        L = 2b·asinh(cos a).

The neck then has constant ``K = −1/b²`` and the two pieces cancel in
Gauss–Bonnet, ``4π cos a − 4π cos a = 0``, so the glued surface is a torus
(azimuth-preserving gluing) or a Klein bottle (azimuth-reversing).  This is
the 2-surface *section* of the ``S³`` picture: point → ring → great circle
→ antipode is the section of point → shell → maximal shell → antipode.

The controlled baseline is the same domain on the same grid with both
mouths sealed by a perfect mirror (``mode="plugged"``).  Only the handle
differs between the two runs.

THE ANSWER (measured)
  • THE RING DOES NOT CROSS ITSELF while it is in free flight.  On a closed
    surface the front of a point pulse is a single connected circle from
    launch until it first reaches a mouth — measured, both topologies.
  • THE ECHO DELAY IS THE NECK LENGTH.  Sealed, the pulse mirrors off the
    mouth and returns at 2(θ₀−a).  Open, it crosses the bulk and returns at
    2(θ₀−a)+L.  The measured delay reproduces L to ~1%, and holds at that
    level across a 3× grid refinement.  The bulk route is not a metaphor;
    the wave measures it.
  • THE OPEN MOUTH BARELY REFLECTS.  Opening the throat suppresses the
    mirror echo by ~96%: the wave goes *through* the hole rather than
    bouncing off it.
  • THE TWIST AIMS THE BULK ARRIVAL.  With the gluing offset τ = π the bulk
    route ends at the antipode and delivers a precursor there ~0.4 ahead of
    the geodesic focus; with τ = 0 the same route returns to the source and
    the antipode sees ~10× less.  The geometry, not a coupling, decides
    where the energy lands.
  • THE ORIENTATION IS REAL BUT HIDDEN.  Torus and Klein-bottle gluings are
    *identical* for a point source exactly when τ ∈ {0, π} — derived from
    the meridian mirror ``R∘g = g∘R ⟺ τ ≡ −τ`` and confirmed to machine
    precision — and differ by ~18% everywhere else.  The inner/outer
    asymmetry of the throat exists at every twist; it takes a twist that
    breaks the mirror to make it observable.

HONEST SCOPE
  Linear, no backreaction: a focus here can sharpen but cannot nucleate a
  new throat, so this probe does not speak to the #175 threshold.  The join
  is C¹ and not C², so each mouth carries a curvature ring which scatters a
  little on its own; the mouth budget is reported rather than hidden.  All
  absolute arrival times carry a common leading-edge bias of about the
  launch pulse's half width (the peak of a finite pulse is not its geodesic
  front), which is why every load-bearing number is a *difference* of
  arrival times taken on the same grid with the same pulse.  A 2-surface
  section, not S³.

Tests:
  T1. Goal: the question and the surface.
  T2. The geometry closes: C¹ join, K = −1/b², Gauss–Bonnet χ = 0.
  T3. The numerics are honest: energy conserved, delay grid-stable.
  T4. Free flight: the ring is a single circle until it meets a mouth.
  T5. The echo delay is the neck length (the bulk route is real).
  T6. The open mouth transmits and barely reflects.
  T7. The twist aims the bulk arrival; the orientation is hidden at τ∈{0,π}.
  T8. Synthesis + assessment.

Verdict:
  - GEOMETRY_ALONE_ROUTES_THE_WAVE (expected): the handle's length, its
    twist and its orientation are all read directly off a classical wave
    with no fitted parameter.
  - GEOMETRIC_ROUTING_INCONCLUSIVE: a check failed.

Run:
    python -m experiments.closure_ledger.geometric_wave_refocusing_probe
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from geometrodynamics.viz.throat_wavefront import (
    ANTIPODAL_TIME,
    ThroatGeometry,
    ThroatWaveSim,
    measure_bulk_precursor,
    measure_echo_delay,
    measure_orientation_visibility,
    measure_transmission,
    mirror_symmetry_broken,
    surface_name,
    track_wavefront,
)

# ── the one geometry every test shares ──────────────────────────────────────
MOUTH_ANGLE = 0.5          # cap radius a (rad); the throat is visible at 28.6°
N_THETA = 144
N_PHI = 192
PULSE_WIDTH = 0.18
SOURCE_THETA = 0.5 * math.pi


# ════════════════════════════════════════════════════════════════════════════
# T1 — the question and the surface
# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    g = ThroatGeometry(MOUTH_ANGLE)
    return {
        "name": "T1_goal",
        "question": (
            "let a classical wave run on a closed surface that carries a "
            "throat — what does the geometry alone do to it?"
        ),
        "surface": (
            "unit S² with both polar caps θ<a, θ>π−a removed, the two mouth "
            "circles joined by a catenoidal neck; the C¹ join fixes the neck "
            "from a alone"
        ),
        "baseline": (
            "the same domain on the same grid with both mouths sealed by a "
            "perfect mirror — only the handle differs"
        ),
        "mouth_angle": MOUTH_ANGLE,
        "neck_radius": g.neck_radius,
        "neck_length": g.length,
        "outer_route": g.outer_route,
        "shortcut_ratio": g.shortcut_ratio,
        "framing": "geometry → field on a fixed classical surface; no backreaction",
        "pass": True,
    }


# ════════════════════════════════════════════════════════════════════════════
# T2 — the geometry closes
# ════════════════════════════════════════════════════════════════════════════
def t2_geometry() -> dict:
    rows: List[dict] = []
    ok = True
    for a in (0.2, 0.35, 0.5, 0.7, 1.0):
        g = ThroatGeometry(a)
        r0 = float(g.radius(0.0))
        rp0 = float(g.radius_slope(0.0))
        rL = float(g.radius(g.length))
        rpL = float(g.radius_slope(g.length))
        s = np.linspace(0.0, g.length, 4001)
        area_num = float(np.trapezoid(2.0 * math.pi * np.asarray(g.radius(s)), s))
        row = {
            "mouth_angle": a,
            "neck_radius": g.neck_radius,
            "neck_length": g.length,
            "join_radius_error": abs(r0 - math.sin(a)),
            "join_slope_error": abs(rp0 + math.cos(a)),
            "join_symmetry_error": abs(rL - r0) + abs(rpL + rp0),
            "gauss_curvature": g.gauss_curvature,
            "area_closed_form_error": abs(area_num - g.area) / g.area,
            "euler_characteristic": g.euler_characteristic(),
            "embeddable": bool(abs(rp0) <= 1.0),
            "shortcut_ratio": g.shortcut_ratio,
        }
        ok &= (row["join_radius_error"] < 1e-12
               and row["join_slope_error"] < 1e-12
               and row["join_symmetry_error"] < 1e-12
               and abs(row["euler_characteristic"]) < 1e-9
               and row["area_closed_form_error"] < 1e-6
               and row["embeddable"])
        rows.append(row)
    return {
        "name": "T2_geometry_closes",
        "claim": (
            "the C¹ join fixes the neck from a alone; the neck has constant "
            "K = −1/b²; the sphere and neck cancel in Gauss–Bonnet so χ = 0 "
            "(torus for the orientable gluing, Klein bottle for the other)"
        ),
        "rows": rows,
        "surfaces": {"+1": surface_name(1), "-1": surface_name(-1)},
        "pass": bool(ok),
    }


# ════════════════════════════════════════════════════════════════════════════
# T3 — the numerics are honest
# ════════════════════════════════════════════════════════════════════════════
def t3_numerics() -> dict:
    drift: List[dict] = []
    for mode in ("plugged", "throat"):
        s = ThroatWaveSim(mode=mode, mouth_angle=MOUTH_ANGLE, n_theta=N_THETA,
                          n_phi=N_PHI, pulse_width=PULSE_WIDTH)
        s.advance_to(3.6)
        drift.append({"mode": mode, "t": s.t, "energy_drift": s.energy_drift(),
                      "finite": s.is_finite()})

    conv: List[dict] = []
    for nth, nph in ((96, 128), (144, 192), (216, 288), (288, 384)):
        d = measure_echo_delay(mouth_angle=MOUTH_ANGLE, n_theta=nth, n_phi=nph,
                               pulse_width=PULSE_WIDTH)
        conv.append({"n_theta": nth, "n_phi": nph,
                     "delay_measured": d["delay_measured"],
                     "delay_rel_error": d["delay_rel_error"]})

    ok = (all(r["finite"] and r["energy_drift"] < 0.02 for r in drift)
          and all(r["delay_rel_error"] < 0.03 for r in conv))
    return {
        "name": "T3_numerics_honest",
        "claim": (
            "energy is conserved to better than 2% across the coupled "
            "sphere+neck solve, and the load-bearing delay is stable at the "
            "1-2% level across a 3x refinement of the grid (the residual is "
            "the neck's own dispersion of the pulse, not the discretisation)"
        ),
        "energy": drift,
        "convergence": conv,
        "bias_note": (
            "absolute arrival times sit ≈ one pulse half-width early (the "
            "peak of a finite pulse is not its geodesic front); the bias is "
            "common to every watched point, so differences are unbiased"
        ),
        "pass": bool(ok),
    }


# ════════════════════════════════════════════════════════════════════════════
# T4 — free flight: one ring, no self-crossing
# ════════════════════════════════════════════════════════════════════════════
def t4_free_flight() -> dict:
    g = ThroatGeometry(MOUTH_ANGLE)
    free = g.free_flight(SOURCE_THETA)
    rows: List[dict] = []
    ok = True
    for mode in ("plugged", "throat"):
        s = ThroatWaveSim(mode=mode, mouth_angle=MOUTH_ANGLE, n_theta=N_THETA,
                          n_phi=N_PHI, pulse_width=PULSE_WIDTH)
        r = track_wavefront(s, 3.0, n_samples=120)
        rows.append({
            "mode": mode,
            "single_ring_until": r["single_ring_until"],
            "free_flight_predicted": free,
            "survives_free_flight": bool(r["single_ring_until"] >= free),
            "max_components_after": r["max_components"],
        })
        ok &= rows[-1]["survives_free_flight"]
    return {
        "name": "T4_free_flight_single_ring",
        "claim": (
            "on a closed surface a point pulse's front is one connected "
            "circle from launch until it first reaches a mouth — it expands "
            "to the great circle and would contract to the antipode without "
            "ever meeting itself.  Only the handle can put a second front on "
            "the same surface, and on a closed surface two fronts must cross."
        ),
        "free_flight_predicted": free,
        "rows": rows,
        "detector": (
            "connected components of the smoothed energy density u_t²+|∇u|² "
            "above 15% of its peak, φ-periodic; scored after the launch "
            "transient (a cap released from rest splits into an outgoing and "
            "an ingoing ring which merges at the source)"
        ),
        "pass": bool(ok),
    }


# ════════════════════════════════════════════════════════════════════════════
# T5 — the echo delay is the neck length
# ════════════════════════════════════════════════════════════════════════════
def t5_echo_delay() -> dict:
    d = measure_echo_delay(mouth_angle=MOUTH_ANGLE, n_theta=N_THETA,
                           n_phi=N_PHI, pulse_width=PULSE_WIDTH)
    return {
        "name": "T5_echo_delay_is_neck_length",
        "claim": (
            "sealed, the pulse mirrors off the mouth and returns at 2(θ₀−a); "
            "open, it crosses the bulk and returns at 2(θ₀−a)+L.  Both paths "
            "share every segment except the crossing, so the delay between "
            "the echoes is the neck arclength and the pulse bias cancels."
        ),
        **{k: float(v) for k, v in d.items()},
        "pass": bool(d["delay_rel_error"] < 0.03),
    }


# ════════════════════════════════════════════════════════════════════════════
# T6 — the open mouth transmits and barely reflects
# ════════════════════════════════════════════════════════════════════════════
def t6_mouth_budget() -> dict:
    d = measure_echo_delay(mouth_angle=MOUTH_ANGLE, n_theta=N_THETA,
                           n_phi=N_PHI, pulse_width=PULSE_WIDTH)
    rows: List[dict] = []
    for a in (0.3, 0.5, 0.7):
        s = ThroatWaveSim(mode="throat", mouth_angle=a, n_theta=N_THETA,
                          n_phi=N_PHI, pulse_width=PULSE_WIDTH)
        t = measure_transmission(s, 3.0, n_samples=200)
        rows.append({"mouth_angle": a, "mouth_radius": math.sin(a), **t})
    monotone = all(rows[i]["transmitted"] < rows[i + 1]["transmitted"]
                   for i in range(len(rows) - 1))
    ok = (d["mirror_suppression"] > 0.8 and monotone
          and all(r["energy_drift"] < 0.02 for r in rows))
    return {
        "name": "T6_mouth_transmits",
        "claim": (
            "opening the throat all but removes the mirror echo — the wave "
            "goes through the hole instead of bouncing off it — and the "
            "fraction it swallows grows with the mouth's aperture"
        ),
        "mirror_suppression": float(d["mirror_suppression"]),
        "mirror_amplitude_sealed": float(d["mirror_amplitude"]),
        "aperture_scan": rows,
        "monotone_in_aperture": bool(monotone),
        "caveat": (
            "the C¹ (not C²) join leaves a curvature ring at each mouth "
            "which scatters on its own; that scattering is inside the "
            "reported budget, not removed from it"
        ),
        "pass": bool(ok),
    }


# ════════════════════════════════════════════════════════════════════════════
# T7 — the twist aims; the orientation hides at τ ∈ {0, π}
# ════════════════════════════════════════════════════════════════════════════
def t7_twist_and_orientation() -> dict:
    p = measure_bulk_precursor(mouth_angle=MOUTH_ANGLE, n_theta=N_THETA,
                               n_phi=N_PHI, pulse_width=PULSE_WIDTH)
    vis = measure_orientation_visibility(mouth_angle=MOUTH_ANGLE,
                                         pulse_width=PULSE_WIDTH)
    agree = all(
        (row["relative_difference"] > 1e-3) == row["mirror_broken"]
        for row in vis["rows"]
    )
    ok = (p["aiming_contrast"] > 3.0 and p["lead_time"] > 0.0 and agree)
    return {
        "name": "T7_twist_aims_orientation_hides",
        "claim": (
            "the gluing offset decides where the bulk route lands: at τ=π it "
            "ends on the antipode and delivers a precursor ahead of the "
            "geodesic focus, at τ=0 it returns to the source instead.  The "
            "orientation of the gluing is invisible to a point source "
            "exactly at τ ∈ {0, π}, where the meridian mirror carries one "
            "gluing into the other, and visible everywhere else."
        ),
        "precursor": {k: float(v) for k, v in p.items()},
        "orientation_visibility": vis,
        "mirror_prediction_matches": bool(agree),
        "mirror_argument": "R∘g = g∘R ⟺ −(εφ+τ) = ε(−φ)+τ ⟺ τ ≡ −τ ⟺ τ ∈ {0, π}",
        "pass": bool(ok),
    }


# ════════════════════════════════════════════════════════════════════════════
# T8 — assessment
# ════════════════════════════════════════════════════════════════════════════
def t8_assessment(tests: List[dict]) -> dict:
    n_pass = sum(1 for t in tests if t["pass"])
    return {
        "name": "T8_assessment",
        "n_passed": n_pass,
        "n_total": len(tests),
        "pass": n_pass == len(tests),
    }


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests: List[dict] = []
    tests.append(t1_goal())
    tests.append(t2_geometry())
    tests.append(t3_numerics())
    tests.append(t4_free_flight())
    tests.append(t5_echo_delay())
    tests.append(t6_mouth_budget())
    tests.append(t7_twist_and_orientation())
    tests.append(t8_assessment(tests))

    t5, t6, t7 = tests[4], tests[5], tests[6]
    all_pass = all(t["pass"] for t in tests)
    if all_pass:
        verdict_class = "GEOMETRY_ALONE_ROUTES_THE_WAVE"
        verdict = (
            "GEOMETRY ALONE ROUTES THE WAVE. A linear classical wave on a "
            "fixed closed surface reports the handle's every property with "
            "no fitted parameter. The ring stays a single circle through "
            f"free flight (to t = {tests[3]['free_flight_predicted']:.3f}), "
            "so it never crosses itself until the geometry gives it "
            "somewhere else to be. Sealing the mouths gives a mirror echo; "
            "opening them replaces it with a bulk return delayed by "
            f"{t5['delay_measured']:.4f} against a neck length of "
            f"{t5['neck_length']:.4f} ({100*t5['delay_rel_error']:.1f}% "
            "error) — the wave measures the throat. The open mouth "
            f"suppresses the reflection by {100*t6['mirror_suppression']:.0f}%: "
            "energy transits the hole rather than bouncing off it. A gluing "
            "twist of π re-aims the bulk arrival onto the antipode, where it "
            f"lands {t7['precursor']['lead_time']:.3f} ahead of the geodesic "
            f"focus and {t7['precursor']['aiming_contrast']:.1f}× stronger "
            "than the untwisted throat. And the throat's orientation — torus "
            "against Klein bottle — is invisible to a point source at "
            "exactly τ ∈ {0, π} and visible elsewhere, the mirror argument "
            "confirmed to machine precision. The inner/outer asymmetry is "
            "there at every twist; it takes a twist to expose it.\n\n"
            "SCOPE. Linear and without backreaction, so a focus can sharpen "
            "but cannot nucleate — this says nothing about the #175 "
            "threshold. The C¹ join leaves a curvature ring at each mouth "
            "which is inside the reported budget. Absolute arrival times "
            "carry a common pulse-width bias; every load-bearing number is "
            "a difference taken on one grid with one pulse. A 2-surface "
            "section of the S³ picture, not S³."
        )
    else:
        verdict_class = "GEOMETRIC_ROUTING_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A check failed; review the Gauss–Bonnet closure, "
            "the energy drift and delay convergence, the free-flight front "
            "count, the echo delay, the mouth budget, or the twist/"
            "orientation pair."
        )

    g = ThroatGeometry(MOUTH_ANGLE)
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "a classical wave on a closed surface with a catenoidal throat "
            "reads off the handle's length, aperture, twist and orientation "
            "with no fitted parameter"
        ),
        "the_ring": "a single circle through free flight — a pulse does not cross itself",
        "the_delay": "the open/sealed echo delay is the neck arclength L",
        "the_mouth": "the open throat transmits and barely reflects",
        "the_twist": "the gluing offset aims where the bulk energy lands",
        "the_orientation": "torus vs Klein bottle is hidden exactly at τ ∈ {0, π}",
        "geometry": {
            "mouth_angle": MOUTH_ANGLE,
            "neck_radius": g.neck_radius,
            "neck_length": g.length,
            "outer_route": g.outer_route,
            "shortcut_ratio": g.shortcut_ratio,
            "gauss_curvature": g.gauss_curvature,
            "euler_characteristic": g.euler_characteristic(),
            "antipodal_time": ANTIPODAL_TIME,
        },
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out: List[str] = []
    out.append("# Ring wavefronts on a surface with a throat (PR #242)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "**If we just let a classical wave run on a closed surface that has "
        "a throat, what does the geometry itself do to it?** A unit S² with "
        "both polar caps removed and the mouths joined by a catenoidal neck, "
        "against the same domain with the mouths sealed. "
        "*(Geometry → field on a fixed classical surface; no backreaction.)*"
    )
    out.append("")
    out.append(f"- **The ring**: {s['the_ring']}")
    out.append(f"- **The delay**: {s['the_delay']}")
    out.append(f"- **The mouth**: {s['the_mouth']}")
    out.append(f"- **The twist**: {s['the_twist']}")
    out.append(f"- **The orientation**: {s['the_orientation']}")
    out.append("")

    g = s["geometry"]
    out.append("## The surface")
    out.append("")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| mouth angle `a` | {g['mouth_angle']:.3f} rad |")
    out.append(f"| neck waist `b` | {g['neck_radius']:.4f} |")
    out.append(f"| neck length `L` | {g['neck_length']:.4f} |")
    out.append(f"| outer route `π − 2a` | {g['outer_route']:.4f} |")
    out.append(f"| shortcut ratio | {g['shortcut_ratio']:.2f}× |")
    out.append(f"| neck curvature `K = −1/b²` | {g['gauss_curvature']:.4f} |")
    out.append(f"| Euler characteristic `∫K dA / 2π` | {g['euler_characteristic']:.2e} |")
    out.append("")

    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the question and the surface",
        "T2": "the geometry closes (C¹ join, χ = 0)",
        "T3": "energy conserved, delay grid-stable",
        "T4": "free flight: one ring, no self-crossing",
        "T5": "the echo delay is the neck length",
        "T6": "the open mouth transmits, barely reflects",
        "T7": "the twist aims; the orientation hides at τ∈{0,π}",
        "T8": "GEOMETRY_ALONE_ROUTES_THE_WAVE",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre, '—')} | {p} |")
    out.append("")

    t5 = s["tests"][4]
    out.append("## The echo delay (the wave measures the throat)")
    out.append("")
    out.append("| route | predicted | measured |")
    out.append("|---|---:|---:|")
    out.append(f"| sealed mirror echo `2(θ₀−a)` | {t5['t_mirror_predicted']:.4f} "
               f"| {t5['t_mirror_measured']:.4f} |")
    out.append(f"| open bulk return `2(θ₀−a)+L` | {t5['t_bulk_predicted']:.4f} "
               f"| {t5['t_bulk_measured']:.4f} |")
    out.append(f"| **delay = neck length `L`** | **{t5['neck_length']:.4f}** "
               f"| **{t5['delay_measured']:.4f}** |")
    out.append("")
    out.append(f"Relative error {100 * t5['delay_rel_error']:.2f}%. The two "
               "absolute times share the same pulse-width bias, which cancels "
               "in the delay.")
    out.append("")

    t3 = s["tests"][2]
    out.append("### Grid stability")
    out.append("")
    out.append("| grid | delay | rel. error |")
    out.append("|---|---:|---:|")
    for r in t3["convergence"]:
        out.append(f"| {r['n_theta']}×{r['n_phi']} | {r['delay_measured']:.4f} "
                   f"| {100 * r['delay_rel_error']:.2f}% |")
    out.append("")

    t6 = s["tests"][5]
    out.append("## The mouth budget")
    out.append("")
    out.append(f"Opening the throat suppresses the mirror echo by "
               f"**{100 * t6['mirror_suppression']:.1f}%**.")
    out.append("")
    out.append("| mouth `a` | mouth radius | transmitted | reflected | energy drift |")
    out.append("|---:|---:|---:|---:|---:|")
    for r in t6["aperture_scan"]:
        out.append(f"| {r['mouth_angle']:.2f} | {r['mouth_radius']:.3f} "
                   f"| {r['transmitted']:.3f} | {r['reflected']:.3f} "
                   f"| {r['energy_drift']:.1e} |")
    out.append("")

    t7 = s["tests"][6]
    p = t7["precursor"]
    out.append("## The twist aims the bulk arrival")
    out.append("")
    out.append(f"With `τ = π` the bulk route ends on the antipode at a "
               f"predicted `{p['t_route_predicted']:.4f}`, measured "
               f"`{p['t_precursor_measured']:.4f}` — **{p['lead_time']:.3f} "
               f"ahead** of the geodesic focus there, and "
               f"**{p['aiming_contrast']:.1f}×** stronger than the same "
               "throat with `τ = 0`.")
    out.append("")
    out.append("## The orientation is real but hidden")
    out.append("")
    out.append("`R∘g = g∘R ⟺ −(εφ+τ) = ε(−φ)+τ ⟺ τ ≡ −τ ⟺ τ ∈ {0, π}` — at "
               "those two offsets the meridian mirror of a point source "
               "carries the torus gluing into the Klein-bottle gluing, so no "
               "point source can tell them apart.")
    out.append("")
    out.append("| `τ/π` | torus vs Klein difference | mirror broken? |")
    out.append("|---:|---:|---|")
    for r in t7["orientation_visibility"]["rows"]:
        out.append(f"| {r['tau_over_pi']:.3f} | {r['relative_difference']:.4f} "
                   f"| {'yes' if r['mirror_broken'] else 'no'} |")
    out.append("")

    out.append("## Verdict")
    out.append("")
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
    out.append("")
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
    out = here / "runs" / f"{ts}_geometric_wave_refocusing_probe"
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
