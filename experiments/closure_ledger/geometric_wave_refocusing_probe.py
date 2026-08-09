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

THREE SURFACES, ONE CLOCK
─────────────────────────
The comparison needs all three terms, not two:

  1. BARE S² — the uncut sphere.  The canonical picture: point → expanding
     ring → great circle → contracting ring → antipodal focus.
  2. PLUGGED — both polar caps cut out, the mouths sealed by a mirror.
     What merely cutting holes does.
  3. THROAT — the same cut sphere, mouths joined by a neck.  What opening a
     second geometric route does.

The neck is a genuine CATENOID — the minimal surface of revolution, H = 0,
r = b·cosh(z/b).  Matching the circumference and its arclength slope to the
sphere at the mouth fixes it from a alone:

    b = sin²a,      L = sin 2a,      z₀ = b·artanh(cos a).

Its curvature VARIES: exactly −1 at each mouth — the sphere's +1 with its
sign flipped, no jump in magnitude — deepening to −1/sin⁴a at the waist.

ON GAUSS–BONNET, HONESTLY
─────────────────────────
For any surface of revolution K = −r''/r and dA = 2πr ds, so ∫K dA =
−2π[r'] depends only on the boundary slopes.  The C¹ join pins those, so
χ = 0 closes for the catenoid and for any other C¹-matched profile alike.
**The closure is a check on the join, never evidence for a particular
neck** — a point the earlier version of this probe overstated.

THE ANSWER (measured)
  • THE BARE FRONT NEVER MEETS ITSELF.  On a closed surface with no
    boundary the front is the geodesic circle of radius t: it sweeps each
    point exactly once in a half period.  Measured per point, the maximum
    arrival count over the whole sphere is EXACTLY 1 — no second front
    anywhere — which is also what calibrates the detector's thresholds.
  • A SEALED MOUTH SENDS ONE BACK; AN OPEN ONE DOES NOT.  Plugged, a second
    front reaches 12.7% of the source hemisphere.  Open, a second front
    covers the surface elsewhere but *none* of the source side.  The same
    fact the echo shows, resolved in space.
  • THE ECHO DELAY IS THE NECK LENGTH.  Sealed, the pulse mirrors and
    returns at 2(θ₀−a); open, it crosses the bulk and returns at
    2(θ₀−a)+L.  The delay reproduces L to well under 1%.
  • THE MOUTH TRANSMITS, MEASURED AS A FLUX.  Of the energy that actually
    reaches the mouth, most crosses into the neck.  This is integrated
    power through two surfaces — not the peak energy that happens to be
    stored in the neck, which the earlier version wrongly called
    "transmitted".
  • THE TWIST AIMS THE BULK ARRIVAL, and the orientation of the gluing is
    invisible to a point source exactly at τ ∈ {0, π}, where the source's
    own meridian mirror carries one gluing into the other.

HONEST SCOPE
  Linear, no backreaction: a focus can sharpen but cannot nucleate, so this
  says nothing about the #175 threshold.  The join is C¹ and not C², so
  each mouth carries a curvature ring which scatters on its own; that is
  inside the reported budget rather than removed from it.  All absolute
  arrival times carry a common leading-edge bias of about the launch
  pulse's half width, which is why every load-bearing number is a
  *difference* of arrival times on one grid with one pulse.  A 2-surface
  section, not S³.

Tests:
  T1. Goal: the question and the three surfaces.
  T2. The neck is a true catenoid; χ = 0 tests the join, not the profile.
  T3. The scheme is conservative: energy to round-off, flux closes.
  T4. Arrival multiplicity: bare vs sealed vs open, per point.
  T5. The echo delay is the neck length.
  T6. The mouth budget, by integrated flux.
  T7. The twist aims; the orientation hides at τ ∈ {0, π}.
  T8. Synthesis + assessment.

Verdict:
  - GEOMETRY_ALONE_ROUTES_THE_WAVE (expected).
  - GEOMETRIC_ROUTING_INCONCLUSIVE: a check failed.

Run:
    python -m experiments.closure_ledger.geometric_wave_refocusing_probe
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.throat_wavefront import (
    ANTIPODAL_TIME,
    BareSphereSim,
    ThroatGeometry,
    ThroatWaveSim,
    measure_arrival_multiplicity,
    measure_bulk_precursor,
    measure_echo_delay,
    measure_mouth_budget,
    measure_orientation_visibility,
    surface_name,
)

# ── the one geometry every test shares ──────────────────────────────────────
MOUTH_ANGLE = 0.75         # cap radius a (rad); keeps 73% of the sphere
N_THETA = 144
N_PHI = 192
PULSE_WIDTH = 0.18
SOURCE_THETA = 0.5 * math.pi
MULTI_N_THETA = 96         # the multiplicity sweep runs the sim twice
MULTI_N_PHI = 128


# ════════════════════════════════════════════════════════════════════════════
# T1 — the question and the three surfaces
# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    g = ThroatGeometry(MOUTH_ANGLE)
    return {
        "name": "T1_goal",
        "question": (
            "let a classical wave run on a closed surface that carries a "
            "throat — what does the geometry alone do to it?"
        ),
        "surfaces": [
            "bare S² — the uncut sphere: point → ring → great circle → "
            "antipodal focus",
            "plugged — both caps cut out, mouths sealed by a mirror: what "
            "cutting holes does",
            "throat — the same cut sphere, mouths joined by a catenoid: what "
            "opening a second route does",
        ],
        "neck": "a true catenoid, H = 0, r = b cosh(z/b), b = sin²a, L = sin 2a",
        "mouth_angle": MOUTH_ANGLE,
        "neck_radius": g.neck_radius,
        "neck_length": g.length,
        "outer_route": g.outer_route,
        "shortcut_ratio": g.shortcut_ratio,
        "framing": "geometry → field on a fixed classical surface; no backreaction",
        "pass": True,
    }


# ════════════════════════════════════════════════════════════════════════════
# T2 — a true catenoid, and what χ = 0 actually tests
# ════════════════════════════════════════════════════════════════════════════
def t2_geometry() -> dict:
    rows: List[dict] = []
    ok = True
    for a in (0.25, 0.5, 0.75, 1.0, 1.3):
        g = ThroatGeometry(a)
        b = g.neck_radius
        s = np.linspace(0.0, g.length, 20001)
        r = np.asarray(g.radius(s))
        z = np.asarray(g.height(s))
        area = float(np.trapezoid(2.0 * math.pi * r, s))
        k_int = float(np.trapezoid(np.asarray(g.curvature(s)) * 2.0 * math.pi * r, s))
        row = {
            "mouth_angle": a,
            "neck_radius": b,
            "waist_closed_form_error": abs(b - math.sin(a) ** 2),
            "length": g.length,
            "length_closed_form_error": abs(g.length - math.sin(2.0 * a)),
            "catenoid_residual": float(np.max(np.abs(r - b * np.cosh(z / b)))),
            "join_radius_error": abs(float(g.radius(0.0)) - math.sin(a)),
            "join_slope_error": abs(float(g.radius_slope(0.0)) + math.cos(a)),
            "curvature_at_mouth": g.curvature_at_mouth,
            "curvature_at_waist": g.curvature_at_waist,
            "curvature_ratio": g.curvature_at_waist / g.curvature_at_mouth,
            "area_closed_form_error": abs(area - g.area) / g.area,
            "neck_total_curvature": g.neck_total_curvature(),
            "neck_total_curvature_numeric": k_int,
            "euler_characteristic": g.euler_characteristic(),
            "shortcut_ratio": g.shortcut_ratio,
        }
        ok &= (row["waist_closed_form_error"] < 1e-14
               and row["length_closed_form_error"] < 1e-14
               and row["catenoid_residual"] < 1e-12
               and row["join_radius_error"] < 1e-13
               and row["join_slope_error"] < 1e-13
               and abs(row["curvature_at_mouth"] + 1.0) < 1e-12
               and row["curvature_ratio"] > 1.0
               and row["area_closed_form_error"] < 1e-6
               and abs(k_int - g.neck_total_curvature()) < 1e-4
               and abs(row["euler_characteristic"]) < 1e-10)
        rows.append(row)
    return {
        "name": "T2_true_catenoid",
        "claim": (
            "the neck is a genuine minimal surface, r = b cosh(z/b) to "
            "machine precision, with b = sin²a and L = sin 2a; its curvature "
            "varies from exactly −1 at each mouth to −1/sin⁴a at the waist"
        ),
        "gauss_bonnet_caveat": (
            "∫K dA = −2π[r'] for any surface of revolution, so χ = 0 follows "
            "from the C¹ join alone and holds for every C¹-matched profile — "
            "it is a check on the join, not evidence for the catenoid"
        ),
        "rows": rows,
        "surfaces": {"+1": surface_name(1), "-1": surface_name(-1)},
        "pass": bool(ok),
    }


# ════════════════════════════════════════════════════════════════════════════
# T3 — the scheme is conservative
# ════════════════════════════════════════════════════════════════════════════
def t3_numerics() -> dict:
    drift: List[dict] = []
    for mode in ("plugged", "throat"):
        s = ThroatWaveSim(mode=mode, mouth_angle=MOUTH_ANGLE, n_theta=N_THETA,
                          n_phi=N_PHI, pulse_width=PULSE_WIDTH)
        s.advance_to(4.0)
        drift.append({"mode": mode, "t": s.t, "energy_drift": s.energy_drift(),
                      "finite": s.is_finite()})

    s = ThroatWaveSim(mode="throat", mouth_angle=MOUTH_ANGLE, n_theta=96,
                      n_phi=128, pulse_width=PULSE_WIDTH)
    acc = 0.0
    while s.t < 2.4:
        p_n, p_s = s.mouth_power()
        acc += (p_n + p_s) * s.dt
        s.step()
    stored = s.neck_energy()
    flux_err = abs(acc - stored) / max(stored, 1e-30)

    conv: List[dict] = []
    for nth, nph in ((96, 128), (144, 192), (216, 288)):
        d = measure_echo_delay(mouth_angle=MOUTH_ANGLE, n_theta=nth, n_phi=nph,
                               pulse_width=PULSE_WIDTH)
        conv.append({"n_theta": nth, "n_phi": nph,
                     "delay_measured": d["delay_measured"],
                     "delay_rel_error": d["delay_rel_error"]})

    ok = (all(r["finite"] and r["energy_drift"] < 1e-10 for r in drift)
          and flux_err < 0.05
          and all(r["delay_rel_error"] < 0.03 for r in conv))
    return {
        "name": "T3_conservative_scheme",
        "claim": (
            "each mouth is one finite-volume face shared by a sphere cell and "
            "a neck cell, evaluated once and handed to both with opposite "
            "signs, so the discrete divergence theorem holds across the mouth "
            "and the discrete energy is conserved to round-off"
        ),
        "energy": drift,
        "mouth_flux_closes_to": flux_err,
        "convergence": conv,
        "notes": (
            "the conserved quantity is the exact leapfrog invariant — the "
            "gradient term is the cross product ⟨∇uⁿ, ∇uⁿ⁻¹⟩, since the "
            "velocity lives at the half step; the launch is a purely outgoing "
            "ring, because a cap released from rest splits into an outgoing "
            "and an ingoing front"
        ),
        "pass": bool(ok),
    }


# ════════════════════════════════════════════════════════════════════════════
# T4 — arrival multiplicity: does a front ever meet another front?
# ════════════════════════════════════════════════════════════════════════════
def t4_arrival_multiplicity() -> dict:
    g = ThroatGeometry(MOUTH_ANGLE)
    kw = dict(n_theta=MULTI_N_THETA, n_phi=MULTI_N_PHI, pulse_width=PULSE_WIDTH)
    sims = [
        ("bare", BareSphereSim(**kw)),
        ("plugged", ThroatWaveSim(mode="plugged", mouth_angle=MOUTH_ANGLE, **kw)),
        ("throat", ThroatWaveSim(mode="throat", mouth_angle=MOUTH_ANGLE, **kw)),
    ]
    rows: List[dict] = []
    for label, sim in sims:
        m = measure_arrival_multiplicity(sim, ANTIPODAL_TIME)
        rows.append({
            "surface": label,
            "max_arrivals": m.max_arrivals,
            "area_fraction_multi": m.area_fraction_multi,
            "source_side_fraction": m.source_side_fraction,
            "first_multi_time": m.first_multi_time,
        })
    # the detector has two free thresholds, so show the answer does not
    # depend on them: calibrate on the bare sphere, whose answer is known
    stability: List[dict] = []
    for hi, lo in ((0.50, 0.15), (0.60, 0.20), (0.70, 0.25), (0.80, 0.30)):
        row = {"hi": hi, "lo": lo}
        for label, mk in (("bare", lambda: BareSphereSim(**kw)),
                          ("plugged", lambda: ThroatWaveSim(
                              mode="plugged", mouth_angle=MOUTH_ANGLE, **kw)),
                          ("throat", lambda: ThroatWaveSim(
                              mode="throat", mouth_angle=MOUTH_ANGLE, **kw))):
            m = measure_arrival_multiplicity(mk(), ANTIPODAL_TIME, hi=hi, lo=lo)
            row[label] = m.source_side_fraction
        stability.append(row)

    bare, plugged, throat = rows
    ok = (bare["area_fraction_multi"] < 1e-9
          and bare["source_side_fraction"] < 1e-9
          and plugged["source_side_fraction"] > 0.01
          and throat["source_side_fraction"] < 1e-9
          and all(r["bare"] < 1e-9 and r["plugged"] > r["throat"]
                  for r in stability))
    return {
        "name": "T4_arrival_multiplicity",
        "claim": (
            "on a closed surface with no boundary the front is the geodesic "
            "circle of radius t and sweeps each point exactly once, so it "
            "never meets itself — the bare sphere returns a maximum arrival "
            "count of exactly 1.  Sealing a mouth puts a second front back "
            "toward the source; opening it puts one downstream of the neck "
            "instead and none back home."
        ),
        "detector": (
            "a per-cell hysteresis trigger on the energy density u_t²+|∇u|², "
            "armed above 50% and re-armed below 15% of that cell's own peak, calibrated on the bare sphere whose answer is known; "
            "plain local-maximum counting fails because a 2+1-dimensional "
            "wave violates Huygens and every front drags a rippling wake"
        ),
        "supersedes": (
            "counting connected components of a level set, which could not "
            "distinguish a hole cutting one ring into arcs from a genuine "
            "second front"
        ),
        "calibration": (
            "the two thresholds are set on the case whose answer is known — "
            "the bare sphere, where the front provably passes once.  Any "
            "hi ≥ 0.5 returns exactly zero second arrivals there and the "
            "sealed/open contrast survives across that range, so the "
            "conclusion is not an artefact of the pair chosen"
        ),
        "threshold_stability": stability,
        "window": ANTIPODAL_TIME,
        "free_flight": g.free_flight(SOURCE_THETA),
        "rows": rows,
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
# T6 — the mouth budget, by integrated flux
# ════════════════════════════════════════════════════════════════════════════
def t6_mouth_budget() -> dict:
    d = measure_echo_delay(mouth_angle=MOUTH_ANGLE, n_theta=N_THETA,
                           n_phi=N_PHI, pulse_width=PULSE_WIDTH)
    rows: List[dict] = []
    for a in (0.4, 0.55, 0.75, 0.9):
        s = ThroatWaveSim(mode="throat", mouth_angle=a, n_theta=N_THETA,
                          n_phi=N_PHI, pulse_width=PULSE_WIDTH)
        rows.append({"mouth_angle": a, "mouth_radius": math.sin(a),
                     **measure_mouth_budget(s)})
    monotone = all(rows[i]["transmission"] < rows[i + 1]["transmission"]
                   for i in range(len(rows) - 1))
    ok = (monotone
          and all(0.0 < r["transmission"] <= 1.0 for r in rows)
          and all(r["energy_drift"] < 1e-10 for r in rows)
          and d["mirror_suppression"] > 0.4)
    return {
        "name": "T6_mouth_budget_by_flux",
        "claim": (
            "of the energy that actually reaches the mouth, most crosses into "
            "the neck, and the fraction rises with the aperture"
        ),
        "method": (
            "integrated power through two surfaces: 'offered' is the energy "
            "crossing a reference circle a few cells inside the mouth, "
            "'through' is the energy crossing the mouth face itself, and "
            "transmission is their ratio.  On a closed surface only part of "
            "the wave ever reaches the mouth, so the total energy is the "
            "wrong denominator."
        ),
        "supersedes": (
            "the peak *stored* energy fraction inside the neck, which is a "
            "snapshot, depends on the neck's length, and ignores everything "
            "that has already passed through — it is not a transmission "
            "coefficient and should not have been reported as one"
        ),
        "aperture_scan": rows,
        "monotone_in_aperture": bool(monotone),
        "mirror_suppression": float(d["mirror_suppression"]),
        "distinction": (
            "mirror suppression is an amplitude ratio at one watched point "
            "and one watched time; transmission is an energy ratio at the "
            "mouth.  They are different measurements of the same fact and "
            "must not be quoted interchangeably."
        ),
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
    tests.append(t4_arrival_multiplicity())
    tests.append(t5_echo_delay())
    tests.append(t6_mouth_budget())
    tests.append(t7_twist_and_orientation())
    tests.append(t8_assessment(tests))

    t3, t4, t5, t6, t7 = tests[2], tests[3], tests[4], tests[5], tests[6]
    bare, plugged, throat = t4["rows"]
    all_pass = all(t["pass"] for t in tests)
    if all_pass:
        verdict_class = "GEOMETRY_ALONE_ROUTES_THE_WAVE"
        verdict = (
            "GEOMETRY ALONE ROUTES THE WAVE. A linear classical wave on a "
            "fixed closed surface reports the handle's every property with "
            "no fitted parameter, and the three-surface comparison separates "
            "what each change of geometry is responsible for. On the bare "
            "sphere the front sweeps each point exactly once — "
            f"{100*bare['area_fraction_multi']:.1f}% of the surface ever sees "
            "a second front and none of the source side does — so a pulse on "
            "a closed surface with no boundary cannot meet itself. Sealing "
            "the mouths puts a second front back toward the source over "
            f"{100*plugged['source_side_fraction']:.1f}% of that hemisphere; "
            "opening the throat puts one downstream of the neck instead and "
            f"{100*throat['source_side_fraction']:.1f}% back home. The echoes "
            "say the same thing in time: the sealed mirror echo and the open "
            f"bulk return differ by {t5['delay_measured']:.4f} against a neck "
            f"length of {t5['neck_length']:.4f} "
            f"({100*t5['delay_rel_error']:.2f}% error) — the wave measures "
            "the throat. Of the energy that actually reaches the mouth, "
            f"{100*t6['aperture_scan'][2]['transmission']:.0f}% crosses into "
            "the neck, rising with the aperture. A gluing twist of π re-aims "
            "the bulk arrival onto the antipode, "
            f"{t7['precursor']['lead_time']:.3f} ahead of the geodesic focus "
            f"and {t7['precursor']['aiming_contrast']:.1f}× stronger than the "
            "untwisted throat. And the throat's orientation — torus against "
            "Klein bottle — is invisible to a point source at exactly "
            "τ ∈ {0, π} and visible elsewhere, the mirror argument confirmed "
            "to machine precision.\n\n"
            "SCHEME. Each mouth is one shared finite-volume face, so the "
            "discrete energy is conserved to round-off "
            f"(drift {t3['energy'][1]['energy_drift']:.1e}) and the mouth "
            "flux closes against the neck's stored energy to "
            f"{100*t3['mouth_flux_closes_to']:.1f}%.\n\n"
            "SCOPE. Linear and without backreaction, so a focus can sharpen "
            "but cannot nucleate — this says nothing about the #175 "
            "threshold. The C¹ join leaves a curvature ring at each mouth, "
            "inside the reported budget. χ = 0 checks the join and not the "
            "profile. Absolute arrival times carry a common pulse-width "
            "bias; every load-bearing number is a difference on one grid "
            "with one pulse. A 2-surface section of the S³ picture, not S³."
        )
    else:
        verdict_class = "GEOMETRIC_ROUTING_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A check failed; review the catenoid identities, "
            "the energy/flux closure, the arrival multiplicity, the echo "
            "delay, the mouth budget, or the twist/orientation pair."
        )

    g = ThroatGeometry(MOUTH_ANGLE)
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "a classical wave on a closed surface with a catenoidal throat "
            "reads off the handle's length, aperture, twist and orientation "
            "with no fitted parameter"
        ),
        "the_bare_front": "sweeps each point once — a pulse cannot meet itself",
        "the_sealed_mouth": "sends a second front back toward the source",
        "the_open_mouth": "sends none back; the second front is downstream of the neck",
        "the_delay": "the open/sealed echo delay is the neck arclength L",
        "the_twist": "the gluing offset aims where the bulk energy lands",
        "the_orientation": "torus vs Klein bottle is hidden exactly at τ ∈ {0, π}",
        "geometry": {
            "mouth_angle": MOUTH_ANGLE,
            "neck_radius": g.neck_radius,
            "neck_length": g.length,
            "outer_route": g.outer_route,
            "shortcut_ratio": g.shortcut_ratio,
            "curvature_at_mouth": g.curvature_at_mouth,
            "curvature_at_waist": g.curvature_at_waist,
            "euler_characteristic": g.euler_characteristic(),
            "free_flight": g.free_flight(SOURCE_THETA),
            "mirror_echo": g.mirror_echo(SOURCE_THETA),
            "throat_loop": g.throat_loop(SOURCE_THETA),
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
        "a throat, what does the geometry itself do to it?** Three surfaces "
        "on one clock: the bare S², the same sphere with both caps cut out "
        "and sealed, and the same cut sphere with the mouths joined by a "
        "catenoid. *(Geometry → field on a fixed classical surface; no "
        "backreaction.)*"
    )
    out.append("")
    for k, label in (("the_bare_front", "The bare front"),
                     ("the_sealed_mouth", "A sealed mouth"),
                     ("the_open_mouth", "An open mouth"),
                     ("the_delay", "The delay"),
                     ("the_twist", "The twist"),
                     ("the_orientation", "The orientation")):
        out.append(f"- **{label}**: {s[k]}")
    out.append("")

    g = s["geometry"]
    out.append("## The surface")
    out.append("")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| mouth angle `a` | {g['mouth_angle']:.3f} rad |")
    out.append(f"| neck waist `b = sin²a` | {g['neck_radius']:.4f} |")
    out.append(f"| neck length `L = sin 2a` | {g['neck_length']:.4f} |")
    out.append(f"| outer route `π − 2a` | {g['outer_route']:.4f} |")
    out.append(f"| shortcut ratio | {g['shortcut_ratio']:.2f}× |")
    out.append(f"| `K` at the mouth | {g['curvature_at_mouth']:.4f} |")
    out.append(f"| `K` at the waist | {g['curvature_at_waist']:.4f} |")
    out.append(f"| `∫K dA / 2π` | {g['euler_characteristic']:.2e} |")
    out.append("")
    out.append("Four geometric times set the clock: the ring reaches the "
               f"mouths at `{g['free_flight']:.3f}`, a sealed echo returns at "
               f"`{g['mirror_echo']:.3f}`, a bulk crossing lands at "
               f"`{g['throat_loop']:.3f}`, and the antipodal focus is at "
               f"`{g['antipodal_time']:.3f}`.")
    out.append("")

    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the question and the three surfaces",
        "T2": "a true catenoid; χ = 0 tests the join",
        "T3": "conservative to round-off; flux closes",
        "T4": "arrival multiplicity per point",
        "T5": "the echo delay is the neck length",
        "T6": "the mouth budget, by integrated flux",
        "T7": "the twist aims; the orientation hides",
        "T8": "GEOMETRY_ALONE_ROUTES_THE_WAVE",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre, '—')} | {p} |")
    out.append("")

    t2 = s["tests"][1]
    out.append("## The neck is a true catenoid")
    out.append("")
    out.append("| `a` | `b = sin²a` | `L = sin 2a` | `\\|r − b cosh(z/b)\\|` | `K` mouth | `K` waist | `χ` |")
    out.append("|---:|---:|---:|---:|---:|---:|---:|")
    for r in t2["rows"]:
        out.append(f"| {r['mouth_angle']:.2f} | {r['neck_radius']:.4f} "
                   f"| {r['length']:.4f} | {r['catenoid_residual']:.1e} "
                   f"| {r['curvature_at_mouth']:.3f} "
                   f"| {r['curvature_at_waist']:.2f} "
                   f"| {r['euler_characteristic']:.1e} |")
    out.append("")
    out.append(f"> {t2['gauss_bonnet_caveat']}.")
    out.append("")

    t4 = s["tests"][3]
    out.append("## Does a front ever meet another front?")
    out.append("")
    out.append(f"Arrivals counted per grid point over `t < π`. {t4['detector']}.")
    out.append("")
    out.append("| surface | max arrivals | area with ≥2 | of the source side |")
    out.append("|---|---:|---:|---:|")
    for r in t4["rows"]:
        out.append(f"| {r['surface']} | {r['max_arrivals']} "
                   f"| {r['area_fraction_multi']:.3f} "
                   f"| {r['source_side_fraction']:.3f} |")
    out.append("")
    out.append("The bare sphere is the case with no second front at all. "
               "Sealing the mouths sends one back toward the source; opening "
               "them sends **none** back — the same fact the echo shows, "
               "resolved in space rather than in time.")
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
    out.append(f"Relative error {100 * t5['delay_rel_error']:.2f}%.")
    out.append("")

    t3 = s["tests"][2]
    out.append("### Scheme and grid stability")
    out.append("")
    out.append(f"Energy drift {t3['energy'][0]['energy_drift']:.1e} (plugged) "
               f"and {t3['energy'][1]['energy_drift']:.1e} (throat); the mouth "
               f"flux closes against the neck's stored energy to "
               f"{100*t3['mouth_flux_closes_to']:.1f}%.")
    out.append("")
    out.append("| grid | delay | rel. error |")
    out.append("|---|---:|---:|")
    for r in t3["convergence"]:
        out.append(f"| {r['n_theta']}×{r['n_phi']} | {r['delay_measured']:.4f} "
                   f"| {100 * r['delay_rel_error']:.2f}% |")
    out.append("")

    t6 = s["tests"][5]
    out.append("## The mouth budget, by integrated flux")
    out.append("")
    out.append(f"{t6['method']}")
    out.append("")
    out.append("| mouth `a` | offered | through | transmission | reflection | peak stored |")
    out.append("|---:|---:|---:|---:|---:|---:|")
    for r in t6["aperture_scan"]:
        out.append(f"| {r['mouth_angle']:.2f} | {r['offered']:.4f} "
                   f"| {r['through']:.4f} | {r['transmission']:.3f} "
                   f"| {r['reflection']:.3f} | {r['peak_stored_fraction']:.3f} |")
    out.append("")
    out.append(f"> **{t6['distinction']}** The sealed echo's amplitude is "
               f"suppressed by {100*t6['mirror_suppression']:.1f}% when the "
               "throat is opened; that is a different measurement from the "
               "transmission column and the two are not interchangeable.")
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
