"""
What a wavefront has to be like to fold: focal geometry across the bulk
(PR #243)

> Framing: ray and wavefront geometry on a *fixed classical* vacuole —
> geometry → field, not quantum gravity.  Nothing here is dynamical.

THE QUESTION
────────────
`geometric_wave_refocusing_probe` runs waves on a surface that already has a
throat.  This one asks the prior question -- what kind of wavefront can fold at
all -- and answers it with the differential geometry of wavefronts.  Every
number is a closed form, and the topology is measured INDEPENDENTLY of those
closed forms, from the front's own area element.

THE VACUOLE
───────────
Two concentric spheres, r = R_inner and r = R_outer, with a FLAT bulk between
them.  Along a straight ray the impact parameter b = r sin α is conserved, so a
ray launched inward reaches radius b and no deeper.

A POINT SOURCE DOES NOT FOLD -- IN THIS FLAT BULK
─────────────────────────────────────────────────
The front of a POINT is the metric sphere |x − P| = t, whose signed area
element is t² sin θ: positive for every t > 0, so it never folds, and the swept
region is the filled ball.

This is a statement about FLAT EUCLIDEAN SPACE, not about point sources.  On a
closed manifold the same source folds: on S² or S³ a point's front converges on
the antipode at t = πR, which is exactly what throat_wavefront.py measures.
Curvature is what gives a point source a focal set; a flat bulk is what denies
it one.

A CIRCLE FOCUSES *COHERENTLY* -- THAT IS WHAT IS SPECIAL
────────────────────────────────────────────────────────
Any curved extended source has a focal set, so "only a ring folds" would be
false.  What singles out the circle is not THAT it folds but HOW.  The signed
area element of its offset tube is t(ρ + t cos v), which vanishes where
cos v = −ρ/t:

  t < ρ   immersed, no fold anywhere;
  t = ρ   the FIRST caustic, infinitely degenerate -- the two roots coincide at
          the ring's centre and the WHOLE RING arrives there at once;
  t > ρ   still singular, at two axis points z = ±√(t² − ρ²) that run outward.

So the fold is not an isolated event: a maximally degenerate first focus,
followed by a persistent singular locus along the symmetry axis.

THE TWO CONDITIONS COINCIDE  (the core result)
──────────────────────────────────────────────
Demand the first caustic land on the inner sphere.  Its centre must sit at
radius R_inner, fixing cos θ₀ = R_inner/R_outer -- and that same ring's rays
leave the outer sphere at sin α = R_inner/R_outer, exactly the grazing ray
tangent to the inner sphere.  One condition, not two.  It forms at

    t = ρ = √(R_outer² − R_inner²).

ACCEPTANCE IS ASYMMETRIC; PROPAGATION IS NOT
────────────────────────────────────────────
Outer → inner accepts only sin α ≤ R_inner/R_outer; inner → outer accepts
everything.  That is an ANGULAR (SOLID-ANGLE) ACCEPTANCE ASYMMETRY, not
nonreciprocal propagation: every individual ray is exactly reversible, and the
probe checks it.  What differs is the measure of launch directions that
connect, because a hemisphere at R_outer and one at R_inner are different sets
of directions.  Nothing about the sphere's symmetry is broken.

HONEST SCOPE
  Ray/front geometry in a FLAT bulk -- exact, and independent of any wave solve.
  A curved bulk works in the effective radius r/√f(r), but only where that is
  MONOTONE on the shell; ShellGeometry validates it and refuses otherwise,
  because a photon sphere admits trapped orbits these closed forms do not
  describe.  Nothing here is dynamical: it says which fronts CAN fold, not what
  happens when one does, and it does not show a throat forming -- that needs
  backreaction, which this programme does not have.

Tests:
  T1. Goal: the vacuole and the two sources.
  T2. A point's front does not fold in a flat bulk -- and why that is a
      statement about the bulk, not about points.
  T3. A circle's front folds at t = ρ, detected from the area element alone.
  T4. The first caustic is infinitely degenerate, and singular after it.
  T5. The critical ring puts its caustic on the inner sphere, and grazes it.
  T6. Solid-angle acceptance asymmetry, with reciprocity verified.
  T7. The result is a statement about the radii, not this one shell.
  T8. Synthesis + assessment.

Verdict:
  - WHOLE_RING_FOCUSES_AT_ONE_POINT (expected).
  - RING_CAUSTIC_INCONCLUSIVE: a check failed.

Run:
    python -m experiments.closure_ledger.ring_caustic_defect_probe
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.radial_caustic import (
    PointSource,
    RingSource,
    ShellGeometry,
    measure_acceptance_asymmetry,
    measure_critical_ring,
    detect_fold,
    measure_front_topology,
)

# the programme's own vacuole: R_MID = 1, ΔR/2 = 0.26 either side
R_INNER = 0.74
R_OUTER = 1.26
N_MC = 400000


def _shell() -> ShellGeometry:
    return ShellGeometry(R_INNER, R_OUTER)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    sh = _shell()
    return {
        "name": "T1_goal",
        "question": "what kind of wavefront can make a throat?",
        "vacuole": {"r_inner": sh.r_inner, "r_outer": sh.r_outer,
                    "gap": sh.gap},
        "method": (
            "focal-set geometry of wavefronts, in closed form, checked "
            "against an independent numerical count — no wave solve"
        ),
        "sources": ["a point on the outer sphere", "a circle of latitude on it"],
        "pass": True,
    }


def t2_point_does_not_fold_in_a_flat_bulk() -> dict:
    sh = _shell()
    pulse = sh.point_source()
    top = measure_front_topology(pulse, t_hi=2.2)
    ok = (top["detected_folds"] is False
          and top["closed_form_fold_time"] is None
          and top["first_caustic_point"] is None
          and top["singular_points_after"] == 0)
    return {
        "name": "T2_point_does_not_fold_in_a_flat_bulk",
        "claim": (
            "the front of a point is the metric sphere |x−P| = t, whose signed "
            "area element t² sin θ is positive for every t > 0; it never folds "
            "and the swept region is the filled ball"
        ),
        "scope_warning": (
            "this is a property of the FLAT bulk, not of point sources.  On a "
            "closed manifold the same source folds: on S²/S³ a point's front "
            "converges on the antipode at t = πR, which throat_wavefront.py "
            "measures directly.  Curvature gives a point a focal set; a flat "
            "bulk denies it one."
        ),
        "detected_folds": top["detected_folds"],
        "scanned_to": top["scanned_to"],
        "crossing_time": pulse.crossing_time,
        "crossing_time_closed_form": sh.gap,
        "pass": bool(ok),
    }


def t3_ring_folds_detected_independently() -> dict:
    """Measured from the area element; the closed form is only the check."""
    sh = _shell()
    ring = sh.critical_ring()
    top = measure_front_topology(ring, t_hi=1.6)
    # and again on rings that have nothing to do with this shell
    cross = []
    for th in (0.3, 0.7, 1.1, 1.4):
        r2 = RingSource(shell=sh, polar_angle=th)
        d = detect_fold(r2, t_hi=1.6)
        cross.append({"polar_angle": th, "rho": r2.radius,
                      "detected": d["fold_time"],
                      "error": abs(d["fold_time"] - r2.radius)})
    ok = (top["detected_folds"]
          and top["detection_error"] < 1e-4
          and all(c["error"] < 1e-4 for c in cross))
    return {
        "name": "T3_ring_folds_detected_independently",
        "claim": (
            "the front of a circle is its offset tube, whose signed area "
            "element t(ρ + t cos v) changes sign first at t = ρ"
        ),
        "method": (
            "the fold time is found by scanning the front's own area element "
            "for a sign change and bisecting — no radius of curvature is "
            "consulted, so the comparison with ρ is a test and not a "
            "tautology.  The test is relative (min J < −tol·max|J|) because a "
            "parametrisation whose area element merely vanishes, as the "
            "direction sphere's poles do, has no caustic; and the orientation "
            "is referenced at small t because (X_u × X_v)·N carries the "
            "handedness of whatever (u,v) ordering a source happens to use."
        ),
        "detected_fold_time": top["detected_fold_time"],
        "closed_form_fold_time": top["closed_form_fold_time"],
        "detection_error": top["detection_error"],
        "cross_check": cross,
        "pass": bool(ok),
    }


def t4_first_caustic_is_degenerate_then_persists() -> dict:
    """Infinitely degenerate at t = ρ, and still singular afterwards."""
    sh = _shell()
    ring = sh.critical_ring()
    d = np.linalg.norm(ring.points(4096) - ring.centre, axis=1)
    spread = float(np.ptp(d))
    mult_at = ring.arrival_multiplicity(ring.centre, ring.fold_time)
    mult_off = ring.arrival_multiplicity(
        ring.centre + np.array([0.05, 0.0, 0.0]), ring.fold_time)
    after = ring.singular_points(1.4 * ring.radius)
    at = ring.singular_points(ring.radius)
    sep = float(abs(after[0, 2] - after[1, 2])) if after.shape[0] == 2 else 0.0
    expect = 2.0 * math.sqrt((1.4 * ring.radius) ** 2 - ring.radius ** 2)
    ok = (spread < 1e-12 and mult_at == -1 and mult_off == 2
          and at.shape[0] == 2 and abs(at[0, 2] - at[1, 2]) < 1e-12
          and after.shape[0] == 2 and abs(sep - expect) < 1e-12)
    return {
        "name": "T4_first_caustic_degenerate_then_persists",
        "claim": (
            "at t = ρ every point of the ring is equidistant from the centre, "
            "so the whole ring arrives simultaneously — an infinitely "
            "degenerate first caustic.  It is not an isolated event: for "
            "t > ρ the tube stays singular at two axis points z = ±√(t²−ρ²) "
            "which separate as the front grows."
        ),
        "ring_distance_spread": spread,
        "multiplicity_at_first_caustic": ("degenerate (whole ring)"
                                          if mult_at == -1 else mult_at),
        "multiplicity_just_off": mult_off,
        "singular_points_at_fold": int(at.shape[0]),
        "singular_points_coincide_at_fold": bool(
            at.shape[0] == 2 and abs(at[0, 2] - at[1, 2]) < 1e-12),
        "singular_separation_at_1p4_rho": sep,
        "separation_closed_form": expect,
        "pass": bool(ok),
    }


def t5_critical_ring() -> dict:
    sh = _shell()
    r = measure_critical_ring(sh)
    ok = (r["caustic_on_inner_error"] < 1e-12
          and r["grazing_error"] < 1e-12
          and abs(r["ray_turning_radius"] - sh.r_inner) < 1e-12
          and abs(r["fold_time"] - r["fold_time_closed_form"]) < 1e-12)
    return {
        "name": "T5_critical_ring_grazes",
        "claim": (
            "the ring whose first caustic lands on the inner sphere is the "
            "ring whose rays graze it: cos θ₀ = R_in/R_out puts the caustic "
            "at radius R_in, and the same ring launches at sin α = "
            "R_in/R_out, the tangent ray.  One condition, not two.  This and "
            "t = √(R_out² − R_in²) are the core result."
        ),
        **{k: float(v) for k, v in r.items()},
        "pass": bool(ok),
    }


def t6_acceptance_asymmetry() -> dict:
    sh = _shell()
    r = measure_acceptance_asymmetry(sh, n_samples=N_MC)
    ok = (abs(r["inward_closed_form"] - r["inward_monte_carlo"]) < 0.01
          and r["outward_closed_form"] == 1.0
          and r["solid_angle_ratio"] > 2.0
          and bool(r["rays_reversible"]))
    return {
        "name": "T6_solid_angle_acceptance_asymmetry",
        "claim": (
            "b = r sin α is conserved and r decreases going in, so only "
            "launch directions with sin α ≤ R_in/R_out reach the inner "
            "sphere, while every direction from the inner sphere reaches the "
            "outer one"
        ),
        "reading": (
            "this is an ANGULAR (solid-angle) acceptance asymmetry, NOT "
            "nonreciprocal propagation.  Every individual ray is exactly "
            "reversible — b is unchanged under reversal, and the accepted "
            "inward rays all have b ≤ R_in, so each one climbs back out along "
            "its own reverse.  What differs is the measure of launch "
            "directions that connect, because a hemisphere at R_out and a "
            "hemisphere at R_in are different sets of directions.  No "
            "symmetry of the sphere is broken; the ordering of the two radii "
            "is the whole of it."
        ),
        **{k: (float(v) if not isinstance(v, bool) else v)
           for k, v in r.items()},
        "pass": bool(ok),
    }


def t7_scales_with_the_ratio() -> dict:
    """The result is about R_in/R_out, not about this one shell."""
    rows: List[dict] = []
    ok = True
    for r_in in (0.2, 0.45, 0.74, 0.95, 1.15):
        sh = ShellGeometry(r_in, R_OUTER)
        ring = sh.critical_ring()
        rows.append({
            "r_inner": r_in,
            "critical_sin": sh.critical_sin,
            "theta0_deg": math.degrees(ring.polar_angle),
            "ring_radius": ring.radius,
            "fold_time": ring.fold_time,
            "caustic_radius_error": abs(ring.centre_radius - r_in),
            "inward_acceptance": sh.acceptance_fraction("inward"),
            "asymmetry_ratio": sh.acceptance_asymmetry,
            "pulse_crossing_time": sh.gap,
        })
        ok &= rows[-1]["caustic_radius_error"] < 1e-12
    # the defect always forms later than the pulse's crossing, and the
    # acceptance tightens as the inner sphere shrinks
    ok &= all(r["fold_time"] > r["pulse_crossing_time"] for r in rows)
    ok &= all(rows[i]["inward_acceptance"] < rows[i + 1]["inward_acceptance"]
              for i in range(len(rows) - 1))
    return {
        "name": "T7_scales_with_the_ratio",
        "claim": (
            "every statement is a function of R_in/R_out alone: the first "
            "caustic lands on the inner sphere at any ratio, always later "
            "than the pulse crosses, and the acceptance tightens as the "
            "inner sphere shrinks"
        ),
        "rows": rows,
        "pass": bool(ok),
    }


def t8_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T8_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_point_does_not_fold_in_a_flat_bulk(),
             t3_ring_folds_detected_independently(),
             t4_first_caustic_is_degenerate_then_persists(),
             t5_critical_ring(), t6_acceptance_asymmetry(),
             t7_scales_with_the_ratio()]
    tests.append(t8_assessment(tests))
    t5, t6 = tests[4], tests[5]
    sh = _shell()

    if all(t["pass"] for t in tests):
        verdict_class = "WHOLE_RING_FOCUSES_AT_ONE_POINT"
        verdict = (
            "THE WHOLE RING FOCUSES AT ONE POINT. Not every source folds the "
            "same way, and the circle's distinction is coherence rather than "
            "the mere existence of a caustic — any curved extended source has "
            "a focal set. What the circle does is arrive all at once.\n\n"
            "In this FLAT bulk a point source does not fold at all: its front "
            "is the metric sphere, whose signed area element stays positive, "
            f"so it crosses at t = ΔR = {sh.gap:.3f} and sweeps the filled "
            "ball behind it. That is a property of the flat bulk and not of "
            "point sources — on a closed manifold the same source converges "
            "on the antipode at t = πR, which the companion wave study "
            "measures directly.\n\n"
            "A circle's offset tube has area element t(ρ + t cos v), and "
            "scanning that element for a sign change — with no radius of "
            "curvature consulted — locates the fold at "
            f"{tests[2]['detected_fold_time']:.6f} against ρ = "
            f"{tests[2]['closed_form_fold_time']:.6f}, an error of "
            f"{tests[2]['detection_error']:.1e}, and reproduces ρ on four "
            "unrelated rings to the same precision. At that first caustic the "
            "ring's points are equidistant from the centre to "
            f"{tests[3]['ring_distance_spread']:.1e}: infinitely degenerate, "
            "the whole ring at once. And it does not end there — for t > ρ "
            "the tube stays singular at two axis points that separate as "
            f"√(t²−ρ²), measured {tests[3]['singular_separation_at_1p4_rho']:.6f} "
            f"against {tests[3]['separation_closed_form']:.6f} at t = 1.4ρ.\n\n"
            "THE CORE RESULT. Requiring that first caustic to land on the "
            f"inner sphere gives θ₀ = {t5['polar_angle_deg']:.2f}°, and that "
            f"same ring launches at sin α = {t5['launch_sin']:.6f} = "
            "R_in/R_out — the grazing ray, tangent to the inner sphere "
            f"(turning radius {t5['ray_turning_radius']:.4f} against R_in = "
            f"{sh.r_inner:.4f}). The ring that focuses on the throat and the "
            "ray that grazes it are the same ring, and it forms at "
            f"t = √(R_out² − R_in²) = {t5['fold_time']:.4f}.\n\n"
            "ACCEPTANCE, NOT NONRECIPROCITY. Outer to inner, "
            f"{100 * t6['inward_closed_form']:.1f}% of the hemisphere arrives "
            f"(Monte-Carlo {100 * t6['inward_monte_carlo']:.1f}%); inner to "
            f"outer, all of it. That is a {t6['solid_angle_ratio']:.1f}× "
            "difference in the *measure of launch directions that connect*, "
            "not in propagation: every ray is reversible and the probe "
            f"verifies it ({t6['rays_reversible']}). No symmetry is broken.\n\n"
            "SCOPE. Ray and front geometry in a flat bulk: exact, and "
            "independent of any wave solve. It says which fronts CAN fold, "
            "not what happens when one does — showing a throat actually form "
            "needs backreaction, which this programme does not have."
        )
    else:
        verdict_class = "RING_CAUSTIC_INCONCLUSIVE"
        verdict = ("INCONCLUSIVE. A check failed; review the focal sets, the "
                   "degeneracy at the focus, the critical ring, the "
                   "acceptance asymmetry, or the ratio scan.")

    ring = sh.critical_ring()
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "a pulse cannot make a wavefront defect and a ring must, because "
            "the focal set of a point is empty and the focal set of a circle "
            "is its centre"
        ),
        "the_pulse": "front stays immersed in a flat bulk; fills its own void",
        "the_ring": "whole ring arrives at its centre at t = ρ, then stays singular on the axis",
        "the_coincidence": "the ring that focuses on the inner sphere grazes it (core result)",
        "the_asymmetry": "solid-angle acceptance ~19% vs 100%; rays remain reversible",
        "geometry": {
            "r_inner": sh.r_inner, "r_outer": sh.r_outer, "gap": sh.gap,
            "critical_sin": sh.critical_sin,
            "critical_angle_deg": math.degrees(sh.critical_angle),
            "ring_theta0_deg": math.degrees(ring.polar_angle),
            "ring_radius": ring.radius,
            "fold_time": ring.fold_time,
            "pulse_crossing_time": sh.gap,
        },
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out: List[str] = []
    g = s["geometry"]
    out.append("# What a wavefront has to be like to fold (PR #243)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "**What kind of wavefront can fold at all?** Answered with the focal "
        "geometry of wavefronts. The closed forms are exact; the topology is "
        "measured *independently* of them, from the front's own area element. "
        "*(Flat bulk, fixed classical vacuole; nothing here is dynamical.)*")
    out.append("")
    for k, lab in (("the_pulse", "The pulse"), ("the_ring", "The ring"),
                   ("the_coincidence", "The coincidence"),
                   ("the_asymmetry", "The asymmetry")):
        out.append(f"- **{lab}**: {s[k]}")
    out.append("")
    out.append("## The vacuole")
    out.append("")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| `R_inner` | {g['r_inner']:.4f} |")
    out.append(f"| `R_outer` | {g['r_outer']:.4f} |")
    out.append(f"| `ΔR` | {g['gap']:.4f} |")
    out.append(f"| grazing `sin α = R_in/R_out` | {g['critical_sin']:.4f} |")
    out.append(f"| grazing angle | {g['critical_angle_deg']:.2f}° |")
    out.append(f"| critical ring `θ₀` | {g['ring_theta0_deg']:.2f}° |")
    out.append(f"| ring radius `ρ` | {g['ring_radius']:.4f} |")
    out.append(f"| **pulse crosses at** | **{g['pulse_crossing_time']:.4f}** |")
    out.append(f"| **ring first caustic at** | **{g['fold_time']:.4f}** |")
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the vacuole and the two sources",
        "T2": "a point does not fold — in this flat bulk",
        "T3": "the ring's fold, detected from the area element alone",
        "T4": "degenerate first caustic, singular afterwards",
        "T5": "the critical ring grazes the inner sphere",
        "T6": "solid-angle acceptance; rays stay reversible",
        "T7": "a statement about the radii, not this shell",
        "T8": "WHOLE_RING_FOCUSES_AT_ONE_POINT",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre, '—')} | {p} |")
    out.append("")

    t2, t3, t4 = s["tests"][1], s["tests"][2], s["tests"][3]
    out.append("## A point does not fold here; a circle folds coherently")
    out.append("")
    out.append("| | folds? | detected | closed form | error |")
    out.append("|---|---|---:|---:|---:|")
    out.append(f"| point source (flat bulk) | **no**, scanned to "
               f"{t2['scanned_to']:.2f} | — | — | — |")
    out.append(f"| ring source | **yes** | {t3['detected_fold_time']:.6f} "
               f"| {t3['closed_form_fold_time']:.6f} "
               f"| {t3['detection_error']:.1e} |")
    out.append("")
    out.append(f"> {t2['scope_warning']}")
    out.append("")
    out.append(f"{t3['method']}")
    out.append("")
    out.append("| ring `θ₀` | `ρ` | detected | error |")
    out.append("|---:|---:|---:|---:|")
    for c in t3["cross_check"]:
        out.append(f"| {c['polar_angle']:.2f} | {c['rho']:.6f} "
                   f"| {c['detected']:.6f} | {c['error']:.1e} |")
    out.append("")
    out.append(f"At the first caustic the ring's points are equidistant to "
               f"`{t4['ring_distance_spread']:.1e}` — the whole ring arrives "
               f"at once ({t4['multiplicity_at_first_caustic']}), against "
               f"multiplicity {t4['multiplicity_just_off']} just off it. And "
               "it does not end there: for `t > ρ` the tube stays singular at "
               "two axis points, separating as `√(t²−ρ²)` — measured "
               f"`{t4['singular_separation_at_1p4_rho']:.6f}` against "
               f"`{t4['separation_closed_form']:.6f}` at `t = 1.4ρ`.")
    out.append("")

    t5 = s["tests"][4]
    out.append("## The critical ring grazes the inner sphere")
    out.append("")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| ring `θ₀` | {t5['polar_angle_deg']:.4f}° |")
    out.append(f"| first caustic at radius | {t5['caustic_radius']:.6f} |")
    out.append(f"| `R_inner` | {t5['r_inner']:.6f} |")
    out.append(f"| error | {t5['caustic_on_inner_error']:.1e} |")
    out.append(f"| launch `sin α` | {t5['launch_sin']:.6f} |")
    out.append(f"| grazing `sin α` | {t5['critical_sin']:.6f} |")
    out.append(f"| error | {t5['grazing_error']:.1e} |")
    out.append(f"| ray turning radius | {t5['ray_turning_radius']:.6f} |")
    out.append("")

    t6 = s["tests"][5]
    out.append("## Solid-angle acceptance is asymmetric; propagation is not")
    out.append("")
    out.append("| direction | closed form | Monte-Carlo |")
    out.append("|---|---:|---:|")
    out.append(f"| outer → inner | {t6['inward_closed_form']:.4f} | "
               f"{t6['inward_monte_carlo']:.4f} |")
    out.append(f"| inner → outer | {t6['outward_closed_form']:.4f} | "
               f"{t6['outward_monte_carlo']:.4f} |")
    out.append("")
    out.append(f"A **{t6['solid_angle_ratio']:.1f}×** difference in accepted "
               f"solid angle. Rays reversible: **{t6['rays_reversible']}**.")
    out.append("")
    out.append(f"> {t6['reading']}")
    out.append("")

    t7 = s["tests"][6]
    out.append("## It scales with the ratio, not the shell")
    out.append("")
    out.append("| `R_in` | `sin α_crit` | `θ₀` | `ρ` | caustic at | pulse at | inward accept |")
    out.append("|---:|---:|---:|---:|---:|---:|---:|")
    for r in t7["rows"]:
        out.append(f"| {r['r_inner']:.2f} | {r['critical_sin']:.4f} "
                   f"| {r['theta0_deg']:.1f}° | {r['ring_radius']:.4f} "
                   f"| {r['fold_time']:.4f} | {r['pulse_crossing_time']:.4f} "
                   f"| {r['inward_acceptance']:.4f} |")
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
    out = here / "runs" / f"{ts}_ring_caustic_defect_probe"
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
