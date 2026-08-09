"""
Why a throat needs a ring: front topology across the bulk (PR #243)

> Framing: ray and wavefront geometry on a *fixed classical* vacuole —
> geometry → field, not quantum gravity.  Nothing here is dynamical.

THE QUESTION
────────────
`geometric_wave_refocusing_probe` runs waves on a surface that already has a
throat.  This one asks the prior question:

    what kind of wavefront can make one?

and answers it with the differential geometry of wavefronts — focal sets —
rather than with a simulation.  Every number below is a closed form, checked
against an independent numerical count.

THE VACUOLE
───────────
Two concentric spheres, r = R_inner and r = R_outer, bulk between them.  Along
a straight ray the impact parameter b = r sin α is conserved, so a ray launched
inward reaches radius b and no deeper.

A PULSE CANNOT MAKE A DEFECT
────────────────────────────
The front of a POINT source is the metric sphere |x − P| = t.  A sphere is
embedded at every t, so the front never touches itself and the swept region is
the filled ball — the pulse "just fills its own void".  The reason is exact:
THE FOCAL SET OF A POINT IS EMPTY, so the front has nothing to fold on.  It
crosses the bulk at t = ΔR and that is all it does.

A RING MUST
───────────
The front of a CIRCLE of radius ρ is the offset tube of that circle, and a
curve does have a focal set — its centres of curvature.  Every point of a
circle shares the SAME centre, so the focal set collapses to one point and the
whole ring arrives there simultaneously at t = ρ.  A degenerate caustic of
infinite multiplicity: the front stops being embedded.  That is the defect.

THE TWO CONDITIONS COINCIDE
───────────────────────────
Demand the defect land on the inner sphere.  Its centre must sit at radius
R_inner, fixing cos θ₀ = R_inner/R_outer — and that same ring's rays leave the
outer sphere at sin α = R_inner/R_outer, exactly the grazing ray tangent to the
inner sphere.  One condition, not two.  The defect forms at

    t = ρ = √(R_outer² − R_inner²).

THE ASYMMETRY
─────────────
Because b = r sin α is conserved and r decreases going in, outer → inner
accepts only sin α ≤ R_inner/R_outer while inner → outer accepts everything.
Same bulk, same rays, opposite directions, ~5× different acceptance.  This is
the inner/outer asymmetry in its plainest form: not a broken symmetry of the
sphere, just the ordering of the two radii.

HONEST SCOPE
  Ray/front geometry in a FLAT bulk — exact, and independent of any wave solve.
  A curved bulk replaces r by the effective radius r/√f(r) and the structure
  carries over, but the closed forms quoted are the flat ones.  Nothing here is
  dynamical: it says which fronts CAN fold, not what happens when one does, and
  in particular it does not show a throat forming — that needs backreaction,
  which this programme does not yet have.

Tests:
  T1. Goal: the vacuole and the two sources.
  T2. A point's focal set is empty; its front never self-intersects.
  T3. A circle's focal set is one point; its front folds at t = ρ.
  T4. The defect is degenerate: the whole ring arrives at once.
  T5. The critical ring puts the defect on the inner sphere, and grazes it.
  T6. Acceptance asymmetry, closed form vs Monte-Carlo.
  T7. The result is a statement about the radii, not this one shell.
  T8. Synthesis + assessment.

Verdict:
  - ONLY_A_RING_FOLDS (expected).
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


def t2_point_never_folds() -> dict:
    sh = _shell()
    pulse = sh.point_source()
    top = measure_front_topology(pulse, np.linspace(0.05, 2.2, 60))
    ok = (top["focal_set_size"] == 0
          and not top["ever_self_intersects"]
          and top["degenerate_times"] == []
          and all(r["off_axis"] in (0, 1) for r in top["rows"]))
    return {
        "name": "T2_point_never_folds",
        "claim": (
            "the front of a point is the metric sphere |x−P| = t, embedded at "
            "every t; the focal set of a point is empty so there is nothing "
            "to fold on, and the swept region is the filled ball — the pulse "
            "fills its own void"
        ),
        "focal_set_size": top["focal_set_size"],
        "ever_self_intersects": top["ever_self_intersects"],
        "max_multiplicity": max(r["off_axis"] for r in top["rows"]),
        "crossing_time": pulse.crossing_time,
        "crossing_time_closed_form": sh.gap,
        "pass": bool(ok),
    }


def t3_ring_must_fold() -> dict:
    sh = _shell()
    ring = sh.critical_ring()
    top = measure_front_topology(ring, np.linspace(0.05, 1.6, 60))
    ok = (top["focal_set_size"] == 1
          and top["ever_self_intersects"]
          and abs(ring.self_intersection_time - ring.radius) < 1e-12
          and len(top["degenerate_times"]) == 1)
    return {
        "name": "T3_ring_must_fold",
        "claim": (
            "the front of a circle is its offset tube; a curve's focal set is "
            "its centres of curvature, and every point of a circle shares one "
            "centre, so the focal set is a single point reached at t = ρ"
        ),
        "focal_set_size": top["focal_set_size"],
        "self_intersection_time": ring.self_intersection_time,
        "ring_radius": ring.radius,
        "degenerate_times": top["degenerate_times"],
        "pass": bool(ok),
    }


def t4_defect_is_degenerate() -> dict:
    """The whole ring arrives at once — infinite multiplicity, not a focus."""
    sh = _shell()
    ring = sh.critical_ring()
    pts = ring.points(4096)
    d = np.linalg.norm(pts - ring.centre, axis=1)
    spread = float(np.ptp(d))
    mult_at_focus = ring.arrival_multiplicity(ring.centre,
                                              ring.self_intersection_time)
    off = ring.centre + np.array([0.05, 0.0, 0.0])
    mult_off = ring.arrival_multiplicity(off, ring.self_intersection_time)
    ok = (spread < 1e-12 and mult_at_focus == -1 and mult_off == 2)
    return {
        "name": "T4_defect_is_degenerate",
        "claim": (
            "at the focal point every point of the ring is equidistant, so "
            "the entire ring arrives simultaneously: a degenerate caustic of "
            "infinite multiplicity, where the front ceases to be embedded — "
            "a codimension-2 defect, not a smooth focus"
        ),
        "ring_distance_spread": spread,
        "multiplicity_at_focus": ("degenerate (whole ring)"
                                  if mult_at_focus == -1 else mult_at_focus),
        "multiplicity_just_off_focus": mult_off,
        "pass": bool(ok),
    }


def t5_critical_ring() -> dict:
    sh = _shell()
    r = measure_critical_ring(sh)
    ok = (r["defect_on_inner_error"] < 1e-12
          and r["grazing_error"] < 1e-12
          and abs(r["ray_turning_radius"] - sh.r_inner) < 1e-12
          and abs(r["defect_time"] - r["defect_time_closed_form"]) < 1e-12)
    return {
        "name": "T5_critical_ring_grazes",
        "claim": (
            "the ring whose defect lands on the inner sphere is the ring "
            "whose rays graze it: cos θ₀ = R_in/R_out puts the focal point at "
            "radius R_in, and the same ring launches at sin α = R_in/R_out, "
            "the tangent ray.  One condition, not two."
        ),
        **{k: float(v) for k, v in r.items()},
        "pass": bool(ok),
    }


def t6_acceptance_asymmetry() -> dict:
    sh = _shell()
    r = measure_acceptance_asymmetry(sh, n_samples=N_MC)
    ok = (abs(r["inward_closed_form"] - r["inward_monte_carlo"]) < 0.01
          and r["outward_closed_form"] == 1.0
          and r["asymmetry_ratio"] > 2.0)
    return {
        "name": "T6_acceptance_asymmetry",
        "claim": (
            "b = r sin α is conserved and r decreases going in, so only rays "
            "with sin α ≤ R_in/R_out reach the inner sphere while every ray "
            "from the inner sphere escapes.  Same bulk, opposite directions, "
            "different acceptance."
        ),
        **{k: float(v) for k, v in r.items()},
        "reading": (
            "not a broken symmetry of the sphere — the ordering of the two "
            "radii is the asymmetry"
        ),
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
            "defect_time": ring.self_intersection_time,
            "defect_radius_error": abs(ring.centre_radius - r_in),
            "inward_acceptance": sh.acceptance_fraction("inward"),
            "asymmetry_ratio": sh.acceptance_asymmetry,
            "pulse_crossing_time": sh.gap,
        })
        ok &= rows[-1]["defect_radius_error"] < 1e-12
    # the defect always forms later than the pulse's crossing, and the
    # acceptance tightens as the inner sphere shrinks
    ok &= all(r["defect_time"] > r["pulse_crossing_time"] for r in rows)
    ok &= all(rows[i]["inward_acceptance"] < rows[i + 1]["inward_acceptance"]
              for i in range(len(rows) - 1))
    return {
        "name": "T7_scales_with_the_ratio",
        "claim": (
            "every statement is a function of R_in/R_out alone: the defect "
            "lands on the inner sphere at any ratio, always later than the "
            "pulse crosses, and the acceptance tightens as the inner sphere "
            "shrinks"
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
    tests = [t1_goal(), t2_point_never_folds(), t3_ring_must_fold(),
             t4_defect_is_degenerate(), t5_critical_ring(),
             t6_acceptance_asymmetry(), t7_scales_with_the_ratio()]
    tests.append(t8_assessment(tests))
    t5, t6 = tests[4], tests[5]
    sh = _shell()

    if all(t["pass"] for t in tests):
        verdict_class = "ONLY_A_RING_FOLDS"
        verdict = (
            "ONLY A RING FOLDS. The difference between a pulse and a ring is "
            "not a matter of degree, it is the focal set. A point has none, "
            "so its front is the metric sphere — embedded at every t, never "
            "touching itself, sweeping the filled ball behind it. It crosses "
            f"the bulk at t = ΔR = {sh.gap:.3f} and does nothing else; it "
            "fills its own void. A circle's focal set is its centres of "
            "curvature, and because every point of a circle shares the same "
            "centre that set collapses to a single point which the entire "
            f"ring reaches simultaneously at t = ρ. Measured, the ring's "
            "points are equidistant from it to "
            f"{tests[3]['ring_distance_spread']:.1e}: a degenerate caustic of "
            "infinite multiplicity, where the front stops being embedded. "
            "That is a codimension-2 defect made by geometry alone.\n\n"
            "The two conditions coincide. Requiring the defect to land on the "
            f"inner sphere gives θ₀ = {t5['polar_angle_deg']:.2f}°, and that "
            f"same ring launches at sin α = {t5['launch_sin']:.6f} = "
            "R_in/R_out — the grazing ray, tangent to the inner sphere "
            f"(turning radius {t5['ray_turning_radius']:.4f} against R_in = "
            f"{sh.r_inner:.4f}). The ring that focuses on the throat and the "
            "ray that grazes it are the same ring.\n\n"
            "And the bulk is not symmetric between its two faces. Outer to "
            f"inner, only {100 * t6['inward_closed_form']:.1f}% of the "
            f"hemisphere arrives (Monte-Carlo "
            f"{100 * t6['inward_monte_carlo']:.1f}%); inner to outer, all of "
            f"it does — a {t6['asymmetry_ratio']:.1f}× asymmetry from the "
            "ordering of two radii, with no symmetry broken anywhere.\n\n"
            "SCOPE. Ray and front geometry in a flat bulk: exact, and "
            "independent of any wave solve. It says which fronts CAN fold, "
            "not what happens when one does — showing a throat actually form "
            "needs backreaction, which this programme does not yet have."
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
        "the_pulse": "front is an embedded sphere; fills its own void",
        "the_ring": "front folds onto its centre at t = ρ — a degenerate caustic",
        "the_coincidence": "the ring that focuses on the inner sphere grazes it",
        "the_asymmetry": "outer→inner ~19% accepted, inner→outer 100%",
        "geometry": {
            "r_inner": sh.r_inner, "r_outer": sh.r_outer, "gap": sh.gap,
            "critical_sin": sh.critical_sin,
            "critical_angle_deg": math.degrees(sh.critical_angle),
            "ring_theta0_deg": math.degrees(ring.polar_angle),
            "ring_radius": ring.radius,
            "defect_time": ring.self_intersection_time,
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
    out.append("# Why a throat needs a ring: front topology across the bulk (PR #243)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "**What kind of wavefront can make a throat?** Answered with the "
        "focal-set geometry of wavefronts, in closed form, checked against an "
        "independent numerical count — no wave solve. *(Fixed classical "
        "vacuole; nothing here is dynamical.)*")
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
    out.append(f"| **ring defect at** | **{g['defect_time']:.4f}** |")
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the vacuole and the two sources",
        "T2": "a point's focal set is empty — never folds",
        "T3": "a circle's focal set is one point — folds at ρ",
        "T4": "the defect is degenerate: the whole ring at once",
        "T5": "the critical ring grazes the inner sphere",
        "T6": "acceptance asymmetry, closed form vs Monte-Carlo",
        "T7": "a statement about the radii, not this shell",
        "T8": "ONLY_A_RING_FOLDS",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre, '—')} | {p} |")
    out.append("")

    t2, t3, t4 = s["tests"][1], s["tests"][2], s["tests"][3]
    out.append("## A pulse cannot; a ring must")
    out.append("")
    out.append("| | focal set | self-intersects | when |")
    out.append("|---|---:|---|---:|")
    out.append(f"| point source | {t2['focal_set_size']} points | "
               f"{'yes' if t2['ever_self_intersects'] else '**never**'} | — |")
    out.append(f"| ring source | {t3['focal_set_size']} point | **yes** | "
               f"`t = ρ = {t3['self_intersection_time']:.4f}` |")
    out.append("")
    out.append(f"At the focus every point of the ring is equidistant to "
               f"`{t4['ring_distance_spread']:.1e}` — the whole ring arrives "
               f"at once ({t4['multiplicity_at_focus']}), against multiplicity "
               f"{t4['multiplicity_just_off_focus']} just off it. That is the "
               "defect: a codimension-2 point where the front stops being "
               "embedded, not a smooth focus.")
    out.append("")

    t5 = s["tests"][4]
    out.append("## The critical ring grazes the inner sphere")
    out.append("")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| ring `θ₀` | {t5['polar_angle_deg']:.4f}° |")
    out.append(f"| defect forms at radius | {t5['defect_radius']:.6f} |")
    out.append(f"| `R_inner` | {t5['r_inner']:.6f} |")
    out.append(f"| error | {t5['defect_on_inner_error']:.1e} |")
    out.append(f"| launch `sin α` | {t5['launch_sin']:.6f} |")
    out.append(f"| grazing `sin α` | {t5['critical_sin']:.6f} |")
    out.append(f"| error | {t5['grazing_error']:.1e} |")
    out.append(f"| ray turning radius | {t5['ray_turning_radius']:.6f} |")
    out.append("")

    t6 = s["tests"][5]
    out.append("## The bulk is not symmetric between its faces")
    out.append("")
    out.append("| direction | closed form | Monte-Carlo |")
    out.append("|---|---:|---:|")
    out.append(f"| outer → inner | {t6['inward_closed_form']:.4f} | "
               f"{t6['inward_monte_carlo']:.4f} |")
    out.append(f"| inner → outer | {t6['outward_closed_form']:.4f} | "
               f"{t6['outward_monte_carlo']:.4f} |")
    out.append("")
    out.append(f"A **{t6['asymmetry_ratio']:.1f}×** asymmetry. {t6['reading'].capitalize()}.")
    out.append("")

    t7 = s["tests"][6]
    out.append("## It scales with the ratio, not the shell")
    out.append("")
    out.append("| `R_in` | `sin α_crit` | `θ₀` | `ρ` | defect at | pulse at | inward accept |")
    out.append("|---:|---:|---:|---:|---:|---:|---:|")
    for r in t7["rows"]:
        out.append(f"| {r['r_inner']:.2f} | {r['critical_sin']:.4f} "
                   f"| {r['theta0_deg']:.1f}° | {r['ring_radius']:.4f} "
                   f"| {r['defect_time']:.4f} | {r['pulse_crossing_time']:.4f} "
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
