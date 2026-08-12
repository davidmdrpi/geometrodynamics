"""
Focus the congruence to the pinch threshold, and let the equations decide

> Framing: still a spin-2 field on a FIXED round S².  What is new is the object
> drawn — a congruence with a deforming cross-section, not a surface and not a
> bright point.

WHY NOT A BRIGHT POINT
──────────────────────
A singularity in GR is a failure of evolution and geodesic completeness;
infinite curvature is one possible signature, not the definition.  Drawing "a
singularity" as a glowing dot assumes the answer.  The object that does not is
a congruence: integrate the geodesic-deviation equation F̈ = ½ḧF from F(0) = I
and watch the Jacobian J = det F, which is the cross-sectional area of the
bundle.  J → 0 is a caustic of the MAP.  It says nothing about the metric.

THREE THINGS THAT HAVE BEEN GETTING CONFLATED
─────────────────────────────────────────────
1  ordinary focus       J dips and recovers; the map stays invertible
2  caustic              J reaches zero; trajectories cross; geometry regular
3  curvature singularity  the geometry itself fails

WHAT IS CHECKED
───────────────
T2  RAYCHAUDHURI IS THE SAME STATEMENT, AND THE RICCI TERM IS ZERO.  The
    residual 6.7e-15 is an ALGEBRAIC IDENTITY holding to round-off, not an
    accuracy check — substituting θ and σ reproduces the deviation equation
    symbolically.  The content is that R_ab u^a u^b is identically 0 because h
    is trace-free, so none of the focusing is matter; all of it is
    shear-squared, second order in the amplitude, which is the shape of every
    threshold below.

T3  ḧ HAS TO COME FROM THE WAVE OPERATOR.  Seeding a three-sample time
    difference with h(−dt) = h(0) injects a spurious impulse ½ḣ(0), and the
    1/dt² and the dt of the update cancel, so THE ERROR DOES NOT SHRINK WITH
    dt.  Both forms converge on the same refinement ladder and converge to
    different answers: J_min = 0.628 (no caustic at all) against −12.84.
    Refinement cannot reveal this; only the operator form can.

T4  THE FRONT IS CAUSAL ONLY IF THE LAUNCH IS COMPACT.  A Gaussian reaches
    d = 2.95 at t = 2.00 against a causal bound of 2.77, and that number is
    GRID-CONVERGED — an analytic tail, not a numerical precursor.  The compact
    bump (1−x²)⁴ arrives at 2.7697 against 2.7699, and its residual earliness
    shrinks under refinement.  This is the constraint wave_constraints found
    for the scalar, and the spin-2 launch needed it too.

T5  THE THREE CASES SEPARATE.  A weak wave only dips (min J = 0.9995 on the
    converging ring).  Raising the amplitude closes the source ring first and
    the converging ring about ten times later.  Case 3 never appears, and
    cannot — see T9.

T6  TWO RINGS, TWO THRESHOLDS, A FACTOR OF TEN APART.  Source ring closes at
    peak strain 0.026, converging ring at 0.247, in a 1.2π window.  The source
    is shaken hardest and from t = 0; the converging ring has to arrive first.

T7  THE NECK IS A RING, AND SPIN WEIGHT IS WHY.  h = sin²d·q vanishes at both
    poles, so the tidal field is identically zero AT the antipode and the
    congruence there is never driven.  The neck sits on a ring of radius
    0.44 w (mean 0.443) — the same ratio across a 3.3× range of pulse
    width, to within the one-cell resolution of the grid.  A spin-2 focus has no centre.

T8  THE CAUSTIC IS A PASSAGE.  Of passage, singular termination, and
    finite-radius reconnection, the equations give passage.  At the source ring
    J crosses zero with slope −17.877, converged to −17.836 under a halved
    timestep (0.2%), and plunges to −471 and stays.  A tangency would have
    the slope going to zero; this one converges to a definite nonzero value.  The
    converging ring only GRAZES: it crosses and returns within a few
    thousandths and the depth of that excursion does not converge — so
    "does focusing drive the neck radius to zero" is answered barely, and only
    at a strain nothing physical would reach.

T9  CASE 3 IS OUT OF SCOPE, AND THAT IS ABOUT THE PROGRAM.  The background is a
    fixed round S² with curvature 1 everywhere at every time.  No Einstein
    equation, no backreaction, nothing that could diverge or terminate.  A
    caustic is the strongest thing available here.  Reconnection is likewise
    not merely unobserved but STRUCTURALLY UNAVAILABLE: each point's F is
    driven only by the external h and never by its neighbours, so the
    congruence cannot act back on anything.

T10 EVERY THRESHOLD CARRIES ITS WINDOW.  The wave on S² is exactly periodic, so
    a material point is driven over and over by the same returning ring and the
    deviation equation is a Hill equation.  The converging-ring threshold falls
    4.07× from a 1.2π window to a 4π one.  A threshold quoted without a window
    is not a number about the wave.

SCOPE
─────
gain is a strength dial and is reported as a peak strain everywhere it matters.
The deviation equation is exact in ξ and linear in h; at the strains where the
converging ring closes (0.25) the field is no longer a weak perturbation, and
that is stated rather than hidden.  What is derived is the separation of the
three cases, the two thresholds, the ring law, and the passage.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.congruence import (
    measure_case_three_is_unreachable,
    measure_raychaudhuri_is_exact,
    measure_the_caustic_is_a_passage,
    measure_the_caustic_thresholds,
    measure_the_front_is_causal,
    measure_the_neck_is_a_ring,
    measure_the_operator_form_matters,
    measure_the_three_cases,
    measure_the_threshold_depends_on_the_window,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "claim": ("visualise geodesic/tidal focusing all the way to the pinch "
                  "threshold, and let the equations say what the threshold "
                  "produces"),
        "object": "a congruence with a deforming cross-section, not a point",
        "quantity": "J = det F, the cross-sectional area of the bundle",
        "the_three_outcomes": ["caustic -> passage",
                               "caustic -> singular termination",
                               "caustic -> finite-radius reconnection"],
        "pass": True,
    }


def t2_raychaudhuri_is_exact() -> dict:
    r = measure_raychaudhuri_is_exact()
    return {"name": "T2_raychaudhuri_is_exact", **r,
            "pass": bool(r["raychaudhuri_is_exact"]
                         and r["the_ricci_term_vanishes"])}


def t3_the_operator_form_matters() -> dict:
    """The error a refinement ladder cannot see."""
    r = measure_the_operator_form_matters()
    return {"name": "T3_the_operator_form_matters", **r,
            "pass": bool(r["both_forms_converge"] and r["but_they_disagree"])}


def t4_the_front_is_causal() -> dict:
    r = measure_the_front_is_causal()
    return {"name": "T4_the_front_is_causal", **r,
            "pass": bool(r["the_gaussian_arrival_is_grid_converged"]
                         and r["the_gaussian_beats_the_bound"]
                         and r["the_compact_launch_respects_it"])}


def t5_the_three_cases_separate() -> dict:
    r = measure_the_three_cases()
    return {"name": "T5_the_three_cases_separate", **r,
            "pass": bool(r["both_cases_appear"] and r["a_weak_wave_only_dips"]
                         and r["case_three_never_appears"])}


def t6_two_rings_two_thresholds() -> dict:
    r = measure_the_caustic_thresholds()
    return {"name": "T6_two_rings_two_thresholds", **r,
            "pass": bool(r["the_source_ring_is_far_easier_to_close"]
                         and r["closing_the_ring_needs_an_enormous_strain"])}


def t7_the_neck_is_a_ring() -> dict:
    r = measure_the_neck_is_a_ring()
    return {"name": "T7_the_neck_is_a_ring", **r,
            "pass": bool(r["the_radius_scales_with_the_width"]
                         and r["the_ratio_converges"]
                         and r["the_neck_is_resolved_off_the_pole"])}


def t8_the_caustic_is_a_passage() -> dict:
    r = measure_the_caustic_is_a_passage()
    return {"name": "T8_the_caustic_is_a_passage", **r,
            "pass": bool(r["the_caustic_is_a_passage"]
                         and r["crossings_are_transversal"]
                         and r["the_source_excursion_is_resolved"])}


def t9_case_three_is_out_of_scope() -> dict:
    r = measure_case_three_is_unreachable()
    return {"name": "T9_case_three_is_out_of_scope", **r,
            "pass": bool(r["case_three_is_out_of_scope"]
                         and r["the_geometry_never_moves"]
                         and not r["evolution_terminated"])}


def t10_every_threshold_carries_its_window() -> dict:
    r = measure_the_threshold_depends_on_the_window()
    return {"name": "T10_every_threshold_carries_its_window", **r,
            "pass": bool(r["threshold_falls_with_the_window"])}


def t11_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T11_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_raychaudhuri_is_exact(),
             t3_the_operator_form_matters(), t4_the_front_is_causal(),
             t5_the_three_cases_separate(), t6_two_rings_two_thresholds(),
             t7_the_neck_is_a_ring(), t8_the_caustic_is_a_passage(),
             t9_case_three_is_out_of_scope(),
             t10_every_threshold_carries_its_window()]
    tests.append(t11_assessment(tests))
    t6, t7, t8 = tests[5], tests[6], tests[7]
    src = t8.get("source_ring", {})
    fine = t8.get("source_ring_at_half_the_timestep", {})
    slope = (src.get("crossing_slopes") or [0.0])[0]
    slope_fine = (fine.get("crossing_slopes") or [0.0])[0]
    depth = src.get("depth_past_zero", 0.0)

    if all(t["pass"] for t in tests):
        verdict_class = "THE_CAUSTIC_IS_A_PASSAGE"
        verdict = (
            "THE CAUSTIC IS A PASSAGE. Focusing was carried to the pinch "
            "threshold without assuming what the threshold would produce, and "
            "the equations chose the first of the three offered outcomes. J "
            f"crosses zero transversally — slope {slope:.3f}, converged to "
            f"{slope_fine:.3f} under a halved timestep, where a tangency would "
            f"have driven it to zero — plunges to {depth:.0f}, and the "
            "evolution continues with the solver's invariant unmoved at "
            f"{t8.get('solver_invariant_drift', 0.0):.1e}. Neither of the other two "
            "outcomes was available, and for different reasons worth keeping "
            "apart. Singular termination needs the geometry to fail, and the "
            "background here is a fixed round S² with curvature 1 at every "
            "time. Finite-radius reconnection needs the congruence to act back "
            "on something, and each point's F is driven only by the external h "
            "and never by its neighbours. So this is not 'we looked and did "
            "not find them'; it is 'this program could not have produced "
            "them', which is a different and more useful statement. What the "
            "program CAN say is where the neck forms and what it costs: on a "
            "ring around the antipode, never at it, because spin weight "
            "forces h to vanish at the pole, at a radius of "
            f"{t7.get('mean_radius_over_width', 0):.2f} w — the same ratio "
            "across a 3.3x range of pulse width, to within one grid cell; and "
            "at peak strain "
            f"{t6.get('converging_ring_threshold_strain', 0):.3f} for the "
            f"converging ring against "
            f"{t6.get('source_ring_threshold_strain', 0):.3f} for the source "
            "ring, a factor of ten apart. Even then the antipodal crossing "
            "only grazes zero and the depth of its excursion does not "
            "converge. A spin-2 focus has no centre, and closing its neck "
            "costs a strain nothing physical would reach.")
    else:
        verdict_class = "INCOMPLETE"
        verdict = ("At least one check did not hold; see the failing test "
                   "before quoting any number from this run.")

    return {
        "probe": "congruence",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "tests": tests,
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    lines = [f"# {s['probe']} probe", "", f"_generated {s['generated_utc']}_",
             ""]
    for t in s["tests"]:
        mark = "PASS" if t["pass"] else "FAIL"
        lines.append(f"## {t['name']} — {mark}")
        for k, v in t.items():
            if k in ("name", "pass"):
                continue
            if isinstance(v, list) and v and isinstance(v[0], dict):
                lines.append(f"- **{k}**:")
                for row in v:
                    lines.append("    - " + ", ".join(
                        f"{a}={_fmt(b)}" for a, b in row.items()))
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
    out = here / "runs" / f"{ts}_congruence_probe"
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
