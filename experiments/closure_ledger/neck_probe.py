"""The neck: the balls actually removed, and the binary question re-asked.

> Scope: a LINEAR conformally coupled scalar on a FIXED round S^3 with two balls
> cut out. The tube has ONE transverse channel, so A is a coupling and the neck
> is a quantum-graph edge rather than a solved geometry. The metric is not
> solved for. NO BACKREACTION.

THE QUESTION
────────────
PR #260 found a growing mode with POINT mouths and gated the roadmap on whether
it survives a finite-radius mouth. PR #261 answered NO -- but with one stated
limitation, and it is the load-bearing one: those mouths were SPHERES IN A FIXED
AMBIENT. The balls were never removed. That round smeared the coupling over
dB_a while still using the WHOLE sphere's Green function, and coupled through
ONE channel per mouth, so only l = 0 talked to the tube.

This round removes the balls. The ambient is Omega = S^3 \\ (B_a(c1) u B_a(c2)),
the tube is glued along the boundary spheres, and the question is re-asked:

    DOES THE NEGATIVE MODE SURVIVE?

THE ANSWER IS STILL NO -- and here it is a THEOREM rather than a sweep.

THE THEOREM
───────────
With the balls removed there is NO SUBTRACTION ANYWHERE, so the composite energy
is a sum of manifestly non-negative terms:

    E[phi,u] = Int_Omega (|grad phi|^2 + phi^2) dV  +  A Int_0^L (|u'|^2 + m^2
    |u|^2) ds

Every term is >= 0, and E = 0 forces phi = 0 on Omega -- whereupon the matching
gives u(0) = u(L) = 0 and Poincare gives A Int |u'|^2 >= (pi/L)^2 A Int |u|^2, so
u = 0 too. Hence

    lam = E / ||.||^2  >  0     for EVERY configuration,
                               ALL MULTIPOLES, NO TRUNCATION, NO SWEEP.

That is the whole argument, and it is why removing the balls is not a refinement
of PR #261 but a CHANGE OF FOOTING. That round had to establish a sign on a
reduced 2x2; this one has positivity of the form itself. The renormalized
g(lam) < 0 that PR #260's mode fed on has nowhere to enter, because nothing is
renormalized.

WHAT IS CHECKED
───────────────
T2  THE EXTERIOR MAP, AND THE INTEGRATOR THAT COMPUTES IT. The exterior
    Dirichlet-to-Neumann map N_l(lam) = -4 pi sin^2 a psi'(a)/psi(a) is obtained
    by shooting v'' + [lam - l(l+2)/sin^2 chi] v = 0 from the far pole, where
    regularity picks the e^{l+1} branch. In the l = 0 channel it agrees with an
    INDEPENDENT closed form to 1.7e-14, which is the check the l >= 1 channels
    inherit. N_l > 0 for every l and every lam < 0, and INCREASES with l, so the
    monopole is the softest channel. N_0 -> 4 pi a, the capacitance of a sphere:
    that is what fixes the normalization as physical rather than conventional.

T3  *** THE THEOREM, AND A CHECK THAT THE OBJECT OBEYS IT. *** Random trial
    configurations are evaluated with the honest measures -- 4 pi sin^2 chi dchi
    outside, A ds inside -- and every Rayleigh quotient is positive and sits
    ABOVE the computed lowest mode, which is what a variational bound must do.
    The degenerate case is checked exactly: the purely interior configuration
    lands on the Poincare floor pi^2/L^2 to 2e-07.

T4  THE BINARY QUESTION, RE-ASKED. Still no. The sweep has the same shape as PR
    #261's and reaches the same place from a DIFFERENT CONSTRUCTION -- the
    diagonal is now 1/N_0, the inverse EXTERIOR map, rather than a smeared
    average over an ambient that still contained the balls. Two independent
    models, one answer.

T5  PR #261's MONOPOLE TRUNCATION WAS NEVER A STABILITY LIMITATION. A tube with
    one transverse channel presents a single number at each mouth, so it drives
    the l = 0 projection and nothing else. The l >= 1 sectors are the exterior's
    own modes with its own boundary condition, their DtN is positive, and they
    are at least 1.62x stiffer than the monopole. They can neither be driven nor
    go unstable. PR #250's (a/d)^l screening bounds missed AMPLITUDE, not the
    stability answer.

T6  WHAT THE FIXED AMBIENT COST, IN NUMBERS. PR #261's f(a)G(a) against this
    round's 1/N_0(a,lam): they agree to 1.3e-04 at a = 0.02, lam = 0, to
    3.8e-03 at the working radius a = 0.05, and part company at large radius and
    deep lam, reaching 11% at a = 0.35, lam = -4. The error tracks the fraction
    of the sphere wrongly left in. A limitation with a measured size is a
    different object from one with a caveat.

    The same test prices the approximation this round has NOT removed: the
    reduced 2x2 keeps a SINGLE-SCATTERING cross term f(a)^2 G(d), where the
    exact two-ball exterior would add the series in which each boundary sphere
    drives the other. Its expansion parameter is cross/self, measured as
    0.80-0.90 times (a/d) across a decade in a -- 0.031 at the working point, so
    the neglected terms enter at 9.5e-04 there and never exceed 0.16 anywhere
    sampled. A correction that small cannot flip a sign, which is the only thing
    the reduced model is asked to decide; and the theorem does not go through the
    reduced model at all.

T7  THE SOFT MODE IS FORCED, NOT FOUND -- AND ONE CORRECTION TO PR #261. The
    same style of argument that kills the growing mode produces this one, from
    the two ends of the gap: F_sym -> +inf as lam -> 0+ (the tube's 2/(A lam L))
    and F_sym -> -inf as lam -> 1-, because the exterior stiffness VANISHES at
    the free ESU threshold,

        N_0  ->  2 pi (pi - a + sin a cos a) (1 - lam)

    reproduced to 3.1e-05 with first-order convergence (error ratio 0.1000 per
    decade). A continuous function running from +inf to -inf has a zero.
    *** AND: PR #261's "exactly one state below the gap" is a statement about
    L < pi, not a structural one. *** The channel functions have POLES at
    lam = (2 pi j/L)^2 and (pi(2j-1)/L)^2; above L = pi these enter the gap and
    each brings a genuine extra state just above it -- three states at L = 8. A
    pole is a SIGN CHANGE WITH NO ZERO, so crossing-counting alone reports
    states that are not there; classifying by the residual separates roots from
    poles by fifteen orders of magnitude. None of this touches the stability
    answer: every one of those states is positive, as the form requires.

T8  THE SOFT MODE SURVIVES THE REMOVAL. lam_0 moves by 1.2e-04 relative at
    a = 0.02 and 7.0e-04 at the working radius, and still tends to
    8 pi a/(A L) -- two mouth capacitances restoring a tube of volume A L. The
    capacitance reading is cleaner here, since N_0 -> 4 pi a is now the exterior
    map's own small-ball limit rather than an identification made after the
    fact.

WHAT THIS UNGATES
─────────────────
PR #260 gated stationary action and backreaction on the negative mode. PR #261
released that gate on a reduced model with the balls left in; this round closes
the gap the release stood on, and upgrades the answer from a sign on a 2x2 to
positivity of the energy itself. What it does NOT settle: the tube still has one
transverse channel, so A is a coupling rather than a solved neck cross-section;
the ambient metric is FIXED; and there is NO BACKREACTION. Those are the next
round's subject, not this one's caveat.

    python -m experiments.closure_ledger.neck_probe
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.neck import (
    WORKING_NECK,
    measure_the_exterior_dtn_is_positive_in_every_channel,
    measure_the_higher_multipoles_decouple,
    measure_the_negative_mode_does_not_survive_the_neck,
    measure_the_quadratic_form_is_positive,
    measure_the_soft_mode_is_forced_by_the_two_ends,
    measure_the_soft_mode_survives_the_removal,
    measure_what_the_fixed_ambient_cost,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("PR #261 answered PR #260's gate with spheres in a FIXED "
                     "ambient -- the balls were never removed, and only l = 0 "
                     "coupled. Remove them: solve on S^3 minus two balls, glue "
                     "the tube along the boundary spheres, and re-ask whether "
                     "the negative mode survives."),
        "scope": ("a linear conformally coupled scalar on a FIXED round S^3 "
                  "with two balls cut out. The tube has ONE transverse "
                  "channel, so A is a coupling and the neck is a quantum-graph "
                  "edge rather than a solved geometry. No backreaction."),
        "pass": True,
    }


def t2_the_exterior_map() -> dict:
    r = measure_the_exterior_dtn_is_positive_in_every_channel()
    return {"name": "T2_the_exterior_map", **r,
            "pass": bool(r["it_is_positive_everywhere"]
                         and r["it_increases_with_ell"]
                         and r["the_capacitance_normalization_holds"]
                         and r["worst_shooting_vs_closed_form"] < 1e-11)}


def t3_the_quadratic_form_is_positive() -> dict:
    r = measure_the_quadratic_form_is_positive()
    return {"name": "T3_the_quadratic_form_is_positive", **r,
            "pass": bool(r["every_quotient_is_positive"]
                         and r["no_trial_beats_the_lowest_mode"]
                         and r["the_interior_only_case_hits_the_poincare_floor"])}


def t4_the_negative_mode_does_not_survive() -> dict:
    r = measure_the_negative_mode_does_not_survive_the_neck()
    return {"name": "T4_the_negative_mode_does_not_survive", **r,
            "pass": bool(r["there_is_no_growing_mode"])}


def t5_the_higher_multipoles_decouple() -> dict:
    r = measure_the_higher_multipoles_decouple()
    return {"name": "T5_the_higher_multipoles_decouple", **r,
            "pass": bool(r["every_channel_is_positive"]
                         and r["every_channel_is_stiffer_than_the_monopole"])}


def t6_what_the_fixed_ambient_cost() -> dict:
    r = measure_what_the_fixed_ambient_cost()
    return {"name": "T6_what_the_fixed_ambient_cost", **r,
            "pass": bool(r["pr_261_was_right_where_it_looked"]
                         and r["the_error_grows_with_the_radius"]
                         and r["cross_over_self_scales_like_a_over_d"]
                         and r["the_neglected_series_cannot_reach_the_sign"])}


def t7_the_soft_mode_is_forced() -> dict:
    r = measure_the_soft_mode_is_forced_by_the_two_ends()
    return {"name": "T7_the_soft_mode_is_forced", **r,
            "pass": bool(
                r["the_symmetric_channel_starts_positive"]
                and r["and_ends_negative"]
                and r["so_a_symmetric_state_is_forced"]
                and r["the_exterior_stiffness_vanishes_linearly_at_the_gap_edge"]
                and r["short_tubes_hold_exactly_one_state"]
                and r["long_tubes_hold_more"]
                and r["every_extra_state_follows_a_tube_harmonic"])}


def t8_the_soft_mode_survives_the_removal() -> dict:
    r = measure_the_soft_mode_survives_the_removal()
    return {"name": "T8_the_soft_mode_survives_the_removal", **r,
            "pass": bool(r["every_one_is_positive_and_below_the_gap"]
                         and r["the_capacitance_formula_survives"]
                         and r["removing_the_balls_barely_moves_it"])}


def t9_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T9_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_exterior_map(),
             t3_the_quadratic_form_is_positive(),
             t4_the_negative_mode_does_not_survive(),
             t5_the_higher_multipoles_decouple(),
             t6_what_the_fixed_ambient_cost(),
             t7_the_soft_mode_is_forced(),
             t8_the_soft_mode_survives_the_removal()]
    tests.append(t9_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8 = tests[1:8]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_NEGATIVE_MODE_DOES_NOT_SURVIVE_THE_NECK"
        verdict = (
            "*** THE ANSWER IS STILL NO, AND IT IS NOW A THEOREM. *** PR #261 "
            "answered PR #260's gate with spheres in a FIXED ambient: the balls "
            "were never removed, and only l = 0 coupled. That was its own "
            "stated limitation and the strongest remaining doubt. This round "
            "removes them -- the ambient is S^3 minus two balls, the tube is "
            "glued along the boundary spheres -- and the answer does not move. "
            "WITH THE BALLS REMOVED THERE IS NO SUBTRACTION ANYWHERE, so the "
            "energy E = Int_Omega (|grad phi|^2 + phi^2) + A Int (|u'|^2 + m^2 "
            "|u|^2) is a sum of manifestly non-negative terms, and E = 0 forces "
            "phi = 0, whence matching gives u(0) = u(L) = 0 and Poincare gives "
            "A Int |u'|^2 >= (pi/L)^2 A Int |u|^2. Hence lam > 0 for EVERY "
            "configuration, ALL MULTIPOLES, NO TRUNCATION, NO SWEEP. That is a "
            "change of footing rather than a refinement: PR #261 had to "
            "establish a SIGN on a reduced 2x2, this round has POSITIVITY OF "
            "THE FORM ITSELF, and the renormalized g(lam) < 0 that PR #260's "
            "mode fed on has nowhere to enter because nothing is renormalized. "
            "The object the theorem is about is checked rather than asserted: "
            "the exterior DtN is computed by shooting the radial equation from "
            "the far pole and agrees with an independent closed form to "
            f"{t2.get('worst_shooting_vs_closed_form', 0):.1e}; it is positive "
            "in every channel and INCREASES with l, so the monopole is the "
            "softest and the higher channels cannot be the first to go soft; "
            "and N_0 -> 4 pi a, the capacitance of a sphere, which fixes the "
            "normalization as physical. Explicit trial configurations give "
            f"Rayleigh quotients from {t3.get('smallest_quotient', 0):.4f} up, "
            "all positive and all above the computed lowest mode "
            f"{t3.get('lowest_computed_mode', 0):.6f}, and the degenerate "
            "purely-interior case lands on the Poincare floor pi^2/L^2 exactly. "
            f"The sweep agrees with the theorem: {t4.get('samples', 0)} samples "
            f"give {t4.get('positives', 0)} positives, worst approach "
            f"{t4.get('worst_approach_to_zero', 0):.1e}, with the tube side "
            "negative and the exterior side positive at every one. *** AND PR "
            "#261's MONOPOLE TRUNCATION WAS NEVER A STABILITY LIMITATION. *** A "
            "one-channel tube presents a single number at each mouth, so it "
            "drives only l = 0; the higher sectors are the exterior's own modes "
            "with positive DtN, at least "
            f"{t5.get('smallest_ratio_to_the_monopole', 0):.2f}x stiffer than "
            "the monopole, and can neither be driven nor go unstable. PR #250's "
            "(a/d)^l screening bounds missed AMPLITUDE, not the answer. WHAT "
            "THE FIXED AMBIENT COST, PRICED: PR #261's f(a)G(a) against this "
            "round's 1/N_0 agrees to "
            f"{t6.get('worst_error_at_the_working_radius', 0):.1e} at the "
            "working radius and departs to "
            f"{t6.get('worst_error_overall', 0):.1%} at a = 0.35, lam = -4, "
            "growing with the fraction of the sphere wrongly left in -- so that "
            "round was right where it looked, and now by a measured margin "
            "rather than a caveat. THE ONE APPROXIMATION LEFT IN THE REDUCED "
            "2x2 IS ALSO PRICED: its cross term f(a)^2 G(d) is "
            "SINGLE-SCATTERING, and the neglected series expands in cross/self "
            "= 0.8 (a/d), which is "
            f"{t6.get('neglected_scattering_at_the_working_point', 0):.1e} at "
            "the working point and at worst "
            f"{t6.get('neglected_scattering_worst', 0):.2f} anywhere sampled -- "
            "too small to flip the sign the reduced model is asked to decide, "
            "and irrelevant to the theorem, which does not go through it. THE "
            "SOFT MODE IS FORCED, NOT FOUND: the same "
            "style of argument that kills the growing mode produces this one "
            "from the two ends of the gap, F_sym -> +inf as lam -> 0+ and "
            "-> -inf as lam -> 1-, the latter because the exterior stiffness "
            "VANISHES at the free ESU threshold, N_0 -> 2 pi (pi - a + sin a "
            f"cos a)(1 - lam), reproduced to {t7.get('gap_edge_slope_error', 0):.1e} "
            "with first-order convergence (error ratio "
            f"{t7.get('gap_edge_error_ratio_per_decade', 0):.4f} per decade). "
            "*** ONE CORRECTION TO PR #261: *** its 'exactly one state below "
            "the gap' is a statement about L < pi, not a structural one. The "
            "channel functions have POLES at lam = (2 pi j/L)^2 and "
            "(pi(2j-1)/L)^2; above L = pi these enter the gap and each brings a "
            "genuine extra state just above it -- three states at L = 8. A pole "
            "is a sign change with NO ZERO, so crossing-counting alone reports "
            f"states that are not there ({t7.get('pole_crossings_that_would_have_been_miscounted', 0)} "
            "of them across the lengths swept here); classifying by the "
            "residual separates roots from poles by fifteen orders of "
            "magnitude. None of it touches the stability answer -- every one of "
            "those states is positive, as the form requires. AND THE SOFT MODE "
            f"SURVIVES THE REMOVAL: lam_0 moves by "
            f"{t8.get('shift_at_the_working_radius', 0):.1e} at the working "
            "radius and still tends to 8 pi a/(A L), matched to "
            f"{t8.get('worst_closed_form_error', 0):.1%}. WHAT IS STILL PUT IN: "
            "the tube has ONE transverse channel, so A is a coupling rather "
            "than a solved neck cross-section; the ambient metric is FIXED; and "
            "there is NO BACKREACTION. PR #260's gate is now closed from both "
            "sides -- resolved mouths and removed balls -- so stationary action "
            "and backreaction can proceed on a background whose positivity is "
            "proved rather than sampled.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "neck", "tests": tests,
            "verdict_class": verdict_class, "verdict": verdict,
            "working_neck": {"separation": WORKING_NECK.separation,
                             "length": WORKING_NECK.length,
                             "area": WORKING_NECK.area,
                             "radius": WORKING_NECK.radius},
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
                for row in v[:30]:
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
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, np.ndarray):
        return np.array2string(v, precision=5)
    if isinstance(v, dict):
        return ", ".join(f"{a}={_fmt(b)}" for a, b in v.items())
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_fmt(x) for x in list(v)[:12]) + "]"
    return str(v)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_neck_probe"
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
