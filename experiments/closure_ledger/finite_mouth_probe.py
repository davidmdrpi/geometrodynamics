"""A finite-radius mouth: does PR #260's negative mode survive?

> Scope: a LINEAR conformally coupled scalar on a FIXED Einstein static universe.
> The mouths are now SPHERES of radius a, coupled to the tube through ONE channel
> each -- monopole matching, so the l >= 1 content on each sphere is dropped.
> They are spheres in a fixed ambient, not a solved neck geometry. NO
> BACKREACTION.

THE QUESTION
────────────
PR #260 built the conservative finite throat and found that with POINT mouths it
carries an exponentially growing mode for every choice of parameters,
generated at the point-mouth/tube INTERFACE -- and, in the sigma L, sigma d >> 1
limit, localizing to a single mouth at a rate sigma* = 2 sqrt(pi/A) that knows
neither the tube's length nor the mouths' separation. That gated the roadmap on
one question:

    DOES THE NEGATIVE MODE SURVIVE A FINITE-RADIUS MOUTH?

THE ANSWER IS NO -- and the reason is sharper than "it goes away".

WHAT CHANGES
────────────
A point mouth needs the 1/(4 pi chi) subtraction and keeps the RENORMALIZED
self-energy g(lam), which is negative. A sphere does not. Smearing the coupling
over a sphere of radius a -- the same operator on both sides, so the composite
stays manifestly self-adjoint -- replaces g by the UNSUBTRACTED Green function:

    G_self(lam) = f(a,lam) G(a,lam),     G_cross(lam) = f(a,lam)^2 G(d,lam)

with f(chi,lam) = sin(w chi)/(w sin chi) the regular radial solution, f(0) = 1.
Both are MEAN-VALUE IDENTITIES, and T2 checks them against direct quadrature on
S^3 rather than asserting them.

WHAT IS CHECKED
───────────────
T2  THE MODEL'S OWN IDENTITIES. <G(.,c2)> over the sphere at distance a from c1
    equals f(a) G(d) to 1.0e-10 -- exact, by the mean-value theorem, computed a
    second way. The self term equals f(a) G(a) to 4.1e-04, grid-limited by the
    Green function's integrable singularity at coincidence and reported as
    quadrature error rather than as a model error.

T3  *** THE NEGATIVE MODE DOES NOT SURVIVE, AND IT IS STRUCTURAL. *** At
    lam = -sigma^2 the two sides of det(A - G) = 0 have OPPOSITE SIGNS for every
    admissible parameter choice:

      the tube      -coth(kL/2)/(A k),  -tanh(kL/2)/(A k)     strictly NEGATIVE
      the ambient   f G(a) +/- f^2 G(d)                       strictly POSITIVE

    the second because f and G are positive on the imaginary axis and
    G(a) > f(a) G(d) whenever a < d/2, which disjoint mouths require anyway. A
    difference of a negative and a positive number has no zero. 3078 samples
    over (a, d, L, A, m, sigma): 0 positives, worst approach -5.1e-04.

T4  AND PR #260's MODE WAS THE LINEARIZATION. That round wrote the mouth as a
    CONSTANT shift 1/(4 pi a), the leading term of
    G(a,lam) = 1/(4 pi a) + g(lam) + O(a). The exact G(a,-k^2) is SCREENED,
    ~ e^{-k a}/(4 pi a), and dies; the constant does not, so it eventually beats
    the tube's -1/(A k) and crosses zero. Measured: the two agree to 0.8% for
    k a <= 0.1 and disagree by 1000% at k a = 3, and the linearized root sits at
    k a = 1.0004, 1.0025, 1.022, 1.127 for a = 0.02, 0.05, 0.15, 0.35 -- i.e. AT
    THE EDGE OF ITS OWN VALIDITY. A mode living where its derivation fails is an
    artifact, and this is the demonstration rather than the suspicion.

T5  WHERE THE MODE WENT: SOFT AND POSITIVE. The composite has exactly one state
    below the free gap lam = 1, in the symmetric channel, and as a -> 0

        lam_0  ->  8 pi a / (A L)

    two mouth capacitances 4 pi a restoring a tube of volume A L, matched to 1%
    by a = 0.02 and 0.2% by a = 0.005. So the point limit drives the mode to zero
    FROM ABOVE. PR #260 did not get a rate slightly wrong; it put a mode that
    approaches 0+ on the other side of zero, at lam ~ -1/a^2.

T6  THE DELAY SURVIVES. Slope 1.0010 in L, saturation at the ambient path d to
    0.0, with the mouth contributing only a sub-leading O(a) shift (measured
    slope -0.39). A first draft predicted -2a from an ambient block missing the
    shell form factor; the measured slope is quoted and the prediction recorded
    as wrong. The contour is also easier: with no growing mode to clear, eps is
    back to PR #259's single requirement, and 0.4 is comfortable where PR #260
    needed eps > 2.

T7  AND SO DO THE STATIC RESULTS. det S -> 0 linearly in lam (rank one), and PR
    #258's W = -beta(lam) holds to 3.6e-12. Both came from the TUBE's zero mode,
    which the mouth does not touch; only the coefficient moves, because the
    ambient diagonal is now +f(a)G(a) instead of the renormalized g_0 < 0.

T8  WHAT IS STILL PUT IN, QUANTIFIED. One channel per mouth couples only l = 0.
    The dropped multipoles obey PR #250's screening law -- the dipole/monopole
    ratio is 0.934 (a/d) across a decade in a -- so the leading omission is the
    dipole at O(a/d), and the dropped power fraction is 6.9e-05 at a = 0.02.

WHAT THIS UNGATES
─────────────────
PR #260 gated stationary action and backreaction on this question. The answer
releases them: the finite-mouth throat is conservative, has a real traversal
delay, and has NO growing mode, so an on-shell action or an A/B/A+B collapse
comparison computed on it measures the physics rather than an instability. What
it does NOT settle is the neck geometry -- these are spheres in a fixed ambient,
with monopole coupling and no backreaction.

    python -m experiments.closure_ledger.finite_mouth_probe
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.finite_mouth import (
    WORKING_MOUTH,
    measure_monopole_matching_is_the_remaining_approximation,
    measure_the_delay_survives_with_a_radius_correction,
    measure_the_instability_was_the_linearization,
    measure_the_mean_value_identities_hold,
    measure_the_mode_became_soft_and_positive,
    measure_the_negative_mode_does_not_survive,
    measure_the_static_results_survive,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("PR #260 gated the roadmap on one question: does its "
                     "negative mode survive a finite-radius mouth? Replace the "
                     "point mouths by spheres of radius a, coupled through one "
                     "channel each, and settle it."),
        "scope": ("a linear conformally coupled scalar on a fixed Einstein "
                  "static universe. The mouths are spheres in a FIXED ambient "
                  "with MONOPOLE coupling -- not a solved neck geometry -- so "
                  "the l >= 1 content on each sphere is dropped, quantified in "
                  "T8. No backreaction."),
        "pass": True,
    }


def t2_the_mean_value_identities() -> dict:
    r = measure_the_mean_value_identities_hold()
    return {"name": "T2_the_mean_value_identities", **r,
            "pass": bool(r["the_cross_identity_is_exact"]
                         and r["the_self_identity_is_grid_limited"])}


def t3_the_negative_mode_does_not_survive() -> dict:
    r = measure_the_negative_mode_does_not_survive()
    return {"name": "T3_the_negative_mode_does_not_survive", **r,
            "pass": bool(r["there_is_no_growing_mode"])}


def t4_the_instability_was_the_linearization() -> dict:
    r = measure_the_instability_was_the_linearization()
    return {"name": "T4_the_instability_was_the_linearization", **r,
            "pass": bool(r["every_linearized_model_has_a_root"]
                         and r["no_exact_model_has_one"]
                         and r["the_root_sits_at_kappa_a_of_order_one"])}


def t5_the_mode_became_soft_and_positive() -> dict:
    r = measure_the_mode_became_soft_and_positive()
    return {"name": "T5_the_mode_became_soft_and_positive", **r,
            "pass": bool(r["every_one_is_positive"]
                         and r["every_one_is_below_the_gap"]
                         and r["the_capacitance_formula_holds"]
                         and r["the_point_limit_approaches_zero_from_above"])}


def t6_the_delay_survives() -> dict:
    r = measure_the_delay_survives_with_a_radius_correction()
    return {"name": "T6_the_delay_survives", **r,
            "pass": bool(r["the_traversal_time_survives"]
                         and r["the_ambient_path_still_takes_over"]
                         and r["the_mouth_shift_is_subleading"])}


def t7_the_static_results_survive() -> dict:
    r = measure_the_static_results_survive()
    return {"name": "T7_the_static_results_survive", **r,
            "pass": bool(r["the_static_response_is_still_rank_one"]
                         and r["det_S_is_linear_in_lambda"]
                         and r["the_defect_is_still_minus_beta"])}


def t8_monopole_matching_is_what_is_left() -> dict:
    r = measure_monopole_matching_is_the_remaining_approximation()
    return {"name": "T8_monopole_matching_is_what_is_left", **r,
            "pass": bool(r["dipole_scales_like_a_over_d"]
                         and r["smallest_dropped_fraction"] < 1e-3)}


def t9_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T9_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_mean_value_identities(),
             t3_the_negative_mode_does_not_survive(),
             t4_the_instability_was_the_linearization(),
             t5_the_mode_became_soft_and_positive(),
             t6_the_delay_survives(),
             t7_the_static_results_survive(),
             t8_monopole_matching_is_what_is_left()]
    tests.append(t9_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8 = tests[1:8]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_NEGATIVE_MODE_DOES_NOT_SURVIVE_A_FINITE_MOUTH"
        verdict = (
            "*** THE ANSWER IS NO. *** PR #260 gated the roadmap on whether its "
            "growing mode survives a finite-radius mouth. It does not, and the "
            "statement is STRUCTURAL rather than parametric. Smearing the "
            "coupling over a sphere of radius a -- the same operator on both "
            "sides, so the composite stays manifestly self-adjoint -- replaces "
            "the renormalized self-energy g(lam), which is negative, by the "
            "UNSUBTRACTED Green function G(a,lam), which is positive. The "
            "mean-value identities that give it, G_self = f(a)G(a) and "
            "G_cross = f(a)^2 G(d), are checked against direct quadrature on "
            f"S^3 to {t2.get('worst_cross_error', 0):.1e} for the cross term "
            f"and {t2.get('worst_self_error', 0):.1e} for the self term, the "
            "latter grid-limited by the integrable singularity at coincidence. "
            "THEN THE SIGNS DECIDE IT: at lam < 0 the tube's channel functions "
            "are strictly NEGATIVE -- a passive interior supplies no restoring "
            "force -- and the ambient's are strictly POSITIVE, because f and G "
            "are positive on the imaginary axis and G(a) > f(a)G(d) whenever "
            "a < d/2, which disjoint mouths require anyway. A difference of a "
            "negative and a positive number has no zero, so no parameter choice "
            f"can produce a growing mode: {t3.get('samples', 0)} samples over "
            "(a, d, L, A, m, sigma) give "
            f"{t3.get('positives', 0)} positives, worst approach "
            f"{t3.get('worst_approach_to_zero', 0):.1e}. *** AND PR #260's MODE "
            "WAS THE LINEARIZATION. *** That round modelled the mouth as a "
            "CONSTANT shift 1/(4 pi a) -- the leading term of "
            "G(a,lam) = 1/(4 pi a) + g(lam) + O(a). The exact G(a,-k^2) is "
            "SCREENED, ~ e^{-k a}/(4 pi a), and dies; the constant does not, so "
            "it eventually beats the tube's -1/(A k) and crosses. The two agree "
            f"to {100*t4.get('worst_gap_below_kappa_a_of_0p1', 0):.1f}% for "
            f"k a <= 0.1 and differ by "
            f"{100*t4.get('worst_gap_at_kappa_a_of_3', 0):.0f}% at k a = 3, and "
            "the linearized root sits at k a = "
            f"{', '.join(f'{v:.3f}' for v in t4.get('linearized_roots_times_radius', []))}"
            " -- AT THE EDGE OF ITS OWN VALIDITY. A mode that lives at the "
            "scale where its derivation fails is an artifact, and that is now "
            "shown rather than suspected. WHERE THE MODE WENT: SOFT AND "
            "POSITIVE. The composite has exactly one state below the free gap "
            "lam = 1, in the symmetric channel, and as a -> 0 it goes to "
            "8 pi a/(A L) -- two mouth capacitances 4 pi a restoring a tube of "
            f"volume A L -- matched to "
            f"{100*t5.get('worst_closed_form_error', 0):.1f}%. So the point "
            "limit drives the mode to zero FROM ABOVE. PR #260 did not get a "
            "rate slightly wrong; it took a mode approaching 0+ and put it on "
            "the other side of zero, at lam ~ -1/a^2. THE GOOD RESULTS SURVIVE: "
            f"the traversal delay has slope {t6.get('slope_in_length', 0):.4f} "
            "in L against a predicted 1, saturating at the ambient path to "
            f"{t6.get('onset_spread_above_d', 0):.1e}, with the mouth adding "
            f"only a sub-leading O(a) shift (slope "
            f"{t6.get('slope_in_radius', 0):.2f}; a first draft predicted -2a "
            "from an ambient block missing the shell form factor, and that "
            "prediction is recorded as wrong). The static response is still "
            "rank one, det S linear in lam, and PR #258's W = -beta(lam) holds "
            f"to {t7.get('worst_defect_error', 0):.1e} -- all three came from "
            "the TUBE's zero mode, which the mouth does not touch. The contour "
            "is easier too: with nothing to clear, eps = 0.4 suffices where PR "
            "#260 needed eps > 2. WHAT IS STILL PUT IN: one channel per mouth, "
            "so only l = 0 couples; the dropped multipoles obey PR #250's "
            "screening law, dipole/monopole = "
            f"{t8.get('rows', [{}])[0].get('ratio_over_a_over_d', 0):.3f} (a/d) "
            "across a decade in a, and the dropped power fraction is "
            f"{t8.get('smallest_dropped_fraction', 0):.1e} at a = 0.02. THE "
            "MOUTHS ARE SPHERES IN A FIXED AMBIENT, NOT A SOLVED NECK, AND "
            "THERE IS NO BACKREACTION. But the gate PR #260 set is answered: "
            "with the mouth resolved there is no growing mode, so the "
            "stationary-action and backreaction rounds can proceed on this "
            "background and measure the physics rather than an artifact.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "finite_mouth", "tests": tests,
            "verdict_class": verdict_class, "verdict": verdict,
            "working_mouth": {"separation": WORKING_MOUTH.separation,
                              "length": WORKING_MOUTH.length,
                              "area": WORKING_MOUTH.area,
                              "radius": WORKING_MOUTH.radius},
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
    out = here / "runs" / f"{ts}_finite_mouth_probe"
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
