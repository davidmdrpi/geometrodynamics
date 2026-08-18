"""A finite conservative throat: a DtN map, a point limit, and a traversal time.

> Scope: a LINEAR scalar field on a FIXED round background (the Einstein static
> universe S3 x R), conformally coupled. The throat now has an INTERIOR -- a
> tube of proper length L, cross-section A and interior mass m -- but the
> interior is ONE-DIMENSIONAL (a single transverse channel) and the mouths are
> still POINTS in the ambient. NO BACKREACTION and no topology change: nothing
> here solves for the geometry that would produce this throat.

WHAT THIS IS FOR
────────────────
Every round from PR #253 to PR #259 carried the same disclaimer: the throat is
point-supported -- no interior, no proper length, no delay -- and what is solved
is a rank-one MOUTH-TRANSFER model: field values only, no normal-derivative
matching, no reflected channel, 1x1 where a conserving junction needs 2x2, and
lossy for kappa < 1. This round pays that debt.

WHAT IS PUT IN
──────────────
Two lines. The tube's Dirichlet-to-Neumann map,

    N(lam) = A k [[cot kL, -csc kL], [-csc kL, cot kL]],   k^2 = lam - m^2,

and the matching to the ambient: value continuity and flux continuity at the
mouths, q = -N Phi. Everything below follows.

The consequence that organizes the round is that the boundary condition is now
FREQUENCY-DEPENDENT. A point throat is a fixed Hermitian A; a finite throat is
A(lam) = -N(lam)^-1, and that dependence IS the interior.

WHAT IS CHECKED
───────────────
T2  THE BOUNDARY CONDITION IS SELF-ADJOINT AT EVERY FREQUENCY. rank[B|C] = 2
    and BC^dag = CB^dag to 0.0 at seven values of lam on both sides of zero.
    The DtN map is checked against the interior it summarizes, by Green's
    identity under quadrature, rather than against itself. The control is PR
    #255's rank-one transfer model, whose defect is 0.3 -- the size of the
    coupling itself. That is the kappa < 1 loss, as an operator property.

T3  *** THE THROAT TRANSMITS AT THE TRAVERSAL TIME. *** The measured object is
    the throat operator's OWN impulse response, R(w) = (A(w) - Gamma(w))^-1
    inverted along the retarded contour, so no source or observer geometry
    enters the answer. Two different predictions, both met:

      r_11 (same mouth in and out)  starts at t = 0    -- a wave that reaches a
                                                          mouth is partly
                                                          reflected INSTANTLY
      r_12 (opposite mouths)        starts at min(L,d)

    d(onset)/dL = 1.007 against a predicted 1 below the ambient path, and the
    onset stops moving above it, to 0.0. THE AMBIENT ALSO CONNECTS THE TWO
    MOUTHS, along a geodesic of length d, and that path is there whether or not
    the mouths are joined -- PR #258's cross-mouth channel and PR #259's beta=0
    control, now SEPARATED IN TIME rather than by rank counting. The point
    throat transmits at t = 0, which is what a point throat is.

T4  AND THE DELAY LEDGER IS A DERIVATION. On the contour cot x = -i - 2i sum
    e^{2ikx} and csc x = -2i sum e^{i(2k+1)x}, verified to 4.5e-16 and 1.7e-15:
    the same-mouth entry carries delays 0, 2L, 4L (an instantaneous reflection
    and its echoes) and the cross-mouth entry L, 3L, 5L. The parities are the
    physics, and the reflected channel is the one the rank-one model does not
    have at all.

T5  THERE IS NO POINT LIMIT -- THERE IS A BAND, AND ITS WIDTH IS 1/L. Freezing
    A at A(lam_0) is exact at lam_0 and nowhere else: the transmission
    amplitude is 4% out at lam = 1.05 lam_0 and 121% out at 3 lam_0. The two
    channels fail differently. The antisymmetric one HAS a limit,
    A_anti -> -L/(2A) with an O(L^2) error measured over a decade in L; the
    symmetric one DIVERGES like 2/(A lam L), because a massless tube holds a
    zero mode and a point cannot.

T6  THE STATIC LIMIT IS RANK ONE, AND PR #258's TOMOGRAPHY BREAKS ON IT. The
    same zero mode empties the symmetric channel at lam = 0: S = Re R collapses
    onto [[1,-1],[-1,1]] to 4e-5, det S -> 0 LINEARLY in lam with coefficient
    149.08, and the disconnection defect W = S_12/det S - G_0 diverges like
    1/lam. A point throat is statically rank two; a massless finite throat is
    rank one. That is a FALSIFIABLE difference from a static measurement.
    Give the tube an interior mass and the rank returns with det S ~ -148.7 m^2,
    and then the closure: off the collapse W = -beta(lam) EXACTLY, to 3.1e-13,
    with beta the tube's own transmission amplitude. PR #258's theorem survives
    the generalization; what it returns is no longer a constant.

T7  AN INTERIOR MASS IS A TRANSMISSION CUTOFF. Below lam = m^2 the tube is
    evanescent: beta -> -csch(kappa L)/(A kappa), matched by its exponential
    asymptote to 7.6e-05, monotone, and with ZERO sign changes against 7 above
    the cutoff. Suppression 3.1e-03. The cutoff is also where the rank
    collapses -- the two statements are one statement read at lam = m^2.

T8  THE MODEL ALWAYS HAS A GROWING MODE, AND IT IS THE MOUTH'S. A(lam)
    decreases and Gamma(lam) increases (PR #257's Gram identity), so A - Gamma
    is strictly monotone between poles and each channel has AT MOST ONE root --
    a count, not a scan. The symmetric channel always has exactly one, at
    lam < 0. Three facts identify it, all as limits with their convergence
    measured: its rate matches sigma* = 2 sqrt(pi/A) to 1.5e-03 with NO L in it;
    two mouth separations agree to 3.9e-09; and the channel splitting is
    1.04 e^{-sigma* d}, the Euclidean propagator between the mouths, so at that
    scale the mouths do not know about each other. Its length scale is
    sqrt(A/4pi) -- the mouth's own radius, exactly where "point mouth" stops
    being an approximation.

T9  SO THE CONTOUR MUST CLEAR IT. Placed 0.03 BELOW sigma*, the inversion
    returns a field with support before its own light cone: a pedestal at 99%
    of the peak, onset 0.0 for an event that cannot begin until t = 0.6. Placed
    above, the pedestal is 1.0e-16. Same species as PR #259's under-resolved
    contour, and reported the same way -- both values and the rule. sigma* has
    a closed form, so the contour can be placed BEFORE the solve.

WHAT IS NOT CLOSED
──────────────────
The mouths are still points, and T8 says exactly what that costs. The interior
is one-dimensional, so A enters as a coupling and not as a geometry. And the
throat is still a fixed background: PR #261's common action and PR #262's A/B/
A+B metric backreaction are the next two steps, and this round hands them an
object with a proper length, a delay, and a conservation law.

    python -m experiments.closure_ledger.finite_throat_probe
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.finite_throat import (
    WORKING_THROAT,
    measure_the_boundary_condition_is_self_adjoint_at_every_frequency,
    measure_the_contour_must_clear_the_growing_mode,
    measure_the_delay_ledger_is_the_bounce_series,
    measure_the_growing_mode_belongs_to_the_mouth,
    measure_the_interior_mass_is_a_transmission_cutoff,
    measure_the_point_throat_is_a_single_frequency_match,
    measure_the_static_limit_is_rank_one_and_the_defect_diverges,
    measure_the_throat_transmits_at_the_traversal_time,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("replace the point-supported throat of PRs #253-#259 with "
                     "a FINITE conservative one -- a tube with a "
                     "Dirichlet-to-Neumann map -- and measure what the point "
                     "idealization was throwing away: a traversal delay, a "
                     "reflected channel, and the rank of the static response."),
        "scope": ("a linear conformally coupled scalar on a fixed Einstein "
                  "static universe. The tube's interior is one-dimensional and "
                  "its mouths are still points in the ambient, which is the "
                  "one approximation the round then measures the cost of. No "
                  "backreaction, no topology change, and L, A and m are chosen "
                  "rather than derived."),
        "pass": True,
    }


def t2_the_boundary_condition_is_self_adjoint() -> dict:
    r = measure_the_boundary_condition_is_self_adjoint_at_every_frequency()
    return {"name": "T2_the_boundary_condition_is_self_adjoint", **r,
            "pass": bool(r["the_finite_throat_is_conservative"]
                         and r["the_control_is_not"]
                         and r["worst_green_identity_residual"] < 1e-5)}


def t3_the_throat_transmits_at_the_traversal_time() -> dict:
    r = measure_the_throat_transmits_at_the_traversal_time()
    return {"name": "T3_the_throat_transmits_at_the_traversal_time", **r,
            "pass": bool(r["the_delay_is_the_traversal_time"]
                         and r["the_ambient_path_takes_over"]
                         and r["reflection_is_instantaneous"]
                         and r["the_point_throat_transmits_instantly"])}


def t4_the_delay_ledger_is_the_bounce_series() -> dict:
    r = measure_the_delay_ledger_is_the_bounce_series()
    return {"name": "T4_the_delay_ledger_is_the_bounce_series", **r,
            "pass": bool(r["the_series_converge_on_the_contour"])}


def t5_the_point_throat_is_a_single_frequency_match() -> dict:
    r = measure_the_point_throat_is_a_single_frequency_match()
    return {"name": "T5_the_point_throat_is_a_single_frequency_match", **r,
            "pass": bool(r["the_band_error_reaches_one"]
                         and r["the_antisymmetric_channel_has_a_limit"]
                         and r["the_symmetric_channel_diverges"])}


def t6_the_static_limit_is_rank_one() -> dict:
    r = measure_the_static_limit_is_rank_one_and_the_defect_diverges()
    return {"name": "T6_the_static_limit_is_rank_one", **r,
            "pass": bool(r["the_static_response_is_rank_one"]
                         and r["det_S_is_linear_in_lambda"]
                         and r["the_defect_diverges"]
                         and r["the_defect_is_still_minus_beta"])}


def t7_the_interior_mass_is_a_cutoff() -> dict:
    r = measure_the_interior_mass_is_a_transmission_cutoff()
    return {"name": "T7_the_interior_mass_is_a_cutoff", **r,
            "pass": bool(r["the_evanescent_side_does_not_oscillate"]
                         and r["the_transmission_is_suppressed_below_cutoff"]
                         and r["everything_stays_real"])}


def t8_the_growing_mode_belongs_to_the_mouth() -> dict:
    r = measure_the_growing_mode_belongs_to_the_mouth()
    return {"name": "T8_the_growing_mode_belongs_to_the_mouth", **r,
            "pass": bool(r["every_throat_has_one"]
                         and r["it_stops_knowing_the_separation"]
                         and r["the_closed_form_holds_once_sigma_L_is_large"]
                         and r["the_split_is_the_euclidean_propagator"])}


def t9_the_contour_must_clear_it() -> dict:
    r = measure_the_contour_must_clear_the_growing_mode()
    return {"name": "T9_the_contour_must_clear_it", **r,
            "pass": bool(r["a_contour_below_the_mode_breaks_causality"])}


def t10_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T10_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_boundary_condition_is_self_adjoint(),
             t3_the_throat_transmits_at_the_traversal_time(),
             t4_the_delay_ledger_is_the_bounce_series(),
             t5_the_point_throat_is_a_single_frequency_match(),
             t6_the_static_limit_is_rank_one(),
             t7_the_interior_mass_is_a_cutoff(),
             t8_the_growing_mode_belongs_to_the_mouth(),
             t9_the_contour_must_clear_it()]
    tests.append(t10_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8, t9 = tests[1:9]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_THROAT_HAS_AN_INTERIOR_AND_A_TRAVERSAL_TIME"
        verdict = (
            "THE THROAT NOW HAS AN INTERIOR, AND THE INTERIOR IS THE DELAY. "
            "PRs #253-#259 all carried the same disclaimer -- point-supported, "
            "no proper length, no delay, a rank-one mouth-transfer model that "
            "is lossy for kappa < 1 -- and this round replaces it with a tube "
            "of length L and cross-section A whose Dirichlet-to-Neumann map is "
            "exact. Two lines are put in: N(lam) = A k [[cot kL, -csc kL], "
            "[-csc kL, cot kL]] and the matching q = -N Phi. Everything else "
            "follows. THE BOUNDARY CONDITION IS SELF-ADJOINT AT EVERY "
            f"FREQUENCY: rank[B|C] = 2 and BC^dag - CB^dag = "
            f"{t2.get('worst_hermiticity_defect', 0):.1e} at seven values of "
            "lam on both sides of zero, with the DtN map checked against the "
            "interior it summarizes by Green's identity to "
            f"{t2.get('worst_green_identity_residual', 0):.1e} rather than "
            "against itself; PR #255's rank-one transfer model, the control, "
            f"has defect {t2.get('the_rank_one_control_defect', 0):.2f} -- the "
            "size of the coupling itself, which is what the kappa < 1 loss was. "
            "*** AND THE THROAT TRANSMITS AT THE TRAVERSAL TIME. *** The "
            "measured object is the throat operator's own impulse response, so "
            "no source or observer geometry enters: r_11, same mouth in and "
            f"out, starts at {t3.get('rows', [{}])[0].get('onset_same_mouth', 0):.1f} "
            "-- a wave that reaches a mouth is partly reflected INSTANTLY -- "
            "and r_12, opposite mouths, starts at min(L, d) with "
            f"d(onset)/dL = {t3.get('slope_below_the_ambient_path', 0):.4f} "
            "against a predicted 1, saturating above the ambient path to "
            f"{t3.get('onset_spread_above_it', 0):.1e}. THE AMBIENT ALSO "
            "CONNECTS THE TWO MOUTHS, along a geodesic of length d, whether or "
            "not they are joined: PR #258's cross-mouth channel and PR #259's "
            "beta = 0 control, now SEPARATED IN TIME instead of by rank "
            "counting. A point throat transmits at "
            f"{t3.get('point_throat_onset_opposite', 0):.1f}, which is what a "
            "point throat is. THE LEDGER IS A DERIVATION, NOT A FIT: on the "
            "contour cot x = -i - 2i sum e^{2ikx} and csc x = -2i sum "
            f"e^{{i(2k+1)x}}, to {t4.get('cot_series_error', 0):.1e} and "
            f"{t4.get('csc_series_error', 0):.1e}, so the same-mouth entry "
            "carries 0, 2L, 4L and the cross-mouth entry L, 3L, 5L. The "
            "parities are the physics, and the reflected channel is the one "
            "the rank-one model does not have at all. THERE IS NO POINT LIMIT "
            "-- THERE IS A BAND. Freezing A at A(lam_0) is exact at lam_0 and "
            f"{100 * t5.get('band', [{}])[-1].get('relative_error', 0):.0f}% "
            "wrong at 3 lam_0, and the two channels fail differently: the "
            "antisymmetric one has a limit, -L/(2A) with an O(L^2) error over a "
            "decade in L, and the symmetric one DIVERGES like 2/(A lam L), "
            "because a massless tube holds a zero mode and a point cannot. "
            "THE SAME ZERO MODE BREAKS PR #258's TOMOGRAPHY. At lam = 0 the "
            "static response collapses onto [[1,-1],[-1,1]] to "
            f"{t6.get('worst_antisymmetry', 0):.1e}, det S goes to zero "
            f"LINEARLY in lam with coefficient "
            f"{t6.get('linear_coefficient', 0):.2f}, and W = S_12/det S - G_0 "
            "diverges like 1/lam. A point throat is statically RANK TWO and a "
            "massless finite throat RANK ONE -- a falsifiable difference from a "
            "static measurement. Give the tube an interior mass and the rank "
            "returns; and off the collapse W = -beta(lam) exactly, to "
            f"{t6.get('worst_defect_error', 0):.1e}, with beta the tube's own "
            "transmission amplitude. PR #258's theorem survives the "
            "generalization and returns the interior's amplitude instead of a "
            "constant. AN INTERIOR MASS IS A TRANSMISSION CUTOFF: below "
            f"lam = m^2 the tube is evanescent, beta matches its exponential "
            f"asymptote to {t7.get('asymptote_error', 0):.1e}, is suppressed by "
            f"{t7.get('suppression_ratio', 0):.1e}, and has "
            f"{t7.get('sign_changes_below', 0)} sign changes against "
            f"{t7.get('sign_changes_above', 0)} above -- monotone decay against "
            "oscillation, which is the discriminator rather than the sign "
            "itself. *** AND THE MODEL ALWAYS HAS A GROWING MODE, WHICH IS THE "
            "MOUTH'S AND NOT THE TUBE'S. *** A(lam) decreases and Gamma(lam) "
            "increases -- PR #257's Gram identity -- so A - Gamma is strictly "
            "monotone between poles and each channel has at most one root, a "
            "count rather than a scan. The symmetric channel always has exactly "
            "one, at lam < 0. Three facts identify what it is, all limits with "
            "their convergence measured: its rate matches sigma* = "
            f"2 sqrt(pi/A) to {t8.get('worst_closed_form_error', 0):.1e} with "
            "NO L in it; two mouth separations agree to "
            f"{t8.get('separation_spread_far', 0):.1e}; and the channel "
            "splitting is 1.04 e^{-sigma* d}, the Euclidean propagator between "
            "the mouths. A mode that ignores the tube's length and the mouths' "
            "separation and does not distinguish the channels is a SINGLE-MOUTH "
            "object, and its scale sqrt(A/4pi) is the mouth's own radius -- "
            "exactly where 'point mouth' stops being an approximation. Every "
            "statement here is therefore made at |lam| << sigma*^2, and that "
            f"band ({t8.get('the_working_band', 0):.2f} at the working point) "
            "is quoted rather than assumed. SO THE CONTOUR MUST CLEAR IT: "
            "placed 0.03 below sigma*, the inversion returns a field with "
            "support before its own light cone -- a pedestal at "
            f"{100 * t9.get('pedestal_below', 0):.0f}% of the peak and an onset "
            f"of {t9.get('onset_below', 0):.1f} for an event that cannot begin "
            f"until {t9.get('true_onset', 0):.1f} -- against "
            f"{t9.get('pedestal_above', 0):.1e} when it is placed above. Same "
            "species as PR #259's under-resolved contour, reported the same "
            "way, and this time the rule is checkable in advance because "
            "sigma* has a closed form. WHAT IS STILL PUT IN: the background, "
            "the mouth positions, and L, A, m -- three numbers with dimensions "
            "instead of four without. NO BACKREACTION: the throat is a fixed "
            "background, and PR #261's common action and PR #262's A/B/A+B "
            "metric backreaction are the next two steps. They now have an "
            "object with a proper length, a delay, and a conservation law.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "finite_throat", "tests": tests,
            "verdict_class": verdict_class, "verdict": verdict,
            "working_throat": {"separation": WORKING_THROAT.separation,
                               "length": WORKING_THROAT.length,
                               "area": WORKING_THROAT.area,
                               "interior_mass": WORKING_THROAT.interior_mass},
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
    if isinstance(v, complex):
        return f"{v.real:.6g}{v.imag:+.6g}j"
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
    if isinstance(o, complex):
        return {"real": o.real, "imag": o.imag}
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_finite_throat_probe"
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
