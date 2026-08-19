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

WHERE THE SELF-ADJOINTNESS LIVES. The conservative object is the ENLARGED
system, ambient (+) tube, with the lam-INDEPENDENT matching above: one
self-adjoint operator on L2(S3) (+) L2([0,L]). Eliminating the tube leaves an
ambient-only problem with a lam-DEPENDENT boundary condition A(lam) = -N(lam)^-1
-- the Weyl (M-) function of that elimination. A(lam) is NOT itself a
self-adjoint operator on the ambient space and is not claimed to be; an
energy-dependent boundary condition never is. What it is, is a matrix Nevanlinna
function, monotone in lam between its poles, and that monotonicity is the
enlarged system's self-adjointness showing through. That lam-dependence IS the
interior.

WHAT IS CHECKED
───────────────
T2  THE ELIMINATION IS FAITHFUL, AND THE ENLARGED SYSTEM IS CONSERVATIVE. At
    each lam the eliminated problem is a maximal self-adjoint boundary condition
    for the ambient problem AT THAT lam: rank[B|C] = 2 and BC^dag = CB^dag to
    0.0, at seven values on both sides of zero. The DtN map is checked against
    the interior it summarizes, by the SESQUILINEAR Green's identity under
    quadrature -- the one that expresses energy -- rather than against itself.
    The control is PR #255's rank-one transfer model, whose defect is 0.3: the
    size of the coupling itself, which is what the kappa < 1 loss was.

T3  *** THE THROAT TRANSMITS AT THE TRAVERSAL TIME. *** The measured object is
    the TWO-MOUTH BLOCK's impulse response, R(w) = (A(w) - Gamma(w))^-1 inverted
    along the retarded contour: the source and observer legs are gone, but Gamma
    is the AMBIENT's own mouth-to-mouth propagator and stays in, so this is the
    coupled ambient+tube response and not the throat alone. That is exactly why
    the second prediction reads min(L,d). Two different predictions, both met:

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

T4  AND THE MASSLESS TUBE'S DELAY LEDGER IS A DERIVATION. On the contour
    cot x = -i - 2i sum e^{2ikx} and csc x = -2i sum e^{i(2k+1)x}, verified to
    4.5e-16 and 1.7e-15: the same-mouth entry carries delays 0, 2L, 4L (an
    instantaneous reflection and its echoes) and the cross-mouth entry L, 3L,
    5L. The parities are the physics, and the reflected channel is the one the
    rank-one model does not have at all. TWO SCOPES: this is the m = 0 kernel,
    where k = w makes e^{ikL} a pure translation -- with m /= 0 the front is
    still at L but e^{ikL} is DISPERSIVE and the echoes are not translated
    copies (T7 is the same physics from the other side) -- and it is the ledger
    of the TUBE KERNEL A(w), not of the coupled R = (A - Gamma)^-1, which also
    carries the ambient's d-paths. That is why T3 reads min(L,d).

T5  THERE IS A POINT LIMIT, AND IT IS NOT A FINITE A. Freezing A at A(lam_0) is
    exact at lam_0 and nowhere else -- 4% out at 1.05 lam_0, 121% out at 3
    lam_0 -- everything in A varies through kL, so the range over which
    freezing it is defensible is an O(1/L) FREQUENCY scale -- and as
    L -> 0 the antisymmetric channel converges to -L/(2A) while the symmetric
    one DIVERGES like 2/(A lam L). (An earlier draft also quoted a universal
    Dlam ~ 2 sqrt(lam)/L; that is only the local linearization Dlam = 2 w Dw,
    it drops the (Dw)^2 term that matters at low frequency, and this test does
    not sweep L to extract a fixed-error bandwidth anyway.) A first draft
    concluded from the divergence that the limit does not exist. It does. A boundary pair is defined up to
    (B,C) -> (MB,MC), so a diverging chart matrix means the limit has LEFT THE
    CHART. Row-scaled, the pair converges to the projector pair

        (B, C) -> (P_anti, -P_sym),   i.e.  Phi_anti = 0 and q_sym = 0

    a MIXED DIRICHLET-NEUMANN stratum: maximal (rank[B|C] = 2 throughout),
    self-adjoint, and reached by no finite Hermitian A because both blocks are
    singular. Convergence is linear in L at rate A lam/2 = 2pi, measured. So the
    correct statement is "no finite-A point limit", and it is exactly the kind
    of stratum PR #257's review said the chart does not cover.

T6  THE STATIC LIMIT IS RANK ONE, AND PR #258's TOMOGRAPHY BREAKS ON IT. The
    same zero mode empties the symmetric channel at lam = 0: S = Re R collapses
    onto [[1,-1],[-1,1]] to 4e-5, det S -> 0 LINEARLY in lam with coefficient
    149.08, and the disconnection defect W = S_12/det S - G_0 diverges like
    1/lam. WHAT THAT FALSIFIES is the generic finite-A family, every member of
    which is rank two -- and NOT point-ness: the tube's own short-tube stratum
    from T5 gives R -> diag(0, -1/Gamma_anti), rank one as well, and the tube
    converges to it. Both are measured here, because the first draft claimed the
    stronger and wrong version.
    Give the tube an interior mass and the rank returns with det S ~ -148.7 m^2,
    and then the closure: off the collapse W = -beta(lam) EXACTLY, to 3.1e-13,
    with beta the tube's own transmission amplitude. PR #258's theorem survives
    the generalization; what it returns is no longer a constant.

T7  AN INTERIOR MASS GIVES THE CHANNEL A MASS GAP. Below lam = m^2 the channel
    is EVANESCENT: beta -> -csch(kappa L)/(A kappa), matched by its exponential
    asymptote to 7.6e-05, monotone, and with ZERO sign changes against 7 above
    the cutoff. Suppression 3.1e-03. (An earlier draft called this "low-pass",
    which is backwards -- it is the LOW frequencies that are suppressed.) The
    cutoff is also where the rank collapses: one statement read at lam = m^2.

T8  *** THE MODEL FAILS THE STABILITY GATE, AND THE FAILURE IS THE MOUTH'S. ***
    A(lam)
    decreases and Gamma(lam) increases (PR #257's Gram identity), so A - Gamma
    is strictly monotone between poles and each channel has AT MOST ONE root --
    a count, not a scan. The symmetric channel always has exactly one, at
    lam < 0. Three facts identify it, all as limits with their convergence
    measured: its rate matches sigma* = 2 sqrt(pi/A) to 1.5e-03 with NO L in it;
    two mouth separations agree to 3.9e-09; and the channel splitting is
    1.04 e^{-sigma* d}, the Euclidean propagator between the mouths, so at that
    scale the mouths do not know about each other. Its length scale is
    sqrt(A/4pi). So the mode is generated at the POINT-MOUTH/TUBE INTERFACE
    and, in the sigma L, sigma d >> 1 limit, LOCALIZES to a single mouth and
    stops knowing L and the other mouth. THE WORKING THROAT IS NOT IN THAT
    LIMIT, and saying so matters: at A = 4 pi the asymptotic form gives
    sigma* = 1 while L = 0.9 gives 1.417, and sigma* runs 1.769 ... 1.152 across
    L = 0.4 ... 3 -- a spread of 0.617, 54% of the smallest value. An earlier
    draft claimed the mode "belongs to the mouth and not the interior", which is
    stronger than the data support. THIS IS THE ROUND'S FALSIFICATION RESULT and nothing
    here cures it: whether a finite-radius mouth or neck geometry removes the
    mode is open, and is the thing to settle before stationary-action or
    backreaction work.

T9  SO THE CONTOUR MUST CLEAR IT -- AND BE RESOLVED, WHICH ARE TWO CONDITIONS.
    Placed 0.03 BELOW sigma*, the inversion returns a field with support before
    its own light cone: pedestal 99% of the peak, onset 0.0 for an event that
    cannot begin until t = 0.6. But an earlier draft stated the rule as
    "eps > sigma*" and its own table contradicted it: at clearance +0.02 the
    contour IS above the mode and the pedestal is still 2.6e-03 with the onset
    at 0. That clearance is 0.95 of the frequency spacing 2 pi/span = 0.0209 --
    the pole is above the contour but UNRESOLVED BY THE GRID, PR #259's lesson
    arriving a second time. The rule has two parts:

        eps > sigma*              analytic:  the contour clears the pole
        eps - sigma* >> 2pi/span  numerical: the grid resolves the clearance

    Both have closed forms, so both are checkable before the solve. At
    clearances of 14.3, 38.2 and 71.6 spacings the pedestal is 1.0e-16 and the
    recovered onset converges -- 0.4646, 0.4601, 0.4555, a spread of 0.0092,
    four time steps. What clearing the contour does NOT do is stabilize
    anything: above sigma* the inversion returns the correct retarded solution
    OF AN UNSTABLE SYSTEM. The delay is read from the causal ONSET, which is
    immune to what the solution does afterwards.

WHICH FREQUENCIES EACH RESULT USES
─────────────────────────────────
The band is not uniform across the round, so: T3 and T4 are statements about the
ANALYTIC STRUCTURE OF THE EXACT MODEL AT ALL FREQUENCIES. A causal onset is a UV
object, and the probe pulse that resolves it (width 0.03) carries content out to
w ~ 30, far above sigma* ~ 1.4. They are exact results ABOUT THIS MODEL, not
predictions about a resolved physical mouth. T2, T5, T6 and T7 are low-frequency
statements and sit inside the band. Relatedly, A is a ONE-DIMENSIONAL COUPLING
STRENGTH and not an area with a mouth radius attached -- reading sqrt(A/4pi) as a
radius at the working point would give one of order the whole unit S3, which is
another way of saying the same thing.

WHAT IS NOT CLOSED
──────────────────
The point-mouth matching is UNSTABLE (T8), and that is the round's closure
result. It GATES the roadmap: an action or a backreaction computed on a
background with a growing mode inherits the mode, so PR #261's common action and
PR #262's A/B/A+B metric backreaction both wait. The next construction is a
FINITE-RADIUS MOUTH OR NECK -- the ambient solved outside two small balls rather
than a point interaction with a radius parameter -- and it has one question to
answer: DOES THE NEGATIVE MODE SURVIVE? If it does, the point-interaction throat
family of PRs #255-#260 is the wrong model of a wormhole mouth. The mouths are points, the interior is
one-dimensional so A is a coupling and not a geometry, and L, A, m are chosen.
What this round does hand the next two is an object with a proper length, a
delay, a conservation law -- and a stated failure mode.

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
    measure_the_contour_must_clear_the_growing_mode,
    measure_the_enlarged_system_is_conservative,
    measure_the_delay_ledger_is_the_bounce_series,
    measure_the_growing_mode_is_interface_localized,
    measure_the_interior_mass_is_a_transmission_cutoff,
    measure_the_short_tube_limit_is_a_mixed_stratum,
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


def t2_the_enlarged_system_is_conservative() -> dict:
    r = measure_the_enlarged_system_is_conservative()
    return {"name": "T2_the_enlarged_system_is_conservative", **r,
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


def t5_the_short_tube_limit_is_a_mixed_stratum() -> dict:
    r = measure_the_short_tube_limit_is_a_mixed_stratum()
    return {"name": "T5_the_short_tube_limit_is_a_mixed_stratum", **r,
            "pass": bool(r["the_limit_exists_and_is_not_a_finite_A"]
                         and r["the_band_error_reaches_one"]
                         and r["convergence_is_linear_in_L"]
                         and r["every_pair_is_maximal"])}


def t6_the_static_limit_is_rank_one() -> dict:
    r = measure_the_static_limit_is_rank_one_and_the_defect_diverges()
    return {"name": "T6_the_static_limit_is_rank_one", **r,
            "pass": bool(r["the_static_response_is_rank_one"]
                         and r["it_falsifies_the_finite_A_family_not_point_ness"]
                         and r["det_S_is_linear_in_lambda"]
                         and r["the_defect_diverges"]
                         and r["the_defect_is_still_minus_beta"])}


def t7_the_interior_mass_is_a_cutoff() -> dict:
    r = measure_the_interior_mass_is_a_transmission_cutoff()
    return {"name": "T7_the_interior_mass_is_a_cutoff", **r,
            "pass": bool(r["the_evanescent_side_does_not_oscillate"]
                         and r["the_transmission_is_suppressed_below_cutoff"]
                         and r["everything_stays_real"])}


def t8_the_growing_mode_is_interface_localized() -> dict:
    r = measure_the_growing_mode_is_interface_localized()
    return {"name": "T8_the_growing_mode_is_interface_localized", **r,
            "pass": bool(r["every_throat_has_one"]
                         and r["it_stops_knowing_the_separation"]
                         and r["the_closed_form_holds_once_sigma_L_is_large"]
                         and r["the_split_is_the_euclidean_propagator"])}


def t9_the_contour_must_clear_it() -> dict:
    r = measure_the_contour_must_clear_the_growing_mode()
    return {"name": "T9_the_contour_must_clear_it", **r,
            "pass": bool(r["a_contour_below_the_mode_breaks_causality"]
                         and r["clearing_the_mode_is_not_enough"]
                         and r["the_resolved_contours_are_clean"]
                         and r["the_recovered_onset_converges"])}


def t10_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T10_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_enlarged_system_is_conservative(),
             t3_the_throat_transmits_at_the_traversal_time(),
             t4_the_delay_ledger_is_the_bounce_series(),
             t5_the_short_tube_limit_is_a_mixed_stratum(),
             t6_the_static_limit_is_rank_one(),
             t7_the_interior_mass_is_a_cutoff(),
             t8_the_growing_mode_is_interface_localized(),
             t9_the_contour_must_clear_it()]
    tests.append(t10_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8, t9 = tests[1:9]

    if all(t["pass"] for t in tests):
        verdict_class = ("THE_INTERIOR_GIVES_A_DELAY_AND_THE_POINT_MOUTH_IS_"
                         "UNSTABLE")
        verdict = (
            "THE THROAT NOW HAS AN INTERIOR, AND THE INTERIOR IS THE DELAY. "
            "PRs #253-#259 all carried the same disclaimer -- point-supported, "
            "no proper length, no delay, a rank-one mouth-transfer model that "
            "is lossy for kappa < 1 -- and this round replaces it with a tube "
            "of length L and cross-section A whose Dirichlet-to-Neumann map is "
            "exact. Two lines are put in: N(lam) = A k [[cot kL, -csc kL], "
            "[-csc kL, cot kL]] and the matching q = -N Phi. Everything else "
            "follows. THE CONSERVATIVE OBJECT IS THE ENLARGED SYSTEM, ambient "
            "(+) tube, with lam-independent matching; eliminating the tube "
            "leaves a lam-DEPENDENT boundary condition -- the Weyl function of "
            "that elimination, which is a Nevanlinna function and not itself a "
            "self-adjoint operator on the ambient space. What is checked is "
            "that the elimination is faithful at each lam: rank[B|C] = 2 and "
            f"BC^dag - CB^dag = "
            f"{t2.get('worst_hermiticity_defect', 0):.1e} at seven values of "
            "lam on both sides of zero, with the DtN map checked against the "
            "interior it summarizes by the SESQUILINEAR Green's identity to "
            f"{t2.get('worst_green_identity_residual', 0):.1e} rather than "
            "against itself; PR #255's rank-one transfer model, the control, "
            f"has defect {t2.get('the_rank_one_control_defect', 0):.2f} -- the "
            "size of the coupling itself, which is what the kappa < 1 loss was. "
            "*** AND THE THROAT TRANSMITS AT THE TRAVERSAL TIME. *** The "
            "measured object is the TWO-MOUTH BLOCK's impulse response -- the "
            "source and observer legs are gone, but Gamma, the ambient's own "
            "mouth-to-mouth propagator, stays in, so it is the coupled "
            "ambient+tube response and not the throat alone: r_11, same mouth in and "
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
            "the rank-one model does not have at all. THERE IS A POINT LIMIT, "
            "AND IT IS NOT A FINITE A. Freezing A at A(lam_0) is exact at "
            f"lam_0 and "
            f"{100 * t5.get('band', [{}])[-1].get('relative_error', 0):.0f}% "
            "wrong at 3 lam_0 -- a band of width ~1/L in omega -- and as L -> 0 "
            "the antisymmetric channel converges to -L/(2A) while the symmetric "
            "one diverges like 2/(A lam L). A first draft concluded that the "
            "limit does not exist. It does: a boundary pair is defined up to "
            "(B,C) -> (MB,MC), so a diverging chart matrix means the limit has "
            "LEFT THE CHART, and row-scaled the pair converges to (P_anti, "
            "-P_sym) -- Phi_anti = 0, q_sym = 0, a MIXED DIRICHLET-NEUMANN "
            f"stratum, maximal throughout and "
            f"{t5.get('distance_to_the_stratum', 0):.2f} away at L = 0.02, "
            "linearly in L. No finite Hermitian A reaches it, which is exactly "
            "the stratum PR #257's review said the chart does not cover. "
            "THE SAME ZERO MODE BREAKS PR #258's TOMOGRAPHY. At lam = 0 the "
            "static response collapses onto [[1,-1],[-1,1]] to "
            f"{t6.get('worst_antisymmetry', 0):.1e}, det S goes to zero "
            f"LINEARLY in lam with coefficient "
            f"{t6.get('linear_coefficient', 0):.2f}, and W = S_12/det S - G_0 "
            "diverges like 1/lam. WHAT THAT FALSIFIES is the generic finite-A "
            "family, every member of which is rank two -- and NOT point-ness: "
            "the tube's own short-tube stratum is rank one as well, and the "
            "tube converges to it. The first draft claimed the stronger and "
            "wrong version. Give the tube an interior mass and the rank "
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
            "itself; the interior has a MASS GAP and the channel below it is "
            "evanescent, which is the opposite of a low-pass filter. *** AND "
            "THE MODEL FAILS THE STABILITY GATE, WITH A GROWING MODE THAT IS "
            "THE MOUTH'S AND NOT THE TUBE'S. *** A(lam) decreases and Gamma(lam) "
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
            "object IN THAT LIMIT: the mode is generated at the point-mouth/"
            "tube INTERFACE and localizes to one mouth only for sigma L, "
            "sigma d >> 1. The working throat is NOT there -- sigma* runs "
            "1.769...1.152 across L = 0.4...3 at fixed A, a spread of "
            f"{t8.get('length_spread_at_the_working_area', 0):.3f} -- so the "
            "honest statement is the interface one, not 'the mouth's alone'. "
            "THIS IS THE ROUND'S FALSIFICATION RESULT, and "
            "nothing here cures it -- whether a finite-radius mouth or neck "
            "geometry removes the mode is open, and is the thing to settle "
            "before stationary-action or backreaction work. SO THE CONTOUR "
            "MUST CLEAR IT: "
            "placed 0.03 below sigma*, the inversion returns a field with "
            "support before its own light cone -- a pedestal at "
            f"{100 * t9.get('pedestal_below', 0):.0f}% of the peak and an onset "
            f"of {t9.get('onset_below', 0):.1f} for an event that cannot begin "
            f"until {t9.get('true_onset', 0):.1f} -- against "
            f"{t9.get('pedestal_resolved', 0):.1e} when it is placed above AND "
            "RESOLVED. Same species as PR #259's under-resolved contour, and "
            "the rule turns out to have TWO parts, which an earlier draft got "
            "wrong: eps > sigma* is the analytic Bromwich condition, but at a "
            "clearance of 0.95 frequency spacings the pedestal is still "
            "2.6e-03 because the grid does not RESOLVE the clearance. Both "
            "eps > sigma* and eps - sigma* >> 2 pi/span are needed, both have "
            "closed forms, and at 14-72 spacings the onset converges to "
            f"{t9.get('onset_spread_across_resolved', 0):.4f}. Clearing the "
            "contour STABILIZES "
            "NOTHING: above sigma* the inversion returns the correct retarded "
            "solution OF AN UNSTABLE SYSTEM, and the delay is read from the "
            "causal onset, which is immune to what happens afterwards. WHICH "
            "FREQUENCIES: the delay and the ledger are statements about the "
            "exact model's analytic structure at ALL frequencies -- a causal "
            "onset is a UV object and the probe pulse carries content to "
            "w ~ 30, far above sigma* -- so they are exact results about this "
            "model rather than predictions about a resolved physical mouth; "
            "the static and low-frequency results sit inside the band. A is a "
            "ONE-DIMENSIONAL COUPLING, not an area with a radius attached. "
            "WHAT IS STILL PUT IN: the background, the mouth positions, and "
            "L, A, m -- three numbers with dimensions where the real-field "
            "point sector has three without them. NO BACKREACTION: the throat "
            "is a fixed background. AND THE INSTABILITY GATES WHAT COMES "
            "NEXT: an action or a backreaction computed on a background with a "
            "growing mode inherits the mode, so PR #261's common action and PR "
            "#262's A/B/A+B metric backreaction both wait on a FINITE-RADIUS "
            "MOUTH OR NECK -- the ambient solved outside two small balls "
            "rather than a point interaction with a radius parameter -- whose "
            "single question is whether the negative mode survives. If it "
            "does, the point-interaction throat family of PRs #255-#260 is the "
            "wrong model of a wormhole mouth. If it does not, the delay and "
            "the conservation law carry over unchanged.")
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
