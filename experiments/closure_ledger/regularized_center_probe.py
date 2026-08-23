#!/usr/bin/env python3
"""
Regularized-center probe — the middle as a finite bearing, not a point
======================================================================

THE PROPOSAL. Every picture in this arc has put a POINT in the middle, because
a point is where the clock-hands story works: two radial arms P_A -> O -> P_B
change direction at O for free, since at a point there is no angular direction
left to change. That is the property the connection is wanted for -- it makes
the link's cost independent of where the mouths sit -- and it is bought with
f = 0, where the geometry stops existing.

Regularise it. Keep dl^2 = ds^2 + f(s)^2 dOmega^2 and replace f(0) = 0 with

    f_min = f0 > 0 .

In the 2-D cross-section the middle becomes a small CIRCLE; in a d-dimensional
bulk it is the space of radial directions, S^(d-1), or RP^(d-1) if the clock
hand is an unoriented axis (n ~ -n). The scalar-flat profile f'^2 = 1 - f0/f is
used as the worked example, because it is the one PR #265 FORCED and so every
closed form below can be checked against physical_throat rather than asserted.

T1  THE ARMS ARE THE REPO'S OWN GEOMETRY, WITH THE SYMMETRY DROPPED. The proper
    distance from the neck to an end of transverse scale F is exactly

        L(F) = sqrt(F(F - f0)) + f0 arcosh sqrt(F/f0) ,

    and one arm's share of the resistance is exactly

        I(F) = (2/f0) sqrt(1 - f0/F) .

    Both check against quadrature, and -- the point of the check -- setting
    f_o = f_i = sin a with f0 = sin^3 a reproduces VacuumThroat.length() and
    .resistance() BIT FOR BIT. The regularised centre is not a new geometry; it
    is the forced throat with the two arms allowed to differ.

T2  *** THE TWO ARMS ARE INDEPENDENT -- INTRINSICALLY. *** L_o != L_i, and
    nothing pathological happens: scanned to L_o/L_i = 437 the resistance, the
    turn cost and the small-angle law all keep their character. Unequal arms
    are honest truncations of ONE scalar-flat profile.
    THE CAVEAT THAT WAS MISSING: they are NOT a fully matched #265 throat with
    the symmetry merely removed. C1 matching to a UNIT ROUND S^3 gives
    F_j = sin a_j and |f'_j| = cos a_j, which with f'^2 = 1 - f0/F_j forces
    f0 = sin^3 a_j = F_j^3 at each end; one common f0 then fixes
    F = f0^(1/3) UNIQUELY, so BOTH MATCHED MOUTHS HAVE THE SAME SCALE.
    Asymmetric arms need asymmetric ambient attachments -- different curvature
    scales, a shell or junction, nonzero K -- and the interior says nothing
    about whether those exist.
    ALSO CORRECTED: R_inner/R_outer does NOT "stop carrying physical
    significance". T3 immediately derives w_i/w_o = f_i/f_o, so the ENDPOINT
    SCALE RATIO is load-bearing. What loses significance is the old vacuole's
    ARBITRARY DRAWING ratio, which was never identified with the endpoint
    scales at all.

T3  SCALE TRANSPORT IS EXPLICIT. A feature of angular width dtheta has physical
    width w(s) = f(s) dtheta, so along the route w_o -> w_min -> w_i, with
    w_i/w_o = f_i/f_o and f0 nowhere in it. A packet does NOT teleport from one
    radius to another at fixed size. (The ANGULAR width is what is transported
    and is constant; the physical width is not.)

T4  *** THE TURN COST IS QUADRATIC, NOT LINEAR. *** This is the round's result,
    and it is a CORRECTION to the proposal it came from.
    The natural guess is that turning through alpha around a bearing of radius
    f0 costs its arc, f0*alpha. That is EXACTLY RIGHT for one route -- down,
    turn on the bearing, up -- and it is an honest upper bound. But it is not
    the geodesic. WHY THE CORNER LOSES IS PYTHAGORAS, NOT LEVERAGE: the corner
    pays its angle as PURE TRANSVERSE MOTION at the neck (f0*alpha, first
    order), while the geodesic TILTS motion that is already radial, and a tilt
    costs sqrt(ds^2 + f^2 dtheta^2) - ds ~ f^2 theta'^2 ds/2 -- SECOND order.
    Integrating Clairaut's f sin(psi) = h gives

        T(alpha) = alpha^2/(2 I2) + O(alpha^4) ,   I2 = int ds/f^2 ,

    LEADING ORDER, not an exact law -- the exact object is the integrated
    turn_cost, and T is 0.9248 of the quadratic value at pi. At q = 2 the same
    I2 is the throat's monopole resistance: T = alpha^2 (4 pi/I2)/(8 pi).
    For two long arms I2 -> 4/f0 and the law reads T -> f0 alpha^2/8.
    Checked against the integrated geodesic, not against the expansion: exact
    to 8e-05 at alpha = 0.1.
    (An earlier draft explained this as turning "where the lever arm f is
    longer, so a given angle costs less arc". BACKWARDS on both counts: an
    angular increment at larger f costs MORE arc, f dtheta, and theta' = h/f^2
    puts the turn where f is SMALLEST -- 76% of it inside f < 2.4 f0.)

T5  SO THE LINEAR GUESS IS AN UPPER BOUND, AND A LOOSE ONE. The geodesic spends
    1.25% of f0*alpha at alpha = 0.1 and 36% at alpha = pi -- and that fraction
    is a function of ALPHA ALONE, the same to three figures across two decades
    of f0, which is itself the signature of the quadratic law.

T6  *** THE HINGE IS NEVER THE EXPENSIVE PART. *** Which is the whole point:
    the property the point centre was wanted for SURVIVES regularisation, and
    survives it more strongly than proposed. At the working point a full
    half-turn costs 8e-04 of the arms -- and pi is the LARGEST SEPARATION
    THERE IS, since bearing_distance reduces any pair of directions to [0,pi].
    So that is the worst case over the whole configuration space, and THERE IS
    NO REACHABLE ORIENTATION at which the hinge costs as much as the journey.
    (An earlier draft reported a "break-even angle" of 104 rad from
    sqrt(2 I2 L). WITHDRAWN: it extrapolates the small-angle law outside its
    own domain AND outside the configuration space. Only the fact that it
    exceeds pi carries content.)

T7  AND THE POINT MODEL IS THE LIMIT, NOT A RIVAL. T is linear in f0 at fixed
    alpha, and the arms' own excess over their naive value obeys
    L(F) - F ~ (f0/2)[ln(4F/f0) - 1]. Both vanish with the bearing, so f0 -> 0
    returns L_o + L_i, independent of where the mouths are.

T8  *** WHAT "INTERSECTION" BECOMES. *** Two fronts no longer collide at r = 0.
    They land on the bearing at angular positions, and the question splits in
    two halves that behave differently: WHETHER they meet is a statement about
    angles, with f0 nowhere in it; HOW BIG the meeting is depends on f0 AND on
    the bearing's dimension -- a LENGTH f0 x (overlap angle) on the drawn S^1
    cross-section, an AREA ~f0^2 x (angular area) on S^2, a VOLUME ~f0^3 on
    S^3. The yes/no criterion is dimension-free; the size law is not.
    THE POINT LIMIT IS NOT WHAT IT LOOKS LIKE: f0 -> 0 does NOT make
    every route meet. It shrinks the overlap AND the gap to zero together, so
    the distinction survives as a yes/no and disappears as a length.

T9  THE DRAWN CIRCLE IS HONEST. Two directions span a totally geodesic 2-plane,
    so the walk between them is a great-circle arc whatever the bearing's
    dimension is -- S^1, S^2, S^3 give the same T. Checked by computing the
    cost from an S^2 dot product and again from the planar angle in the frame
    the pair spans. The one thing that DOES change the answer is the projective
    identification n ~ -n, which replaces alpha by min(alpha, pi - alpha) and
    is worth a factor of 418 near alpha = pi.

T10 WHAT IS GENERAL, AND WHAT IS ONLY ABOUT THIS NECK. The derivation of
    T = alpha^2/(2I) never used the profile: minimising int f^2 theta'^2/2 ds
    at fixed int theta' ds = alpha gives theta' = h/f^2, so alpha = h I and the
    excess is h^2 I / 2 = alpha^2/(2I) for ANY rotationally symmetric f(s) with
    a finite minimum, with its own I. Checked on a second, unrelated profile --
    f = sqrt(f0^2 + s^2), not scalar-flat and not glued to anything -- where it
    holds to 1.3e-04 across six decades of f0.
    What does NOT carry over is the correction beyond leading order: the shape
    T/(alpha^2/2I) at alpha = pi is 0.9250 on the scalar-flat profile and
    0.8886 on the hyperbolic one. SO THE QUADRATIC LAW IS A STATEMENT ABOUT
    NECKS; THE ~8-11% DEFICIT AT A HALF TURN IS A STATEMENT ABOUT A PARTICULAR
    ONE.

T11 *** THE MOMENT HIERARCHY: I2 IS UNIVERSAL, I4 IS LOCAL. *** Expanding the
    geodesic exactly, dl/ds - 1 = (1 - h^2/f^2)^(-1/2) - 1, gives
    T = h^2 I2/2 + 3 h^4 I4/8 and alpha = h I2 + h^3 I4/2 with I_n = int ds/f^n.
    Eliminating h:

        T(alpha) = alpha^2/(2 I2)  -  alpha^4 I4/(8 I2^4)  +  O(alpha^6) ,
        shape = T/(alpha^2/2 I2) = 1 - alpha^2 I4/(4 I2^3) .

    THAT IS THE DIVISION OF LABOUR. I2 is the only moment in the leading term,
    and it is the SAME I2 that sets the monopole conductance -- so the
    quadratic hinge is universal in the strong sense that a quantity the
    throat already had fixes it. I4 is the first moment shared with nothing,
    and it is where the neck's SHAPE is first felt:
        scalar-flat  I2 = 4/f0,   I4 = 32/(15 f0^3)  ->  shape = 1 - a^2/120
        hyperbolic   I2 = pi/f0,  I4 = pi/(2 f0^3)   ->  shape = 1 - a^2/(8 pi^2)
    1/120 against 1/79 sets HOW FAST they separate: the shape difference is
    alpha^2 (1/79 - 1/120) = alpha^2 x 4.33e-03, so 4.3e-05 at alpha = 0.1
    growing to 3.5e-02 at pi. (An earlier draft said the two "agree to eight
    digits at alpha = 0.1". WRONG AND WITHDRAWN: 8 digits is how well EACH
    profile matches ITS OWN quartic law, 7.9e-09; the two match EACH OTHER
    only to 4.3e-05. Two quantities, and the smaller was reported for the
    larger.) NOTE ALSO that I4 is not the ENTIRE profile dependence -- I6 and
    beyond enter at O(alpha^6), and at pi the quartic prediction misses by
    7.2e-03. And I2 is not itself universal (4/f0 here, pi/f0 there); what is
    universal is the LEADING FUNCTIONAL FORM alpha^2/(2 I2), with I4 the first
    ADDITIONAL INDEPENDENT MOMENT to enter.
    The moments are closed forms checked against quadrature; the shapes are
    checked against the integrated geodesic.

T12 *** THE IDENTITY UNDERNEATH: ONE DIRICHLET FORM. *** Static monopole flux
    and infinitesimal rotation of the throat are not two problems sharing a
    number. They are ONE variational problem on the tube,

        minimise  E[phi] = int f^2 phi'^2 ds  at fixed increment ,

    with Euler-Lagrange (w phi')' = 0: the current w phi' is conserved, the
    increment is c int ds/w, the minimum is (increment)^2 / int ds/w.
    THE TWO WEIGHTS ARE NOT AUTOMATICALLY THE SAME, and an earlier draft said
    they were. The azimuth's weight is the METRIC coefficient f^2, for every
    bearing dimension. The monopole's is the VOLUME element f^q on an S^q
    cross-section, so its resistance is int ds/f^q. THEY COINCIDE EXACTLY AT
    q = 2 -- which is the case in play, a 3-D spatial throat with S^2
    cross-sections, so the identity holds where it is used; but it is a fact
    about that dimension, not about the construction, and unlike the
    great-circle reduction of the hinge it is NOT dimension-free.
        phi = u     : the current IS the flux 4 pi f^2 u'; conductance 4 pi/I2.
        phi = theta : the current IS Clairaut's h = f^2 theta'; cost a^2/(2 I2).
    The sharpest form is not about the two numbers but the two PROFILES:
    normalised to their own total, the monopole potential and the geodesic
    azimuth are THE SAME FUNCTION of position along the tube. The deviation
    falls as alpha^2 -- 4.9e-03 at alpha = 1 down to 4.9e-07 at alpha = 0.01,
    a clean factor of 100 per decade. Which is why "INFINITESIMAL" is the right
    word: at finite alpha the geodesic feels I4 and the potential does not.

WHAT IS PUT IN, AND WHAT IS NOT CLAIMED
───────────────────────────────────────
This is GEOMETRY: a metric, its geodesics, and the transport of an angular
width along them. NO FIELD EQUATION is solved on it. Nothing here says which of
the three candidate readings of the bulk -- finite bearing, finite caustic ring,
finite wormhole neck with moving attachments -- is the right one; it works out
the first far enough to be compared with the others. The scalar-flat profile is
a worked example chosen because it is checkable, not a claim that the bearing
must be that profile. And the arms are static: nothing here evolves.

    python -m experiments.closure_ledger.regularized_center_probe
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.regularized_center import (
    WORKING_CENTER,
    measure_the_arm_length_is_the_repos_own_formula,
    measure_the_bearing_replaces_collision_with_overlap,
    measure_the_bearing_shrinks_back_to_the_point,
    measure_the_corner_route_is_only_an_upper_bound,
    measure_the_fourth_moment_is_where_the_neck_shape_enters,
    measure_the_hinge_and_the_monopole_are_one_dirichlet_form,
    measure_the_law_does_not_depend_on_the_profile,
    measure_where_the_turning_happens,
    measure_the_hinge_is_never_the_expensive_part,
    measure_the_turn_cost_does_not_care_which_sphere,
    measure_the_turn_cost_is_quadratic_not_linear,
    measure_the_two_arms_are_independent,
    measure_the_width_rescales_with_the_profile,
)


def run_probe() -> dict:
    checks: List[dict] = []

    arms = measure_the_arm_length_is_the_repos_own_formula()
    checks.append({
        "id": "T1", "name": "the arms are the repo's geometry, symmetry dropped",
        "detail": arms,
        "pass": bool(arms["the_closed_forms_are_the_integrals"]
                     and arms["it_contains_the_symmetric_throat_exactly"])})

    indep = measure_the_two_arms_are_independent()
    checks.append({
        "id": "T2", "name": "*** the arms are independent INTRINSICALLY (attachment caveat) ***",
        "detail": indep,
        "pass": bool(indep["the_arms_can_be_very_unequal"]
                     and indep["and_nothing_about_the_hinge_changes"]
                     and indep["so_asymmetric_arms_need_asymmetric_attachments"])})

    width = measure_the_width_rescales_with_the_profile()
    checks.append({
        "id": "T3", "name": "scale transport is explicit",
        "detail": width,
        "pass": bool(width["the_ratio_is_the_scale_ratio_to_the_last_ulp"]
                     and width["so_a_packet_does_not_teleport_unchanged"])})

    quad = measure_the_turn_cost_is_quadratic_not_linear()
    checks.append({
        "id": "T4", "name": "*** the turn cost is quadratic, not linear ***",
        "detail": quad,
        "pass": bool(quad["the_small_angle_law_holds"]
                     and quad["the_prefactor_is_one_eighth"]
                     and quad["the_resistance_is_the_repos"])})

    bound = measure_the_corner_route_is_only_an_upper_bound()
    checks.append({
        "id": "T5", "name": "the linear guess is an upper bound, and loose",
        "detail": bound,
        "pass": bool(bound["and_is_an_upper_bound_on_the_geodesic"]
                     and bound["the_fraction_is_a_function_of_alpha_alone"])})

    cheap = measure_the_hinge_is_never_the_expensive_part()
    checks.append({
        "id": "T6", "name": "*** the hinge is never the expensive part ***",
        "detail": cheap,
        "pass": bool(cheap["the_hinge_is_always_cheap"]
                     and cheap["no_reachable_orientation_breaks_even"])})

    limit = measure_the_bearing_shrinks_back_to_the_point()
    checks.append({
        "id": "T7", "name": "the point model is the limit, not a rival",
        "detail": limit,
        "pass": bool(limit["turn_cost_is_linear_in_the_neck"]
                     and limit["the_arm_asymptotic_holds"])})

    meet = measure_the_bearing_replaces_collision_with_overlap()
    checks.append({
        "id": "T8", "name": "*** intersection becomes overlap on the bearing ***",
        "detail": meet,
        "pass": bool(meet["the_verdict_does_not_depend_on_the_neck"]
                     and meet["no_route_passes_through_a_singular_point"])})

    sphere = measure_the_turn_cost_does_not_care_which_sphere()
    checks.append({
        "id": "T9", "name": "the drawn circle is honest; the identification bites",
        "detail": sphere,
        "pass": bool(sphere["the_sphere_and_the_drawn_circle_agree"]
                     and sphere["the_projective_identification_matters"])})

    general = measure_the_law_does_not_depend_on_the_profile()
    checks.append({
        "id": "T10",
        "name": "the law is about necks; the large-angle deficit is not",
        "detail": general,
        "pass": bool(general["the_law_holds_on_the_second_profile"]
                     and general["the_correction_is_profile_dependent"])})

    moments = measure_the_fourth_moment_is_where_the_neck_shape_enters()
    checks.append({
        "id": "T11",
        "name": "*** the leading form is universal; I4 is the first extra moment ***",
        "detail": moments,
        "pass": bool(moments["the_shape_law_holds_at_small_angle"]
                     and moments["the_fourth_order_beats_the_second"]
                     and moments["the_second_moment_is_shared_the_fourth_is_not"]
                     and moments["the_two_profiles_do_not_agree_to_eight_digits"]
                     and moments["the_separation_grows_as_alpha_squared"])})

    identity = measure_the_hinge_and_the_monopole_are_one_dirichlet_form()
    checks.append({
        "id": "T12",
        "name": "*** one Dirichlet form -- and it is a q = 2 statement ***",
        "detail": identity,
        "pass": bool(identity["the_profiles_coincide_as_alpha_vanishes"]
                     and identity["the_deviation_is_second_order_in_alpha"]
                     and identity["the_two_readings_agree"]
                     and identity["the_identity_is_a_q_equals_two_statement"])})

    turning = measure_where_the_turning_happens()
    checks.append({
        "id": "T13",
        "name": "*** the turn sits AT the neck -- a corrected explanation ***",
        "detail": turning,
        "pass": bool(turning["the_turn_is_concentrated_at_the_neck"]
                     and turning["the_geodesic_is_cheaper"])})

    c = WORKING_CENTER
    return {
        "probe": "regularized_center",
        "question": "can the centre keep its special role -- the clock-hand "
                    "hinge whose cost does not care where the mouths are -- "
                    "without being a point?",
        "answer": "yes, and the hinge is cheaper than the proposal assumed: "
                  "the geodesic turn cost is quadratic, T(alpha) = "
                  "alpha^2/(2I), not the bearing's arc f0*alpha",
        "bearing": {"neck": c.neck, "outer_scale": c.outer,
                    "inner_scale": c.inner,
                    "outer_length": c.outer_length(),
                    "inner_length": c.inner_length(),
                    "arm_ratio": c.outer_length() / c.inner_length(),
                    "resistance": c.resistance(),
                    "conductance": c.conductance()},
        "turn_cost_law": "T(alpha) = alpha^2 / (2I) = alpha^2 * "
                         "conductance / (8 pi);  I -> 4/f0 gives f0 alpha^2/8",
        "at_a_half_turn": {
            "alpha": math.pi,
            "geodesic_turn_cost": c.turn_cost(math.pi),
            "linear_guess": c.neck * math.pi,
            "arm_length_sum": c.arm_length_sum(),
            "turn_over_arms": c.turn_cost(math.pi) / c.arm_length_sum()},
        "checks": checks,
        "passed": sum(1 for x in checks if x["pass"]),
        "total": len(checks),
    }


def render_markdown(summary: dict) -> str:
    b = summary["bearing"]
    h = summary["at_a_half_turn"]
    lines = [
        "# Regularized-center probe — the middle as a finite bearing",
        "",
        f"**Question.** {summary['question']}",
        "",
        f"**Answer.** {summary['answer']}",
        "",
        "| the bearing | |",
        "|--|--|",
        f"| neck `f₀` | `{b['neck']:.3e}` |",
        f"| outer scale / arm | `{b['outer_scale']:.4g}` / `{b['outer_length']:.6f}` |",
        f"| inner scale / arm | `{b['inner_scale']:.4g}` / `{b['inner_length']:.6f}` |",
        f"| arm ratio `L_o/L_i` | `{b['arm_ratio']:.4f}` — **not one, and need not be** |",
        f"| resistance `I` | `{b['resistance']:.6e}` |",
        f"| conductance `4π/I` | `{b['conductance']:.6e}` |",
        "",
        f"**The law.** `{summary['turn_cost_law']}`",
        "",
        f"At a half turn (`α = π`): geodesic hinge cost "
        f"`{h['geodesic_turn_cost']:.4e}` against the linear guess "
        f"`{h['linear_guess']:.4e}` — and `{h['turn_over_arms']:.2e}` of the "
        f"arms' own `{h['arm_length_sum']:.6f}`.",
        "",
        f"**{summary['passed']}/{summary['total']} checks pass.**",
        "",
        "| id | check | result |",
        "|----|-------|--------|",
    ]
    for c in summary["checks"]:
        lines.append(f"| {c['id']} | {c['name']} | "
                     f"{'PASS' if c['pass'] else 'FAIL'} |")

    quad = next(c for c in summary["checks"] if c["id"] == "T4")["detail"]
    lines += ["", "## The turn cost, against the integrated geodesic", "",
              f"`{quad['law']}` — and `{quad['conductance_form']}`, "
              f"so the geometric hinge and the monopole channel are one "
              f"integral.  `{quad['long_arm_limit']}`", "",
              "| `f₀` | `α` | `T(α)` | `α²/(2I)` | ratio | `f₀α` | `T/(f₀α)` |",
              "|--|--|--|--|--|--|--|"]
    for r in quad["rows"]:
        lines.append(
            f"| `{r['neck']:.0e}` | `{r['alpha']:.4f}` | `{r['turn_cost']:.5e}` "
            f"| `{r['small_angle_law']:.5e}` | `{r['law_ratio']:.6f}` "
            f"| `{r['linear_guess']:.2e}` "
            f"| `{r['fraction_of_the_linear_guess']:.5f}` |")
    lines += ["",
              f"**Proposed:** {quad['what_was_proposed']}.",
              "",
              f"**Measured:** {quad['what_is_measured']}.",
              ""]

    same = quad["the_same_I_as_physical_throat"]
    lines += [
        f"The resistance is not a new quantity: at `a = {same['mouth_radius']}` "
        f"this module gives `{same['here']:.10e}` and "
        f"`physical_throat.resistance()` gives "
        f"`{same['physical_throat']:.10e}` — difference "
        f"`{same['difference']:.1e}`.", ""]

    indep = next(c for c in summary["checks"] if c["id"] == "T2")["detail"]
    lines += ["## The arms do not have to match", "",
              "| `f_i` | `L_o` | `L_i` | `L_o/L_i` | `T(1 rad)` | `α²/(2I)` |",
              "|--|--|--|--|--|--|"]
    for r in indep["rows"]:
        lines.append(f"| `{r['inner_scale']:.4g}` | `{r['outer_length']:.6f}` "
                     f"| `{r['inner_length']:.6f}` | `{r['length_ratio']:.2f}` "
                     f"| `{r['turn_cost_one_radian']:.4e}` "
                     f"| `{r['small_angle_law']:.4e}` |")
    lines += ["", indep["why_it_matters"].capitalize() + ".", ""]

    meet = next(c for c in summary["checks"] if c["id"] == "T8")["detail"]
    lines += ["## What intersection becomes", "",
              f"**Criterion.** {meet['criterion']} — {meet['whether_is_an_angular_question']}.",
              "",
              "| `f₀` | separation | `w_a` | `w_b` | meet? | overlap on the bearing | gap |",
              "|--|--|--|--|--|--|--|"]
    for r in meet["rows"]:
        lines.append(
            f"| `{r['neck']:.0e}` | `{r['angular_separation']:.2f}` "
            f"| `{r['angular_width_a']:.2f}` | `{r['angular_width_b']:.2f}` "
            f"| {'**yes**' if r['they_meet'] else 'no'} "
            f"| `{r['overlap_length_on_the_bearing']:.2e}` "
            f"| `{r['gap_length_on_the_bearing']:.2e}` |")
    lines += ["", f"**The point limit, correctly stated.** "
                  f"{meet['the_point_limit_correctly_stated']}.", ""]
    return "\n".join(lines)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = (Path(__file__).resolve().parent / "runs"
           / f"{stamp}_regularized_center_probe")
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0 if summary["passed"] == summary["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
