"""The center as a finite bearing.

The discipline here is the one the arc settled into: every closed form is
checked against an integral that does not use it, and the two claims that the
round is *for* — that the turn cost is quadratic and that the hinge stays
cheap — are checked against the integrated geodesic rather than against the
expansion that predicts them.
"""

import math

import numpy as np
import pytest

from geometrodynamics.viz import regularized_center as rc
from geometrodynamics.waves import physical_throat as pt


# ── the closed forms are the integrals ──────────────────────────────────────
@pytest.mark.parametrize("neck", [1e-4, 1e-2, 0.1])
@pytest.mark.parametrize("scale", [0.5, 1.0])
def test_the_arm_length_is_the_proper_distance(neck, scale):
    c = rc.RegularizedCenter(neck=neck, outer=scale, inner=scale)
    assert abs(c.outer_length() / c.length_by_quadrature() - 1.0) < 1e-9


@pytest.mark.parametrize("neck", [1e-4, 1e-2, 0.1])
@pytest.mark.parametrize("scale", [0.5, 1.0])
def test_the_arm_resistance_is_the_integral_of_ds_over_f_squared(neck, scale):
    c = rc.RegularizedCenter(neck=neck, outer=scale, inner=scale)
    closed = rc.arm_resistance(scale, neck)
    assert abs(closed / c.resistance_by_quadrature() - 1.0) < 1e-9


@pytest.mark.parametrize("a", [0.02, 0.05, 0.10, 0.20])
def test_the_symmetric_case_is_the_repos_throat_bit_for_bit(a):
    """The regularised centre has to *contain* the geometry PR #265 forced.

    Not "agree with to a tolerance" — the arm formula is the throat's own
    formula with the symmetry assumption dropped, so with the symmetry put
    back it must return the same floats.
    """
    t = pt.VacuumThroat(mouth_radius=a)
    c = rc.RegularizedCenter(neck=t.neck_radius(), outer=t.mouth_f(),
                             inner=t.mouth_f())
    assert abs(c.arm_length_sum() - t.length()) < 1e-15
    assert abs(c.resistance() - t.resistance()) < 1e-8
    assert abs(c.conductance() / t.conductance() - 1.0) < 1e-12


def test_the_closed_form_measurement_reports_both_checks():
    got = rc.measure_the_arm_length_is_the_repos_own_formula()
    assert got["the_closed_forms_are_the_integrals"]
    assert got["it_contains_the_symmetric_throat_exactly"]
    assert got["worst_length_quadrature_error"] < 1e-9
    assert got["worst_resistance_quadrature_error"] < 1e-9


def test_an_arm_narrower_than_the_neck_is_refused():
    with pytest.raises(ValueError):
        rc.arm_length(1e-4, 1e-3)
    with pytest.raises(ValueError):
        rc.arm_resistance(1e-4, 1e-3)


def test_the_singular_limit_of_the_arm_length_is_the_scale_itself():
    assert rc.arm_length(0.7, 0.0) == 0.7
    for f0 in (1e-4, 1e-6, 1e-8):
        assert abs(rc.arm_length(1.0, f0) - 1.0) < 1e-3
        assert rc.arm_length(1.0, f0) > 1.0, "the arm is longer, never shorter"


# ── the arms are independent ────────────────────────────────────────────────
def test_the_two_arms_need_not_match():
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.01)
    assert c.outer_length() > 80.0 * c.inner_length()
    got = rc.measure_the_two_arms_are_independent()
    assert got["the_arms_can_be_very_unequal"]
    assert got["and_nothing_about_the_hinge_changes"]
    assert got["largest_arm_ratio"] > 400.0


def test_the_profile_is_monotone_on_each_arm_and_pinches_at_the_neck():
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    s, f = c.profile(n=300)
    assert np.all(np.diff(s) > 0.0), "s must be a coordinate, not fold back"
    assert abs(f.min() - c.neck) < 1e-15
    j = int(np.argmin(f))
    assert np.all(np.diff(f[:j]) < 0.0) and np.all(np.diff(f[j:]) > 0.0)
    # the two arms really are different lengths, and in the drawn object
    assert abs(-s[0] - c.outer_length()) < 1e-12
    assert abs(s[-1] - c.inner_length()) < 1e-12
    assert -s[0] > s[-1], "the outer arm is the longer one here"


# ── the width rescales ──────────────────────────────────────────────────────
def test_the_physical_width_follows_the_profile():
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    dtheta = 0.2
    wo = c.feature_width(c.outer, dtheta)
    wm = c.feature_width(c.neck, dtheta)
    wi = c.feature_width(c.inner, dtheta)
    assert wo > wi > wm
    assert abs(wo / wm - c.outer / c.neck) < 1e-9


def test_the_end_to_end_width_ratio_is_the_scale_ratio_to_one_ulp():
    """An identity in the field, one rounding in the drawn numbers.

    `w_i/w_o` and `f_i/f_o` are the same real number, but the drawn ratio
    multiplies by `Δθ` twice first.  Asserting bitwise equality here would be
    asserting something about float multiplication, not about the geometry —
    the arc has made that mistake before, so this pins the ulp instead.
    """
    got = rc.measure_the_width_rescales_with_the_profile()
    assert got["the_ratio_is_the_scale_ratio_to_the_last_ulp"]
    assert got["worst_residue_in_ulps"] <= 1.0
    assert got["worst_residue_in_ulps"] > 0.0, "and it is not bitwise"
    for row in got["rows"]:
        assert abs(row["end_to_end_ratio"] / row["scale_ratio"] - 1.0) < 1e-15


def test_the_width_ratio_does_not_involve_the_neck():
    a = rc.RegularizedCenter(neck=1e-2, outer=1.0, inner=0.35)
    b = rc.RegularizedCenter(neck=1e-6, outer=1.0, inner=0.35)
    ra = a.feature_width(a.inner, 0.3) / a.feature_width(a.outer, 0.3)
    rb = b.feature_width(b.inner, 0.3) / b.feature_width(b.outer, 0.3)
    assert ra == rb


# ── the result: the turn cost is quadratic ──────────────────────────────────
@pytest.mark.parametrize("alpha", [0.02, 0.05, 0.1, 0.3])
def test_the_turn_cost_matches_the_small_angle_law(alpha):
    """`T(α) = α²/(2I)` against the *integrated* geodesic, not against itself."""
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    assert abs(c.turn_cost(alpha) / c.turn_cost_small_angle(alpha) - 1.0) < 1e-3


def test_the_turn_cost_is_quadratic_and_not_linear():
    """Doubling the angle must quadruple the cost, not double it."""
    c = rc.RegularizedCenter(neck=1e-4, outer=1.0, inner=0.35)
    for base in (0.02, 0.05, 0.1):
        ratio = c.turn_cost(2.0 * base) / c.turn_cost(base)
        assert abs(ratio - 4.0) < 0.02, "quadratic"
        assert abs(ratio - 2.0) > 1.0, "and emphatically not linear"


def test_the_turn_cost_is_linear_in_the_neck_radius():
    alpha = 1.0
    scaled = [rc.RegularizedCenter(neck=f0, outer=1.0,
                                   inner=0.35).turn_cost(alpha) / f0
              for f0 in (1e-3, 1e-4, 1e-5)]
    assert max(scaled) / min(scaled) - 1.0 < 1e-2


def test_the_long_arm_prefactor_is_one_eighth():
    """`I → 4/f₀` for two long arms, so `T → f₀α²/8`."""
    c = rc.RegularizedCenter(neck=1e-6, outer=1.0, inner=1.0)
    assert abs(c.resistance() * c.neck - 4.0) < 1e-4
    for alpha in (0.05, 0.1, 0.2):
        assert abs(c.turn_cost(alpha) / (c.neck * alpha ** 2 / 8.0) - 1.0) < 1e-3


def test_the_turn_cost_uses_the_same_resistance_as_the_physical_throat():
    """The geometric hinge cost and the monopole conductance are one integral."""
    t = pt.VacuumThroat(mouth_radius=0.05)
    c = rc.RegularizedCenter(neck=t.neck_radius(), outer=t.mouth_f(),
                             inner=t.mouth_f())
    assert abs(c.resistance() - t.resistance()) < 1e-8
    alpha = 0.1
    via_conductance = alpha ** 2 * t.conductance() / (8.0 * math.pi)
    assert abs(c.turn_cost_small_angle(alpha) / via_conductance - 1.0) < 1e-12
    assert abs(c.turn_cost(alpha) / via_conductance - 1.0) < 1e-3


def test_the_quadratic_measurement_reports_the_correction_it_makes():
    got = rc.measure_the_turn_cost_is_quadratic_not_linear()
    assert got["the_small_angle_law_holds"]
    assert got["the_prefactor_is_one_eighth"]
    assert got["the_resistance_is_the_repos"]
    assert got["it_is_still_close_at_pi"]
    assert got["worst_small_angle_error"] < 1e-3
    assert "f0 alpha" in got["what_was_proposed"]
    assert "alpha^2" in got["what_is_measured"]


# ── the corner route is an upper bound ──────────────────────────────────────
@pytest.mark.parametrize("alpha", [0.1, 1.0, math.pi])
@pytest.mark.parametrize("neck", [1e-3, 1e-5])
def test_the_corner_route_bounds_the_geodesic_from_above(alpha, neck):
    c = rc.RegularizedCenter(neck=neck, outer=1.0, inner=0.35)
    assert c.corner_route_length(alpha) > c.route_length(alpha)
    assert c.route_length(alpha) > c.arm_length_sum()


def test_the_corner_route_is_exactly_arms_plus_arc():
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    for alpha in (0.3, 1.0, math.pi):
        assert abs(c.corner_route_length(alpha)
                   - (c.arm_length_sum() + c.neck * alpha)) < 1e-15


def test_how_loose_the_bound_is_depends_on_the_angle_and_not_the_neck():
    got = rc.measure_the_corner_route_is_only_an_upper_bound()
    assert got["and_is_an_upper_bound_on_the_geodesic"]
    assert got["the_fraction_is_a_function_of_alpha_alone"]
    assert got["fraction_drift_between_two_decades_of_neck"] < 1e-2
    # the geodesic spends a small fraction of the arc at small angle
    assert got["fraction_of_the_arc_actually_spent_at_alpha_0p1"] < 0.02
    assert got["fraction_at_pi"] < 0.4


# ── the hinge stays cheap ───────────────────────────────────────────────────
def test_turning_is_never_the_expensive_part_of_the_route():
    """And `α = π` is the worst case, not a sample of one.

    `bearing_distance` reduces any pair of directions to `[0, π]`, so there is
    no separation past it to check. The withdrawn "break-even angle" of 104 rad
    extrapolated the small-angle law outside both its domain and the bearing's
    configuration space.
    """
    got = rc.measure_the_hinge_is_never_the_expensive_part()
    assert got["the_hinge_is_always_cheap"]
    assert got["no_reachable_orientation_breaks_even"]
    assert got["worst_case_turn_over_arms"] < 1e-2
    assert "withdrawn" in " ".join(got.keys())
    assert "not a physical angle" in got["the_break_even_angle_is_withdrawn"]
    for row in got["rows"]:
        assert row["turn_over_arms"] < 1e-2
    # pi really is the ceiling the measurement claims it is
    for alpha in (3.5, 5.0, 2.0 * math.pi + 1.0):
        assert rc.bearing_distance(alpha) <= math.pi + 1e-12


def test_the_cheapness_is_stated_without_claiming_it_is_free():
    got = rc.measure_the_hinge_is_never_the_expensive_part()
    assert "not say" in " ".join(got.keys()).replace("_", " ")
    caveat = got["what_this_does_not_say"].lower()
    assert "zero" in caveat and "linear in f0" in caveat


# ── the point limit ─────────────────────────────────────────────────────────
def test_the_bearing_shrinks_back_to_the_point():
    got = rc.measure_the_bearing_shrinks_back_to_the_point()
    assert got["turn_cost_is_linear_in_the_neck"]
    assert got["the_arm_asymptotic_holds"]
    rows = got["rows"]
    assert rows[0]["turn_cost"] > 100.0 * rows[-1]["turn_cost"]
    assert rows[-1]["arm_length_sum"] < rows[0]["arm_length_sum"]


def test_the_arm_excess_asymptotic():
    """`L(F) − F ≃ (f₀/2)[ln(4F/f₀) − 1]` — an independent closed form."""
    for f0 in (1e-4, 1e-6, 1e-8):
        excess = rc.arm_length(1.0, f0) - 1.0
        predicted = 0.5 * f0 * (math.log(4.0 / f0) - 1.0)
        assert abs(excess / predicted - 1.0) < 2e-2


def test_route_length_at_zero_angle_is_just_the_arms():
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    assert c.route_length(0.0) == c.arm_length_sum()
    assert c.turn_cost(0.0) == 0.0
    assert c.clairaut_constant(0.0) == 0.0


# ── intersection becomes overlap ────────────────────────────────────────────
def test_whether_two_fronts_meet_is_an_angular_question():
    got = rc.measure_the_bearing_replaces_collision_with_overlap()
    assert got["the_verdict_does_not_depend_on_the_neck"]
    assert got["both_the_overlap_and_the_gap_scale_with_the_neck"]
    assert got["no_route_passes_through_a_singular_point"]
    # the same angular configuration must give the same verdict at either neck
    by_config = {}
    for r in got["rows"]:
        key = (r["angular_separation"], r["angular_width_a"],
               r["angular_width_b"])
        by_config.setdefault(key, set()).add(r["they_meet"])
    assert all(len(v) == 1 for v in by_config.values())
    assert {True, False} == {r["they_meet"] for r in got["rows"]}


def test_the_point_limit_does_not_make_everything_meet():
    """The correction to the obvious reading of the `f₀ → 0` limit.

    Shrinking the bearing does not merge every angular position into one
    meeting; it shrinks the overlap *and* the gap together. So a pair that
    misses at `f₀ = 1e-2` still misses at `1e-6` — its separation just becomes
    unmeasurably small.
    """
    got = rc.measure_the_bearing_replaces_collision_with_overlap()
    misses = [r for r in got["rows"] if not r["they_meet"]]
    assert misses, "the scan has to contain a pair that misses"
    for r in misses:
        assert r["gap_angle"] > 0.0
        assert r["overlap_length_on_the_bearing"] == 0.0
        assert r["gap_length_on_the_bearing"] == r["neck"] * r["gap_angle"]
    text = got["the_point_limit_correctly_stated"].lower()
    assert "does not" in text and "yes/no" in text


# ── the bearing's dimension does not matter, but the identification does ────
def test_the_turn_cost_reduces_to_a_great_circle_on_any_sphere():
    got = rc.measure_the_turn_cost_does_not_care_which_sphere()
    assert got["the_sphere_and_the_drawn_circle_agree"]
    assert got["worst_great_circle_cost_difference"] < 1e-9
    for row in got["great_circle_reduction"]:
        assert abs(row["angle_from_the_dot_product"]
                   - row["angle_drawn_on_the_2d_circle"]) < 1e-12


def test_the_projective_identification_is_what_changes_the_answer():
    got = rc.measure_the_turn_cost_does_not_care_which_sphere()
    assert got["the_projective_identification_matters"]
    assert got["the_projective_saving_near_pi"] > 100.0
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    # below pi/2 the identification does nothing at all
    for alpha in (0.2, 1.0, math.pi / 2):
        assert c.turn_cost(alpha, projective=True) == c.turn_cost(alpha)
    # and above it, reversing an unoriented axis is nearly free
    assert c.turn_cost(3.0, projective=True) < 0.01 * c.turn_cost(3.0)


@pytest.mark.parametrize("alpha,expected", [
    (0.0, 0.0), (0.5, 0.5), (math.pi / 2, math.pi / 2),
    (math.pi, 0.0), (2.0 * math.pi, 0.0), (-0.5, 0.5),
])
def test_the_bearing_distance_wraps_and_identifies(alpha, expected):
    assert abs(rc.bearing_distance(alpha, projective=True) - expected) < 1e-12


def test_the_oriented_bearing_distance_wraps_but_does_not_identify():
    assert abs(rc.bearing_distance(math.pi) - math.pi) < 1e-12
    assert abs(rc.bearing_distance(1.5 * math.pi) - 0.5 * math.pi) < 1e-12
    assert abs(rc.bearing_distance(-2.0) - 2.0) < 1e-12


# ── the geodesic machinery itself ───────────────────────────────────────────
def test_only_geodesics_below_the_neck_scale_cross_at_all():
    """Clairaut's `κ = h/f₀` must land inside `(0, 1)`: `h > f₀` turns back."""
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    previous = 0.0
    for alpha in (0.1, 0.5, 1.0, 2.0, math.pi):
        k = c.clairaut_constant(alpha)
        assert 0.0 < k < 1.0
        assert k > previous, "a wider sweep needs a larger kappa"
        previous = k


def test_the_swept_angle_diverges_as_kappa_approaches_one():
    """Which is why every angle is reachable and the bracket is safe."""
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    swept = [c._sweep(k) for k in (0.5, 0.9, 0.99, 0.9999)]
    assert all(b > a for a, b in zip(swept, swept[1:]))
    assert swept[-1] > 3.0 * math.pi


def test_the_arm_extent_in_x_is_the_throats_own_half_length():
    """`X = arcosh√(F/f₀)`, and at the symmetric point `arcosh(1/sin a)`."""
    for a in (0.02, 0.05, 0.20):
        t = pt.VacuumThroat(mouth_radius=a)
        c = rc.RegularizedCenter(neck=t.neck_radius(), outer=t.mouth_f(),
                                 inner=t.mouth_f())
        assert abs(c.half_length_in_x(t.mouth_f())
                   - math.acosh(1.0 / math.sin(a))) < 1e-13
        # and PR #267's identity: tanh X = cos a, e^-X = tan(a/2)
        x = c.half_length_in_x(t.mouth_f())
        assert abs(math.tanh(x) - math.cos(a)) < 1e-13
        assert abs(math.exp(-x) - math.tan(0.5 * a)) < 1e-13


def test_tanh_of_the_arm_extent_is_the_resistance_factor():
    """`tanh X = √(1 − f₀/F)`, which is what ties `T` to `I` analytically."""
    f0 = 1e-3
    for F in (0.01, 0.35, 1.0):
        c = rc.RegularizedCenter(neck=f0, outer=F, inner=F)
        x = c.half_length_in_x(F)
        assert abs(math.tanh(x) - math.sqrt(1.0 - f0 / F)) < 1e-13
        assert abs(rc.arm_resistance(F, f0) - (2.0 / f0) * math.tanh(x)) < 1e-9


@pytest.mark.parametrize("alpha", [0.02, 0.1, 1.0, math.pi])
def test_the_turn_cost_is_stable_across_seven_decades_of_neck(alpha):
    """The formulation check that the first draft failed.

    `T` is as small as `1e-15` while the arms are `O(1)`, so building it as a
    difference is all cancellation, and integrating in `f` hides a spike of
    width `√f₀` from an adaptive quadrature. Written as a direct integral in
    `x`, the shape `T/(α²/2I)` settles to a limit and stays there. Both failed
    drafts break this test.

    The shape is *not* exactly `f₀`-independent, and the tolerance says so.
    It depends on the arm extents `X = arcosh√(F/f₀)`, which grow as the neck
    shrinks, so it converges rather than being constant: at `α = π` it runs
    `0.92309 → 0.92478 → 0.92495 → 0.92497` over `f₀ = 1e-2 … 1e-5` and then
    holds. The decades checked here are the converged ones.
    """
    shapes = []
    for f0 in (1e-5, 1e-7, 1e-10):
        c = rc.RegularizedCenter(neck=f0, outer=1.0, inner=0.35)
        shapes.append(c.turn_cost(alpha) / c.turn_cost_small_angle(alpha))
    assert max(shapes) / min(shapes) - 1.0 < 1e-5
    assert all(0.9 < s <= 1.0 + 1e-9 for s in shapes)


@pytest.mark.parametrize("alpha", [0.1, 1.0, math.pi])
def test_the_shape_converges_as_the_arms_outgrow_the_neck(alpha):
    """And the approach is monotone, which is what says it is a limit."""
    shapes = [rc.RegularizedCenter(neck=f0, outer=1.0, inner=0.35).turn_cost(alpha)
              / rc.RegularizedCenter(neck=f0, outer=1.0,
                                     inner=0.35).turn_cost_small_angle(alpha)
              for f0 in (1e-2, 1e-3, 1e-4, 1e-5)]
    assert all(b > a for a, b in zip(shapes, shapes[1:])), "monotone"
    steps = [b - a for a, b in zip(shapes, shapes[1:])]
    assert all(b < a for a, b in zip(steps, steps[1:])), "and decelerating"


def test_the_turn_cost_is_never_built_by_subtraction():
    """`route_length` is arms + hinge, not hinge = route − arms."""
    c = rc.RegularizedCenter(neck=1e-9, outer=1.0, inner=0.35)
    alpha = 0.05
    t = c.turn_cost(alpha)
    assert t > 0.0
    assert abs(t / c.turn_cost_small_angle(alpha) - 1.0) < 1e-3
    # a subtraction at this neck would be pure noise: check it would have been
    naive = (c.arm_length_sum() + t) - c.arm_length_sum()
    assert abs(naive / t - 1.0) > 1e-4, (
        "if the difference happened to be accurate here the test proves "
        "nothing -- it must not be")
    assert c.route_length(alpha) == c.arm_length_sum() + t


def test_the_route_length_grows_with_the_angle():
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    lengths = [c.route_length(a) for a in (0.0, 0.3, 1.0, 2.0, math.pi)]
    assert all(b > a for a, b in zip(lengths, lengths[1:]))


# ── what is general, and what is only about this neck ───────────────────────
def test_the_quadratic_law_holds_on_a_second_unrelated_profile():
    """`T = α²/(2I)` never used the profile, so it must survive changing it.

    `f = √(f₀² + s²)` is not scalar-flat and is not glued to anything. The law
    holds there to `1.3e-04`, across six decades of `f₀`.
    """
    got = rc.measure_the_law_does_not_depend_on_the_profile()
    assert got["the_law_holds_on_the_second_profile"]
    assert got["worst_small_angle_error"] < 1e-3
    for row in got["rows"]:
        if row["alpha"] <= 0.1:
            assert abs(row["shape"] - 1.0) < 1e-3


def test_but_the_large_angle_correction_is_profile_dependent():
    """The leading law is about necks; the deficit at a half turn is not.

    `0.925` on the scalar-flat profile against `0.889` on the hyperbolic one —
    same `α²/(2I)`, different `O(α⁴)`. Worth pinning because it marks exactly
    how far the generality claim reaches.
    """
    got = rc.measure_the_law_does_not_depend_on_the_profile()
    assert got["the_correction_is_profile_dependent"]
    here = got["shape_at_pi_here"]
    flat = got["shape_at_pi_scalar_flat"]
    assert 0.88 < here < 0.90
    assert 0.92 < flat < 0.93
    assert abs(here - flat) > 0.02


def test_the_second_profile_has_its_own_resistance():
    """`I = (1/f₀)∫dx/cosh x` for `f = √(f₀²+s²)`, and `→ π/f₀` for long arms."""
    got = rc.hyperbolic_neck(1e-6, 1.0, 1.0, 0.1)
    assert abs(got["resistance"] * 1e-6 - math.pi) < 1e-4
    assert 0.0 < got["kappa"] < 1.0
    # and it is a different number from the scalar-flat neck's 4/f0
    flat = rc.RegularizedCenter(neck=1e-6, outer=1.0, inner=1.0)
    assert abs(flat.resistance() * 1e-6 - 4.0) < 1e-4
    assert abs(got["resistance"] / flat.resistance() - math.pi / 4.0) < 1e-3


# ── the moment hierarchy ────────────────────────────────────────────────────
@pytest.mark.parametrize("neck", [1e-2, 1e-4])
@pytest.mark.parametrize("scale", [0.5, 1.0])
@pytest.mark.parametrize("order", [2, 4, 6])
def test_the_moments_are_the_integrals(neck, scale, order):
    """`I_n = (2/f₀^{n−1})∫₀^T(1−t²)^{n−2}dt` against direct quadrature."""
    from scipy.integrate import quad
    closed = rc.arm_moment(scale, neck, order)
    got = quad(lambda t: 2.0 * math.sqrt(neck + t * t)
               / (neck + t * t) ** order,
               0.0, math.sqrt(scale - neck), limit=600)[0]
    assert abs(closed / got - 1.0) < 1e-9


def test_the_second_moment_is_the_resistance():
    c = rc.RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    assert c.moment(2) == c.resistance()
    assert rc.arm_moment(1.0, 1e-3, 2) == rc.arm_resistance(1.0, 1e-3)


def test_odd_and_degenerate_moment_orders_are_refused():
    for bad in (0, 1, 3, -2):
        with pytest.raises(ValueError):
            rc.arm_moment(1.0, 1e-3, bad)


def test_the_long_arm_moments_have_their_closed_forms():
    """`I₂ → 4/f₀` and `I₄ → 32/(15f₀³)` for two long scalar-flat arms."""
    c = rc.RegularizedCenter(neck=1e-6, outer=1.0, inner=1.0)
    assert abs(c.resistance() * c.neck - 4.0) < 1e-5
    assert abs(c.fourth_moment() * c.neck ** 3 - 32.0 / 15.0) < 1e-9


@pytest.mark.parametrize("alpha", [0.1, 0.3, 0.6])
def test_the_fourth_order_expansion_beats_the_second(alpha):
    """`α²/(2I₂) − α⁴I₄/(8I₂⁴)` must be closer to the geodesic than `α²/(2I₂)`."""
    c = rc.RegularizedCenter(neck=1e-6, outer=1.0, inner=1.0)
    exact = c.turn_cost(alpha)
    second = abs(exact / c.turn_cost_small_angle(alpha) - 1.0)
    fourth = abs(exact / c.turn_cost_to_fourth_order(alpha) - 1.0)
    assert fourth < second
    assert fourth < 1e-4


def test_the_shape_coefficient_is_one_over_one_hundred_and_twenty():
    """Scalar-flat long arms: shape `= 1 − α²/120`, an exact rational.

    The residual is not asserted against a flat tolerance but against its own
    scaling: `shape = 1 − α²/120 + O(α⁴)`, so the miss must grow as `α⁴` — a
    factor of 16 between `α = 0.3` and `α = 0.6`. A flat bound would pass for
    the wrong reason at small angle and fail for the right one at large.
    """
    c = rc.RegularizedCenter(neck=1e-6, outer=1.0, inner=1.0)
    coeff = c.fourth_moment() / (4.0 * c.resistance() ** 3)
    assert abs(coeff - 1.0 / 120.0) < 1e-7
    misses = []
    for alpha in (0.15, 0.3, 0.6):
        miss = abs(c.turn_cost(alpha) / c.turn_cost_small_angle(alpha)
                   - (1.0 - alpha ** 2 / 120.0))
        assert miss < 2e-5
        misses.append(miss / alpha ** 4)
    assert max(misses) / min(misses) - 1.0 < 0.05, "the residual is O(alpha^4)"


def test_the_two_shapes_do_not_agree_to_eight_digits():
    """A withdrawn claim, pinned so it cannot come back.

    `8` digits is how well each profile matches *its own* quartic law; the two
    profiles match *each other* only to `4.3e-05` at `α = 0.1`. The separation
    is `α²(1/79 − 1/120)` and grows quadratically.
    """
    got = rc.measure_the_fourth_moment_is_where_the_neck_shape_enters()
    assert got["the_two_profiles_do_not_agree_to_eight_digits"]
    assert got["the_separation_grows_as_alpha_squared"]
    d = got["how_far_apart_at_alpha_0p1"]
    assert 4.0e-5 < d < 5.0e-5, "not eight digits, and not nothing either"
    coeff = 1.0 / (8.0 * math.pi ** 2) - 1.0 / 120.0
    for row in got["shape_difference"]:
        assert row["difference"] > 0.0, "the hyperbolic neck is always lower"
        if row["alpha"] <= 1.0:
            assert abs(row["difference"] / (row["alpha"] ** 2 * coeff)
                       - 1.0) < 0.05
    # and the thing that *is* good to eight digits is a different quantity
    flat = rc.RegularizedCenter(neck=1e-6, outer=1.0, inner=1.0)
    own = abs(flat.turn_cost(0.1) / flat.turn_cost_small_angle(0.1)
              - (1.0 - 0.01 / 120.0))
    assert own < 1e-8
    assert own < 1e-3 * d, "the two must not be confusable"


def test_i4_is_not_the_entire_profile_dependence():
    """`I₆` and beyond enter at `O(α⁶)` and matter by `α = π`."""
    got = rc.measure_the_fourth_moment_is_where_the_neck_shape_enters()
    assert "I6" in got["i4_is_not_the_entire_profile_dependence"]
    c = rc.RegularizedCenter(neck=1e-6, outer=1.0, inner=1.0)
    near = abs(c.turn_cost(0.1) / c.turn_cost_to_fourth_order(0.1) - 1.0)
    far = abs(c.turn_cost(math.pi) / c.turn_cost_to_fourth_order(math.pi) - 1.0)
    assert near < 1e-8
    assert far > 1e-3, "the quartic truncation visibly fails at pi"


def test_the_universality_claim_is_about_the_form_not_the_number():
    """`I₂` itself is profile-dependent; the functional form is what is shared."""
    got = rc.measure_the_fourth_moment_is_where_the_neck_shape_enters()
    text = got["the_division_of_labour"]
    assert "LEADING FUNCTIONAL FORM" in text
    assert "NOT itself universal" in text
    flat = rc.RegularizedCenter(neck=1e-6, outer=1.0, inner=1.0).resistance()
    hyper = rc.hyperbolic_neck(1e-6, 1.0, 1.0, 0.1)["resistance"]
    assert abs(flat * 1e-6 - 4.0) < 1e-4
    assert abs(hyper * 1e-6 - math.pi) < 1e-4
    assert abs(flat / hyper - 4.0 / math.pi) < 1e-4, "different numbers"


def test_the_two_profiles_differ_first_at_the_fourth_moment():
    got = rc.measure_the_fourth_moment_is_where_the_neck_shape_enters()
    assert got["the_shape_law_holds_at_small_angle"]
    assert got["the_fourth_order_beats_the_second"]
    assert got["the_second_moment_is_shared_the_fourth_is_not"]
    # I2 differs only by the profile's overall scale (4/f0 vs pi/f0);
    # the SHAPE coefficient is where they part company
    flat = got["scalar_flat"]["shape_coefficient"]
    hyper = got["hyperbolic"]["shape_coefficient"]
    assert abs(flat - 1.0 / 120.0) < 1e-6
    assert abs(hyper - 1.0 / (8.0 * math.pi ** 2)) < 1e-6
    assert hyper > 1.4 * flat, "and they are not close"
    for row in got["hyperbolic_rows"]:
        if row["alpha"] <= 0.1:
            assert abs(row["shape_measured"] - row["shape_predicted"]) < 1e-6


# ── the identity underneath ─────────────────────────────────────────────────
def test_the_dirichlet_identity_is_a_two_dimensional_cross_section_statement():
    """The two weights are not automatically the same.

    The azimuth's weight is the metric coefficient `f²` for every bearing
    dimension; the monopole's is the volume element `f^q`. They coincide at
    `q = 2` and nowhere else — which is the physical case here, but is a fact
    about that dimension rather than about the construction.
    """
    got = rc.measure_the_hinge_and_the_monopole_are_one_dirichlet_form()
    assert got["the_identity_is_a_q_equals_two_statement"]
    by_dim = {d["angular_dimension"]: d for d in got["by_angular_dimension"]}
    assert by_dim[2]["they_coincide"]
    assert not by_dim[1]["they_coincide"]
    assert not by_dim[3]["they_coincide"]
    assert by_dim[1]["monopole_resistance"] < by_dim[2]["monopole_resistance"]
    assert by_dim[3]["monopole_resistance"] > by_dim[2]["monopole_resistance"]
    assert "not dimension-free" in got["and_q_equals_two_is_the_physical_case"]


@pytest.mark.parametrize("q", [1, 2, 3, 4])
def test_the_monopole_resistance_is_the_integral_of_ds_over_f_to_the_q(q):
    from scipy.integrate import quad
    f0, F = 1e-3, 0.5
    closed = rc.monopole_resistance(F, f0, q)
    got = quad(lambda t: 2.0 * math.sqrt(f0 + t * t) / (f0 + t * t) ** q,
               0.0, math.sqrt(F - f0), limit=600)[0]
    assert abs(closed / got - 1.0) < 1e-8


def test_a_zero_dimensional_cross_section_is_refused():
    with pytest.raises(ValueError):
        rc.monopole_resistance(0.5, 1e-3, 0)


def test_the_overlap_size_law_is_scoped_by_dimension():
    got = rc.measure_the_bearing_replaces_collision_with_overlap()
    text = got["how_big_depends_on_the_bearing_dimension"]
    assert "LENGTH" in text and "AREA" in text and "VOLUME" in text
    assert "criterion is dimension-free" in text
    sizes = {r["angular_dimension"]: r for r in got["overlap_size_by_dimension"]}
    assert sizes[1]["example_at_neck_1e-3_overlap_0p225"] > \
        sizes[2]["example_at_neck_1e-3_overlap_0p225"] > \
        sizes[3]["example_at_neck_1e-3_overlap_0p225"]
    # the yes/no verdict is unaffected by any of this
    assert got["the_verdict_does_not_depend_on_the_neck"]


# ── where the turning actually happens ──────────────────────────────────────
def test_the_turn_is_concentrated_at_the_neck_not_up_the_arm():
    """Correcting this module's own first explanation.

    `θ′ = h/f²` puts the angular rate highest where `f` is smallest, so the
    geodesic hugs the neck — it does not turn "up the arm where the lever is
    longer", and an angular increment at larger `f` costs *more* arc, not less.
    """
    got = rc.measure_where_the_turning_happens()
    assert got["the_turn_is_concentrated_at_the_neck"]
    assert got["the_geodesic_is_cheaper"]
    assert got["fraction_done_by_2p4_necks"] > 0.7
    rows = got["rows"]
    assert all(b["fraction_of_the_arms_turn_done"]
               >= a["fraction_of_the_arms_turn_done"]
               for a, b in zip(rows, rows[1:])), "monotone outward"
    assert rows[0]["fraction_of_the_arms_turn_done"] > 0.4, \
        "nearly half the turn is done within 1.3 neck radii"


def test_the_stated_reason_is_pythagoras_not_leverage():
    got = rc.measure_where_the_turning_happens()
    assert "Pythagoras" in got["the_corrected_reason"]
    assert "PURE transverse" in got["the_corrected_reason"]
    assert "backwards" in got["what_the_first_draft_said"]
    # and the arithmetic behind it: the corner is first order, the tilt second
    c = rc.RegularizedCenter(neck=1e-5, outer=1.0, inner=1.0)
    corner = [c.neck * a for a in (0.05, 0.1)]
    tilt = [c.turn_cost(a) for a in (0.05, 0.1)]
    assert abs(corner[1] / corner[0] - 2.0) < 1e-9, "the arc is linear"
    assert abs(tilt[1] / tilt[0] - 4.0) < 1e-3, "the tilt is quadratic"


def test_the_monopole_and_the_azimuth_are_one_dirichlet_form():
    """Normalised, the static potential and the geodesic azimuth coincide."""
    got = rc.measure_the_hinge_and_the_monopole_are_one_dirichlet_form()
    assert got["the_profiles_coincide_as_alpha_vanishes"]
    assert got["the_deviation_is_second_order_in_alpha"]
    assert got["the_two_readings_agree"]
    assert got["conductance_difference"] < 1e-12
    rows = got["profile_convergence"]
    assert all(b["worst_profile_deviation"] < a["worst_profile_deviation"]
               for a, b in zip(rows, rows[1:])), "monotone in alpha"
    assert rows[-1]["worst_profile_deviation"] < 1e-7


def test_the_deviation_between_the_two_profiles_falls_as_alpha_squared():
    """A clean factor of 100 per decade — which is what makes it `O(α²)`."""
    got = rc.measure_the_hinge_and_the_monopole_are_one_dirichlet_form()
    ratios = [r["over_alpha_squared"] for r in got["profile_convergence"]]
    assert max(ratios) / min(ratios) - 1.0 < 1e-2


def test_the_hinge_can_be_read_off_the_conductance_alone():
    """`T(α) = α²·(4π/I₂)/(8π)` — the same number by two routes."""
    t = pt.VacuumThroat(mouth_radius=0.05)
    c = rc.RegularizedCenter(neck=t.neck_radius(), outer=t.mouth_f(),
                             inner=t.mouth_f())
    for alpha in (0.01, 0.05, 0.2):
        via_conductance = alpha ** 2 * t.conductance() / (8.0 * math.pi)
        assert abs(via_conductance / c.turn_cost_small_angle(alpha)
                   - 1.0) < 1e-12


def test_the_identity_says_what_it_does_not_cover():
    got = rc.measure_the_hinge_and_the_monopole_are_one_dirichlet_form()
    caveat = got["what_is_not_claimed"].lower()
    assert "beyond leading order" in caveat
    assert "i4" in caveat and "static potential does not" in caveat


# ── the point limit, in the sharper words ───────────────────────────────────
def test_the_point_limit_separates_incidence_from_interaction_region():
    got = rc.measure_the_bearing_replaces_collision_with_overlap()
    assert "angular incidence" in got["what_survives_the_limit"]
    assert "f0" in got["what_survives_the_limit"]
    assert "interaction region" in got["what_collapses"]
    text = got["the_point_limit_correctly_stated"].lower()
    assert "angular incidence survives" in text
    assert "interaction region collapses" in text
    # and the arithmetic behind the words
    for row in got["rows"]:
        reach = 0.5 * (row["angular_width_a"] + row["angular_width_b"])
        assert row["they_meet"] == (row["angular_separation"] < reach)
        assert row["overlap_length_on_the_bearing"] <= row["neck"] * math.pi
