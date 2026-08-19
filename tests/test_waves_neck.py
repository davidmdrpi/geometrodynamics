"""Tests for the neck: the balls actually removed from the ambient.

The load-bearing ones here are not the sweeps.  They are the tests that pin the
*constructions* against each other — shooting against a closed form, the
quadratic form against Poincaré, roots against poles — because the round's
claim is a theorem and a theorem is only as good as the object it is about.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.finite_mouth import (
    FOUR_PI, FiniteMouthThroat, mouth_green, regular_radial,
)
from geometrodynamics.waves.neck import (
    NeckThroat, WORKING_NECK,
    exterior_dtn, exterior_dtn_monopole, exterior_log_derivative,
    measure_the_exterior_dtn_is_positive_in_every_channel,
    measure_the_higher_multipoles_decouple,
    measure_the_negative_mode_does_not_survive_the_neck,
    measure_the_quadratic_form_is_positive,
    measure_the_soft_mode_is_forced_by_the_two_ends,
    measure_the_soft_mode_survives_the_removal,
    measure_what_the_fixed_ambient_cost,
    rayleigh_quotient,
)


# ── the exterior map ───────────────────────────────────────────────────────
def test_the_shooting_integrator_reproduces_the_monopole_closed_form():
    """Two independent constructions of ``N₀``, agreeing to ``1e-13``.

    This is the test the rest of the module leans on.  The ``ℓ ≥ 1`` channels
    have no closed form and are known only through the integrator, so the
    integrator has to be checked somewhere it can be — and ``ℓ = 0`` is that
    place.
    """
    for a in (0.02, 0.2, 0.5, 1.0):
        for sigma in (0.05, 1.0, 5.0, 30.0):
            shot = -(FOUR_PI * math.sin(a) ** 2) * exterior_log_derivative(
                a, -sigma ** 2, 0)
            closed = exterior_dtn_monopole(a, -sigma ** 2)
            assert shot == pytest.approx(closed, rel=1e-11)


def test_the_exterior_map_is_positive_below_zero():
    """The ambient half of the answer, in every channel.

    ``v'' + [λ − ℓ(ℓ+2)/sin²χ]v = 0`` has a strictly negative bracket at
    ``λ < 0``, so ``v`` cannot oscillate and the map keeps a sign.
    """
    for a in (0.02, 0.2, 0.5):
        for sigma in (0.01, 1.0, 10.0, 100.0):
            for ell in (0, 1, 2, 5):
                assert exterior_dtn(a, -sigma ** 2, ell) > 0.0


def test_the_exterior_map_is_positive_inside_the_gap_too():
    """``N₀ > 0`` on ``0 < λ < 1`` as well, which is where the bound states are.

    Without this the crossing found in `bound_states` could be an artifact of
    the ambient block changing sign rather than a genuine state.
    """
    for a in (0.02, 0.05, 0.2, 0.5):
        for lam in np.linspace(1e-6, 1.0 - 1e-6, 200):
            assert exterior_dtn_monopole(a, float(lam)) > 0.0


def test_the_map_is_stiffer_at_higher_multipoles():
    for a in (0.05, 0.3):
        for sigma in (0.1, 3.0):
            values = [exterior_dtn(a, -sigma ** 2, ell)
                      for ell in (0, 1, 2, 3, 5)]
            assert all(b > c for c, b in zip(values, values[1:]))


def test_the_small_ball_limit_is_the_capacitance():
    """``N₀ → 4πa``, which is what fixes the normalization as physical."""
    ratios = [exterior_dtn_monopole(a, 0.0) / (FOUR_PI * a)
              for a in (0.02, 0.01, 0.005)]
    assert all(abs(r - 1.0) < 0.01 for r in ratios)
    # and the approach is linear in a, not merely close
    assert ratios[1] - 1.0 == pytest.approx(0.5 * (ratios[0] - 1.0), rel=0.05)


def test_the_gap_edge_slope_is_the_closed_form():
    """``N₀ → 2π(π − a + sin a cos a)(1 − λ)`` as ``λ → 1⁻``.

    This is why the soft mode exists: the exterior stiffness vanishes at the
    free ESU threshold, so the ambient diagonal ``1/N₀`` diverges there.
    """
    for a in (0.05, 0.2, 0.5, 1.0):
        slope = 2.0 * math.pi * (math.pi - a + math.sin(a) * math.cos(a))
        measured = exterior_dtn_monopole(a, 1.0 - 1e-6) / 1e-6
        assert measured == pytest.approx(slope, rel=1e-3)


def test_a_ball_outside_the_sphere_is_rejected():
    with pytest.raises(ValueError):
        exterior_log_derivative(4.0, -1.0, 0)


# ── the quadratic form ─────────────────────────────────────────────────────
def test_the_interior_only_configuration_hits_the_poincare_floor():
    """The degenerate case of the theorem, checked exactly.

    ``φ ≡ 0`` outside forces ``u(0) = u(L) = 0``, and the lowest such ``u`` is
    ``sin(πs/L)`` with quotient exactly ``π²/L²``.  If this did not land, the
    measures in `rayleigh_quotient` would be wrong and the theorem's degenerate
    case would be untested.
    """
    for length in (0.4, 0.9, 2.5):
        t = NeckThroat(separation=1.3, length=length, area=FOUR_PI,
                       radius=0.05)
        r = rayleigh_quotient(t, lambda chi: 0.0,
                              lambda s: math.sin(math.pi * s / length),
                              n=4001)
        assert r["quotient"] == pytest.approx((math.pi / length) ** 2,
                                              rel=1e-5)


def test_the_energy_and_norm_are_both_positive():
    rng = np.random.default_rng(7)
    for _ in range(12):
        c = rng.normal(size=3)

        def profile(chi, c=c):
            return float(c[0] + c[1] * chi + c[2] * chi ** 2)

        r = rayleigh_quotient(WORKING_NECK, profile,
                              lambda s: profile(WORKING_NECK.radius), n=1001)
        assert r["energy"] > 0.0 and r["norm"] > 0.0
        assert r["quotient"] > 0.0


def test_no_trial_configuration_beats_the_lowest_mode():
    """A variational bound has to bound.  If a trial quotient fell below the
    computed lowest mode, one of the two would be wrong."""
    m = measure_the_quadratic_form_is_positive(trials=25)
    assert m["smallest_quotient"] > m["lowest_computed_mode"]


# ── the answer ─────────────────────────────────────────────────────────────
def test_the_tube_side_is_negative_and_the_exterior_side_positive():
    """The whole answer in one line, and neither half is scanned for."""
    for a in (0.01, 0.05, 0.35):
        for d in (0.8, 1.3, 3.0):
            if a >= 0.5 * d:
                continue
            t = NeckThroat(separation=d, length=0.9, area=FOUR_PI, radius=a)
            for sigma in (1e-3, 0.5, 10.0, 200.0):
                (at, ab), (gt, gb) = t.signed_parts(sigma)
                assert at < 0.0 and ab < 0.0
                assert gt > 0.0 and gb > 0.0


def test_there_is_no_growing_mode():
    for t in (WORKING_NECK,
              NeckThroat(separation=0.8, length=0.3, area=0.2, radius=0.01),
              NeckThroat(separation=3.0, length=3.0, area=math.pi,
                         radius=0.35)):
        modes = t.growing_modes(sigma_max=300.0, n=200)
        assert modes["symmetric"] is None
        assert modes["antisymmetric"] is None
        assert modes["worst_symmetric"] < 0.0


def test_disjointness_is_enforced():
    with pytest.raises(ValueError):
        NeckThroat(separation=1.0, radius=0.6)


# ── the soft mode, and the correction to PR #261 ───────────────────────────
def test_the_symmetric_channel_runs_from_plus_to_minus_across_the_gap():
    """Existence without a scan — the counterpart of the sign argument that
    kills the growing mode."""
    for length in (0.3, 0.9, 4.0, 8.0):
        t = NeckThroat(separation=1.3, length=length, area=FOUR_PI,
                       radius=0.05)
        assert t.channel_functions(1e-9)[0] > 0.0
        assert t.channel_functions(1.0 - 1e-9)[0] < 0.0


def test_the_soft_mode_tends_to_the_capacitance_formula():
    rows = measure_the_soft_mode_survives_the_removal()["rows"]
    small = [r for r in rows if r["radius"] <= 0.02]
    assert all(abs(1.0 - r["ratio_to_closed_form"]) < 0.015 for r in small)


def test_removing_the_balls_barely_moves_the_soft_mode():
    """PR #261's approximation, priced.  ``7e-04`` at the working radius."""
    m = measure_the_soft_mode_survives_the_removal()
    assert m["shift_at_the_working_radius"] < 1e-3
    assert m["worst_shift_from_removal"] < 1e-2


def test_a_pole_is_not_a_bound_state():
    """The correction to PR #261, pinned.

    At ``L = 8`` the tube's own harmonics are inside the gap, and each is a sign
    change of the channel function with **no zero**.  Counting crossings gives
    five states; two of them are poles sitting exactly on ``(2πj/L)²`` and
    ``(π(2j−1)/L)²``.
    """
    t = NeckThroat(separation=1.3, length=8.0, area=FOUR_PI, radius=0.05)
    b = t.bound_states(n=6000)
    poles = t.tube_poles_in_the_gap()
    assert len(b["symmetric_poles"]) + len(b["antisymmetric_poles"]) == 2
    assert b["states_below_the_gap"] == 3
    for got, want in zip(b["symmetric_poles"], poles["symmetric"]):
        assert got == pytest.approx(want, rel=1e-6)
    for got, want in zip(b["antisymmetric_poles"], poles["antisymmetric"]):
        assert got == pytest.approx(want, rel=1e-6)


def test_the_root_pole_separation_is_not_a_tuned_threshold():
    """Fifteen orders of magnitude between the two residuals."""
    t = NeckThroat(separation=1.3, length=8.0, area=FOUR_PI, radius=0.05)
    b = t.bound_states(n=6000)
    for root in b["symmetric_roots"]:
        assert abs(t.channel_functions(root)[0]) < 1e-8
    for pole in b["symmetric_poles"]:
        assert abs(t.channel_functions(pole)[0]) > 1e8


def test_short_tubes_hold_exactly_one_state():
    """PR #261's claim, restated with its condition attached."""
    for length in (0.3, 0.9, 2.0, 3.0):
        t = NeckThroat(separation=1.3, length=length, area=FOUR_PI,
                       radius=0.05)
        b = t.bound_states(n=4000)
        assert b["states_below_the_gap"] == 1
        assert b["antisymmetric"] is None
        assert 0.0 < b["symmetric"] < 1.0
        assert not t.tube_poles_in_the_gap()["symmetric"]
        assert not t.tube_poles_in_the_gap()["antisymmetric"]


def test_the_first_tube_harmonic_enters_the_gap_at_pi():
    """``L = π`` is the boundary, exactly."""
    just_below = NeckThroat(separation=1.3, length=math.pi - 1e-3,
                            area=FOUR_PI, radius=0.05)
    just_above = NeckThroat(separation=1.3, length=math.pi + 1e-3,
                            area=FOUR_PI, radius=0.05)
    assert not just_below.tube_poles_in_the_gap()["antisymmetric"]
    assert just_above.tube_poles_in_the_gap()["antisymmetric"]


# ── what PR #261's fixed ambient cost ──────────────────────────────────────
def test_the_fixed_ambient_error_grows_with_the_radius():
    rows = measure_what_the_fixed_ambient_cost()["rows"]
    at_zero = [r for r in rows if r["lambda"] == 0.0]
    errors = [r["relative_error"] for r in at_zero]
    assert all(b > c for c, b in zip(errors, errors[1:]))


def test_the_two_ambient_constructions_agree_at_small_radius():
    """``f(a)G(a)`` against ``1/N₀(a)``: ``1.3e-04`` at ``a = 0.02``, ``λ = 0``."""
    for a in (0.02, 0.01):
        kept = regular_radial(a, 0.0) * mouth_green(a, 0.0)
        removed = 1.0 / exterior_dtn_monopole(a, 0.0)
        assert kept == pytest.approx(removed, rel=2e-4)


def test_the_cross_term_is_single_scattering_and_the_remainder_is_bounded():
    """The one approximation left in the reduced 2×2, priced.

    ``f(a)²G(d)`` is the single-scattering cross term; the exact two-ball
    exterior adds the series in which each boundary sphere drives the other.
    Its expansion parameter is ``cross/self``, which comes out as ``0.8·(a/d)``,
    so the neglected terms are ``~1e-03`` at the working point and cannot reach
    the sign the reduced model is being asked to decide.
    """
    m = measure_what_the_fixed_ambient_cost()
    assert m["cross_over_self_scales_like_a_over_d"]
    assert m["neglected_scattering_at_the_working_point"] < 1e-2
    assert m["the_neglected_series_cannot_reach_the_sign"]
    at_zero = [r for r in m["single_scattering"] if r["lambda"] == 0.0]
    assert all(0.7 < r["over_a_over_d"] < 1.0 for r in at_zero)


def test_the_neck_and_the_fixed_ambient_agree_on_the_soft_mode():
    neck = NeckThroat(separation=1.3, length=0.9, area=FOUR_PI, radius=0.02)
    fixed = FiniteMouthThroat(separation=1.3, length=0.9, area=FOUR_PI,
                              radius=0.02)
    assert neck.bound_states()["symmetric"] == pytest.approx(
        fixed.bound_states()["symmetric"], rel=1e-3)


# ── the measurements report what they claim ────────────────────────────────
def test_the_dtn_measurement_passes():
    m = measure_the_exterior_dtn_is_positive_in_every_channel()
    assert m["it_is_positive_everywhere"]
    assert m["it_increases_with_ell"]
    assert m["the_capacitance_normalization_holds"]
    assert m["worst_shooting_vs_closed_form"] < 1e-11


def test_the_quadratic_form_measurement_passes():
    m = measure_the_quadratic_form_is_positive(trials=15)
    assert m["every_quotient_is_positive"]
    assert m["no_trial_beats_the_lowest_mode"]
    assert m["the_interior_only_case_hits_the_poincare_floor"]


def test_the_stability_measurement_passes():
    m = measure_the_negative_mode_does_not_survive_the_neck(
        radii=(0.05, 0.35), separations=(1.3, 3.0), lengths=(0.9,),
        areas=(math.pi, FOUR_PI), sigmas=(1e-3, 1.0, 50.0))
    assert m["there_is_no_growing_mode"]
    assert m["positives"] == 0
    assert m["worst_approach_to_zero"] < 0.0


def test_the_multipole_measurement_passes():
    m = measure_the_higher_multipoles_decouple()
    assert m["every_channel_is_positive"]
    assert m["every_channel_is_stiffer_than_the_monopole"]
    assert m["smallest_ratio_to_the_monopole"] > 1.0


def test_the_fixed_ambient_measurement_passes():
    m = measure_what_the_fixed_ambient_cost()
    assert m["pr_261_was_right_where_it_looked"]
    assert m["the_error_grows_with_the_radius"]
    assert m["worst_error_overall"] > 0.1


def test_the_existence_measurement_passes():
    m = measure_the_soft_mode_is_forced_by_the_two_ends()
    assert m["the_symmetric_channel_starts_positive"]
    assert m["and_ends_negative"]
    assert m["so_a_symmetric_state_is_forced"]
    assert m["short_tubes_hold_exactly_one_state"]
    assert m["long_tubes_hold_more"]
    assert m["every_extra_state_follows_a_tube_harmonic"]
    assert m["the_exterior_stiffness_vanishes_linearly_at_the_gap_edge"]
    assert m["pole_crossings_that_would_have_been_miscounted"] > 0


def test_the_soft_mode_measurement_passes():
    m = measure_the_soft_mode_survives_the_removal()
    assert m["every_one_is_positive_and_below_the_gap"]
    assert m["the_capacitance_formula_survives"]
    assert m["removing_the_balls_barely_moves_it"]
