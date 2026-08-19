"""Tests for the finite-radius mouth, and PR #260's question."""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.finite_mouth import (
    FOUR_PI, FiniteMouthThroat, WORKING_MOUTH,
    measure_monopole_matching_is_the_remaining_approximation,
    measure_the_delay_survives_with_a_radius_correction,
    measure_the_instability_was_the_linearization,
    measure_the_mean_value_identities_hold,
    measure_the_mode_became_soft_and_positive,
    measure_the_negative_mode_does_not_survive,
    measure_the_static_results_survive,
    mouth_green, regular_radial, screened_products,
    shell_average_cross, shell_average_self,
)
from geometrodynamics.waves.throat_operator import gamma_at


# ── the ambient block ──────────────────────────────────────────────────────
def test_the_regular_solution_is_one_at_the_origin():
    for lam in (-4.0, -0.3, 0.0, 0.5, 3.0):
        assert regular_radial(1e-9, lam) == pytest.approx(1.0, abs=1e-9)


def test_the_unsubtracted_green_function_is_never_negative_below_zero():
    """Half of the answer: ``G(χ,λ) > 0`` for every ``λ < 0``.

    Asserted as ``≥ 0``, because ``e^{−κχ}`` underflows to exactly zero once
    ``κχ ≳ 745``.  That floor is harmless to the argument — the tube's side is
    *strictly* negative, so an ambient term that has underflowed to zero still
    cannot cancel it — and `test_the_difference_stays_negative_through_underflow`
    checks precisely that.
    """
    for chi in (0.02, 0.4, 1.3, 3.0):
        for sigma in (0.01, 1.0, 20.0, 300.0):
            assert mouth_green(chi, -sigma ** 2) >= 0.0
            if sigma * chi < 700.0:
                assert mouth_green(chi, -sigma ** 2) > 0.0


def test_the_green_function_matches_the_point_kernel_off_coincidence():
    """Away from the mouth it is the same object PRs #257–#260 used."""
    for lam in (-2.0, 0.0, 0.4):
        assert mouth_green(1.3, lam) == pytest.approx(
            float(gamma_at(lam, 1.3)[0, 1].real), rel=1e-12)


def test_the_screened_products_agree_with_the_direct_ones():
    for a, d, s in ((0.05, 1.3, 0.3), (0.15, 1.3, 2.0), (0.35, 2.2, 7.0)):
        lam = -s ** 2
        f = regular_radial(a, lam)
        got = screened_products(a, d, s)
        assert got[0] == pytest.approx(f * mouth_green(a, lam), rel=1e-11)
        assert got[1] == pytest.approx(f * f * mouth_green(d, lam), rel=1e-11)


def test_the_screened_products_survive_where_the_factors_overflow():
    """The factors overflow; the products do not."""
    with pytest.raises(OverflowError):
        regular_radial(0.35, -(1e4 ** 2))
    both = screened_products(0.35, 2.2, 1e4)
    assert all(math.isfinite(v) and v >= 0.0 for v in both)
    # at a merely large σ they are still strictly positive and finite
    both = screened_products(0.35, 2.2, 200.0)
    assert all(math.isfinite(v) and v > 0.0 for v in both)


def test_the_difference_stays_negative_through_underflow():
    """Where the ambient term underflows to zero the conclusion still holds,
    because the tube's side is strictly negative on its own."""
    t = FiniteMouthThroat(separation=2.2, length=3.0, area=0.2, radius=0.35)
    for sigma in (300.0, 3e3, 1e4):
        (a_sym, a_anti), (g_sym, g_anti) = t.signed_parts(sigma)
        assert a_sym < 0.0 and a_anti < 0.0
        assert g_sym >= 0.0 and g_anti >= 0.0
        assert a_sym - g_sym < 0.0 and a_anti - g_anti < 0.0


def test_the_mean_value_identity_holds_for_the_cross_term():
    for a, d, lam in ((0.05, 1.3, 0.16), (0.15, 2.2, 2.89)):
        got = shell_average_cross(a, d, lam)
        assert got == pytest.approx(regular_radial(a, lam) * mouth_green(d, lam),
                                    rel=1e-9)


def test_the_self_identity_holds_to_quadrature_accuracy():
    for a, lam in ((0.05, 0.16), (0.35, 2.89)):
        got = shell_average_self(a, lam)
        assert got == pytest.approx(regular_radial(a, lam) * mouth_green(a, lam),
                                    rel=1e-3)


# ── the throat ─────────────────────────────────────────────────────────────
def test_the_mouths_must_be_disjoint():
    with pytest.raises(ValueError):
        FiniteMouthThroat(separation=1.0, radius=0.5)
    with pytest.raises(ValueError):
        FiniteMouthThroat(separation=1.0, radius=0.0)


def test_the_tube_block_is_pr_260s_with_no_mouth_term():
    t = WORKING_MOUTH
    assert t.tube().mouth_radius is None
    assert t.tube().length == t.length and t.tube().area == t.area


def test_the_ambient_matrix_is_symmetric_and_positive_below_zero():
    for t in (WORKING_MOUTH,
              FiniteMouthThroat(separation=2.2, length=3.0, area=0.2,
                                radius=0.3)):
        for sigma in (0.05, 1.0, 12.0, 200.0):
            _, (g_sym, g_anti) = t.signed_parts(sigma)
            assert g_sym > 0.0 and g_anti > 0.0


def test_the_tube_side_is_negative_below_zero():
    for t in (WORKING_MOUTH,
              FiniteMouthThroat(separation=2.2, length=3.0, area=0.2,
                                radius=0.3)):
        for sigma in (0.05, 1.0, 12.0, 200.0):
            (a_sym, a_anti), _ = t.signed_parts(sigma)
            assert a_sym < 0.0 and a_anti < 0.0


def test_there_is_no_growing_mode():
    """**The answer.**"""
    for t in (WORKING_MOUTH,
              FiniteMouthThroat(separation=0.8, length=0.3, area=0.2,
                                radius=0.01),
              FiniteMouthThroat(separation=3.0, length=3.0, area=FOUR_PI,
                                radius=0.6, interior_mass=1.5)):
        modes = t.growing_modes()
        assert modes["symmetric"] is None
        assert modes["antisymmetric"] is None
        assert modes["worst_symmetric"] < 0.0
        assert modes["the_tube_side_is_negative"]
        assert modes["the_ambient_side_is_positive"]


def test_the_point_mouth_of_pr_260_does_have_one():
    """The contrast, run through the same object with the screening switched
    off — so the difference is the linearization and nothing else."""
    t = WORKING_MOUTH
    sig = np.geomspace(1e-2, 1e3, 400)
    exact = [t.negative_lambda_channels(float(s))[0] for s in sig]
    linear = [t.linearized_channels(float(s))[0] for s in sig]
    assert max(exact) < 0.0
    assert max(linear) > 0.0


def test_the_linearized_root_sits_at_kappa_a_of_order_one():
    r = measure_the_instability_was_the_linearization()
    for value in r["linearized_roots_times_radius"]:
        assert 0.5 < value < 2.0


def test_there_is_exactly_one_soft_bound_state():
    t = WORKING_MOUTH
    states = t.bound_states()
    assert states["antisymmetric"] is None
    assert 0.0 < states["symmetric"] < 1.0


def test_the_soft_mode_matches_the_capacitance_formula():
    t = FiniteMouthThroat(separation=1.3, length=0.9, area=FOUR_PI,
                          radius=0.005)
    lam0 = t.bound_states()["symmetric"]
    assert lam0 == pytest.approx(t.soft_mode_closed_form(), rel=5e-3)
    assert t.soft_mode_closed_form() == pytest.approx(
        8.0 * math.pi * 0.005 / (FOUR_PI * 0.9), rel=1e-12)


def test_the_soft_mode_goes_to_zero_with_the_radius():
    big = FiniteMouthThroat(radius=0.2).bound_states()["symmetric"]
    small = FiniteMouthThroat(radius=0.01).bound_states()["symmetric"]
    assert 0.0 < small < big < 1.0


def test_the_response_spectrum_inverts_the_channel_matrix():
    t = WORKING_MOUTH
    om = np.array([0.7 + 0.4j, 2.3 + 0.4j])
    r = t.response_spectrum(om)
    for i in range(2):
        assert abs(np.linalg.det(r[i])) > 0.0
    # and on the real axis it reproduces the channel functions
    lam = 1.7
    s, a = t.channel_functions(lam)
    om_real = np.array([math.sqrt(lam) + 0j])
    inv = np.linalg.inv(t.response_spectrum(om_real)[0]).real
    v_s = np.array([1.0, 1.0]) / math.sqrt(2.0)
    v_a = np.array([1.0, -1.0]) / math.sqrt(2.0)
    assert float(v_s @ inv @ v_s) == pytest.approx(s, rel=1e-9)
    assert float(v_a @ inv @ v_a) == pytest.approx(a, rel=1e-9)


def test_the_static_response_is_still_rank_one():
    s = WORKING_MOUTH.static_response(1e-8)
    assert abs(s[0, 0] + s[0, 1]) < 1e-6 * abs(s[0, 0])


# ── the measurements ───────────────────────────────────────────────────────
def test_measure_the_mean_value_identities_hold():
    r = measure_the_mean_value_identities_hold()
    assert r["the_cross_identity_is_exact"]
    assert r["the_self_identity_is_grid_limited"]


def test_measure_the_negative_mode_does_not_survive():
    r = measure_the_negative_mode_does_not_survive()
    assert r["there_is_no_growing_mode"]
    assert r["positives"] == 0
    assert r["samples"] > 1000
    assert r["worst_approach_to_zero"] < 0.0


def test_measure_the_instability_was_the_linearization():
    r = measure_the_instability_was_the_linearization()
    assert r["every_linearized_model_has_a_root"]
    assert r["no_exact_model_has_one"]
    assert r["the_root_sits_at_kappa_a_of_order_one"]
    assert r["worst_gap_below_kappa_a_of_0p1"] < 0.02
    assert r["worst_gap_at_kappa_a_of_3"] > 1.0


def test_measure_the_mode_became_soft_and_positive():
    r = measure_the_mode_became_soft_and_positive()
    assert r["every_one_is_positive"]
    assert r["every_one_is_below_the_gap"]
    assert r["no_antisymmetric_bound_state"]
    assert r["the_capacitance_formula_holds"]
    assert r["the_point_limit_approaches_zero_from_above"]


def test_measure_the_delay_survives_with_a_radius_correction():
    r = measure_the_delay_survives_with_a_radius_correction()
    assert r["the_traversal_time_survives"]
    assert r["the_ambient_path_still_takes_over"]
    assert r["the_mouth_shift_is_subleading"]
    assert r["slope_in_length"] == pytest.approx(1.0, abs=0.05)


def test_measure_the_static_results_survive():
    r = measure_the_static_results_survive()
    assert r["the_static_response_is_still_rank_one"]
    assert r["det_S_is_linear_in_lambda"]
    assert r["the_defect_is_still_minus_beta"]


def test_measure_monopole_matching_is_the_remaining_approximation():
    r = measure_monopole_matching_is_the_remaining_approximation()
    assert r["dipole_scales_like_a_over_d"]
    assert r["smallest_dropped_fraction"] < 1e-3
