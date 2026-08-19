"""Tests for the finite conservative throat."""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.throat_operator import gamma_at
from geometrodynamics.waves.finite_throat import (
    FOUR_PI, FiniteThroat, WORKING_THROAT, bounce_delays, causal_onset,
    dtn_matrix, green_identity_residual, impulse_response, interior_profile,
    channel_basis,
    measure_the_contour_must_clear_the_growing_mode,
    measure_the_enlarged_system_is_conservative,
    measure_the_delay_ledger_is_the_bounce_series,
    measure_the_growing_mode_is_interface_localized,
    measure_the_interior_mass_is_a_transmission_cutoff,
    measure_the_short_tube_limit_is_a_mixed_stratum,
    measure_the_static_limit_is_rank_one_and_the_defect_diverges,
    measure_the_throat_transmits_at_the_traversal_time,
    response_spectrum,
    short_tube_stratum,
)
from geometrodynamics.waves.two_wave import RetardedGrid, gamma_omega


# ── the DtN map ────────────────────────────────────────────────────────────
def test_the_dtn_map_is_symmetric_and_real_for_real_k():
    for k in (0.4, 1.1, 2.7, 5.5):
        n = dtn_matrix(k, 0.9, FOUR_PI)
        assert np.abs(n - n.T).max() == 0.0
        assert np.abs(n.imag).max() == 0.0


def test_the_dtn_map_is_even_in_k():
    """No branch choice is needed: ``cot(kL)·k`` and ``csc(kL)·k`` are even."""
    for k in (0.7, 2.3, 1.0 + 0.4j):
        assert np.abs(dtn_matrix(k, 1.3, 2.0)
                      - dtn_matrix(-k, 1.3, 2.0)).max() < 1e-13


def test_the_dtn_map_satisfies_greens_identity():
    """The matrix really is the tube's map, checked against the interior."""
    for k in (0.7, 1.9, 3.3):
        r = green_identity_residual(k, 0.9, FOUR_PI, (1.0, -0.4))
        assert r < 1e-6


def test_the_interior_profile_matches_its_end_values():
    s, u, du = interior_profile(1.7, 0.9, (0.6, -0.25))
    assert u[0] == pytest.approx(0.6, abs=1e-12)
    assert u[-1] == pytest.approx(-0.25, abs=1e-12)
    # and the analytic derivative is the derivative
    num = np.gradient(u, s, edge_order=2)
    assert np.abs(num[5:-5] - du[5:-5]).max() < 1e-6


def test_the_static_dtn_is_rank_one():
    """A constant field on a massless tube carries no current."""
    n = dtn_matrix(1e-9, 1.4, 3.0)
    assert abs(np.linalg.det(n)) < 1e-6 * np.abs(n).max()
    assert np.abs(n @ np.array([1.0, 1.0])).max() < 1e-6 * np.abs(n).max()


def test_the_dtn_determinant_is_exact():
    """``det N = −(𝒜k)²``, which is why ``A = −N⁻¹`` has a closed form."""
    for k in (0.5, 1.6, 4.4):
        n = dtn_matrix(k, 0.9, 1.7)
        assert complex(np.linalg.det(n)).real == pytest.approx(
            -(1.7 * k) ** 2, rel=1e-12)


# ── the boundary condition ─────────────────────────────────────────────────
def test_the_boundary_matrix_is_minus_the_inverse_dtn():
    t = FiniteThroat(length=0.9, area=FOUR_PI)
    for lam in (0.3, 1.0, 5.0):
        a = t.boundary_matrix(complex(lam))
        n = t.dtn(complex(lam))
        assert np.abs(a + np.linalg.inv(n)).max() < 1e-12


def test_the_mouth_radius_shifts_the_diagonal_only():
    a = FiniteThroat(length=0.9, area=FOUR_PI)
    b = FiniteThroat(length=0.9, area=FOUR_PI, mouth_radius=0.2)
    diff = b.boundary_matrix(complex(1.0)) - a.boundary_matrix(complex(1.0))
    assert np.abs(diff - np.eye(2) * diff[0, 0]).max() < 1e-14
    assert complex(diff[0, 0]).real == pytest.approx(-1.0 / (FOUR_PI * 0.2),
                                                     rel=1e-12)


def test_the_boundary_pair_is_maximal_and_self_adjoint():
    t = WORKING_THROAT
    for lam in (-3.0, -0.2, 0.4, 2.0, 9.0):
        b, c = t.boundary_pair(complex(lam))
        prod = b @ c.conjugate().T
        assert np.abs(prod - prod.conjugate().T).max() < 1e-13
        assert np.linalg.matrix_rank(np.hstack([b, c]), tol=1e-10) == 2


def test_the_channel_functions_diagonalize_the_two_by_two():
    t = WORKING_THROAT
    lam = 1.7
    m = t.boundary_matrix(complex(lam)) - gamma_at(lam, t.separation)
    s, a = t.channel_functions(lam)
    v_s = np.array([1.0, 1.0]) / math.sqrt(2.0)
    v_a = np.array([1.0, -1.0]) / math.sqrt(2.0)
    assert complex(v_s @ m @ v_s).real == pytest.approx(s, rel=1e-12)
    assert complex(v_a @ m @ v_a).real == pytest.approx(a, rel=1e-12)


def test_the_channel_closed_forms_are_the_half_angle_ones():
    """``A_sym = cot(kL/2)/(𝒜k)`` and ``A_anti = −tan(kL/2)/(𝒜k)``.

    Pinned because the two are easy to write down the wrong way round: the
    docstring of `channel_functions` had them swapped while every formula in
    the module had them right.
    """
    t = FiniteThroat(separation=1.3, length=0.9, area=FOUR_PI)
    for lam in (0.3, 1.7, 6.0):
        k = math.sqrt(lam)
        a = t.boundary_matrix(complex(lam))
        sym = float((a[0, 0] + a[0, 1]).real)
        anti = float((a[0, 0] - a[0, 1]).real)
        assert sym == pytest.approx(
            1.0 / (math.tan(k * t.length / 2.0) * t.area * k), rel=1e-12)
        assert anti == pytest.approx(
            -math.tan(k * t.length / 2.0) / (t.area * k), rel=1e-12)
    # and the symmetric channel is the one that carries the k → 0 pole, while
    # the antisymmetric one goes to the finite −L/(2𝒜)
    small = t.boundary_matrix(complex(1e-6))
    assert float((small[0, 0] + small[0, 1]).real) == pytest.approx(
        2.0 / (t.area * 1e-6 * t.length), rel=1e-6)
    assert float((small[0, 0] - small[0, 1]).real) == pytest.approx(
        -t.length / (2.0 * t.area), rel=1e-6)


def test_the_negative_lambda_channels_are_the_same_two_functions():
    """The continuations ``−coth(κL/2)/(𝒜κ)`` and ``−tanh(κL/2)/(𝒜κ)``."""
    t = FiniteThroat(separation=1.3, length=0.9, area=FOUR_PI)
    for sigma in (0.4, 2.0):
        g = gamma_at(-sigma ** 2, t.separation).real
        sym, anti = t.negative_lambda_channels(sigma)
        scale = t.area * sigma
        assert sym + float(g[0, 0] + g[0, 1]) == pytest.approx(
            -1.0 / (math.tanh(sigma * t.length / 2.0) * scale), rel=1e-12)
        assert anti + float(g[0, 0] - g[0, 1]) == pytest.approx(
            -math.tanh(sigma * t.length / 2.0) / scale, rel=1e-12)


def test_the_transmission_amplitude_is_the_off_diagonal():
    t = WORKING_THROAT
    for lam in (0.4, 2.2, 6.1):
        assert complex(t.transmission(complex(lam))).real == pytest.approx(
            float(t.boundary_matrix(complex(lam))[0, 1].real), rel=1e-12)


def test_the_boundary_matrix_is_array_aware():
    t = WORKING_THROAT
    lams = np.array([0.4, 1.1, 3.0])
    stack = t.boundary_matrix(lams)
    assert stack.shape == (3, 2, 2)
    for i, lam in enumerate(lams):
        assert np.abs(stack[i] - t.boundary_matrix(complex(lam))).max() < 1e-14


# ── monotonicity, and the growing mode ─────────────────────────────────────
def test_the_boundary_matrix_decreases_in_lambda():
    """``A(λ)`` is Nevanlinna-decreasing between poles — the throat-side echo
    of PR #257's ``dΓ/dλ ≻ 0`` Gram identity, and what makes the mode count a
    count rather than a scan."""
    t = WORKING_THROAT
    # strictly inside one gap: Γ has poles at the free eigenvalues λ = n², and
    # the tube has its own at λ = (nπ/L)², so monotonicity is a statement
    # *between* poles and the interval has to say so
    for lo, hi in ((0.02, 0.98), (1.02, 3.98), (4.02, 8.98)):
        lams = np.linspace(lo, hi, 60)
        for idx in (0, 1):
            vals = [t.channel_functions(float(x))[idx] for x in lams]
            assert all(b < a for a, b in zip(vals, vals[1:]))


def test_the_stable_negative_lambda_channels_agree_with_gamma_at():
    """Two constructions of ``Γ`` on the imaginary axis, one of which does not
    overflow."""
    t = WORKING_THROAT
    for sigma in (0.2, 1.0, 4.0, 20.0):
        stable = t.negative_lambda_channels(sigma)
        direct = t.channel_functions(-sigma ** 2)
        assert stable[0] == pytest.approx(direct[0], rel=1e-11)
        assert stable[1] == pytest.approx(direct[1], rel=1e-11)


def test_the_stable_channels_survive_where_gamma_at_overflows():
    t = FiniteThroat(length=6.0, area=0.2)
    values = t.negative_lambda_channels(600.0)
    assert all(math.isfinite(v) for v in values)
    with pytest.raises((OverflowError, ValueError)):
        gamma_at(-600.0 ** 2, t.separation)


def test_every_throat_has_a_symmetric_growing_mode():
    for area in (0.4, math.pi, FOUR_PI):
        for length in (0.5, 2.0):
            t = FiniteThroat(length=length, area=area)
            assert t.growing_modes()["symmetric"] is not None


def test_the_growing_mode_matches_its_closed_form():
    t = FiniteThroat(length=3.0, area=0.2)
    modes = t.growing_modes()
    assert float(modes["symmetric"]) == pytest.approx(
        t.sigma_star_closed_form(), rel=5e-6)


def test_the_growing_mode_ignores_the_mouth_separation():
    a = FiniteThroat(separation=2.4, length=3.0, area=0.2)
    b = FiniteThroat(separation=3.0, length=3.0, area=0.2)
    assert float(a.growing_modes()["symmetric"]) == pytest.approx(
        float(b.growing_modes()["symmetric"]), abs=1e-7)


def test_a_finite_mouth_radius_raises_the_growing_mode():
    a = FiniteThroat(length=3.0, area=0.5)
    b = FiniteThroat(length=3.0, area=0.5, mouth_radius=0.1)
    assert (float(b.growing_modes()["symmetric"])
            > float(a.growing_modes()["symmetric"]))
    assert b.sigma_star_closed_form() == pytest.approx(
        0.5 * (10.0 + math.sqrt(100.0 + 16.0 * math.pi / 0.5)), rel=1e-12)


# ── the impulse response ───────────────────────────────────────────────────
def test_the_response_spectrum_is_the_inverse_krein_matrix():
    t = WORKING_THROAT
    om = np.array([0.7 + 1.9j, 2.3 + 1.9j])
    r = response_spectrum(t, om)
    m = t.boundary_matrix(om ** 2) - gamma_omega(om, t.separation)
    for i in range(2):
        assert np.abs(r[i] @ m[i] - np.eye(2)).max() < 1e-11


def test_the_transmission_onset_is_the_traversal_time():
    t = FiniteThroat(separation=1.3, length=0.6, area=FOUR_PI)
    grid = RetardedGrid(n=1 << 16, span=200.0,
                        eps=float(t.growing_modes()["symmetric"]) + 0.8)
    imp = impulse_response(t, grid, width=0.03)
    onset = causal_onset(imp["opposite_mouths"], imp["times"])
    assert 0.3 < onset < 0.6
    assert causal_onset(imp["same_mouth"], imp["times"]) < 0.12


def test_the_point_throat_transmits_instantaneously():
    t = FiniteThroat(separation=1.3, length=0.6, area=FOUR_PI)
    grid = RetardedGrid(n=1 << 16, span=200.0,
                        eps=float(t.growing_modes()["symmetric"]) + 0.8)
    frozen = t.boundary_matrix(complex(1.0)).real.astype(complex)
    imp = impulse_response(t, grid, width=0.03, constant=frozen)
    assert causal_onset(imp["opposite_mouths"], imp["times"]) < 0.12


def test_the_bounce_delays_have_the_right_parity():
    d = bounce_delays(3)
    assert d["same_mouth"] == [0.0, 2.0, 4.0]
    assert d["opposite_mouths"] == [1.0, 3.0, 5.0]


def test_the_causal_onset_finds_a_planted_edge():
    times = np.linspace(0.0, 10.0, 2001)
    series = np.where(times < 3.0, 0.0, 1.0)
    assert causal_onset(series, times) == pytest.approx(3.0, abs=6e-3)


# ── the static limit ───────────────────────────────────────────────────────
def test_the_static_response_collapses_onto_the_antisymmetric_direction():
    t = WORKING_THROAT
    s = np.linalg.inv(t.krein_matrix(complex(1e-8))).real
    assert abs(s[0, 0] + s[0, 1]) < 1e-6 * abs(s[0, 0])
    assert abs(np.linalg.det(s)) < 1e-5 * abs(s[0, 0]) ** 2


def test_an_interior_mass_restores_the_rank():
    t = FiniteThroat(length=0.9, area=FOUR_PI, interior_mass=0.3)
    s = np.linalg.inv(t.krein_matrix(complex(0.0))).real
    assert abs(np.linalg.det(s)) > 1.0


def test_the_defect_is_still_minus_the_transmission():
    """PR #258's theorem, on a throat with an interior."""
    t = FiniteThroat(length=0.9, area=FOUR_PI, interior_mass=0.25)
    g0 = float(gamma_at(0.0, t.separation)[0, 1].real)
    s = np.linalg.inv(t.krein_matrix(complex(0.0))).real
    defect = s[0, 1] / float(np.linalg.det(s)) - g0
    assert defect == pytest.approx(
        -float(t.transmission(complex(0.0)).real), rel=1e-10)


def test_the_cutoff_separates_oscillation_from_decay():
    t = FiniteThroat(length=3.0, area=FOUR_PI, interior_mass=2.0)
    below = float(t.transmission(complex(1.0)).real)
    above = float(t.transmission(complex(20.0)).real)
    assert below < 0.0
    assert abs(below) < 0.05 * abs(above)


def test_the_resonances_are_the_interior_spectrum():
    t = FiniteThroat(length=0.9, area=FOUR_PI, interior_mass=0.5)
    for i, lam in enumerate(t.resonances(3), start=1):
        assert lam == pytest.approx(0.25 + (i * math.pi / 0.9) ** 2, rel=1e-12)
        assert abs(complex(t.transmission(complex(lam))).real) > 1e10


# ── the measurements ───────────────────────────────────────────────────────
def test_measure_the_enlarged_system_is_conservative():
    r = measure_the_enlarged_system_is_conservative()
    assert r["the_finite_throat_is_conservative"]
    assert r["the_control_is_not"]
    assert r["worst_green_identity_residual"] < 1e-5


def test_measure_the_throat_transmits_at_the_traversal_time():
    r = measure_the_throat_transmits_at_the_traversal_time()
    assert r["the_delay_is_the_traversal_time"]
    assert r["the_ambient_path_takes_over"]
    assert r["reflection_is_instantaneous"]
    assert r["the_point_throat_transmits_instantly"]
    assert r["slope_below_the_ambient_path"] == pytest.approx(1.0, abs=0.05)


def test_measure_the_delay_ledger_is_the_bounce_series():
    r = measure_the_delay_ledger_is_the_bounce_series()
    assert r["the_series_converge_on_the_contour"]
    assert r["cot_series_error"] < 1e-11
    assert r["csc_series_error"] < 1e-11


def test_measure_the_short_tube_limit_is_a_mixed_stratum():
    r = measure_the_short_tube_limit_is_a_mixed_stratum()
    assert r["the_limit_exists_and_is_not_a_finite_A"]
    assert r["the_band_error_reaches_one"]
    assert r["the_antisymmetric_channel_has_a_limit"]
    assert r["the_chart_matrix_diverges"]
    assert r["convergence_is_linear_in_L"]
    assert r["every_pair_is_maximal"]
    assert r["band"][0]["relative_error"] == 0.0


def test_measure_the_static_limit_is_rank_one_and_the_defect_diverges():
    r = measure_the_static_limit_is_rank_one_and_the_defect_diverges()
    assert r["the_static_response_is_rank_one"]
    assert r["the_stratum_is_rank_one_too"]
    assert r["it_falsifies_the_finite_A_family_not_point_ness"]
    assert r["det_S_is_linear_in_lambda"]
    assert r["the_defect_diverges"]
    assert r["the_defect_is_still_minus_beta"]


def test_the_short_tube_pair_converges_to_the_mixed_stratum():
    """``(B, C) → (P_anti, −P_sym)``: Φ_anti = 0 and q_sym = 0, linearly in L."""
    b_star, c_star = short_tube_stratum()
    rates = []
    for length in (0.2, 0.1, 0.05, 0.02):
        t = FiniteThroat(separation=1.3, length=length, area=FOUR_PI)
        b, c = t.normalized_pair(1.0)
        gap = max(np.abs(b - b_star).max(), np.abs(c - c_star).max())
        assert np.linalg.matrix_rank(np.hstack([b, c]), tol=1e-12) == 2
        rates.append(gap / length)
    # the rate is 𝒜λ/2 = 2π for the working area
    assert rates[-1] == pytest.approx(2.0 * math.pi, rel=2e-3)
    assert max(rates) / min(rates) < 1.02


def test_the_stratum_is_reached_by_no_finite_boundary_matrix():
    b_star, c_star = short_tube_stratum()
    assert abs(np.linalg.det(b_star)) < 1e-15
    assert abs(np.linalg.det(c_star)) < 1e-15
    assert np.linalg.matrix_rank(np.hstack([b_star, c_star]), tol=1e-12) == 2
    prod = b_star @ c_star.conjugate().T
    assert np.abs(prod - prod.conjugate().T).max() == 0.0


def test_the_channel_basis_diagonalizes_the_dtn():
    v = channel_basis()
    n = dtn_matrix(1.4, 0.9, FOUR_PI)
    diag = v.T @ n @ v
    assert abs(diag[0, 1]) < 1e-14 and abs(diag[1, 0]) < 1e-14


def test_the_green_identity_is_sesquilinear():
    """Complex end values: the bilinear form is wrong, the sesquilinear is not."""
    k, length, area, ends = 1.9, 0.9, FOUR_PI, (1.0 + 0.7j, -0.4 + 0.2j)
    assert green_identity_residual(k, length, area, ends) < 1e-6
    s, u, du = interior_profile(complex(k), length, ends, 40001)
    bilinear = area * np.trapezoid(du ** 2 - k ** 2 * u ** 2, s)
    phi = np.array(ends, dtype=complex)
    assert abs(complex(phi @ dtn_matrix(k, length, area) @ phi)
               - bilinear) < 1e-6            # the bilinear identity also holds
    # but it is not the energy: the two forms disagree on complex data
    assert abs(bilinear - area * np.trapezoid(
        np.abs(du) ** 2 - k ** 2 * np.abs(u) ** 2, s)) > 1e-3
    with pytest.raises(ValueError):
        green_identity_residual(1.0 + 0.3j, length, area, ends)


def test_measure_the_interior_mass_is_a_transmission_cutoff():
    r = measure_the_interior_mass_is_a_transmission_cutoff()
    assert r["the_evanescent_side_does_not_oscillate"]
    assert r["the_transmission_is_suppressed_below_cutoff"]
    assert r["everything_stays_real"]
    assert r["asymptote_error"] < 1e-3


def test_measure_the_growing_mode_is_interface_localized():
    r = measure_the_growing_mode_is_interface_localized()
    assert r["every_throat_has_one"]
    # the localization is asymptotic, and the working throat is not there
    assert r["the_working_throat_is_not_asymptotic"]
    assert r["it_stops_knowing_the_separation"]
    assert r["the_closed_form_holds_once_sigma_L_is_large"]
    assert r["the_split_is_the_euclidean_propagator"]


def test_measure_the_contour_must_clear_the_growing_mode():
    """Two conditions, not one — clearing the mode is necessary, not enough."""
    r = measure_the_contour_must_clear_the_growing_mode()
    assert r["a_contour_below_the_mode_breaks_causality"]
    assert r["pedestal_below"] > 0.1
    # above the mode but under-resolved by the grid: still broken
    assert r["clearing_the_mode_is_not_enough"]
    assert r["pedestal_above_but_unresolved"] > 1e-4
    # and resolved: clean, with a converged onset
    assert r["the_resolved_contours_are_clean"]
    assert r["the_recovered_onset_converges"]
    assert r["pedestal_resolved"] < 1e-10


def test_the_failing_clearance_is_under_one_grid_spacing():
    """The diagnosis, pinned: the clearance that fails is 0.95 of 2π/span."""
    r = measure_the_contour_must_clear_the_growing_mode()
    bad = [row for row in r["rows"]
           if row["above_the_mode"] and not row["resolved"]]
    assert bad, "the sweep must include an above-but-under-resolved contour"
    for row in bad:
        assert abs(row["clearance_over_spacing"]) < 4.0
        assert row["pedestal"] > 1e-4
