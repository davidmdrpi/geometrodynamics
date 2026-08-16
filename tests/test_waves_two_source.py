"""Tests for the two-source invariant of a point-supported throat."""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.throat_operator import MouthPair, gamma_at
from geometrodynamics.waves.throat_positivity import (
    boundary_pair_from_unitary, positivity_defect)
from geometrodynamics.waves.two_source import (
    WORKING_BOUNDARY, WORKING_SEPARATION, cross_matrix, defect_of_pair,
    disconnection_defect, energy_functional, free_interaction_energy, geodesic,
    green_at, interaction_energy, invisible_partner, is_real_field_compatible,
    isotropy_profile, measure_a_real_field_forces_beta_real,
    measure_anisotropy_is_not_the_signature,
    measure_phase_sensitive_sources_need_only_one_spectral_parameter,
    measure_the_antipodal_endpoint_on_its_own,
    measure_the_blind_spot_of_a_single_frequency_test,
    measure_the_defect_is_the_mouth_mixing_amplitude,
    measure_the_invariant_is_recoverable_from_observations,
    measure_the_invariant_vanishes_when_a_source_is_removed,
    measure_the_throat_channel_has_the_rank_of_the_boundary_condition,
    measure_two_disconnected_scatterers_lie_on_a_surface,
    measure_two_spectral_parameters_reconstruct_the_boundary_matrix,
    mouth_channel_invariant, mouth_positions, random_points, recover_boundary,
    recover_complex_response, recover_response, response_matrix,
    response_of_pair, ring_points, source_vector, static_response)


def _working() -> MouthPair:
    a1, a2, b = WORKING_BOUNDARY
    return MouthPair(WORKING_SEPARATION, a1, a2, b)


# ── the background ──────────────────────────────────────────────────────────
def test_the_green_function_is_the_gamma_off_diagonal():
    """One implementation, read at two kinds of separation."""
    for chi in (0.3, 1.1, 2.4, math.pi):
        for lam in (0.0, -1.0, -9.0, 0.5):
            assert green_at(chi, lam) == pytest.approx(
                float(gamma_at(lam, chi)[0, 1].real), abs=0.0)


def test_the_static_green_function_is_the_coulomb_tail_plus_g0():
    """``G₀(χ) = 1/(4πχ) + g₀ + O(χ)`` — the subtraction the throat lives on."""
    g0 = -1.0 / (4.0 * math.pi ** 2)
    for chi in (1e-2, 3e-3, 1e-3):
        residual = green_at(chi) - 1.0 / (4.0 * math.pi * chi)
        assert residual == pytest.approx(g0, rel=1e-2)
        # and the next term is the one the expansion predicts.  Not pushed
        # below χ ~ 1e-3: the closed form is written in e = π − χ so that the
        # *antipode* is accurate, which costs relative precision in sin e at
        # the opposite corner — a deliberate trade, and the coincidence limit
        # is the excluded point anyway.
        assert residual - g0 == pytest.approx(chi / (24.0 * math.pi), rel=2e-2)


def test_the_mouths_sit_at_the_requested_separation():
    for d in (0.2, 1.3, 3.0, math.pi):
        c1, c2 = mouth_positions(d)
        assert geodesic(c1, c2) == pytest.approx(d, abs=1e-12)


def test_a_ring_is_at_constant_distance_from_its_centre():
    centre = random_points(1, seed=5)[0]
    ring = ring_points(centre, 1.0, 40, seed=6)
    for p in ring:
        assert np.linalg.norm(p) == pytest.approx(1.0, abs=1e-12)
        assert geodesic(centre, p) == pytest.approx(1.0, abs=1e-12)


# ── the invariant ───────────────────────────────────────────────────────────
def test_the_invariant_is_the_krein_resolvent():
    """The branch decomposition sums to ``Re G_A(y_A, y_B)``, independently."""
    pair = _working()
    pts = random_points(6, seed=11)
    r = np.linalg.inv(pair.krein_matrix(0.0))
    for k in range(3):
        y_a, y_b = pts[2 * k], pts[2 * k + 1]
        v_a = source_vector(y_a, pair.separation)
        v_b = source_vector(y_b, pair.separation)
        direct = free_interaction_energy(y_a, y_b)
        expect = float(direct + (v_a @ r @ v_b).real)
        parts = mouth_channel_invariant(pair, y_a, y_b)
        assert parts["total"] == pytest.approx(expect, rel=1e-12)
        assert interaction_energy(pair, y_a, y_b) == pytest.approx(expect,
                                                                   rel=1e-12)


def test_the_invariant_is_bilinear_in_the_two_sources():
    """Scaling either source scales it, so removing one zeroes it."""
    pair = _working()
    pts = random_points(2, seed=12)
    base = interaction_energy(pair, pts[0], pts[1]) - free_interaction_energy(
        pts[0], pts[1])
    v_a = source_vector(pts[0], pair.separation)
    v_b = source_vector(pts[1], pair.separation)
    s = static_response(pair)
    assert (3.0 * v_a) @ s @ (2.0 * v_b) == pytest.approx(6.0 * base, rel=1e-12)
    assert (0.0 * v_a) @ s @ v_b == 0.0


def test_the_free_interaction_depends_only_on_the_separation():
    centre = random_points(1, seed=13)[0]
    ring = ring_points(centre, 0.9, 60, seed=14)
    vals = [free_interaction_energy(centre, p) for p in ring]
    assert max(vals) - min(vals) < 1e-14


def test_the_throat_breaks_that_and_so_does_a_disconnected_pair():
    pair = _working()
    prof = isotropy_profile(pair, chi=1.0, n=60)
    assert prof["the_free_field_is_isotropic"]
    assert prof["throat_spread"] > 1e-3
    assert prof["disconnected_spread"] > 1e-3


# ── the rank statement ──────────────────────────────────────────────────────
def test_the_throat_table_is_rank_two_at_any_source_count():
    pair = _working()
    for n in (4, 8, 16):
        out = cross_matrix(pair, random_points(n, seed=17))
        assert out["throat_rank"] == 2
        assert out["direct_rank"] == n


def test_the_response_is_hermitian_for_every_self_adjoint_pair():
    rng = np.random.default_rng(19)
    gam = gamma_at(0.0, WORKING_SEPARATION)
    for _ in range(12):
        x = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        q, r = np.linalg.qr(x)
        u = q @ np.diag(np.diag(r) / np.abs(np.diag(r)))
        b, c = boundary_pair_from_unitary(u)
        resp = response_matrix(b, c, gam)
        assert np.abs(resp - resp.conjugate().T).max() < 1e-11


def test_the_chart_response_is_the_general_one():
    pair = _working()
    cond = pair.boundary_condition()
    general = response_matrix(cond.B, cond.C, pair.gamma(0.0))
    assert np.abs(general - response_of_pair(pair)).max() < 1e-12


def test_a_real_dirichlet_direction_drops_the_static_table_to_rank_one():
    """A rank-one boundary condition is rank one to static sources only when
    its Dirichlet direction is real; a complex one still fills both channels."""
    gam = gamma_at(0.0, WORKING_SEPARATION)
    pts = random_points(10, seed=21)
    vv = np.array([source_vector(p, WORKING_SEPARATION) for p in pts]).T

    def table_rank(direction: np.ndarray) -> int:
        v = direction / np.linalg.norm(direction)
        w = np.array([-v[1].conjugate(), v[0].conjugate()])
        u = np.outer(v, v.conjugate()) + np.exp(1.1j) * np.outer(
            w, w.conjugate())
        b, c = boundary_pair_from_unitary(u)
        assert abs(np.linalg.det(b)) < 1e-12
        resp = response_matrix(b, c, gam)
        assert np.linalg.matrix_rank(resp, tol=1e-10) == 1
        return int(np.linalg.matrix_rank(vv.T @ resp.real @ vv, tol=1e-12))

    assert table_rank(np.array([0.6, -0.8], dtype=complex)) == 1
    assert table_rank(np.array([0.6, -0.8 + 0.5j])) == 2


# ── the discriminator ───────────────────────────────────────────────────────
def test_two_disconnected_scatterers_satisfy_the_surface_equation():
    rng = np.random.default_rng(23)
    for _ in range(40):
        d = float(rng.uniform(0.25, math.pi))
        a1, a2 = rng.uniform(0.12, 0.7, size=2)
        assert abs(defect_of_pair(MouthPair(d, a1, a2, 0.0))) < 1e-12


def test_the_defect_is_minus_the_real_mouth_mixing_amplitude():
    rng = np.random.default_rng(29)
    for _ in range(40):
        d = float(rng.uniform(0.25, math.pi))
        a1, a2 = rng.uniform(0.12, 0.7, size=2)
        b = float(rng.uniform(-0.3, 0.3))
        assert defect_of_pair(MouthPair(d, a1, a2, b)) == pytest.approx(
            -b, abs=1e-12)


def test_the_defect_does_not_move_with_the_loewner_margin():
    """The raw invariant diverges at the cone's boundary; ``𝒲`` does not."""
    gam = gamma_at(0.0, WORKING_SEPARATION).real
    beta, pts, vals, defects = 0.06, random_points(2, seed=31), [], []
    for eps in (0.4, 0.02, 0.004):
        a = float(gam[0, 0] + abs(beta - gam[0, 1]) + eps)
        pair = MouthPair(WORKING_SEPARATION, a, a, beta)
        margin = positivity_defect(pair.boundary_matrix(),
                                   pair.separation)["min_eigenvalue"]
        assert margin == pytest.approx(eps, abs=1e-12)
        vals.append(abs(interaction_energy(pair, pts[0], pts[1])))
        defects.append(defect_of_pair(pair))
    assert vals[-1] > 3.0 * vals[0]
    assert max(abs(w + beta) for w in defects) < 1e-12


def test_the_defect_is_undefined_when_the_response_is_singular():
    with pytest.raises(ValueError):
        disconnection_defect(np.zeros((2, 2)), 0.05)


# ── the blind spot, and the repair ──────────────────────────────────────────
def test_a_connected_throat_can_be_invisible_at_one_frequency():
    a1, a2, rb = 0.30, 0.35, -0.05
    ib = invisible_partner(a1, a2, rb, WORKING_SEPARATION)
    assert ib is not None and ib > 0.1
    pair = MouthPair(WORKING_SEPARATION, a1, a2, complex(rb, ib))
    assert abs(defect_of_pair(pair)) < 1e-12
    assert abs(complex(pair.beta)) > 0.2                # not a small coupling


def test_the_blind_family_is_inside_the_stable_cone():
    a1, a2, rb = 0.30, 0.35, -0.05
    ib = invisible_partner(a1, a2, rb, WORKING_SEPARATION)
    pair = MouthPair(WORKING_SEPARATION, a1, a2, complex(rb, ib))
    margin = positivity_defect(pair.boundary_matrix(),
                               pair.separation)["min_eigenvalue"]
    assert margin > 0.05


def test_the_blind_family_is_visible_at_a_second_spectral_parameter():
    a1, a2, rb = 0.30, 0.35, -0.05
    ib = invisible_partner(a1, a2, rb, WORKING_SEPARATION)
    pair = MouthPair(WORKING_SEPARATION, a1, a2, complex(rb, ib))
    assert abs(defect_of_pair(pair, 0.8)) > 1e-3
    assert abs(defect_of_pair(pair, 0.3)) > 1e-3


def test_invisible_partner_returns_none_where_there_is_no_root():
    """Between ``0`` and ``G_d`` both factors have the same sign, so the
    invisibility equation has no real root."""
    gd = float(gamma_at(0.0, WORKING_SEPARATION).real[0, 1])
    assert invisible_partner(0.3, 0.35, 0.5 * gd, WORKING_SEPARATION) is None


def test_the_other_blind_branch_is_excluded_by_the_stability_gate():
    """``Re β > G_d`` is blind too, and unstable — PR #257 removes it."""
    gd = float(gamma_at(0.0, WORKING_SEPARATION).real[0, 1])
    rb = gd + 0.02
    ib = invisible_partner(0.3, 0.35, rb, WORKING_SEPARATION)
    assert ib is not None
    pair = MouthPair(WORKING_SEPARATION, 0.3, 0.35, complex(rb, ib))
    assert abs(defect_of_pair(pair)) < 1e-10
    margin = positivity_defect(pair.boundary_matrix(),
                               pair.separation)["min_eigenvalue"]
    assert margin < 0.0


def test_two_spectral_parameters_reconstruct_the_boundary_matrix():
    rng = np.random.default_rng(37)
    for _ in range(3):
        a1, a2 = rng.uniform(0.2, 0.5, size=2)
        rb, ib = rng.uniform(-0.15, 0.15), rng.uniform(-0.25, 0.25)
        pair = MouthPair(WORKING_SEPARATION, float(a1), float(a2),
                         complex(rb, ib))
        out = recover_boundary(pair)
        assert out["max_parameter_error"] < 1e-6
        assert out["residual"] < 1e-9


# ── the protocol ────────────────────────────────────────────────────────────
def test_the_response_is_recoverable_from_measured_invariants():
    pair = _working()
    pts = random_points(20, seed=41)
    obs = [(pts[2 * k], pts[2 * k + 1],
            interaction_energy(pair, pts[2 * k], pts[2 * k + 1]))
           for k in range(10)]
    rec = recover_response(obs, pair.separation)
    assert np.abs(rec["response"] - static_response(pair)).max() < 1e-10
    gam = pair.gamma(0.0).real
    assert disconnection_defect(rec["response"],
                                float(gam[0, 1])) == pytest.approx(
        -complex(pair.beta).real, abs=1e-9)


def test_three_placements_are_enough():
    pair = _working()
    pts = random_points(6, seed=43)
    obs = [(pts[2 * k], pts[2 * k + 1],
            interaction_energy(pair, pts[2 * k], pts[2 * k + 1]))
           for k in range(3)]
    rec = recover_response(obs, pair.separation)
    assert rec["n_observations"] == 3
    assert np.abs(rec["response"] - static_response(pair)).max() < 1e-8


# ── the antipodal endpoint, tested as itself ────────────────────────────────
def test_at_the_antipode_the_invariant_diverges_as_the_boundary_data_vanishes():
    pts = random_points(2, seed=47)
    prev = None
    for eps in (1e-3, 1e-4, 1e-5):
        pair = MouthPair(math.pi, eps, eps, 0.0)
        val = abs(interaction_energy(pair, pts[0], pts[1]) * eps)
        if prev is not None:
            assert val == pytest.approx(prev, rel=2e-2)
        prev = val


def test_at_the_antipode_the_defect_stays_exactly_zero():
    for eps in (1e-2, 1e-4, 1e-6):
        assert abs(defect_of_pair(MouthPair(math.pi, eps, eps, 0.0))) < 1e-12


def test_the_identity_survives_the_antipodal_endpoint():
    for beta in (-0.2, 0.06, 0.25):
        pair = MouthPair(math.pi, 0.30, 0.35, beta)
        assert defect_of_pair(pair) == pytest.approx(-beta, abs=1e-12)


def test_the_working_point_is_strictly_inside_the_cone():
    pair = _working()
    out = positivity_defect(pair.boundary_matrix(), pair.separation)
    assert out["non_negative"]
    assert out["min_eigenvalue"] > 0.3


# ── the measurements ────────────────────────────────────────────────────────
def test_measure_the_invariant_vanishes_when_a_source_is_removed():
    r = measure_the_invariant_vanishes_when_a_source_is_removed()
    assert r["it_vanishes_exactly"]
    assert r["it_is_not_vacuous"]
    assert r["is_bilinear"]
    assert r["inside_the_cone"]


def test_measure_the_throat_channel_has_the_rank_of_the_boundary_condition():
    r = measure_the_throat_channel_has_the_rank_of_the_boundary_condition()
    assert r["the_throat_table_is_rank_two_in_the_chart"]
    assert r["the_direct_table_has_full_rank"]
    assert r["the_complex_response_has_the_rank_of_B"]
    assert r["a_complex_dirichlet_direction_still_fills_both_channels"]
    assert r["a_real_dirichlet_direction_drops_the_table_to_one"]
    assert r["the_bound_is_two_at_every_source_count"]


def test_measure_anisotropy_is_not_the_signature():
    r = measure_anisotropy_is_not_the_signature()
    assert r["the_free_field_is_isotropic"]
    assert r["both_break_it"]
    assert r["anisotropy_does_not_discriminate"]


def test_measure_two_disconnected_scatterers_lie_on_a_surface():
    r = measure_two_disconnected_scatterers_lie_on_a_surface()
    assert r["the_disconnected_family_is_a_surface"]
    assert r["connected_throats_are_off_it"]


def test_measure_the_defect_is_the_mouth_mixing_amplitude():
    r = measure_the_defect_is_the_mouth_mixing_amplitude()
    assert r["W_is_minus_beta"]
    assert r["every_row_is_inside_the_cone"]
    assert r["the_discriminator_does_not_see_the_resonance"]


def test_measure_the_invariant_is_recoverable_from_observations():
    r = measure_the_invariant_is_recoverable_from_observations()
    assert r["the_protocol_closes"]
    assert r["W_error"] < 1e-9


def test_measure_the_blind_spot_of_a_single_frequency_test():
    r = measure_the_blind_spot_of_a_single_frequency_test()
    assert r["the_blind_family_is_not_empty"]
    assert r["the_upper_branch_is_excluded_by_the_stability_gate"]
    assert r["the_lower_branch_survives_the_stability_gate"]
    assert r["but_no_blind_point_is_real_field_compatible"]
    assert r["they_are_invisible_at_lambda_zero"]
    assert r["they_are_visible_at_a_second_spectral_parameter"]
    assert r["largest_stable_invisible_coupling"] > 0.2
    assert r["every_stable_coupling_is_smaller_than_its_self_energies"]


def test_measure_two_spectral_parameters_reconstruct_the_boundary_matrix():
    r = measure_two_spectral_parameters_reconstruct_the_boundary_matrix()
    assert r["the_boundary_matrix_is_reconstructed"]
    assert r["even_the_blind_family_is_reconstructed"]


def test_measure_the_antipodal_endpoint_on_its_own():
    r = measure_the_antipodal_endpoint_on_its_own()
    assert r["the_antipodal_value_is_minus_g0"]
    assert r["the_invariant_diverges_like_one_over_epsilon"]
    assert r["the_defect_stays_zero"]
    assert r["the_identity_survives_the_endpoint"]


# ── the corrections from the #258 review ────────────────────────────────────
def test_the_cross_term_comes_from_a_real_quadratic_functional():
    """Not a multiplication by zero: three evaluations of a functional that
    carries its own self-energy terms."""
    pair = _working()
    pts = random_points(6, seed=53)
    a, b = 1.7, -0.9
    for k in range(3):
        y_a, y_b = pts[2 * k], pts[2 * k + 1]
        cross = (energy_functional(pair, y_a, y_b, a, b)
                 - energy_functional(pair, y_a, y_b, a, 0.0)
                 - energy_functional(pair, y_a, y_b, 0.0, b))
        assert cross == pytest.approx(a * b * interaction_energy(pair, y_a,
                                                                 y_b),
                                      abs=1e-14)
        # the self-energies are actually there and are not small
        assert abs(energy_functional(pair, y_a, y_b, a, 0.0)) > 1e-3
    y_a, y_b = pts[0], pts[1]
    assert energy_functional(pair, y_a, y_b, 0.0, 0.0) == 0.0


def test_a_real_field_needs_a_real_boundary_matrix():
    """Conjugation-invariance of the domain is ``A = A*``, and a complex ``β``
    makes a real static source produce a complex field."""
    d = WORKING_SEPARATION
    pts = random_points(2, seed=59)
    v_a = source_vector(pts[0], d)
    v_b = source_vector(pts[1], d)
    for beta, real_ok in ((0.06, True), (complex(0.06, 0.20), False)):
        pair = MouthPair(d, 0.30, 0.35, beta)
        assert is_real_field_compatible(pair.boundary_matrix()) is real_ok
        field = complex(free_interaction_energy(pts[0], pts[1])
                        + v_b @ response_of_pair(pair) @ v_a)
        assert (abs(field.imag) < 1e-15) is real_ok


def test_the_blind_family_is_outside_the_real_field_sector():
    """Every invisible point needs Im β != 0, so there is no blind family for a
    real scalar."""
    d = WORKING_SEPARATION
    for a1, a2, rb in ((0.30, 0.35, -0.05), (0.50, 0.40, -0.02),
                       (0.25, 0.25, -0.10)):
        ib = invisible_partner(a1, a2, rb, d)
        assert ib is not None and ib > 1e-3
        pair = MouthPair(d, a1, a2, complex(rb, ib))
        assert not is_real_field_compatible(pair.boundary_matrix())
        assert abs(complex(pair.beta)) < min(a1, a2)      # comparable, smaller


def test_phase_sensitive_sources_need_only_one_spectral_parameter():
    """``A = Γ + R⁻¹`` once both quadratures are measured."""
    d = WORKING_SEPARATION
    pair = MouthPair(d, 0.30, 0.35, complex(-0.05, 0.24))
    quad = recover_complex_response(pair, *random_points(2, seed=61))
    assert quad["the_quadratures_give_the_kernel"]
    r = response_of_pair(pair)
    assert np.abs(gamma_at(0.0, d) + np.linalg.inv(r)
                  - pair.boundary_matrix()).max() < 1e-12


def test_the_reconstruction_uses_positive_spectral_parameters():
    """λ = ω², so a negative λ is an imaginary frequency; the default pair is
    positive and below the free ground state."""
    pair = MouthPair(WORKING_SEPARATION, 0.3, 0.4, complex(0.05, 0.12))
    out = recover_boundary(pair)
    assert all(0.0 < lam < 1.0 for lam in out["lambdas"])
    assert out["max_parameter_error"] < 1e-9
    assert out["residual"] < 1e-9


def test_the_reconstruction_reports_a_failure_to_converge():
    """A single bad start does land in a local minimum; the residual is what
    catches it, so it has to be reported and small."""
    rng = np.random.default_rng(20260821)
    for _ in range(6):
        a1, a2 = rng.uniform(0.15, 0.6, size=2)
        rb, ib = rng.uniform(-0.2, 0.2), rng.uniform(-0.3, 0.3)
        out = recover_boundary(MouthPair(WORKING_SEPARATION, float(a1),
                                         float(a2), complex(rb, ib)))
        assert out["residual"] < 1e-9
        assert out["max_parameter_error"] < 1e-9


def test_the_off_diagonal_block_is_not_a_throat_signature():
    """A disconnected pair fills the cross-mouth channel too."""
    d = WORKING_SEPARATION
    pts = random_points(2, seed=67)
    disc = MouthPair(d, 0.30, 0.35, 0.0)
    parts = mouth_channel_invariant(disc, pts[0], pts[1])
    assert abs(parts["cross_mouth"]) > 1e-6
    assert abs(defect_of_pair(disc)) < 1e-12       # …and W still says "no"


def test_measure_a_real_field_forces_beta_real():
    r = measure_a_real_field_forces_beta_real()
    assert r["a_real_beta_gives_a_real_field"]
    assert r["a_complex_beta_does_not"]
    assert r["so_for_PR254s_field_there_is_no_blind_family"]


def test_measure_phase_sensitive_sources_need_only_one_spectral_parameter():
    r = measure_phase_sensitive_sources_need_only_one_spectral_parameter()
    assert r["the_quadratures_give_the_kernel"]
    assert r["one_spectral_parameter_suffices"]
    assert r["worst_boundary_error"] < 1e-9
