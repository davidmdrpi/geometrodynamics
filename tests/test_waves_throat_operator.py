"""
Tests for the flux-conserving two-mouth throat operator.

The load-bearing claims, after the review:

* the free ESU Green function has a **closed form** whose short-distance
  expansion is ``1/(4πχ) + g(ω)``, so a point interaction is definable at all;
* a linear two-mouth boundary condition is a pair ``B φ^reg = C q``; it is
  self-adjoint iff ``rank[B|C] = 2`` and ``BC†`` is Hermitian, and that is
  **exactly** the statement that no flux is created at the mouths;
* self-adjointness makes ``λ = ω²`` real and **nothing more** — ``λ`` can be
  negative, and then ``ω`` is imaginary and a mode grows.  Positivity is a
  separate condition with a closed form for the exchange-symmetric family;
* ``det(C − BΓ) = 0`` is the **rank-two mouth-active sector**: ``(n+1)² − 2``
  modes per level never move, and the sector has a piece below the free ground
  state as well as the interlacing pieces;
* and PR #255's relation embeds exactly as a **non-self-adjoint** boundary
  condition with singular ``B`` — checked against that round's own ``1 − L``.
"""

from __future__ import annotations

import cmath
import math

import numpy as np
import pytest

from geometrodynamics.waves.throat_operator import (
    BoundaryCondition,
    DirectionalThroat,
    MouthPair,
    boundary_from_scattering,
    boundary_mixing,
    channel_endpoints,
    esu_free_spectrum,
    free_green,
    gamma_at,
    is_stable,
    measure_hermiticity_is_flux_conservation,
    measure_self_adjointness_makes_lambda_real_not_omega,
    measure_the_green_function_has_a_finite_part,
    measure_the_mouth_active_sector_is_rank_two,
    measure_the_pr255_boundary_condition_is_not_self_adjoint,
    measure_the_spectrum_is_conjugation_symmetric_in_beta,
    measure_the_stability_region_in_the_boundary_family,
    mode_charges,
    mouth_active_spectrum,
    mouth_flux,
    negative_lambda_modes,
    regularized_green,
    spectrum_by_channel,
    stability_thresholds,
    untouched_multiplicity,
)


# ── the free Green function ─────────────────────────────────────────────────
@pytest.mark.parametrize("chi", [0.4, 1.3, 2.6])
@pytest.mark.parametrize("omega", [0.37, 1.63, 2.4, 5.21])
def test_the_closed_form_is_the_branch_series_limit(chi, omega):
    """An independent construction: PR #255 summed images, this sums nothing."""
    from geometrodynamics.waves.branch_coupling import mouth_transfer
    got = mouth_transfer(omega, chi, 1e-6)
    assert got.real == pytest.approx(free_green(chi, omega), abs=1e-10)
    assert abs(got.imag) < 1e-5


def test_the_green_function_is_real_on_the_real_axis():
    """A closed universe has no radiation condition to make it complex."""
    for w in (0.4, 1.7, 3.3):
        assert isinstance(free_green(1.3, w), float)


def test_the_green_function_is_finite_at_the_antipode():
    """``sin(ω(π−χ))`` vanishes where ``sin χ`` does, so the antipodal focus is
    finite — ``G → ω/(4π sin πω)`` — unlike the coincidence point."""
    for w in (0.7, 2.3, 4.4):
        want = w / (4.0 * math.pi * math.sin(math.pi * w))
        errs = [abs(free_green(math.pi - e, w) - want) / abs(want)
                for e in (1e-3, 1e-5)]
        assert errs[0] < 1e-4 and errs[1] < 1e-8
        assert errs[0] > errs[1]


@pytest.mark.parametrize("bad", [0.0, math.pi])
def test_the_green_function_is_singular_at_the_poles(bad):
    with pytest.raises(ValueError):
        free_green(bad, 1.3)


@pytest.mark.parametrize("m", [1, 2, 3])
def test_the_green_function_has_a_pole_at_every_free_eigenfrequency(m):
    with pytest.raises(ValueError):
        free_green(1.3, float(m))
    with pytest.raises(ValueError):
        regularized_green(float(m))


def test_the_finite_part_is_the_regularized_green_function():
    """``G(χ) − 1/(4πχ) → g(ω)``, first order in ``χ``."""
    for w in (0.37, 1.63, 2.4):
        errs = []
        for chi in (1e-2, 1e-3, 1e-4):
            fp = free_green(chi, w) - 1.0 / (4.0 * math.pi * chi)
            errs.append(abs(fp - regularized_green(w)))
        assert errs[0] > errs[1] > errs[2]
        assert errs[0] / errs[1] == pytest.approx(10.0, rel=0.05)


def test_the_free_spectrum_is_the_integers():
    assert esu_free_spectrum(4) == [1.0, 2.0, 3.0, 4.0]


# ── Γ on both sides of λ = 0 ────────────────────────────────────────────────
def test_gamma_is_real_symmetric_for_every_sign_of_lambda():
    """Which is why the eigenvalues are real — on both sides of threshold."""
    for lam in (-9.0, -0.4, 0.0, 0.5, 5.3, 12.0):
        g = gamma_at(lam, 1.3)
        assert np.abs(g.imag).max() == 0.0
        assert np.abs(g - g.T).max() == 0.0


def test_gamma_is_continuous_through_threshold():
    below = gamma_at(-1e-9, 1.3)
    at = gamma_at(0.0, 1.3)
    above = gamma_at(1e-9, 1.3)
    assert np.abs(below - at).max() < 1e-7
    assert np.abs(above - at).max() < 1e-7


def test_the_threshold_values_are_the_closed_form():
    th = stability_thresholds(1.3)
    assert th["g_at_zero"] == pytest.approx(-1.0 / (4.0 * math.pi ** 2))
    assert th["G_d_at_zero"] == pytest.approx(
        (math.pi - 1.3) / (4.0 * math.pi ** 2 * math.sin(1.3)))
    at = gamma_at(0.0, 1.3)
    assert float(at[0, 0].real) == pytest.approx(th["g_at_zero"])
    assert float(at[0, 1].real) == pytest.approx(th["G_d_at_zero"])


def test_the_channel_functions_fall_monotonically_along_the_imaginary_axis():
    """Which is why the stability condition is a threshold and not a scan."""
    p = MouthPair(1.3, 0.0, 0.0, 0.0)
    lams = -np.geomspace(1e-6, 900.0, 3000)
    for k in (0, 1):
        v = np.array([p.channel_functions(float(x))[k] for x in lams])
        assert (np.diff(v) < 0).all()


# ── boundary conditions ─────────────────────────────────────────────────────
def test_a_hermitian_A_is_a_self_adjoint_boundary_condition():
    bc = MouthPair(1.3, 0.2, -0.13, 0.15 + 0.07j).boundary_condition()
    assert bc.is_maximal()
    assert bc.is_self_adjoint()
    assert bc.self_adjointness_defect() == 0.0


def test_a_non_hermitian_A_is_maximal_but_not_self_adjoint():
    bc = BoundaryCondition(np.eye(2, dtype=complex),
                           np.array([[0.2, 0.3], [0.1, -0.1]], dtype=complex))
    assert bc.is_maximal()
    assert not bc.is_self_adjoint()
    assert bc.self_adjointness_defect() > 0.0


def test_exchange_symmetry_is_detected():
    assert MouthPair(1.3, 0.05, 0.05, 0.03).is_exchange_symmetric()
    assert not MouthPair(1.3, 0.05, 0.06, 0.03).is_exchange_symmetric()
    assert not MouthPair(1.3, 0.05, 0.05, 0.03j).is_exchange_symmetric()


def test_the_secular_function_is_real_in_lambda_on_both_sides():
    p = MouthPair(1.3, 0.2, -0.13, 0.15 + 0.07j)
    for lam in (-12.0, -0.3, 0.5, 5.3, 11.0):
        det = complex(np.linalg.det(p.krein_matrix(lam)))
        assert abs(det.imag) < 1e-13 * max(abs(det), 1.0)


# ── the Cayley form is a parametrization, not a scattering matrix ───────────
@pytest.mark.parametrize("scale", [0.05, 0.1, 0.2])
def test_the_cayley_transform_is_unitary(scale):
    m = boundary_mixing(MouthPair(1.3, 0.2, -0.13, 0.15 + 0.07j)
                        .boundary_matrix(), scale)
    assert m["unitarity_defect"] < 1e-13
    for c in m["column_norms"]:
        assert c == pytest.approx(1.0, abs=1e-13)


@pytest.mark.parametrize("scale", [0.05, 0.2])
def test_the_cayley_transform_inverts(scale):
    p = MouthPair(1.3, -0.4, 0.07, -0.09 + 0.31j)
    back = boundary_from_scattering(
        boundary_mixing(p.boundary_matrix(), scale)["S"], scale)
    assert np.abs(back - p.boundary_matrix()).max() < 1e-12
    assert np.abs(back - back.conjugate().T).max() < 1e-12


def test_the_cayley_entries_depend_on_the_reference_scale():
    """So they are boundary-mixing coefficients, not reflection and
    transmission — the correction the review forced."""
    a = MouthPair(1.3, 0.2, -0.13, 0.15 + 0.07j).boundary_matrix()
    diag = [boundary_mixing(a, c)["diagonal_mixing"][0]
            for c in (0.05, 0.1, 0.2)]
    assert max(diag) - min(diag) > 0.1
    assert "not physical r and t" in boundary_mixing(a)["caveat"]


# ── flux ────────────────────────────────────────────────────────────────────
def test_hermitian_boundary_data_conserves_flux_identically():
    rng = np.random.default_rng(7)
    p = MouthPair(1.3, 0.31, -0.17, -0.2 + 0.4j)
    for _ in range(50):
        q = rng.normal(0, 1, 2) + 1j * rng.normal(0, 1, 2)
        f = mouth_flux(q, p.boundary_matrix())
        assert abs(f["net"]) / f["scale"] < 1e-14


def test_a_pure_transmission_throat_moves_flux_between_the_mouths():
    rng = np.random.default_rng(11)
    p = MouthPair(1.3, 0.0, 0.0, 0.25 + 0.1j)
    q = rng.normal(0, 1, 2) + 1j * rng.normal(0, 1, 2)
    f = mouth_flux(q, p.boundary_matrix())
    assert f["per_mouth"][0] == pytest.approx(-f["per_mouth"][1], abs=1e-15)


def test_a_non_hermitian_boundary_matrix_does_not_conserve_flux():
    a = np.array([[0.0, 0.0], [3.33, 0.0]], dtype=complex)
    rng = np.random.default_rng(3)
    nets = [abs(mouth_flux(rng.normal(0, 1, 2) + 1j * rng.normal(0, 1, 2),
                           a)["net"]) for _ in range(20)]
    assert float(np.median(nets)) > 1e-3


# ── positivity is a separate condition ──────────────────────────────────────
def test_the_advertised_examples_have_growing_modes():
    """The claim this round originally got wrong, pinned as a test."""
    got = {}
    for (a1, a2, b, want_sigma) in (
            (0.2, -0.13, 0.15 + 0.07j, 2.470532),
            (-0.4, 0.07, -0.09 + 0.31j, 7.090982)):
        modes = negative_lambda_modes(MouthPair(1.3, a1, a2, b))
        assert len(modes) == 1
        assert modes[0]["sigma"] == pytest.approx(want_sigma, abs=1e-5)
        assert modes[0]["lmbda"] < 0.0
        got[a1] = modes[0]["growth_rate"]
    assert got[0.2] > 0.0 and got[-0.4] > 0.0


def test_the_default_exchange_symmetric_pair_is_stable():
    st = is_stable(MouthPair(1.3, 0.05, 0.05, 0.03))
    assert st["stable"] and st["n_growing_modes"] == 0
    assert st["closed_form_agrees_with_the_scan"]


@pytest.mark.parametrize("alpha,beta,stable", [
    (0.05, 0.03, True), (0.2, 0.0, True),
    (0.0, 0.0, False), (-0.05, 0.0, False), (0.05, 0.30, False)])
def test_the_closed_form_stability_test_matches_the_scan(alpha, beta, stable):
    st = is_stable(MouthPair(1.3, alpha, alpha, beta))
    assert st["stable"] is stable
    assert st["closed_form_stable"] is stable
    assert st["closed_form_agrees_with_the_scan"]


def test_the_stability_wedge_is_where_the_thresholds_say():
    th = stability_thresholds(1.3)
    assert th["symmetric_threshold"] == pytest.approx(0.02308202, abs=1e-7)
    assert th["antisymmetric_threshold"] == pytest.approx(-0.07374262,
                                                          abs=1e-7)
    # just inside and just outside the symmetric edge
    eps = 1e-4
    a = th["symmetric_threshold"]
    assert is_stable(MouthPair(1.3, a / 2 + eps, a / 2 + eps, a / 2))["stable"]
    assert not is_stable(MouthPair(1.3, a / 2 - eps, a / 2 - eps,
                                   a / 2))["stable"]


# ── the mouth-active sector, scoped ─────────────────────────────────────────
def test_only_two_modes_per_level_can_move():
    assert untouched_multiplicity(0) == 0
    assert untouched_multiplicity(1) == 2
    assert untouched_multiplicity(4) == 23
    for n in range(1, 6):
        assert untouched_multiplicity(n) == (n + 1) ** 2 - 2


def test_the_sector_has_a_piece_below_the_free_ground_state():
    """A piece an ω-scan starting above 1 cannot see."""
    spec = mouth_active_spectrum(MouthPair(1.3, 0.05, 0.05, 0.03), 4)
    below = [r for r in spec if r["sector"] == "below the free ground state"]
    assert len(below) == 1
    assert 0.0 < below[0]["lmbda"] < 1.0


def test_the_interlacing_gaps_carry_two_roots_each():
    spec = mouth_active_spectrum(MouthPair(1.3, 0.05, 0.05, 0.03), 5)
    per = {}
    for r in spec:
        if r["sector"] == "interlacing":
            per[r["gap"]] = per.get(r["gap"], 0) + 1
            assert r["gap"] ** 2 < r["lmbda"] < (r["gap"] + 1) ** 2
    assert set(per.values()) == {2}


def test_the_antisymmetric_channel_does_not_span_the_first_gap():
    """The ``m = 1`` pole cancels: the constant mode is equal at both mouths."""
    p = MouthPair(1.3, 0.05, 0.05, 0.03)
    first = channel_endpoints(p, 1)
    assert first["symmetric_spans_the_whole_line"]
    assert not first["antisymmetric_spans_the_whole_line"]
    assert abs(first["antisymmetric_at_lower"]) < 1.0
    assert channel_endpoints(p, 2)["antisymmetric_spans_the_whole_line"]


def test_a_root_in_the_first_gap_is_conditional_on_alpha_minus_beta():
    got = []
    for b in (0.03, 0.10, 0.14):
        r = spectrum_by_channel(MouthPair(1.3, 0.05, 0.05, b), 1)["rows"][0]
        got.append(r["antisymmetric"] is not None)
    assert len(set(got)) > 1


def test_channel_splitting_needs_an_exchange_symmetric_pair():
    with pytest.raises(ValueError):
        spectrum_by_channel(MouthPair(1.3, 0.05, 0.06, 0.03))


def test_mode_charges_annihilate_the_krein_matrix():
    p = MouthPair(1.3, 0.05, 0.05, 0.03)
    for r in mouth_active_spectrum(p, 3):
        if r["sector"] == "growing":
            continue
        q = mode_charges(r["lmbda"], p)
        assert np.linalg.norm(p.krein_matrix(r["lmbda"]) @ q) < 1e-5


# ── the PR #255 embedding ───────────────────────────────────────────────────
@pytest.mark.parametrize("kappa", [0.3, 1.0])
@pytest.mark.parametrize("omega", [1.3 + 0.2j, 2.77 - 0.4j, 4.11 + 0.05j])
def test_the_embedding_reproduces_pr255s_own_pole_condition(kappa, omega):
    ctl = DirectionalThroat(1.3, 1.0, +1, kappa)
    assert abs(ctl.secular(omega) - ctl.pr255_pole_condition(omega)) < 1e-13


def test_the_embedding_is_maximal_but_not_self_adjoint():
    bc = DirectionalThroat(1.3, 1.0, +1, 0.3).boundary_condition(2.3)
    assert bc.is_maximal()
    assert not bc.is_self_adjoint()
    assert bc.self_adjointness_defect() > 1e-3


def test_the_old_control_was_a_different_function():
    """Recorded rather than dropped: ``A = [[0,0],[1/gain,0]]`` is not it."""
    ctl = DirectionalThroat(1.3, 1.0, +1, 0.3)
    w = 2.77 - 0.4j
    sp = cmath.sin(math.pi * w)
    g = -w * cmath.cos(math.pi * w) / (4.0 * math.pi * sp)
    gd = cmath.sin(w * (math.pi - 1.3)) / (4.0 * math.pi * math.sin(1.3) * sp)
    old = complex(np.linalg.det(
        np.array([[0.0, 0.0], [1.0 / ctl.gain(w), 0.0]], dtype=complex)
        - np.array([[g, gd], [gd, g]], dtype=complex)))
    assert abs(old - ctl.pr255_pole_condition(w)) > 1e-2


# ── the phase of β ──────────────────────────────────────────────────────────
def test_the_spectrum_depends_on_the_phase_of_beta():
    """Γ is not diagonal, so the relative phase is physical."""
    dets = []
    for b in (0.1655 + 0j, 0.15 + 0.07j, 0.07 + 0.15j):
        dets.append(MouthPair(1.3, 0.2, -0.13, b).secular(5.3))
    assert max(dets) - min(dets) > 1e-6


def test_the_spectrum_is_invariant_under_conjugating_beta():
    """Time reversal — and the reason 'non-reciprocal' was the wrong word."""
    for b in (0.15 + 0.07j, -0.09 + 0.31j):
        a = MouthPair(1.3, 0.2, -0.13, b).secular(5.3)
        c = MouthPair(1.3, 0.2, -0.13, b.conjugate()).secular(5.3)
        assert a == pytest.approx(c, abs=1e-15)


# ── the measurements ────────────────────────────────────────────────────────
def test_measure_the_green_function_has_a_finite_part():
    r = measure_the_green_function_has_a_finite_part()
    assert r["the_closed_form_is_the_branch_series"]
    assert r["the_remainder_is_first_order_in_chi"]
    assert r["the_antipodal_focus_is_finite"]


def test_measure_hermiticity_is_flux_conservation():
    r = measure_hermiticity_is_flux_conservation()
    assert r["flux_is_conserved_identically"]
    assert r["what_one_mouth_absorbs_the_other_emits"]
    assert r["the_control_does_not_conserve"]
    assert r["the_cayley_transform_is_unitary"]
    assert r["the_cayley_entries_are_not_physical_amplitudes"]


def test_measure_self_adjointness_makes_lambda_real_not_omega():
    """The correction: Hermiticity buys real λ, not positive λ."""
    r = measure_self_adjointness_makes_lambda_real_not_omega()
    assert r["the_secular_function_is_real_in_lambda"]
    assert r["hermiticity_gives_real_lambda"]
    assert r["hermiticity_does_not_give_positivity"]
    assert r["n_unstable_examples"] == 2
    assert "withdrawn" in " ".join(r.keys())
    assert "imaginary axis" in r["why_it_was_missed"]


def test_measure_the_stability_region_in_the_boundary_family():
    r = measure_the_stability_region_in_the_boundary_family()
    assert r["the_channel_functions_are_monotone_on_the_imaginary_axis"]
    assert r["the_closed_form_agrees_with_every_probe"]
    assert r["the_closed_form_matches_everywhere"]
    assert r["grid_mismatches"] == 0
    assert r["both_signs_are_represented"]


def test_measure_the_mouth_active_sector_is_rank_two():
    r = measure_the_mouth_active_sector_is_rank_two()
    assert r["at_most_two_modes_move_per_level"]
    assert r["untouched_modes_at_level_4"] == 23
    assert r["there_is_a_sector_below_the_ground_state"]
    assert r["two_per_interlacing_gap"]
    assert r["the_first_gap_antisymmetric_endpoints_are_finite"]
    assert r["the_symmetric_channel_does_span_it"]
    assert r["existence_in_the_first_gap_is_conditional"]


def test_measure_the_pr255_boundary_condition_is_not_self_adjoint():
    r = measure_the_pr255_boundary_condition_is_not_self_adjoint()
    assert r["the_embedding_is_exact"]
    assert r["worst_embedding_error"] < 1e-12
    assert r["the_old_control_was_a_different_model"]
    assert r["every_embedding_is_maximal"]
    assert r["none_is_self_adjoint"]
    assert "classification, not a diagnosis" in r["what_is_not_concluded"]


def test_measure_the_spectrum_is_conjugation_symmetric_in_beta():
    r = measure_the_spectrum_is_conjugation_symmetric_in_beta()
    assert r["the_phase_of_beta_is_physical"]
    assert r["the_spectrum_is_conjugation_symmetric"]
    assert "non-reciprocal" in r["the_withdrawn_claim"]


def test_the_module_says_what_is_still_put_in():
    from geometrodynamics.waves import throat_operator
    doc = throat_operator.__doc__ or ""
    assert "point-supported" in doc
    assert "no backreaction" in doc.lower()
    assert "It does not buy" in doc
    assert "mouth-active" in doc
