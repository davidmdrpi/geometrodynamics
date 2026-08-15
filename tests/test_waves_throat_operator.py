"""
Tests for the flux-conserving two-mouth throat operator.

The load-bearing claims:

* the free ESU Green function has a **closed form** whose short-distance
  expansion is ``1/(4πχ) + g(ω)``, so a point interaction is definable at all;
* self-adjoint two-mouth boundary conditions are ``U(2)`` — the Cayley transform
  of a Hermitian ``A`` is unitary and inverts back, with **reflection** on the
  diagonal and **transmission** off it, ``|r|² + |t|² = 1``;
* ``A`` Hermitian **is** flux conservation, ``Im(q† A q) = 0`` identically;
* and therefore the coupled spectrum is **real for every coupling** — which
  retires PR #255's off-axis poles as an artefact of its directional,
  non-conserving relation, measured here against that exact control;
* the coupled spectrum **interlaces** the free one two per gap, and returns to
  ``ω = n+1`` as ``‖A‖ → ∞`` (which is the decoupling limit, not ``A → 0``).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.throat_operator import (
    DirectionalThroat,
    MouthPair,
    boundary_from_scattering,
    complex_root_search,
    coupled_spectrum,
    esu_free_spectrum,
    free_green,
    measure_self_adjointness_makes_the_spectrum_real,
    measure_the_boundary_operator_is_unitary_with_both_channels,
    measure_the_coupled_spectrum_interlaces_the_free_one,
    measure_the_directional_model_is_what_pr255_solved,
    measure_the_flux_balance_is_exactly_hermiticity,
    measure_the_green_function_has_a_finite_part,
    mode_charges,
    mouth_flux,
    regularized_green,
    spectrum_by_channel,
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
    """``sin(ω(π−χ))`` vanishes exactly where ``sin χ`` does, so the antipodal
    focus is *finite* — ``G → ω/(4π sin πω)`` — unlike the coincidence point."""
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


# ── the operator is self-adjoint ────────────────────────────────────────────
def test_the_boundary_matrix_is_hermitian_by_construction():
    p = MouthPair(1.3, 0.2, -0.13, 0.15 + 0.07j)
    a = p.boundary_matrix()
    assert p.is_self_adjoint()
    assert np.abs(a - a.conjugate().T).max() == 0.0


def test_gamma_is_real_symmetric_on_the_real_axis():
    g = MouthPair(1.3).gamma(2.3)
    assert np.abs(g.imag).max() == 0.0
    assert np.abs(g - g.T).max() == 0.0


def test_the_secular_function_is_real_for_hermitian_boundary_data():
    p = MouthPair(1.3, 0.2, -0.13, 0.15 + 0.07j)
    for w in (1.3, 2.77, 4.11):
        det = complex(np.linalg.det(p.krein_matrix(w)))
        assert abs(det.imag) < 1e-14 * max(abs(det), 1.0)


def test_exchange_symmetry_is_detected():
    assert MouthPair(1.3, 0.05, 0.05, 0.03).is_exchange_symmetric()
    assert not MouthPair(1.3, 0.05, 0.06, 0.03).is_exchange_symmetric()
    assert not MouthPair(1.3, 0.05, 0.05, 0.03j).is_exchange_symmetric()


# ── the unitary boundary operator ───────────────────────────────────────────
@pytest.mark.parametrize("scale", [0.05, 0.1, 0.2])
def test_the_cayley_transform_is_unitary(scale):
    p = MouthPair(1.3, 0.2, -0.13, 0.15 + 0.07j)
    s = p.scattering_matrix(scale)
    assert np.abs(s.conjugate().T @ s - np.eye(2)).max() < 1e-13


@pytest.mark.parametrize("scale", [0.05, 0.2])
def test_the_cayley_transform_inverts(scale):
    p = MouthPair(1.3, -0.4, 0.07, -0.09 + 0.31j)
    back = boundary_from_scattering(p.scattering_matrix(scale), scale)
    assert np.abs(back - p.boundary_matrix()).max() < 1e-12
    assert np.abs(back - back.conjugate().T).max() < 1e-12


def test_reflection_and_transmission_are_both_present_and_conserve():
    """The channel PR #255 had none of, and the sum rule it could not satisfy."""
    ch = MouthPair(1.3, 0.2, 0.2, 0.15).channels(0.1)
    assert abs(ch["reflection_1"]) > 1e-3
    assert abs(ch["transmission_12"]) > 1e-3
    assert ch["sum_of_squares_mouth_1"] == pytest.approx(1.0, abs=1e-13)
    assert ch["sum_of_squares_mouth_2"] == pytest.approx(1.0, abs=1e-13)


def test_a_real_beta_is_reciprocal_and_a_complex_one_is_not():
    assert MouthPair(1.3, 0.2, 0.2, 0.15).channels()["reciprocal"]
    assert not MouthPair(1.3, 0.2, 0.2, 0.15 + 0.05j).channels()["reciprocal"]


def test_a_decoupled_mouth_pair_transmits_nothing():
    ch = MouthPair(1.3, 0.2, 0.2, 0.0).channels(0.1)
    assert abs(ch["transmission_12"]) < 1e-15
    assert abs(abs(ch["reflection_1"]) - 1.0) < 1e-13


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


def test_the_directional_control_does_not_conserve_flux():
    ctl = DirectionalThroat(1.3, 1.0, +1, 0.3)
    rng = np.random.default_rng(3)
    nets = []
    for _ in range(20):
        q = rng.normal(0, 1, 2) + 1j * rng.normal(0, 1, 2)
        f = mouth_flux(q, ctl.boundary_matrix(2.3))
        nets.append(abs(f["net"]) / f["scale"])
    assert float(np.median(nets)) > 1e-3
    assert ctl.anti_hermitian_part(2.3) > 0.0


# ── the spectrum ────────────────────────────────────────────────────────────
def test_the_coupled_spectrum_has_exactly_two_roots_per_gap():
    spec = coupled_spectrum(MouthPair(1.3, 0.05, 0.05, 0.03), 6)
    per = {m: 0 for m in range(1, 7)}
    for r in spec:
        per[r["gap"]] += 1
        assert r["gap"] < r["omega"] < r["gap"] + 1
        assert abs(r["secular"]) < 1e-6
    assert set(per.values()) == {2}


def test_the_channel_split_reproduces_the_determinant_roots():
    p = MouthPair(1.3, 0.05, 0.05, 0.03)
    by_ch = spectrum_by_channel(p, 5)["rows"]
    det = sorted(r["omega"] for r in coupled_spectrum(p, 5))
    chan = sorted(w for r in by_ch
                  for w in (r["symmetric"], r["antisymmetric"])
                  if w is not None)
    assert len(det) == len(chan)
    for a, b in zip(det, chan):
        assert a == pytest.approx(b, abs=1e-8)


def test_the_channel_functions_are_monotone_across_every_gap():
    """Why there is exactly one root per channel per gap."""
    p = MouthPair(1.3, 0.05, 0.05, 0.03)
    for m in (1, 2, 3):
        ws = np.linspace(m + 1e-6, m + 1 - 1e-6, 1500)
        for k in (0, 1):
            v = np.array([p.channel_functions(float(w))[k] for w in ws])
            assert (np.diff(v) > 0).all()


def test_channel_splitting_needs_an_exchange_symmetric_pair():
    with pytest.raises(ValueError):
        spectrum_by_channel(MouthPair(1.3, 0.05, 0.06, 0.03))


def test_the_free_spectrum_returns_as_the_boundary_norm_grows():
    """Decoupling is ``‖A‖ → ∞`` — the diagonal is an *inverse* length."""
    worst = []
    for t in (1e2, 1e3, 1e4):
        p = MouthPair(1.3, 0.05 * t, 0.05 * t, 0.03 * t)
        got = [w for r in spectrum_by_channel(p, 4)["rows"]
               for w in (r["symmetric"], r["antisymmetric"]) if w is not None]
        worst.append(max(abs(w - round(w)) for w in got))
    assert worst[0] > worst[1] > worst[2]
    assert worst[2] < 1e-3
    assert worst[1] / worst[2] == pytest.approx(10.0, rel=0.05)


def test_mode_charges_annihilate_the_krein_matrix():
    p = MouthPair(1.3, 0.05, 0.05, 0.03)
    for r in coupled_spectrum(p, 3):
        q = mode_charges(r["omega"], p)
        assert np.linalg.norm(p.krein_matrix(r["omega"]) @ q) < 1e-6


# ── the headline: no off-axis roots ─────────────────────────────────────────
def test_a_self_adjoint_throat_has_no_root_off_the_real_axis():
    p = MouthPair(1.3, 0.2, -0.13, 0.15 + 0.07j)

    def sec(z):
        return complex(np.linalg.det(_krein(z, p)))

    r = complex_root_search(sec, (1.1, 6.9))
    assert r["n_roots"] > 0
    assert r["n_off_axis"] == 0
    assert r["worst_abs_imaginary"] < 1e-12


@pytest.mark.parametrize("kappa", [0.3, 1.0])
def test_the_directional_control_puts_roots_off_the_axis(kappa):
    """Including at unit transmission — it is the directionality, not the loss."""
    r = complex_root_search(
        DirectionalThroat(1.3, 1.0, +1, kappa).secular, (1.1, 6.9))
    assert r["n_off_axis"] == r["n_roots"] > 0
    assert r["worst_abs_imaginary"] > 1e-3
    assert r["n_growing"] > 0


def _krein(z, pair):
    from geometrodynamics.waves.throat_operator import _krein_complex
    return _krein_complex(z, pair)


# ── the measurements ────────────────────────────────────────────────────────
def test_measure_the_green_function_has_a_finite_part():
    r = measure_the_green_function_has_a_finite_part()
    assert r["the_closed_form_is_the_branch_series"]
    assert r["the_remainder_is_first_order_in_chi"]


def test_measure_the_boundary_operator_is_unitary_with_both_channels():
    r = measure_the_boundary_operator_is_unitary_with_both_channels()
    assert r["the_cayley_transform_is_unitary"]
    assert r["every_mouth_conserves"]
    assert r["both_channels_are_present"]
    assert r["pr255_is_outside_U2_unless_kappa_is_one"]
    assert r["worst_unitarity_defect"] < 1e-13


def test_measure_self_adjointness_makes_the_spectrum_real():
    """The result this round exists for."""
    r = measure_self_adjointness_makes_the_spectrum_real()
    assert r["the_secular_function_is_real_on_the_axis"]
    assert r["nothing_off_the_axis"]
    assert r["worst_abs_imaginary_over_all_seeds"] < 1e-12
    assert r["the_control_fails_both_tests"]
    assert r["and_the_control_is_unstable_even_at_unit_transmission"]
    assert r["control_directional"]["n_off_axis"] > 0
    assert r["control_directional"]["n_growing"] > 0


def test_measure_the_coupled_spectrum_interlaces_the_free_one():
    r = measure_the_coupled_spectrum_interlaces_the_free_one()
    assert r["exactly_two_per_gap"]
    assert r["every_root_strictly_between_free_ones"]
    assert r["both_channel_functions_are_monotone"]
    assert r["the_shift_goes_like_one_over_the_boundary_norm"]
    assert r["free_spectrum_recovered"]


def test_measure_the_flux_balance_is_exactly_hermiticity():
    r = measure_the_flux_balance_is_exactly_hermiticity()
    assert r["flux_is_conserved_identically"]
    assert r["what_one_mouth_absorbs_the_other_emits"]
    assert r["the_control_does_not_conserve"]
    assert r["worst_relative_net_flux"] < 1e-14


def test_measure_the_directional_model_is_what_pr255_solved():
    r = measure_the_directional_model_is_what_pr255_solved()
    assert r["no_reflection_channel"]
    assert r["not_hermitian_at_any_frequency"]
    assert r["anti_hermitian_part_is_the_whole_matrix"]
    assert r["the_boundary_matrix_depends_on_frequency"]
    assert r["only_through_its_phase"]
    assert "exact for the model it posed" in r["what_is_not_wrong_with_it"]


def test_the_module_says_what_is_still_put_in():
    from geometrodynamics.waves import throat_operator
    doc = throat_operator.__doc__ or ""
    assert "point-supported" in doc
    assert "no backreaction" in doc.lower()
    assert "is a choice" in doc
