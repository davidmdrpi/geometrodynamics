"""
Tests for the ℓ-resolved coupling between two concentric shells.

The load-bearing claim is that the mutual stiffness carries the Laplacian
eigenvalue ``ℓ(ℓ+1)``, so the previous round's decoupling is the ``ℓ = 0`` case
of a multipole fact rather than a separate theorem about spheres.  That is
checked against brute-force integration which never expands in multipoles.

Also held down, because both were got wrong first:

* a pure ``P₁`` deformation is **not** a translation past linear order, so its
  area cost is ``4/3`` and not zero — a naive translation-invariance test would
  report a stiffness that is not there;
* translation invariance nevertheless holds exactly, along the family that
  carries the induced ``ℓ = 0`` and ``ℓ = 2`` pieces.

And the constitutive gap is asserted rather than hidden: a perfect fluid shell
has no shear modulus, so ``μ_s`` is an input.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.shells.multipole import (
    ShellPair,
    area_second_variation,
    measure_a_fluid_shell_has_no_shear_modulus,
    measure_the_ell_zero_decoupling_is_a_zero_eigenvalue,
    measure_the_translation_mode_does_not_couple,
    rigid_pair_mutual_energy,
    measure_the_area_cost_of_a_deformation,
    measure_the_coupling_is_screened_by_separation,
    measure_the_mutual_coupling_is_the_laplacian_eigenvalue,
    measure_the_pure_dipole_is_not_a_translation,
    measure_translation_invariance_is_exact,
    mutual_stiffness,
    transfer_exponent,
    translation_family,
)


# ── the central result ──────────────────────────────────────────────────────
def test_the_mutual_stiffness_vanishes_only_at_ell_zero():
    for b, a in ((2.0, 5.0), (1.0, 100.0), (4.9, 5.0)):
        assert mutual_stiffness(0, b, a) == 0.0
        for ell in (1, 2, 3, 8):
            assert mutual_stiffness(ell, b, a) > 0.0


def test_the_prefactor_is_the_laplacian_eigenvalue():
    b, a = 2.0, 5.0
    for ell in (1, 2, 3, 5):
        expected = (ell * (ell + 1) * (b / a) ** ell
                    / (a * (2 * ell + 1) ** 2))
        assert mutual_stiffness(ell, b, a) == pytest.approx(expected,
                                                            rel=1e-14)


def test_the_closed_form_matches_brute_force_integration():
    """The control never expands in multipoles, so it can disagree."""
    r = measure_the_mutual_coupling_is_the_laplacian_eigenvalue()
    assert r["the_closed_form_is_confirmed"] is True
    assert r["worst_relative_error"] < 1e-4
    zero = [row for row in r["rows"] if row["ell"] == 0][0]
    assert zero["closed_form"] == 0.0
    assert abs(zero["brute_force"]) < 1e-9


def test_the_ell_zero_decoupling_is_a_zero_eigenvalue():
    r = measure_the_ell_zero_decoupling_is_a_zero_eigenvalue()
    assert r["ell_zero_is_exactly_zero_at_every_separation"] is True
    assert r["every_other_ell_couples"] is True
    for row in r["rows"]:
        assert row["ell_0"] == 0.0
        assert row["ell_2"] > 0.0


def test_the_ell_zero_result_is_not_claimed_to_be_birkhoff():
    """Birkhoff is a GR theorem; this is its static Newtonian analogue."""
    r = measure_the_ell_zero_decoupling_is_a_zero_eigenvalue()
    note = r["this_is_not_birkhoff"]
    assert "GR theorem" in note
    assert "Newtonian analogue" in note


# ── the ℓ = 1 control, done on the energy ───────────────────────────────────
def test_the_shell_theorem_holds_for_the_mutual_energy():
    """A rigid displacement leaves the mutual energy at ``−G m_b m_a / a``."""
    b, a = 2.0, 5.0
    ref = rigid_pair_mutual_energy(b, a)
    assert ref == pytest.approx(-1.0 / a, abs=1e-12)
    for d in (0.1, 0.8, 2.5):
        assert rigid_pair_mutual_energy(b, a, d_b=d) == pytest.approx(
            ref, rel=1e-12)


def test_the_translation_mode_does_not_couple_but_the_shape_mode_does():
    """The correction: ``ℓ = 1`` couples as a *shape*, not as a translation."""
    r = measure_the_translation_mode_does_not_couple()
    assert r["shell_theorem_holds"] is True
    assert abs(r["rigid_translation_cross_derivative"]) < 1e-9
    assert r["pure_P1_shape_coupling"] == pytest.approx(1.7778e-02, rel=1e-3)
    assert r["the_translation_mode_does_not_couple"] is True
    assert r["the_shape_mode_does"] is True
    assert r["so_coupling_starts_at_ell_two"] is True


def test_the_area_test_is_recorded_as_insufficient():
    r = measure_the_translation_mode_does_not_couple()
    why = r["why_the_area_test_was_not_enough"]
    assert "not about the mutual gravitational energy" in why
    assert "it does not" in why


def test_the_mutual_stiffness_scales_with_the_masses():
    base = mutual_stiffness(2, 2.0, 5.0, m_b=1.0, m_a=1.0)
    assert mutual_stiffness(2, 2.0, 5.0, m_b=3.0, m_a=1.0) == \
        pytest.approx(3.0 * base, rel=1e-14)
    assert mutual_stiffness(2, 2.0, 5.0, m_b=2.0, m_a=5.0) == \
        pytest.approx(10.0 * base, rel=1e-14)


# ── the trap ────────────────────────────────────────────────────────────────
def test_a_pure_dipole_deformation_costs_area():
    """``4/3``, not zero — which is why the naive zero-mode test is wrong."""
    assert area_second_variation(1) == pytest.approx(4.0 / 3.0, rel=1e-14)
    r = measure_the_pure_dipole_is_not_a_translation()
    assert r["the_dipole_area_cost_is_not_zero"] is True
    assert r["dipole_second_variation"] == pytest.approx(4.0 / 3.0, rel=1e-6)
    assert r["so_a_pure_P1_test_would_be_wrong"] is True


def test_the_area_second_variation_is_the_exact_rational():
    for ell, val in ((0, 2.0), (1, 4 / 3), (2, 8 / 5), (3, 2.0),
                     (4, 22 / 9), (5, 32 / 11)):
        assert area_second_variation(ell) == pytest.approx(val, rel=1e-14)
    r = measure_the_area_cost_of_a_deformation()
    assert r["every_value_is_the_exact_rational"] is True


def test_the_area_second_variation_matches_direct_integration():
    r = measure_the_pure_dipole_is_not_a_translation()
    assert r["the_closed_form_is_confirmed"] is True
    assert r["worst_relative_error"] < 1e-6


# ── and translation invariance holds anyway ─────────────────────────────────
def test_the_translation_family_carries_induced_monopole_and_quadrupole():
    f = translation_family(0.1, radius=1.0)
    assert set(f) == {0, 1, 2}
    assert f[1] == pytest.approx(0.1)
    assert f[0] == pytest.approx(-f[2])            # equal and opposite
    assert f[2] == pytest.approx(0.01 / 3.0, rel=1e-12)


def test_translation_invariance_is_exact_along_the_right_family():
    r = measure_translation_invariance_is_exact()
    assert r["the_exact_sphere_is_area_preserving"] is True
    assert r["the_truncated_family_is_preserving_to_order_d4"] is True
    assert r["the_pure_dipole_is_not"] is True
    for row in r["rows"]:
        assert row["exact_sphere"] < 1e-12
        assert row["second_order_family"] < row["pure_P1"]


def test_the_family_error_falls_faster_than_the_dipole_error():
    """`O(d⁴)` against `O(d²)`, so the gap widens as `d` shrinks."""
    r = measure_translation_invariance_is_exact(displacements=(0.02, 0.10))
    small, large = r["rows"][0], r["rows"][1]
    assert small["family_beats_pure_by"] > large["family_beats_pure_by"]


# ── the second half of the answer ───────────────────────────────────────────
def test_the_coupling_is_screened_geometrically():
    r = measure_the_coupling_is_screened_by_separation()
    assert r["every_shape_mode_couples"] is True
    assert r["but_it_falls_geometrically"] is True
    assert r["suppression_from_ell_1_to_ell_8"] > 100.0


def test_the_suppression_is_b_over_a_to_the_ell():
    pair = ShellPair(b=2.0, a=5.0)
    for ell in (0, 1, 2, 5):
        assert pair.coupling_ratio(ell) == pytest.approx(0.4 ** ell, rel=1e-14)


# ── the constitutive gap ────────────────────────────────────────────────────
def test_a_fluid_shell_has_no_shear_response():
    r = measure_a_fluid_shell_has_no_shear_modulus()
    assert r["a_fluid_shell_has_no_shear_response"] is True
    assert r["an_elastic_one_needs_an_extra_input"] is True
    assert r["the_shear_modulus_is_never_fitted"] is True
    for row in r["rows"]:
        assert row["fluid_shear_stiffness"] == 0.0
        assert row["elastic_shear_stiffness"] > 0.0


def test_the_shear_stiffness_vanishes_at_the_modes_it_should():
    """Shape response cannot restore a dilation or a translation."""
    elastic = ShellPair(shear_modulus=1.0)
    assert elastic.shear_stiffness(1, 1.0) == 0.0        # translation
    assert elastic.shear_stiffness(2, 1.0) > 0.0
    assert elastic.shear_stiffness(0, 1.0) < 0.0         # (ℓ−1)(ℓ+2) = −2


# ── the radial structure ────────────────────────────────────────────────────
def test_the_interior_exponent_is_ell_and_the_exterior_carries_the_dimension():
    for d in (4, 5, 6):
        assert transfer_exponent(3, dim=d) == 3.0
        assert transfer_exponent(3, dim=d, outward=True) == -(3 + d - 3)


def test_bad_arguments_are_refused():
    with pytest.raises(ValueError):
        transfer_exponent(-1)
    with pytest.raises(ValueError):
        transfer_exponent(2, dim=3)
    with pytest.raises(ValueError):
        area_second_variation(-2)
    with pytest.raises(ValueError):
        mutual_stiffness(2, 5.0, 2.0)          # b must be inside a
    with pytest.raises(ValueError):
        ShellPair(b=5.0, a=2.0)
