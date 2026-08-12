"""
Tests for Israel junctions, and for keeping three observables apart.

The module exists because "does the detached shell replace the exotic matter?"
is three questions wearing one coat, so most of these tests are about the three
answers being allowed to disagree:

* the shell's own ``σ`` — and the theorem that an anti-aligned one is always
  negative, which is an identity and is asserted as one;
* the force it puts on the throat, which is positive and comes from ordinary
  matter;
* the stiffness, which can be positive for both modes and still requires the
  throat to carry a pathological equation of state.

Also held down: the controls reproduce published shells before anything new is
asked, and the vanishing off-diagonal is recorded as structural rather than
measured, because Birkhoff is imported the moment the region between the
surfaces is written as Schwarzschild with a constant mass.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.shells.junction import (
    ALIGNED,
    ANTI_ALIGNED,
    DetachedShell,
    Region,
    ThroatShellSystem,
    Z2Throat,
    measure_any_minimal_surface_is_exotic,
    measure_the_detached_shell_can_be_ordinary,
    measure_the_hessian_has_no_flat_direction,
    measure_the_junction_reproduces_known_shells,
    measure_the_shell_force_on_the_throat,
    measure_the_stability_window,
    measure_the_three_observables_are_independent,
    measure_the_throat_and_shell_are_decoupled,
    surface_stress,
)


# ── the controls ────────────────────────────────────────────────────────────
def test_a_bubble_carries_ordinary_matter_of_the_right_mass():
    r = 10.0
    s = surface_stress(Region(), Region(mass=0.001), r,
                       eps_inner=+1, eps_outer=+1)
    assert s["sigma"] > 0.0
    assert s["is_exotic"] is False
    assert s["rest_mass"] == pytest.approx(0.001, rel=1e-3)


def test_a_z2_throat_reproduces_visser():
    reg = Region(mass=0.5)
    s = surface_stress(reg, reg, 3.0, eps_inner=-1, eps_outer=+1)
    ref = -math.sqrt(reg.f(3.0)) / (2.0 * math.pi * 3.0)
    assert s["sigma"] == pytest.approx(ref, abs=1e-15)
    assert s["is_exotic"] is True


def test_the_controls_measurement_agrees():
    r = measure_the_junction_reproduces_known_shells()
    assert r["visser_is_reproduced"] is True
    assert r["the_bubble_is_ordinary"] is True
    assert r["worst_visser_error"] < 1e-12


def test_a_shell_inside_a_horizon_is_rejected():
    with pytest.raises(ValueError):
        Region(mass=10.0).beta(3.0)


# ── observable 1: the exoticity theorem ─────────────────────────────────────
def test_an_anti_aligned_shell_is_exotic_for_any_bulk():
    """``σ = −(β₊+β₋)/4πGR`` — the two terms add, so the sign is forced."""
    cases = [(Region(mass=0.3), Region(mass=0.9)),
             (Region(mass=0.9), Region(mass=0.3)),
             (Region(), Region()),
             (Region(lambda_=0.001), Region(mass=0.2, charge=0.4))]
    for inner, outer in cases:
        s = DetachedShell(inner, outer, 8.0, ANTI_ALIGNED).stress()
        assert s["sigma"] < 0.0
        assert s["is_exotic"] is True


def test_the_orientation_is_not_a_relabelling():
    """Same regions, same radius; only the gluing differs, and so does σ."""
    inner, outer = Region(mass=0.7), Region(mass=1.0)
    al = DetachedShell(inner, outer, 20.0, ALIGNED).stress()
    anti = DetachedShell(inner, outer, 20.0, ANTI_ALIGNED).stress()
    assert al["sigma"] > 0.0 > anti["sigma"]
    assert al["orientation"] == +1 and anti["orientation"] == -1


def test_the_sweep_finds_no_counterexample_and_could_have():
    r = measure_any_minimal_surface_is_exotic(samples=20_000)
    assert r["every_minimal_surface_is_exotic"] is True
    assert r["anti_aligned_positive_sigma"] == 0
    assert r["worst_anti_aligned_sigma"] < 0.0
    # the control: the same sweep does find positive σ when one exists
    assert r["the_sweep_can_find_positive_sigma"] is True
    assert 0.3 < r["aligned_positive_fraction"] < 0.7


def test_a_throat_cannot_be_configured_non_exotic():
    """A throat is a minimal surface, so the theorem applies to it too."""
    for mass in (0.1, 0.7, 1.0):
        for b0 in (3.0, 5.0, 12.0):
            assert Z2Throat(mass=mass, b0=b0).static()["is_exotic"] is True


def test_the_aligned_shell_can_be_ordinary():
    r = measure_the_detached_shell_can_be_ordinary()
    assert r["the_aligned_shell_is_ordinary"] is True
    assert r["the_anti_aligned_shell_is_exotic"] is True


# ── observable 2: the force ─────────────────────────────────────────────────
def test_the_shell_pushes_the_throat_outward():
    r = measure_the_shell_force_on_the_throat()
    assert r["the_force_opposes_closure"] is True
    assert r["it_grows_with_the_screened_mass"] is True
    assert r["it_matches_2_G_dM_over_b_squared"] is True


def test_no_shell_means_no_force():
    sysm = ThroatShellSystem(adm_mass=1.0, screened=0.0, b0=5.0)
    assert abs(sysm.shell_force()) < 1e-12
    assert abs(sysm.shell_potential_shift(5.0)) < 1e-14


# ── observable 3: the stiffness ─────────────────────────────────────────────
def test_the_throat_stiffness_is_linear_in_beta_squared():
    t = Z2Throat(mass=0.7, b0=5.0)
    a, b, c = t.stiffness(-1.0), t.stiffness(0.0), t.stiffness(1.0)
    assert (b - a) == pytest.approx(c - b, rel=1e-10)


def test_screening_enlarges_the_window_but_never_reaches_ordinary_matter():
    r = measure_the_stability_window()
    assert r["screening_raises_the_critical_beta2"] is True
    assert r["the_window_always_needs_negative_beta2"] is True
    assert r["the_throat_is_always_exotic"] is True
    for row in r["rows"]:
        assert row["beta2_critical"] < 0.0


def test_both_normal_modes_can_be_positive_at_once():
    sysm = ThroatShellSystem(adm_mass=1.0, screened=0.5, b0=5.0, a0=12.0)
    modes = sysm.normal_modes(beta2_throat=-2.0, beta2_shell=0.5)
    assert modes["both_positive"] is True
    assert all(v > 0.0 for v in modes["eigenvalues"])


# ── the three answers are allowed to disagree ───────────────────────────────
def test_the_three_observables_disagree_on_one_system():
    r = measure_the_three_observables_are_independent()
    assert r["observable_1_shell_is_exotic"] is False    # shell is ordinary
    assert r["observable_2_supports_the_throat"] is True  # and it supports
    assert r["the_throat_is_still_exotic"] is True        # and yet
    assert r["they_do_not_agree"] is True


# ── the decoupling, and what it is worth ────────────────────────────────────
def test_moving_the_shell_changes_it_a_lot_and_the_throat_not_at_all():
    r = measure_the_throat_and_shell_are_decoupled()
    assert r["shell_sigma_varies_by"] > 2.0        # genuinely different shells
    assert r["throat_sigma_spread"] == 0.0         # to the last bit
    assert r["the_off_diagonal_vanishes"] is True


def test_the_vanishing_off_diagonal_is_labelled_structural():
    """It follows from Birkhoff, which is imported — so it is not evidence."""
    r = measure_the_throat_and_shell_are_decoupled()
    note = r["but_that_is_structural_not_measured"]
    assert "Birkhoff is assumed" in note
    assert "proves nothing more" in note


def test_the_spectrum_has_no_flat_direction_under_dilation():
    r = measure_the_hessian_has_no_flat_direction()
    assert r["eigenvalues_scale_as_inverse_length_squared"] is True
    assert r["no_flat_direction"] is True
    assert r["smallest_absolute_eigenvalue"] > 1e-10


def test_the_hessian_is_symmetric():
    sysm = ThroatShellSystem(adm_mass=1.0, screened=0.3, b0=5.0, a0=20.0)
    H = sysm.hessian(-1.0, 0.5)
    assert H[0, 1] == pytest.approx(H[1, 0], abs=1e-15)
    assert np.all(np.isfinite(H))


# ── guards ──────────────────────────────────────────────────────────────────
def test_an_unknown_orientation_is_refused():
    with pytest.raises(ValueError):
        DetachedShell(Region(), Region(mass=1.0), 10.0, orientation=0)


def test_a_throat_inside_its_horizon_is_refused():
    with pytest.raises(ValueError):
        Z2Throat(mass=5.0, b0=2.0)
