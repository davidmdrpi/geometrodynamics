"""
Tests for Darmois–Israel junctions, and for the scope of what they force.

Everything here is inside **Einstein gravity, thin shells, spherical symmetry,
vacuum bulk**, with the dimension a parameter.  The tests are mostly about
three things not being allowed to blur together:

* the shell's own ``σ``, whose sign is *forced* by the retained branches for two
  of the four gluings and free for the other two;
* the shell's effect on the throat, which is a potential **gradient** and not an
  equilibrium-consistent force;
* the stiffness, which is a stiffness and not a normal-mode frequency, because
  no kinetic metric has been derived.

Also held down: ``ε`` is derived from the gluing rather than supplied, the
forced signs survive ``D = 4, 5, 6``, and the vanishing off-diagonal is recorded
as structural rather than measured.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.shells.junction import (
    ANTI_BUBBLE,
    GLUINGS,
    INNER,
    MAXIMAL_SURFACE,
    MINIMAL_SURFACE,
    ORDINARY_BUBBLE,
    OUTER,
    DetachedShell,
    Gluing,
    Region,
    ThroatShellSystem,
    Z2Throat,
    measure_the_detached_shell_can_be_ordinary,
    measure_the_forced_signs_hold_in_any_dimension,
    measure_the_gluing_fixes_the_sign,
    measure_the_junction_reproduces_known_shells,
    measure_the_shell_potential_gradient,
    measure_the_stability_window,
    measure_the_stiffnesses_scale_dimensionally,
    measure_the_three_observables_are_independent,
    measure_the_throat_and_shell_are_decoupled,
    surface_stress,
)


# ── the controls ────────────────────────────────────────────────────────────
def test_a_bubble_carries_ordinary_matter_of_the_right_mass():
    r = 10.0
    s = surface_stress(Region.from_mass(0.0), Region.from_mass(0.001),
                       ORDINARY_BUBBLE, r)
    assert s["sigma"] > 0.0 and s["is_exotic"] is False
    assert 4.0 * math.pi * r * r * s["sigma"] == pytest.approx(0.001, rel=1e-3)


def test_a_z2_throat_reproduces_visser():
    reg = Region.from_mass(0.5)
    s = surface_stress(reg, reg, MINIMAL_SURFACE, 3.0)
    ref = -math.sqrt(reg.f(3.0)) / (2.0 * math.pi * 3.0)
    assert s["sigma"] == pytest.approx(ref, abs=1e-15)


def test_the_controls_measurement_agrees():
    r = measure_the_junction_reproduces_known_shells()
    assert r["visser_is_reproduced"] is True
    assert r["worst_visser_error"] < 1e-12


# ── ε is derived, not supplied ──────────────────────────────────────────────
def test_epsilon_follows_from_the_retained_branches():
    assert (MINIMAL_SURFACE.eps_minus, MINIMAL_SURFACE.eps_plus) == (-1, +1)
    assert (MAXIMAL_SURFACE.eps_minus, MAXIMAL_SURFACE.eps_plus) == (+1, -1)
    assert (ORDINARY_BUBBLE.eps_minus, ORDINARY_BUBBLE.eps_plus) == (+1, +1)
    assert (ANTI_BUBBLE.eps_minus, ANTI_BUBBLE.eps_plus) == (-1, -1)


def test_eta_alone_does_not_decide_the_sign():
    """The correction: ``η = −1`` covers two gluings with opposite signs."""
    assert MINIMAL_SURFACE.eta == MAXIMAL_SURFACE.eta == -1
    assert MINIMAL_SURFACE.forced_sign == -1
    assert MAXIMAL_SURFACE.forced_sign == +1
    assert ORDINARY_BUBBLE.forced_sign is None
    assert ANTI_BUBBLE.forced_sign is None
    r = measure_the_gluing_fixes_the_sign()
    assert r["eta_alone_does_not_decide"] is True
    assert r["every_forced_sign_is_realised"] is True


def test_there_are_four_gluings_and_a_bad_branch_is_refused():
    assert len(GLUINGS) == 4
    assert len({(g.minus, g.plus) for g in GLUINGS}) == 4
    with pytest.raises(ValueError):
        Gluing("SIDEWAYS", OUTER)


# ── the forced signs ────────────────────────────────────────────────────────
def test_a_minimal_surface_is_exotic_and_a_maximal_one_is_not():
    for minus, plus in ((Region(mu=0.3), Region(mu=0.9)),
                        (Region(mu=0.9), Region(mu=0.3)),
                        (Region(), Region()),
                        (Region(lambda_=0.001), Region(mu=0.2, charge=0.4))):
        assert DetachedShell(minus, plus, 8.0,
                             MINIMAL_SURFACE).stress()["sigma"] < 0.0
        assert DetachedShell(minus, plus, 8.0,
                             MAXIMAL_SURFACE).stress()["sigma"] > 0.0


def test_the_forced_signs_survive_five_dimensions():
    r = measure_the_forced_signs_hold_in_any_dimension(dims=(4, 5),
                                                       samples=6000)
    assert r["no_violations_in_any_dimension"] is True
    for row in r["rows"]:
        assert row["minimal_surface_violations"] == 0
        assert row["maximal_surface_violations"] == 0
        assert row["worst_minimal_sigma"] < 0.0 < row["worst_maximal_sigma"]


def test_a_throat_cannot_be_configured_non_exotic():
    for mu in (0.2, 1.4, 2.0):
        for b0 in (3.0, 5.0, 12.0):
            assert Z2Throat(mu=mu, b0=b0).static()["is_exotic"] is True
    assert Z2Throat(mu=4.0, b0=5.0, dim=5).static()["is_exotic"] is True


def test_the_non_exotic_gluing_is_the_disconnected_one():
    r = measure_the_detached_shell_can_be_ordinary()
    assert r["the_minimal_surface_is_exotic"] is True
    assert r["the_maximal_surface_is_ordinary"] is True
    assert r["but_it_is_disconnected"] is True
    assert MINIMAL_SURFACE.connects_to_infinity_on_both_sides is True
    assert MAXIMAL_SURFACE.connects_to_infinity_on_both_sides is False


# ── the dimension is a parameter ────────────────────────────────────────────
def test_the_stress_carries_the_dimension_in_the_trace_weight():
    reg = Region(mu=0.0)
    ref = surface_stress(reg, reg, MINIMAL_SURFACE, 5.0)["sigma"]
    for d in (5, 6):
        flat = Region(mu=0.0, dim=d)
        s = surface_stress(flat, flat, MINIMAL_SURFACE, 5.0)["sigma"]
        assert s == pytest.approx(ref * (d - 2) / 2.0, rel=1e-12)


def test_mixed_dimensions_are_refused():
    with pytest.raises(ValueError):
        surface_stress(Region(dim=4), Region(dim=5), ORDINARY_BUBBLE, 5.0)
    with pytest.raises(ValueError):
        Region(dim=3)
    with pytest.raises(ValueError):
        Region.from_mass(1.0, dim=5)


# ── the static solve is generic across gluings ──────────────────────────────
def test_every_gluing_gets_a_static_solution_with_zero_residual():
    minus, plus = Region(mu=1.0), Region(mu=2.0)
    for g in GLUINGS:
        st = DetachedShell(minus, plus, 12.0, g).static()
        assert abs(st["residual_Vprime"]) < 1e-8
        assert np.isfinite(st["sigma"]) and np.isfinite(st["p"])


def test_the_stiffness_matches_the_independently_derived_value():
    """Regression against the hand-derived, RK4-verified shell stiffness."""
    sh = DetachedShell(Region.from_mass(0.5), Region.from_mass(1.0), 12.0,
                       ORDINARY_BUBBLE)
    assert sh.static()["p"] == pytest.approx(1.0647e-05, rel=1e-3)
    assert sh.stiffness(0.5) == pytest.approx(0.0223298686, rel=1e-6)


def test_the_throat_stiffness_is_linear_in_beta_squared():
    t = Z2Throat(mu=1.4, b0=5.0)
    a, b, c = t.stiffness(-1.0), t.stiffness(0.0), t.stiffness(1.0)
    assert (b - a) == pytest.approx(c - b, rel=1e-6)


# ── observable 2, as a gradient ─────────────────────────────────────────────
def test_the_screening_gradient_opposes_closure_and_is_labelled_as_such():
    r = measure_the_shell_potential_gradient()
    assert r["the_gradient_opposes_closure"] is True
    assert r["it_grows_with_the_screened_mass"] is True
    assert "not an equilibrium-consistent force" in \
        r["it_is_not_an_equilibrium_consistent_force"] or \
        "no equation-of-state response" in \
        r["it_is_not_an_equilibrium_consistent_force"]
    for row in r["rows"]:
        assert row["acceleration_contribution"] == pytest.approx(
            0.5 * row["minus_dV_db"], rel=1e-9, abs=1e-15)


def test_no_shell_means_no_gradient():
    sysm = ThroatShellSystem(mu_outer=2.0, screened=0.0, b0=5.0)
    assert abs(sysm.shell_potential_gradient()["minus_dV_db"]) < 1e-12


# ── observable 3, as a stiffness ────────────────────────────────────────────
def test_screening_enlarges_the_window_but_never_reaches_positive_beta2():
    r = measure_the_stability_window()
    assert r["screening_raises_the_critical_beta2"] is True
    assert r["the_window_always_needs_negative_beta2"] is True
    assert r["beta2_is_an_eos_derivative_not_a_sound_speed"] is True
    for row in r["rows"]:
        assert row["beta2_critical"] < 0.0


def test_the_stiffnesses_are_not_called_normal_modes():
    sysm = ThroatShellSystem(mu_outer=2.0, screened=1.0, b0=5.0, a0=12.0)
    st = sysm.stiffnesses(beta2_throat=-2.0, beta2_shell=0.5)
    assert "not normal-mode frequencies" in st["note"]
    assert "kinetic metric" in st["note"]
    assert st["both_stiffnesses_positive"] is True


def test_the_stiffness_matrix_is_symmetric_and_finite():
    sysm = ThroatShellSystem(mu_outer=2.0, screened=0.6, b0=5.0, a0=20.0)
    m = sysm.stiffness_matrix(-1.0, 0.5)
    assert m[0, 1] == pytest.approx(m[1, 0], abs=1e-15)
    assert np.all(np.isfinite(m))


# ── the three answers are allowed to disagree ───────────────────────────────
def test_the_three_observables_disagree_on_one_system():
    r = measure_the_three_observables_are_independent()
    assert r["observable_1_shell_is_exotic"] is False
    assert r["observable_2_opposes_closure"] is True
    assert r["the_throat_is_still_exotic"] is True
    assert r["they_do_not_agree"] is True


# ── the decoupling, and what it is worth ────────────────────────────────────
def test_moving_the_shell_changes_it_a_lot_and_the_throat_not_at_all():
    r = measure_the_throat_and_shell_are_decoupled()
    assert r["shell_sigma_varies_by"] > 2.0
    assert r["throat_sigma_spread"] == 0.0
    assert r["the_off_diagonal_vanishes"] is True


def test_the_vanishing_off_diagonal_is_labelled_structural():
    r = measure_the_throat_and_shell_are_decoupled()
    assert "Birkhoff is assumed" in r["but_that_is_structural_not_measured"]
    assert "proves nothing more" in r["but_that_is_structural_not_measured"]
    assert "in this vacuum spherical model" in r["what_it_establishes"]
    assert "not that every" in r["what_it_establishes"]


def test_the_scaling_check_is_labelled_as_units_only():
    r = measure_the_stiffnesses_scale_dimensionally()
    assert r["stiffnesses_scale_as_inverse_length_squared"] is True
    assert "does not show" in r["this_is_a_units_check_only"]
    assert "dilation mode" in r["this_is_a_units_check_only"]


# ── guards ──────────────────────────────────────────────────────────────────
def test_a_throat_inside_its_horizon_is_refused():
    with pytest.raises(ValueError):
        Z2Throat(mu=10.0, b0=2.0)


def test_a_shell_inside_a_horizon_is_refused():
    with pytest.raises(ValueError):
        Region(mu=20.0).beta(3.0)
