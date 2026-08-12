"""
Tests for the constraint audit of the height-field representation.

Three objections, three different fates, and the tests are arranged to keep
them apart:

* **causality** — the front is causal and the antipode is genuinely dark, so
  this objection does not land on the solve;
* **the permanent offset** — it does land, one step over. ``ℓ = 0`` has
  ``ω = 0``, so a Gaussian's ``w²/4`` monopole is frozen into a closed surface
  forever, and the fix has to stay inside the pulse or it just moves the
  problem to the far side;
* **height is not energy** — the constant offset displaces every point and
  carries exactly zero energy, which is the part no choice of source repairs.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.viz.wave_constraints import (
    CleanSphereSim,
    compact_monopole_free_profile,
    measure_a_localised_correction_is_the_only_one_that_works,
    measure_launch_and_focus_are_the_same,
    measure_the_clean_source_leaves_the_antipode_alone,
    measure_the_constant_monopole,
    measure_the_front_is_causal,
    measure_the_growth_is_invisible_for_most_of_the_trip,
    measure_the_monopole_carries_no_energy,
    measure_the_offset_never_leaves,
    measure_the_ring_bookkeeping,
    monopole_free_profile,
    sphere_monopole,
)


# ── causality: this objection does not land ─────────────────────────────────
def test_the_front_is_causal():
    r = measure_the_front_is_causal()
    assert r["nothing_outruns_the_front"] is True
    assert r["the_antipode_stays_dark_until_the_front"] is True
    assert abs(r["antipode_at_launch"]) < 1e-100
    assert abs(r["antipode_at_two"]) < 1e-12
    assert abs(r["antipode_at_pi"]) > 0.5


def test_the_antipode_climbs_monotonically_toward_arrival():
    """130 orders of magnitude, in step with ``t → π`` — not a leak."""
    r = measure_the_front_is_causal()
    amps = [abs(row["antipode"]) for row in r["rows"]]
    assert all(b >= a for a, b in zip(amps, amps[1:]))


# ── the permanent offset: this one does ─────────────────────────────────────
def test_a_gaussian_carries_a_monopole_of_w_squared_over_four():
    r = measure_the_constant_monopole()
    assert r["the_monopole_is_w_squared_over_four"] is True
    assert r["monopole_at_launch"] == pytest.approx(
        r["predicted_w_squared_over_four"], rel=0.03)


def test_the_monopole_never_changes():
    """``ℓ = 0`` has ``ω = 0``: nowhere to propagate, nothing to decay into."""
    r = measure_the_constant_monopole(frames=40)
    assert r["it_never_changes"] is True
    assert r["monopole_drift_over_a_full_return"] < 1e-3


def test_the_surface_never_returns_home():
    r = measure_the_offset_never_leaves(frames=200, returns=1.0)
    assert r["every_point_is_offset"] is True
    assert r["the_surface_never_returns_home"] is True
    assert r["quietest_moment"] > r["predicted_offset"]


def test_the_offset_carries_no_energy():
    """The measurable core of "a deforming surface is the wrong picture"."""
    r = measure_the_monopole_carries_no_energy()
    assert r["it_carries_no_energy"] is True
    assert r["energy_of_the_offset"] == pytest.approx(0.0, abs=1e-24)
    assert r["energy_of_the_pulse"] > 0.1


# ── the fix has to stay inside the pulse ────────────────────────────────────
def test_subtracting_the_mean_just_moves_the_problem():
    r = measure_a_localised_correction_is_the_only_one_that_works()
    assert r["subtracting_the_mean_just_moves_it"] is True
    # it removes exactly what it deposits on the far side
    assert abs(r["mean_subtracted_antipode"]) == pytest.approx(
        r["gaussian_monopole"], rel=1e-6)


def test_a_localised_difference_of_bumps_is_monopole_free():
    r = measure_a_localised_correction_is_the_only_one_that_works()
    assert r["the_localised_source_is_clean"] is True
    assert abs(r["localised_monopole"]) < 1e-12
    assert abs(r["localised_antipode"]) < 1e-12


def test_the_compact_profile_is_exactly_zero_outside_its_support():
    """Not small — zero.  A Gaussian can never say that."""
    d = np.linspace(0.0, math.pi, 4000)
    f, c, grad = compact_monopole_free_profile(d, a=0.30, k=2.0)
    assert abs(sphere_monopole(f, d)) < 1e-12
    outside = d > 0.60 + 1e-9
    assert np.all(f[outside] == 0.0)
    assert np.all(grad[outside] == 0.0)
    assert 0.0 < c < 1.0


def test_the_compact_profile_rejects_a_support_that_does_not_fit():
    d = np.linspace(0.0, math.pi, 500)
    with pytest.raises(ValueError):
        compact_monopole_free_profile(d, a=2.0, k=2.0)
    with pytest.raises(ValueError):
        compact_monopole_free_profile(d, a=0.3, k=0.5)


def test_the_compact_source_wins_on_both_counts():
    """Monopole-free *and* quieter on the far side than the original."""
    r = measure_the_clean_source_leaves_the_antipode_alone(frames=40)
    assert r["the_compact_source_wins_on_both"] is True
    assert abs(r["compact_monopole"]) < 1e-3 * abs(r["gaussian_monopole"])
    assert r["compact_far_side"] <= r["gaussian_far_side"]


def test_a_wider_gaussian_corrector_costs_far_side_quiet():
    """The honest trade the compact profile exists to avoid."""
    r = measure_the_clean_source_leaves_the_antipode_alone(frames=40)
    assert r["a_wider_corrector_costs_far_side_quiet"] is True
    assert r["wide_gaussian_far_side"] > 1e3 * r["gaussian_far_side"]


def test_the_clean_launch_reproduces_the_ordinary_one_when_asked():
    a = CleanSphereSim(n_radial=1200, monopole_free=False)
    assert abs(a.monopole() - 0.18 ** 2 / 4.0) < 0.02 * 0.18 ** 2 / 4.0
    b = CleanSphereSim(n_radial=1200, monopole_free=True)
    assert abs(b.monopole()) < 1e-4


# ── the ring: right physics, invisible in the picture ───────────────────────
def test_the_ring_conserves_its_energy():
    """``A²·(circumference) = const``, i.e. ``A ∝ 1/√(sin d)``."""
    r = measure_the_ring_bookkeeping(frames=200)
    assert r["energy_on_the_ring_is_conserved"] is True
    assert r["spread"] < 0.25
    assert r["n_samples"] > 20


def test_the_growth_is_invisible_for_most_of_the_trip():
    r = measure_the_growth_is_invisible_for_most_of_the_trip(frames=200)
    assert r["most_of_the_trip_shows_nothing"] is True
    assert r["flat_fraction_of_the_trip"] > 0.6
    assert r["growth_from_the_equator"] > 3.0


def test_the_focus_never_beats_the_launch():
    """On a compact surface the source is a focus too, so it cannot.

    This is the correction to "height should grow as the ring compresses": it
    does grow, relative to the *equator*, and still does not reach the launch.
    """
    r = measure_launch_and_focus_are_the_same(widths=(0.30, 0.12))
    assert r["the_focus_does_not_beat_the_launch"] is True
    for row in r["rows"]:
        assert row["ratio"] < 1.1
