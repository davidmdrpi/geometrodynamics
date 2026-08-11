"""
Tests for the great-circle slice and the glued bulk it wraps through.

The load-bearing claim here is a **negative** one, so it is tested hardest:
gluing ``R_outer`` to ``R_inner`` makes the radial direction a circle and the
curve's home a torus, but it never lets the curve wind.  The curve is a graph
``r = f(σ)`` with ``f`` single-valued, so every outward crossing of the seam is
paid for by an inward one.

Around that:

* the slice is a *cut of the 2-D solve*, not a separate 1-D wave — so it
  inherits the sphere's caustic;
* the wrap threshold is exactly the half-gap over the run's peak;
* the curve closes, which is what makes the winding number well posed;
* and different pulses differ in the *arc* they carry onto the far sheet,
  which is just their own width, and not in any topological number.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.viz.circle_slice import (
    RETURN_TIME,
    TWO_PI,
    BulkAnnulus,
    CircleSlice,
    measure_different_waves_wrap_differently,
    measure_the_curve_closes,
    measure_the_lobes_meet_at_the_antipode,
    measure_the_seam_is_crossed_in_pairs,
    measure_the_slice_is_the_sphere,
    measure_the_wrap_threshold,
)
from geometrodynamics.viz.warped_sphere import NestedShells


@pytest.fixture(scope="module")
def slc():
    return CircleSlice(pulse_width=0.18, n_sigma=721, n_radial=900)


# ── the gluing itself ───────────────────────────────────────────────────────
def test_the_bulk_glues_its_two_boundaries():
    b = BulkAnnulus()
    assert b.gap == pytest.approx(b.r_outer - b.r_inner)
    # just inside the outer boundary stays on sheet 0...
    drawn, sheet = b.wrap(np.array([b.r_outer - 1e-9]))
    assert int(sheet[0]) == 0
    # ...and just past it reappears at the inner one, on sheet 1
    drawn, sheet = b.wrap(np.array([b.r_outer + 1e-9]))
    assert int(sheet[0]) == 1
    assert float(drawn[0]) == pytest.approx(b.r_inner, abs=1e-8)
    # below the inner boundary lands on sheet −1, near the outer one
    drawn, sheet = b.wrap(np.array([b.r_inner - 1e-9]))
    assert int(sheet[0]) == -1
    assert float(drawn[0]) == pytest.approx(b.r_outer, abs=1e-8)


def test_wrapping_is_periodic_in_the_gap():
    b = BulkAnnulus()
    r = np.array([1.03, 0.81, 1.19])
    d0, s0 = b.wrap(r)
    d1, s1 = b.wrap(r + 3.0 * b.gap)
    assert np.allclose(d0, d1)
    assert np.all(s1 - s0 == 3)


def test_the_drawn_radius_always_lands_in_the_annulus():
    b = BulkAnnulus(NestedShells(r_mid=1.0, delta=0.3))
    r = np.linspace(-4.0, 6.0, 5000)
    drawn, _ = b.wrap(r)
    assert np.all(drawn >= b.r_inner - 1e-12)
    assert np.all(drawn < b.r_outer + 1e-12)


# ── the slice is the sphere ─────────────────────────────────────────────────
def test_the_slice_is_a_cut_of_the_two_dimensional_solve(slc):
    """Not a 1-D wave on a circle — so the caustic is the sphere's."""
    r = measure_the_slice_is_the_sphere(slc)
    assert r["is_a_cut_of_the_sphere"] is True
    assert r["worst_mismatch_against_the_sphere"] < 1e-9
    assert r["both_arms_are_the_same"] is True


def test_the_distance_along_the_slice_is_the_absolute_angle(slc):
    d = slc.distance()
    assert np.allclose(d, np.abs(slc.sigma))
    assert float(d[0]) == pytest.approx(math.pi)      # σ = −π is the antipode
    assert float(d[len(d) // 2]) == pytest.approx(0.0, abs=1e-12)


def test_the_lobes_meet_at_the_antipode(slc):
    r = measure_the_lobes_meet_at_the_antipode(slc)
    assert r["lobes_are_symmetric"] is True
    assert r["arrives_on_time"] is True
    assert r["final_left_lobe"] == pytest.approx(-r["final_right_lobe"], abs=1e-9)


# ── the threshold ───────────────────────────────────────────────────────────
def test_the_wrap_threshold_is_the_half_gap_over_the_peak(slc):
    r = measure_the_wrap_threshold(slc, frames=120)
    assert r["threshold_is_the_half_gap_over_the_peak"] is True
    assert r["relative_error"] < 1e-3
    assert r["display_gain_is_below_threshold"] is True


def test_below_the_threshold_nothing_wraps(slc):
    slc.reset()
    slc.advance_to(2.9)
    assert slc.seam_crossings(gain=slc.gain)["unsigned"] == 0
    assert slc.excursion(gain=slc.gain)["wraps"] is False


# ── the topology ────────────────────────────────────────────────────────────
def test_the_curve_closes_wrapped_or_not(slc):
    """What makes the winding number a well-posed integer in the first place."""
    r = measure_the_curve_closes(slc)
    assert r["closes_exactly"] is True
    assert r["below_endpoint_gap"] < 1e-12
    assert r["above_endpoint_gap"] < 1e-12


def test_the_seam_is_crossed_in_pairs():
    """The load-bearing negative result: amplitude never buys charge."""
    r = measure_the_seam_is_crossed_in_pairs(
        CircleSlice(pulse_width=0.18, n_sigma=721, n_radial=900),
        gains=(1.6, 3.6, 5.0), frames=100)
    assert r["amplitude_buys_crossings"] is True
    assert r["most_unsigned_crossings"] >= 2
    assert r["the_signed_count_is_always_zero"] is True
    assert r["a_height_field_cannot_wind"] is True
    for row in r["rows"]:
        assert row["signed"] == 0
        assert row["winding"] == 0
        assert row["unsigned"] % 2 == 0          # in and out, always


def test_a_single_valued_height_has_degree_zero(slc):
    """Stated directly, at absurd gains, so the reason is on the record.

    ``ρ(σ)`` comes from a function on the circle, so its degree as a map
    ``S¹ → S¹`` is zero however violently it is driven.
    """
    slc.reset()
    slc.advance_to(2.98)
    for gain in (0.2, 1.0, 5.0, 25.0, 200.0):
        assert slc.winding_number(gain=gain) == 0
        assert slc.seam_crossings(gain=gain)["signed"] == 0


def test_the_crossing_ledger_and_the_degree_agree(slc):
    """Two independent computations of the same integer."""
    slc.reset()
    for i in range(24):
        slc.advance_to((i + 1) * RETURN_TIME / 24)
        for gain in (0.4, 1.5, 6.0):
            assert (slc.seam_crossings(gain=gain)["signed"]
                    == slc.winding_number(gain=gain))


# ── different waves ─────────────────────────────────────────────────────────
def test_different_waves_differ_in_arc_not_in_winding():
    r = measure_different_waves_wrap_differently(widths=(0.36, 0.14, 0.08),
                                                 gain=0.80, frames=120)
    assert r["none_of_them_winds"] is True
    assert r["wide_pulses_carry_a_broader_arc"] is True
    assert r["arc_ratio"] > 2.0
    assert r["the_far_sheet_arc_is_the_pulse"] is True


def test_the_far_sheet_arc_scales_with_the_pulse_width():
    r = measure_different_waves_wrap_differently(widths=(0.36, 0.14, 0.08),
                                                 gain=0.80, frames=120)
    ratios = [row["arc_over_pulse_width"] for row in r["rows"]]
    assert max(ratios) - min(ratios) < 0.35 * (sum(ratios) / len(ratios))


# ── drawing ─────────────────────────────────────────────────────────────────
def test_segments_cover_the_curve_and_break_only_at_the_seam(slc):
    slc.reset()
    slc.advance_to(2.98)
    thr = (slc.bulk.r_outer - slc.bulk.r_mid) / float(np.max(np.abs(slc.field())))
    gain = 3.6 * thr
    segs = slc.segments(gain=gain)
    crossings = slc.seam_crossings(gain=gain)["unsigned"]
    assert len(segs) == crossings + 1
    assert sum(len(s) for s in segs) == slc.n_sigma
    for s in segs:
        assert np.all(np.isfinite(s))
        rad = np.linalg.norm(s, axis=1)
        assert np.all(rad >= slc.bulk.r_inner - 1e-9)
        assert np.all(rad <= slc.bulk.r_outer + 1e-9)


def test_an_unwrapped_slice_is_a_single_segment(slc):
    slc.reset()
    slc.advance_to(1.2)
    segs = slc.segments(gain=slc.gain)
    assert len(segs) == 1
    assert len(segs[0]) == slc.n_sigma


def test_the_unrolled_coordinates_are_the_drawn_ones(slc):
    slc.reset()
    slc.advance_to(1.9)
    sigma, rho = slc.unrolled(gain=slc.gain)
    assert np.allclose(sigma, slc.sigma)
    assert np.all((rho >= 0.0) & (rho < 1.0 + 1e-12))
    drawn = slc.bulk.r_inner + rho * slc.bulk.gap
    assert np.allclose(np.linalg.norm(slc.points(gain=slc.gain), axis=1), drawn)


def test_the_deformation_scales_with_the_gain(slc):
    slc.reset()
    slc.advance_to(1.5)
    base = slc.radius(gain=0.0)
    assert np.allclose(base, slc.bulk.r_mid)
    one = slc.radius(gain=1e-3) - base
    two = slc.radius(gain=2e-3) - base
    assert np.allclose(two, 2.0 * one, rtol=1e-9)
