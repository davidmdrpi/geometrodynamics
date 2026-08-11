"""
Tests for the converging ring, the gap it can reach across, and the fold.

Two questions that look like one, kept apart here because they answer
differently:

* **reaching** is set by the gap — threshold ``δ / max|u|``, and shrinking the
  gap buys the converging ring a *lead* on the focal pulse;
* **crossing** is not, and cannot be.  A curve ``r = f(σ)`` with ``f``
  single-valued is embedded, so it never self-intersects at any gap or any
  gain.  The intersection detector is validated against shapes that do cross,
  so the negative result is not a broken counter.

What crosses is a **fold**, which needs tangential freedom.  Its threshold is
set by the front's curvature, is independent of the gap, and scales as the
square of the pulse width.  Folding is necessary but not sufficient: a crossing
always comes with a fold, a fold need not cross.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.viz.circle_slice import ANTIPODAL_TIME, TWO_PI
from geometrodynamics.viz.slice_folding import (
    LagrangianSlice,
    measure_a_graph_never_intersects,
    measure_folding_ignores_the_gap,
    measure_folding_is_necessary_for_crossing,
    measure_the_fold_threshold,
    measure_the_fold_threshold_converges,
    measure_the_fold_threshold_scales_with_the_pulse,
    measure_the_intersection_test_can_see_one,
    measure_the_ring_grows_as_it_converges,
    measure_when_the_ring_spans_the_gap,
    self_intersections,
)
from geometrodynamics.viz.warped_sphere import NestedShells


# ── the detector, before anything is concluded from it ──────────────────────
def test_the_intersection_detector_works():
    r = measure_the_intersection_test_can_see_one()
    assert r["circle"] == 0
    assert r["limacon_with_an_inner_loop"] > 0
    assert r["figure_eight"] > 0
    assert r["a_folded_loop"] > 0
    assert r["the_detector_works"] is True


def test_the_detector_ignores_neighbouring_segments():
    """Adjacent segments share an endpoint; that is not a crossing."""
    sig = np.linspace(0.0, TWO_PI, 40, endpoint=False)
    square = np.stack([np.cos(sig), np.sin(sig)], axis=-1)
    assert self_intersections(square) == 0
    assert self_intersections(square, closed=False) == 0


# ── the ring ────────────────────────────────────────────────────────────────
def test_the_ring_thins_then_grows():
    r = measure_the_ring_grows_as_it_converges(frames=160, n_radial=1200)
    assert r["the_ring_thins_then_grows"] is True
    assert r["growth_from_the_equator"] > 3.0
    assert r["equator_distance"] == pytest.approx(0.5 * math.pi, abs=0.25)


def test_the_ring_follows_the_focusing_law():
    """``1/√(sin d)`` — a closing ring, not a 1-D pulse."""
    r = measure_the_ring_grows_as_it_converges(frames=160, n_radial=1200)
    assert r["follows_one_over_root_sin"] is True
    assert r["law_ratio_mean"] == pytest.approx(1.0, abs=0.25)
    assert r["law_ratio_spread"] < 0.25


# ── reaching across ─────────────────────────────────────────────────────────
def test_the_spanning_threshold_is_the_gap_over_the_peak():
    r = measure_when_the_ring_spans_the_gap(
        deltas=(0.26, 0.09), gains=(0.05, 0.20, 0.80), frames=120)
    assert r["threshold_is_delta_over_peak"] is True


def test_shrinking_the_gap_buys_dwell_and_lead():
    """The reaching stops being an instant and becomes a state."""
    r = measure_when_the_ring_spans_the_gap(
        deltas=(0.26, 0.09), gains=(0.40, 0.80), frames=120)
    assert r["shrinking_the_gap_buys_dwell"] is True
    assert r["shrinking_the_gap_buys_lead"] is True
    assert r["the_converging_ring_spans_well_before_the_pulse"] is True
    assert r["longest_lead_before_the_focus"] > 0.5


# ── ...and still never crossing ─────────────────────────────────────────────
def test_a_graph_never_intersects_at_any_gap_or_gain():
    r = measure_a_graph_never_intersects(
        deltas=(0.26, 0.05), gains=(1.0, 8.0), frames=30, n_sigma=601)
    assert r["a_graph_is_embedded"] is True
    assert r["worst_self_intersections"] == 0
    assert r["and_it_was_genuinely_driven"] is True
    assert r["most_seam_crossings"] > 5
    for row in r["rows"]:
        assert row["self_intersections"] == 0


# ── the fold ────────────────────────────────────────────────────────────────
def test_the_lagrangian_slice_is_a_graph_when_the_mixing_is_zero():
    """``λ = 0`` is exactly the height field — embedded however hard it is driven."""
    s = LagrangianSlice(mixing=0.0, n_sigma=601, n_radial=1200)
    s.reset()
    s.advance_to(3.0)
    for gain in (0.1, 2.0, 30.0):
        assert s.is_folded(gain=gain) is False
        assert self_intersections(s.points(gain=gain)) == 0
        assert np.allclose(s.sigma(gain=gain), s.sigma0)


def test_the_displacement_is_periodic():
    """Otherwise the closing chord fakes a crossing where the map never folded.

    ``σ = −π`` and ``σ = +π`` are the same point, so the sample set excludes one
    of them and the derivatives wrap.  Getting this wrong reported an
    intersection at times when the Jacobian was positive everywhere.
    """
    s = LagrangianSlice(mixing=1.0, n_sigma=801, n_radial=1200)
    s.reset()
    s.advance_to(3.0)
    assert len(s.sigma0) == 801
    assert float(s.sigma0[-1]) < math.pi
    du = s._du()
    # ``u`` is even in σ, so ∂_σu is odd.  Index 0 sits exactly on σ = −π (the
    # antipode) and reflection σ → −σ maps index k to −k mod n; σ = 0 itself is
    # not a sample point, which is why the check is oddness rather than a value.
    assert abs(float(du[0])) < 1e-8
    n = len(du)
    k = np.arange(1, n)
    assert np.allclose(du[k], -du[(-k) % n], atol=1e-9)
    total = float(np.sum(du)) * s.dsigma
    assert abs(total) < 1e-9              # a periodic derivative integrates to 0


def test_the_fold_threshold_is_front_curvature():
    r = measure_the_fold_threshold(frames=90, n_sigma=801, n_radial=1200)
    assert r["the_threshold_is_front_curvature"] is True
    assert r["relative_error"] < 1e-6
    assert r["past_it_the_curve_self_intersects"] is True
    assert r["below_it_the_curve_does_not"] is True


def test_it_folds_on_the_converging_ring():
    """Where ``∂²_σu`` peaks — which is exactly the steep converging front."""
    r = measure_the_fold_threshold(frames=90, n_sigma=801, n_radial=1200)
    assert r["it_folds_at_the_converging_ring"] is True
    assert r["distance_to_the_antipode"] < 0.3
    assert r["folds_first_at_time"] == pytest.approx(ANTIPODAL_TIME, abs=0.3)


def test_the_fold_threshold_converges_under_refinement():
    """``∂²_σu`` at a focusing front is what a coarse grid gets wrong."""
    r = measure_the_fold_threshold_converges(
        grids=((801, 1000), (1601, 2000)), frames=80)
    assert r["it_converges"] is True
    assert r["last_step_drift"] < 0.05


def test_a_crossing_always_folds_but_a_fold_need_not_cross():
    r = measure_folding_is_necessary_for_crossing(frames=45, n_sigma=601)
    assert r["a_crossing_always_folds"] is True
    assert r["crossing_without_fold"] == 0
    assert r["nothing_crosses_below_threshold"] is True


# ── the two knobs ───────────────────────────────────────────────────────────
def test_the_fold_threshold_does_not_know_the_gap_exists():
    """The load-bearing separation: shrinking the vacuole cannot buy a crossing."""
    r = measure_folding_ignores_the_gap(deltas=(0.26, 0.12, 0.05), frames=60)
    assert r["fold_threshold_is_gap_independent"] is True
    assert r["fold_threshold_spread"] < 1e-9
    assert r["span_threshold_scales_with_the_gap"] is True


def test_the_fold_threshold_scales_with_the_pulse_squared():
    r = measure_the_fold_threshold_scales_with_the_pulse(
        widths=(0.36, 0.18, 0.12), frames=80, n_sigma=1201, n_radial=1500)
    assert r["narrower_folds_sooner"] is True
    assert r["scales_as_the_width_squared"] is True
    ratios = [row["product_over_width_squared"] for row in r["rows"]]
    assert max(ratios) - min(ratios) < 0.12 * (sum(ratios) / len(ratios))


# ── mechanics ───────────────────────────────────────────────────────────────
def test_the_jacobian_matches_the_predicted_threshold():
    s = LagrangianSlice(mixing=1.0, n_sigma=801, n_radial=1200)
    s.reset()
    s.advance_to(3.0)
    p = s.predicted_fold_product()
    assert s.fold_margin(gain=0.98 * p) > 0.0
    assert s.fold_margin(gain=1.02 * p) < 0.0


def test_points_are_finite_and_scale_with_the_gain():
    s = LagrangianSlice(mixing=1.0, n_sigma=401, n_radial=900,
                        shells=NestedShells(r_mid=1.0, delta=0.2))
    s.reset()
    s.advance_to(1.6)
    assert np.all(np.isfinite(s.points(gain=0.3)))
    assert np.all(s.radius(gain=50.0) > 0.0)        # multiplicative law
    base = s.radius(gain=0.0)
    assert np.allclose(base, s.shells.r_mid)
