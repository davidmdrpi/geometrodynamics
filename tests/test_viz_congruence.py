"""
Tests for the congruence, its Jacobian, and the three cases it separates.

The module exists to stop three different things being called by one name, so
the tests are mostly about keeping them apart: a dip that recovers, a caustic of
the map, and a curvature singularity — which is not reachable here and is
asserted to be out of scope rather than quietly absent.

Held down here:

* the Ricci term is *identically* zero, which is what makes the focusing purely
  shear — the Raychaudhuri residual beside it is an algebraic identity and so
  proves nothing about the integration, which the grid ladders cover instead;
* ``ḧ`` comes from the wave operator, because the seeded time difference is
  wrong by an amount that does not shrink with ``dt``;
* the launch is compactly supported, because a Gaussian tail deforms the far
  side before the front arrives and refining the grid cannot fix it;
* the neck is a ring around the antipode and never the antipode itself, which
  is spin weight and not an accident of the sampling;
* and every threshold carries the window it was measured in.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.viz.congruence import (
    ANTIPODAL_TIME,
    NECK_CAP,
    TidalCongruence,
    compact_launch,
    measure_case_three_is_unreachable,
    measure_raychaudhuri_is_exact,
    measure_the_caustic_is_a_passage,
    measure_the_caustic_thresholds,
    measure_the_front_is_causal,
    measure_the_neck_is_a_ring,
    measure_the_operator_form_matters,
    measure_the_three_cases,
    measure_the_threshold_depends_on_the_window,
    spacetime_history,
)
from geometrodynamics.viz.spin2_tidal import Spin2WaveSim


# ── the state itself ────────────────────────────────────────────────────────
def test_the_congruence_starts_undeformed():
    c = TidalCongruence(n_sigma=181)
    assert np.allclose(c.a, 1.0)
    assert np.allclose(c.b, 1.0)
    assert np.allclose(c.jacobian(), 1.0)
    assert not c.has_caustic()


def test_the_duplicated_endpoint_is_excluded():
    """``σ = ±π`` is one point of the slice; sampling both double-counts it."""
    c = TidalCongruence(n_sigma=181)
    assert c.sigma[0] == pytest.approx(-math.pi)
    assert c.sigma[-1] < math.pi - 1e-6
    assert len(np.unique(np.round(c.sigma, 12))) == c.n_sigma


def test_zero_gain_never_deforms_anything():
    c = TidalCongruence(gain=0.0, n_sigma=181)
    c.advance_to(0.6 * ANTIPODAL_TIME)
    assert np.allclose(c.jacobian(), 1.0, atol=1e-12)
    assert np.allclose(c.expansion(), 0.0, atol=1e-12)


def test_the_antipode_is_never_driven():
    """``h = sin²d·q`` vanishes at the poles, so ``ḧ`` does too — spin weight."""
    c = TidalCongruence(gain=80.0, n_sigma=181)
    k = int(np.argmax(c.distance))            # the station at σ = −π
    assert c.distance[k] == pytest.approx(math.pi, abs=1e-12)
    c.advance_to(1.1 * ANTIPODAL_TIME)
    assert abs(float(c.h_ddot()[k])) < 1e-10
    assert float(c.jacobian()[k]) == pytest.approx(1.0, abs=1e-9)


# ── the launch ──────────────────────────────────────────────────────────────
def test_the_compact_launch_has_finite_support_and_goes_outward():
    sim = Spin2WaveSim(n=800, pulse_width=0.18)
    compact_launch(sim)
    assert np.all(sim.q[sim.d >= 0.18] == 0.0)
    sim.advance_to(1.0)
    d, _ = sim.peak()
    assert d == pytest.approx(1.0, abs=0.05)      # the front travels at one


def test_a_gaussian_launch_reaches_the_far_side_early_and_a_compact_one_does_not():
    r = measure_the_front_is_causal()
    assert r["the_gaussian_beats_the_bound"] is True
    assert r["the_gaussian_arrival_is_grid_converged"] is True
    assert r["the_compact_launch_respects_it"] is True
    for row in r["rows"]:
        assert row["gaussian_arrival"] < row["compact_arrival"]


# ── the equation ────────────────────────────────────────────────────────────
def test_raychaudhuri_holds_to_round_off_with_no_ricci_term():
    """The residual is an identity, so round-off is the *expected* answer.

    Substituting ``θ = A + B`` and ``σ = (A − B)/2`` into ``−θ²/2 − 2σ²`` gives
    ``−A² − B²`` exactly, which is what the deviation equation produces.  The
    test that carries content is the Ricci term being *identically* zero — not
    small, zero — which is what makes the focusing purely shear-driven.
    """
    r = measure_raychaudhuri_is_exact(frames=200)
    assert r["raychaudhuri_is_exact"] is True
    assert r["worst_raychaudhuri_residual"] < 1e-9
    assert r["worst_ricci_term"] == 0.0
    assert r["so_the_focusing_is_all_shear"] is True


def test_the_operator_form_disagrees_with_a_seeded_time_difference():
    """Both converge; they converge to different answers.

    The seeded history injects a spurious ``½ḣ(0)`` impulse whose ``1/dt²``
    cancels against the ``dt`` of the update, so a refinement ladder reproduces
    the wrong number instead of exposing it.
    """
    r = measure_the_operator_form_matters(n_radials=(600, 1200))
    assert r["both_forms_converge"] is True
    assert r["but_they_disagree"] is True
    ops = [row["operator_min_J"] for row in r["rows"]]
    fds = [row["seeded_difference_min_J"] for row in r["rows"]]
    assert min(ops) < 0.0 < min(fds)          # one caustics, the other does not


# ── the three cases ─────────────────────────────────────────────────────────
def test_a_weak_wave_only_dips_and_a_strong_one_caustics():
    r = measure_the_three_cases(gains=(1.0, 150.0), n_sigma=181)
    assert r["both_cases_appear"] is True
    assert r["a_weak_wave_only_dips"] is True
    assert r["rows"][0]["case"] == "ordinary focus"
    assert r["rows"][-1]["case"] == "caustic"


def test_case_three_is_declared_out_of_scope_rather_than_absent():
    r = measure_case_three_is_unreachable(n_sigma=121)
    assert r["case_three_is_out_of_scope"] is True
    assert r["the_geometry_never_moves"] is True
    assert r["evolution_terminated"] is False
    assert r["background_curvature_range"] == 0.0
    assert r["the_field_stayed_finite"] is True
    assert "no Einstein equation" in r["reason"]


# ── the two rings ───────────────────────────────────────────────────────────
def test_the_source_ring_closes_long_before_the_converging_one():
    r = measure_the_caustic_thresholds(n_sigma=181)
    a = r["source_ring_threshold_strain"]
    b = r["converging_ring_threshold_strain"]
    assert a is not None and b is not None
    assert b > 2.0 * a
    assert r["closing_the_ring_needs_an_enormous_strain"] is True
    assert r["window_in_units_of_pi"] == pytest.approx(1.2)


def test_a_threshold_without_a_window_would_be_meaningless():
    """The wave is periodic, so waiting longer lets a weaker wave suffice."""
    r = measure_the_threshold_depends_on_the_window(
        windows=(1.2, 2.0, 3.0), n_sigma=121, n_radial=700)
    assert r["threshold_falls_with_the_window"] is True
    assert r["ratio_first_to_last"] > 1.5
    for row in r["rows"]:
        assert row["threshold_gain"] is not None


# ── the neck ────────────────────────────────────────────────────────────────
def test_the_neck_radius_scales_with_the_pulse_width():
    r = measure_the_neck_is_a_ring(widths=(0.12, 0.24, 0.40), n_sigma=721,
                                   grids=(361, 721))
    assert r["the_radius_scales_with_the_width"] is True
    assert r["the_neck_is_resolved_off_the_pole"] is True
    assert 0.2 < r["mean_radius_over_width"] < 0.8
    # the spread is quantisation, not disagreement
    assert r["relative_spread"] <= 1.05 * r["one_cell_in_the_same_units"]


def test_the_neck_is_searched_in_the_cap_not_at_the_region_cut():
    """A ``d > π/2`` search finds the truncation, not a neck.

    Before the refocus the far half's minimum sits exactly on the cut, which is
    an artefact of where the region was severed; the cap avoids it.
    """
    c = TidalCongruence(gain=60.0, n_sigma=361)
    c.advance_to(0.9 * ANTIPODAL_TIME)
    far = c.far_mask()
    J = c.jacobian()
    at_cut = float(c.distance[far][int(np.argmin(J[far]))])
    assert at_cut == pytest.approx(0.5 * math.pi, abs=0.02)   # the artefact
    cap = (math.pi - c.distance) < NECK_CAP
    assert not np.any(cap & far & (np.abs(c.distance - 0.5 * math.pi) < 0.02))


# ── the passage ─────────────────────────────────────────────────────────────
def test_the_caustic_is_a_passage_and_the_crossing_is_transversal():
    r = measure_the_caustic_is_a_passage(n_sigma=361)
    assert r["the_caustic_is_a_passage"] is True
    assert r["crossings_are_transversal"] is True
    assert r["the_source_excursion_is_resolved"] is True
    assert r["everything_stayed_finite"] is True
    assert abs(r["crossing_slope"]) > 1.0
    assert r["crossing_slope_drift"] < 0.05
    assert abs(r["solver_invariant_drift"]) < 1e-10


def test_the_converging_ring_only_grazes_zero():
    """Reported because it is the answer, not hidden because it is untidy."""
    r = measure_the_caustic_is_a_passage(n_sigma=361)
    assert r["the_converging_ring_only_grazes"] is True
    assert r["converging_ring_depth_drift"] > r[
        "source_depth_drift_under_refinement"]


def test_reconnection_is_recorded_as_structurally_unavailable():
    r = measure_the_caustic_is_a_passage(n_sigma=181)
    assert "never by its neighbours" in r["reconnection_was_never_available"]
    assert r["not_a_termination"] is True


# ── the history panel's data ────────────────────────────────────────────────
def test_the_spacetime_history_has_the_right_shape_and_a_causal_bound():
    h = spacetime_history(gain=60.0, n_sigma=121, frames=40,
                          n_radial=700, t_end=1.1 * ANTIPODAL_TIME)
    assert h["J"].shape == (121, len(h["t"]))
    assert h["pinched"].shape == h["J"].shape
    assert np.all(np.diff(h["t"]) > 0)
    # the bound is d − w, floored at zero, and the antipode is the latest
    assert h["causal_bound"].min() == 0.0
    assert h["causal_bound"].max() == pytest.approx(math.pi - 0.18, abs=1e-9)


def test_nothing_deforms_before_the_front_arrives():
    """The far side is untouched until the causal bound, to round-off."""
    c = TidalCongruence(gain=100.0, n_sigma=361)
    probe = int(np.argmin(np.abs(c.distance - 2.6)))
    bound = float(c.causal_bound(c.distance[probe]))
    c.advance_to(0.75 * bound)
    assert float(c.jacobian()[probe]) == pytest.approx(1.0, abs=1e-9)
