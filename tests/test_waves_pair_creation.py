"""
Tests for the two-wave collision kinematics behind pair creation.

The load-bearing claims are three, and none of them is that the caustic is
bright:

* the threshold is a condition on the **opening angle**, so focusing is neither
  sufficient (collinear beams have ``s = 0`` at any amplitude) nor necessary
  (crossed beams need no focus at all);
* ``1 − cos θ = (1 − cos δ)/sin²t`` is an identity of geodesic triangles, so it
  is the same on ``S²`` and ``S³`` — checked against embedded tangent vectors
  that never use the law of cosines;
* therefore ``s(t)`` is **U-shaped**, the threshold cuts **two** windows, and
  only the antipodal one is a collision of independently propagated waves.

Also held down, because the figure would otherwise lie: the drawn momentum
arrows are a *projection*, and a projection does not preserve angles — by tens
of degrees, and differently at the two crossing points, whose true opening angle
is identical.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.pair_creation import (
    TWO_M,
    WavePair,
    breit_wheeler_cross_section,
    crossing_locus,
    crossing_window,
    invariant_along_the_crossing,
    mandelstam_s,
    measure_focusing_alone_creates_no_invariant_mass,
    measure_only_the_far_window_is_independent,
    measure_the_collision_is_head_on_twice,
    measure_the_cross_section_is_imported_and_checked,
    measure_the_invariant_is_a_triangle_identity,
    measure_the_pair_conserves_orientation,
    measure_the_projected_angle_is_not_the_opening_angle,
    measure_the_threshold_opens_two_windows,
    opening_angle,
    outgoing_momentum,
    threshold_windows,
)


# ── the invariant, which is the whole threshold ─────────────────────────────
@pytest.mark.parametrize("amp", [1.0, 1e3, 1e6, 1e12])
def test_collinear_momenta_have_no_invariant_mass_at_any_amplitude(amp):
    """Focusing raises E. It never reaches threshold, at any brightness."""
    assert mandelstam_s(0.0, amp, amp) == 0.0


def test_head_on_is_the_maximum_and_sets_the_energy_threshold():
    assert mandelstam_s(math.pi, 1.0, 1.0) == pytest.approx(4.0)
    for theta in (0.3, 1.0, 2.0, 3.0):
        assert mandelstam_s(theta) < mandelstam_s(math.pi)


def test_crossed_beams_need_no_focusing():
    """Not necessary either — two beams, no caustic anywhere."""
    assert mandelstam_s(math.pi, 1.5, 1.5) > (TWO_M * 1.0) ** 2
    assert mandelstam_s(math.pi, 0.5, 0.5) < (TWO_M * 1.0) ** 2


# ── the triangle identity ───────────────────────────────────────────────────
def test_crossing_window_endpoints():
    lo, hi = crossing_window(0.42)
    assert lo == pytest.approx(0.21)
    assert hi == pytest.approx(math.pi - 0.21)


@pytest.mark.parametrize("delta", [0.05, 0.42, 1.0, 3.0])
def test_the_window_ends_are_exactly_head_on(delta):
    """Tested on the invariant: arccos near −1 loses eight digits."""
    lo, hi = crossing_window(delta)
    assert invariant_along_the_crossing(delta, lo) == pytest.approx(4.0,
                                                                    abs=1e-12)
    assert invariant_along_the_crossing(delta, hi) == pytest.approx(4.0,
                                                                    abs=1e-12)


@pytest.mark.parametrize("delta", [0.15, 0.42, 1.0, 2.0])
def test_the_equator_angle_is_exactly_the_separation(delta):
    assert opening_angle(delta, math.pi / 2) == pytest.approx(delta, abs=1e-12)


@pytest.mark.parametrize("delta", [0.2, 0.42, 1.1])
def test_the_invariant_is_u_shaped(delta):
    lo, hi = crossing_window(delta)
    ts = np.linspace(lo, hi, 401)
    s = np.array([invariant_along_the_crossing(delta, float(t)) for t in ts])
    i = int(np.argmin(s))
    assert ts[i] == pytest.approx(math.pi / 2, abs=1e-2)
    assert s[0] == pytest.approx(s[-1], abs=1e-12)     # symmetric
    assert s[0] > s[i] and s[-1] > s[i]


def test_the_two_closed_forms_agree():
    for delta in (0.3, 1.2, 2.4):
        lo, hi = crossing_window(delta)
        for t in np.linspace(lo + 1e-6, hi - 1e-6, 50):
            a = invariant_along_the_crossing(delta, float(t))
            b = mandelstam_s(opening_angle(delta, float(t)))
            assert a == pytest.approx(b, rel=1e-9, abs=1e-12)


def test_outside_the_window_the_fronts_do_not_meet():
    lo, hi = crossing_window(0.42)
    with pytest.raises(ValueError):
        opening_angle(0.42, lo - 0.05)
    with pytest.raises(ValueError):
        invariant_along_the_crossing(0.42, hi + 0.05)


def test_separation_must_be_a_real_separation():
    for bad in (0.0, -0.2, math.pi, 4.0):
        with pytest.raises(ValueError):
            crossing_window(bad)


# ── the embedded control ────────────────────────────────────────────────────
@pytest.mark.parametrize("dim", [3, 4])
def test_momenta_are_unit_and_tangent_to_the_sphere(dim):
    pair = WavePair(delta=0.42, energy=1.4, dim=dim)
    lo, hi = pair.window
    for t in (lo + 0.05, math.pi / 2, hi - 0.05):
        for x in pair.locus_at(t, samples=2):
            x = x / np.linalg.norm(x)
            for src in (pair.source_a, pair.source_b):
                p = outgoing_momentum(src, x, t)
                assert np.linalg.norm(p) == pytest.approx(1.0, abs=1e-12)
                assert float(np.dot(p, x)) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("dim", [3, 4])
def test_the_embedded_angle_matches_the_closed_form(dim):
    """The control that never uses the law of cosines."""
    pair = WavePair(delta=0.7, energy=1.4, dim=dim)
    lo, hi = pair.window
    for t in np.linspace(lo + 1e-3, hi - 1e-3, 25):
        x = pair.locus_at(float(t), samples=1)[0]
        x = x / np.linalg.norm(x)
        pa = outgoing_momentum(pair.source_a, x, float(t))
        pb = outgoing_momentum(pair.source_b, x, float(t))
        assert 1.0 - float(np.dot(pa, pb)) == pytest.approx(
            1.0 - math.cos(opening_angle(0.7, float(t))), abs=1e-11)


def test_the_crossing_locus_really_is_equidistant():
    a = np.array([0.0, 0.0, 1.0])
    b = np.array([math.sin(0.5), 0.0, math.cos(0.5)])
    for t in (0.4, 1.2, 2.5):
        for x in crossing_locus(a, b, t, samples=2):
            x = x / np.linalg.norm(x)
            da = math.acos(np.clip(float(np.dot(x, a)), -1, 1))
            db = math.acos(np.clip(float(np.dot(x, b / np.linalg.norm(b))),
                                   -1, 1))
            assert da == pytest.approx(t, abs=1e-10)
            assert db == pytest.approx(t, abs=1e-10)


def test_the_locus_is_two_points_on_S2_and_a_circle_on_S3():
    p2 = WavePair(delta=0.42, dim=3)
    p3 = WavePair(delta=0.42, dim=4)
    assert len(p2.locus_at(1.2, samples=2)) == 2
    ring = p3.locus_at(1.2, samples=24)
    assert len(ring) == 24
    # a genuine circle: every point equidistant from both sources
    for x in ring:
        x = x / np.linalg.norm(x)
        assert math.acos(np.clip(float(np.dot(x, p3.source_a)), -1, 1)) == \
            pytest.approx(1.2, abs=1e-10)


# ── the two windows ─────────────────────────────────────────────────────────
def test_below_the_energy_threshold_there_is_no_window():
    assert threshold_windows(0.42, 0.5, 1.0) == []
    assert threshold_windows(0.42, 0.999, 1.0) == []


def test_at_the_energy_threshold_the_window_has_zero_width():
    assert threshold_windows(0.42, 1.0, 1.0) == []
    # and that is exactly the head-on invariant touching (2m)²
    assert invariant_along_the_crossing(0.42, 0.21) == pytest.approx(4.0)


def test_there_are_two_windows_and_they_are_mirror_images():
    w = threshold_windows(0.42, 1.4, 1.0)
    assert len(w) == 2
    (n0, n1), (f0, f1) = w
    assert n0 + f1 == pytest.approx(math.pi, abs=1e-12)
    assert n1 + f0 == pytest.approx(math.pi, abs=1e-12)


def test_the_windows_merge_above_the_merge_energy():
    merge = 1.0 / math.sin(0.21)
    assert len(threshold_windows(0.42, 0.99 * merge, 1.0)) == 2
    assert len(threshold_windows(0.42, 1.01 * merge, 1.0)) == 1


def test_the_far_window_is_reached_only_after_a_quarter_turn():
    _, (f0, _) = threshold_windows(0.42, 1.4, 1.0)
    assert f0 > math.pi / 2


# ── the imported QED ────────────────────────────────────────────────────────
def test_cross_section_vanishes_at_and_below_threshold():
    assert breit_wheeler_cross_section(4.0, 1.0) == 0.0
    assert breit_wheeler_cross_section(3.0, 1.0) == 0.0


def test_cross_section_peak_matches_the_textbook():
    ss = np.linspace(4.0 + 1e-6, 40.0, 200_000)
    v = np.array([breit_wheeler_cross_section(float(x)) for x in ss])
    i = int(np.argmax(v))
    assert math.sqrt(ss[i]) / 2.0 == pytest.approx(1.40, abs=0.01)
    assert math.sqrt(1.0 - 4.0 / ss[i]) == pytest.approx(0.701, abs=0.005)
    assert v[i] / (8.0 / 3.0) == pytest.approx(0.256, abs=0.005)


def test_cross_section_falls_at_large_s():
    peak = max(breit_wheeler_cross_section(float(x))
               for x in np.linspace(4.1, 40.0, 5000))
    assert breit_wheeler_cross_section(1000.0) < 0.1 * peak


# ── the pair ────────────────────────────────────────────────────────────────
def test_the_pair_must_be_one_up_and_one_down():
    assert WavePair(orientations=(+1, -1)).net_orientation() == 0
    with pytest.raises(ValueError):
        WavePair(orientations=(+1, +1))


def test_the_pair_geometry_is_what_was_asked_for():
    pair = WavePair(delta=0.42, energy=1.4)
    assert math.acos(float(np.dot(pair.source_a, pair.source_b))) == \
        pytest.approx(0.42, abs=1e-12)
    anti_a, anti_b = pair.antipodes
    assert math.acos(np.clip(float(np.dot(anti_a, anti_b)), -1, 1)) == \
        pytest.approx(0.42, abs=1e-12)


def test_only_S2_and_S3_are_accepted():
    for bad in (2, 5):
        with pytest.raises(ValueError):
            WavePair(dim=bad)


# ── the measurements ────────────────────────────────────────────────────────
def test_measure_focusing_alone_creates_no_invariant_mass():
    r = measure_focusing_alone_creates_no_invariant_mass()
    assert r["focusing_is_not_sufficient"]
    assert r["focusing_is_not_necessary"]
    assert r["threshold_is_at_energy_equal_mass"]
    # and the complication is stated rather than buried
    assert r["a_converging_shell_does_contain_opposed_rays"]


def test_measure_the_invariant_is_a_triangle_identity():
    r = measure_the_invariant_is_a_triangle_identity(n_random=60)
    assert r["the_closed_form_is_confirmed"]
    assert r["and_it_is_dimension_independent"]
    assert r["worst_over_all_dimensions"] < 1e-11


def test_measure_the_collision_is_head_on_twice():
    r = measure_the_collision_is_head_on_twice()
    assert r["both_ends_are_head_on"]
    assert r["the_equator_angle_is_exactly_the_separation"]
    assert r["the_minimum_is_at_the_equator"]
    assert r["worst_head_on_error"] < 1e-12


def test_measure_the_threshold_opens_two_windows():
    r = measure_the_threshold_opens_two_windows(scan=4000)
    assert r["the_scan_agrees_with_the_closed_form"]
    assert r["below_E_equals_m_there_is_no_window"]
    assert r["there_are_exactly_two_windows_in_between"]
    assert r["and_they_merge_above"]


def test_measure_only_the_far_window_is_independent():
    r = measure_only_the_far_window_is_independent()
    assert r["near_collision_is_within_the_source_region"]
    assert r["far_collision_is_past_a_quarter_turn"]
    assert r["ratio_of_path_lengths"] > 5.0


def test_measure_the_projected_angle_is_not_the_opening_angle():
    """The figure cannot carry the claim, and this is by how much."""
    r = measure_the_projected_angle_is_not_the_opening_angle(n_t=20)
    assert r["momenta_are_perpendicular_to_the_sphere"]
    assert r["the_projection_misreads_the_angle"]
    assert r["and_it_misreads_the_two_crossings_differently"]
    assert r["worst_projected_error_deg"] > 5.0


def test_measure_the_cross_section_is_imported_and_checked():
    r = measure_the_cross_section_is_imported_and_checked()
    assert r["matches_the_textbook_peak"]
    assert r["zero_at_threshold"] and r["zero_below_threshold"]
    assert r["falls_at_large_s"]


def test_measure_the_pair_conserves_orientation():
    r = measure_the_pair_conserves_orientation()
    assert r["the_labels_cancel"]
    assert r["above_threshold_there"]
    assert r["crossing_locus_size_on_S2"] == 2
    assert r["on_S3_the_locus_is_a_circle"] == 16
    # the throat is flagged as this program's reading, with the bill inherited
    note = r["but_the_throat_is_interpretation"]
    assert "reading" in note and "inherited, not paid" in note
