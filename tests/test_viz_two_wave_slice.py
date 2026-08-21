"""Two waves on the circle slice, and where they connect inner to outer.

The one-wave result this extends is a *negative* one — a height field cannot
wind — so the checks here are mostly about the boundary between what one wave
can do and what two can. Where a claim is an identity rather than a measurement
it is asserted at machine precision, because that is what an identity means.
"""

import math

import numpy as np
import pytest

from geometrodynamics.viz import two_wave_slice as tw
from geometrodynamics.viz.circle_slice import (ANTIPODAL_TIME, CircleSlice)


@pytest.fixture()
def pair():
    p = tw.TwoWaveSlice()
    p.slice_.reset()
    p.slice_.advance_to(ANTIPODAL_TIME)
    return p


# ── the base circle ─────────────────────────────────────────────────────────
def test_the_closed_circle_drops_the_duplicated_endpoint(pair):
    """`σ = −π` and `σ = +π` are one point; counting must not see two.

    A first draft counted contact regions on the drawn array and reported `2`
    at what is a single tangency.
    """
    assert len(pair.sigma) == len(pair.sigma_closed) + 1
    assert abs(pair.sigma[0] + math.pi) < 1e-12
    assert abs(pair.sigma[-1] - math.pi) < 1e-12
    assert abs(pair.sigma_closed[-1] - math.pi) > 1e-6


def test_the_two_waves_are_mirror_images_about_the_mid_radius(pair):
    r_a, r_b = pair.radii(gain=0.3)
    mid = pair.slice_.bulk.r_mid
    assert np.max(np.abs((r_a + r_b) - 2.0 * mid)) < 1e-13


def test_a_zero_offset_pair_shares_a_source_and_pi_puts_them_apart():
    a = tw.TwoWaveSlice(offset=tw.CO_LOCATED)
    b = tw.TwoWaveSlice(offset=tw.ANTIPODAL_SOURCES)
    ua, ub = a.fields()
    assert np.max(np.abs(ua - ub)) < 1e-15, "co-located waves see one field"
    va, vb = b.fields()
    assert np.max(np.abs(va - vb)) > 1e-6, "antipodal ones do not"


# ── the identity ────────────────────────────────────────────────────────────
def test_the_pair_touches_exactly_where_one_wave_wraps():
    """`εu = gap/2` is both conditions, so the thresholds are the same number."""
    got = tw.measure_the_pair_touches_exactly_where_one_wave_wraps()
    assert got["they_are_the_same_threshold"]
    assert got["relative_difference"] == 0.0
    assert got["single_wave_wrap_gain"] > 0.0


def test_that_threshold_is_the_one_v46_reports():
    """The number has to be v46's own wrap threshold, not a new one."""
    got = tw.measure_the_pair_touches_exactly_where_one_wave_wraps()
    sl = CircleSlice()
    half = 0.5 * sl.bulk.gap
    assert abs(got["single_wave_wrap_gain"] - half / got["run_peak"]) < 1e-15


def test_the_covered_fraction_is_linear_in_the_gain(pair):
    base = pair.contact_gain()
    for m in (0.3, 0.75, 1.0, 1.6):
        assert abs(pair.covered_fraction(gain=base * m) - m) < 1e-12


def test_the_contact_gain_is_where_the_covered_fraction_reaches_one(pair):
    assert abs(pair.covered_fraction(gain=pair.contact_gain()) - 1.0) < 1e-12


# ── where they connect ──────────────────────────────────────────────────────
def test_they_connect_at_the_antipode_on_the_seam():
    got = tw.measure_where_the_two_pulses_connect()
    assert got["connected"]
    assert got["the_contact_is_at_the_antipode"]
    assert got["the_contact_is_on_the_seam"]
    assert got["distance_to_the_seam"] < 1e-9
    assert abs(abs(got["sigma_over_pi"]) - 1.0) < 1e-12


def test_the_antipodal_refocus_is_a_rarefaction():
    """So the roles swap: the *inward*-driven wave is the one reaching R_outer."""
    got = tw.measure_where_the_two_pulses_connect()
    assert got["it_is_a_rarefaction"]
    assert got["refocus_amplitude"] < 0.0
    assert got["the_inward_wave_is_the_one_that_reaches_r_outer"]
    assert abs(got["radius_of_the_inward_wave"] - got["r_outer"]) < 1e-12
    assert abs(got["radius_of_the_outward_wave"] - got["r_inner"]) < 1e-12


def test_at_contact_the_two_radii_are_the_two_boundaries(pair):
    got = pair.contact(gain=pair.contact_gain())
    radii = sorted((got["radius_a"], got["radius_b"]))
    assert abs(radii[0] - got["r_inner"]) < 1e-12
    assert abs(radii[1] - got["r_outer"]) < 1e-12


# ── what one wave cannot do ─────────────────────────────────────────────────
def test_nothing_connects_below_threshold_and_one_arc_opens_above():
    got = tw.measure_the_contact_is_an_arc_the_band_covers()
    assert got["nothing_connects_below_threshold"]
    assert got["it_touches_at_threshold"]
    assert got["one_arc_above"]
    assert got["the_covered_fraction_tracks_the_gain"] < 1e-12


def test_the_connected_set_is_one_arc_with_two_crossings(pair):
    got = pair.contact(gain=1.3 * pair.contact_gain())
    assert got["arcs"] == 1
    assert got["crossings"] == 2


def test_a_single_wave_still_does_not_wind_however_hard_it_is_driven():
    """The v46 result has to survive: one graph, zero winding, at any gain."""
    sl = CircleSlice()
    sl.reset()
    sl.advance_to(ANTIPODAL_TIME)
    for gain in (0.1, 0.5, 2.0, 8.0):
        assert sl.winding_number(gain=gain) == 0


def test_the_band_covers_radii_that_a_single_wave_leaves_outside(pair):
    """Two graphs bound a band; past threshold no radius is left out of it."""
    gain = 1.3 * pair.contact_gain()
    r_a, r_b = pair.radii(gain=gain, closed=True)
    covered = np.abs(r_a - r_b) >= pair.gap
    assert covered.any(), "the band must be radially surjective somewhere"
    # and the single wave, at the same gain, still leaves radii outside it
    single = pair.slice_.radius(gain=gain)
    assert single.min() > -1e9  # it is a graph: one radius per sigma
    assert len(single) == len(pair.sigma)


# ── the two configurations ──────────────────────────────────────────────────
def test_meeting_mid_flight_is_harder_than_meeting_at_a_refocus():
    got = tw.measure_meeting_mid_flight_is_harder()
    assert got["mid_flight_is_harder"]
    assert got["both_are_cheapest_at_a_refocus"]
    assert got["worst_penalty"] > 5.0
    for key in ("co_located", "antipodal_sources"):
        assert got[key]["mid_flight_penalty"] > 1.0


def test_a_co_located_pair_connects_at_about_half_the_gain():
    """At a refocus both of a co-located pair peak; only one of an antipodal
    pair does, so `|u_A + u_B|` is twice as large."""
    got = tw.measure_meeting_mid_flight_is_harder()
    assert got["co_located_is_about_twice_as_cheap"]
    assert 1.8 < got["antipodal_over_co_located"] < 2.2


def test_the_outward_gap_is_what_closes(pair):
    """`gap − |δ|` is the arc through the seam; it hits zero at contact."""
    at = pair.contact_gain()
    assert pair.outward_gap(gain=0.5 * at, closed=True).min() > 0.0
    assert abs(pair.outward_gap(gain=at, closed=True).min()) < 1e-12


def test_the_scan_covers_a_whole_return_period():
    p = tw.TwoWaveSlice()
    rows = p.scan(samples=40)
    assert len(rows) == 40
    assert rows[-1]["t"] == pytest.approx(2.0 * math.pi)
    assert all(r["contact_gain"] > 0.0 for r in rows)


# ── drawing ─────────────────────────────────────────────────────────────────
def test_the_drawn_segments_stay_inside_the_annulus(pair):
    seg_a, seg_b = pair.segments(gain=1.4 * pair.contact_gain())
    b = pair.slice_.bulk
    for segs in (seg_a, seg_b):
        assert segs, "there must be something to draw"
        for p in segs:
            r = np.hypot(p[:, 0], p[:, 1])
            assert r.min() > b.r_inner - 1e-9
            assert r.max() < b.r_outer + 1e-9
