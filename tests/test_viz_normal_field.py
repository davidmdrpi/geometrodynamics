"""
Tests for drawing the wave as vectors rather than as the graph of their tips.

The whole point is a contrast: the same field, the same instant, two objects.
The graph is embedded — that is four rounds of negative results — and the normal
field crosses itself freely, because neighbouring normals meet at the centre of
curvature.

Held down here:

* the frame is right (outward normals, correct curvature on shapes with known
  answers);
* the crossing threshold is the radius of curvature, and the focusing drives it
  down;
* the reset from ``R_outer`` to ``R_inner`` is a *separate* mechanism;
* and the duplicated ``σ = ±π`` sample is excluded, because including it puts
  two coincident normals in the bundle and scores a crossing at every length.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.viz.normal_field import (
    NormalField,
    measure_normals_intersect_where_a_graph_cannot,
    measure_the_gap_matters_again,
    measure_the_reset_adds_intersections,
    measure_the_threshold_is_the_radius_of_curvature,
    segment_crossings,
)


# ── the crossing counter ────────────────────────────────────────────────────
def test_the_counter_sees_a_crossing_and_ignores_a_miss():
    a = np.array([[-1.0, 0.0], [0.0, -1.0]])
    b = np.array([[1.0, 0.0], [0.0, 1.0]])
    assert segment_crossings(a, b) == 1
    apart = np.array([[-1.0, 0.0], [-1.0, 5.0]])
    bpart = np.array([[1.0, 0.0], [1.0, 5.0]])
    assert segment_crossings(apart, bpart) == 0


# ── the frame ───────────────────────────────────────────────────────────────
def test_the_normal_points_outward_and_a_circle_has_unit_curvature():
    """A round slice is the control: κ = 1/R everywhere, normals radial."""
    nf = NormalField(delta=0.26, gain=0.0, n_sigma=2001, stride=8)
    nf.reset()
    X, N, kappa = nf.frame()
    r = np.linalg.norm(X, axis=1)
    assert np.allclose(r, nf.shells.r_mid, atol=1e-9)
    assert np.allclose(kappa, 1.0 / nf.shells.r_mid, atol=1e-6)
    # the normal is radial, and outward
    assert np.allclose(N, X / r[:, None], atol=1e-9)
    assert np.all(np.sum(N * X, axis=1) > 0.0)


def test_a_round_slice_has_a_curvature_radius_equal_to_its_radius():
    nf = NormalField(delta=0.26, gain=0.0, n_sigma=2001)
    nf.reset()
    assert nf.envelope_distance() == pytest.approx(nf.shells.r_mid, rel=1e-6)


def test_outward_normals_of_a_circle_never_cross():
    """They diverge — so any crossing later is the wave's doing, not the frame's."""
    nf = NormalField(delta=0.26, gain=0.0, n_sigma=2001, stride=8)
    nf.reset()
    for L in (0.1, 0.5, 0.9):
        assert nf.crossings(L) == 0


# ── the duplicated endpoint ─────────────────────────────────────────────────
def test_the_duplicated_endpoint_is_excluded():
    """``σ = −π`` and ``σ = +π`` are one point; sampling both fakes a crossing.

    With both in the bundle the two coincident normals score a crossing at every
    length, including lengths at which nothing can possibly have met.
    """
    nf = NormalField(gain=0.30, n_sigma=4001, stride=8)
    nf.reset()
    nf.advance_to(2.6)
    idx = nf._sample()
    assert idx[-1] < nf.slice.n_sigma - 1
    a, _ = nf.vectors(1e-4)
    seps = np.linalg.norm(a[:, None, :] - a[None, :, :], axis=-1)
    np.fill_diagonal(seps, np.inf)
    assert seps.min() > 1e-6                       # no coincident bases
    assert nf.crossings(1e-4) == 0                 # nothing crosses at 1e-4


# ── the contrast that matters ───────────────────────────────────────────────
def test_the_graph_never_crosses_but_the_normals_do():
    r = measure_normals_intersect_where_a_graph_cannot(t=3.0)
    assert r["the_graph_never_crosses"] is True
    assert r["graph_self_intersections"] == 0
    assert r["the_normals_do"] is True
    assert r["most_normal_crossings"] > 10


# ── the threshold ───────────────────────────────────────────────────────────
def test_the_threshold_is_the_radius_of_curvature():
    r = measure_the_threshold_is_the_radius_of_curvature(times=(1.2, 3.0))
    assert r["the_envelope_bounds_the_drawn_crossing"] is True
    assert r["the_focus_sharpens_the_surface"] is True
    assert r["rho_at_the_focus"] < r["rho_at_the_start"]


def test_the_drawn_crossing_never_precedes_the_envelope():
    """A finite stride can only lag the continuous envelope, never beat it."""
    r = measure_the_threshold_is_the_radius_of_curvature(times=(2.0, 2.6, 3.0))
    for row in r["rows"]:
        assert row["first_drawn_crossing"] >= row["rho_min"]


def test_crossings_are_monotone_in_the_vector_length():
    nf = NormalField(gain=0.30, n_sigma=2001, stride=8)
    nf.reset()
    nf.advance_to(3.0)
    counts = [nf.crossings(L) for L in (0.10, 0.25, 0.45, 0.70)]
    assert all(b >= a for a, b in zip(counts, counts[1:]))
    assert counts[-1] > counts[0]


# ── the reset ───────────────────────────────────────────────────────────────
def test_the_reset_adds_crossings_and_grows_with_the_stub():
    r = measure_the_reset_adds_intersections(t=3.0)
    assert r["the_reset_adds_crossings"] is True
    assert r["it_grows_with_the_stub"] is True
    for row in r["rows"]:
        assert row["with_reset"] >= row["without_reset"]


def test_the_stub_re_enters_at_the_inner_shell_at_the_exit_angle():
    nf = NormalField(gain=0.30, n_sigma=2001, stride=16)
    nf.reset()
    nf.advance_to(3.0)
    _, clipped, base, stub = nf.vectors_with_reset(0.45)
    assert len(base) > 0
    # every re-entry sits exactly on the inner shell
    assert np.allclose(np.linalg.norm(base, axis=1), nf.shells.r_inner,
                       atol=1e-9)
    # and nothing drawn on the outward leg passes the outer shell
    assert np.max(np.linalg.norm(clipped, axis=1)) <= nf.shells.r_outer + 1e-9
    assert len(stub) == len(base)


def test_nothing_wraps_when_the_vectors_are_short():
    nf = NormalField(delta=0.26, gain=0.30, n_sigma=2001, stride=16)
    nf.reset()
    nf.advance_to(1.0)
    assert nf.n_wrapped(0.01) == 0
    assert nf.crossings(0.01, with_reset=True) == nf.crossings(0.01)


# ── the gap ─────────────────────────────────────────────────────────────────
def test_the_gap_changes_the_count_now():
    """What ``slice_folding`` had severed: the gap is a knob on intersections again."""
    r = measure_the_gap_matters_again(deltas=(0.40, 0.16, 0.09), t=3.0)
    assert r["the_gap_changes_the_count"] is True
    assert r["every_gap_crosses"] is True
    counts = {row["delta"]: row["with_reset"] for row in r["rows"]}
    assert counts[0.09] > 0
