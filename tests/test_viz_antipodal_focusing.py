"""
Smoke + physics tests for the antipodal-focusing wave study.

These verify the headless wave sims reproduce the load-bearing claim — an
open (plane) geometry disperses while a closed (sphere) geometry refocuses
at the antipode near ``t = πR`` (the PR #166 ``ψ = f/sinχ`` result in
analogue) — and that the ``draw_*`` panels render under the ``Agg``
backend.
"""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")  # must precede pyplot import

import numpy as np
import matplotlib.pyplot as plt
import pytest

from geometrodynamics.viz import (
    PlaneWaveSim,
    SphereWaveSim,
    draw_contrast,
    draw_focus_object,
    draw_focusing_factor,
    focusing_factor,
    measure_refocus,
    plot_contrast_panel,
    plot_focusing_panel,
    plot_object_panel,
)
from geometrodynamics.viz.antipodal_focusing import HALF_PERIOD, PERIOD


@pytest.fixture(autouse=True)
def _cleanup_figs():
    yield
    plt.close("all")


# ── the geometric focusing factor ───────────────────────────────────────────
def test_focusing_factor_unity_at_equator_diverges_at_poles():
    assert focusing_factor(math.pi / 2) == pytest.approx(1.0)
    assert focusing_factor(0.02) > 40.0            # 1/sinθ → ∞ at the pole
    assert focusing_factor(math.pi - 0.02) > 40.0


# ── the sims step cleanly and stay bounded ──────────────────────────────────
def test_plane_sim_steps_and_stays_finite():
    sim = PlaneWaveSim()
    sim.run(400)
    assert sim.t > 0.0
    assert np.all(np.isfinite(sim.u))
    assert np.max(np.abs(sim.u)) < 5.0


def test_sphere_sim_steps_and_stays_finite():
    sim = SphereWaveSim()
    sim.run(400)
    assert sim.t > 0.0
    assert np.all(np.isfinite(sim.u))
    assert np.max(np.abs(sim.u)) < 5.0


# ── the central contrast: open disperses, closed refocuses ──────────────────
def test_plane_peak_decays_open_geometry():
    sim = PlaneWaveSim()
    p0, _ = sim.peak()
    sim.advance_to(HALF_PERIOD)
    p1, _ = sim.peak()
    # the ring has run to the rim and drained: no return path
    assert p1 < 0.4 * p0


def test_sphere_refocuses_at_the_antipode_near_half_period():
    r = measure_refocus(n_periods=0.62)
    # refocus lands at the antipode (θ ≈ π) at the great-circle half-period
    assert r.refocus_theta > 0.9 * math.pi
    assert r.refocus_time == pytest.approx(HALF_PERIOD, abs=0.35)
    # and by then the open plane has dispersed
    assert r.plane_peak_ratio < 0.4


def test_sphere_energy_is_approximately_conserved_closed():
    sim = SphereWaveSim()
    sim.run(10)
    e0 = sim.energy()
    sim.run(600)
    e1 = sim.energy()
    # closed surface: no boundary loss (leapfrog + pole limit ⇒ small drift)
    assert abs(e1 - e0) / e0 < 0.1


def test_sphere_wave_leaves_the_pole_then_returns():
    # peak marches away from the launch pole, then the far-pole caustic
    # brings it back — a signature the plane cannot produce
    sim = SphereWaveSim()
    locs = []
    for _ in range(12):
        sim.run(120)
        locs.append(sim.peak()[1])
    assert max(locs) > 0.85 * math.pi      # reaches the antipode
    # after the refocus the peak location turns back toward the launch pole
    i_max = int(np.argmax(locs))
    assert i_max < len(locs) - 1
    assert locs[-1] < max(locs)


# ── rendering under Agg ─────────────────────────────────────────────────────
def test_contrast_panel_draws_four_axes():
    fig = plot_contrast_panel(t=HALF_PERIOD)
    # two field images + two peak strips
    assert len(fig.axes) == 4


def test_contrast_panel_accepts_history_tracks():
    fig = plt.figure()
    plane, sphere = PlaneWaveSim(), SphereWaveSim()
    hist = {"t": [], "plane": [], "sphere": []}
    for _ in range(5):
        plane.run(20)
        sphere.run(20)
        hist["t"].append(sphere.t)
        hist["plane"].append(plane.peak()[0])
        hist["sphere"].append(sphere.peak()[0])
    axes = draw_contrast(fig, plane, sphere, history=hist)
    assert set(axes) == {"ax_plane", "ax_sphere", "ax_plane_strip", "ax_sphere_strip"}
    assert len(fig.axes) == 4


def test_focusing_panel_renders():
    fig = plot_focusing_panel()
    assert len(fig.axes) == 1
    # the plane, sphere and 1/sinθ curves are all present
    assert len(fig.axes[0].get_lines()) >= 3


def test_object_panel_both_states_render():
    fig_n = plot_object_panel(nucleated=True)
    fig_d = plot_object_panel(nucleated=False)
    assert len(fig_n.axes) == 1
    assert len(fig_d.axes) == 1


def test_draw_functions_rebind_same_figure():
    fig = plt.figure()
    draw_focusing_factor(fig)
    n1 = len(fig.axes)
    draw_focus_object(fig, nucleated=True)
    n2 = len(fig.axes)
    assert n1 == 1
    assert n2 == 1  # object panel uses one axis; old axes cleared


def test_period_constants_consistent():
    assert HALF_PERIOD == pytest.approx(math.pi)
    assert PERIOD == pytest.approx(2 * math.pi)
