"""
Geometry, physics and rendering tests for the throat-wavefront study.

These pin the load-bearing claims of PR #242 — the catenoidal neck is fixed
by a ``C¹`` join and closes Gauss–Bonnet at ``χ = 0``; the coupled
sphere+neck solve conserves energy; a point pulse's front is a single ring
through free flight; the open-versus-sealed echo delay is the neck
arclength; and the throat's orientation is invisible exactly at the twist
offsets that keep a point source's meridian mirror — and check that every
``draw_*`` panel renders under ``Agg``.
"""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")  # must precede pyplot import

import numpy as np
import matplotlib.pyplot as plt
import pytest

from geometrodynamics.viz import (
    ThroatGeometry,
    ThroatWaveSim,
    draw_geometry,
    draw_hemispheres,
    draw_map,
    draw_neck_strip,
    export_frames,
    measure_echo_delay,
    measure_transmission,
    mirror_symmetry_broken,
    orientation_parity,
    peak_in_window,
    plot_wavefront_panel,
    sphere_image,
    surface_name,
    track_wavefront,
    watch_point,
    wavefront_components,
)
from geometrodynamics.viz.throat_wavefront import ANTIPODAL_TIME


@pytest.fixture(autouse=True)
def _cleanup_figs():
    yield
    plt.close("all")


# ── geometry: the C¹ join fixes the neck ────────────────────────────────────
@pytest.mark.parametrize("a", [0.2, 0.5, 0.9, 1.3])
def test_neck_joins_the_sphere_to_first_order(a):
    g = ThroatGeometry(a)
    assert g.radius(0.0) == pytest.approx(math.sin(a), abs=1e-13)
    assert g.radius_slope(0.0) == pytest.approx(-math.cos(a), abs=1e-13)
    # and symmetrically at the far mouth
    assert g.radius(g.length) == pytest.approx(math.sin(a), abs=1e-12)
    assert g.radius_slope(g.length) == pytest.approx(math.cos(a), abs=1e-12)


@pytest.mark.parametrize("a", [0.2, 0.5, 0.9, 1.3])
def test_gauss_bonnet_closes_at_chi_zero(a):
    """Sphere ``+4π cos a`` and neck ``−4π cos a`` cancel: a handle, χ = 0."""
    g = ThroatGeometry(a)
    assert g.sphere_area == pytest.approx(4.0 * math.pi * math.cos(a))
    assert g.gauss_curvature * g.area == pytest.approx(-4.0 * math.pi * math.cos(a))
    assert g.euler_characteristic() == pytest.approx(0.0, abs=1e-10)


def test_neck_area_matches_the_numerical_surface_of_revolution():
    g = ThroatGeometry(0.5)
    s = np.linspace(0.0, g.length, 20001)
    area = np.trapezoid(2.0 * math.pi * np.asarray(g.radius(s)), s)
    assert area == pytest.approx(g.area, rel=1e-8)


@pytest.mark.parametrize("a", [0.2, 0.5, 0.9, 1.3])
def test_neck_is_embeddable_as_a_surface_of_revolution(a):
    """``|r'| ≤ cos a ≤ 1`` everywhere, so ``z' = √(1−r'²)`` stays real."""
    g = ThroatGeometry(a)
    s = np.linspace(0.0, g.length, 501)
    assert np.max(np.abs(g.radius_slope(s))) <= 1.0 + 1e-12
    z = g.height(s)
    assert np.all(np.isfinite(z))
    assert np.all(np.diff(z) >= -1e-12)          # monotone along the neck


def test_bulk_route_is_a_shortcut_and_the_loop_exceeds_the_echo_by_L():
    g = ThroatGeometry(0.5)
    assert g.inner_route < g.outer_route
    assert g.shortcut_ratio > 3.0
    th0 = 0.5 * math.pi
    assert g.throat_loop(th0) - g.mirror_echo(th0) == pytest.approx(g.length)


def test_mouth_angle_is_validated():
    with pytest.raises(ValueError):
        ThroatGeometry(0.0)
    with pytest.raises(ValueError):
        ThroatGeometry(0.5 * math.pi)


# ── orientation bookkeeping ─────────────────────────────────────────────────
def test_orientation_parity_names_the_surface():
    assert orientation_parity(1) == 1
    assert orientation_parity(-1) == -1
    assert surface_name(1) == "torus"
    assert surface_name(-1) == "klein_bottle"
    with pytest.raises(ValueError):
        orientation_parity(0)


def test_mirror_survives_exactly_at_tau_zero_and_pi():
    n = 128
    assert not mirror_symmetry_broken(0, n)          # τ = 0
    assert not mirror_symmetry_broken(n // 2, n)     # τ = π
    assert mirror_symmetry_broken(n // 4, n)
    assert mirror_symmetry_broken(3 * n // 8, n)


# ── the solve runs, stays finite, conserves energy ──────────────────────────
@pytest.mark.parametrize("mode", ["plugged", "throat"])
def test_sim_steps_and_conserves_energy(mode):
    sim = ThroatWaveSim(mode=mode, n_theta=64, n_phi=96, pulse_width=0.2)
    sim.advance_to(2.5)
    assert sim.t >= 2.5
    assert sim.is_finite()
    assert sim.energy_drift() < 0.03


def test_sim_rejects_bad_configuration():
    with pytest.raises(ValueError):
        ThroatWaveSim(mode="open")
    with pytest.raises(ValueError):
        ThroatWaveSim(n_phi=97)


def test_plugged_run_keeps_no_energy_in_the_neck():
    sim = ThroatWaveSim(mode="plugged", n_theta=64, n_phi=96)
    sim.advance_to(2.0)
    assert sim.neck_energy_fraction() == 0.0


def test_open_throat_fills_the_neck_and_the_plugged_control_does_not():
    kw = dict(n_theta=64, n_phi=96, pulse_width=0.2)
    throat = ThroatWaveSim(mode="throat", **kw)
    t = measure_transmission(throat, 2.5, n_samples=120)
    assert t["transmitted"] > 0.1
    assert t["energy_drift"] < 0.03
    assert t["transmitted"] + t["reflected"] == pytest.approx(1.0)


def test_transmission_grows_with_the_mouth_aperture():
    fracs = []
    for a in (0.3, 0.6):
        sim = ThroatWaveSim(mode="throat", mouth_angle=a, n_theta=64,
                            n_phi=96, pulse_width=0.2)
        fracs.append(measure_transmission(sim, 2.2, n_samples=100)["transmitted"])
    assert fracs[0] < fracs[1]


# ── the wavefront is a single ring while it is in free flight ───────────────
@pytest.mark.parametrize("mode", ["plugged", "throat"])
def test_front_is_one_ring_until_it_reaches_a_mouth(mode):
    sim = ThroatWaveSim(mode=mode, n_theta=96, n_phi=128, pulse_width=0.18)
    free = sim.geom.free_flight(sim.source_theta)
    r = track_wavefront(sim, 2.4, n_samples=96)
    assert r["single_ring_until"] >= free


def test_wavefront_components_counts_one_early():
    sim = ThroatWaveSim(mode="throat", n_theta=96, n_phi=128, pulse_width=0.18)
    sim.advance_to(0.8)                       # well inside free flight
    assert wavefront_components(sim) == 1


# ── the echo delay reads the neck length ────────────────────────────────────
@pytest.mark.slow
def test_echo_delay_reproduces_the_neck_length():
    d = measure_echo_delay(n_theta=96, n_phi=128, pulse_width=0.18)
    assert d["delay_rel_error"] < 0.05
    # and the open throat lets the pulse through instead of bouncing it
    assert d["mirror_suppression"] > 0.7


# ── watched-point helpers ───────────────────────────────────────────────────
def test_watch_point_and_peak_window_agree_on_a_synthetic_peak():
    t = np.linspace(0.0, 4.0, 801)
    a = np.exp(-((t - 2.3) / 0.15) ** 2)
    p = peak_in_window(t, a, 1.5, 3.0)
    assert p.time == pytest.approx(2.3, abs=2e-3)
    assert p.amplitude == pytest.approx(1.0, rel=1e-2)
    assert math.isnan(peak_in_window(t, a, 9.0, 9.5).time)


def test_watch_point_returns_a_series_at_the_source():
    sim = ThroatWaveSim(mode="plugged", n_theta=48, n_phi=64)
    ts, a = watch_point(sim, 0.5 * math.pi, 0.0, 1.0, n_samples=40)
    assert ts.shape == a.shape == (40,)
    assert np.all(np.diff(ts) > 0)
    assert abs(a[0]) > abs(a[-1])              # the cap leaves the source


# ── the antipodal focus still lands where the bare sphere puts it ───────────
@pytest.mark.slow
def test_antipodal_focus_is_near_the_bare_sphere_half_period():
    """The caps and the neck perturb, they do not relocate, the focus."""
    sim = ThroatWaveSim(mode="plugged", n_theta=96, n_phi=128, pulse_width=0.22)
    ts, a = watch_point(sim, 0.5 * math.pi, math.pi, 3.6, n_samples=600)
    p = peak_in_window(ts, a, 2.4, 3.6)
    assert abs(p.time - ANTIPODAL_TIME) < 0.35


# ── orientation is hidden at τ ∈ {0, π} ─────────────────────────────────────
@pytest.mark.parametrize("steps,visible", [(0, False), (32, True), (64, False)])
def test_orientation_visibility_matches_the_mirror_argument(steps, visible):
    fields = []
    for parity in (1, -1):
        s = ThroatWaveSim(mode="throat", n_theta=64, n_phi=128,
                          twist_parity=parity, twist_steps=steps,
                          pulse_width=0.2)
        s.advance_to(2.6)
        fields.append(s.u.copy())
    rel = np.max(np.abs(fields[0] - fields[1])) / np.max(np.abs(fields[0]))
    assert mirror_symmetry_broken(steps, 128) is visible
    if visible:
        assert rel > 1e-2
    else:
        assert rel < 1e-9


# ── rendering ───────────────────────────────────────────────────────────────
def _small(mode="throat"):
    s = ThroatWaveSim(mode=mode, n_theta=48, n_phi=64, pulse_width=0.22)
    s.advance_to(1.3)
    return s


def test_sphere_image_masks_the_caps_and_the_disc_exterior():
    sim = _small()
    img, mask = sphere_image(sim, n_px=64)
    assert img.shape == (64, 64)
    assert np.isnan(img[0, 0])                 # outside the disc
    assert np.any(mask)
    assert np.all(np.isfinite(img[mask]))


def test_draw_functions_render():
    sim = _small()
    fig, axes = plt.subplots(2, 2, figsize=(6, 5))
    draw_hemispheres(axes[0][0], axes[0][1], sim, n_px=64)
    draw_neck_strip(axes[1][0], sim)
    draw_map(axes[1][1], sim)
    assert fig.axes


def test_draw_neck_strip_annotates_the_sealed_control():
    sim = _small(mode="plugged")
    fig, ax = plt.subplots()
    draw_neck_strip(ax, sim)
    assert any("sealed" in t.get_text() for t in ax.texts)


def test_draw_geometry_renders_both_pieces():
    fig, ax = plt.subplots()
    draw_geometry(ax, ThroatGeometry(0.5))
    assert len(ax.lines) >= 4                  # two sphere arcs, two neck arcs


def test_plot_wavefront_panel_builds_a_full_figure():
    fig = plot_wavefront_panel(_small(), figsize=(8.0, 5.0))
    assert len(fig.axes) == 5


# ── export ──────────────────────────────────────────────────────────────────
def test_export_frames_round_trips_to_the_declared_shape():
    import base64

    sim = ThroatWaveSim(mode="throat", n_theta=48, n_phi=64, pulse_width=0.22)
    data = export_frames(sim, t_end=1.2, frames=6, rows=16, cols=24,
                         neck_rows=6)
    raw = base64.b64decode(data["sphere_b64"])
    assert len(raw) == 6 * 16 * 24
    neck = base64.b64decode(data["neck_b64"])
    assert len(neck) == 6 * 6 * 24
    assert data["scale"] > 0.0
    assert len(data["times"]) == 6
    assert data["surface"] == "torus"
    arr = np.frombuffer(raw, dtype=np.uint8)
    assert arr.min() >= 0 and arr.max() <= 255
    assert data["compand"] == "linear"


def test_signed_sqrt_companding_preserves_sign_and_lifts_weak_structure():
    import base64

    sim = ThroatWaveSim(mode="throat", n_theta=48, n_phi=64, pulse_width=0.22)
    kw = dict(t_end=2.4, frames=10, rows=16, cols=24, neck_rows=6)
    lin = export_frames(sim, compand="linear", **kw)
    sqr = export_frames(sim, compand="signed_sqrt", **kw)
    a = np.frombuffer(base64.b64decode(lin["sphere_b64"]), dtype=np.uint8).astype(int) - 128
    b = np.frombuffer(base64.b64decode(sqr["sphere_b64"]), dtype=np.uint8).astype(int) - 128
    assert sqr["compand"] == "signed_sqrt"
    # the sign is untouched wherever the linear encoding resolved one at all
    live = (a != 0) & (np.abs(a) < 127)
    assert np.all(np.sign(a[live]) == np.sign(b[live]))
    # and the weak tail is lifted out of the floor
    assert np.mean(np.abs(b) > 8) > np.mean(np.abs(a) > 8)


def test_export_rejects_an_unknown_companding():
    sim = ThroatWaveSim(mode="throat", n_theta=32, n_phi=48)
    with pytest.raises(ValueError):
        export_frames(sim, t_end=0.3, frames=2, compand="log")
