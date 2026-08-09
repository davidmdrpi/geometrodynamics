"""
Geometry, physics and rendering tests for the throat-wavefront study.

These pin the load-bearing claims of PR #242:

* the neck is a genuine **catenoid** — ``r = b cosh(z/b)`` to machine
  precision, with ``b = sin²a`` and ``L = sin 2a`` — whose curvature varies
  from ``−1`` at each mouth to ``−1/sin⁴a`` at the waist;
* ``∫K dA`` over a surface of revolution depends only on the boundary
  slopes, so ``χ = 0`` tests the ``C¹`` join and *not* the profile;
* the mouths are a single shared finite-volume face, so the coupled solve
  conserves its discrete energy to round-off;
* a front on a closed surface with no boundary crosses each point once,
  while a sealed mouth sends a second front back toward the source and an
  open one does not;
* the open/sealed echo delay is the neck arclength;
* the throat's orientation is invisible exactly where a point source's own
  meridian mirror survives the gluing;

and that every ``draw_*`` panel renders under ``Agg``.
"""

from __future__ import annotations

import base64
import math

import matplotlib

matplotlib.use("Agg")  # must precede pyplot import

import numpy as np
import matplotlib.pyplot as plt
import pytest

from geometrodynamics.viz import (
    BareSphereSim,
    ThroatGeometry,
    ThroatWaveSim,
    draw_geometry,
    draw_hemispheres,
    draw_map,
    draw_neck_strip,
    export_frames,
    measure_arrival_multiplicity,
    measure_echo_delay,
    measure_mouth_budget,
    mirror_symmetry_broken,
    orientation_parity,
    peak_in_window,
    plot_wavefront_panel,
    sphere_image,
    surface_name,
    watch_point,
)
from geometrodynamics.viz.throat_wavefront import ANTIPODAL_TIME

MOUTH = 0.75


@pytest.fixture(autouse=True)
def _cleanup_figs():
    yield
    plt.close("all")


# ── the neck is a real catenoid ─────────────────────────────────────────────
@pytest.mark.parametrize("a", [0.2, 0.5, 0.75, 1.1, 1.4])
def test_neck_is_exactly_a_catenoid(a):
    """``r = b cosh(z/b)`` in the axial gauge, to machine precision."""
    g = ThroatGeometry(a)
    b = g.neck_radius
    s = np.linspace(0.0, g.length, 4001)
    r = np.asarray(g.radius(s))
    z = np.asarray(g.height(s))
    assert np.max(np.abs(r - b * np.cosh(z / b))) < 1e-12


@pytest.mark.parametrize("a", [0.2, 0.5, 0.75, 1.1, 1.4])
def test_closed_forms_for_waist_and_length(a):
    g = ThroatGeometry(a)
    assert g.neck_radius == pytest.approx(math.sin(a) ** 2)
    assert g.length == pytest.approx(math.sin(2.0 * a))
    assert g.radius(0.0) == pytest.approx(math.sin(a), abs=1e-13)
    assert g.radius_slope(0.0) == pytest.approx(-math.cos(a), abs=1e-13)
    assert g.radius(g.length) == pytest.approx(math.sin(a), abs=1e-12)
    assert g.radius_slope(g.length) == pytest.approx(math.cos(a), abs=1e-12)


@pytest.mark.parametrize("a", [0.2, 0.5, 0.75, 1.1])
def test_curvature_varies_and_is_minus_one_at_the_mouth(a):
    """The sphere's ``+1`` flips to exactly ``−1``, then deepens to the waist."""
    g = ThroatGeometry(a)
    assert g.curvature_at_mouth == pytest.approx(-1.0, abs=1e-12)
    assert g.curvature_at_waist == pytest.approx(-1.0 / math.sin(a) ** 4)
    assert g.curvature_at_waist < g.curvature_at_mouth      # strictly deeper
    s = np.linspace(0.0, g.length, 501)
    k = np.asarray(g.curvature(s))
    assert k.min() == pytest.approx(g.curvature_at_waist, rel=1e-6)
    assert np.ptp(k) > 0.1                                  # genuinely varying


@pytest.mark.parametrize("a", [0.2, 0.5, 0.75, 1.1])
def test_arclength_parametrisation_is_unit_speed(a):
    g = ThroatGeometry(a)
    s = np.linspace(0.0, g.length, 20001)
    dr = np.gradient(np.asarray(g.radius(s)), s)
    dz = np.gradient(np.asarray(g.height(s)), s)
    assert np.max(np.abs(dr[5:-5] ** 2 + dz[5:-5] ** 2 - 1.0)) < 1e-5
    assert np.max(np.abs(g.radius_slope(s))) <= 1.0         # embeddable


@pytest.mark.parametrize("a", [0.2, 0.5, 0.75, 1.1])
def test_neck_area_matches_the_closed_form(a):
    g = ThroatGeometry(a)
    s = np.linspace(0.0, g.length, 40001)
    area = np.trapezoid(2.0 * math.pi * np.asarray(g.radius(s)), s)
    assert area == pytest.approx(g.area, rel=1e-7)


# ── Gauss–Bonnet closes on the join, not on the profile ─────────────────────
@pytest.mark.parametrize("a", [0.2, 0.5, 0.75, 1.1, 1.4])
def test_total_curvature_depends_only_on_the_boundary_slopes(a):
    g = ThroatGeometry(a)
    s = np.linspace(0.0, g.length, 40001)
    r = np.asarray(g.radius(s))
    numeric = np.trapezoid(np.asarray(g.curvature(s)) * 2.0 * math.pi * r, s)
    assert numeric == pytest.approx(g.neck_total_curvature(), rel=1e-6)
    assert g.neck_total_curvature() == pytest.approx(-4.0 * math.pi * math.cos(a))
    assert g.euler_characteristic() == pytest.approx(0.0, abs=1e-12)


def test_chi_zero_is_insensitive_to_the_profile():
    """A different C¹-matched profile closes χ = 0 just as exactly.

    Pins the honest reading of the closure: it constrains the join, so it
    cannot be quoted as evidence for the catenoid in particular.
    """
    a = 0.75
    g = ThroatGeometry(a)
    s = np.linspace(0.0, g.length, 40001)
    # a cubic Hermite profile matching r and r' at both ends
    x = s / g.length
    h = (2 * x ** 3 - 3 * x ** 2 + 1) * math.sin(a) + \
        (x ** 3 - 2 * x ** 2 + x) * (-math.cos(a) * g.length) + \
        (-2 * x ** 3 + 3 * x ** 2) * math.sin(a) + \
        (x ** 3 - x ** 2) * (math.cos(a) * g.length)
    d2 = np.gradient(np.gradient(h, s), s)
    total = np.trapezoid((-d2 / h) * 2.0 * math.pi * h, s)
    assert total == pytest.approx(-4.0 * math.pi * math.cos(a), rel=1e-3)


def test_bulk_route_is_a_shortcut_and_the_loop_exceeds_the_echo_by_L():
    g = ThroatGeometry(MOUTH)
    assert g.inner_route < g.outer_route
    assert g.shortcut_ratio > 1.0
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
    assert not mirror_symmetry_broken(0, n)
    assert not mirror_symmetry_broken(n // 2, n)
    assert mirror_symmetry_broken(n // 4, n)
    assert mirror_symmetry_broken(3 * n // 8, n)


# ── the shared mouth face makes the scheme conservative ─────────────────────
@pytest.mark.parametrize("mode", ["plugged", "throat"])
def test_discrete_energy_is_conserved_to_round_off(mode):
    sim = ThroatWaveSim(mode=mode, mouth_angle=MOUTH, n_theta=64, n_phi=96,
                        pulse_width=0.22)
    sim.advance_to(4.0)
    assert sim.is_finite()
    assert sim.energy_drift() < 1e-10


def test_mouth_power_integrates_to_the_neck_energy():
    """The coupling is conservative, not merely plausible."""
    sim = ThroatWaveSim(mode="throat", mouth_angle=MOUTH, n_theta=64,
                        n_phi=96, pulse_width=0.22)
    acc = 0.0
    for _ in range(900):
        p_n, p_s = sim.mouth_power()
        acc += (p_n + p_s) * sim.dt
        sim.step()
    stored = sim.neck_energy()
    assert stored > 1e-3
    assert acc == pytest.approx(stored, rel=0.05)


def test_sim_rejects_bad_configuration():
    with pytest.raises(ValueError):
        ThroatWaveSim(mode="open")
    with pytest.raises(ValueError):
        ThroatWaveSim(n_phi=97)


def test_plugged_run_keeps_no_energy_in_the_neck():
    sim = ThroatWaveSim(mode="plugged", mouth_angle=MOUTH, n_theta=64, n_phi=96)
    sim.advance_to(2.0)
    assert sim.neck_energy_fraction() == 0.0


# ── the mouth budget is a flux, not a snapshot ──────────────────────────────
def test_mouth_budget_is_a_transmission_coefficient():
    sim = ThroatWaveSim(mode="throat", mouth_angle=MOUTH, n_theta=96,
                        n_phi=128, pulse_width=0.18)
    b = measure_mouth_budget(sim)
    assert 0.0 < b["transmission"] <= 1.0
    assert b["transmission"] + b["reflection"] == pytest.approx(1.0)
    assert b["through"] <= b["offered"]
    assert b["energy_drift"] < 1e-10


def test_transmission_grows_with_the_mouth_aperture():
    vals = []
    for a in (0.4, 0.8):
        sim = ThroatWaveSim(mode="throat", mouth_angle=a, n_theta=96,
                            n_phi=128, pulse_width=0.18)
        vals.append(measure_mouth_budget(sim)["transmission"])
    assert vals[0] < vals[1]


def test_mouth_budget_needs_an_open_throat():
    with pytest.raises(ValueError):
        measure_mouth_budget(ThroatWaveSim(mode="plugged", n_theta=48, n_phi=64))


# ── the bare sphere ─────────────────────────────────────────────────────────
def test_bare_sphere_field_depends_only_on_geodesic_distance():
    sim = BareSphereSim(n_theta=64, n_phi=96, pulse_width=0.2)
    sim.advance_to(1.1)
    # two points equidistant from the source must carry the same value
    assert sim.sample(0.5 * math.pi, 0.8) == pytest.approx(
        sim.sample(0.5 * math.pi, -0.8), abs=1e-12)
    assert sim.u.shape == (64, 96)
    assert sim.energy_drift() < 1e-2
    assert sim.neck_energy_fraction() == 0.0
    assert sim.v.shape[0] == 0


def test_bare_sphere_front_reaches_the_antipode_near_half_period():
    sim = BareSphereSim(n_theta=48, n_phi=64, pulse_width=0.2)
    ts, a = watch_point(sim, 0.5 * math.pi, math.pi, 3.6, n_samples=500)
    p = peak_in_window(ts, a, 2.2, 3.6)
    assert abs(p.time - ANTIPODAL_TIME) < 0.4


# ── arrival multiplicity: the honest self-intersection diagnostic ───────────
@pytest.mark.slow
def test_bare_front_crosses_each_point_once_and_never_returns_to_the_source():
    """No boundary, no second front: the ring cannot meet itself."""
    sim = BareSphereSim(n_theta=64, n_phi=96, pulse_width=0.2, n_radial=500)
    m = measure_arrival_multiplicity(sim, math.pi)
    assert m.area_fraction_multi < 0.10
    assert m.source_side_fraction < 1e-9


@pytest.mark.slow
def test_a_sealed_mouth_sends_a_front_back_and_an_open_one_does_not():
    """At the resolution the claim is made at — the mouth needs resolving.

    On a coarser grid the numerically under-resolved mouth reflects enough
    to put some second fronts back on the source side even when the throat
    is open, and the separation degrades to a factor of two.
    """
    kw = dict(mouth_angle=MOUTH, n_theta=96, n_phi=128, pulse_width=0.18)
    plugged = measure_arrival_multiplicity(
        ThroatWaveSim(mode="plugged", **kw), math.pi)
    throat = measure_arrival_multiplicity(
        ThroatWaveSim(mode="throat", **kw), math.pi)
    # both put a second front somewhere — only the mirror puts one back home
    assert plugged.area_fraction_multi > 0.05
    assert throat.area_fraction_multi > 0.05
    assert plugged.source_side_fraction > 0.01
    assert throat.source_side_fraction < plugged.source_side_fraction / 4.0


def test_multiplicity_counts_have_the_grid_shape():
    sim = ThroatWaveSim(mode="plugged", mouth_angle=MOUTH, n_theta=32,
                        n_phi=48, pulse_width=0.25)
    m = measure_arrival_multiplicity(sim, 1.4)
    assert m.counts.shape == (32, 48)
    assert m.max_arrivals >= 1
    assert 0.0 <= m.area_fraction_multi <= 1.0


# ── the echo delay reads the neck length ────────────────────────────────────
@pytest.mark.slow
def test_echo_delay_reproduces_the_neck_length():
    d = measure_echo_delay(mouth_angle=MOUTH, n_theta=96, n_phi=128,
                           pulse_width=0.18)
    assert d["delay_rel_error"] < 0.03
    assert d["mirror_suppression"] > 0.4


# ── watched-point helpers ───────────────────────────────────────────────────
def test_peak_window_finds_a_synthetic_peak_to_sub_sample_accuracy():
    t = np.linspace(0.0, 4.0, 801)
    a = np.exp(-((t - 2.3) / 0.15) ** 2)
    p = peak_in_window(t, a, 1.5, 3.0)
    assert p.time == pytest.approx(2.3, abs=2e-3)
    assert p.amplitude == pytest.approx(1.0, rel=1e-2)
    assert math.isnan(peak_in_window(t, a, 9.0, 9.5).time)


def test_watch_point_returns_a_series_at_the_source():
    sim = ThroatWaveSim(mode="plugged", mouth_angle=MOUTH, n_theta=48, n_phi=64)
    ts, a = watch_point(sim, 0.5 * math.pi, 0.0, 1.0, n_samples=40)
    assert ts.shape == a.shape == (40,)
    assert np.all(np.diff(ts) > 0)
    assert abs(a[0]) > abs(a[-1])


# ── orientation is hidden at τ ∈ {0, π} ─────────────────────────────────────
@pytest.mark.parametrize("steps,visible", [(0, False), (32, True), (64, False)])
def test_orientation_visibility_matches_the_mirror_argument(steps, visible):
    fields = []
    for parity in (1, -1):
        s = ThroatWaveSim(mode="throat", mouth_angle=MOUTH, n_theta=64,
                          n_phi=128, twist_parity=parity, twist_steps=steps,
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
    s = ThroatWaveSim(mode=mode, mouth_angle=MOUTH, n_theta=48, n_phi=64,
                      pulse_width=0.22)
    s.advance_to(1.1)
    return s


def test_sphere_image_masks_the_caps_and_the_disc_exterior():
    sim = _small()
    img, mask = sphere_image(sim, n_px=64)
    assert img.shape == (64, 64)
    assert np.isnan(img[0, 0])
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
    draw_geometry(ax, ThroatGeometry(MOUTH))
    assert len(ax.lines) >= 4


def test_plot_wavefront_panel_builds_a_full_figure():
    fig = plot_wavefront_panel(_small(), figsize=(8.0, 5.0))
    assert len(fig.axes) == 5


# ── export ──────────────────────────────────────────────────────────────────
def test_export_frames_round_trips_to_the_declared_shape():
    sim = ThroatWaveSim(mode="throat", mouth_angle=MOUTH, n_theta=48,
                        n_phi=64, pulse_width=0.22)
    data = export_frames(sim, t_end=1.2, frames=6, rows=16, cols=24,
                         neck_rows=6)
    raw = base64.b64decode(data["sphere_b64"])
    assert len(raw) == 6 * 16 * 24
    assert len(base64.b64decode(data["neck_b64"])) == 6 * 6 * 24
    assert data["scale"] > 0.0
    assert len(data["times"]) == 6
    assert data["surface"] == "torus"
    assert data["compand"] == "linear"


def test_export_frames_handles_a_sim_without_a_neck():
    sim = BareSphereSim(n_theta=48, n_phi=64, pulse_width=0.22)
    data = export_frames(sim, t_end=1.0, frames=4, rows=16, cols=24,
                         neck_rows=6)
    assert len(base64.b64decode(data["sphere_b64"])) == 4 * 16 * 24
    assert data["mode"] == "bare"
    assert data["neck_length"] == 0.0
    assert all(f == 0.0 for f in data["neck_energy_fraction"])


def test_signed_sqrt_companding_preserves_sign_and_lifts_weak_structure():
    sim = ThroatWaveSim(mode="throat", mouth_angle=MOUTH, n_theta=48,
                        n_phi=64, pulse_width=0.22)
    kw = dict(t_end=2.4, frames=10, rows=16, cols=24, neck_rows=6)
    lin = export_frames(sim, compand="linear", **kw)
    sqr = export_frames(sim, compand="signed_sqrt", **kw)
    a = np.frombuffer(base64.b64decode(lin["sphere_b64"]), dtype=np.uint8).astype(int) - 128
    b = np.frombuffer(base64.b64decode(sqr["sphere_b64"]), dtype=np.uint8).astype(int) - 128
    assert sqr["compand"] == "signed_sqrt"
    live = (a != 0) & (np.abs(a) < 127)
    assert np.all(np.sign(a[live]) == np.sign(b[live]))
    assert np.mean(np.abs(b) > 8) > np.mean(np.abs(a) > 8)


def test_export_rejects_an_unknown_companding():
    sim = ThroatWaveSim(mode="throat", n_theta=32, n_phi=48)
    with pytest.raises(ValueError):
        export_frames(sim, t_end=0.3, frames=2, compand="log")
