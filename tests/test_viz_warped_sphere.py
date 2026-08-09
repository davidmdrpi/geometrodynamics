"""
Tests for the restored geometry: one continuous S², warped by a solved wave.

These pin the properties that make the picture mean something:

* the surface is **one closed manifold** — poles single-valued, seam matched;
* it stays strictly **between the two dolls**, at every time, by construction;
* the displacement is an **order-preserving** map of the solved field, so no
  feature is invented and no sign is flipped;
* the largest deformation happens **at the antipode**, at ``t ≈ π``, with the
  residual error explained by the pulse's own width and shrinking with it;
* and the rings drawn on the surface really are **on** it.
"""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")  # must precede pyplot import

import numpy as np
import pytest

from geometrodynamics.viz.warped_sphere import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    NestedShells,
    WarpedSphere,
    measure_containment,
    measure_focus,
)

# small but honest: the field is 1-D, so a coarse *mesh* costs no accuracy
COARSE = dict(n_theta=41, n_phi=61)


@pytest.fixture(scope="module")
def surface():
    return WarpedSphere(**COARSE)


# ── the dolls ───────────────────────────────────────────────────────────────
def test_shells_default_to_the_radial_caustic_vacuole():
    sh = NestedShells()
    assert sh.r_inner == pytest.approx(0.74)
    assert sh.r_outer == pytest.approx(1.26)
    assert sh.r_mid == pytest.approx(1.0)


def test_shells_validate_their_gap():
    with pytest.raises(ValueError):
        NestedShells(r_mid=1.0, delta=1.0)      # inner doll would collapse
    with pytest.raises(ValueError):
        NestedShells(r_mid=1.0, delta=0.0)
    with pytest.raises(ValueError):
        NestedShells(r_mid=-1.0, delta=0.2)


def test_shell_containment_is_strict():
    sh = NestedShells()
    assert sh.contains([0.75, 1.0, 1.25]) is True
    assert sh.contains([0.74]) is False          # touching is not inside
    assert sh.contains([1.26]) is False


def test_unit_sphere_mesh_is_closed():
    X, Y, Z = NestedShells().unit_sphere(21, 33)
    assert np.allclose(X[:, 0], X[:, -1])
    assert np.allclose(Y[:, 0], Y[:, -1])
    assert np.ptp(Z[0, :]) == pytest.approx(0.0)
    assert np.allclose(X ** 2 + Y ** 2 + Z ** 2, 1.0)


# ── one continuous surface ──────────────────────────────────────────────────
def test_the_warped_surface_is_closed(surface):
    """The whole point: no caps cut out, no seam, no hole at the poles."""
    surface.reset()
    surface.advance_to(1.0)
    assert surface.is_closed() is True


def test_the_surface_carries_its_poles(surface):
    """A grid of cell centres would miss them; this mesh includes them."""
    assert surface.theta[0] == pytest.approx(0.0)
    assert surface.theta[-1] == pytest.approx(math.pi)
    assert surface.phi[0] == pytest.approx(0.0)
    assert surface.phi[-1] == pytest.approx(2.0 * math.pi)


def test_zero_field_is_an_exact_sphere():
    s = WarpedSphere(**COARSE)
    s.reset()
    s.sim._sim.u[:] = 0.0
    s.sim._sim.u_prev[:] = 0.0
    r = s.radius()
    assert np.allclose(r, s.shells.r_mid)


# ── the warp is a display, and an honest one ────────────────────────────────
def test_displacement_preserves_sign_everywhere(surface):
    surface.reset()
    surface.advance_to(1.3)
    u = surface.field()
    d = surface.displacement()
    assert np.all(np.sign(d) == np.sign(u))


def test_displacement_is_strictly_increasing_in_the_field(surface):
    """tanh is monotone, so amplitude *ordering* survives the display."""
    surface.reset()
    surface.advance_to(1.3)
    u = surface.field().ravel()
    d = surface.displacement().ravel()
    order = np.argsort(u)
    diffs = np.diff(d[order])
    assert np.all(diffs >= -1e-15)


def test_displacement_is_bounded_by_the_gap(surface):
    surface.reset()
    surface.advance_to(0.4)
    assert np.max(np.abs(surface.displacement())) < surface.shells.delta


def test_gain_changes_visibility_not_sign(surface):
    a = WarpedSphere(gain=0.5, **COARSE)
    b = WarpedSphere(gain=3.0, **COARSE)
    for s in (a, b):
        s.reset()
        s.advance_to(0.9)
    da, db = a.displacement(), b.displacement()
    assert np.all(np.sign(da) == np.sign(db))
    assert np.max(np.abs(db)) > np.max(np.abs(da))


def test_gain_and_resolution_are_validated():
    with pytest.raises(ValueError):
        WarpedSphere(gain=0.0, **COARSE)
    with pytest.raises(ValueError):
        WarpedSphere(n_theta=4, n_phi=61)


def test_calibration_uses_the_run_s_own_peak():
    """u_ref is a property of the solve, not a hand-set number."""
    fresh = WarpedSphere(**COARSE)
    assert fresh.reference_amplitude > 0.9          # the launch peaks at 1
    assert fresh.reference_amplitude < 3.0
    assert fresh.t == pytest.approx(0.0)            # calibration rewinds


# ── it never touches a doll ─────────────────────────────────────────────────
def test_the_surface_stays_between_the_dolls(surface):
    r = measure_containment(surface, t_end=RETURN_TIME, frames=40)
    assert r["contained"] is True
    assert r["closest_approach_inner"] > 0.0
    assert r["closest_approach_outer"] > 0.0
    assert r["r_min"] > r["r_inner"]
    assert r["r_max"] < r["r_outer"]


# ── the focus ───────────────────────────────────────────────────────────────
def test_the_deepest_warp_is_at_the_antipode(surface):
    r = measure_focus(surface, frames=90)
    assert r["distance_error"] < 1e-9               # exactly the antipode
    assert r["time_error"] < 0.25                   # early by ~ the pulse width


@pytest.mark.parametrize("width", [0.24, 0.12])
def test_the_focus_time_bias_is_the_pulse_width(width):
    """Narrower pulse, earlier-arrival bias shrinks — it is not a solver error."""
    wide = measure_focus(WarpedSphere(pulse_width=0.24, **COARSE), frames=120)
    narrow = measure_focus(WarpedSphere(pulse_width=width, **COARSE), frames=120)
    assert narrow["time_error"] <= wide["time_error"] + 1e-12
    assert narrow["peak_time"] < ANTIPODAL_TIME


def test_the_focus_reverses_sign(surface):
    """It arrives as a mound and immediately inverts — the surface rings."""
    surface.reset()
    surface.advance_to(3.05)
    out = surface.excursion()["outward_fraction"]
    surface.advance_to(ANTIPODAL_TIME)
    inward = surface.excursion()["inward_fraction"]
    assert out > 0.4
    assert inward > 0.4


# ── the rings are on the surface ────────────────────────────────────────────
def test_front_distances_track_the_front(surface):
    surface.reset()
    surface.advance_to(0.8)
    early = surface.front_distances()["crest_distance"]
    surface.advance_to(1.6)
    later = surface.front_distances()["crest_distance"]
    assert later > early
    assert 0.0 < early < math.pi


def test_a_geodesic_circle_lies_on_the_warped_surface(surface):
    surface.reset()
    surface.advance_to(1.1)
    d = 0.9
    ring = surface.geodesic_circle(d, n=64)
    radii = np.linalg.norm(ring, axis=1)
    assert np.ptp(radii) < 1e-12                    # one radius, as it must be
    assert radii[0] == pytest.approx(float(surface.radius_at_distance(d)))
    # ...and that radius is the surface's own there
    cos_d = ring @ surface.source_direction / radii
    assert np.allclose(np.arccos(np.clip(cos_d, -1, 1)), d)


def test_marker_sits_on_the_surface(surface):
    surface.reset()
    surface.advance_to(0.6)
    p = surface.marker(math.pi)
    assert np.linalg.norm(p) == pytest.approx(
        float(surface.radius_at_distance(math.pi)))


# ── the field sampler the renderer leans on ─────────────────────────────────
def test_field_at_matches_the_grid_field(surface):
    """field_at must agree with BareSphereSim's own grid, exactly."""
    sim = surface.sim
    sim.reset()
    sim.advance_to(0.7)
    TH, PH = np.meshgrid(sim.theta, sim.phi, indexing="ij")
    assert np.allclose(sim.field_at(TH, PH), sim.u)


def test_field_at_distance_is_the_radial_solve(surface):
    sim = surface.sim
    sim.reset()
    sim.advance_to(0.5)
    d = np.array([0.0, 0.3, 1.2, math.pi])
    assert np.allclose(sim.field_at_distance(d), sim._sim.sample(d))
    assert sim.sample(0.5 * math.pi, 0.0) == pytest.approx(
        float(sim.field_at(0.5 * math.pi, 0.0)))
