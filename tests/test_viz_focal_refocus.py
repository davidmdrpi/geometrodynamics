"""
Tests for the antipodal refocus of a continuous trace-free deformation.

The claim under test is not "a wave refocuses" — that is the scalar's story
too — but what a **spin-2** field can and cannot do when it does:

* it cannot sit on its own focal point, because ``h = sin²d·q`` vanishes at
  both poles for every ``q``, so the density piles into a *ring*;
* the ring's radius is the pulse width, not a numerical floor;
* the finite amplification is **shared with the scalar**, so it is tested here
  as a *negative* result — it is not a spin-2 protection mechanism;
* a material patch changes shape without changing size, and reports the local
  stretch axis in the limit of a small patch;
* and the area law's second-order residual is largest exactly on the focal
  ring, which is where the linear description runs out.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.viz.embedded_wave import (
    EmbeddedTidalSurface,
    MaterialPatch,
    measure_patch_area_invariance,
    measure_patch_shape_history,
    measure_patch_size_convergence,
    measure_where_the_area_law_fails,
)
from geometrodynamics.viz.spin2_tidal import (
    ANTIPODAL_TIME,
    Spin2WaveSim,
    TidalField,
    measure_focal_energy,
)

PULSE_WIDTH = 0.18


@pytest.fixture(scope="module")
def surface():
    return EmbeddedTidalSurface(
        sim=Spin2WaveSim(n=900, pulse_width=PULSE_WIDTH),
        n_theta=61, n_phi=91)


@pytest.fixture(scope="module")
def focal():
    return measure_focal_energy(TidalField(sim=Spin2WaveSim(
        n=900, pulse_width=PULSE_WIDTH)))


# ── the energy density ──────────────────────────────────────────────────────
def test_the_density_is_the_kinetic_invariant_not_the_conserved_one():
    """``∫ρ_E dA`` is *half* the energy and oscillates — the invariant is not.

    Worth pinning down: reading the kinetic integral as a conservation check
    reports a drift of order one where there is none.
    """
    sim = Spin2WaveSim(n=900, pulse_width=PULSE_WIDTH)
    totals = []
    for i in range(60):
        sim.advance_to((i + 1) * ANTIPODAL_TIME / 60)
        totals.append(sim.total_energy_measure())
    swing = (max(totals) - min(totals)) / max(totals)
    assert swing > 0.1                        # the kinetic half really swings
    assert abs(sim.energy_drift()) < 1e-8     # ...while the invariant holds


def test_the_focus_is_a_node(focal):
    """A spin-2 field cannot sit on its own focal point."""
    assert focal["antipode_over_peak"] < 1e-4
    assert focal["ring_radius"] > 1e-3
    assert focal["concentrates_in_a_ring"] is True


def test_the_focal_time_is_the_antipodal_transit(focal):
    assert focal["focal_time"] == pytest.approx(ANTIPODAL_TIME, abs=0.25)
    assert focal["focal_distance"] < ANTIPODAL_TIME


def test_the_ring_radius_tracks_the_pulse_width():
    """The hole is the wave's own length scale, not a grid artefact."""
    ratios = []
    for w in (0.24, 0.12):
        r = measure_focal_energy(TidalField(sim=Spin2WaveSim(
            n=900, pulse_width=w)))
        ratios.append(r["ring_radius"] / w)
    assert all(0.7 < x < 1.2 for x in ratios)
    assert abs(ratios[0] - ratios[1]) < 0.15


def test_the_amplification_is_finite(focal):
    assert 1.2 < focal["amplification"] < 4.0
    assert abs(focal["invariant_drift"]) < 1e-8


def test_the_amplification_is_not_a_spin_two_effect():
    """A negative result, kept as a test so it cannot quietly become a claim.

    The scalar refocuses with the same ``O(1)`` gain, so the factor says
    nothing about the spin — only the node and the ring do.
    """
    from geometrodynamics.viz.spin2_tidal import (
        measure_amplification_is_not_protection,
    )

    r = measure_amplification_is_not_protection(widths=(0.24, 0.12), frames=200)
    t_lo, t_hi = r["tensor_amplification_range"]
    s_lo, s_hi = r["scalar_amplification_range"]
    assert 1.2 < t_lo and t_hi < 4.0
    assert 1.2 < s_lo and s_hi < 4.0
    assert abs(t_hi - s_hi) < 1.0
    assert r["amplification_is_not_a_spin_2_effect"] is True
    assert r["but_the_focal_node_is"] is True


# ── the material patch ──────────────────────────────────────────────────────
def test_the_patch_shape_is_read_in_its_own_plane(surface):
    """The drawn patch is a curved cap; its shape only means anything in-plane.

    Projecting against the *radial* direction instead leaves part of the
    surface's tilt in the answer and rotates the measured long axis out of the
    surface — which is exactly what it did before this was fixed.
    """
    surface.reset()
    surface.advance_to(3.08)
    p = MaterialPatch(surface, centre_distance=math.pi - 0.94 * PULSE_WIDTH,
                      centre_azimuth=1.5 * math.pi, radius=0.10)
    sh = p.shape()
    u, v, n = sh["long_axis"], sh["short_axis"], sh["normal"]
    for a, b in ((u, v), (u, n), (v, n)):
        assert float(a @ b) == pytest.approx(0.0, abs=1e-9)
    assert sh["flatness"] < 0.25              # genuinely a flattened disc
    assert sh["aspect_ratio"] > 1.0


def test_the_patch_changes_shape_without_changing_size():
    r = measure_patch_area_invariance(gain=1e-2, frames=60)
    assert r["area_is_invariant"] is True
    assert r["relative_area_swing"] < 1e-5
    assert r["max_aspect_ratio_at_display_gain"] > 1.02


def test_the_patch_long_axis_is_the_stretch_eigenvector(surface):
    r = measure_patch_shape_history(surface, centre_distance=1.20, frames=80)
    assert r["aligns_with_the_stretch_axis"] is True
    assert r["long_axis_alignment"] > 0.98
    assert r["max_aspect_ratio"] > 1.0


def test_near_the_focus_alignment_is_a_question_of_patch_size():
    """A patch straddling the focal ring averages a sign change; a small one does not."""
    r = measure_patch_size_convergence(
        centre_distance=math.pi - 0.94 * PULSE_WIDTH,
        radii=(0.24, 0.12, 0.05), frames=90)
    assert r["improves_as_the_patch_shrinks"] is True
    assert r["converges_to_the_eigenvector"] is True
    assert r["smallest_patch_alignment"] > r["worst_alignment"]


def test_the_focal_patch_distorts_harder_than_the_smooth_one():
    r = measure_where_the_area_law_fails(frames=90)
    assert r["the_focus_distorts_harder"] is True
    assert r["focal_max_aspect_ratio"] > r["smooth_max_aspect_ratio"]
    assert r["focal_over_smooth_distortion"] > 3.0


def test_the_area_residual_is_second_order_in_the_gain():
    """So the failure at the focus is the linear description running out.

    Not a numerical defect: halve the gain and the residual quarters.
    """
    r = measure_where_the_area_law_fails(frames=60)
    assert r["residual_is_second_order"] is True
    assert r["residual_exponent"] == pytest.approx(2.0, abs=0.35)
    assert r["the_area_law_fails_first_at_the_focus"] is True
    residuals = r["focal_area_residuals"]
    assert all(b > a for a, b in zip(residuals, residuals[1:]))


def test_patch_geometry_is_finite_and_closed(surface):
    surface.reset()
    surface.advance_to(1.4)
    p = MaterialPatch(surface, centre_distance=1.2, radius=0.12,
                      n_rings=5, n_spokes=32)
    pts, bnd, tris = p.points(), p.boundary(), p.triangles()
    for a in (pts, bnd, tris):
        assert np.all(np.isfinite(a))
    assert bnd.shape == (33, 3)                        # closed ring
    assert np.allclose(bnd[0], bnd[-1])
    assert tris.shape[1:] == (3, 3)
    assert p.area() > 0.0
    # Flat triangles chord the curved cap, so the undeformed area is slightly
    # *under* the spherical one; refining the tiling closes the gap.
    coarse = p.area(gain=0.0) / p.reference_area
    fine = MaterialPatch(surface, centre_distance=1.2, radius=0.12,
                         n_rings=20, n_spokes=160).area(gain=0.0) \
        / p.reference_area
    assert coarse < fine < 1.0
    assert fine == pytest.approx(1.0, rel=1e-3)


def test_a_patch_radius_must_be_sane(surface):
    with pytest.raises(ValueError):
        MaterialPatch(surface, centre_distance=1.0, radius=0.0)
    with pytest.raises(ValueError):
        MaterialPatch(surface, centre_distance=1.0, radius=2.0)
