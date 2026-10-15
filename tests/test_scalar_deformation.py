"""Scalar-to-shape fidelity, gluing, and external-time sampling."""

import numpy as np
import pytest

from geometrodynamics.viz.klein_handles import build_surface, default_handles, HandleWaveSim
from geometrodynamics.viz.scalar_deformation import ScalarDeformation


@pytest.mark.parametrize("kinds", [(), (True,), (False,), (True, False)])
def test_signed_radial_graph_and_shared_rims(kinds):
    surf = build_surface(handles=default_handles(len(kinds), kinds), n_sphere=1600, n_ring=24)
    cut = surf.cut()
    warp = ScalarDeformation(surf, cut)
    u = np.linspace(-0.4, 0.4, surf.n_vertices)[cut.orig]
    points = warp.coordinates(u, 1.0, 3.0)
    np.testing.assert_array_equal(warp.coordinates(u, 1.0, 0.0), warp.reference)
    np.testing.assert_allclose(warp.coordinates(-u, 1.0, 3.0) - warp.reference,
                               -(points-warp.reference), atol=3e-16)
    sphere = surf.on_sphere[cut.orig]
    np.testing.assert_allclose(np.linalg.norm(points[sphere], axis=1),
                               1 + warp.excursion*np.tanh(3*u[sphere]), atol=5e-16)
    for ids in surf.tube_ids:
        # Both glued ends share precisely the sphere's radial displacement,
        # including the far rim's gluing permutation.
        for ring in (ids[0], ids[-1]):
            np.testing.assert_allclose(warp.vectors[ring], warp.excursion*surf.points[ring])
    np.testing.assert_allclose(points[surf.n_vertices:], points[cut.orig[surf.n_vertices:]])
    assert np.max(np.linalg.norm(points-warp.reference, axis=1)) <= warp.excursion


def test_breathing_at_waist_and_fixed_amplitude_scale():
    surf = build_surface(handles=default_handles(1), n_sphere=1600, n_ring=24)
    warp = ScalarDeformation(surf, surf.cut())
    ids = surf.tube_ids[0][surf.waist]
    np.testing.assert_allclose(np.linalg.norm(warp.vectors[ids], axis=1),
                               warp.excursion*surf.neck.b)
    u = np.ones(len(warp.reference))*.1
    big = warp.coordinates(u, 1, 3)-warp.reference
    small = warp.coordinates(u/10, 1, 3)-warp.reference
    np.testing.assert_allclose(small, big*np.tanh(.03)/np.tanh(.3), atol=3e-16)


def test_clock_cut_faces_are_not_welded_or_averaged():
    surf = build_surface(handles=default_handles(1), n_sphere=1600, n_ring=24)
    cut = surf.cut()
    warp = ScalarDeformation(surf, cut)
    u = np.zeros(len(cut.orig))
    copies = np.flatnonzero(cut.handle == 0)
    u[copies] = .2  # samples of the same intrinsic vertices at different times
    p = warp.coordinates(u, 1, 4)
    np.testing.assert_array_equal(p[cut.orig[copies]], warp.reference[copies])
    assert np.min(np.linalg.norm(p[copies]-p[cut.orig[copies]], axis=1)) > 0


def test_real_external_frames_drive_shape_without_changing_the_solver():
    surf = build_surface(handles=default_handles(1), n_sphere=1600, n_ring=24)
    sim = HandleWaveSim(surf, offsets=[-.4], n_frames=51, t_max=1.8).run()
    warp = ScalarDeformation(surf, sim.cut)
    original_points = surf.points.copy()
    original_frames = sim.frames.copy()
    original_K = sim.K.copy()
    k = sim.frame_at(.8)
    raw = sim.frames[k].astype(float)/sim.amplitude
    p = warp.coordinates(sim.frames[k], sim.amplitude, 6)
    sphere = surf.on_sphere[sim.cut.orig]
    np.testing.assert_allclose(np.linalg.norm(p[sphere], axis=1),
                               1+warp.excursion*np.tanh(6*raw[sphere]))
    warp.coordinates(sim.frames[k], sim.amplitude, 20)
    np.testing.assert_array_equal(surf.points, original_points)
    np.testing.assert_array_equal(sim.frames, original_frames)
    assert (sim.K != original_K).nnz == 0
    assert sim.energy_drift < 1e-10


def test_invalid_display_input_is_rejected():
    surf = build_surface(handles=(), n_sphere=100)
    warp = ScalarDeformation(surf, surf.cut())
    u = np.zeros(len(warp.reference))
    for amplitude, gain in [(0, 1), (1, -1), (np.nan, 1), (1, np.inf)]:
        with pytest.raises(ValueError):
            warp.coordinates(u, amplitude, gain)
    with pytest.raises(ValueError):
        warp.coordinates(u[:-1], 1)
