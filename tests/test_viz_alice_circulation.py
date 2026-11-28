"""Charge without charge on one antipodally closed sheet with two Alice handles: geometry of the schematic."""
import math

import numpy as np
import pytest

from geometrodynamics.viz import alice_circulation as ac
from geometrodynamics.viz.klein_handles import Neck


def test_identification_is_the_antipodal_map():
    p = np.array([0.3, -1.1])
    assert np.allclose(ac.identify(p), -p)
    assert np.allclose(ac.identify(ac.identify(p)), p)


def test_each_handle_joins_antipodal_points_of_one_surface():
    mouths = ac.mouths_on_single_surface()
    assert len(mouths) == 2
    flat = [m for pair in mouths for m in pair]
    for a, b in mouths:
        assert np.allclose(a, -b)                     # an antipodal pair
    for i in range(4):                                # four distinct mouths
        for j in range(i + 1, 4):
            assert np.linalg.norm(flat[i] - flat[j]) > 2 * ac.RIM


def test_throat_is_the_finite_mouth_neck():
    nk = Neck(ac.A_MOUTH)
    assert ac.column_radius(0.0) == pytest.approx(nk.b)
    assert ac.column_radius(ac.HEIGHT / 2 - ac.DEPTH) == pytest.approx(float(nk.r(nk.S)))


def test_charges_are_sourceless():
    q = ac.mouth_charges()
    for sheet in ("top", "bottom"):
        assert q[(sheet, 1)] + q[(sheet, 2)] == 0     # Gauss: no net charge on either sheet
    for handle in (1, 2):
        assert q[("top", handle)] == -q[("bottom", handle)]   # flux through a throat: opposite ends


def test_worldline_closes_continuously():
    s = np.linspace(0, 1, 4001)
    P = np.array([ac.worldline(x).xyz for x in s])
    assert np.allclose(P[0], P[-1], atol=1e-9)
    steps = np.linalg.norm(np.diff(P, axis=0), axis=1)
    assert steps.max() < 0.05                          # no jump between segments


def test_both_time_orientations_on_one_circuit():
    assert ac.orientation_flips() == 2                 # one reversal at each Alice handle
    s = np.linspace(0, 1, 2000, endpoint=False)
    pts = [ac.worldline(x) for x in s]
    assert {p.orientation for p in pts} == {+1, -1}
    assert all(p.orientation == +1 for p in pts if p.sheet == "top")
    assert all(p.orientation == -1 for p in pts if p.sheet == "bottom")


def test_render_smoke(tmp_path):
    ac.render(path_png=tmp_path / "still.png", frames=4, balls=3, dpi=40)
    assert (tmp_path / "still.png").stat().st_size > 10_000


def test_handles_are_the_seam_type_that_admits_threading_flux():
    for s_t, s_f in ac.seam_signs().values():
        assert s_t == -1                               # an Alice handle reverses time orientation
        assert ac.flux_through_tube_allowed(s_t, s_f)  # ... and must conjugate charge to carry the flux
    assert not ac.flux_through_tube_allowed(-1, +1)


def test_orientation_flips_at_the_waists():
    s = np.linspace(0, 1, 8000, endpoint=False)
    pts = [ac.worldline(x) for x in s]
    flips = [pts[i].xyz for i in range(len(pts)) if pts[i].orientation != pts[i - 1].orientation]
    assert len(flips) == 2
    for xyz in flips:
        assert abs(xyz[2]) < 0.01                      # at the narrowest point of a throat
