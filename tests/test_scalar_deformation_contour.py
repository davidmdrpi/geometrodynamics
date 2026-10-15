"""Optional PyVista regression: unknown arrivals must never create NaN geometry.

No render window is created. Skips on installations without PyVista.
"""
import numpy as np
import pytest

pv = pytest.importorskip('pyvista')
from scripts.klein_wormhole_explorer import reached_contour


def test_contour_excludes_unreached_cells():
    mesh = pv.Sphere(theta_resolution=24, phi_resolution=24)
    values = mesh.points[:, 2].copy()
    values[mesh.points[:, 0] > .1] = np.nan
    mesh['arrival'] = values
    front = reached_contour(mesh, 0.0)
    assert front.n_points > 0
    assert np.isfinite(front.points).all()
    assert np.isfinite(front.bounds).all()
    assert np.max(front.points[:, 0]) <= .1 + 1e-7
    mesh['arrival'][:] = np.nan
    assert reached_contour(mesh, 0.0).n_points == 0
