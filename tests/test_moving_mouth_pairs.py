"""Independent force, Noether, constraint, and geometry checks."""
from dataclasses import replace

import numpy as np
import pytest

from geometrodynamics.waves.antipodal_kick import Parameters
from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs
from geometrodynamics.viz.klein_handles import build_surface, is_orientable, euler_characteristic
from geometrodynamics.viz.moving_mouth_pairs import MovingPairGeometry, axis_rotation


def rotate(x, axis, angle):
    x = np.asarray(x)
    return x*np.cos(angle)+np.cross(axis, x)*np.sin(angle)+axis*(axis@x)*(1-np.cos(angle))


def test_pair_torque_is_negative_hamiltonian_gradient():
    model = MovingMouthPairs(Parameters(cutoff=10, mass=.1))
    s = model.initial()
    model.drift(s, 2.9)
    s.q[:] = [[.03, -.01], [.02, .04]]
    torque = model.interaction(s)[3].sum(axis=1)
    for h in range(2):
        for i, axis in enumerate(np.eye(3)):
            plus, minus = s.copy(), s.copy()
            plus.x[h] = rotate(s.x[h], axis, 1e-6)
            minus.x[h] = rotate(s.x[h], axis, -1e-6)
            derivative = (model.energy(plus)-model.energy(minus))/2e-6
            assert torque[h, i] == pytest.approx(-derivative, abs=2e-10)


def test_reciprocal_kick_conserves_three_momenta_and_mouth_ledger():
    model = MovingMouthPairs(Parameters(cutoff=12, mass=.1))
    s = model.initial()
    model.drift(s, 2.9)
    s.x[0] = np.array([.2, .3, np.sqrt(.87)])
    s.x[1] = np.array([np.sqrt(.8), -.4, .2])
    before = model.wave_angular(s)+s.angular.sum(axis=0)
    model.kick(s, .17)
    np.testing.assert_allclose(model.wave_angular(s)+s.angular.sum(axis=0), before, atol=1e-16, rtol=0)
    np.testing.assert_allclose(s.angular, s.mouth_impulse.sum(axis=1), atol=1e-17)
    np.testing.assert_allclose(np.sum(s.x*s.angular, axis=1), 0, atol=1e-17)
    assert np.linalg.norm(s.angular.sum(axis=0)) > 1e-6


def test_independent_free_geodesics_and_reversible_integrator():
    model = MovingMouthPairs(Parameters(cutoff=8, mass=.1, coupling=0))
    s = model.initial()
    s.angular = np.array([[.02, 0, 0], [0, 0, .03]])
    model.drift(s, .5)
    np.testing.assert_allclose(s.x[0], [0, -np.sin(.1), np.cos(.1)], atol=1e-15)
    np.testing.assert_allclose(s.x[1], [np.cos(.15), np.sin(.15), 0], atol=1e-15)
    model = MovingMouthPairs(Parameters(cutoff=8, mass=.1))
    s = model.initial()
    before = s.copy()
    for _ in range(20):
        model.step(s, .01)
    for _ in range(20):
        model.step(s, -.01)
    for key in before.__dataclass_fields__:
        np.testing.assert_allclose(getattr(s, key), getattr(before, key), atol=2e-14, rtol=0)


def test_mirror_null_controls_and_time_refinement():
    par = Parameters(cutoff=8, mass=.1)
    model = MovingMouthPairs(par)
    plus, coarse = model.evolve(dt=.012, end=.8)
    minus, _ = model.evolve(model.initial(-.7), dt=.012, end=.8)
    np.testing.assert_allclose(minus['x'], plus['x']*np.array([1, -1, 1]), atol=1e-13)
    np.testing.assert_allclose(minus['angular'], plus['angular']*np.array([-1, 1, -1]), atol=1e-13)
    _, fine = model.evolve(dt=.006, end=.8)
    assert fine['max_energy_relative'] < .4*coarse['max_energy_relative']
    for obj, state in [(model, model.initial(empty=True)),
                       (MovingMouthPairs(replace(par, coupling=0)), model.initial())]:
        data, _ = obj.evolve(state, dt=.01, end=.3)
        np.testing.assert_array_equal(data['angular'], 0.)
        np.testing.assert_array_equal(data['q'], 0.)


def test_moving_mesh_keeps_rims_gluing_and_fixed_far_landmarks():
    surf = build_surface(n_sphere=600, n_ring=16)
    original_tris = surf.tris.copy()
    geom = MovingPairGeometry(surf)
    axes = np.array([rotate([0, 0, 1], np.array([1, 0, 0]), .05),
                     rotate([1, 0, 0], np.array([0, 0, 1]), .08)])
    points = geom.coordinates(axes)
    for h, ids in enumerate(surf.tube_ids):
        rot = axis_rotation(geom.axes[h], axes[h])
        np.testing.assert_allclose(points[ids], surf.points[ids]@rot.T, atol=5e-16, rtol=0)
        for sign, rim in zip([1, -1], surf.rims[h]):
            np.testing.assert_allclose(points[rim].mean(axis=0), sign*surf.neck.ca*axes[h], atol=5e-16)
    far = geom.sphere[np.all(geom.weights == 0, axis=1)]
    np.testing.assert_allclose(points[far], surf.points[far], atol=3e-16, rtol=0)
    np.testing.assert_array_equal(surf.tris, original_tris)
    assert not is_orientable(surf.tris)
    assert euler_characteristic(surf.tris) == -2
    assert geom.minimum_sphere_orientation_ratio(points) > 0
    deformed = geom.deform(points, points[geom.sphere, 1], amplitude=1)
    for ids in surf.tube_ids:
        for rim in [ids[0], ids[-1]]:
            expected = points[rim]*(1+.15*np.tanh(points[rim, 1]))[:, None]
            np.testing.assert_allclose(deformed[rim], expected, atol=3e-16)
    # Moving coordinates must be used to sample an Eulerian field u(n)=n_y.
    assert np.max(abs(points[geom.sphere, 1]-surf.points[geom.sphere, 1])) > .01
    with pytest.raises(ValueError):
        geom.coordinates(np.array([[0., 0., 1.], [0., 0., 1.]]))
