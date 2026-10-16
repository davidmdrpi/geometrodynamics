"""Independent invariants of the scalar/attachment coupling and focus controls."""
from dataclasses import replace

import numpy as np
import pytest

from geometrodynamics.waves.antipodal_kick import AntipodalKick, Parameters


def rotation(x, axis, angle):
    return (x*np.cos(angle)+np.cross(axis, x)*np.sin(angle)
            + axis*(axis@x)*(1-np.cos(angle)))


def test_free_phase_space_refocus_and_energy():
    model = AntipodalKick(Parameters(cutoff=12, coupling=0))
    for eta in [0., -.7, .7]:
        state = model.initial(eta)
        target_a, target_p = model.target_state(eta)
        assert model.wave_energy(state.a, state.p) == pytest.approx(.1)
        model.drift(state, np.pi)
        np.testing.assert_allclose(state.a, target_a, atol=2e-16, rtol=0)
        np.testing.assert_allclose(state.p, target_p, atol=2e-15, rtol=0)


def test_rotational_derivatives_against_cartesian_finite_difference():
    model = AntipodalKick(Parameters(cutoff=10))
    a, _ = model.target_state(.7)
    x = np.array([.8, .4, np.sqrt(.2)])
    h = 1e-6
    w = model.footprint(x)
    for axis, d in zip(np.eye(3), model.generators):
        values = [np.vdot(model.footprint(rotation(x, axis, sign*h)), a).real
                  for sign in [1, -1]]
        assert np.vdot(w, d@a).real == pytest.approx((values[0]-values[1])/(2*h), abs=1e-9)


def test_torque_is_negative_energy_gradient():
    model = AntipodalKick(Parameters(cutoff=10))
    s = model.initial()
    model.drift(s, np.pi-.17)
    s.q = .03
    h = 1e-6
    torque = model.interaction(s)[3]
    for i, axis in enumerate(np.eye(3)):
        states = [s.copy(), s.copy()]
        for state, sign in zip(states, [1, -1]):
            state.x = rotation(s.x, axis, sign*h)
        gradient = (model.energy(states[0])-model.energy(states[1]))/(2*h)
        assert torque[i] == pytest.approx(-gradient, abs=1e-10)


def test_interaction_kick_exchanges_all_three_momenta():
    model = AntipodalKick(Parameters(cutoff=12))
    s = model.initial()
    model.drift(s, 2.9)
    s.x = np.array([.8, .4, np.sqrt(.2)])
    j = model.wave_angular(s)+s.angular
    model.kick(s, .13)
    assert np.linalg.norm(s.angular) > 1e-5
    np.testing.assert_allclose(model.wave_angular(s)+s.angular, j, atol=1e-16, rtol=0)
    assert abs(s.x@s.angular) < 1e-16


def test_dephasing_preserves_modal_energy_and_angular_momentum():
    model = AntipodalKick(Parameters(cutoff=12))
    s, d = model.initial(), model.initial(dephased=True)
    for l in range(13):
        mask = model.ell == l
        es = abs(s.p[mask])**2+model.omega[mask]**2*abs(s.a[mask])**2
        ed = abs(d.p[mask])**2+model.omega[mask]**2*abs(d.a[mask])**2
        np.testing.assert_allclose(es, ed, atol=1e-17, rtol=1e-13)
    np.testing.assert_allclose(model.wave_angular(s), model.wave_angular(d), atol=1e-16)


def test_uncoupled_attachment_stays_at_rest():
    model = AntipodalKick(Parameters(cutoff=8, coupling=0))
    a, stats = model.evolve(dt=.02)
    np.testing.assert_array_equal(a['angular'], 0)
    np.testing.assert_array_equal(a['q'], 0)
    np.testing.assert_allclose(a['x'], np.tile([1, 0, 0], (len(a['t']), 1)), atol=0)
    assert stats['max_energy_relative'] < 1e-12


def test_symmetric_and_mirrored_controls_and_second_order_energy():
    model = AntipodalKick(Parameters(cutoff=10))
    _, symmetric = model.evolve(model.initial(0), dt=.008)
    assert np.linalg.norm(symmetric['final_impulse']) < 1e-14
    _, plus = model.evolve(dt=.008)
    _, minus = model.evolve(model.initial(-.7), dt=.008)
    np.testing.assert_allclose(plus['final_impulse'], -np.array(minus['final_impulse']), atol=1e-13)
    _, fine = model.evolve(dt=.004)
    assert fine['max_energy_relative'] < .4*plus['max_energy_relative']
    assert fine['max_angular_relative'] < 1e-12
    assert fine['max_sphere_error'] < 1e-12


def test_invalid_parameters():
    with pytest.raises(ValueError):
        AntipodalKick(replace(Parameters(), mass=0))
