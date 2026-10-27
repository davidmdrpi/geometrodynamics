import numpy as np
import pytest
from scipy.integrate import simpson

from geometrodynamics.waves.directional_flash_response import (
    coasting_focus_audit, electric, electric_integral, power_audit,
    rim_geometry, rim_modes,
)
from geometrodynamics.waves.pair_birth import PairBirth, Parameters, rotate


@pytest.fixture
def flash():
    model = PairBirth(Parameters(cutoff=12))
    state = model.initial()
    q, p = model.flash(state)
    return model, state, q.sum(axis=0), p.sum(axis=0)


def test_electric_is_minus_time_derivative_and_is_tangent(flash):
    model, state, q, p = flash
    t, h = 2.7, 1e-5
    _, pt = model.free_wave(q, p, t)
    e = electric(model, pt, state.x)
    np.testing.assert_allclose(e, electric_integral(model, q, p, state.x, t-h, t+h)/(2*h),
                               atol=2e-11, rtol=1e-8)
    np.testing.assert_allclose(np.sum(e*state.x, axis=1), 0, atol=1e-16)


def test_integral_matches_independent_time_quadrature(flash):
    model, state, q, p = flash
    times = np.linspace(2.4, 3.8, 501)
    values = [electric(model, model.free_wave(q, p, t)[1], -state.x) for t in times]
    actual = simpson(values, x=times, axis=0)
    np.testing.assert_allclose(actual, electric_integral(model, q, p, -state.x, times[0], times[-1]),
                               atol=2e-11, rtol=1e-8)
    np.testing.assert_allclose(electric_integral(model, q, p, -state.x, times[-1], times[0]),
                               -actual, atol=2e-11, rtol=1e-8)


def test_source_superposition_sign_and_time_reversal(flash):
    model, state, q, p = flash
    copies, momenta = model.flash(state)
    for t in [np.pi, 2*np.pi, 3.8]:
        total = electric(model, model.free_wave(q, p, t)[1], -state.x)
        pieces = sum(electric(model, model.free_wave(qq, pp, t)[1], -state.x)
                     for qq, pp in zip(copies, momenta))
        np.testing.assert_allclose(total, pieces, atol=2e-16)
        np.testing.assert_allclose(electric(model, model.free_wave(-q, -p, t)[1], -state.x), -total,
                                   atol=2e-16)
        np.testing.assert_allclose(electric(model, model.free_wave(q, p, -t)[1], -state.x), -total,
                                   atol=2e-16)


def test_same_footprint_gives_reciprocal_power(flash):
    model, state, q, p = flash
    rng = np.random.default_rng(349)
    for _ in range(8):
        x = rng.normal(size=(2, 3)); x /= np.linalg.norm(x, axis=1)[:, None]
        v = rng.normal(size=(2, 3)); v -= np.sum(x*v, axis=1)[:, None]*x
        pt = model.free_wave(q, p, rng.uniform(-8, 8))[1]
        audit = power_audit(model, pt, x, v)
        assert abs(audit['residual']) < 2e-16


def test_rim_m1_recovers_uniform_transported_vector():
    _, er, ep, phi, _ = rim_geometry([0, 0, 1], .3, first=[1, 0, 0])
    # Parallel transport of the vector (2,-.7) from the centre.
    radial = 2*np.cos(phi)-.7*np.sin(phi)
    azimuthal = -2*np.sin(phi)-.7*np.cos(phi)
    fields = radial[:, None]*er + azimuthal[:, None]*ep
    modes = rim_modes(fields, er, ep, phi)
    np.testing.assert_allclose(modes['transported_mean'], [2, -.7], atol=8e-16)
    assert abs(modes['radial_m0']) < 3e-16
    assert abs(modes['azimuthal_m0']) < 3e-16
    # Centred radial m=0 stress cannot generate a tangential mean.
    np.testing.assert_allclose(rim_modes(er, er, ep, phi)['transported_mean'], 0, atol=2e-16)


def test_vector_field_projection_rotates_covariantly():
    model = PairBirth(Parameters(cutoff=12))
    axis = np.array([1., 2., 3.]); axis /= np.linalg.norm(axis)
    plane = rotate(np.array([[0., 0., 1.], [1., 0., 0.]]), axis, .73)
    original, turned = model.initial(), model.initial(plane=plane)
    responses = []
    for state in [original, turned]:
        q, p = model.flash(state)
        pt = model.free_wave(q.sum(axis=0), p.sum(axis=0), np.pi)[1]
        responses.append(electric(model, pt, -state.x))
    np.testing.assert_allclose(responses[1], rotate(responses[0], axis, .73), atol=2e-15)


def test_no_phantom_capture_in_standard_first_flight():
    result = coasting_focus_audit(.3, 1/3)
    assert len(result['events']) == 4
    assert all(e['event'] == 'empty_focus' for e in result['events'])
    assert all(e['time'] < result['contact_time'] for e in result['events'])


def test_capture_stops_free_surface_route_and_contact_is_excluded():
    result = coasting_focus_audit(.3, .8)
    assert [e['target'] for e in result['events']] == [1, 0]
    assert all(e['return_number'] == 1 and e['offset'] < .3 for e in result['events'])
    # Contact at exactly pi is not an open-mouth capture.
    result = coasting_focus_audit(.3, (np.pi-.6)/np.pi)
    assert result['events'] == []


@pytest.mark.parametrize('angle', [.05, .16, .3, .5])
def test_first_return_capture_band_from_independent_geometry(angle):
    low, high = 1-3*angle/np.pi, 1-2*angle/np.pi
    for speed in [low-1e-5, (low+high)/2, high+1e-5]:
        events = coasting_focus_audit(angle, speed)['events']
        captures = [e for e in events if e['return_number']==1 and e['target'] is not None]
        assert bool(captures) == (low <= speed < high)


def test_slow_mouth_can_recapture_its_own_second_return():
    result = coasting_focus_audit(.3, .03)
    captured = [e for e in result['events'] if e['target'] is not None]
    assert len(captured) == 2
    assert all(e['return_number']==2 and e['target']==e['source'] for e in captured)
