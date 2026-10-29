import numpy as np
import pytest
from numpy.testing import assert_allclose
from geometrodynamics.waves.maxwell_throat_response import Bridge
from geometrodynamics.waves.einstein_throat_evolution import (
    initial_acceleration, central_history, integrate_central, optical_length, group_delays,
)


@pytest.mark.parametrize('fraction', [0., .5, .9, .99])
def test_parametric_history_matches_independent_second_order_evolution(fraction):
    b = Bridge(charge_fraction=fraction)
    exact = central_history(b, np.linspace(0, .97*np.pi, 101))
    numeric = integrate_central(b, exact['tau'], rtol=3e-12)
    assert_allclose(numeric['radius'], exact['radius'], atol=2e-10*b.b, rtol=3e-8)
    assert_allclose(numeric['velocity'], exact['velocity'], atol=2e-8, rtol=3e-8)
    r, u = exact['radius'], exact['velocity']
    assert_allclose(u*u+1-2*b.mass/r+b.geometric_charge**2/r**2, 0, atol=3e-12)
    assert_allclose(exact['acceleration'][0], initial_acceleration(b, 0), rtol=1e-13)
    assert np.all(exact['null_expansion_plus'][1:] < 0)
    assert_allclose(exact['null_expansion_plus'], exact['null_expansion_minus'])


def test_adm_angular_equation_independently_reproduces_acceleration():
    for fraction in [0., .7]:
        b = Bridge(charge_fraction=fraction)
        x = np.linspace(-b.limit, b.limit, 45)
        r = b.profile(x)['radius']
        first = b.discriminant*np.sinh(2*x)/r
        second = b.mass/r**2-b.geometric_charge**2/r**3
        rho = b.geometric_charge**2/(8*np.pi*r**4)
        for pressure in [-2., 0., 2.]:
            radial, tangent = -rho+pressure, rho
            ricci_theta = 1-first**2-r*second
            trace_stress = radial+2*tangent
            kdot_theta = ricci_theta-8*np.pi*r*r*tangent+4*np.pi*r*r*(trace_stress-rho)
            assert_allclose(-kdot_theta/r, initial_acceleration(b, x, pressure), rtol=2e-12, atol=1e-12)


def test_static_support_exactly_cancels_einstein_radius_acceleration():
    for fraction in [0., .5, .9]:
        for depth in [0., .8, .99]:
            b = Bridge(charge_fraction=fraction, lapse_depth=depth)
            x = np.linspace(-b.limit, b.limit, 91)
            r = b.profile(x)['radius']
            support_pressure = (b.support(x)[1]+b.geometric_charge**2/r**4)/(8*np.pi)
            assert_allclose(initial_acceleration(b, x, support_pressure, static_lapse=True), 0, atol=3e-12)
            # At the waist the lapse-gradient term vanishes for every depth.
            assert_allclose(initial_acceleration(b, 0, static_lapse=True), initial_acceleration(b, 0))


def test_same_constraint_data_do_not_fix_acceleration_without_pressure_law():
    b = Bridge(charge_fraction=.9)
    rho = b.geometric_charge**2/(8*np.pi*b.b**4)
    support = b.support(0)[1]/(8*np.pi)+rho
    free = initial_acceleration(b, 0)
    assert free < 0
    assert_allclose(initial_acceleration(b, 0, support), 0, atol=1e-12)
    assert_allclose(initial_acceleration(b, 0, 2*support), -free, rtol=1e-13)


def test_neutral_singularity_and_charged_inner_horizon_are_different_endpoints():
    b = Bridge()
    with pytest.raises(ValueError):
        central_history(b, np.pi)
    near = central_history(b, np.array([.8, .9, .95, .99])*np.pi)
    assert np.all(np.diff(near['kretschmann']) > 0)
    for f in [.5, .9]:
        b = Bridge(charge_fraction=f)
        end = central_history(b, np.pi)
        assert_allclose(end['radius'], b.geometric_charge**2/b.b, rtol=1e-14)
        assert_allclose(end['tau'], np.pi*b.mass, rtol=1e-14)
        assert end['radius'] > 0 and end['acceleration'] > 0
        with pytest.raises(ValueError):
            central_history(b, 1.01*np.pi)


@pytest.mark.parametrize('depth,omega', [(0., 10.), (.8, 10.), (.99, 8.)])
def test_wigner_smith_delays_are_hermitian_and_probability_weighted(depth, omega):
    b = Bridge(lapse_depth=depth)
    a = group_delays(b, omega, step=1e-4)
    fine = group_delays(b, omega, step=5e-5)
    assert a['hermiticity_error'] < 3e-6
    expected = a['reflection']*a['reflection_delay']+a['transmission']*a['transmission_delay']
    assert_allclose(a['mean_delay'], expected, atol=3e-12)
    assert_allclose(a['mean_delay'], fine['mean_delay'], rtol=2e-6, atol=2e-7)


def test_high_frequency_delay_tends_to_optical_length_and_quadrature_converges():
    b = Bridge(lapse_depth=.8)
    optical = optical_length(b)
    assert_allclose(optical, optical_length(b, order=768), rtol=1e-12)
    a = group_delays(b, 300.)['transmission_delay']
    c = group_delays(b, 600.)['transmission_delay']
    assert abs(c-optical) < abs(a-optical)/3
    assert abs(c-optical)/optical < .001
