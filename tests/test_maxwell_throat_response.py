import numpy as np
import pytest
from numpy.testing import assert_allclose

from geometrodynamics.waves.maxwell_throat_response import (
    Bridge, solve_port, force_matrices, scattering_force_derivative,
    boundary_dtn, field_basis, realify, seam_map,
    static_flux,
)
from geometrodynamics.waves.physical_throat import VacuumThroat


def test_existing_bridge_geometry_and_required_support():
    for a in [.2, .3, .45]:
        b = Bridge(a)
        old = VacuumThroat(a)
        assert_allclose(b.b, old.neck_radius(), rtol=2e-15)
        assert_allclose(2*b.half_length, old.length(), rtol=2e-15)
        assert_allclose(b.profile(b.limit)['radius'], np.sin(a), atol=1e-15)
        x = np.linspace(-b.limit, b.limit, 51)
        stress = b.support(x)
        r = b.profile(x)['radius']
        assert_allclose(stress[0], 0, atol=2e-11)
        assert_allclose(stress[1], -b.b/r**3, rtol=1e-13)
        assert_allclose(stress[2], b.b/(2*r**3), rtol=1e-13)
        for depth in [0., .4, .8]:
            rho, pr, pt = Bridge(a, depth).support(0.)
            assert_allclose(rho+pr, -1/b.b**2, rtol=2e-14)


@pytest.mark.parametrize('depth', [0., .4, .8])
def test_port_flux_reciprocity_and_independent_dtn(depth):
    b = Bridge(lapse_depth=depth)
    for w in [5., 17., 60.]:
        result = solve_port(b, w)
        s = result['scattering']
        assert_allclose(s.T, s, atol=2e-14)
        assert_allclose(s.conj().T@s, np.eye(2), atol=3e-8)
        d = boundary_dtn(b, w)
        independent = -np.linalg.solve(d-1j*w*np.eye(2), d+1j*w*np.eye(2))
        assert_allclose(s, independent, atol=3e-8)
        u, j = field_basis(result, np.array([-b.limit, b.limit]))
        assert_allclose(u[0, 0], 1+s[0, 0], atol=1e-12)
        assert_allclose(j[0, 0], 1j*w*(1-s[0, 0]), atol=2e-10)
        assert_allclose(j[1, 0], 1j*w*s[1, 0], atol=2e-10)


@pytest.mark.parametrize('parameter', ['waist', 'stretch'])
def test_stress_action_force_equals_scattering_phase_derivative(parameter):
    for depth, w in [(0., 10.), (0., 65.), (.6, 17.)]:
        b = Bridge(lapse_depth=depth)
        value = force_matrices(b, w, rtol=2e-12)
        integral = value[parameter+'_force']
        phase = scattering_force_derivative(b, w, parameter=parameter, step=3e-6)
        assert_allclose(integral, phase, rtol=3e-6, atol=2e-8)


def test_coherent_radius_load_is_nonnegative_for_all_phases():
    b = Bridge()
    for w in [4., 12., 40., 70.]:
        result = force_matrices(b, w)
        matrix = result['waist_force']
        assert np.min(np.linalg.eigvalsh(matrix)) > -2e-12
        for phase in np.linspace(0, 2*np.pi, 25):
            incoming = np.array([1., np.exp(1j*phase)])/np.sqrt(2)
            assert np.vdot(incoming, matrix@incoming).real >= -2e-12


def test_deeper_lapse_changes_scattering_at_fixed_spatial_slice():
    a = solve_port(Bridge(), 20.)
    b = solve_port(Bridge(lapse_depth=.8), 20.)
    assert abs(a['transmission']-b['transmission']) > .01


def test_time_reversal_is_real_linear_not_a_complex_phase():
    spatial = np.diag([-1., 1.])
    reverse = seam_map(spatial, True)
    u = np.array([.2+.7j, -.4+.3j])
    actual = reverse @ np.r_[u.real, u.imag]
    expected = spatial @ u.conj()
    assert_allclose(actual, np.r_[expected.real, expected.imag])
    assert_allclose(reverse@reverse, np.eye(4))
    # Reversing time swaps incoming/outgoing: S* S = I for reciprocal lossless S.
    s = solve_port(Bridge(), 35.)['scattering']
    assert_allclose(s.conj()@s, np.eye(2), atol=3e-9)
    real = realify(s)
    assert_allclose(real.T@real, np.eye(4), atol=3e-9)


def test_no_radiative_monopole_or_vacuum_lapse_substitution():
    with pytest.raises(ValueError):
        solve_port(Bridge(), 3., ell=0)
    with pytest.raises(ValueError):
        Bridge(lapse_depth=1.)


def test_static_flux_energy_and_both_shape_forces():
    from dataclasses import replace
    for depth in [0., .6]:
        b = Bridge(lapse_depth=depth)
        result = static_flux(b, charge=.7)
        if depth == 0:
            exact = .7**2/(8*np.pi)*VacuumThroat(.3).resistance()
            assert_allclose(result['energy'], exact, rtol=3e-14)
        for parameter in ['waist', 'stretch']:
            step = 1e-5
            plus = static_flux(replace(b, **{parameter: getattr(b, parameter)+step}), .7)['energy']
            minus = static_flux(replace(b, **{parameter: getattr(b, parameter)-step}), .7)['energy']
            assert_allclose(result[parameter+'_force'], -(plus-minus)/(2*step), rtol=2e-9)
        assert result['waist_force'] > 0
        assert result['stretch_force'] < 0


def test_charged_slice_satisfies_maxwell_hamiltonian_constraint_and_rim_match():
    for fraction in [.3, .8, .98]:
        b = Bridge(charge_fraction=fraction)
        x = np.linspace(-b.limit, b.limit, 101)
        r = b.profile(x)['radius']
        first = b.discriminant*np.sinh(2*x)/r
        second = b.mass/r**2-b.geometric_charge**2/r**3
        scalar_curvature = 2/r**2-4*second/r-2*first**2/r**2
        assert_allclose(scalar_curvature, 2*b.geometric_charge**2/r**4, rtol=2e-12, atol=1e-11)
        assert_allclose(r[-1], np.sin(b.angle), rtol=1e-14)
        assert_allclose(first[-1], np.cos(b.angle), rtol=1e-14)
        stress = b.support(x)
        assert_allclose(stress[0], b.geometric_charge**2/r**4, rtol=2e-12, atol=1e-11)
        # The fixed radial Maxwell field has rho+p_s=0; additional static support violates NEC.
        rho, pr, pt = b.support(0.)
        assert rho+pr < 0
        flux = static_flux(b, np.sqrt(4*np.pi)*b.geometric_charge)
        assert flux['energy'] > 0 and flux['stretch_force'] < 0
