from dataclasses import replace

import numpy as np
import pytest

from geometrodynamics.waves.coherent_mouth_recoil import (
    PathParameters, closed_paths, finite_mass_scatter, focus_overlap,
    junction, path_sum, recoil_operator, scatter, translated_junction,
)


def test_interface_matches_value_and_flux_continuity():
    # Physical (not energy-normalized) travelling amplitudes satisfy
    # 1+r=t_phys and Ys*(1-r)=Yt*t_phys for incidence from the surface.
    for z in np.geomspace(.001, 1000, 37):
        s = junction(z)
        r, t = s[:, 0]
        assert abs(1+r-t/np.sqrt(z)) < 1e-12
        assert abs(1-r-np.sqrt(z)*t) < 1e-12
        np.testing.assert_allclose(s.T@s, np.eye(2), atol=6e-16)


@pytest.mark.parametrize('z', [.001, .03, 3-2*np.sqrt(2), 1., 7.])
def test_no_lone_wave_sign_force_and_analytic_interference(z):
    r, t = junction(z)[0]
    for phi in np.linspace(0, 2*np.pi, 17):
        solo = scatter([np.exp(1j*phi), 0], z)
        assert solo['impulse'] == pytest.approx(2*r*r)
        for ratio in [.01, .3, 1, 4]:
            a = np.array([1, np.sqrt(ratio)*np.exp(1j*phi)])
            out = scatter(a, z)
            expected = 2*r*r*(1-ratio)+4*r*t*np.sqrt(ratio)*np.cos(phi)
            assert out['impulse'] == pytest.approx(expected, abs=4e-15)
            assert abs(out['energy_residual']) < 3e-15
            assert scatter(-a, z)['impulse'] == pytest.approx(out['impulse'])


def test_phase_average_and_balanced_optimum():
    phases = 2*np.pi*np.arange(128)/128
    for z in [.001, .05, 3-2*np.sqrt(2), .8, 1.]:
        values = [scatter([1, np.exp(1j*p)], z)['impulse']/2 for p in phases]
        assert abs(np.mean(values)) < 1e-15
        assert max(abs(np.asarray(values))) <= 1+1e-14
    z = (np.sqrt(2)-1)**2
    assert scatter([1, -1], z)['impulse']/2 == pytest.approx(-1)
    assert scatter([1, 1], z)['impulse']/2 == pytest.approx(1)
    assert np.linalg.eigvalsh(recoil_operator(z))[0] == pytest.approx(-np.sqrt(2))


def test_recoil_operator_bound_for_arbitrary_coherent_inputs():
    rng = np.random.default_rng(8341)
    for z in np.geomspace(.001, 100, 23):
        r = junction(z)[0, 0]
        for _ in range(20):
            a = rng.normal(size=2)+1j*rng.normal(size=2)
            result = scatter(a, z)
            assert result['impulse'] == pytest.approx(np.vdot(a, recoil_operator(z)@a).real)
            assert abs(result['impulse']) <= 2*abs(r)*result['energy']+1e-12


@pytest.mark.parametrize('velocity', [-.7, 0., .7])
def test_finite_mass_energy_momentum_and_mass_shell(velocity):
    rng = np.random.default_rng(307)
    mass = 1.3; momentum = mass*velocity/np.sqrt(1-velocity**2)
    for energy in [1e-5, .03, 1., 4.]:
        for z in [.01, .171572875, 1., 5.]:
            a = rng.normal(size=2)+1j*rng.normal(size=2)
            a *= np.sqrt(energy/np.vdot(a, a).real)
            result = finite_mass_scatter(a, z, mass, momentum)
            assert abs(result['energy_residual']) < 1e-13
            assert abs(result['momentum_residual']) < 1e-13
            assert result['mouth_energy']**2-result['momentum']**2 == pytest.approx(mass**2)
            assert result['rest_energy_scale'] > 0
            if z == 1:
                assert result['impulse'] == pytest.approx(0, abs=1e-13)


def test_finite_mass_reflection_has_recoil_redshift():
    energy, mass = .3, 2.
    result = finite_mass_scatter([np.sqrt(energy), 0], 1e-18, mass)
    # Independent elastic photon-mirror formula in c=1 units.
    expected = energy*mass/(mass+2*energy)
    assert abs(result['outgoing'][0])**2 == pytest.approx(expected)
    assert result['impulse'] == pytest.approx(energy+expected)


def test_closed_paths_preserve_seed_energy_and_do_not_create_a_source():
    for handle in ['klein', 'torus']:
        for coherent in [False, True]:
            p = PathParameters(surface_bins=8, throat_bins=3, periods=12, handle=handle)
            result = closed_paths(p, coherent=coherent, sample=True)
            assert result['max_energy_error'] < 3e-14
            assert result['final_energy'] == pytest.approx(1.)
            assert all(abs(row[-1]-1) < 3e-14 for row in result['sample'])
            times = [row[0] for row in result['sample']]
            assert times[1] == pytest.approx(3*np.pi/8)  # first throat exit


def test_incoherent_network_cannot_see_compression_sign():
    p = PathParameters(surface_bins=8, throat_bins=3, periods=12)
    a = closed_paths(p, coherent=False)
    b = closed_paths(replace(p, handle='torus', seed_phase=.73), coherent=False)
    np.testing.assert_allclose(a['impulse'], b['impulse'], atol=2e-14)


def test_resolvent_matches_explicit_all_walks_and_tracks_escape():
    for handle in ['klein', 'torus']:
        result = path_sum(1.13, .19, .565, survival=.8, handle=handle)
        assert result['walk_sum_relative_error'] < 1e-14
        assert result['unitary_error'] < 1e-15
        assert result['escape_power'] == pytest.approx((1-.8**2)*result['exposure'])
        assert abs(result['impulse_per_exposure']) <= 2*abs(junction(.19)[0, 0])+1e-14
    with pytest.raises(ValueError):
        path_sum(1, .1, .5, survival=1.)


def test_aperture_overlap_scales_energy_not_phase():
    assert focus_overlap(0) == 0
    assert focus_overlap(1) == pytest.approx(1-np.exp(-.5))
    assert focus_overlap(20) == pytest.approx(1)
    assert focus_overlap(.001) == pytest.approx(.001**2/2, rel=1e-6)


def test_stress_equals_translation_virtual_work():
    for z in [.01, .171572875, 1., 3.]:
        w, h = 1.37, 1e-6
        derivative = (translated_junction(h,w,z)-translated_junction(-h,w,z))/(2*h)
        force = -1j/w*junction(z).T@derivative
        np.testing.assert_allclose(force,recoil_operator(z),atol=3e-12)


def test_advanced_clock_delays_match_history_map_without_energy_gain():
    from geometrodynamics.waves.antipodal_throat_transport import inertial_mouth
    length, offset = .565, 2.
    a = inertial_mouth(0., horizon=20., anchor=0.)
    b = inertial_mouth(np.pi, horizon=20., anchor=offset)
    result = path_sum(1.137,.17,length,clock_offset=offset,survival=.8)
    assert result['throat_delays'][0] == pytest.approx(b.external(a.sigma(5.)+length)-5.)
    assert result['throat_delays'][1] == pytest.approx(a.external(b.sigma(5.)+length)-5.)
    assert result['advanced_boundary']
    assert result['unitary_error'] < 1e-15
    assert result['walk_sum_relative_error'] < 1e-14
