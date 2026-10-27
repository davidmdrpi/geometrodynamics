import math
import numpy as np
import pytest

from geometrodynamics.bulk.focus_regularity import RegularFocus, round_caustic_control
from geometrodynamics.bulk.gauss_bonnet import _riemann


@pytest.mark.parametrize('epsilon', [0., .1])
def test_stress_against_independent_coordinate_riemann(epsilon):
    model = RegularFocus(epsilon=epsilon)
    for fraction in [0., .31, -.67]:
        s = fraction*model.geometry.half_length
        point = np.array([0., s, 1.1, 1.2, .4])
        exact = model.invariants(s)
        errors = []
        for h in [2e-4, 1e-4]:
            metric = model.metric(point)
            inv = np.linalg.inv(metric)
            rup = _riemann(model.metric, point, h)
            ricci = np.einsum('mamb->ab', rup)
            scalar = np.einsum('ab,ab->', inv, ricci)
            einstein = ricci-.5*metric*scalar
            stress = np.diag(einstein)/abs(np.diag(metric))
            target = np.array([exact['kappa_rho'], exact['kappa_p_radial']]+[exact['kappa_p_tangential']]*3)
            errors.append(np.max(abs(stress-target)))
            lower = np.einsum('am,mbcd->abcd', metric, rup)
            raised = np.einsum('abcd,ai,bj,ck,dl->ijkl', lower, inv, inv, inv, inv)
            K = np.einsum('abcd,abcd->', lower, raised)
            assert np.isclose(K, exact['kretschmann'], rtol=2e-5)
        assert errors[1] < 3e-5
        assert errors[1] < .4*errors[0]


def test_equal_mouth_data_and_neck_size_do_not_select_transit():
    first, second = RegularFocus(), RegularFocus(epsilon=.1)
    for s in [-first.geometry.half_length, first.geometry.half_length]:
        np.testing.assert_allclose(first.fields(s), second.fields(s), atol=1e-14)
    assert first.fields(0)[0] == second.fields(0)[0]
    assert math.isclose(first.transit_time(), first.geometry.proper_length, rel_tol=1e-13)
    assert second.transit_time() < .98*first.transit_time()
    assert np.max(abs(first.transfer(2/first.geometry.half_length)-second.transfer(2/second.geometry.half_length))) > .01


@pytest.mark.parametrize('epsilon', [0., .1])
def test_constraint_conservation_nec_and_expansions(epsilon):
    m = RegularFocus(epsilon=epsilon)
    grid = np.linspace(-m.geometry.half_length, m.geometry.half_length, 101)
    d = m.invariants(grid)
    assert np.max(abs(d['kappa_rho'])) < 1e-12  # K_ij=0, spatial R=0
    assert np.all(d['radius'] >= m.geometry.neck_radius)
    assert np.all(d['lapse'] >= 1.)
    assert np.all(d['theta_plus']*d['theta_minus'] <= 0.)
    assert m.invariants(0)['kappa_radial_nec'] == pytest.approx(-3/m.geometry.neck_radius**2)
    for u in [-.6, -.2, .3, .7]:
        s = u*m.geometry.half_length
        h = 1e-5*m.geometry.neck_radius
        pprime = (m.invariants(s+h)['kappa_p_radial']-m.invariants(s-h)['kappa_p_radial'])/(2*h)
        r, rp, _, _, lp, _ = m.fields(s)
        d = m.invariants(s)
        residual = pprime+lp*(d['kappa_rho']+d['kappa_p_radial'])+3*rp/r*(d['kappa_p_radial']-d['kappa_p_tangential'])
        assert abs(residual) < 2e-6


def test_static_monopole_integral_and_flux_preservation():
    from scipy.integrate import quad
    m = RegularFocus(epsilon=.1)
    b, S = m.geometry.neck_radius, m.geometry.half_length
    integral = quad(lambda s: b*b/(m.fields(s)[3]*m.fields(s)[0]**3), -S, S)[0]
    np.testing.assert_allclose(m.transfer(0), [[1, integral], [0, 1]], rtol=1e-10, atol=1e-11)
    state = np.array([1+.3j, -.2+.7j])
    for ell in [0, 1]:
        for wS in [.5, 2., 5.]:
            M = m.transfer(wS/S, ell)
            tighter = m.transfer(wS/S, ell, rtol=2e-13)
            np.testing.assert_allclose(M, tighter, rtol=1e-9, atol=1e-9)
            assert abs(np.linalg.det(M)-1) < 1e-9
            out = M@state
            assert abs(np.imag(np.conj(out[0])*out[1])-np.imag(np.conj(state[0])*state[1])) < 1e-9


def test_scale_covariance_and_regular_round_caustic():
    a, b = RegularFocus(epsilon=.1), RegularFocus(radius=3, epsilon=.1)
    assert b.transit_time() == pytest.approx(3*a.transit_time())
    assert b.invariants(0)['kretschmann'] == pytest.approx(a.invariants(0)['kretschmann']/81)
    np.testing.assert_allclose(a.transfer(2/a.geometry.half_length), b.transfer(2/b.geometry.half_length), atol=1e-10)
    control = round_caustic_control()
    assert abs(control['jacobi_factor']) < 1e-14
    assert control['kretschmann'] == 24.


@pytest.mark.parametrize('kwargs', [{'epsilon':-.1}, {'epsilon':float('nan')}, {'radius':0}, {'mouth_angle':2}])
def test_invalid_parameters_rejected(kwargs):
    with pytest.raises(ValueError):
        RegularFocus(**kwargs)
