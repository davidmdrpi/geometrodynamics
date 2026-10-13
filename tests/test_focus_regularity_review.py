import math

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from geometrodynamics.bulk.focus_regularity_review import SmoothRegularFocus, vacuum_neck_control
from geometrodynamics.bulk.gauss_bonnet import _riemann


def test_vacuum_geodesic_and_domain_condition():
    d = vacuum_neck_control()
    b = d['neck_radius']
    # Independent proper-time quadrature; then integrate the geodesic
    # acceleration r''=-b^2/r^3 from rest, before its singular endpoint.
    integral = quad(lambda r: r/math.sqrt(b*b-r*r), 0, b)[0]
    assert integral == pytest.approx(d['waist_proper_time_to_singularity'], rel=1e-11)
    t = np.linspace(0, .95*b, 100)
    sol = solve_ivp(lambda _, y: [y[1], -b*b/y[0]**3], (0, t[-1]), [b, 0.],
                    t_eval=t, rtol=2e-12, atol=1e-14)
    assert sol.success
    np.testing.assert_allclose(sol.y[0], np.sqrt(b*b-t*t), rtol=2e-10)
    r_m = math.hypot(b, d['half_length'])
    x_m = math.exp(r_m/b)*math.sqrt((r_m-b)/(r_m+b))
    assert math.exp(d['log_kruskal_mouth_x']) == pytest.approx(x_m)
    assert d['exterior_causally_excluded'] and x_m > 1
    assert not vacuum_neck_control(mouth_angle=1.3)['exterior_causally_excluded']


def test_vacuum_curvature_coefficient_independently():
    b = vacuum_neck_control()['neck_radius']
    def metric(p):
        _, r, chi, theta, _ = p
        F = 1-b*b/r**2
        return np.diag([-F, 1/F, r*r, r*r*np.sin(chi)**2,
                        r*r*np.sin(chi)**2*np.sin(theta)**2])
    for radius in [.7*b, 1.5*b]:  # interior and exterior; avoid horizon chart
        point = np.array([0., radius, 1.1, 1.2, .4])
        g = metric(point)
        inv = np.linalg.inv(g)
        rup = _riemann(metric, point, 1e-5)
        ricci = np.einsum('mamb->ab', rup)
        assert np.max(abs(np.diag(ricci)/np.diag(g))) < .001
        lower = np.einsum('am,mbcd->abcd', g, rup)
        raised = np.einsum('abcd,ai,bj,ck,dl->ijkl', lower, inv, inv, inv, inv)
        assert np.einsum('abcd,abcd->', lower, raised) == pytest.approx(72*b**4/radius**8, rel=1e-5)


@pytest.mark.parametrize('epsilon', [0., .1])
def test_smooth_collar_stress_against_coordinate_curvature(epsilon):
    m = SmoothRegularFocus(epsilon=epsilon)
    S = m.geometry.half_length
    for fraction in [.4, .94, .97, 1.01]:
        s = fraction*S
        p = np.array([0., s, 1.1, 1.2, .4])
        g = m.metric(p)
        inv = np.linalg.inv(g)
        d = m.invariants(s)
        target = np.array([d['kappa_rho'], d['kappa_p_radial']]+[d['kappa_p_tangential']]*3)
        errors = []
        for h in [1e-4, 5e-5]:
            rup = _riemann(m.metric, p, h)
            ricci = np.einsum('mamb->ab', rup)
            scalar = np.einsum('ab,ab->', inv, ricci)
            stress = np.diag(ricci-.5*g*scalar)/abs(np.diag(g))
            errors.append(np.max(abs(stress-target)))
        assert errors[1] < 1e-3
        assert errors[1] < .4*errors[0]
        h = 1e-6*S
        dp = (m.invariants(s+h)['kappa_p_radial']-m.invariants(s-h)['kappa_p_radial'])/(2*h)
        r, rp, _, _, lp, _ = m.fields(s)
        assert abs(dp+lp*(d['kappa_rho']+d['kappa_p_radial'])+
                   3*rp/r*(d['kappa_p_radial']-d['kappa_p_tangential'])) < 1e-5


@pytest.mark.parametrize('collar', [.1, .2])
def test_smooth_matching_and_distinct_conservative_transport(collar):
    a = SmoothRegularFocus(collar_fraction=collar)
    b = SmoothRegularFocus(epsilon=.1, collar_fraction=collar)
    S = a.geometry.half_length
    for sign in [-1, 1]:
        expected = [math.sin(.5), sign*math.cos(.5), -math.sin(.5), 1., 0., 0.]
        np.testing.assert_allclose(a.fields(sign*S), expected, atol=1e-13)
        np.testing.assert_allclose(b.fields(sign*S), expected, atol=1e-13)
    grid = np.linspace(-S, S, 1001)
    da, db = a.invariants(grid), b.invariants(grid)
    np.testing.assert_array_equal(da['kappa_rho'], db['kappa_rho'])
    assert np.max(abs(da['kappa_rho'])) > 1.  # disclosed, nonvacuum collar
    assert np.all(db['radius'] >= b.geometry.neck_radius)
    assert np.all(np.isfinite(db['kretschmann']))
    assert np.all(db['theta_plus']*db['theta_minus'] <= 0)
    assert b.invariants(0)['kappa_radial_nec'] == pytest.approx(-3/b.geometry.neck_radius**2)
    assert b.transit_time() < .98*a.transit_time()
    omega = 2/S
    A, B = a.transfer(omega), b.transfer(omega)
    assert np.max(abs(A-B)) > .01
    np.testing.assert_allclose(B, b.transfer(omega, rtol=2e-13), atol=1e-9, rtol=1e-9)
    assert abs(np.linalg.det(B)-1) < 1e-9
    state = np.array([1+.3j, -.2+.7j])
    out = B@state
    assert np.imag(np.conj(out[0])*out[1]) == pytest.approx(np.imag(np.conj(state[0])*state[1]), abs=1e-9)


@pytest.mark.parametrize('fraction', [0, -.1, .3, float('nan')])
def test_invalid_collar(fraction):
    with pytest.raises(ValueError):
        SmoothRegularFocus(collar_fraction=fraction)


def test_high_gain_collar_transfer_meets_absolute_precision_gate():
    m = SmoothRegularFocus(epsilon=.1, collar_fraction=.2)
    omega = .5/m.geometry.half_length
    a = m.transfer(omega, 1, rtol=2e-13)
    b = m.transfer(omega, 1, rtol=5e-14)
    assert np.max(abs(a-b)) < 1e-8
    assert abs(np.linalg.det(a)-1) < 1e-8
