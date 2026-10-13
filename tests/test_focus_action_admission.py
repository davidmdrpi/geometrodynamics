import math
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from geometrodynamics.bulk.focus_action_admission import EllisGhost, scalar_boundary_obstruction
from geometrodynamics.bulk.gauss_bonnet import _riemann


def test_original_neck_scalar_lapse_obstruction_by_integration():
    d = scalar_boundary_obstruction()
    b, S, rm = d['neck_radius'], d['half_length'], d['mouth_radius']
    sol = solve_ivp(lambda s, y: [y[1], -2*s/(s*s+b*b)*y[1]-b*b/(s*s+b*b)**2*y[0]],
                    (0, S), [rm/b, 0.], rtol=1e-12, atol=1e-13)
    assert sol.success
    assert sol.y[0, -1] == pytest.approx(1., abs=2e-12)
    assert sol.y[1, -1] == pytest.approx(d['forced_right_mouth_lapse_derivative'], abs=2e-11)
    assert d['unit_lapse_positive_integral'] == pytest.approx(d['integral_exact'], abs=1e-13)
    assert d['unit_lapse_positive_integral'] > 0
    assert d['forced_right_mouth_lapse_derivative'] < -1


@pytest.mark.parametrize('s', [0., .31, -1.1])
def test_action_stress_and_curvature_independently(s):
    m = EllisGhost()
    point = np.array([0., s, 1.1, 1.2, .4])
    g = m.metric(point)
    inv = np.linalg.inv(g)
    rup = _riemann(m.metric, point, 1e-4)
    ricci = np.einsum('mamb->ab', rup)
    scalar = np.einsum('ab,ab->', inv, ricci)
    stress = np.diag(ricci-.5*g*scalar)/abs(np.diag(g))
    rho, ps, pt = m.stress(s)
    np.testing.assert_allclose(stress, [rho, ps, pt, pt, pt], atol=2e-6)
    lower = np.einsum('am,mbcd->abcd', g, rup)
    raised = np.einsum('abcd,ai,bj,ck,dl->ijkl', lower, inv, inv, inv, inv)
    r, rp, rpp, _, cp = m.fields(s)
    assert np.einsum('abcd,abcd->', lower, raised) == pytest.approx(60/r**12, abs=1e-5)
    assert rp*rp == pytest.approx(1-1/r**4, abs=1e-12)
    h = 1e-5
    cpp = (m.fields(s+h)[4]-m.fields(s-h)[4])/(2*h)
    assert abs(cpp+3*rp/r*cp) < 1e-8  # scalar equation, not only Einstein


def test_master_potential_darboux_identity_away_from_gauge_pole():
    m = EllisGhost()
    for x in [.3, .8, 2.]:
        def log_derivative(s):
            r, rp, rpp, _, _ = m.fields(s)
            return -.5*rp/r-rpp/rp
        h = 1e-5
        S = log_derivative(x)
        Sp = (log_derivative(x+h)-log_derivative(x-h))/(2*h)
        assert S*S-Sp == pytest.approx(m.potential(x), abs=3e-8)
    assert m.potential(0) == -7


def test_radial_instability_convergence_and_variational_bound():
    m = EllisGhost()
    modes = [m.radial_mode(spacing=h) for h in [.04, .02, .01]]
    a, b, c = [v['lambda0'] for v in modes]
    assert abs(c-b) < .3*abs(b-a)
    assert abs((4*c-b)/3-(4*b-a)/3) < 1e-5
    assert (4*c-b)/3 == pytest.approx(-2.98495893027790, abs=1e-5)
    assert all(v['residual_l2'] < 1e-8 and v['lambda0'] < -.01 for v in modes)
    assert abs(m.radial_mode(8)['lambda0']-m.radial_mode(16)['lambda0']) < 1e-5
    assert m.gaussian_rayleigh()['quotient'] < -1


def test_scale_and_scalar_flux_matching():
    a, b = EllisGhost(), EllisGhost(b=3.)
    r, rp, rpp, chi, cp = a.fields(.8)
    R, Rp, Rpp, Chi, Cp = b.fields(2.4)
    np.testing.assert_allclose([R,Rp,Rpp,Chi,Cp], [3*r,rp,rpp/3,chi,cp/3])
    assert r**3*cp == pytest.approx(math.sqrt(6))
    assert R**3*Cp == pytest.approx(9*math.sqrt(6))
    assert a.stress(0)[0] == pytest.approx(-3, abs=1e-14)


@pytest.mark.parametrize('kwargs', [{'b':0}, {'extent':0}, {'b':float('nan')}])
def test_invalid_background(kwargs):
    with pytest.raises(ValueError):
        EllisGhost(**kwargs)
