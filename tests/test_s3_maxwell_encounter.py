import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp

from geometrodynamics.waves.s3_maxwell_encounter import (
    DipolePulse, Worldtube, radial_basis, quadrature, surface_rates, stored, encounter,
    ConductingReceiver,
)


def test_radial_normalization_and_independent_eigen_ode():
    n = np.arange(2, 15)
    chi, w = quadrature(0, np.pi, 160)
    h, d = radial_basis(n, chi)
    f, fp = h*np.sin(chi)**2, d*np.sin(chi)
    assert_allclose((f*w) @ f.T, np.eye(len(n)), atol=2e-14)
    assert_allclose((fp*w) @ fp.T+(h*(2*w*np.sin(chi)**2)) @ h.T,
                    np.diag(n*n), atol=4e-12)
    # Integrate the differential equation, independently of the polynomial recurrence.
    for j in [2, 5, 11]:
        x = np.linspace(.2, 2.8, 91)
        hh, dd = radial_basis([j], x)
        exact = hh[0]*np.sin(x)**2
        sol = solve_ivp(lambda t, y: [y[1], (2/np.sin(t)**2-j*j)*y[0]],
                        [x[0], x[-1]], [exact[0], dd[0, 0]*np.sin(x[0])],
                        t_eval=x, rtol=2e-11, atol=2e-12)
        assert_allclose(sol.y[0], exact, atol=3e-10)


def test_spatial_refocus_is_full_field_identity_for_arbitrary_phase():
    rng = np.random.default_rng(34)
    points = rng.normal(size=(121, 4))
    points /= np.linalg.norm(points, axis=1)[:, None]
    for phase in [0, .37, np.pi/2]:
        pulse = DipolePulse(32, .2, phase, duality=.42)
        a, b = pulse.prepare(points), pulse.prepare(-points)
        for time in [.31, 1.2]:
            e, m, _ = pulse.fields(time+np.pi, a)
            ei, mi, _ = pulse.fields(time, b)
            # Antipodal pullback on embedded tangent vectors supplies -I.
            assert_allclose(e, -ei, atol=2e-12)
            assert_allclose(m, -mi, atol=2e-12)


def test_local_maxwell_equations_from_coordinate_differences():
    pulse = DipolePulse(32, .2, phase=.27)
    chi = np.linspace(.3, 2.7, 17)
    theta, time, step = .8, 1.1, 2e-6
    def components(ch, th, t):
        n = np.array([np.sin(th), 0, np.cos(th)])
        x = np.column_stack([np.cos(ch), np.sin(ch)[:, None]*n])
        ec = np.column_stack([-np.sin(ch), np.cos(ch)[:, None]*n])
        et = np.array([0., np.cos(th), 0., -np.sin(th)])
        ef = np.array([0., 0., 1., 0.])
        e, b, _ = pulse.fields(t, pulse.prepare(x))
        return e @ ef, np.sum(b*ec, axis=1), b @ et
    e, br, bt = components(chi, theta, time)
    dt = (np.array(components(chi, theta, time+step))-
          np.array(components(chi, theta, time-step)))/(2*step)
    cp, cm = components(chi+step, theta, time), components(chi-step, theta, time)
    tp, tm = components(chi, theta+step, time), components(chi, theta-step, time)
    curl_b = ((np.sin(chi+step)*cp[2]-np.sin(chi-step)*cm[2])/(2*step)
              -(tp[1]-tm[1])/(2*step))/np.sin(chi)
    curl_e_r = (np.sin(theta+step)*tp[0]-np.sin(theta-step)*tm[0])/(2*step*np.sin(theta)*np.sin(chi))
    curl_e_t = -(np.sin(chi+step)*cp[0]-np.sin(chi-step)*cm[0])/(2*step*np.sin(chi))
    assert_allclose(dt[0], curl_b, atol=3e-8)
    assert_allclose(dt[1], -curl_e_r, atol=3e-8)
    assert_allclose(dt[2], -curl_e_t, atol=3e-8)


@pytest.mark.parametrize('deformation', [0., .2])
def test_deformed_worldtube_local_energy_and_killing_momentum_balance(deformation):
    pulse = DipolePulse(32, .25, .3)
    tube = Worldtube(.3, .16, deformation, 20)
    time, step = np.pi-.22, 2e-5
    endpoints = stored(pulse, tube, [time-step, time+step], radial_order=20)
    derivative = (endpoints[1]-endpoints[0])/(2*step)
    rates = surface_rates(pulse, tube, time)
    assert_allclose(derivative, [-rates['energy_flux'], -rates['momentum_flux']], atol=3e-7)
    assert_allclose(rates['outgoing_characteristic']-rates['incoming_characteristic'],
                    rates['energy_flux'], atol=2e-15)


def test_finite_window_stored_momentum_cancels_apparent_impulse():
    result = encounter(DipolePulse(32, .25), Worldtube(.3, .16, .2, 20),
                       np.pi-.7, np.pi-.12, 64, 20)
    assert abs(result['incident_only_impulse']) > .01
    assert abs(result['transparent_matter_impulse']) < 1e-9
    assert abs(result['energy_balance_residual']) < 1e-9


def test_surface_geometry_and_duality_invariance():
    tube = Worldtube(.3, .12, .2, 18)
    assert_allclose(np.sum(tube.points*tube.normal, axis=1), 0, atol=2e-15)
    assert_allclose(np.sum(tube.normal*tube.normal, axis=1), 1, atol=2e-15)
    a = surface_rates(DipolePulse(), tube, 2.9)
    b = surface_rates(DipolePulse(duality=np.pi/2), tube, 2.9)
    assert_allclose(list(a.values()), list(b.values()), atol=2e-14)
    sphere = Worldtube(.3, angular_order=12)
    assert_allclose(sphere.area.sum(), 4*np.pi*np.sin(.3)**2, atol=2e-15)


def test_round_s3_phase_closure_depends_on_internal_sign():
    for phase in [0, .3, np.pi/2]:
        pulse = DipolePulse(32, .2, phase)
        q, p = pulse.state(np.pi)
        parity = (-1.)**pulse.n
        for eta in [-1, 1]:
            defect = .5*np.sum((pulse.n*(q-eta*parity*pulse.q0))**2+
                                (p+eta*parity*pulse.p0)**2)
            expected = 4*pulse.energy*(np.sin(phase)**2 if eta == 1 else np.cos(phase)**2)
            assert_allclose(defect, expected, atol=3e-14)


def test_conductor_spectrum_against_continuum_boundary_root():
    from scipy.optimize import brentq
    receiver = ConductingReceiver(cells=1024)
    length = np.pi-receiver.radius
    def boundary(w):
        return w*np.cos(w*length)-np.sin(w*length)/np.tan(length)
    for j in range(4):
        approx = receiver.omega[j]
        exact = brentq(boundary, approx-.08, approx+.08)
        assert abs(approx-exact)/exact < 2e-5


def test_conductor_pressure_agrees_with_independent_shape_work():
    from scipy.linalg import eigh_tridiagonal
    radius, step, cells = .3, 1e-4, 2048
    receiver = ConductingReceiver(radius, cells=cells)
    def lowest_frequency(a):
        h = (np.pi-a)/(cells+1)
        chi = h*np.arange(1, cells+1)
        value = eigh_tridiagonal(2/h**2+2/np.sin(chi)**2, np.full(cells-1, -1/h**2),
                                select='i', select_range=(0, 0), eigvals_only=True)[0]
        return np.sqrt(value)
    derivative = (lowest_frequency(radius+step)-lowest_frequency(radius-step))/(2*step)
    # Unit energy single mode: action=1/omega; average stress=-f_b'^2/(2 omega^2).
    load_stress = -receiver.boundary_derivative[0]**2/(2*receiver.omega[0]**2)
    load_shape_work = -derivative/receiver.omega[0]
    assert_allclose(load_stress, load_shape_work, rtol=1e-4)


def test_conductor_discrete_energy_and_negative_radius_load():
    r = ConductingReceiver(cells=256)
    for t in [.0, 2.8, 3.5]:
        q = r.coefficients*np.cos(r.omega*t)
        p = -r.omega*r.coefficients*np.sin(r.omega*t)
        # Reconstruct the energy from the local difference operator, not modal squares.
        y, v = r.modes @ q, r.modes @ p
        energy = .5*(v@v+np.sum(r.diagonal*y*y)+2*np.sum(r.off*y[:-1]*y[1:]))
        assert_allclose(energy, 1, atol=3e-13)
    result = r.encounter(np.pi-1, np.pi+1)
    assert result['radius_impulse'] < -.1
    assert result['vector_impulse'] == 0
