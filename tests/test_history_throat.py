from dataclasses import replace

import numpy as np
import pytest

from geometrodynamics.waves.history_throat import (
    EndpointHistory, HistoryThroat, Parameters,
)
from geometrodynamics.waves.coupled_mouth_throat import real_basis


def model(**kwargs):
    return HistoryThroat(replace(Parameters(lmax=3, nodes=7, horizon=1.2,
        neck_horizon=1.2, spacing=.02), **kwargs))


def test_fourier_footprint_and_rotation_derivative():
    m = model(); q = np.array([-.4, 0., .7, 3.4])
    f, fq = m.footprint(q)
    expected = real_basis(m.par.lmax, np.arccos(np.cos(q)),
                         np.where(np.sin(q) >= 0, 0., np.pi))*m.smooth
    assert np.max(abs(f-expected)) < 1e-12
    assert np.max(abs(fq-f@m.d)) < 1e-12


@pytest.mark.parametrize('orientation', [-1, 1])
def test_full_canonical_gradient_and_proper_clocks(orientation):
    m = model(time_orientation=orientation)
    t = m.sigma
    neck = EndpointHistory(t, np.c_[.01*np.sin(t), .02*np.cos(t)],
                              np.c_[.01*np.cos(t), -.02*np.sin(t)])
    y = m.initial(); y[2*m.n+4:2*m.n+6] = [.2, .3]
    y[2*m.n+6:] = [.01, -.02]
    rng = np.random.default_rng(17); y[:2*m.n] += rng.normal(0, .002, 2*m.n)
    h = 1e-6; grad = np.empty_like(y)
    for i in range(len(y)):
        plus, minus = y.copy(), y.copy(); plus[i] += h; minus[i] -= h
        grad[i] = (m.hamiltonian(plus, neck)-m.hamiltonian(minus, neck))/(2*h)
    n = m.n
    expected = np.r_[grad[n:2*n], -grad[:n], grad[2*n+2:2*n+4],
                     -grad[2*n:2*n+2], grad[2*n+6:2*n+8], -grad[2*n+4:2*n+6]]
    rhs = m.rhs(0., y, neck)
    assert np.max(abs(rhs-expected)) < 1e-8
    assert np.max(abs(rhs[2*n+4:2*n+6]**2+rhs[2*n:2*n+2]**2-1)) < 1e-14
    assert abs(grad@rhs) < 1e-9
    jgrad = np.zeros_like(y)
    jgrad[:n] = m.d@y[n:2*n]; jgrad[n:2*n] = -m.d@y[:n]
    jgrad[2*n+2:2*n+4] = 1.
    assert abs(jgrad@rhs) < 1e-12


def test_clock_energy_changes_inertia():
    m = model(); neck = EndpointHistory(m.sigma, np.zeros((len(m.sigma), 2)),
                                       np.zeros((len(m.sigma), 2)))
    a = m.initial(); b = a.copy(); b[2*m.n+6] += .1
    assert m.rhs(0., b, neck)[2*m.n] < m.rhs(0., a, neck)[2*m.n]


def test_free_clocks_and_exact_transfer_map():
    m = model(coupling=0., charge=0., packet_energy=0., birth_a=.4, birth_b=-.2)
    r = m.solve(); t = r['t']; clocks = r['states'][:, 2*m.n+4:2*m.n+6]
    assert np.max(abs(clocks-t[:, None]*np.sqrt(1-np.array([.4, .2])**2))) < 1e-12
    event = m.transfer_event(r, .1, 0, 1, .2)
    expected = (.1*np.sqrt(1-.4**2)+.2)/np.sqrt(1-.2**2)
    assert event['arrival'] == pytest.approx(expected, abs=1e-12)
    assert m.transfer_event(r, 1.2, 0, 1, 4.) is None


def test_time_reversal_reverses_event_map_not_field_sign():
    m = model(coupling=0., charge=0., packet_energy=0., birth_a=.1, birth_b=-.1,
              time_orientation=-1, mouth_angle=.2)
    r = m.solve()
    a = m.transfer_event(r, .3, 0, 1, .2)
    b = m.transfer_event(r, .4, 0, 1, .2)
    rate = np.sqrt(1-.1**2)
    assert a['arrival'] == pytest.approx((1.2-.2)/rate-.3)
    assert b['arrival'] == pytest.approx((1.2-.2)/rate-.4)
    assert a['time_map_orientation'] == -1
    assert m.transfer_event(r, .8, 0, 1, .2)['elapsed_exterior'] < 0


def test_history_lookup_rejects_extrapolation():
    h = EndpointHistory([0., 1.], np.zeros((2, 2)), np.zeros((2, 2)))
    with pytest.raises(ValueError, match='outside'):
        h.sample([-.1, .2])


def test_free_surface_returns_sign_flipped_without_relaunch():
    m = model(coupling=0., charge=0., horizon=2*np.pi, neck_horizon=2*np.pi)
    r = m.solve()
    assert r['summary']['stop'] == 'horizon'
    assert np.max(abs(r['states'][-1, :2*m.n]+r['states'][0, :2*m.n])) < 1e-7


@pytest.mark.parametrize('orientation', [-1, 1])
def test_reciprocal_neck_clock_work_and_exterior_ledger(orientation):
    m = model(coupling=.1, time_orientation=orientation, spacing=.01,
              fixed_point_tolerance=2e-9)
    r = m.solve()
    assert r['summary']['converged']
    assert r['summary']['energy_relative_residual'] < 1e-7
    assert r['summary']['angular_absolute_residual'] < 1e-7
    assert abs(r['summary']['neck_clock_work_residual']) < 2e-7
    assert np.max(abs(r['neck_states'])) > 1e-5


def test_failed_iteration_is_not_accepted():
    r = model(max_iterations=1).solve()
    assert not r['summary']['converged']


def test_nonpositive_clock_mass_is_rejected():
    m = model(); y = m.initial(); y[2*m.n+6:] = -2.
    zero = np.zeros((len(m.sigma), 2)); neck = EndpointHistory(m.sigma, zero, zero)
    with pytest.raises(ValueError, match='Nonpositive'):
        m.rhs(0., y, neck)
