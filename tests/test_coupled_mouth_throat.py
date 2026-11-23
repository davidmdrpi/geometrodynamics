"""Independent action, conservation and event checks for the scalar candidate."""
from dataclasses import replace
import numpy as np
import pytest
from geometrodynamics.waves.coupled_mouth_throat import CoupledMouthThroat, Parameters, real_basis


def model(**changes):
    return CoupledMouthThroat(replace(Parameters(lmax=3, nodes=7), **changes))


def test_harmonic_derivatives():
    theta = np.array([.4, 1.2, 2.7]); phi = np.array([.2, -.3, 1.1]); h = 1e-6
    _, dt, dp = real_basis(4, theta, phi, True)
    assert np.max(abs(dt-(real_basis(4, theta+h, phi)-real_basis(4, theta-h, phi))/(2*h))) < 1e-8
    assert np.max(abs(dp-(real_basis(4, theta, phi+h)-real_basis(4, theta, phi-h))/(2*h))) < 1e-8


@pytest.mark.parametrize('offset', [0., -1.2])
def test_rhs_is_hamiltonian_gradient(offset):
    m = model(offset=offset); y = m.initial(); rng = np.random.default_rng(29)
    y[:2*m.n] += rng.normal(0, .001, 2*m.n)
    y[2*m.n:2*m.n+3] += [.012, -.017, -.003]
    y[2*m.n+3:2*m.n+6] = [.01, -.02, .004]
    h = 1e-7; grad = []
    for i in range(len(y)-1):
        yp, ym = y.copy(), y.copy(); yp[i] += h; ym[i] -= h
        grad.append((sum(m.evaluate(yp, True)[1].values())-sum(m.evaluate(ym, True)[1].values()))/(2*h))
    grad = np.array(grad)
    expected = np.r_[grad[m.n:2*m.n], -grad[:m.n], grad[2*m.n+3:], -grad[2*m.n:2*m.n+3]]
    assert np.max(abs(m.evaluate(y)[:-1]-expected)) < 2e-8


def test_no_wave_is_exact_rest():
    m = model(energy=0.); y = m.initial()
    assert np.array_equal(m.evaluate(y), np.zeros_like(y))


def test_energy_and_zero_offset_noether():
    m = model(offset=0.); _, ys, _ = m.evolve(horizon=1.)
    es = [sum(m.evaluate(y, True)[1].values()) for y in ys]
    js = [m.angular(y) for y in ys]
    assert np.ptp(es) < 1e-8
    assert np.ptp(js) < 1e-8
    assert np.linalg.norm(ys[-1][2*m.n:2*m.n+2]-ys[0][2*m.n:2*m.n+2]) > 1e-6


def test_rejects_initial_clock_margin():
    with pytest.raises(ValueError, match='clock tilt'):
        model(offset=-2.)


def test_fixed_constraints_conserve_energy():
    m = model(fixed_mouths=True, fixed_radius=True); _, ys, _ = m.evolve(horizon=.6)
    assert np.max(abs(ys[:, 2*m.n:2*m.n+3]-ys[0, 2*m.n:2*m.n+3])) == 0
    es = [sum(m.evaluate(y, True)[1].values()) for y in ys]
    assert np.ptp(es) < 1e-8


def test_clock_reaction_completes_angular_ledger():
    for fixed in [False, True]:
        m = model(fixed_mouths=fixed); _, ys, _ = m.evolve(horizon=1.)
        j = np.array([m.angular(y) for y in ys])
        assert np.ptp(j+ys[:, -1]) < 1e-8
        assert np.ptp(j) > 1e-5


def test_cutoff_stops_without_inventing_release():
    m = model(energy=0.); y = m.initial(); y[2*m.n+5] = -.12
    _, ys, stop = m.evolve(horizon=2., initial=y)
    assert stop == 'radius_cutoff'
    assert abs(ys[-1, 2*m.n+2]-.25*m.par.radius0) < 1e-12
    assert np.max(abs(ys[:, :2*m.n])) == 0
    es = [sum(m.evaluate(y, True)[1].values()) for y in ys]
    assert np.max(abs(np.array(es)-.036)) < 1e-8


def test_uncoupled_full_return_sign_without_relaunch():
    m = model(offset=0., coupling=0.); _, ys, stop = m.evolve(horizon=2*np.pi)
    assert stop == 'horizon'
    assert np.max(abs(ys[-1, :2*m.n]+ys[0, :2*m.n])) < 1e-7
