"""Antipodal focusing on S3 against S2: refocusing, closure, absorption by a finite mouth, router."""
import numpy as np
import pytest
from scipy.optimize import brentq

from geometrodynamics.waves.antipodal_pair_birth import AntipodalPairBirth, closure_band
from geometrodynamics.waves.pair_birth import Parameters
from geometrodynamics.waves.physical_throat import VacuumThroat
from geometrodynamics.waves.s3_focusing import (absorbed_by_mouth, closure_correction, midpoint_grid,
                                                neck_length_s3, s3_birth_phase, s3_zonal_evolve, zonal_energy,
                                                zonal_fd, zonal_modes_s2, zonal_modes_s3, zonal_momentum)

COARSE = Parameters(cutoff=12, footprint=.2)


def flash(n, width=.25):
    psi = midpoint_grid(n)
    return np.exp(-(psi/width)**2), np.zeros(n)


def test_free_evolution_for_pi_on_s3_is_minus_the_antipodal_map():
    errors = []
    for n in (200, 400):
        u0, p0 = flash(n)
        u, p = s3_zonal_evolve(u0, p0, np.pi)
        np.testing.assert_allclose(u, -u0[::-1], atol=1e-12)        # mirror image, sign flipped
        np.testing.assert_allclose(p, 0., atol=1e-12)
        uf, pf = zonal_fd(u0, p0, np.pi, sphere=3)                   # independent finite volumes
        s = np.sin(midpoint_grid(n))
        errors.append(np.sqrt(np.sum((uf + u0[::-1])**2*s**2)/np.sum(u0**2*s**2)))
        assert abs(zonal_energy(uf, pf) - zonal_energy(u0, p0)) < 1e-5*zonal_energy(u0, p0)
    assert errors[1] < errors[0]/8 and errors[1] < 1e-4


def test_on_s2_the_focus_turns_displacement_into_velocity():
    norms = []
    for n in (200, 400):
        u0, p0 = flash(n)
        u, p = zonal_fd(u0, p0, np.pi, sphere=2)
        s = np.sin(midpoint_grid(n))
        norms.append(np.sqrt(np.sum(u**2*s)/np.sum(u0**2*s)))
        assert np.sqrt(np.sum((u + u0[::-1])**2*s)/np.sum(u0**2*s)) > .9   # not a mirror image
    assert norms[1] < norms[0]/3 and norms[1] < 5e-3                     # displacement -> 0


@pytest.mark.parametrize('k', [0, 1, 2, 3])
def test_a_flash_released_at_rest_closes_exactly_on_s3_and_is_half_off_on_s2(k):
    u0, p0 = flash(400)
    m, q, p = zonal_modes_s3(u0, p0)
    assert closure_correction(m, (-1.)**m, q, p, (2*k + 1)*np.pi) < 1e-25
    l, q2, p2 = zonal_modes_s2(u0, p0, lmax=150)
    assert closure_correction(l + .5, (-1.)**(l + 1), q2, p2, (2*k + 1)*np.pi) == pytest.approx(.5, abs=1e-12)


@pytest.mark.parametrize('k', [0, 1, 2, 3])
def test_on_s3_the_sign_of_the_identification_picks_which_half_of_the_data_closes(k):
    # S = (-1)^m is the sign U(pi) carries for the conformal scalar; for Maxwell U(pi) is the plain
    # pullback (#349). Data at rest close with the first and are wholly off with the opposite sign;
    # data with Q = 0 do the reverse; a uniform 45-degree birth phase is half off either way.
    T = (2*k + 1)*np.pi
    for theta, same, opposite in ((0., 0., 1.), (np.pi/2, 1., 0.), (np.pi/4, .5, .5)):
        m, q, p = zonal_modes_s3(*s3_birth_phase(*flash(400), theta))
        assert closure_correction(m, (-1.)**m, q, p, T) == pytest.approx(same, abs=1e-12)
        assert closure_correction(m, -(-1.)**m, q, p, T) == pytest.approx(opposite, abs=1e-12)


def test_a_flash_that_closes_carries_no_momentum_at_its_foci():
    # Closing data have P = 0 or Q = 0 at every odd multiple of pi, so p du/dpsi vanishes there;
    # mixed data carry momentum density at the focus, and do not close.
    for theta, expected in ((0., 0.), (np.pi/2, 0.), (np.pi/4, .77)):
        u0, p0 = s3_birth_phase(*flash(400), theta)
        energy = zonal_energy(u0, p0)
        for t in (0., np.pi, 3*np.pi):
            assert zonal_momentum(*s3_zonal_evolve(u0, p0, t))/energy == pytest.approx(expected, abs=1e-2 if expected else 1e-13)


@pytest.mark.parametrize('sphere', [3, 2])
@pytest.mark.parametrize('across', [False, True])
def test_a_resting_emitter_gives_an_absorbing_mouth_no_momentum_at_any_offset(sphere, across):
    # Every circle through the emitter that meets the mouth carries two rays with opposite angular
    # momentum, and the mouth takes both.
    for delta in (0., .2, .31, .6, 1.):
        r = absorbed_by_mouth(.3, delta, 0., across=across, sphere=sphere, nodes=120)
        assert abs(r['along']) < 1e-14 and r['across'] < 1e-14 and r['tangent_check'] < 1e-14


@pytest.mark.parametrize('sphere', [3, 2])
def test_a_mouth_holding_the_focus_takes_the_whole_flash_and_absorption_is_continuous(sphere):
    a, u = .3, .27
    for delta in (0., .15, .29):
        along = absorbed_by_mouth(a, delta, u, sphere=sphere, nodes=120)
        across = absorbed_by_mouth(a, delta, u, across=True, sphere=sphere, nodes=120)
        assert along['energy_fraction'] == pytest.approx(1, abs=1e-12)
        assert along['along'] == pytest.approx(1, abs=1e-12)             # the flash's full momentum
        assert across['along'] == pytest.approx(np.cos(delta), abs=1e-12)  # the rest becomes spin
    rest = [absorbed_by_mouth(a, d, 0., sphere=sphere, nodes=200)['energy_fraction'] for d in (a + 1e-6, .45, .9)]
    alpha = np.arcsin(np.sin(a)/np.sin(np.array([a + 1e-6, .45, .9])))
    expected = 1 - np.cos(alpha) if sphere == 3 else 2*alpha/np.pi
    np.testing.assert_allclose(rest, expected, atol=2e-4)              # no jump at the rim


def test_s3_neck_router_closes_in_its_band_and_never_captures_in_quadrature():
    a = .3
    L = neck_length_s3(a)
    assert L == pytest.approx(VacuumThroat(a).length_by_quadrature(), rel=1e-8)
    m = AntipodalPairBirth(COARSE, reflection_phase=-1, neck_length=L)
    low, uk = closure_band(a, 1, length=L)
    assert low < closure_band(a, 1)[0]                                  # a longer neck, a wider band

    def drifts(u):
        v = brentq(lambda v: np.linalg.norm(m.velocity(m.initial(v)), axis=1)[0] - u, .02, .95, xtol=1e-14)
        h, info = m.run(m.initial(v), 3*((np.pi - 2*a)/u + np.pi) - .5, cycles=3, dt=.01)
        return max(abs(e['speed'][0] - u) for e in info['events'] if e['kind'] == 'birth')

    assert drifts(low*(1 + 3e-3)) < 1e-12 and drifts(uk) < 1e-12 and drifts(low*(1 - 3e-3)) > 1e-3
    ms = AntipodalPairBirth(COARSE, reflection_phase=-1, neck_length=L, exchange='signed')
    h, info = ms.run(ms.initial(.3), 40., cycles=4, dt=.01)
    couplings = {e['coupling'] for e in info['events'] if e['kind'] == 'capture'}
    assert couplings and couplings <= {-1., 1.}
    assert info['energy_drift'] < 1e-14 and info['angular_drift'] < 1e-13


@pytest.mark.parametrize('phase', [-1., 1.])
def test_klein_and_torus_agree_on_s3_only_with_minus_one_per_focus(phase):
    # -1 per focus is the field's sign in embedded components; +1 is Maxwell's plain pullback.
    L = neck_length_s3(.3)
    runs = {}
    for handle in ('klein', 'torus'):
        m = AntipodalPairBirth(COARSE, reflection_phase=phase, neck_length=L, exchange='signed', handle=handle)
        h, info = m.run(m.initial(.3), 40., cycles=4, dt=.01)
        exits = sum(abs(e['drate']) for e in info['events'] if e['kind'] == 'throat_exit')
        runs[handle] = h['x'], exits
        assert info['energy_drift'] < 1e-14
    (xk, ek), (xt, et) = runs['klein'], runs['torus']
    if phase < 0:
        assert np.array_equal(xk, xt) and ek == et == 0
    else:
        assert xk.shape != xt.shape or not np.array_equal(xk, xt)
        assert ek > 1e-3 and et > 1e-3
