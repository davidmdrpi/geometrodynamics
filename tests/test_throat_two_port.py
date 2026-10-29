"""The matched bridge as a Maxwell two-port, and the strict router that routes fronts through it."""
import numpy as np
import pytest
from scipy.optimize import brentq

from geometrodynamics.waves.antipodal_pair_birth import AntipodalPairBirth, closure_band, two_port_closure_band
from geometrodynamics.waves.pair_birth import Parameters
from geometrodynamics.waves.physical_throat import VacuumThroat
from geometrodynamics.waves.throat_two_port import ThroatTwoPort, split

COARSE = Parameters(cutoff=12, footprint=.2)


@pytest.mark.parametrize('depth, expected', [(0., (2.585e-7, 5.679e-5, .851877)),
                                             (.4, (3.854e-6, .005030, .999993)),
                                             (.8, (.204059, .999668, .999991))])
def test_two_port_reproduces_the_published_transmission_table(depth, expected):
    # #349 (3f72d9a), a = 0.3, l = 1, its own integrator; this one shares no code with it.
    tp = ThroatTwoPort(.3, depth)
    assert tp.length == pytest.approx(VacuumThroat(.3).length(), rel=1e-14)
    for omega, T in zip((10., 20., 60.), expected):
        assert abs(tp.scattering(omega)[0])**2 == pytest.approx(T, rel=2e-3)


def test_two_port_reproduces_the_published_delay_table():
    # #349 (dc9c689): its width-0.2 TE-dipole S3 pulse, modes n = 2..48 at omega = n, weighted by
    # energy; transmitted and reflected delays weighted by output. The pulse weights are rebuilt
    # here from f_n = sin^2 chi C_(n-2)^(2)(cos chi)/sqrt(pi (n^2 - 1)/8).
    from numpy.polynomial.legendre import leggauss
    from scipy.special import eval_gegenbauer
    x, wq = leggauss(256)
    chi, wq = np.pi*(x + 1)/2, wq*np.pi/2
    n = np.arange(2, 49)
    f = np.array([np.sin(chi)**2*eval_gegenbauer(k - 2, 2, np.cos(chi)) for k in n])/np.sqrt(np.pi*(n*n - 1)/8)[:, None]
    q = f @ (wq*np.sin(chi)**2*np.exp((np.cos(chi) - 1)/.2**2))
    w = (n*q)**2/np.sum((n*q)**2)
    table = {0.: (6.3588e-7, .66218, .60848, .47126), .4: (2.1436e-5, .83802, .86973, .64425),
             .8: (.27958, 1.40704, 2.00190, 1.76187), .9: (.78995, 1.96375, 2.53000, 2.45260),
             .95: (.92556, 2.75065, 3.77024, 4.32810), .99: (.75704, 6.07056, 5.93775, 4.34448)}
    for depth, (T, optical, tau_t, tau_r) in table.items():
        tp = ThroatTwoPort(.3, depth)
        r = tp.response(n.astype(float))
        wt, wr = w*r['T'], w*r['R']
        assert wt.sum() == pytest.approx(T, rel=1e-4)
        assert tp.optical_length == pytest.approx(optical, abs=1e-5)
        assert wt @ r['raw_tau_t']/wt.sum() == pytest.approx(tau_t, abs=2e-5)
        assert wr @ r['raw_tau_r']/wr.sum() == pytest.approx(tau_r, abs=2e-5)


def test_two_port_is_lossless_symmetric_and_converged():
    coarse, fine = ThroatTwoPort(.3, .6, cells=2**12), ThroatTwoPort(.3, .6, cells=2**16)
    for omega in (1.4, 7., 25., 70.):
        t, r = fine.scattering(omega)
        assert abs(t)**2 + abs(r)**2 == pytest.approx(1, abs=1e-12)
        assert abs((r/t).real) < 1e-9*abs(r/t)                    # symmetric barrier: r and t in quadrature
        assert abs(coarse.scattering(omega)[0]/t - 1) < 1e-5
    resp = fine.response([3., 70.])
    np.testing.assert_allclose(resp['raw_tau_r'], resp['raw_tau_t'], rtol=1e-6)   # equal Wigner delays
    assert resp['raw_tau_t'][0] < fine.optical_length            # tunnelling: held to the light cone
    assert resp['tau_t'][0] == fine.optical_length
    assert resp['tau_t'][1] == pytest.approx(fine.optical_length, rel=.05)       # high frequency: optics


def test_split_conserves_energy_and_filters_the_spectrum():
    m = AntipodalPairBirth(COARSE, throat=ThroatTwoPort(.3, .9))
    w = m.initial(.3).fronts[0].spectrum
    assert w.sum() == pytest.approx(1) and np.all(w >= 0)
    parts = split(w, m.response)
    assert parts['T'] + parts['R'] == pytest.approx(1, abs=1e-14)
    assert parts['transmitted'].sum() == pytest.approx(1) and parts['reflected'].sum() == pytest.approx(1)
    mean = m.spectral_omega @ w
    assert m.spectral_omega @ parts['transmitted'] > mean > m.spectral_omega @ parts['reflected']   # high-pass


@pytest.mark.parametrize('exchange', ['rest', 'signed'])
@pytest.mark.parametrize('depth', [0., .95])
def test_two_port_router_conserves_and_reflects_back_out_of_the_catching_mouth(exchange, depth):
    m = AntipodalPairBirth(COARSE, reflection_phase=-1, exchange=exchange, throat=ThroatTwoPort(.3, depth))
    h, info = m.run(m.initial(.3), 40., cycles=4, dt=.01)
    assert info['energy_drift'] < 1e-13 and info['angular_drift'] < 1e-13
    ev = info['events']
    captures = [e for e in ev if e['kind'] == 'capture']
    assert captures and all(e['transmitted']['fraction'] + e['reflected']['fraction'] == pytest.approx(1)
                            for e in captures)
    exits = [e for e in ev if e['kind'] == 'throat_exit']
    assert any(e['route'] == 'reflect' for e in exits)
    for e in exits:
        # Reflected parts leave the mouth that caught them, transmitted ones the other mouth.
        assert (e['mouth'] == e['source']) == (e['route'] == 'reflect')
        if exchange == 'rest':
            assert np.allclose(e['recoil'], 0)


def test_with_unit_lapse_the_handle_no_longer_matters():
    # With +1 per focus and a whole crossing in the neck's length the handles differ; through the
    # two-port at N = 1 almost nothing crosses, and Klein and torus agree.
    rates = {}
    for throat in (None, ThroatTwoPort(.3, 0.)):
        for handle in ('klein', 'torus'):
            m = AntipodalPairBirth(COARSE, reflection_phase=1, exchange='signed', handle=handle, throat=throat,
                                   neck_length=None if throat else VacuumThroat(.3).length())
            h, info = m.run(m.initial(.3), 40., cycles=4, dt=.01)
            rates[throat is None, handle] = sum(e['drate'] for e in info['events']
                                                if e['kind'] in ('capture', 'throat_exit'))
    assert abs(rates[True, 'klein'] - rates[True, 'torus']) > 1e-2
    assert abs(rates[False, 'klein'] - rates[False, 'torus']) < 1e-7


def test_the_closure_band_narrows_to_the_reflection_delay():
    a = .3
    m = AntipodalPairBirth(COARSE, reflection_phase=-1, throat=ThroatTwoPort(a, 0.))
    low, uk = two_port_closure_band(m, 1)
    assert closure_band(a, 1, length=m.length)[0] < low < uk      # narrower than with the neck's length

    def drifts(u):
        v = brentq(lambda v: np.linalg.norm(m.velocity(m.initial(v)), axis=1)[0] - u, .02, .95, xtol=1e-14)
        h, info = m.run(m.initial(v), 3*((np.pi - 2*a)/u + np.pi) - .5, cycles=3, dt=.01)
        return max(abs(e['speed'][0] - u) for e in info['events'] if e['kind'] == 'birth')

    assert drifts(low*(1 + 3e-3)) < 1e-12 and drifts(uk) < 1e-12 and drifts(low*(1 - 3e-3)) > 1e-3


def test_a_reflected_front_retraces_its_path_and_the_mouth_takes_the_difference():
    # Signed exchange: the throat holds (1 - c) J for the reflected part, which leaves the same mouth
    # with that momentum reversed in the mouth's frame, so (unsaturated) the mouth recoils by about
    # twice what the throat held: a mirror.
    m = AntipodalPairBirth(COARSE, reflection_phase=1, exchange='signed', throat=ThroatTwoPort(.3, 0.))
    h, info = m.run(m.initial(.3), 40., cycles=4, dt=.01)
    out = [e for e in info['events'] if e['kind'] == 'throat_exit' and e['route'] == 'reflect'
           and not e['saturated'] and np.linalg.norm(e['J_out']) > 1e-6]
    assert out
    for e in out:
        held = np.array(e['recoil']) + np.array(e['J_out'])
        assert np.linalg.norm(e['J_out']) == pytest.approx(np.linalg.norm(held), rel=1e-12)
        assert np.dot(e['J_out'], held) < 0                  # reversed, up to the mouth's motion meanwhile
