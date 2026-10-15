import numpy as np

from geometrodynamics.waves import r3_half_clock as h
from geometrodynamics.waves import r3_return_map as rm
from experiments.closure_ledger import r3_half_clock_probe as probe


def test_clock_involution_covariance():
    rng = np.random.default_rng(21)
    for _ in range(12):
        y = np.array([1., 0., .2, -1., .03, .01])+rng.normal(0, .02, 6)
        np.testing.assert_allclose(rm.rhs(h.clock_sign(y)), h.clock_sign(rm.rhs(y)),
                                   atol=1e-14, rtol=1e-14)
        np.testing.assert_array_equal(h.clock_sign(h.clock_sign(y)), y)


def test_lift_denominators():
    assert str(h.lifted_rotation(2, 5)) == '7/10'
    assert str(h.lifted_rotation(3, 7)) == '5/7'
    assert str(h.lifted_rotation(3, 8)) == '11/16'
    assert str(h.lifted_rotation(1, 2)) == '3/4'


def test_square_against_full_equations_off_resonance():
    z = np.array([1., 0., .01, 0.])
    half = h.half_map(z, full=True)
    assert half['state'][3] < 0 and half['constraint'] < 1e-10
    full = rm.esu_map(z, tol=(1e-13, 1e-15))[0]
    np.testing.assert_allclose(h.half_map(half['z']), full, atol=1e-10, rtol=0)


def records(harmonic=7, amplitude=2e-8):
    points = [dict(j=j, ok=True, lam=amplitude*np.sin(harmonic*(2*np.pi*j/84+.013)),
                   residual=1e-14, constraint=1e-14, cond=100.) for j in range(84)]
    noise = [dict(j=j, ok=True, primary=points[j]['lam'], alternatives={
        method: dict(lam=points[j]['lam'], residual=1e-14)
        for method in ('matrix_DOP853', 'matrix_Radau')}) for j in probe.NOISE_PHASES]
    known = dict(identity=[dict(square_error=1e-12, correct_lift_error=1e-12, wrong_lift_error=.02)])
    return points, noise, known


def test_score_can_fail_for_wrong_harmonic_and_small_signal():
    assert probe.score(*records())['label'] == 'SEVENTH_HARMONIC_PREDICTION_SUPPORTED'
    assert probe.score(*records(harmonic=14))['label'] == 'SEVENTH_HARMONIC_PREDICTION_FAILED'
    assert probe.score(*records(amplitude=1e-13))['label'] == 'SEVENTH_HARMONIC_PREDICTION_FAILED'


def test_missing_noisy_and_incomplete_scans_are_unresolved():
    points, noise, known = records()
    assert probe.score(points[:-1], noise, known)['label'] == 'NUMERICALLY_UNRESOLVED'
    assert probe.score(points, noise[:-1], known)['label'] == 'NUMERICALLY_UNRESOLVED'
    noise[0]['alternatives']['matrix_Radau']['lam'] += 1e-8
    assert probe.score(points, noise, known)['label'] == 'NUMERICALLY_UNRESOLVED'


def test_identity_failure_is_reported():
    points, noise, known = records()
    known['identity'][0]['square_error'] = 1e-4
    assert probe.score(points, noise, known)['label'] == 'HALF_CLOCK_IDENTITY_FAILED'
