"""Review controls preserve the unresolved archive and reject invalid evidence."""
import copy
import hashlib
import json
import subprocess
import sys

import numpy as np
import pytest

from experiments.closure_ledger import r3_review_controls as review
from experiments.closure_ledger.r3_prefreeze import r2_degree_scratch as scratch
from geometrodynamics.waves import r3_resonance as r3
from geometrodynamics.waves import nonlinear_supported_tt as dynamics


def test_historical_replay_preserves_verdict_and_does_not_certify_shadowing():
    out = review.replay()
    assert out['registered_verdict'] == 'UNRESOLVED'
    assert out['diagnostic_replay'] == 'VERIFIED'
    assert out['trajectory_shadowing'] == 'NOT_ESTABLISHED'
    assert out['largest_later_A_kick'] > 1e-7


@pytest.mark.parametrize('damage', ['missing', 'duplicate', 'short', 'nan', 'negative', 'windows'])
def test_incomplete_or_nonfinite_diagnostics_are_rejected(damage):
    raw = copy.deepcopy(json.loads(review.ARCHIVE.read_text())['raw'])
    if damage == 'missing':
        raw['runs'].pop()
    elif damage == 'duplicate':
        raw['runs'][-1] = copy.deepcopy(raw['runs'][0])
    elif damage == 'short':
        raw['runs'][0]['increments'].pop()
    elif damage == 'nan':
        raw['runs'][0]['residuals'][0] = float('nan')
    elif damage == 'negative':
        raw['runs'][0]['A_dev'][0] = -1.
    else:
        raw['runs'][0]['windows'].pop()
    with pytest.raises(ValueError):
        review.validate_raw(raw)


def test_jointly_altered_diagnostics_and_result_cannot_authenticate(tmp_path):
    rec = json.loads(review.ARCHIVE.read_text())
    rec['raw']['runs'][0]['increments'][0] += .01
    rec['result'] = review.probe.score(rec['raw'])
    path = tmp_path/'altered.json'
    path.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match='fingerprint'):
        review.replay(path)


def test_replay_tolerates_reversed_summation_order(monkeypatch):
    def reversed_sum(inc):
        inc = np.asarray(inc, float)
        t = (np.arange(len(inc))+.5)/len(inc)
        w = np.exp(-1/(t*(1-t)))
        return float(sum(w[::-1]*inc[::-1])/(2*np.pi*w.sum()))
    monkeypatch.setattr(r3, 'birkhoff', reversed_sum)
    assert review.replay()['registered_verdict'] == 'UNRESOLVED'


@pytest.mark.parametrize('damage', ['verdict', 'resolved', 'label', 'number', 'nan'])
def test_tolerant_replay_still_rejects_changed_decisions_and_material_errors(monkeypatch, damage):
    score = review.probe.score

    def altered(raw):
        fresh = score(raw)
        if damage == 'verdict':
            fresh['verdict'] = 'PASS'
        elif damage == 'resolved':
            fresh['rows'][0]['resolved'] = not fresh['rows'][0]['resolved']
        elif damage == 'label':
            fresh['failures'][0] = 'different failure'
        elif damage == 'number':
            fresh['rows'][0]['c'] += 1e-5
        else:
            fresh['rows'][0]['c'] = float('nan')
        return fresh

    monkeypatch.setattr(review.probe, 'score', altered)
    with pytest.raises(ValueError):
        review.replay()


def test_exact_static_event_time_and_hopf_countercontrol():
    eps = np.array([.001, .01, .05])
    lam = 1.3
    tau = review.exact_event_time(eps, lam)
    np.testing.assert_allclose(-np.sqrt(3)/2*np.sin(2*tau)+eps*lam, 0., atol=1e-16)
    assert np.ptp(tau/eps) > 1e-4  # the historical equality was imposed by linearisation
    for value in (np.nan, np.sqrt(3)/2, 1.):
        with pytest.raises(ValueError):
            review.exact_event_time(value, 1.)
    X, h = scratch.grid(12)
    F = .3*X+.02*scratch.hopf(X)
    np.testing.assert_allclose(np.linalg.norm(F, axis=-1), np.hypot(.3, .02), atol=1e-15)
    with pytest.raises(ValueError):
        scratch.degree(np.zeros_like(X), h)


def test_constraints_do_not_detect_a_state_reset():
    before = r3.section_state(.01, 1)
    after = before.copy()
    after[0] += 1e-5
    after = r3.resolve_clock_velocity(after)
    assert abs(dynamics.constraints(before)['residual'][0]) < 1e-13
    assert abs(dynamics.constraints(after)['residual'][0]) < 1e-13
    assert np.linalg.norm(after-before) > 1e-5


def test_linear_reference_exposes_finite_window_bias_and_angle_dependence():
    primary = review.linear_control()
    secondary = review.linear_control('RK45', 512)
    assert abs(primary['rho0']-1.4846664084) < 1e-10
    assert abs(primary['rho']['48']-secondary['rho']['48']) < 1e-10
    assert primary['half_window_difference'] > 3e-4
    assert abs(primary['bias48']) > 2e-5
    assert abs(primary['alternative_angle_rho48']-primary['rho']['48']) > 2e-4
    # Authenticate the original record and its original producer binding,
    # retained unchanged after the replay-only portability fix.
    blob = (review.probe.ROOT/'experiments/closure_ledger/runs/20260929_r3_review/review.json').read_bytes()
    assert hashlib.sha256(blob).hexdigest() == '73185d30bdc008ef8f7e63a0988aa604b8d900fbe326ac759ce8daed4dc2e583'
    rec = json.loads(blob)
    assert rec['source_sha256'] == {'experiments/closure_ledger/r3_review_controls.py':
                                  '9c362d5c5065afbf9b9d0ffabdacf4c32b1c8812ec9a618f7007d8db1e1fa64c'}
    np.testing.assert_allclose(primary['increments'], rec['linear'][0]['increments'], atol=1e-10, rtol=0)


def test_interval_command_runs_without_missing_deg_file():
    result = subprocess.run([sys.executable, '-m',
                             'experiments.closure_ledger.r3_prefreeze.r2_degree_intervals',
                             '--grid', '6', '--seeds', '8'],
                            cwd=review.probe.ROOT, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert 'KINEMATIC_ONLY' in result.stdout
