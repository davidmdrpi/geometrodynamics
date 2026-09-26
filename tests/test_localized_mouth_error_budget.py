"""Independent geometry calibration, evidence replay and false-pass controls."""
import copy
import gzip
import json
from pathlib import Path
import subprocess
import sys
import mpmath as mp
import pytest
from geometrodynamics.waves import coordinate_budget as c
from experiments.closure_ledger import localized_mouth_error_budget_probe as p

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT/'experiments/closure_ledger/runs/20260916_localized_mouth'
NEW = ROOT/'experiments/closure_ledger/runs/20260922_localized_mouth_error_budget'


@pytest.fixture(scope='module')
def inputs():
    return p.load_inputs(OLD)


@pytest.fixture(scope='module')
def evidence():
    return json.loads((NEW/'error_budget.json').read_text())


def test_coordinate_engine_analytic_geometry_and_nonzero_divergence():
    with mp.workdps(80):
        assert p.calibration_pass(p.calibrate())


def test_exact_round_quartet_solves_jordan_constraints():
    with mp.workdps(80):
        def values(s, t, azimuth):
            q = mp.sqrt(3)/2
            factor = 1/mp.cosh(s)**2
            n = [mp.sin(t)*mp.cos(azimuth), mp.sin(t)*mp.sin(azimuth), mp.cos(t)]
            return ([factor, factor, factor*mp.sin(t)**2]+[mp.mpf(0)]*3
                    +[q*mp.tanh(s)]+[q/mp.cosh(s)*v for v in n]+[mp.mpf(0)]*4)
        r = c.coordinate_check(values, ['.7', '.91', '.37'])
        assert abs(r['R']-6) < mp.mpf('1e-70')
        assert abs(r['spatial']-mp.mpf(9)/4) < mp.mpf('1e-70')
        assert abs(r['H']) < mp.mpf('1e-70')
        assert c.norm(r['M']) < mp.mpf('1e-70')


def test_nonzero_scalar_current_sign_and_geometric_weight():
    with mp.workdps(80):
        def values(s, t, azimuth):
            v = c.calibration_values(tensor=True)(s, t, azimuth)
            v[6], v[10] = s, mp.mpf(1)/5
            return v
        r = c.coordinate_check(values, ['.7', '.91', '.37'])
        expected_H = mp.mpf(7)/8*(2+2*mp.mpf('.7')**2)-mp.mpf(1)/25-1-3
        assert abs(r['H']-expected_H) < mp.mpf('1e-70')
        assert abs(r['M'][0]-(-mp.mpf(7)/4+mp.mpf(1)/5)) < mp.mpf('1e-70')


def test_full_archived_evidence_replays_and_preserves_history(inputs, evidence):
    result = p.score(evidence, inputs)
    assert result['passed'] == 8
    assert all(result['verdicts'].values())
    assert result['historical'] == dict(original='6/8', stable_reconstruction='7/8', unchanged=True)
    assert inputs[-1]['passed'] == 7
    assert not any(inputs[-1]['verdicts'].values())


def test_local_stencils_stay_inside_all_saved_polynomials(inputs, evidence):
    with mp.workdps(80):
        model = c.SavedData(inputs[0], inputs[1])
        for row in evidence['rows']:
            s, h = mp.mpf(row['point'][0]), mp.mpf(row['h0'])
            for poly in (model.psi, model.theta, model.tensor):
                assert all(not s-2*h <= knot <= s+2*h for knot in poly.x)


@pytest.mark.parametrize('damage', ['K', 'Pi', 'H', 'step', 'point', 'missing', 'source', 'nan'])
def test_evidence_rejects_changed_measurements(inputs, evidence, monkeypatch, damage):
    # The separate full replay test recomputes these measurements. Here cache
    # that baseline to exercise corruption rejection without eight full runs.
    lookup = {tuple(row['point']): row for row in evidence['rows']}
    monkeypatch.setattr(p, 'measure', lambda record, profiles, point: lookup[tuple(p.packed(list(map(mp.mpf, point))))])
    bad = copy.deepcopy(evidence)
    row = bad['rows'][0]
    if damage in ('K', 'Pi'):
        row['reference80'][damage][0] = '1'
    elif damage == 'H':
        row['local'][0]['values']['H'] = '1'
    elif damage == 'step':
        row['local'][0]['h'] = '0.1'
    elif damage == 'point':
        row['point'][0] = '0.1'
    elif damage == 'missing':
        bad['rows'].pop()
    elif damage == 'source':
        bad['sources'][p.SOURCES[0]] = '0'*64
    else:
        row['reference80']['R'] = 'nan'
    result = p.score(bad, inputs)
    assert not result['gates']['evidence']
    assert not any(result['verdicts'].values())


def test_convergent_difference_cannot_hide_large_constant_residual(evidence):
    with mp.workdps(80):
        row = copy.deepcopy(evidence['rows'][0])
        old = mp.mpf(row['reference80']['H'])
        shift = mp.mpf(row['reference80']['Hscale'])/100-old
        for key in ('reference60', 'reference80'):
            row[key]['H'] = p.packed(mp.mpf(row[key]['H'])+shift)
        for item in row['local']:
            item['values']['H'] = p.packed(mp.mpf(item['values']['H'])+shift)
        result = p.row_score(row)
        assert result['checks']['convergence']
        assert result['checks']['differentiation_accuracy']
        assert not result['checks']['absolute']


def test_hash_guard_rejects_changed_saved_polynomial(monkeypatch):
    original_reader = p.read_bytes
    def corrupted(path):
        blob = original_reader(path)
        if path.name == 'refinement.json.gz':
            data = json.loads(gzip.decompress(blob))
            data['solutions'][-1]['solution']['c'][0][0][0] += 1
            return gzip.compress(json.dumps(data).encode())
        return blob
    monkeypatch.setattr(p, 'read_bytes', corrupted)
    with pytest.raises(ValueError, match='hash mismatch'):
        p.load_inputs(OLD)


def test_failed_cli_clears_previous_affirmative_verdict(tmp_path):
    (tmp_path/'error_budget_verdict.json').write_text(json.dumps({'verdicts': {k: True for k in p.FALSE}}))
    command = [sys.executable, '-m', 'experiments.closure_ledger.localized_mouth_error_budget_probe',
               '--input-dir', str(tmp_path/'missing'), '--output-dir', str(tmp_path)]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode != 0
    assert not any(json.loads((tmp_path/'error_budget_verdict.json').read_text())['verdicts'].values())
