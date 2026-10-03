import copy
import json
import shutil
import numpy as np
import pytest
from scipy.linalg import block_diag, expm
from geometrodynamics.waves import r3_normal_quotient as nq
from geometrodynamics.waves.r3_extension import E
from experiments.closure_ledger import r3_normal_quotient_probe as probe


def rotation(theta=.8):
    return np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])


def test_jordan_quotient_retains_conjugate_neutral_directions():
    J = np.array([[1., 2.], [0., 1.]])
    M = block_diag(J, J, J, J, rotation())
    generators = np.eye(10)[:, [0, 2, 4, 6]]
    B, Q, R, leak = nq.quotient(M, generators)
    out = nq.modal_check(B)
    assert leak < 1e-14 and out['modal_ok']
    assert out['neutral_dimension'] == 4
    assert max(out['power_norms']) < 1.00000001
    assert np.linalg.norm(np.linalg.matrix_power(M, 100), 2) > 100


def test_unremoved_jordan_and_hyperbolicity_fail():
    J = np.array([[1., 2.], [0., 1.]])
    assert not nq.modal_check(block_diag(J, np.eye(2), rotation()))['modal_ok']
    assert not nq.modal_check(block_diag(np.eye(4), np.diag([1.2, 1/1.2])))['modal_ok']


def test_noninvariant_removed_direction_is_detected():
    M = np.eye(3); M[1, 0] = .2
    assert nq.quotient(M, np.eye(3)[:, :1])[3] == pytest.approx(.2)


def test_orthogonal_coordinates_preserve_power_envelope():
    C = block_diag(np.eye(4), rotation())
    O = np.linalg.qr(np.random.default_rng(320).normal(size=(6, 6)))[0]
    a, b = nq.modal_check(C), nq.modal_check(O.T@C@O)
    assert b['modal_ok']
    assert np.allclose(a['power_norms'], b['power_norms'], atol=1e-10)


def test_analytic_rotation_generators():
    z = np.r_[1., .01, np.arange(10)*.03]
    analytic = nq.rotation_generators(z)
    for k, (i, j) in enumerate(((0, 1), (0, 2), (1, 2))):
        O = np.zeros((3, 3)); O[i, j] = 1; O[j, i] = -1
        values = []
        for sign in (-1, 1):
            R = expm(sign*1e-6*O); zz = z.copy()
            for sl in (slice(2, None, 2), slice(3, None, 2)):
                tensor = np.einsum('k,kij->ij', z[sl], E)
                zz[sl] = np.einsum('kij,ij->k', E, R@tensor@R.T)
            values.append(zz)
        assert np.max(abs((values[1]-values[0])/2e-6-analytic[:, k])) < 1e-9

from experiments.closure_ledger import r3_normal_quotient_replay as replay


def evidence():
    return tuple(json.loads((probe.RUN/name).read_text()) for name in ('raw.json', 'result.json'))


def test_authenticated_replay():
    out = replay.replay()
    assert out['verdict'] == 'BOUNDED_CENTER_QUOTIENT_NUMERICALLY'
    assert out['unreduced'] == 'UNREDUCED_HYPERBOLIC_PAIR_PRESENT'
    assert all(x['modal']['neutral_dimension'] == 4 for x in out['reduced'])


@pytest.mark.parametrize('damage', ['matrix', 'sample', 'direction', 'nan', 'frame', 'schedule'])
def test_invalid_evidence_rejected(damage):
    raw, saved = evidence()
    if damage == 'matrix': raw['samples'][0]['M2'][0][0] += .01
    if damage == 'sample': raw['samples'].pop()
    if damage == 'direction': raw['perturbations'].pop()
    if damage == 'nan': raw['samples'][0]['M2'][1][1] = float('nan')
    if damage == 'frame': saved['reduced'][0]['U'][0][0] += .1
    if damage == 'schedule': raw['perturbations'].reverse()
    with pytest.raises(ValueError): replay.score(raw, saved)


def test_changed_endpoint_fails_physical_gate():
    raw, saved = evidence()
    raw['perturbations'][0]['plus']['final'][1] += .001
    out = replay.score(raw, saved)
    assert not out['checks']['Q6']
    assert out['verdict'] == 'NORMAL_RESPONSE_UNRESOLVED'


@pytest.mark.parametrize('filename', ['raw.json', 'result.json', 'manifest.json'])
def test_changed_archive_or_label_rejected(tmp_path, filename):
    for name in ('raw.json', 'result.json', 'manifest.json'):
        shutil.copy(probe.RUN/name, tmp_path/name)
    p = tmp_path/filename
    if filename == 'result.json':
        r=json.loads(p.read_text());r['verdict']='UNSUPPORTED_SUCCESS';p.write_text(json.dumps(r))
    else: p.write_text(p.read_text()+' ')
    with pytest.raises(ValueError, match='fingerprint'): replay.replay(tmp_path)


def test_equivalent_recorded_quotient_frame_replays():
    raw, saved = evidence(); b=saved['reduced'][0]
    O=np.linalg.qr(np.random.default_rng(320).normal(size=(8,8)))[0]
    b['Q']=(np.array(b['Q'])@O).tolist()
    b['U']=(O.T@np.array(b['U'])).tolist()
    out=replay.score(raw,saved)
    assert all(out['checks'].values())
