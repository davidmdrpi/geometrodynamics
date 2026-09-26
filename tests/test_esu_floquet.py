import json
from pathlib import Path
import numpy as np
import pytest
from geometrodynamics.waves import esu_floquet as fl
from experiments.closure_ledger import esu_floquet_probe as probe

ARCHIVE = Path(__file__).resolve().parents[1]/'experiments/closure_ledger/runs/20260926_esu_floquet/esu_floquet.json'


def test_free_conformal_field_refocuses_exactly():
    for n in (2, 3, 7, 30):
        o = fl.observables('V', n, fl.monodromy('C', n))
        assert o['defect'] < 1e-10 and abs(o['fidelity']-1) < 1e-10


def test_uncoupled_vector_is_free_field():
    for n in (2, 5):
        M = fl.monodromy('V', n, coupled=False)
        assert np.allclose(M, fl.monodromy('C', n), atol=1e-11)


def test_tensor_n2_reproduces_prior_published_map():
    tr = np.trace(fl.monodromy('T', 2, t1=np.pi/2))
    assert abs(tr-probe.PRIOR_T2_HALF_TRACE) < 1e-8


def test_bare_tensor_phase_is_analytic():
    for n in (2, 9):
        o = fl.observables('T', n, fl.monodromy('B', n))
        a = abs(np.pi*(np.sqrt(n*(n+2))-(n+1))) % (2*np.pi)
        assert abs(o['theta']-min(a, 2*np.pi-a)) < 1e-8


def test_parity_conventions():
    assert fl.parity('S', 2) == -1 and fl.parity('S', 3) == 1
    assert fl.parity('V', 2) == 1 and fl.parity('T', 3) == -1


def test_scalar_unused_trace_equation_and_vector_ij_hold():
    rng = np.random.default_rng(1)
    for n in (2, 6):
        for t in rng.uniform(0, np.pi, 5):
            assert fl.scalar_trace_residual(t, rng.uniform(-1, 1, 4), n) < 1e-10
            assert fl.vector_constraint_residual(t, rng.uniform(-1, 1, 2), n) < 1e-10


def test_wkb_average_is_analytic():
    w = fl.wkb_masses()
    assert abs(w['mean_2R2_over_f']-w['analytic_mean_2R2_over_f']) < 1e-9


def test_classifier_recognises_exact_and_absent():
    exact = [dict(n=n, defect=0., fidelity=1.) for n in range(2, 81)]
    assert probe.classify_refocusing(exact)[0] == 'EXACT'
    absent = [dict(n=n, defect=1., fidelity=.2*(-1)**n) for n in range(2, 81)]
    assert probe.classify_refocusing(absent)[0] == 'ABSENT'


@pytest.mark.skipif(not ARCHIVE.exists(), reason='archive not generated')
def test_archive_replays_and_rejects_tampering():
    data = json.loads(ARCHIVE.read_text())
    assert probe.replay(data, degrees=[2, 17, 80])
    bad = json.loads(ARCHIVE.read_text())
    bad['sectors']['S']['rows'][5]['matrix'][0][0] += 1e-6
    assert not probe.replay(bad, degrees=[])
    bad = json.loads(ARCHIVE.read_text())
    bad['sectors']['T']['refocusing'] = 'EXACT'
    assert not probe.replay(bad, degrees=[])


def test_failed_run_withdraws_verdicts(tmp_path, monkeypatch):
    def boom():
        raise ArithmeticError('forced')
    monkeypatch.setattr(probe, 'run', boom)
    out = tmp_path/'r.json'
    monkeypatch.setattr('sys.argv', ['x', '--output', str(out)])
    with pytest.raises(ArithmeticError):
        probe.main()
    data = json.loads(out.read_text())
    assert set(data['verdicts'].values()) == {'UNRESOLVED'}
