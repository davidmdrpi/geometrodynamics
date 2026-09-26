import copy
import hashlib
import json
from pathlib import Path
import numpy as np
import pytest
from geometrodynamics.waves import esu_floquet as fl
from experiments.closure_ledger import esu_floquet_probe as probe
from experiments.closure_ledger import esu_floquet_g3_extension as ext_mod
from experiments.closure_ledger import esu_floquet_symbolic as symbolic

RUN = Path(__file__).resolve().parents[1]/'experiments/closure_ledger/runs/20260926_esu_floquet'
ARCHIVE, EXT, G1 = RUN/'esu_floquet.json', RUN/'g3_extension.json', RUN/'g1_symbolic.json'
have_archive = pytest.mark.skipif(not (ARCHIVE.exists() and G1.exists()), reason='archive not generated')


def test_free_conformal_field_refocuses_exactly():
    for n in (2, 3, 7, 30):
        o = fl.observables('V', n, fl.monodromy('C', n))
        assert o['defect'] < 1e-10 and abs(o['fidelity']-1) < 1e-10


def test_uncoupled_vector_is_free_field():
    for n in (2, 5):
        assert np.allclose(fl.monodromy('V', n, coupled=False), fl.monodromy('C', n), atol=1e-11)


def test_tensor_n2_reproduces_prior_published_map():
    assert abs(np.trace(fl.monodromy('T', 2, t1=np.pi/2))-probe.PRIOR_T2_HALF_TRACE) < 1e-8


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
    assert probe.classify_refocusing([dict(n=n, defect=0., fidelity=1.) for n in range(2, 81)])[0] == 'EXACT'
    absent = [dict(n=n, defect=1., fidelity=.2*(-1)**n) for n in range(2, 81)]
    assert probe.classify_refocusing(absent)[0] == 'ABSENT'


# ---------------------------------------------------------------- G1 records
@have_archive
def test_g1_record_is_bound_complete_and_below_threshold():
    rec = json.loads(G1.read_text())
    assert symbolic.g1_valid(rec)
    for mutate in (lambda r: r['scalar'].__setitem__('3', 1e-6),
                   lambda r: r['vector'].pop('2'),
                   lambda r: r['sources'].__setitem__('geometrodynamics/waves/esu_floquet.py', '0'*64),
                   lambda r: r.__setitem__('threshold', 1.)):
        bad = copy.deepcopy(rec)
        mutate(bad)
        assert not symbolic.g1_valid(bad)


@have_archive
def test_failed_or_missing_g1_withholds_every_certified_verdict():
    data = json.loads(ARCHIVE.read_text())
    bad = json.loads(G1.read_text())
    bad['tensor']['2'] = 1e-3
    res = probe.score(data['raw'], bad)
    assert set(res['verdicts'].values()) == {'UNRESOLVED'}
    assert set(res['verdicts_odd_sector'].values()) == {'UNRESOLVED'}
    assert not probe.replay(data, bad, degrees=[])
    assert set(probe.score(data['raw'], {})['verdicts'].values()) == {'UNRESOLVED'}


# ------------------------------------------------------------ probe replay
@have_archive
def test_archive_replays():
    assert probe.replay(json.loads(ARCHIVE.read_text()), degrees=[2, 17, 80])


def _tampered(edit):
    data = json.loads(ARCHIVE.read_text())
    edit(data)
    return probe.replay(data, degrees=[])


@have_archive
@pytest.mark.parametrize('edit', [
    lambda d: d['result']['verdicts'].__setitem__('S_STABILITY', 'ELLIPTIC_2_TO_80'),
    lambda d: d['result']['verdicts'].__setitem__('T_REFOCUSING', 'EXACT'),
    lambda d: d['result']['verdicts_odd_sector'].__setitem__('T_stability', 'HYPERBOLIC_AT[2]'),
    lambda d: d['result']['sectors']['S'].__setitem__('stability', 'ELLIPTIC_2_TO_80'),
    lambda d: d['result']['sectors']['T']['wkb'].__setitem__('verdict', 'FAIL'),
    lambda d: d['result']['sectors']['V']['gates'].__setitem__('G3', True),
    lambda d: d['result']['sectors']['T']['gates'].__setitem__('G2', False),
    lambda d: d['result']['sectors']['S']['gates'].__setitem__('G4', False),
    lambda d: d['result']['controls'].__setitem__('pass', False),
    lambda d: d['result']['sectors']['S']['rows'][0].__setitem__('stability', 'ELLIPTIC'),
    lambda d: d['raw']['sectors']['S']['rows'][0]['half_period_matrix'][0].__setitem__(0, 0.5),
    lambda d: d['raw']['sectors']['S']['rows'][5]['matrix'][0].__setitem__(0, 1.0),
    lambda d: d['raw']['controls']['C1_defects'].__setitem__(3, 1.0),
    lambda d: d['raw']['sectors']['V']['rows'][70].__setitem__('rk4_16_diff', 1e-6),
])
def test_replay_rejects_altered_verdicts_gates_and_evidence(edit):
    assert not _tampered(edit)


@have_archive
def test_consistently_forged_failed_raw_evidence_withdraws_verdicts():
    """Raw failure plus a self-consistent result still cannot certify."""
    data = json.loads(ARCHIVE.read_text())
    g1 = json.loads(G1.read_text())
    data['raw']['controls']['C1_defects'][3] = 1.0
    data['result'] = probe.score(data['raw'], g1)
    assert set(data['result']['verdicts'].values()) == {'UNRESOLVED'}
    assert not probe.close(data['raw']['controls'], probe.measure_controls(), 1e-9)  # remeasurement exposes it


# --------------------------------------------------------------- extension
@have_archive
@pytest.mark.skipif(not EXT.exists(), reason='extension not generated')
def test_extension_replays_and_binds_archive_bytes():
    raw = ARCHIVE.read_bytes()
    ext = json.loads(EXT.read_text())
    assert ext_mod.replay(ext, raw)
    assert ext['archive_sha256'] == hashlib.sha256(raw).hexdigest()
    assert not ext_mod.replay(ext, raw+b' ')
    bad = copy.deepcopy(ext)
    bad['verdicts_ext']['S_STABILITY_EXT'] = 'ELLIPTIC_2_TO_80'
    assert not ext_mod.replay(bad, raw)


@have_archive
@pytest.mark.skipif(not EXT.exists(), reason='extension not generated')
def test_extension_cannot_certify_when_other_gates_fail():
    data = json.loads(ARCHIVE.read_text())
    g1 = json.loads(G1.read_text())
    errors = json.loads(EXT.read_text())['errors']
    good = ext_mod.score(data, g1, errors)
    assert good['verdicts_ext']['S_STABILITY_EXT'] == 'HYPERBOLIC_AT[2]'
    for forge in ('controls', 'G2', 'G4', 'G1'):
        d, g = copy.deepcopy(data), copy.deepcopy(g1)
        if forge == 'controls':
            d['raw']['controls']['C3_T2_half_trace'] = 0.
        elif forge == 'G2':
            for X in fl.SECTORS:
                d['raw']['sectors'][X]['rows'][4]['constraint_residual'] = 1e-3
        elif forge == 'G4':
            for X in fl.SECTORS:
                M = np.asarray(d['raw']['sectors'][X]['rows'][4]['matrix'])
                d['raw']['sectors'][X]['rows'][4]['matrix'] = (1.01*M).tolist()
        else:
            g['scalar']['2'] = 1.
        res = ext_mod.score(d, g, errors)
        assert not any(k.endswith(('STABILITY_EXT', 'REFOCUSING_EXT')) for k in res['verdicts_ext']), forge


def test_failed_run_withdraws_verdicts(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise ArithmeticError('forced')
    monkeypatch.setattr(probe, 'run', boom)
    out = tmp_path/'r.json'
    monkeypatch.setattr('sys.argv', ['x', '--output', str(out)])
    with pytest.raises(ArithmeticError):
        probe.main()
    assert set(json.loads(out.read_text())['result']['verdicts'].values()) == {'UNRESOLVED'}
