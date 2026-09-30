import json
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from geometrodynamics.waves import esu_floquet as fl
from geometrodynamics.waves import nonlinear_supported_tt as d
from geometrodynamics.waves import r3_resonance as r3
from experiments.closure_ledger import r3_resonance_probe as probe
from experiments.closure_ledger.esu_floquet_probe import close

ARCHIVE = probe.RUN_DIR/'r3_resonance.json'


def _flow(y, t):
    return solve_ivp(lambda _, z: d.conformal_rhs(z), (0, t), y, method='DOP853', rtol=1e-12, atol=1e-14).y[:, -1]


def test_section_data_satisfy_the_constraint():
    for eps in (0., .04, .16):
        y = r3.section_state(eps, -1)
        assert abs(d.constraints(y)['residual'][0]) < 1e-13 and y[7] < 0


def test_lrs_subsystem_is_invariant():
    y = _flow(r3.section_state(.1, 1), 2.)
    _, _, q, qp, M, L = d.unpack(y)
    assert np.allclose(q[1:], 0, atol=1e-14) and np.allclose(qp[1:], 0, atol=1e-14)
    assert np.allclose(M, np.diag(np.diag(M)), atol=1e-14) and abs(M[0, 0]-M[1, 1]) < 1e-13
    assert np.allclose(L, np.diag(np.diag(L)), atol=1e-14)


def test_linearised_section_map_is_the_310_tensor_n2_map():
    h, T = 1e-6, np.zeros((2, 2))
    for j in range(2):
        y = r3.section_state(0., 1)
        if j == 0:
            y = r3.section_state(h, 1)
        else:
            _, _, _, _, M, L = d.unpack(y)
            y = d.pack(y[0], y[1], y[3:7], y[7:11], M, h*r3.B0)
            y = r3.resolve_clock_velocity(y)
        T[:, j] = np.array(r3.tensor(_flow(y, np.pi)))/h
    assert abs(np.trace(T)-np.trace(fl.monodromy('T', 2))) < 1e-5


def test_birkhoff_average_recovers_a_quasiperiodic_rotation_number():
    k = np.arange(48)
    inc = 2*np.pi*1.4 + .3*np.cos(2*np.pi*.37*k) + .1*np.sin(2*np.pi*.74*k+1)
    assert abs(r3.birkhoff(inc)-1.4) < 1e-6


@pytest.mark.skipif(not ARCHIVE.exists(), reason='archive not generated')
def test_archive_rescores_and_binds_sources():
    rec = json.loads(ARCHIVE.read_text())
    assert rec['freeze'] == probe.FREEZE and rec['correction'] == probe.CORRECTION
    assert rec['sources'] == probe.sources()
    again = json.loads(json.dumps(probe.score(rec['raw'])))
    # Labels exactly; floats to 1e-9 (BLAS summation order differs across platforms).
    assert again['verdict'] == rec['result']['verdict'] and close(again, rec['result'], 1e-9)


@pytest.mark.skipif(not ARCHIVE.exists(), reason='archive not generated')
def test_tampered_increments_change_the_score():
    rec = json.loads(ARCHIVE.read_text())
    run = next(r for r in rec['raw']['runs'] if r['eps'] == .01 and r['pol'] == 1 and r['integrator'] == 'secondary')
    run['increments'] = [x+1e-3 for x in run['increments']]
    assert not close(json.loads(json.dumps(probe.score(rec['raw']))), rec['result'], 1e-9)


@pytest.mark.skipif(not ARCHIVE.exists(), reason='archive not generated')
def test_rescore_tolerates_summation_order(monkeypatch):
    rec = json.loads(ARCHIVE.read_text())

    def reversed_sum(inc):
        inc = np.asarray(inc, float)
        t = (np.arange(len(inc))+.5)/len(inc)
        w = np.exp(-1/(t*(1-t)))
        return float(sum(w[::-1]*inc[::-1])/(2*np.pi*w.sum()))
    monkeypatch.setattr(r3, 'birkhoff', reversed_sum)
    again = json.loads(json.dumps(probe.score(rec['raw'])))
    assert again != rec['result'] and close(again, rec['result'], 1e-9)
