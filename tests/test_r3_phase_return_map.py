import numpy as np
import pytest
import copy
import json
from scipy.integrate import solve_ivp
from geometrodynamics.waves.taylor_jets import ring
from geometrodynamics.waves import r3_phase_return_map as rm


def synthetic(twist, warped=False):
    algebra = ring(4)
    a, b, q, p = [algebra.variable(i) for i in range(4)]
    angle = 2*np.pi*.37
    R = np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])
    center = rm.transform(R, [q, p])
    action = (q*q+p*p)/2
    F = [2*a, b/2, center[0]+2*np.pi*twist*action*center[1],
         center[1]-2*np.pi*twist*action*center[0]]
    if warped:
        # Two exact canonical shears, so the centre graph has a nonconstant
        # induced area density and the linear action normalization is fixed.
        def C(v):
            a, b, q, p = v
            a, b, q, p = a+.3*q*q, b, q, p-.6*q*b
            return [a, b+.2*p*p, q-.4*a*p, p]
        def inverse(v):
            a, b, q, p = v
            a, b, q, p = a, b-.2*p*p, q+.4*a*p, p
            return [a-.3*q*q, b, q, p+.6*q*b]
        F = C([v.compose(inverse([a, b, q, p])) for v in F])
    return np.array([v.c for v in F])


@pytest.mark.parametrize('twist', [-.4, 0., .3])
@pytest.mark.parametrize('warped', [False, True])
def test_known_twist_survives_canonical_embedding(twist, warped):
    coefficients = synthetic(twist, warped)
    result = rm.reduce_map(coefficients)
    assert abs(result['nu']-twist) < 1e-11
    assert abs(result['radial_resonant']) < 1e-11
    assert result['graph_residual'] < 1e-11
    assert rm.symplectic_residual(coefficients)[0] < 1e-11


def test_taylor_arithmetic_and_degree_four_remainder():
    algebra = ring(2)
    x, y = [algebra.variable(i) for i in range(2)]
    f = (x+2*y).exp()/((2+x-y)**.5)
    errors = []
    for h in [.02, .01, .005]:
        point = np.array([h, -.3*h])
        exact = np.exp(point[0]+2*point[1])/np.sqrt(2+point[0]-point[1])
        errors.append(abs(f.evaluate(point)-exact))
    assert min(np.log2(np.array(errors[:-1])/errors[1:])) > 3.9
    np.testing.assert_allclose((f*(2+x-y)**.5).c, (x+2*y).exp().c, atol=1e-14)


def test_phase_reduction_matches_unmodified_full_vector_field():
    y0 = rm.initial_full(np.array([.002, .003, .025, -.02]))
    sol = solve_ivp(lambda t, y: rm.full.conformal_rhs(y), (0, .3), y0,
                    method='DOP853', rtol=1e-12, atol=1e-14)
    y = sol.y[:, -1]
    A, u, q, p, M, L = rm.full.unpack(y)
    x = np.log(np.diag(M))@np.diag(rm.B0)/2
    v = np.trace(L@rm.B0)
    phase = np.arctan2(-p[0]/2, q[0])
    full_rhs = rm.full.conformal_rhs(y)
    clock = (p[0]**2-q[0]*full_rhs[7])/(2*(q[0]**2+p[0]**2/4))
    reduced = np.array(rm.phase_rhs(phase, [A, u, x, v, y[2]]))
    expected = [u, full_rhs[1], v, np.trace(full_rhs[20:29].reshape(3, 3)@rm.B0), 1.]
    np.testing.assert_allclose(clock*reduced, expected, atol=2e-12, rtol=2e-12)


@pytest.fixture(scope='module')
def raw_record():
    from experiments.closure_ledger import r3_phase_return_map_probe as probe
    return probe.load_raw(probe.RUN/'raw.json.gz.b64')


def test_archive_replay_authenticates_sources_and_frozen_decision():
    from experiments.closure_ledger.r3_phase_return_map_replay import replay
    result = replay()
    assert result['replay'] == 'VERIFIED'
    assert result['verdict'] == 'SHIFT_AWAY_FROM_TARGET'
    assert all(result['gates'].values())


@pytest.mark.parametrize('damage', ['missing', 'nan', 'preparation', 'endpoint'])
def test_incomplete_or_altered_raw_evidence_is_rejected(raw_record, damage):
    from experiments.closure_ledger import r3_phase_return_map_probe as probe
    raw = copy.deepcopy(raw_record)
    if damage == 'missing':
        raw['validation'].pop()
    elif damage == 'nan':
        raw['jets'][0]['coefficients'][0][1] = float('nan')
    elif damage == 'preparation':
        raw['validation'][0]['z'][0] += .001
    else:
        raw['validation'][0]['history']['returned'][0] += .001
    with pytest.raises(ValueError):
        probe.analyze(raw)


def test_damaged_interior_state_fails_numerical_gate(raw_record):
    from experiments.closure_ledger import r3_phase_return_map_probe as probe
    raw = copy.deepcopy(raw_record)
    raw['validation'][0]['history']['samples'][100][1] += .1
    result = probe.analyze(raw)
    assert not result['gates']['G5'] and result['verdict'] == 'UNRESOLVED'


def test_unresolved_coefficient_does_not_become_a_physical_verdict(raw_record, monkeypatch):
    from experiments.closure_ledger import r3_phase_return_map_probe as probe
    original = rm.reduce_map
    def zero_twist(coefficients):
        result = original(coefficients)
        result['nu'] = 0.
        return result
    monkeypatch.setattr(rm, 'reduce_map', zero_twist)
    result = probe.analyze(raw_record)
    assert not result['gates']['G4'] and result['verdict'] == 'UNRESOLVED'


def test_modified_verdict_cannot_authenticate(tmp_path):
    from experiments.closure_ledger import r3_phase_return_map_probe as probe
    from experiments.closure_ledger.r3_phase_return_map_replay import replay
    (tmp_path/'raw.json.gz.b64').write_bytes((probe.RUN/'raw.json.gz.b64').read_bytes())
    report = json.loads((probe.RUN/'report.json').read_text())
    report['verdict'] = 'SHIFT_TOWARD_TARGET'
    (tmp_path/'report.json').write_text(json.dumps(report))
    with pytest.raises(ValueError, match='fingerprint'):
        replay(tmp_path)
