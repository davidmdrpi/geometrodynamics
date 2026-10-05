import json
import numpy as np
import pytest
from geometrodynamics.waves import r3_breaking as b

R = .2   # toy resonant circle: I = .02


def _c(th):
    return np.array([0, 0, R*np.cos(th), -R*np.sin(th)])


def _dc(th):
    return np.array([0, 0, -R*np.sin(th), -R*np.cos(th)])


def _scan(eps, n=60, seed_error=0.):
    P = b.toy_map(eps)
    lams = []
    for j in range(n):
        ph = 2*np.pi*j/n
        seeds = [_c(ph+i*4*np.pi/5)*(1+seed_error) for i in range(5)]
        lams.append(b.scan_point(P, seeds, _c(ph), _dc(ph))['lam'])
    return np.array(lams)


def test_unbroken_toy_gives_zero_lambda_everywhere():
    lam = _scan(0., seed_error=1e-3)
    out = b.classify(lam, 1e-14, 5)
    assert out['Lambda'] < 1e-13 and out['label'] == 'UNBROKEN_LOOP'


@pytest.mark.parametrize('eps', [1e-9, 1e-6])
def test_broken_toy_is_a_chain_with_the_analytic_amplitude(eps):
    lam = _scan(eps)
    out = b.classify(lam, 1e-14, 5)
    # five kicks of eps*sin(5 phi) in I; dI = R dr, so the radial obstruction is 5 eps/R
    assert out['label'] == 'BROKEN_CHAIN' and out['sign_changes'] == 10 and out['dominant_harmonic'] == 5
    assert abs(out['Lambda']/(5*eps/R)-1) < 1e-3


def test_classification_rule_bands():
    assert b.classify(np.full(60, 1e-12), 1e-13, 5)['label'] == 'UNBROKEN_LOOP'
    assert b.classify(1e-9*np.sin(5*np.linspace(0, 2*np.pi, 60, endpoint=False)), 1e-13, 5)['label'] == 'BROKEN_CHAIN'
    assert b.classify(5e-11*np.sin(5*np.linspace(0, 2*np.pi, 60, endpoint=False)), 1e-13, 5)['label'] == 'INDETERMINATE'
    # structure without the q-fold sign pattern is not a chain
    assert b.classify(1e-6*np.sin(np.linspace(0, 2*np.pi, 60, endpoint=False)), 1e-13, 5)['label'] == 'INDETERMINATE'
    # a resolution above the cap cannot certify an unbroken loop
    assert b.classify(np.zeros(60), 1e-9, 5)['label'] == 'INDETERMINATE'


def test_centre_block_recovers_product_multipliers_without_hyperbolic_products():
    rng = np.random.default_rng(1)
    for a in (.3, 1.1, None):
        Cm = np.array([[1, .7], [0, 1.]]) if a is None else np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
        D = np.zeros((4, 4))
        D[0, 0], D[1:3, 1:3], D[3, 3] = 85., Cm, 1/85
        S = [np.linalg.qr(rng.normal(size=(4, 4)))[0] @ np.diag([1, 2, .5, 1.5]) for _ in range(5)]
        blocks = [S[(i+1) % 5] @ D @ np.linalg.inv(S[i]) for i in range(5)]
        C = b.centre_block(blocks)
        assert abs(np.trace(C)-np.trace(np.linalg.matrix_power(Cm, 5))) < 1e-10 and abs(np.linalg.det(C)-1) < 1e-10


def test_chain_orbits_have_opposite_residues():
    P = b.toy_map(1e-6)
    res = []
    for ph in (0., np.pi/5):
        o = b.periodic_orbit(P, [_c(ph+.02+i*4*np.pi/5) for i in range(5)])
        assert o['residual'] < 1e-12
        res.append((2-np.trace(b.centre_block([b.fd_jacobian(P, np.array(z)) for z in o['nodes']])))/4)
    # R = -(1/4) tr-shift: |R| = (10 pi |nu|)(25 eps)/4 to leading order
    assert res[0]*res[1] < 0 and abs(abs(res[0])/(10*np.pi*25e-6/4)-1) < .02


def test_curves_interpolate_and_close():
    th = 2*np.pi*np.arange(63)/63
    K = np.c_[1+.1*np.cos(th), .2*np.sin(2*th), .3*np.cos(th), -.3*np.sin(th)]
    c, dc = b.trig_curve(K)
    x = .7
    assert np.abs(c(x)-[1+.1*np.cos(x), .2*np.sin(2*x), .3*np.cos(x), -.3*np.sin(x)]).max() < 1e-13
    assert np.abs(dc(x)-[-.1*np.sin(x), .4*np.cos(2*x), -.3*np.sin(x), -.3*np.cos(x)]).max() < 1e-12
    s, ds = b.spline_curve(K)
    assert np.abs(s(0.)-K[0]).max() < 1e-14 and np.abs(s(2*np.pi-1e-12)-K[0]).max() < 1e-10


_TOY = b.toy_map(1e-7)


def _toy_alt(z):
    return _TOY(z)+1e-15


def _zero(z):
    return 0.


def _fake_residue(nodes, steps):
    return 0., [[1., 1.], [0., 0.]]


def test_probe_pipeline_on_toy(tmp_path, monkeypatch):
    from experiments.closure_ledger import r3_breaking_probe as probe
    P = _TOY
    monkeypatch.setattr(probe, 'RUN_DIR', tmp_path)
    monkeypatch.setattr(probe, 'N_PHASE', 30)
    monkeypatch.setattr(probe, 'NOISE_PHASES', (0, 15))
    for name in ('P4', 'P6'):
        monkeypatch.setattr(probe, name, P)
    for name in ('P4_radau', 'P6_radau'):
        monkeypatch.setattr(probe, name, _toy_alt)
    for name in ('constraint4', 'constraint6'):
        monkeypatch.setattr(probe, name, _zero)
    seeds = lambda ph: ([_c(ph+i*4*np.pi/5) for i in range(5)], _c(ph), _dc(ph))
    monkeypatch.setattr(probe, 'main_seeds', seeds)
    monkeypatch.setattr(probe, 'control_seeds', seeds)
    monkeypatch.setattr(probe, 'lrs_circle', lambda: (None, .5, 2.6, 2.4))
    monkeypatch.setattr(probe, '_jet_residue', _fake_residue)
    probe.main(['control', 'main', 'noise', 'orbits'])
    out = probe.score()
    assert out['main']['label'] == 'BROKEN_CHAIN' and out['main']['orbit_gate'] and out['main']['distinct_orbits'] == 2
    rec = json.loads((tmp_path/'result.json').read_text())
    assert rec['sources'] == probe.sources()


_RUN = __import__('pathlib').Path(__file__).resolve().parents[1]/'experiments/closure_ledger/runs/20261005_r3_breaking'


@pytest.mark.skipif(not (_RUN/'result.json').exists(), reason='breaking scan not run')
def test_authenticated_replay_verifies():
    from experiments.closure_ledger import r3_breaking_replay as replay
    out = replay.replay()
    rec = json.loads((_RUN/'result.json').read_text())['result']
    assert out['replay'] == 'VERIFIED' and out['main'] == rec['main']['label'] and out['control'] == rec['control']['label']


@pytest.mark.skipif(not (_RUN/'result.json').exists(), reason='breaking scan not run')
@pytest.mark.parametrize('damage', ['bytes', 'residual', 'phase', 'shape', 'nan'])
def test_replay_rejects_damaged_archives(tmp_path, damage):
    import shutil
    from experiments.closure_ledger import r3_breaking_replay as replay
    for name in replay.SHA256:
        shutil.copy(_RUN/name, tmp_path/name)
    if damage == 'bytes':
        p = tmp_path/'scan_main.json'
        p.write_text(p.read_text().replace('"ok": true', '"ok": false', 1))
        with pytest.raises(ValueError, match='fingerprint'):
            replay.replay(directory=tmp_path)
        return
    rec = json.loads((tmp_path/'scan_main.json').read_text())
    p = next(x for x in rec['points'] if x['ok'])
    if damage == 'residual':
        p['residual'] = 1e-6
    elif damage == 'phase':
        p['phi'] += 1e-3
    elif damage == 'shape':
        p['nodes'] = p['nodes'][:-1]
    else:
        p['lam'] = float('nan')
    with pytest.raises(ValueError):
        replay.validate_scan(rec, 5, 4)


@pytest.mark.skipif(not (_RUN/'result.json').exists(), reason='breaking scan not run')
def test_full_reevaluation_rejects_altered_lambda():
    from experiments.closure_ledger import r3_breaking_probe as probe
    from experiments.closure_ledger import r3_breaking_replay as replay
    rec = json.loads((_RUN/'scan_main.json').read_text())
    rec['points'] = [p for p in rec['points'] if p['ok']][:2]
    assert replay.reevaluate(rec, probe.P4) < 1e-10
    rec['points'][0]['lam'] += 1e-8
    with pytest.raises(ValueError, match='re-close'):
        replay.reevaluate(rec, probe.P4)
