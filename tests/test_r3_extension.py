import numpy as np
import pytest
from scipy.linalg import expm
from geometrodynamics.waves import r3_extension as rx
from geometrodynamics.waves import nonlinear_supported_tt as d
from geometrodynamics.waves.jets import Jet, variables

LAM, A0 = 3., 2.5
B2, B3, G, DL = .7, -.4, .3, .25


def _trig(v, n):
    return getattr(v, n)() if isinstance(v, Jet) else getattr(np, n)(v)


def toy(z):
    """Hyperbolic x two 1:1-resonant elliptic modes; normal-form quartic
    H4 = B2 I2^2 + B3 I3^2 + G I2 I3 + DL (Q2Q3 + P2P3)^2, conjugated by a shear."""
    def R(z):
        Q1, P1, Q2, P2, Q3, P3 = z
        I2, I3 = (Q2*Q2+P2*P2)/2, (Q3*Q3+P3*P3)/2
        rot = lambda Q, P, al: (Q*_trig(al, 'cos')+P*_trig(al, 'sin'), -Q*_trig(al, 'sin')+P*_trig(al, 'cos'))
        return [LAM*Q1, P1/LAM, *rot(Q2, P2, A0+2*B2*I2+G*I3), *rot(Q3, P3, A0+2*B3*I3+G*I2)]

    def PhiG(z):
        Q1, P1, Q2, P2, Q3, P3 = z
        t = 2*DL*(Q2*Q3+P2*P3)
        c, s = _trig(t, 'cos'), _trig(t, 'sin')
        return [Q1, P1, c*Q2+s*P3, c*P2-s*Q3, c*Q3+s*P2, c*P3-s*Q2]

    def SP(z, sg):
        Q1, P1, Q2, P2, Q3, P3 = z
        return [Q1, P1+sg*(-.4*Q2*Q3+.2*Q3*Q3), Q2, P2+sg*(.9*Q2*Q2-.4*Q1*Q3), Q3, P3+sg*(-.4*Q1*Q2+.4*Q3*Q1)]
    return SP(PhiG(R(SP(z, -1))), 1)


def test_multimode_normal_form_recovers_known_quartic():
    nf = rx.multimode_normal_form(toy(variables(np.zeros(6), 3)), np.zeros(6))
    pred = lambda s, t, D: (B2*s*s+B3*t*t+G*s*t+4*DL*s*t*np.cos(D)**2)/np.pi
    for w, args in (((1, 0), (1, 0, 0)), ((0, 1), (0, 1, 0)), ((1, 1), (.5, .5, 0)), ((1, 1j), (.5, .5, np.pi/2)),
                    ((.6, .8*np.exp(.7j)), (.36, .64, .7))):
        assert abs(rx.nu_of(nf['T'], np.array(w), nf['c'])-pred(*args)) < 1e-11
    sup = rx.sup_nu(nf['T'], nf['c'], restarts=16)
    assert abs(sup['max']-B2/np.pi) < 1e-9 and abs(sup['min']-B3/np.pi) < 1e-9
    assert sup['sym2_upper_bound'] >= sup['max']-1e-12 and sup['dissipative_max'] < 1e-12


def test_matrix_equations_match_full_system():
    rng = np.random.default_rng(5)
    for _ in range(5):
        A, Ap, q, qp = 1+.1*rng.normal(), .1*rng.normal(), .3*rng.normal(), rng.normal()
        M = expm(2*sum(.1*rng.normal()*e for e in rx.E))
        L = rng.normal(size=(3, 3))*.2
        f = d.conformal_rhs(d.pack(A, Ap, np.r_[q, 0, 0, 0], np.r_[qp, 0, 0, 0], M, L))
        _, _, _, _, Mp, Lp = d.unpack(f)
        r = rx.rhs([A, Ap, q, qp, M.tolist(), L.tolist()])
        assert abs(r[1]-f[1]) < 1e-13 and abs(r[3]-f[7]) < 1e-13
        assert np.abs(np.array(r[4])-Mp).max() < 1e-13 and np.abs(np.array(r[5])-Lp).max() < 1e-13


def test_section_parametrisation_round_trips_through_order_three():
    z = variables(np.r_[1., np.zeros(11)], 3)
    back = rx.state_to_section(rx.section_to_state(z))
    assert max(np.abs((b-a).c).max() for a, b in zip(z, back)) < 1e-13
    assert np.abs(rx.constraint(rx.section_to_state(z)).c).max() < 1e-13


def test_so3_representation_is_orthogonal_and_multiplicative():
    rng = np.random.default_rng(2)
    R1, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    R2, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    D1, D2 = rx.so3_rep(R1), rx.so3_rep(R2)
    assert np.allclose(D1.T @ D1, np.eye(5)) and np.allclose(rx.so3_rep(R1 @ R2), D1 @ D2)


RESULT = __import__('pathlib').Path(__file__).resolve().parents[1]/'experiments/closure_ledger/runs/20260929_r3_extension/result.json'


@pytest.mark.skipif(not RESULT.exists(), reason='extension not run')
def test_extension_archive_rescores_and_binds_sources():
    import json
    from experiments.closure_ledger import r3_extension_probe as probe
    from experiments.closure_ledger.esu_floquet_probe import close
    rec = json.loads(RESULT.read_text())
    A = json.loads((probe.RUN_DIR/'part_A.json').read_text())
    B = json.loads((probe.RUN_DIR/'part_B.json').read_text())
    assert A['sources'] == B['sources'] == probe.sources()
    a, b = probe.score_a(A, rec['theta0']), probe.score_b(B, rec['theta0'])
    assert a['label'] == rec['A']['label'] and b['label'] == rec['B']['label'] and b['checks'] == rec['B']['checks']
    assert close(json.loads(json.dumps(a)), rec['A'], 1e-9)
    assert abs(b['nu_max']-rec['B']['nu_max']) < 1e-7 and abs(b['nu_lrs']-rec['B']['nu_lrs']) < 1e-9
    A['circles'][3]['omega'] += 1e-3        # a manufactured turn must be detected
    assert probe.score_a(A, rec['theta0'])['label'] == 'TURN_IN_FAMILY'


@pytest.mark.skipif(not RESULT.exists(), reason='extension not run')
def test_authenticated_replay_verifies():
    from experiments.closure_ledger import r3_extension_replay as replay
    out = replay.replay()
    assert out['replay'] == 'VERIFIED' and out['A'] == 'NO_TURN_IN_FAMILY'
    assert out['B'] == 'SOME_POLARISATION_SHIFTS_TOWARD' and out['accepted_circles'] == 27


def _accepted(A):
    return next(i for i, c in enumerate(A['circles']) if c['ok'])


@pytest.mark.skipif(not RESULT.exists(), reason='extension not run')
@pytest.mark.parametrize('damage', ['residual', 'delete_K', 'nan_K', 'tail', 'action', 'flip_failed_to_ok',
                                    'flip_ok_to_failed', 'off_ladder', 'shape_K'])
def test_replay_validation_rejects_damaged_part_a(damage):
    import copy
    import json
    from experiments.closure_ledger import r3_extension_probe as probe
    from experiments.closure_ledger import r3_extension_replay as replay
    A = copy.deepcopy(json.loads((probe.RUN_DIR/'part_A.json').read_text()))
    i = _accepted(A)
    if damage == 'residual':
        A['circles'][i]['residual'] = 1.0
    elif damage == 'delete_K':
        del A['circles'][i]['K']
    elif damage == 'nan_K':
        A['circles'][i]['K'][0][0] = float('nan')
    elif damage == 'tail':
        A['circles'][i]['fourier_tail'] = 5e-11    # passes the threshold, contradicts K
    elif damage == 'action':
        A['circles'][i]['action'] *= 1.001
    elif damage == 'flip_failed_to_ok':
        j = next(k for k, c in enumerate(A['circles']) if not c['ok'])
        A['circles'][j]['ok'] = True
    elif damage == 'flip_ok_to_failed':
        A['circles'][i]['ok'] = False
    elif damage == 'off_ladder':
        A['circles'][i+1]['a'] *= 1.01
    else:
        A['circles'][i]['K'] = A['circles'][i]['K'][:-1]
    with pytest.raises(ValueError):
        replay.validate_a(A)


@pytest.mark.skipif(not RESULT.exists(), reason='extension not run')
def test_replay_rejects_altered_archive_bytes(tmp_path):
    import shutil
    from experiments.closure_ledger import r3_extension_probe as probe
    from experiments.closure_ledger import r3_extension_replay as replay
    for name in replay.SHA256:
        shutil.copy(probe.RUN_DIR/name, tmp_path/name)
    (tmp_path/'part_A.json').write_text((tmp_path/'part_A.json').read_text().replace('"ok": true', '"ok": false', 1))
    with pytest.raises(ValueError, match='fingerprint'):
        replay.replay(directory=tmp_path)
