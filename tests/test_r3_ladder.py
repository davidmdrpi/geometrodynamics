import json
import numpy as np
import pytest
gmpy2 = pytest.importorskip('gmpy2')
from gmpy2 import mpfr
import mpmath
from geometrodynamics.waves import lrs_taylor as lt
from geometrodynamics.waves import r3_breaking as b
from geometrodynamics.waves import r3_return_map as rm
from experiments.closure_ledger import r3_ladder_probe as lp

PTS = ([1.0, 0.0, 0.05, 0.0], [0.999, 0.003, 0.08, -0.04])


def _f(v):
    return np.array([float(x) for x in v])


# ---------- the high-precision map ----------
def test_isotropic_sector_is_exact():
    """x = p_x = 0: clock q'' = -4q (return time pi) and A'' = -A + A^3, checked against mpmath's own
    Taylor integrator."""
    z = [1.001, 0.002, 0.0, 0.0]
    with gmpy2.context(gmpy2.get_context(), precision=160):
        P, con, T = lt.hp_map(z)
        assert abs(T-gmpy2.const_pi()) < 1e-36 and con < 1e-38 and abs(P[2]) == 0 and abs(P[3]) == 0
    mpmath.mp.dps = 45
    f = mpmath.odefun(lambda t, y: [y[1], -y[0]+y[0]**3], 0, [mpmath.mpf(1.001), -mpmath.mpf(0.002)/6])
    A, Ap = f(mpmath.pi)
    assert abs(mpmath.mpf(str(P[0]))-A) < 1e-30 and abs(mpmath.mpf(str(P[1]))+6*Ap) < 1e-30
    mpmath.mp.dps = 15


@pytest.mark.parametrize('z', PTS)
def test_hp_map_agrees_with_esu_map_and_with_itself(z):
    a, ca, _ = lt.hp_map(z, **lt.C1)
    c, cc, _ = lt.hp_map(z, **lt.C2)
    assert np.abs(_f(a)-rm.esu_map(np.array(z))[0]).max() < 1e-11
    with gmpy2.context(gmpy2.get_context(), precision=200):
        assert max(abs(x-y) for x, y in zip(a, c)) < 1e-37
    assert ca < 1e-38 and cc < 1e-48


@pytest.mark.parametrize('z', PTS)
def test_section_map_is_the_square_of_the_half_map(z):
    with gmpy2.context(gmpy2.get_context(), precision=160):
        P, _, T = lt.hp_map(z)
        h1, _, t1 = lt.half_map(z)
        h2, _, t2 = lt.half_map(h1)
        assert max(abs(x-y) for x, y in zip(P, h2)) < 1e-37 and abs(T-t1-t2) < 1e-38


def test_hp_map_is_symplectic():
    z = PTS[1]
    with gmpy2.context(gmpy2.get_context(), precision=160):
        h = mpfr(1e-20)
        cols = []
        for k in range(4):
            zp, zm = [mpfr(v) for v in z], [mpfr(v) for v in z]
            zp[k] += h
            zm[k] -= h
            cols.append(_f([(a-c)/(2*h) for a, c in zip(lt.hp_map(zp)[0], lt.hp_map(zm)[0])]))
    J = np.array(cols).T
    assert np.abs(J.T @ rm.OMEGA @ J-rm.OMEGA).max() < 1e-12


# ---------- high-precision chord on the toy map ----------
def toy_hp(eps, rho0=.4, nu=-1., lam_u=85.):
    def P(z):
        u, s, x, p = z
        I, phi = (x*x+p*p)/2, gmpy2.atan2(-p, x)
        I2 = I+mpfr(eps)*gmpy2.sin(5*phi)
        phi2 = phi+2*gmpy2.const_pi()*(mpfr(rho0)+mpfr(nu)*(I2-mpfr(.02)))
        r = gmpy2.sqrt(2*I2)
        return [mpfr(lam_u)*u, s/mpfr(lam_u), r*gmpy2.cos(phi2), -r*gmpy2.sin(phi2)]
    return P


R = .2


def _c(th):
    return np.array([0, 0, R*np.cos(th), -R*np.sin(th)])


def _dc(th):
    return np.array([0, 0, -R*np.sin(th), -R*np.cos(th)])


def _point(eps, ph=.3):
    return b.scan_point(b.toy_map(eps), [_c(ph+i*4*np.pi/5)*(1+1e-3) for i in range(5)], _c(ph), _dc(ph))


def test_hp_chord_unbroken_toy_reaches_zero():
    pt = _point(0.)
    r = lp.hp_chord(pt, lt.C1, 1e-38, P=toy_hp(0.))
    assert r['residual'] < 1e-38 and abs(r['lam_float']) < 1e-38 and abs(pt['lam']) < 1e-13


def test_hp_chord_broken_toy_matches_double_and_converges_across_precisions():
    pt = _point(1e-7)
    r1 = lp.hp_chord(pt, lt.C1, 1e-38, P=toy_hp(1e-7))
    r2 = lp.hp_chord(pt, lt.C2, 1e-48, start=r1, P=toy_hp(1e-7))
    assert abs(r1['lam_float']-pt['lam']) < 1e-13 and abs(r1['lam_float']) > 1e-7
    with gmpy2.context(gmpy2.get_context(), precision=200):
        assert abs(mpfr(r1['lam'])-mpfr(r2['lam'])) < 1e-37
    assert r1['jac'] == 'double'


def test_hp_chord_falls_back_to_hp_jacobian_when_chord_stalls():
    pt = _point(1e-7)
    pt['J'] = (np.array(pt['J'])*1.6).tolist()
    r = lp.hp_chord(pt, lt.C1, 1e-38, P=toy_hp(1e-7))
    assert r['jac'] == 'hp' and r['residual'] < 1e-38


# ---------- registered rules ----------
def test_half_map_orders():
    assert {lp.tag(p, q): lp.h_order(p, q) for p, q in lp.RUNGS} == {
        '5_11': 11, '4_9': 18, '3_7': 7, '5_12': 24, '2_5': 10, '3_8': 16, '4_11': 22}


def test_brackets_and_anchors():
    from experiments.closure_ledger import r3_breaking_probe as bp
    for p, q in lp.RUNGS:
        br = lp.bracket(p, q)
        assert 0 <= br['s'] <= 1 and br['a_bracket'][0] <= br['a'] <= br['a_bracket'][1]
    assert abs(lp.bracket(2, 5)['s']-bp.lrs_circle()[1]) < 1e-15
    pr = lp.predictions()
    assert abs(pr['rungs']['2_5']['M1']/1.57e-11-1) < 1e-12 and abs(pr['rungs']['2_5']['M2']/1.57e-11-1) < 1e-12
    lo, hi = .01*pr['rungs']['3_7']['M1'], 100*pr['rungs']['3_7']['M2']
    assert 2e-12 < lo < 3e-12 and 6e-7 < hi < 7e-7


def _r(Lam, resolved=True, r=1e-30, **kw):
    return dict(status='RESOLVED' if resolved else 'UNRESOLVED', resolved=resolved, Lambda=Lam,
                upper=Lam if resolved else 10*r, **kw)


def test_primary_label_bands():
    pred = dict(M1=2.4e-10, M2=6.5e-9)
    assert lp.primary(_r(1e-9), pred)[0] == 'ORDINARY_BREAKING'
    assert lp.primary(_r(1e-14), pred)[0] == 'ANOMALOUS_SUPPRESSION'
    assert lp.primary(_r(1e-31, resolved=False), pred)[0] == 'ANOMALOUS_SUPPRESSION'
    assert lp.primary(_r(1e-5), pred)[0] == 'ENHANCED_BREAKING'
    assert lp.primary(dict(status='INDETERMINATE'), pred)[0] == 'INCONCLUSIVE'


def test_secondary_labels():
    assert lp.signal_25(_r(1.5e-11, double_error=1e-12, dominant_harmonic=10)) == 'SIGNAL_CONFIRMED'
    assert lp.signal_25(_r(1.5e-11, double_error=1e-12, dominant_harmonic=5)) == 'OTHER'
    assert lp.signal_25(_r(1e-31, resolved=False)) == 'SIGNAL_ARTEFACT'
    rs = {'a': _r(1., Q_h=10, dominant_harmonic=10), 'b': _r(1., Q_h=7, dominant_harmonic=7)}
    assert lp.harmonic_selection(rs) == 'HALF_MAP_SELECTION'
    rs['b']['dominant_harmonic'] = 14
    assert lp.harmonic_selection(rs) == 'EXTRA_SELECTION'
    rs['b']['dominant_harmonic'] = 5
    assert lp.harmonic_selection(rs) == 'VIOLATED'


def test_exponent_fit_prefers_the_generating_exponent():
    rng = np.random.default_rng(0)
    rungs = {}
    for p, q in lp.RUNGS:
        a, Q = lp.bracket(p, q)['a'], lp.h_order(p, q)
        rungs[lp.tag(p, q)] = _r(1e-3*(a/3)**Q*10**(.2*rng.normal()), Q_h=Q, q=q, a=a)
    out = lp.exponent_fit(rungs)
    assert out['label'] == 'EXPONENT_QH' and abs(np.log10(out['Q_h']['R']/3)) < .2
    for p, q in lp.RUNGS:
        r = rungs[lp.tag(p, q)]
        r['Lambda'] = 1e-3*(r['a']/3)**q
    assert lp.exponent_fit(rungs)['label'] == 'EXPONENT_Q'


# ---------- pipeline on the toy ----------
def _toy_bracket(p, q):
    th = 2*np.pi*np.arange(63)/63
    return dict(K=np.array([_c(t) for t in th]), s=.5, a=R, action=R*R/2, a_bracket=[R, R],
                omega_bracket=[4.*np.pi/5, 4.*np.pi/5], omega_star=4*np.pi/5)


def _zero(z):
    return 0.


def test_probe_pipeline_on_toy(tmp_path, monkeypatch):
    from experiments.closure_ledger import r3_breaking_probe as bp
    monkeypatch.setattr(lp, 'RUN_DIR', tmp_path)
    monkeypatch.setattr(lp, 'N_PHASE', 20)
    monkeypatch.setattr(lp, 'RUNGS', ((2, 5),))
    monkeypatch.setattr(lp, 'DECISIVE', (2, 5))
    monkeypatch.setattr(lp, 'bracket', _toy_bracket)
    monkeypatch.setattr(bp, 'P4', b.toy_map(1e-7))
    monkeypatch.setattr(bp, 'constraint4', _zero)
    monkeypatch.setattr(lp, 'hp_P', lambda cfg: toy_hp(1e-7))
    lp.main(['scan', 'hp', 'hpnoise'])
    out = lp.score()
    r = out['rungs']['2_5']
    assert r['status'] == 'RESOLVED' and r['usable'] == 20 and r['noise'] < 1e-34
    assert abs(r['Lambda']/(5e-7/R)-1) < 1e-2 and r['dominant_harmonic'] == 5
    rec = json.loads((tmp_path/'result.json').read_text())
    assert rec['sources'] == lp.sources()
    with pytest.raises(FileExistsError):
        lp.score()


# ---------- authenticated replay (after the run) ----------
_RUN = __import__('pathlib').Path(__file__).resolve().parents[1]/'experiments/closure_ledger/runs/20261010_r3_ladder'


@pytest.mark.skipif(not (_RUN/'result.json').exists(), reason='ladder not run')
def test_ladder_replay_verifies():
    from experiments.closure_ledger import r3_ladder_replay as replay
    out = replay.replay()
    rec = json.loads((_RUN/'result.json').read_text())['result']
    assert out['replay'] == 'VERIFIED' and out['primary'] == rec['primary']


@pytest.mark.skipif(not (_RUN/'result.json').exists(), reason='ladder not run')
@pytest.mark.parametrize('damage', ['bytes', 'residual', 'lam', 'nodes', 'order'])
def test_ladder_replay_rejects_damaged_archives(tmp_path, damage):
    import shutil
    from experiments.closure_ledger import r3_ladder_replay as replay
    for name in replay.SHA256:
        shutil.copy(_RUN/name, tmp_path/name)
    if damage == 'bytes':
        p = tmp_path/'hp_3_7.json'
        p.write_text(p.read_text().replace('"ok": true', '"ok": false', 1))
        with pytest.raises(ValueError, match='fingerprint'):
            replay.replay(directory=tmp_path)
        return
    scan = json.loads((tmp_path/'scan_3_7.json').read_text())
    rows = json.loads((tmp_path/'hp_3_7.json').read_text())['rows']
    r = next(x for x in rows if x['ok'])
    if damage == 'residual':
        r['residual'] = 1e-20
    elif damage == 'lam':
        r['lam_float'] += 1e-9
    elif damage == 'nodes':
        r['nodes'] = r['nodes'][:-1]
    else:
        rows = rows[1:]+rows[:1]
    with pytest.raises(ValueError):
        replay.validate_rows(scan, rows, 1e-35)


@pytest.mark.skipif(not (_RUN/'result.json').exists(), reason='ladder not run')
def test_ladder_reclosure_rejects_altered_lambda():
    from experiments.closure_ledger import r3_ladder_replay as replay
    scan = json.loads((_RUN/'scan_3_7.json').read_text())
    rows = json.loads((_RUN/'hp_3_7.json').read_text())['rows']
    assert replay.reclose(scan, rows, every=30) < 1e-33
    j = next(k for k in range(0, 60, 30) if rows[k]['ok'])
    with gmpy2.context(gmpy2.get_context(), precision=160):
        rows[j]['lam'] = str(gmpy2.mpfr(rows[j]['lam'])+gmpy2.mpfr('1e-30'))
    with pytest.raises(ValueError, match='re-close'):
        replay.reclose(scan, rows, every=30)
