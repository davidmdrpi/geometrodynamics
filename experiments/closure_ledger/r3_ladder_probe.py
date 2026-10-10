"""Resonance-breaking ladder across the LRS circle family at high precision (docs/r3_ladder_prereg.md).

The section map is a square, P = h o h (h: half-return map, lrs_taylor.half_map), so the P-resonance
p/q is the h-resonance (p+q)/(2q) of order Q_h = denominator of (p+q)/(2q). Stages:
  scan     double-precision phase scan of every rung (frozen r3_breaking machinery, DOP853 map)
  hp       chord re-solve of every converged scan point with the 48-digit Taylor map (config C1)
  hpnoise  continue each C1 solution with the 60-digit map (config C2); noise = |dlambda| + residuals
  score    registered labels
Usage: python -m experiments.closure_ledger.r3_ladder_probe [stage ...]
"""
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path
import gmpy2
from gmpy2 import mpfr
import numpy as np
from geometrodynamics.waves import r3_breaking as b
from geometrodynamics.waves import lrs_taylor as lt
from experiments.closure_ledger import r3_breaking_probe as bp

ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = ROOT/'experiments/closure_ledger/runs/20261010_r3_ladder'
SOURCES = bp.SOURCES+('geometrodynamics/waves/lrs_taylor.py', 'experiments/closure_ledger/r3_ladder_probe.py',
                      'docs/r3_ladder_prereg.md')
EXT = ROOT/'experiments/closure_ledger/runs/20260929_r3_extension/part_A.json'
RUNGS = ((5, 11), (4, 9), (3, 7), (5, 12), (2, 5), (3, 8), (4, 11))
N_PHASE = 60
TOL, ACCEPT = 1e-12, 1e-11            # double-precision scan, as r3_breaking
MAX_FAILED = 6
HP_TOL = {'C1': 1e-35, 'C2': 1e-45}
HP_MAXIT = 14
HP_FLOOR = 1e-30
FD_HP = 1e-20
# anchors (archived, docs/r3_breaking.md): Lambda at the LRS 2/5 crossing, and its h-order
LAMBDA_25 = 1.57e-11
DECISIVE = (3, 7)
WINDOW = (1e-2, 1e2)                  # ORDINARY window around the anchored predictions


def sources():
    return {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in SOURCES}


def utc():
    return datetime.now(timezone.utc).isoformat()


def tag(p, q):
    return f'{p}_{q}'


def h_order(p, q):
    return Fraction(p+q, 2*q).denominator


def _write(name, rec):
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    path = RUN_DIR/name
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(dict(sources=sources(), **rec), indent=1))


def _read(name):
    return json.loads((RUN_DIR/name).read_text())


def bracket(p, q):
    """Consecutive accepted archived circles whose rotation numbers bracket p/q; linear interpolation in omega."""
    circles = sorted([c for c in json.loads(EXT.read_text())['circles'] if c.get('ok')], key=lambda c: c['a'])
    target = 2*np.pi*p/q
    for c1, c2 in zip(circles, circles[1:]):
        if (c1['omega']-target)*(c2['omega']-target) <= 0:
            s = (c1['omega']-target)/(c1['omega']-c2['omega'])
            K = np.array(c1['K'])+s*(np.array(c2['K'])-np.array(c1['K']))
            return dict(K=K, s=float(s), a=float(c1['a']+s*(c2['a']-c1['a'])),
                        action=float(c1['action']+s*(c2['action']-c1['action'])),
                        a_bracket=[c1['a'], c2['a']], omega_bracket=[c1['omega'], c2['omega']], omega_star=target)
    raise ValueError(f'no archived bracket for {p}/{q}')


def predictions():
    """Anchored at 2/5 with exponent Q_h: M1 uniform coefficient eps a^Q_h, M2 uniform radius (a/R)^Q_h."""
    a25 = bracket(2, 5)['a']
    eps, R = LAMBDA_25/a25**10, a25/LAMBDA_25**.1
    out = {}
    for p, q in RUNGS:
        a, Q = bracket(p, q)['a'], h_order(p, q)
        out[tag(p, q)] = dict(Q_h=Q, a=a, M1=eps*a**Q, M2=(a/R)**Q)
    return dict(eps=eps, R=R, rungs=out)


# ---------- double-precision scan ----------
def _scan_one(args):
    p, q, j = args
    br = bracket(p, q)
    c, dc = b.trig_curve(br['K'])
    phi = 2*np.pi*j/N_PHASE
    seeds = [c(phi+i*br['omega_star']) for i in range(q)]
    t = time.time()
    try:
        pt = b.scan_point(bp.P4, seeds, c(phi), dc(phi), tol=TOL)
        pt['constraint_max'] = max(bp.constraint4(z) for z in pt['nodes'])
        pt['seed_distance'] = float(max(np.abs(np.array(z)-s).max() for z, s in zip(pt['nodes'], seeds)))
        pt['ok'] = bool(pt['residual'] <= ACCEPT)
    except (ArithmeticError, ValueError, np.linalg.LinAlgError) as e:
        pt = dict(ok=False, error=f'{type(e).__name__}: {e}')
    pt.update(j=j, phi=phi, seconds=time.time()-t)
    return pt


def stage_scan():
    for p, q in RUNGS:
        start = utc()
        with Pool(4) as pool:
            pts = pool.map(_scan_one, [(p, q, j) for j in range(N_PHASE)])
        br = bracket(p, q)
        _write(f'scan_{tag(p, q)}.json', dict(p=p, q=q, started_utc=start, finished_utc=utc(), points=pts,
                                               **{k: v for k, v in br.items() if k != 'K'}))
        print(f'scan {p}/{q} done', flush=True)


# ---------- high-precision chord ----------
def hp_P(cfg):
    return lambda z: lt.hp_map(z, **cfg)[0]


def hp_residual(P, Z, lam, g, t0, c0):
    """The scan system of r3_breaking._residual in mpfr."""
    q, n = len(Z), len(g)
    F = []
    for i in range(q):
        F.extend(a-c for a, c in zip(P(Z[i]), Z[(i+1) % q]))
    for k in range(n):
        F[(q-1)*n+k] -= lam*g[k]
    F.append(gmpy2.fsum(t0[k]*(Z[0][k]-c0[k]) for k in range(n)))
    return F


def hp_jacobian(P, Z, g, t0):
    blocks = []
    h = mpfr(FD_HP)
    for z in Z:
        cols = []
        for k in range(len(z)):
            zp, zm = list(z), list(z)
            zp[k] += h
            zm[k] -= h
            cols.append([float((a-c)/(2*h)) for a, c in zip(P(zp), P(zm))])
        blocks.append(np.array(cols).T)
    return b._jacobian(blocks, np.array([float(v) for v in g]), np.array([float(v) for v in t0]))


def hp_chord(point, cfg, tol, maxit=HP_MAXIT, start=None, P=None):
    """Mixed-precision chord Newton: residual in mpfr at cfg['bits'], correction with the archived double
    Jacobian; switches once to a high-precision FD Jacobian if an iteration contracts by less than 1e3."""
    P = P or hp_P(cfg)
    with gmpy2.context(gmpy2.get_context(), precision=cfg['bits']):
        src = start or point
        Z = [[mpfr(v) for v in z] for z in src['nodes']]
        lam = mpfr(src['lam'])
        g, t0, c0 = ([mpfr(v) for v in point[k]] for k in ('g', 't0', 'c0'))
        J, jac, hist = np.array(point['J']), 'double', []
        n = len(Z[0])
        for it in range(maxit):
            F = hp_residual(P, Z, lam, g, t0, c0)
            hist.append(max(abs(f) for f in F))
            if hist[-1] < tol:
                break
            if jac == 'double' and len(hist) >= 2 and hist[-1] > hist[-2]*1e-3:
                J, jac = hp_jacobian(P, Z, g, t0), 'hp'
            dX = np.linalg.solve(J, -np.array([float(f) for f in F]))
            Z = [[z[k]+mpfr(dX[i*n+k]) for k in range(n)] for i, z in enumerate(Z)]
            lam = lam+mpfr(dX[-1])
        return dict(nodes=[[str(v) for v in z] for z in Z], lam=str(lam), lam_float=float(lam),
                    residual=float(hist[-1]), history=[float(h) for h in hist], jac=jac, iterations=it+1)


def _hp_one(args):
    p, q, j = args
    pt = _read(f'scan_{tag(p, q)}.json')['points'][j]
    if not pt['ok']:
        return dict(j=j, ok=False)
    t = time.time()
    try:
        r = hp_chord(pt, lt.C1, HP_TOL['C1'])
        r['ok'] = bool(r['residual'] < HP_TOL['C1'])
    except (ArithmeticError, ValueError, np.linalg.LinAlgError) as e:
        r = dict(ok=False, error=f'{type(e).__name__}: {e}')
    r.update(j=j, seconds=time.time()-t)
    return r


def _hp2_one(args):
    p, q, j = args
    pt = _read(f'scan_{tag(p, q)}.json')['points'][j]
    r1 = _read(f'hp_{tag(p, q)}.json')['rows'][j]
    if not r1['ok']:
        return dict(j=j, ok=False)
    t = time.time()
    try:
        r = hp_chord(pt, lt.C2, HP_TOL['C2'], start=r1)
        r['ok'] = bool(r['residual'] < HP_TOL['C2'])
    except (ArithmeticError, ValueError, np.linalg.LinAlgError) as e:
        r = dict(ok=False, error=f'{type(e).__name__}: {e}')
    r.update(j=j, seconds=time.time()-t)
    return r


def _stage_hp(name, fn):
    for p, q in RUNGS:
        start = utc()
        with Pool(4) as pool:
            rows = pool.map(fn, [(p, q, j) for j in range(N_PHASE)])
        _write(f'{name}_{tag(p, q)}.json', dict(p=p, q=q, started_utc=start, finished_utc=utc(), rows=rows))
        print(f'{name} {p}/{q} done', flush=True)


# ---------- scoring ----------
def score_rung(scan, hp1, hp2):
    p, q = scan['p'], scan['q']
    Q = h_order(p, q)
    pts = scan['points']
    good = [j for j in range(N_PHASE) if pts[j]['ok'] and hp1[j]['ok'] and hp2[j]['ok']]
    out = dict(p=p, q=q, Q_h=Q, a=scan['a'], action=scan['action'], usable=len(good))
    if len(good) < N_PHASE-MAX_FAILED:
        out.update(status='INDETERMINATE', reason='too many failed phases')
        return out
    with gmpy2.context(gmpy2.get_context(), precision=lt.C2['bits']):
        d12 = max(abs(mpfr(hp2[j]['lam'])-mpfr(hp1[j]['lam'])) for j in good)
    nu = float(d12)+max(max(hp1[j]['residual'], hp2[j]['residual']) for j in good)
    r = max(10*nu, HP_FLOOR)
    lam = np.array([hp1[j]['lam_float'] for j in good])
    Lam = float(np.abs(lam).max())
    resolved = bool(Lam >= 10*r)
    out.update(noise=nu, resolution=r, Lambda=Lam, resolved=resolved, upper=Lam if resolved else 10*r,
               status='RESOLVED' if resolved else 'UNRESOLVED',
               double_error=float(max(abs(pts[j]['lam']-hp1[j]['lam_float']) for j in good)),
               Lambda_double=float(max(abs(pts[j]['lam']) for j in good)),
               sign_changes=b.sign_changes(lam))
    if len(good) == N_PHASE:
        spec = np.abs(np.fft.rfft(lam))/N_PHASE
        out['harmonics'] = [float(x) for x in spec]
        out['dominant_harmonic'] = int(np.argmax(spec[1:])+1)
    return out


def primary(r37, pred):
    lo = WINDOW[0]*min(pred['M1'], pred['M2'])
    hi = WINDOW[1]*max(pred['M1'], pred['M2'])
    if r37['status'] == 'INDETERMINATE':
        return 'INCONCLUSIVE', lo, hi
    if r37['resolved'] and lo <= r37['Lambda'] <= hi:
        return 'ORDINARY_BREAKING', lo, hi
    if r37['upper'] < lo:
        return 'ANOMALOUS_SUPPRESSION', lo, hi
    if r37['resolved'] and r37['Lambda'] > hi:
        return 'ENHANCED_BREAKING', lo, hi
    return 'INCONCLUSIVE', lo, hi


def signal_25(r25):
    if r25['status'] == 'INDETERMINATE':
        return 'OTHER'
    if (r25['resolved'] and 5e-12 <= r25['Lambda'] <= 5e-11 and r25['double_error'] <= 5e-12
            and r25.get('dominant_harmonic') == 10):
        return 'SIGNAL_CONFIRMED'
    if r25['upper'] <= 1.6e-13:
        return 'SIGNAL_ARTEFACT'
    return 'OTHER'


def harmonic_selection(rungs):
    res = [r for r in rungs.values() if r.get('resolved') and 'dominant_harmonic' in r]
    if len(res) < 2:
        return 'UNTESTED'
    k = [(r['dominant_harmonic'], r['Q_h']) for r in res]
    if all(d == Q for d, Q in k):
        return 'HALF_MAP_SELECTION'
    if all(d % Q == 0 for d, Q in k):
        return 'EXTRA_SELECTION'
    return 'VIOLATED'


def exponent_fit(rungs):
    """y = log10 Lambda - n log10 a = alpha - n log10 R over resolved rungs, n = Q_h or q; RMS in decades."""
    res = [r for r in rungs.values() if r.get('resolved')]
    if len(res) < 4:
        return dict(label='UNTESTED', resolved=len(res))
    out = dict(resolved=len(res))
    for key, n in (('Q_h', [r['Q_h'] for r in res]), ('q', [r['q'] for r in res])):
        n = np.array(n, float)
        y = np.array([np.log10(r['Lambda'])-m*np.log10(r['a']) for r, m in zip(res, n)])
        X = np.c_[np.ones_like(n), -n]
        coef = np.linalg.lstsq(X, y, rcond=None)[0]
        rms = float(np.sqrt(np.mean((X @ coef-y)**2)))
        out[key] = dict(alpha=float(coef[0]), R=float(10**coef[1]), rms=rms)
    h, q = out['Q_h']['rms'], out['q']['rms']
    if h <= 1 and h <= q/2:
        out['label'] = 'EXPONENT_QH'
    elif q <= 1 and q <= h/2:
        out['label'] = 'EXPONENT_Q'
    elif h > 1 and q > 1:
        out['label'] = 'NEITHER_FITS'
    else:
        out['label'] = 'UNDISCRIMINATED'
    return out


def score():
    pred = predictions()
    rungs = {}
    for p, q in RUNGS:
        t = tag(p, q)
        r = score_rung(_read(f'scan_{t}.json'), _read(f'hp_{t}.json')['rows'], _read(f'hpnoise_{t}.json')['rows'])
        r['prediction'] = pred['rungs'][t]
        if r.get('status') != 'INDETERMINATE':
            for m in ('M1', 'M2'):
                r[f'log10_dev_{m}'] = float(np.log10(r['upper']/pred['rungs'][t][m]))
        rungs[t] = r
    label, lo, hi = primary(rungs[tag(*DECISIVE)], pred['rungs'][tag(*DECISIVE)])
    out = dict(primary=label, window=[lo, hi], signal_2_5=signal_25(rungs[tag(2, 5)]),
               harmonic_selection=harmonic_selection(rungs), exponent_fit=exponent_fit(rungs),
               exact_integrability=('INTEGRABLE_TO_HP_RESOLUTION' if not any(r.get('resolved') for r in rungs.values())
                                    else 'BREAKING_RESOLVED'),
               anchors=dict(eps=pred['eps'], R=pred['R']), rungs=rungs)
    _write('result.json', dict(result=out))
    return out


def main(argv):
    stages = argv or ['scan', 'hp', 'hpnoise', 'score']
    if 'scan' in stages:
        _write('started.json', dict(started_utc=utc(), stages=stages))
    for s in stages:
        t = time.time()
        if s == 'scan':
            stage_scan()
        elif s == 'hp':
            _stage_hp('hp', _hp_one)
        elif s == 'hpnoise':
            _stage_hp('hpnoise', _hp2_one)
        elif s == 'score':
            r = score()
            print(json.dumps({k: v for k, v in r.items() if k != 'rungs'}, indent=1))
            for k, v in r['rungs'].items():
                print(k, v.get('status'), {x: v.get(x) for x in ('Q_h', 'Lambda', 'resolution', 'dominant_harmonic',
                                                                    'sign_changes', 'double_error')})
        print(f'stage {s}: {time.time()-t:.0f} s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
