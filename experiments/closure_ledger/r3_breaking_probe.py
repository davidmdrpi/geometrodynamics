"""Resonance-breaking scan: LRS 2/5 crossing (main) and the #319 two-return loop (control).

Freeze: docs/r3_breaking_prereg.md. Stages (each writes one JSON binding the source hashes):
  control  phase scan of the #319 diagonal two-return loop (q = 2, 6D), 60 phases
  main     phase scan of the interpolated LRS circle at rotation 2/5 (q = 5, 4D), 60 phases
  noise    Radau chord re-solve at 6 phases of each scan
  orbits   isolated period-5 orbits from every sign-change bracket of the main scan; residues
  score    registered labels
Usage: python -m experiments.closure_ledger.r3_breaking_probe [stage ...]
"""
import hashlib
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from geometrodynamics.waves import r3_breaking as b
from geometrodynamics.waves import r3_family as f
from geometrodynamics.waves import r3_return_map as rm

ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = ROOT/'experiments/closure_ledger/runs/20261005_r3_breaking'
SOURCES = ('geometrodynamics/waves/r3_breaking.py', 'geometrodynamics/waves/r3_family.py',
           'geometrodynamics/waves/r3_return_map.py', 'geometrodynamics/waves/r3_extension.py',
           'geometrodynamics/waves/jets.py', 'geometrodynamics/waves/nonlinear_supported_tt.py',
           'experiments/closure_ledger/r3_breaking_probe.py')
EXT = ROOT/'experiments/closure_ledger/runs/20260929_r3_extension/part_A.json'
FAM = ROOT/'experiments/closure_ledger/runs/20261001_r3_family/stage_F.json'
N_PHASE = 60
NOISE_PHASES = (0, 10, 20, 30, 40, 50)
P_RES, Q_RES = 2, 5
OMEGA_STAR = 2*np.pi*P_RES/Q_RES
A_BRACKET = (.23475303506039574, .2791699795623058)
TOL, ACCEPT = 1e-12, 1e-11
MAX_FAILED = 6


def sources():
    return {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in SOURCES}


def _write(name, rec):
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    rec = dict(sources=sources(), **rec)
    (RUN_DIR/name).write_text(json.dumps(rec, indent=1))


def _read(name):
    return json.loads((RUN_DIR/name).read_text())


# ---------- maps ----------
def P4(z, method='DOP853'):
    return rm.esu_map(np.asarray(z, float), method=method)[0]


def P4_radau(z):
    return P4(z, 'Radau')


def P6(z, method='DOP853'):
    return f.P(np.r_[np.asarray(z, float), np.zeros(6)], method=method)[:6]


def P6_radau(z):
    return P6(z, 'Radau')


def constraint4(z):
    return rm.esu_map(np.asarray(z, float))[1]


def constraint6(z):
    return f.P(np.r_[np.asarray(z, float), np.zeros(6)], full=True)[3]


# ---------- curves ----------
def lrs_circle():
    """Linear interpolation (in omega) of the two archived LRS circles bracketing omega = 4 pi/5."""
    circ = {c['a']: c for c in json.loads(EXT.read_text())['circles'] if c.get('ok')}
    c1, c2 = circ[A_BRACKET[0]], circ[A_BRACKET[1]]
    s = (c1['omega']-OMEGA_STAR)/(c1['omega']-c2['omega'])
    K = np.array(c1['K'])+s*(np.array(c2['K'])-np.array(c1['K']))
    return K, float(s), float(c1['omega']), float(c2['omega'])


def main_seeds(phi):
    K = lrs_circle()[0]
    c, dc = b.trig_curve(K)
    return [c(phi+i*OMEGA_STAR) for i in range(Q_RES)], c(phi), dc(phi)


def control_seeds(phi):
    V = np.array([p['v'] for p in json.loads(FAM.read_text())['points']])
    c0, dc0 = b.spline_curve(V[:, :6])
    c1, _ = b.spline_curve(V[:, 6:])
    return [c0(phi), c1(phi)], c0(phi), dc0(phi)


def phases():
    return [2*np.pi*j/N_PHASE for j in range(N_PHASE)]


# ---------- stages ----------
def _scan_one(args):
    which, j = args
    phi = phases()[j]
    seeds, c0, t0 = (main_seeds if which == 'main' else control_seeds)(phi)
    P, con = (P4, constraint4) if which == 'main' else (P6, constraint6)
    t = time.time()
    try:
        pt = b.scan_point(P, seeds, c0, t0, tol=TOL)
        pt['constraint_max'] = max(con(z) for z in pt['nodes'])
        pt['seed_distance'] = float(max(np.abs(np.array(z)-s).max() for z, s in zip(pt['nodes'], seeds)))
        pt['ok'] = bool(pt['residual'] <= ACCEPT)
    except (ArithmeticError, ValueError, np.linalg.LinAlgError) as e:
        pt = dict(ok=False, error=f'{type(e).__name__}: {e}')
    pt.update(j=j, phi=phi, seconds=time.time()-t)
    return pt


def stage_scan(which):
    with Pool(4) as pool:
        pts = pool.map(_scan_one, [(which, j) for j in range(N_PHASE)])
    extra = {}
    if which == 'main':
        _, s, w1, w2 = lrs_circle()
        extra = dict(interp_s=s, omega_bracket=[w1, w2], a_bracket=list(A_BRACKET), omega_star=OMEGA_STAR)
    _write(f'scan_{which}.json', dict(which=which, q=Q_RES if which == 'main' else 2, points=pts, **extra))


def _noise_one(args):
    which, j = args
    pt = _read(f'scan_{which}.json')['points'][j]
    if not pt['ok']:
        return dict(j=j, ok=False)
    r = b.chord_resolve(P4_radau if which == 'main' else P6_radau, pt, tol=TOL)
    return dict(j=j, ok=True, lam_dop=pt['lam'], lam_radau=r['lam'], radau_residual=r['residual'],
                history=r['history'])


def stage_noise():
    with Pool(4) as pool:
        out = {w: pool.map(_noise_one, [(w, j) for j in NOISE_PHASES]) for w in ('control', 'main')}
    _write('noise.json', out)


def _noise_value(rows):
    rows = [r for r in rows if r['ok']]
    return float(max(abs(r['lam_radau']-r['lam_dop']) for r in rows)+max(r['radau_residual'] for r in rows))


def _brackets(pts):
    """Sign-change brackets (j, j+1 cyclic) of lambda over converged consecutive phases."""
    out = []
    for k in range(N_PHASE):
        p, n = pts[k], pts[(k+1) % N_PHASE]
        if p['ok'] and n['ok'] and np.sign(p['lam']) != np.sign(n['lam']) and p['lam'] != 0:
            out.append((k, 'up' if n['lam'] > p['lam'] else 'down'))
    return out


def _orbit_one(k):
    pts = _read('scan_main.json')['points']
    p, n = pts[k], pts[(k+1) % N_PHASE]
    w = p['lam']/(p['lam']-n['lam'])
    seeds = [(1-w)*np.array(a)+w*np.array(c) for a, c in zip(p['nodes'], n['nodes'])]
    try:
        o = b.periodic_orbit(P4, seeds, tol=TOL)
        o['ok'] = bool(o['residual'] <= ACCEPT)
        blocks = [b.fd_jacobian(P4, np.array(z)) for z in o['nodes']]
        C = b.centre_block(blocks)
        o['residue_fd'] = float((2-np.trace(C))/4)
        o['constraint_max'] = max(constraint4(z) for z in o['nodes'])
    except (ArithmeticError, ValueError, np.linalg.LinAlgError) as e:
        o = dict(ok=False, error=f'{type(e).__name__}: {e}')
    o['bracket'] = k
    return o


def _same_orbit(A, B, tol=1e-8):
    A, B = np.array(A), np.array(B)
    return any(np.abs(np.roll(A, s, axis=0)-B).max() < tol for s in range(len(A)))


def _jet_residue(nodes, steps):
    blocks = [f.DP(np.r_[z, np.zeros(8)], dims=4, steps=steps)[1] for z in nodes]
    C = b.centre_block(blocks)
    return float((2-np.trace(C))/4), [list(map(float, np.linalg.eigvals(C).real)), list(map(float, np.linalg.eigvals(C).imag))]


def stage_orbits():
    pts = _read('scan_main.json')['points']
    br = _brackets(pts)
    with Pool(4) as pool:
        orbs = pool.map(_orbit_one, [k for k, _ in br])
    for o, (_, d) in zip(orbs, br):
        o['direction'] = d
    distinct = []
    for i, o in enumerate(orbs):
        if o['ok'] and not any(_same_orbit(o['nodes'], orbs[j]['nodes']) for j in distinct):
            distinct.append(i)
    with Pool(4) as pool:
        jets = pool.starmap(_jet_residue, [(orbs[i]['nodes'], s) for i in distinct for s in ((1024, 2048), (2048, 4096))])
    for n, i in enumerate(distinct):
        (r1, _), (r2, mult) = jets[2*n], jets[2*n+1]
        orbs[i].update(residue_jet=r2, residue_jet_coarse=r1, centre_multipliers=mult)
    _write('orbits.json', dict(brackets=[[k, d] for k, d in br], orbits=orbs, distinct=distinct))


def score():
    out = {}
    noise = _read('noise.json')
    for which, q in (('control', 2), ('main', Q_RES)):
        pts = _read(f'scan_{which}.json')['points']
        ok = [p for p in pts if p['ok']]
        nu = _noise_value(noise[which])
        cl = b.classify([p['lam'] for p in ok], nu, q)
        cl.update(converged=len(ok), noise=nu, failed=N_PHASE-len(ok),
                  max_residual=max(p['residual'] for p in ok), max_constraint=max(p['constraint_max'] for p in ok),
                  max_cond=max(p['cond'] for p in ok), max_seed_distance=max(p['seed_distance'] for p in ok))
        if len(ok) < N_PHASE-MAX_FAILED:
            cl['label'] = 'INDETERMINATE'
            cl['reason'] = 'scan failed at too many phases'
        out[which] = cl
    m = out['main']
    if m['label'] == 'BROKEN_CHAIN':
        orb = _read('orbits.json')
        dirs = {orb['orbits'][i]['direction'] for i in orb['distinct']}
        m['orbit_gate'] = bool(len(orb['distinct']) >= 2 and dirs == {'up', 'down'})
        m['distinct_orbits'] = len(orb['distinct'])
        if not m['orbit_gate']:
            m['label'] = 'INDETERMINATE'
            m['reason'] = 'isolated-orbit gate failed'
    _write('result.json', dict(result=out))
    return out


def main(argv):
    stages = argv or ['control', 'main', 'noise', 'orbits', 'score']
    for s in stages:
        t = time.time()
        if s in ('control', 'main'):
            stage_scan(s)
        elif s == 'noise':
            stage_noise()
        elif s == 'orbits':
            if score_main_needs_orbits():
                stage_orbits()
        elif s == 'score':
            print(json.dumps(score(), indent=1))
        print(f'stage {s}: {time.time()-t:.0f} s', flush=True)


def score_main_needs_orbits():
    pts = _read('scan_main.json')['points']
    return len(_brackets(pts)) > 0


if __name__ == '__main__':
    main(sys.argv[1:])
