"""Authenticated replay of the family/stability archives (added 2026-10-01, #319 review).

The frozen scorer trusts saved residuals and does not check the sample
schedule or the link between samples and stage F. This replay fails closed.
It:
- authenticates the three archives by pinned SHA-256;
- validates finite values and required shapes;
- recomputes the registered 12-point arclength schedule from stage F and
  requires the samples to be exactly those points, in order;
- re-evaluates every stage-F node pair on the full system, checking
  nontriviality and two-node closure against the saved residual;
- checks each sample's saved M2 against the closure geometry: det, block
  decoupling, and a unit multiplier along the node-pair tangent;
- only then re-scores with the frozen probe functions.

With --full it also recomputes every sample's M2 by exact jets and compares.
The frozen probe and archives are unchanged.
"""
import argparse
import hashlib
import json
import numpy as np
from multiprocessing import Pool
from experiments.closure_ledger import r3_family_probe as probe
from geometrodynamics.waves import r3_family as rf

SHA256 = {'stage_F.json': '2f8c1e57aab3ded64fba82c0889f6a09abf60f92b92cad95e0a357adde6a53dd',
          'stage_S.json': '8265820aa6ee2eceb844d53bcaabf23a8a1dc453619062341e06f86770a5b545',
          'result.json': '40070b9854c5ea3f7759c92de0871999a3f3f5da551eef9f62db892a89206cf8'}
ESU = np.r_[1., np.zeros(5)]


def _arr(x, shape, what):
    a = np.asarray(x, dtype=float)
    if a.shape != shape or not np.isfinite(a).all():
        raise ValueError(f'{what}: expected finite array of shape {shape}')
    return a


def _closure(v):
    p0, p1 = probe._P6(v[:6]), probe._P6(v[6:])
    return float(max(np.abs(p0-v[6:]).max(), np.abs(p1-v[:6]).max()))


def schedule(F):
    pts = [np.asarray(F['start']['v'])]+[np.asarray(p['v']) for p in F['points']]
    Z = np.array([p[:6] for p in pts])
    s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(Z, axis=0), axis=1))]
    return pts, [int(np.argmin(abs(s-x))) for x in np.linspace(0, s[-1], probe.N_SAMPLE, endpoint=False)]


def validate_f(F, pool=None):
    if F.get('stop', {}).get('reason') != 'RETURNED_TO_START':
        raise ValueError('stage F did not return to its start')
    recs = [F['start']]+list(F['points'])
    if len(recs) < 20:
        raise ValueError('stage F: too few loop points')
    V = np.array([_arr(r['v'], (12,), 'stage F node pair') for r in recs])
    for r in recs:
        if not np.isfinite(r['residual']) or r['residual'] < 0:
            raise ValueError('stage F: invalid saved residual')
    sep = np.linalg.norm(V[:, :6]-V[:, 6:], axis=1)
    dist = np.minimum(np.linalg.norm(V[:, :6]-ESU, axis=1), np.linalg.norm(V[:, 6:]-ESU, axis=1))
    if sep.min() <= 1e-5 or dist.min() <= 1e-3:
        raise ValueError('stage F: trivial or one-return node pair')
    steps = np.linalg.norm(np.diff(V, axis=0), axis=1)
    if steps.max() > 3*probe.DS or steps.min() < 1e-6:
        raise ValueError('stage F: node pairs are not a continuation sequence')
    if pool is not None:
        re = pool.map(_closure, list(V))
        bad = [i for i, (x, r) in enumerate(zip(re, recs)) if x > 1e-10 or abs(x-r['residual']) > 1e-10]
        if bad:
            raise ValueError(f'stage F: saved residuals not reproduced at {bad[:5]}')
        return float(max(re))
    return None


def validate_s(S, F):
    pts, picks = schedule(F)
    smp = S.get('samples')
    if not isinstance(smp, list) or len(smp) != probe.N_SAMPLE:
        raise ValueError('stage S: wrong sample count')
    if [s.get('index') for s in smp] != picks:
        raise ValueError('stage S: samples are not the registered arclength schedule')
    for s, j in zip(smp, picks):
        v = _arr(s['v'], (12,), 'sample node pair')
        if np.abs(v-pts[j]).max() > 0:
            raise ValueError(f'stage S: sample {j} does not match stage F')
        M2 = _arr(s['M2'], (12, 12), 'M2')
        if abs(np.linalg.det(M2)-1) > 1e-6:
            raise ValueError('stage S: det M2 != 1')
        if max(np.abs(M2[:6, 6:]).max(), np.abs(M2[6:, :6]).max()) > 1e-8*np.abs(M2).max():
            raise ValueError('stage S: diagonal/off-diagonal blocks couple')
        t = pts[(j+1) % len(pts)][:6]-pts[j-1][:6]           # chord tangent of the loop from stage F
        t /= np.linalg.norm(t)
        lam, vec = np.linalg.eig(M2[:6, :6])
        k = int(np.argmin(abs(lam-1)))
        e = np.real(vec[:, k])/np.linalg.norm(np.real(vec[:, k]))
        if abs(lam[k]-1) > 1e-5 or np.sqrt(max(0., 1-(e @ t)**2)) > 1e-2:
            raise ValueError('stage S: unit-multiplier direction of M2 is not the stage-F loop tangent')
        for k in ('closure_jet', 'radau_closure', 'tangent_defect', 'smin_exact', 's2_exact'):
            if not np.isfinite(s[k]) or s[k] < 0:
                raise ValueError('stage S: invalid '+k)
    d = _arr(S.get('direct_perturbation'), (6,), 'direct perturbation')
    if (d < 0).any():
        raise ValueError('stage S: invalid direct-perturbation record')


def _jetM2(v):
    p0, J0 = rf.DP(np.r_[v[:6], np.zeros(6)], 12)
    p1, J1 = rf.DP(np.r_[v[6:], np.zeros(6)], 12)
    return J1 @ J0


def replay(directory=probe.RUN_DIR, full=False):
    blobs = {}
    for name, digest in SHA256.items():
        blob = (directory/name).read_bytes()
        if hashlib.sha256(blob).hexdigest() != digest:
            raise ValueError('archive fingerprint mismatch: '+name)
        blobs[name] = json.loads(blob)
    F, S, R = blobs['stage_F.json'], blobs['stage_S.json'], blobs['result.json']
    src = probe.sources()
    if not (F['sources'] == S['sources'] == R['sources'] == src):
        raise ValueError('archive sources do not match the committed code')
    validate_s(S, F)
    with Pool(4) as pool:
        worst = validate_f(F, pool)
        out = dict(replay='VERIFIED', max_reevaluated_loop_closure=worst)
        if full:
            M2s = pool.map(_jetM2, [np.asarray(s['v']) for s in S['samples']])
            rel = max(np.abs(a-np.asarray(s['M2'])).max()/np.abs(a).max() for a, s in zip(M2s, S['samples']))
            if rel > 1e-8:
                raise ValueError('recomputed M2 disagrees with the archive')
            out['max_recomputed_M2_difference'] = float(rel)
    fresh = probe.score(F, S)
    for k in ('FAMILY', 'ACTION_CONSISTENCY', 'DIAGONAL_TRANSVERSE', 'OFF_DIAGONAL'):
        if fresh[k] != R[k]:
            raise ValueError('categorical replay disagreement: '+k)
    if {k: bool(v) for k, v in fresh['checks'].items()} != R['checks']:
        raise ValueError('check replay disagreement')
    if abs(fresh['loop_action']-R['loop_action']) > 1e-12:
        raise ValueError('loop action disagreement')
    out.update({k: R[k] for k in ('FAMILY', 'ACTION_CONSISTENCY', 'DIAGONAL_TRANSVERSE', 'OFF_DIAGONAL', 'loop_action')})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--full', action='store_true')
    print(json.dumps(replay(full=ap.parse_args().full), indent=1))


if __name__ == '__main__':
    main()
