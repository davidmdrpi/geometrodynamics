"""Authenticated replay of the R3 extension archives (added after the #317 review).

The frozen probe's score_a trusts each circle's stored `ok` flag. This replay
does not. It authenticates the three archive files by SHA-256, validates the
required data, recomputes every acceptance decision from the stored numbers
and from K itself against the registered thresholds, and checks the ladder
schedule. Only then does it re-score with the frozen probe functions. With
--full it also re-evaluates the invariance residual of every accepted circle
on the full system. The frozen probe and archives are unchanged.
"""
import argparse
import hashlib
import json
from multiprocessing import Pool
import numpy as np
from experiments.closure_ledger import r3_extension_probe as probe
from experiments.closure_ledger.esu_floquet_probe import close
from geometrodynamics.waves import r3_return_map as rm

SHA256 = {
    'part_A.json': '3a0cff0ea10342dc660558f006db7a75abd21e6b542697a269a5d0bced8aa547',
    'part_B.json': 'c0e8e04e14f473b9275a9d69a53a570086a6d364a154a6604b0d5b61da16b4a0',
    'result.json': '124f9064019bf7a138c28b1ff805a49e18f4f59db4a58453ea03b4a4d8790386',
}
THRESH = dict(residual=1e-10, fourier_tail=1e-10, constraint_max=1e-10)
REQUIRED_OK = ('a', 'omega', 'action', 'residual', 'fourier_tail', 'constraint_max', 'iterations', 'K')


def _finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and bool(np.isfinite(x))


def _tail(K):
    M = len(K)
    Ck = np.abs(np.fft.fft(K-K.mean(0), axis=0))/M
    kk = np.abs(np.fft.fftfreq(M, 1/M))
    return float(Ck[kk >= M//2-2].max())


def validate_a(A):
    circles = A.get('circles')
    if not isinstance(circles, list) or not circles:
        raise ValueError('part A: no circles')
    if A.get('status') not in ('FAMILY_ENDED', 'CAP_REACHED', 'REACHED_PI'):
        raise ValueError('part A: invalid status')
    a_prev, expect = None, probe.A_START
    for i, c in enumerate(circles):
        if not isinstance(c.get('ok'), bool) or not _finite(c.get('a')):
            raise ValueError(f'circle {i}: missing ok flag or amplitude')
        if abs(c['a']/expect-1) > 1e-12:
            raise ValueError(f'circle {i}: amplitude off the registered ladder')
        if c['ok']:
            for k in REQUIRED_OK:
                if k not in c:
                    raise ValueError(f'circle {i}: missing {k}')
            K = np.asarray(c['K'], dtype=float)
            if K.shape != (probe.GRID, 4) or not np.isfinite(K).all():
                raise ValueError(f'circle {i}: malformed K')
            for k in ('omega', 'action', 'residual', 'fourier_tail', 'constraint_max'):
                if not _finite(c[k]) or (k != 'omega' and c[k] < 0):
                    raise ValueError(f'circle {i}: invalid {k}')
            if any(c[k] > t for k, t in THRESH.items()):
                raise ValueError(f'circle {i}: stored numbers fail the registered acceptance thresholds')
            if abs(abs(rm.action(K))-c['action']) > 1e-12*max(1., c['action']):
                raise ValueError(f'circle {i}: action inconsistent with K')
            if abs(_tail(K)-c['fourier_tail']) > 1e-13:
                raise ValueError(f'circle {i}: Fourier tail inconsistent with K')
            c1 = np.exp(-2j*np.pi*np.arange(probe.GRID)/probe.GRID) @ K[:, 2]/probe.GRID
            if abs(c1.real-c['a']/2) > 1e-9 or abs(c1.imag) > 1e-9:
                raise ValueError(f'circle {i}: amplitude condition violated by K')
            a_prev, expect = c['a'], c['a']*probe.A_FACTOR
        else:
            recomputed_ok = c.get('error') is None and all(
                _finite(c.get(k)) and c[k] <= t for k, t in THRESH.items())
            if recomputed_ok:
                raise ValueError(f'circle {i}: stored failure contradicts its own numbers')
            if a_prev is None:
                raise ValueError(f'circle {i}: first rung failed')
            ratio = c['a']/a_prev
            expect = a_prev*2**.125 if abs(ratio-probe.A_FACTOR) < 1e-12 else None
            if expect is None and i != len(circles)-1:
                raise ValueError(f'circle {i}: continued after a failed retry')
    if A['status'] == 'FAMILY_ENDED' and circles[-1]['ok']:
        raise ValueError('part A: ended on an accepted circle')


def validate_b(B):
    for key in ('normal_form', 'T', 'T_coarse', 'c_coarse', 'constraint_max', 'section_q_max', 'coef_scale'):
        if key not in B:
            raise ValueError('part B: missing '+key)
    for key in ('T', 'T_coarse'):
        for part in ('re', 'im'):
            arr = np.asarray(B[key][part], dtype=float)
            if arr.shape != (5, 5, 5, 5) or not np.isfinite(arr).all():
                raise ValueError(f'part B: malformed {key}')


def invariance_residuals(A, workers=4):
    out = []
    with Pool(workers) as pool:
        for c in A['circles']:
            if not c['ok']:
                continue
            K = np.asarray(c['K'])
            PK = np.array(pool.map(rm.esu_map, list(K)), dtype=object)
            PK = np.array([p[0] for p in PK])
            T, _ = rm.shift_matrix(len(K), c['omega'])
            out.append(float(np.max(abs(PK-T @ K))))
    return out


def replay(full=False, directory=probe.RUN_DIR):
    blobs = {}
    for name, digest in SHA256.items():
        blob = (directory/name).read_bytes()
        if hashlib.sha256(blob).hexdigest() != digest:
            raise ValueError('archive fingerprint mismatch: '+name)
        blobs[name] = json.loads(blob)
    A, B, R = blobs['part_A.json'], blobs['part_B.json'], blobs['result.json']
    src = probe.sources()
    if A['sources'] != src or B['sources'] != src:
        raise ValueError('archive sources do not match the committed code')
    validate_a(A)
    validate_b(B)
    a, b = probe.score_a(A, R['theta0']), probe.score_b(B, R['theta0'])
    if a['label'] != R['A']['label'] or b['label'] != R['B']['label'] or b['checks'] != R['B']['checks']:
        raise ValueError('categorical replay disagreement')
    if not close(json.loads(json.dumps(a)), R['A'], 1e-9) or abs(b['nu_max']-R['B']['nu_max']) > 1e-7:
        raise ValueError('numerical replay disagreement')
    out = dict(replay='VERIFIED', A=a['label'], B=b['label'], nu_max=b['nu_max'],
               accepted_circles=a['n_circles'], a_max=a['a_max'],
               continuation_beyond_a_max='UNRESOLVED (numerical failure, not a physical endpoint)',
               closure_rationals='CANDIDATES ONLY; periodic histories not verified')
    if full:
        res = invariance_residuals(A)
        if max(res) > THRESH['residual']:
            raise ValueError('re-evaluated invariance residual exceeds the registered threshold')
        out['max_reevaluated_invariance_residual'] = max(res)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--full', action='store_true')
    args = ap.parse_args()
    print(json.dumps(replay(args.full), indent=1))


if __name__ == '__main__':
    main()
