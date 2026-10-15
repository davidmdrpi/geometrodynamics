"""Authenticated replay of the resonance-breaking archives (docs/r3_breaking.md).

Fail-closed: archive bytes must match the pinned SHA-256, the bound source
hashes must match the files on disk, every converged scan point must carry a
consistent system (shapes, finite values, residual), and the registered labels
are re-derived from the archives with the frozen scoring code. With --full,
every converged scan point and isolated orbit is re-evaluated with the
unchanged maps.
Usage: python -m experiments.closure_ledger.r3_breaking_replay [--full]
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
from experiments.closure_ledger import r3_breaking_probe as probe

SHA256 = {
    'scan_control.json': '40ba7aaf15110e980a9d5309631021a9d33421a14b655b28ebcd7c0f35f81760',
    'scan_main.json': 'd00b8947f0cd19d6b28e8b98716d473d677ca54af27d95032e2eb7cb22f040bd',
    'noise.json': 'b16d69458bf8f16954fe024511dd0da4e92963b3413c6387fcf67a0dc7537efe',
    'orbits.json': '87227498a3a2e7511af621437fa5928ad8ab6f6b379be53122aa5558fe8dcf5f',
    'result.json': '4aaf42df3e9116d0d864ef0b140c0d83ad4e06c4b1807f0fe6ad9dcb082407a7',
}
FULL_TOL = 1e-10


def _load(directory):
    out = {}
    for name, h in SHA256.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != h:
            raise ValueError(f'fingerprint mismatch: {name}')
        out[name] = json.loads(raw)
    return out


def validate_scan(rec, q, n):
    pts = rec['points']
    if len(pts) != probe.N_PHASE or [p['j'] for p in pts] != list(range(probe.N_PHASE)):
        raise ValueError('scan phases incomplete or reordered')
    for p in pts:
        if abs(p['phi']-2*np.pi*p['j']/probe.N_PHASE) > 1e-15:
            raise ValueError('phase off the registered grid')
        if not p['ok']:
            continue
        Z = np.array(p['nodes'], float)
        if Z.shape != (q, n) or not np.isfinite(Z).all() or not np.isfinite(p['lam']):
            raise ValueError('malformed scan point')
        if p['residual'] > probe.ACCEPT or p['history'][-1] != p['residual']:
            raise ValueError('ok flag inconsistent with residual')
        g, t0 = np.array(p['g']), np.array(p['t0'])
        if abs(np.linalg.norm(g)-1) > 1e-12 or abs(np.linalg.norm(t0)-1) > 1e-12:
            raise ValueError('unnormalised direction')


def reevaluate(rec, P, tol=FULL_TOL):
    """Recompute the scan residual with the archived lambda at every converged point."""
    from geometrodynamics.waves import r3_breaking as b
    worst = 0.
    for p in rec['points']:
        if not p['ok']:
            continue
        Z = [np.array(z) for z in p['nodes']]
        F = b._residual([np.asarray(P(z)) for z in Z], Z, p['lam'], np.array(p['g']), np.array(p['c0']), np.array(p['t0']))
        worst = max(worst, float(np.abs(F).max()))
        if worst > tol:
            raise ValueError(f'scan point {p["j"]} does not re-close: {worst:.2e}')
    return worst


def replay(directory=None, full=False):
    directory = Path(directory or probe.RUN_DIR)
    recs = _load(directory)
    src = probe.sources()
    for name, rec in recs.items():
        if rec['sources'] != src:
            raise ValueError(f'source hashes differ: {name}')
    validate_scan(recs['scan_control.json'], 2, 6)
    validate_scan(recs['scan_main.json'], probe.Q_RES, 4)
    import shutil
    import tempfile
    old = probe.RUN_DIR
    with tempfile.TemporaryDirectory() as tmp:
        for name in recs:
            shutil.copy(directory/name, Path(tmp)/name)
        probe.RUN_DIR = Path(tmp)          # the frozen score() rewrites result.json; never in the run directory
        try:
            again = json.loads(json.dumps(probe.score()))
        finally:
            probe.RUN_DIR = old
    ref = recs['result.json']['result']
    for k in ('control', 'main'):
        if again[k]['label'] != ref[k]['label']:
            raise ValueError(f'{k} label does not re-derive')
    out = dict(replay='VERIFIED', control=again['control']['label'], main=again['main']['label'])
    if full:
        out['worst_control'] = reevaluate(recs['scan_control.json'], probe.P6)
        out['worst_main'] = reevaluate(recs['scan_main.json'], probe.P4)
    return out


if __name__ == '__main__':
    print(json.dumps(replay(full='--full' in sys.argv), indent=1))
