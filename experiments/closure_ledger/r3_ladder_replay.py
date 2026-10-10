"""Authenticated replay of the ladder archives (docs/r3_ladder_prereg.md).

Fail-closed: archive bytes must match the pinned SHA-256, the bound source hashes must match the files
on disk, every usable high-precision row must be well formed, and the registered labels are re-derived
with the frozen scoring code in a temporary directory. With --full, the high-precision solution at
every 10th phase of every rung is re-closed with the unchanged 160-bit map.
Usage: python -m experiments.closure_ledger.r3_ladder_replay [--full]
"""
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path
import gmpy2
from gmpy2 import mpfr
import numpy as np
from geometrodynamics.waves import lrs_taylor as lt
from experiments.closure_ledger import r3_ladder_probe as probe

SHA256 = {}   # pinned after the run
FULL_TOL = 1e-33
LABELS = ('primary', 'signal_2_5', 'harmonic_selection', 'exact_integrability')


def names():
    out = ['started.json', 'result.json']
    for p, q in probe.RUNGS:
        out += [f'{s}_{probe.tag(p, q)}.json' for s in ('scan', 'hp', 'hpnoise')]
    return out


def _load(directory):
    if set(SHA256) != set(names()):
        raise ValueError('pinned archive list incomplete')
    out = {}
    for name, h in SHA256.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != h:
            raise ValueError(f'fingerprint mismatch: {name}')
        out[name] = json.loads(raw)
    return out


def validate_rows(scan, rows, tol):
    if len(rows) != probe.N_PHASE or [r['j'] for r in rows] != list(range(probe.N_PHASE)):
        raise ValueError('rows incomplete or reordered')
    q = scan['q']
    for r, pt in zip(rows, scan['points']):
        if not r['ok']:
            continue
        if not pt['ok']:
            raise ValueError('high-precision row without a converged scan point')
        if len(r['nodes']) != q or any(len(z) != 4 for z in r['nodes']):
            raise ValueError('malformed node set')
        if not (r['residual'] < tol and r['history'][-1] == r['residual']):
            raise ValueError('ok flag inconsistent with residual')
        if abs(float(mpfr(r['lam']))-r['lam_float']) > 1e-15*max(1., abs(r['lam_float'])):
            raise ValueError('lambda string and float disagree')


def reclose(scan, rows, every=10):
    worst = 0.
    P = probe.hp_P(lt.C1)
    with gmpy2.context(gmpy2.get_context(), precision=lt.C1['bits']):
        for j in range(0, probe.N_PHASE, every):
            r, pt = rows[j], scan['points'][j]
            if not r['ok']:
                continue
            Z = [[mpfr(v) for v in z] for z in r['nodes']]
            g, t0, c0 = ([mpfr(v) for v in pt[k]] for k in ('g', 't0', 'c0'))
            F = probe.hp_residual(P, Z, mpfr(r['lam']), g, t0, c0)
            worst = max(worst, float(max(abs(f) for f in F)))
            if worst > FULL_TOL:
                raise ValueError(f'row {j} does not re-close: {worst:.2e}')
    return worst


def replay(directory=None, full=False):
    directory = Path(directory or probe.RUN_DIR)
    recs = _load(directory)
    src = probe.sources()
    for name, rec in recs.items():
        if rec['sources'] != src:
            raise ValueError(f'source hashes differ: {name}')
    for p, q in probe.RUNGS:
        t = probe.tag(p, q)
        validate_rows(recs[f'scan_{t}.json'], recs[f'hp_{t}.json']['rows'], probe.HP_TOL['C1'])
        validate_rows(recs[f'scan_{t}.json'], recs[f'hpnoise_{t}.json']['rows'], probe.HP_TOL['C2'])
    old = probe.RUN_DIR
    with tempfile.TemporaryDirectory() as tmp:
        for name in recs:
            if name != 'result.json':
                shutil.copy(directory/name, Path(tmp)/name)
        probe.RUN_DIR = Path(tmp)
        try:
            again = json.loads(json.dumps(probe.score()))
        finally:
            probe.RUN_DIR = old
    ref = recs['result.json']['result']
    for k in LABELS:
        if again[k] != ref[k]:
            raise ValueError(f'{k} does not re-derive')
    if again['exponent_fit']['label'] != ref['exponent_fit']['label']:
        raise ValueError('exponent_fit does not re-derive')
    out = dict(replay='VERIFIED', **{k: again[k] for k in LABELS}, exponent_fit=again['exponent_fit']['label'])
    if full:
        out['worst_reclosure'] = max(reclose(recs[f'scan_{probe.tag(p, q)}.json'], recs[f'hp_{probe.tag(p, q)}.json']['rows'])
                                     for p, q in probe.RUNGS)
    return out


if __name__ == '__main__':
    print(json.dumps(replay(full='--full' in sys.argv), indent=1))
