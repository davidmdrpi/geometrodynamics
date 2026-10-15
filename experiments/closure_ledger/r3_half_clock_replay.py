"""Authenticate and independently reevaluate the half-clock experiment.

    python -m experiments.closure_ledger.r3_half_clock_replay [--full]

--full checks every held-out phase against the full homogeneous matrix
equations with direct full-clock event detection, rather than the producer's
scalar half-clock equations. No output archive is modified.
"""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from geometrodynamics.waves import r3_breaking as b
from experiments.closure_ledger import r3_half_clock_probe as probe
from experiments.closure_ledger import r3_half_clock_validation as validation
from experiments.closure_ledger import r3_half_clock_recovery as recovery
from experiments.closure_ledger import r3_half_clock_validation_parallel as parallel

MANIFEST_SHA256 = 'e63f3717c53c89c264fee817cfb124292ddc86074f83bfcc0a288530a5f7f69d'
FILES = ('known.json', 'scan.json', 'noise.json', 'result.json', 'validation.json', 'execution.json')
FREEZE = '82c832f633a9bada476598df8b4806442d1f3b81'


def load(directory):
    raw = (directory/'manifest.json').read_bytes()
    if hashlib.sha256(raw).hexdigest() != MANIFEST_SHA256:
        raise ValueError('manifest fingerprint mismatch')
    manifest = json.loads(raw)
    if manifest['freeze_commit'] != FREEZE or set(manifest['sha256']) != set(FILES):
        raise ValueError('manifest provenance or inventory mismatch')
    records = {}
    for name in FILES:
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != manifest['sha256'][name]:
            raise ValueError(f'archive fingerprint mismatch: {name}')
        rec = json.loads(raw)
        expected = probe.bindings()
        if name == 'validation.json':
            for path in ('experiments/closure_ledger/r3_half_clock_validation.py',
                         'docs/r3_half_clock_validation_note.md', parallel.PATH):
                expected[path] = hashlib.sha256((probe.ROOT/path).read_bytes()).hexdigest()
            if rec['inputs'] != {n: manifest['sha256'][n] for n in ('scan.json', 'result.json')}:
                raise ValueError('supplemental input hashes differ')
        if name == 'execution.json':
            for path in recovery.EXTRA:
                expected[path] = hashlib.sha256((probe.ROOT/path).read_bytes()).hexdigest()
        if rec['bindings'] != expected:
            raise ValueError(f'source or input fingerprint mismatch: {name}')
        records[name] = rec
    return records


def compare(a, c):
    if isinstance(a, dict):
        if not isinstance(c, dict) or a.keys() != c.keys():
            raise ValueError('score structure differs')
        for key in a:
            compare(a[key], c[key])
    elif isinstance(a, list):
        if not isinstance(c, list) or len(a) != len(c):
            raise ValueError('score length differs')
        for x, y in zip(a, c):
            compare(x, y)
    elif isinstance(a, float):
        if not np.isfinite(a) or not np.isfinite(c) or not np.isclose(a, c, rtol=1e-10, atol=1e-18):
            raise ValueError('score numeric value differs')
    elif a != c:
        raise ValueError('score discrete value differs')


def replay(directory=None, full=False):
    directory = Path(directory or probe.DIRECTORY)
    recs = load(directory)
    points = recs['scan.json']['points']
    recovery.verify_identical_body()
    parallel.verify_identical_body()
    execution = recs['execution.json']
    if (execution['completed_phase_count'] != probe.N_PHASE
            or execution['numerical_body_identical'] is not True):
        raise ValueError('incomplete execution record')
    # Every omitted checkpoint is exactly reconstructible from the aggregate.
    reconstructed = {
        f'phase_{p["j"]:02d}.json': hashlib.sha256(json.dumps(
            dict(bindings=execution['bindings'], point=p), indent=1, allow_nan=False).encode()).hexdigest()
        for p in points}
    if reconstructed != execution['checkpoint_sha256']:
        raise ValueError('aggregate differs from saved phase checkpoints')
    if len(points) != probe.N_PHASE:
        raise ValueError('incomplete scan')
    for j, point in enumerate(points):
        if point['j'] != j or abs(point['phase']-2*np.pi*j/probe.N_PHASE) > 1e-15:
            raise ValueError('wrong phase grid or ordering')
        if not point['ok']:
            continue
        for key, shape in (('nodes', (7, 4)), ('J', (29, 29)), ('g', (4,)),
                           ('c0', (4,)), ('t0', (4,))):
            val = np.asarray(point[key])
            if val.shape != shape or not np.isfinite(val).all():
                raise ValueError(f'malformed {key} at phase {j}')
        if point['residual'] > probe.ACCEPT or not np.isfinite(point['lam']):
            raise ValueError('invalid accepted point')
        if abs(np.linalg.norm(point['t0'])-1) > 1e-12:
            raise ValueError('unnormalized phase tangent')
        if np.max(abs(np.array(point['g'])-b.symplectic(4).T @ point['t0'])) > 1e-12:
            raise ValueError('incorrect obstruction direction')
    calculated = probe.score(points, recs['noise.json']['rows'], recs['known.json'])
    archived = {k: v for k, v in recs['result.json'].items() if k not in ('bindings', 'written_utc')}
    compare(calculated, archived)
    supplement = recs['validation.json']
    for row in supplement['rows']:
        if row['ok'] and row['primary'] != points[row['j']]['lam']:
            raise ValueError('supplemental primary value differs')
    compare(validation.summarize(supplement['rows'], archived), supplement['result'])
    out = dict(replay='VERIFIED', freeze_commit=FREEZE, label=calculated['label'],
               independently_confirmed=supplement['result']['confirmed'])
    if full:
        worst = 0.
        for point in points:
            if not point['ok']:
                raise ValueError('cannot fully reevaluate an incomplete scan')
            Z = np.array(point['nodes'])
            residual = b._residual([probe.full_matrix_map(z) for z in Z], Z, point['lam'],
                                  np.array(point['g']), np.array(point['c0']), np.array(point['t0']))
            worst = max(worst, float(np.max(abs(residual))))
            if not np.isfinite(worst) or worst > 1e-10:
                raise ValueError(f'full-equation reclosure failed at phase {point["j"]}: {worst}')
        out.update(full_equation_phases=len(points), worst_full_equation_residual=worst)
    return out


if __name__ == '__main__':
    print(json.dumps(replay(full='--full' in sys.argv), indent=1))
