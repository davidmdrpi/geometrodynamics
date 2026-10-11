"""Checkpointed scheduling of the unchanged frozen per-phase calculation."""
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import inspect
import json
from pathlib import Path
import time
import numpy as np
from geometrodynamics.waves import r3_breaking as b, r3_half_clock as h
from experiments.closure_ledger import r3_half_clock_probe as probe
from experiments.closure_ledger.r3_half_clock_probe import N_PHASE, P_RES, Q_RES, TOL, ACCEPT

CHECKPOINTS = probe.ROOT.parent/'half-clock-checkpoints'
EXTRA = ('docs/r3_half_clock_recovery_note.md', 'experiments/closure_ledger/r3_half_clock_recovery.py')


def point(j):
    (c, dc), _ = probe.curve()
    phase = 2 * np.pi * j / N_PHASE
    start = time.time()
    try:
        seeds = [c(phase + 2 * np.pi * P_RES * i / Q_RES) for i in range(Q_RES)]
        row = b.scan_point(h.squared_map, seeds, c(phase), dc(phase), tol=TOL, maxit=12)
        Z = np.array(row['nodes'])
        residual = b._residual([h.squared_map(z) for z in Z], Z, row['lam'], np.array(row['g']), np.array(row['c0']), np.array(row['t0']))
        row['residual'] = float(np.max(abs(residual)))
        row['ok'] = bool(row['residual'] <= ACCEPT)
        row['constraint'] = max((h.half_map(z, full=True)['constraint'] for z in Z))
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        row = dict(ok=False, error=f'{type(error).__name__}: {error}')
    row.update(j=j, phase=phase, seconds=time.time() - start)
    return row


def verify_identical_body():
    frozen = ast.parse(inspect.getsource(probe.scan)).body[0]
    loop = next(n for n in frozen.body if isinstance(n, ast.For))
    recovered = ast.parse(inspect.getsource(point)).body[0]
    normalize = lambda nodes: ast.dump(ast.Module(body=nodes, type_ignores=[]))
    if normalize(loop.body[:-2]) != normalize(recovered.body[1:-1]):
        raise ValueError('per-phase numerical code differs from frozen producer')


def main():
    verify_identical_body()
    known = json.loads((probe.DIRECTORY/'known.json').read_text())
    if known['bindings'] != probe.bindings():
        raise ValueError('frozen source/input changed')
    bound = probe.bindings()
    bound.update({p: hashlib.sha256((probe.ROOT/p).read_bytes()).hexdigest() for p in EXTRA})
    CHECKPOINTS.mkdir(exist_ok=True)
    rows = {}
    for file in sorted(CHECKPOINTS.glob('phase_*.json')):
        rec = json.loads(file.read_text())
        if rec['bindings'] != bound:
            raise ValueError('checkpoint binding differs')
        j = rec['point']['j']
        if file.name != f'phase_{j:02d}.json' or j in rows or j not in range(N_PHASE):
            raise ValueError('checkpoint index differs')
        rows[j] = rec['point']
    started = datetime.now(timezone.utc).isoformat()
    if not (probe.DIRECTORY/'scan.json').exists():
        with ProcessPoolExecutor(max_workers=4) as pool:
            pending = {pool.submit(point, j): j for j in range(N_PHASE) if j not in rows}
            for future in as_completed(pending):
                row = future.result()
                j = row['j']
                with (CHECKPOINTS/f'phase_{j:02d}.json').open('x') as stream:
                    json.dump(dict(bindings=bound, point=row), stream, indent=1, allow_nan=False)
                rows[j] = row
                print(f'checkpoint {len(rows)}/{N_PHASE}, phase {j}, ok={row["ok"]}', flush=True)
        _, bracket = probe.curve()
        probe.write('scan.json', dict(bracket=bracket, points=[rows[j] for j in range(N_PHASE)]))
    scan = json.loads((probe.DIRECTORY/'scan.json').read_text())
    if scan['bindings'] != probe.bindings():
        raise ValueError('scan binding differs')
    if not (probe.DIRECTORY/'noise.json').exists():
        checks = probe.noise(scan['points'])
        probe.write('noise.json', dict(rows=checks))
    checks = json.loads((probe.DIRECTORY/'noise.json').read_text())
    if checks['bindings'] != probe.bindings():
        raise ValueError('noise binding differs')
    result = probe.score(scan['points'], checks['rows'], known)
    if not (probe.DIRECTORY/'result.json').exists():
        probe.write('result.json', result)
    execution = dict(bindings=bound, started_utc=started,
                     completed_utc=datetime.now(timezone.utc).isoformat(),
                     numerical_body_identical=True, workers=4,
                     prior_attempt='Interrupted after 63 printed convergence flags; scan not persisted',
                     completed_phase_count=len(scan['points']),
                     checkpoint_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                        for p in sorted(CHECKPOINTS.glob('phase_*.json'))})
    with (probe.DIRECTORY/'execution.json').open('x') as stream:
        json.dump(execution, stream, indent=1, allow_nan=False)
    print(json.dumps(result, indent=1))


if __name__ == '__main__':
    main()
