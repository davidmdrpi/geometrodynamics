"""Checkpointed scheduling of the already-specified independent checks.

Written after observing the primary 3/7 spectrum, to run the pre-specified
checks concurrently while the registered noise stage runs. No numerical
parameter, phase, integrator, or scoring condition changes. The numerical
body is checked structurally against the previously published implementation.
"""
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import inspect
import json
import numpy as np
from geometrodynamics.waves import r3_breaking as b
from experiments.closure_ledger import r3_half_clock_probe as probe
from experiments.closure_ledger import r3_half_clock_validation as validation

PATH = 'experiments/closure_ledger/r3_half_clock_validation_parallel.py'
CHECKPOINTS = probe.ROOT.parent/'half-clock-validation-checkpoints'


def check(job):
    method, j = job
    mapping = probe.full_matrix_map if method == 'DOP853' else probe.radau_matrix_map
    points = json.loads((probe.DIRECTORY/'scan.json').read_text())['points']
    try:
        if not points[j]['ok']:
            raise ValueError('primary point failed')
        rec = b.chord_resolve(mapping, points[j], tol=probe.TOL, maxit=6)
        row = dict(method=method, j=j, primary=points[j]['lam'], ok=True, **rec)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        row = dict(method=method, j=j, ok=False, error=f'{type(error).__name__}: {error}')
    return row


def verify_identical_body():
    frozen = ast.parse(inspect.getsource(validation.main)).body[0]
    outer = next(n for n in frozen.body if isinstance(n, ast.For))
    inner = next(n for n in outer.body if isinstance(n, ast.For))
    current = ast.parse(inspect.getsource(check)).body[0]
    normalize = lambda nodes: ast.dump(ast.Module(body=nodes, type_ignores=[]))
    if normalize(inner.body[:-2]) != normalize(current.body[3:-1]):
        raise ValueError('independent-check numerical body differs')


def main():
    verify_identical_body()
    scan = json.loads((probe.DIRECTORY/'scan.json').read_text())
    if scan['bindings'] != probe.bindings():
        raise ValueError('primary source/input differs')
    bound = probe.bindings()
    for path in ('experiments/closure_ledger/r3_half_clock_validation.py',
                 'docs/r3_half_clock_validation_note.md', PATH):
        bound[path] = hashlib.sha256((probe.ROOT/path).read_bytes()).hexdigest()
    scan_hash = hashlib.sha256((probe.DIRECTORY/'scan.json').read_bytes()).hexdigest()
    jobs = [('DOP853', j) for j in range(84)]+[('Radau', j) for j in validation.RADAU_PHASES]
    CHECKPOINTS.mkdir(exist_ok=True)
    rows = {}
    for file in sorted(CHECKPOINTS.glob('*.json')):
        rec = json.loads(file.read_text())
        if rec['bindings'] != bound or rec['scan_sha256'] != scan_hash:
            raise ValueError('supplement checkpoint provenance differs')
        row = rec['row']; key = (row['method'], row['j'])
        if key not in jobs or key in rows:
            raise ValueError('supplement checkpoint index differs')
        rows[key] = row
    with ProcessPoolExecutor(max_workers=4) as pool:
        pending = {pool.submit(check, job): job for job in jobs if job not in rows}
        for future in as_completed(pending):
            row = future.result(); key = (row['method'], row['j'])
            with (CHECKPOINTS/f'{key[0]}_{key[1]:02d}.json').open('x') as stream:
                json.dump(dict(bindings=bound, scan_sha256=scan_hash, row=row), stream,
                          indent=1, allow_nan=False)
            rows[key] = row
            print(f'independent checkpoint {len(rows)}/91, {key}, ok={row["ok"]}', flush=True)
    if not (probe.DIRECTORY/'result.json').exists():
        print('All checks saved; awaiting registered result for scoring', flush=True)
        return
    ordered = [rows[job] for job in jobs]
    result = json.loads((probe.DIRECTORY/'result.json').read_text())
    inputs = {n: hashlib.sha256((probe.DIRECTORY/n).read_bytes()).hexdigest()
              for n in ('scan.json', 'result.json')}
    output = dict(bindings=bound, inputs=inputs, written_utc=datetime.now(timezone.utc).isoformat(),
                  execution='four workers; exact frozen numerical body; checkpointed',
                  rows=ordered, result=validation.summarize(ordered, result))
    with (probe.DIRECTORY/'validation.json').open('x') as stream:
        json.dump(output, stream, indent=1, allow_nan=False)
    print(json.dumps(output['result'], indent=1))


if __name__ == '__main__':
    main()
