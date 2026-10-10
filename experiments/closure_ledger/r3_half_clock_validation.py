"""Supplemental checks specified during production; see dated validation note."""
from datetime import datetime, timezone
import hashlib
import json

import numpy as np

from geometrodynamics.waves import r3_breaking as b
from experiments.closure_ledger import r3_half_clock_probe as probe

RADAU_PHASES = (0, 14, 28, 42, 56, 70, 83)


def summarize(rows, primary):
    expected = [('DOP853', j) for j in range(84)]+[('Radau', j) for j in RADAU_PHASES]
    if [(r['method'], r['j']) for r in rows] != expected or not all(r['ok'] for r in rows):
        return dict(confirmed=False, reason='incomplete supplemental checks')
    nu = max(abs(r['lam']-r['primary'])+r['residual'] for r in rows)
    resolution = max(primary['resolution'], 10*nu)
    alt = np.array([r['lam'] for r in rows[:84]])
    ref = np.array([r['primary'] for r in rows[:84]])
    spec = abs(np.fft.rfft(alt))/84
    dominant = int(np.argmax(spec[1:22])+1)
    original_retained = bool(primary['label'] == 'SEVENTH_HARMONIC_PREDICTION_SUPPORTED'
                             and nu <= 1e-10 and primary['harmonics'][7] >= 10*resolution)
    alternative_passes = bool(dominant == 7 and spec[7] >= 10*resolution and b.sign_changes(alt) == 14)
    return dict(noise=nu, combined_resolution=resolution, original_retained=original_retained,
                alternative_passes=alternative_passes, confirmed=original_retained and alternative_passes,
                alternative_harmonics=spec.tolist(), alternative_dominant=dominant,
                alternative_sign_changes=b.sign_changes(alt), correlation=float(np.corrcoef(ref, alt)[0, 1]),
                max_difference=float(max(abs(alt-ref))))


def main():
    output = probe.DIRECTORY/'validation.json'
    if output.exists():
        raise FileExistsError('supplement already exists')
    points = json.loads((probe.DIRECTORY/'scan.json').read_text())['points']
    result = json.loads((probe.DIRECTORY/'result.json').read_text())
    rows = []
    for method, indices, mapping in (('DOP853', range(84), probe.full_matrix_map),
                                      ('Radau', RADAU_PHASES, probe.radau_matrix_map)):
        for j in indices:
            try:
                if not points[j]['ok']:
                    raise ValueError('primary point failed')
                rec = b.chord_resolve(mapping, points[j], tol=probe.TOL, maxit=6)
                row = dict(method=method, j=j, primary=points[j]['lam'], ok=True, **rec)
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
                row = dict(method=method, j=j, ok=False, error=f'{type(error).__name__}: {error}')
            rows.append(row)
            print(f'supplement {method} phase {j}: ok={row["ok"]}', flush=True)
    paths = ('experiments/closure_ledger/r3_half_clock_validation.py',
             'docs/r3_half_clock_validation_note.md')
    hashes = probe.bindings()
    hashes.update({p: hashlib.sha256((probe.ROOT/p).read_bytes()).hexdigest() for p in paths})
    inputs = {n: hashlib.sha256((probe.DIRECTORY/n).read_bytes()).hexdigest()
              for n in ('scan.json', 'result.json')}
    out = dict(bindings=hashes, inputs=inputs, written_utc=datetime.now(timezone.utc).isoformat(),
               rows=rows, result=summarize(rows, result))
    with output.open('x') as stream:
        json.dump(out, stream, indent=1, allow_nan=False)
    print(json.dumps(out['result'], indent=1))


if __name__ == '__main__':
    main()
