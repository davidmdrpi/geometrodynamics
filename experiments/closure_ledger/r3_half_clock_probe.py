"""Frozen half-clock explanation and held-out LRS 3/7 experiment.

Run after publishing the protocol and this producer:
  OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.r3_half_clock_probe
No archived input or historical source is modified. Output creation is exclusive.
"""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import time

import numpy as np

from geometrodynamics.waves import r3_breaking as b
from geometrodynamics.waves import r3_half_clock as h
from geometrodynamics.waves import r3_return_map as rm

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT/'experiments/closure_ledger/runs/20261010_r3_half_clock'
EXT = 'experiments/closure_ledger/runs/20260929_r3_extension/part_A.json'
KNOWN = 'experiments/closure_ledger/runs/20261005_r3_breaking/scan_main.json'
SOURCES = (
    'docs/r3_half_clock_prereg.md',
    'geometrodynamics/waves/r3_half_clock.py',
    'experiments/closure_ledger/r3_half_clock_probe.py',
    'geometrodynamics/waves/r3_breaking.py',
    'geometrodynamics/waves/r3_return_map.py',
    'geometrodynamics/waves/nonlinear_supported_tt.py',
    'geometrodynamics/waves/jets.py',
    'tests/test_r3_half_clock.py',
)
N_PHASE = 84
P_RES, Q_RES = 3, 7
NOISE_PHASES = (1, 13, 25, 37, 49, 61, 73)
TOL, ACCEPT = 2e-13, 2e-12


def bindings():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in (*SOURCES, EXT, KNOWN)}


def write(name, rec):
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    with (DIRECTORY/name).open('x') as stream:
        json.dump(dict(bindings=bindings(), written_utc=datetime.now(timezone.utc).isoformat(),
                       **rec), stream, indent=1, allow_nan=False)


def curve():
    circles = [c for c in json.loads((ROOT/EXT).read_text())['circles'] if c.get('ok')]
    target = 2*np.pi*P_RES/Q_RES
    hi = min((c for c in circles if c['omega'] > target), key=lambda c: c['omega'])
    lo = max((c for c in circles if c['omega'] < target), key=lambda c: c['omega'])
    s = (hi['omega']-target)/(hi['omega']-lo['omega'])
    K = (1-s)*np.array(hi['K'])+s*np.array(lo['K'])
    return b.trig_curve(K), dict(amplitude_bracket=[hi['a'], lo['a']], weight=float(s))


def full_matrix_map(z):
    return rm.esu_map(np.asarray(z), tol=(1e-13, 1e-15))[0]


def radau_matrix_map(z):
    return rm.esu_map(np.asarray(z), method='Radau', tol=(1e-12, 1e-14))[0]


def validate_known():
    """Retrospective explanation only; no held-out inference from these cases."""
    rows = []
    circles = [c for c in json.loads((ROOT/EXT).read_text())['circles'] if c.get('ok')]
    for index in (0, 20, 24):
        cc = circles[index]
        c, _ = b.trig_curve(cc['K'])
        for phase in (.37, 1.21, 2.43):
            z = c(phase)
            hz = h.half_map(z, full=True)
            hh = h.half_map(hz['z'])
            matrix, constraint, time_full = rm.esu_map(z, tol=(1e-13, 1e-15))
            rows.append(dict(a=cc['a'], phase=phase,
                square_error=float(np.max(abs(hh-matrix))),
                correct_lift_error=float(np.max(abs(hz['z']-c(phase+cc['omega']/2+np.pi)))),
                wrong_lift_error=float(np.max(abs(hz['z']-c(phase+cc['omega']/2)))),
                constraint=max(hz['constraint'], constraint)))
    points = json.loads((ROOT/KNOWN).read_text())['points']
    replay = []
    for j in range(0, 60, 5):
        pt = points[j]
        rr = b.chord_resolve(h.squared_map, pt, tol=TOL, maxit=6)
        replay.append(dict(j=j, old_lambda=pt['lam'], **rr))
    return dict(identity=rows, known_2_5=replay)


def scan():
    (c, dc), bracket = curve()
    rows = []
    for j in range(N_PHASE):
        phase = 2*np.pi*j/N_PHASE
        start = time.time()
        try:
            seeds = [c(phase+2*np.pi*P_RES*i/Q_RES) for i in range(Q_RES)]
            row = b.scan_point(h.squared_map, seeds, c(phase), dc(phase),
                               tol=TOL, maxit=12)
            # The frozen scanner's final history can precede its final update;
            # always recompute the residual at the actual returned nodes here.
            Z = np.array(row['nodes'])
            residual = b._residual([h.squared_map(z) for z in Z], Z, row['lam'],
                                  np.array(row['g']), np.array(row['c0']), np.array(row['t0']))
            row['residual'] = float(np.max(abs(residual)))
            row['ok'] = bool(row['residual'] <= ACCEPT)
            row['constraint'] = max(h.half_map(z, full=True)['constraint'] for z in Z)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            row = dict(ok=False, error=f'{type(error).__name__}: {error}')
        row.update(j=j, phase=phase, seconds=time.time()-start)
        rows.append(row)
        print(f'3/7 phase {j+1}/{N_PHASE}: ok={row["ok"]}', flush=True)
    return dict(bracket=bracket, points=rows)


def noise(points):
    rows = []
    for j in NOISE_PHASES:
        if not points[j]['ok']:
            rows.append(dict(j=j, ok=False))
            continue
        alternatives = {}
        try:
            for label, mapping in (('matrix_DOP853', full_matrix_map), ('matrix_Radau', radau_matrix_map)):
                rr = b.chord_resolve(mapping, points[j], tol=TOL, maxit=6)
                alternatives[label] = rr
            rows.append(dict(j=j, ok=True, primary=points[j]['lam'], alternatives=alternatives))
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            rows.append(dict(j=j, ok=False, error=f'{type(error).__name__}: {error}'))
        print(f'independent noise phase {j}: ok={rows[-1]["ok"]}', flush=True)
    return rows


def score(points, noise_rows, known):
    """Finite-resolution prediction, not a theorem that allowed coefficients are nonzero."""
    if (len(points) != N_PHASE or [p['j'] for p in points] != list(range(N_PHASE))
            or not all(p['ok'] for p in points)
            or [r['j'] for r in noise_rows] != list(NOISE_PHASES)
            or not all(r['ok'] for r in noise_rows)):
        return dict(label='NUMERICALLY_UNRESOLVED', reason='incomplete scan or noise checks')
    lam = np.array([p['lam'] for p in points])
    if not np.isfinite(lam).all():
        return dict(label='NUMERICALLY_UNRESOLVED', reason='nonfinite lambda')
    nu = max(abs(a['lam']-r['primary'])+a['residual'] for r in noise_rows
             for a in r['alternatives'].values())
    radius = max(10*nu, 1e-11)
    spec = abs(np.fft.rfft(lam))/N_PHASE  # complex Fourier coefficient magnitude, NOT doubled
    dominant = int(np.argmax(spec[1:22])+1)
    numerical = (nu <= 1e-10 and max(p['residual'] for p in points) <= ACCEPT
                 and max(p['constraint'] for p in points) <= 1e-10
                 and max(p['cond'] for p in points) <= 1e8)
    identity = (max(r['square_error'] for r in known['identity']) <= 1e-9
                and max(r['correct_lift_error'] for r in known['identity']) <= 1e-8
                and min(r['wrong_lift_error'] for r in known['identity']) >= 1e-3)
    predicted = dominant == 7 and spec[7] >= 10*radius and b.sign_changes(lam) == 14
    label = ('NUMERICALLY_UNRESOLVED' if not numerical else
             'HALF_CLOCK_IDENTITY_FAILED' if not identity else
             'SEVENTH_HARMONIC_PREDICTION_SUPPORTED' if predicted else
             'SEVENTH_HARMONIC_PREDICTION_FAILED')
    return dict(label=label, noise=nu, resolution=radius, harmonics=spec.tolist(),
                dominant_harmonic=dominant, sign_changes=b.sign_changes(lam),
                max_lambda=float(max(abs(lam))), numerical_gate=bool(numerical),
                identity_gate=bool(identity), prediction_gate=bool(predicted))


def main():
    if DIRECTORY.exists():
        raise FileExistsError('run directory already exists; no silent rerun')
    known = validate_known()
    write('known.json', known)
    production = scan()
    write('scan.json', production)
    checks = noise(production['points'])
    write('noise.json', dict(rows=checks))
    result = score(production['points'], checks, known)
    write('result.json', result)
    print(json.dumps(result, indent=1))


if __name__ == '__main__':
    main()
