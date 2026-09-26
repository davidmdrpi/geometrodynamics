"""Prospectively frozen coordinate error budget; retains historical failures."""
import argparse
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import mpmath as mp
from geometrodynamics.waves import coordinate_budget as c
from . import localized_mouth_refinement_probe as refinement
from .evidence_archive import read_bytes

FREEZE = 'f9025668b2d6ad4bd94a2f704aa699e8849831a7'
REFINEMENT_HASH = 'b9fd7d2e652ccdae2012c636c7fc249fc8b5cca05cb872df5b32576f591cb3ae'
ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('geometrodynamics/waves/coordinate_budget.py',
           'experiments/closure_ledger/localized_mouth_error_budget_probe.py',
           'docs/localized_mouth_error_budget_prereg.md')
FALSE = dict(FOUR_SCALAR_HANDLE_CONSTRAINT_DATA=False, LOCALIZED_BULK_MOUTH_INITIAL_DATA=False)


def packed(value):
    if isinstance(value, dict):
        return {k: packed(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [packed(v) for v in value]
    return mp.nstr(value, 65)


def agreement(a, b):
    """Strict structure, finite numeric evidence, scaled 1e-40 agreement."""
    if isinstance(b, dict):
        return isinstance(a, dict) and a.keys() == b.keys() and all(agreement(a[k], b[k]) for k in b)
    if isinstance(b, list):
        return isinstance(a, list) and len(a) == len(b) and all(agreement(x, y) for x, y in zip(a, b))
    try:
        x, y = mp.mpf(a), mp.mpf(b)
        return bool(mp.isfinite(x) and mp.isfinite(y) and abs(x-y) <= mp.mpf('1e-40')*max(1, abs(y)))
    except (ValueError, TypeError):
        return False


def load_inputs(directory):
    original = refinement.read_original(directory/'probe.json.gz')
    raw = gzip.decompress(read_bytes(directory/'refinement.json.gz'))
    if hashlib.sha256(raw).hexdigest() != REFINEMENT_HASH:
        raise ValueError('refinement archive hash mismatch')
    refined = json.loads(raw)
    historical = refinement.score(refined, original)
    if historical['gates'] != {k: k != 'physical' for k in refinement.p.GATES}:
        raise ValueError('historical 7/8 result changed')
    record = next(r for r in refined['solutions'] if (r['L'], r['eta'], r['n_initial']) == (5.5, .3, 513))
    points = [v['geometry']['point'] for v in refined['physical'][0]['cases']]
    return record, refined['profiles']['5.5'], points, historical


def calibrate():
    with mp.workdps(80):
        point = [mp.mpf('.7'), mp.mpf('.91'), mp.mpf('.37')]
        rows = {}
        for name, sphere, tensor in [('sphere', True, False), ('cylinder', False, False), ('tensor', False, True)]:
            values = c.calibration_values(sphere, tensor)
            direct = c.coordinate_check(values, point, f=1, cosmological=0)
            expected = dict(R=mp.mpf(6 if sphere else 2),
                            trace_K=2*point[0] if tensor else mp.mpf(0),
                            K2=2*point[0]**2 if tensor else mp.mpf(0),
                            divergence=[mp.mpf(-2 if tensor else 0), mp.mpf(0), mp.mpf(0)])
            row = dict(reference=packed(direct), expected=packed(expected))
            if sphere:
                row['finite_differences'] = [packed(c.coordinate_check(values, point, mp.mpf(h), f=1, cosmological=0))
                                             for h in ('.008', '.004', '.002')]
            rows[name] = row
        return rows


def calibration_pass(rows):
    try:
        if set(rows) != {'sphere', 'cylinder', 'tensor'}:
            return False
        for row in rows.values():
            if not agreement({k: row['reference'][k] for k in row['expected']}, row['expected']):
                return False
        errors = [abs(mp.mpf(v['R'])-6) for v in rows['sphere']['finite_differences']]
        return len(errors) == 3 and all(8 <= a/b <= 24 for a, b in zip(errors, errors[1:]))
    except (KeyError, ValueError, TypeError, ZeroDivisionError):
        return False


def measure(record, profiles, point):
    refs = []
    for digits in (60, 80):
        with mp.workdps(digits):
            data = c.SavedData(record, profiles)
            x = [mp.mpf(v) for v in point]
            refs.append(packed(c.coordinate_check(data.values, x)))
    with mp.workdps(80):
        data = c.SavedData(record, profiles)
        x = [mp.mpf(v) for v in point]
        h = data.step(x)
        return dict(point=packed(x), h0=packed(h), reference60=refs[0], reference80=refs[1],
                    local=[dict(h=packed(h/2**j), values=packed(c.coordinate_check(data.values, x, h/2**j))) for j in range(3)],
                    historical_steps=[dict(h=step, values=packed(c.coordinate_check(data.values, x, mp.mpf(step))))
                                      for step in ('.008', '.004', '.002')])


def row_score(row):
    ref = row['reference80']
    hs, ms = mp.mpf(ref['Hscale']), mp.mpf(ref['Mscale'])
    if hs < 1 or ms < 1:
        raise ValueError('invalid normalization')
    H, M = mp.mpf(ref['H']), list(map(mp.mpf, ref['M']))
    errors = [[], []]
    absolute = [abs(H)/hs, c.norm(M)/ms]
    for item in row['local']:
        values = item['values']
        errors[0].append(abs(mp.mpf(values['H'])-H)/hs)
        errors[1].append(c.norm([mp.mpf(v)-r for v, r in zip(values['M'], M)])/ms)
    last = row['local'][-1]['values']
    absolute += [abs(mp.mpf(last['H']))/hs, c.norm(list(map(mp.mpf, last['M'])))/ms]
    ratios, active, converges = [], [], True
    for sequence in errors:
        rr, count = [], 0
        for a, b in zip(sequence, sequence[1:]):
            ratio = a/b if b else None
            rr.append(None if ratio is None else packed(ratio))
            if a > mp.mpf('1e-40'):
                count += 1
                converges &= ratio is not None and 8 <= ratio <= 24
        ratios.append(rr)
        active.append(count)
    checks = dict(precision=agreement(row['reference60'], ref),
                  absolute=all(mp.isfinite(v) and v < mp.mpf('1e-5') for v in absolute),
                  convergence=bool(converges),
                  differentiation_accuracy=all(v[-1] < mp.mpf('1e-8') for v in errors),
                  wrong_f=abs(mp.mpf(ref['wrong_f_H'])) > 10*max(absolute[:2]))
    historical = []
    for item in row['historical_steps']:
        v = item['values']
        historical.append(dict(h=item['h'], H=packed(mp.mpf(v['H'])/hs),
                               H_error=packed((mp.mpf(v['H'])-H)/hs),
                               M_error=packed(c.norm([mp.mpf(z)-r for z, r in zip(v['M'], M)])/ms)))
    return dict(checks=checks, absolute=packed(absolute), errors=packed(errors),
                ratios=ratios, active=active, historical_steps=historical)


def provenance():
    return dict(freeze=FREEZE, original_sha256=refinement.ORIGINAL_HASH,
                refinement_sha256=REFINEMENT_HASH,
                sources={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCES})


def run(inputs):
    calibration = calibrate()
    with mp.workdps(80):
        if not calibration_pass(calibration):
            raise ArithmeticError('coordinate calibration failed; data measurement not started')
    record, profiles, points, _ = inputs
    rows = []
    for i, point in enumerate(points):
        rows.append(measure(record, profiles, point))
        print(f'physical point {i+1}/{len(points)} measured', flush=True)
    return dict(**provenance(), dependencies={name: importlib.metadata.version(name) for name in ('mpmath', 'numpy', 'scipy', 'sympy')},
                calibration=calibration, rows=rows)


def score(data, inputs, replay=True):
    record, profiles, points, historical = inputs
    gates = dict(historical['gates'])
    gates['physical'] = False
    evidence = all(data.get(k) == v for k, v in provenance().items())
    summaries = []
    with mp.workdps(80):
        calibration = calibrate()
        evidence &= agreement(data['calibration'], calibration)
        evidence &= len(data['rows']) == len(points) == 20
        for row, point in zip(data['rows'], points):
            model = c.SavedData(record, profiles)
            x = [mp.mpf(v) for v in point]
            h = model.step(x)
            evidence &= row['point'] == packed(x) and row['h0'] == packed(h)
            evidence &= [v['h'] for v in row['local']] == [packed(h/2**j) for j in range(3)]
            evidence &= [v['h'] for v in row['historical_steps']] == ['.008', '.004', '.002']
            if replay:
                evidence &= agreement(row, measure(record, profiles, point))
            summaries.append(row_score(row))
        physical = (calibration_pass(calibration) and len(summaries) == 20
                    and all(all(row['checks'].values()) for row in summaries)
                    and all(sum(row['active'][k] for row in summaries) > 0 for k in (0, 1)))
    gates['evidence'] &= bool(evidence)
    gates['physical'] = bool(physical)
    passed = all(v for k, v in gates.items() if k != 'localization')
    return dict(freeze=FREEZE, gates=gates, passed=sum(gates.values()), total=8,
                verdicts=dict(FOUR_SCALAR_HANDLE_CONSTRAINT_DATA=passed,
                              LOCALIZED_BULK_MOUTH_INITIAL_DATA=passed and gates['localization']),
                historical=dict(original='6/8', stable_reconstruction='7/8', unchanged=True),
                calibration_pass=calibration_pass(calibration), points=summaries,
                unestablished=historical['unestablished'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--rescore', type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    verdict_path = args.output_dir/'error_budget_verdict.json'
    report_path = args.output_dir/'error_budget.md'
    verdict_path.write_text(json.dumps(dict(verdicts=FALSE, error='incomplete run'))+'\n')
    report_path.write_text('# Physical error budget incomplete\n\nNo affirmative verdict.\n')
    try:
        inputs = load_inputs(args.input_dir)
        data = json.loads(args.rescore.read_text()) if args.rescore else run(inputs)
        # Replay is mandatory for archived input. A fresh run already measured
        # every row from the hash-verified inputs in this process.
        result = score(data, inputs, replay=bool(args.rescore))
        (args.output_dir/'error_budget.json').write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
        verdict_path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
        report_path.write_text('# Independent physical error budget\n\n'
                              f"Extension: **{result['passed']}/8**. Original 6/8 and stable reconstruction 7/8 remain unchanged.\n\n"
                              f'Public freeze: `{FREEZE}`.\n\n'
                              '```json\n'+json.dumps(result['verdicts'], indent=2)+'\n```\n')
        print(json.dumps({k: result[k] for k in ('passed', 'gates', 'verdicts')}, indent=2))
        return 0 if all(result['gates'].values()) else 1
    except Exception as error:
        verdict_path.write_text(json.dumps(dict(verdicts=FALSE, error=str(error)))+'\n')
        report_path.write_text('# Physical error budget failed\n\nNo affirmative verdict.\n\n'+str(error)+'\n')
        raise


if __name__ == '__main__':
    raise SystemExit(main())
