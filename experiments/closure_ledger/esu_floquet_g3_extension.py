"""Prospective G3 extension (docs/esu_floquet_g3_extension_prereg.md, addendum 9a8ef99).

The extension replaces only the G3 criterion. Every other gate (G1, G2, G4,
controls) is recomputed from the archive's raw evidence and a validated G1
record, and must pass. The recorded archive hash is that of the exact bytes
validated.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves import esu_floquet as fl
from . import esu_floquet_probe as probe

ADDENDUM = '9a8ef99'
ARCHIVE = probe.RUN_DIR/'esu_floquet.json'
STEPS = (2**10, 2**11, 2**12)
FLOOR = 1e-9


def measure(data):
    """RK4 trace errors against the archived primary maps (raw extension evidence)."""
    errors = {}
    for X in fl.SECTORS:
        rows = data['raw']['sectors'][X]['rows']
        mats = {s: probe.rk4_all(X, probe.DEGREES, s) for s in STEPS}
        errors[X] = [[abs(float(np.trace(mats[s][i]))-float(np.trace(np.asarray(r['matrix'])))) for s in STEPS]
                     for i, r in enumerate(rows)]
    return errors


def score(data, g1_record, errors):
    """Derive extension verdicts; affirmative only if all non-G3 gates pass."""
    base = probe.score(data['raw'], g1_record)
    out = dict(sectors={}, verdicts_ext={})
    for X in fl.SECTORS:
        sec_base = base['sectors'][X]
        rows = sec_base['rows']
        per, ok, nonvacuous = [], True, False
        for r, e in zip(rows, errors[X]):
            ratios = [a/b if b > 0 else None for a, b in zip(e, e[1:]) if a > FLOOR]
            passed = all(q is not None and 8 <= q <= 32 for q in ratios)
            agree = r['integrators']['rk4_17_diff'] < 1e-7 and r['integrators']['secondary_trace_diff'] < 1e-7
            ok &= passed and agree
            nonvacuous |= bool(ratios) and r['n'] >= 20
            per.append(dict(n=r['n'], errors=e, ratios=ratios, passed=bool(passed and agree)))
        g3_ext = bool(ok and nonvacuous)
        other = sec_base['gates']['non_G3_pass']
        out['sectors'][X] = dict(rows=per, g3_ext=g3_ext, nonvacuous=bool(nonvacuous), non_G3_gates=bool(other),
                                 passed=bool(g3_ext and other))
        if g3_ext and other:
            odd = probe.odd_subset(X, rows)
            out['verdicts_ext'][f'{X}_STABILITY_EXT'] = probe.stability_label(rows)
            out['verdicts_ext'][f'{X}_REFOCUSING_EXT'] = probe.classify_refocusing(rows)[0]
            out['verdicts_ext'][f'{X}_STABILITY_ODD_EXT'] = probe.stability_label(odd)
            out['verdicts_ext'][f'{X}_REFOCUSING_ODD_EXT'] = probe.classify_refocusing(odd)[0]
            if X == 'V':
                w = sec_base['wkb']
                out['verdicts_ext']['V_WKB_PREDICTION_EXT'] = 'PASS' if w['relative_error'] <= .1 else 'FAIL'
        else:
            out['verdicts_ext'][f'{X}_EXT'] = 'UNRESOLVED'
    return out


def run(archive_bytes, g1_record=None):
    g1_record = probe.load_g1() if g1_record is None else g1_record
    data = json.loads(archive_bytes)
    if not probe.replay(data, g1_record, degrees=[]):
        raise ValueError('archive evidence does not replay')
    errors = measure(data)
    return dict(addendum=ADDENDUM, archive_sha256=hashlib.sha256(archive_bytes).hexdigest(),
                steps=list(STEPS), floor=FLOOR, errors=errors, **score(data, g1_record, errors))


def replay(ext, archive_bytes, g1_record=None, full=False, tol=1e-9):
    """The extension record must bind these archive bytes and follow from them."""
    try:
        g1_record = probe.load_g1() if g1_record is None else g1_record
        data = json.loads(archive_bytes)
        ok = ext.get('addendum') == ADDENDUM and ext.get('archive_sha256') == hashlib.sha256(archive_bytes).hexdigest()
        ok &= ext.get('steps') == list(STEPS) and ext.get('floor') == FLOOR
        ok &= probe.replay(data, g1_record, degrees=[])
        derived = score(data, g1_record, ext['errors'])
        ok &= probe.close({k: ext[k] for k in ('sectors', 'verdicts_ext')}, derived, tol)
        if full:
            ok &= probe.close(ext['errors'], measure(data), tol)
        return bool(ok)
    except (KeyError, TypeError, ValueError, AttributeError):
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--archive', type=Path, default=ARCHIVE)
    ap.add_argument('--replay', type=Path)
    args = ap.parse_args()
    archive_bytes = args.archive.read_bytes()
    if args.replay:
        ok = replay(json.loads(args.replay.read_text()), archive_bytes, full=True)
        print('extension replay evidence:', ok)
        raise SystemExit(0 if ok else 1)
    args.output.write_text(json.dumps(dict(verdicts_ext={}, error='incomplete'))+'\n')
    try:
        result = run(archive_bytes)
    except Exception as err:
        args.output.write_text(json.dumps(dict(verdicts_ext={}, error=str(err)))+'\n')
        raise
    args.output.write_text(json.dumps(result, indent=1, allow_nan=False)+'\n')
    print(json.dumps(result['verdicts_ext'], indent=1))


if __name__ == '__main__':
    main()
