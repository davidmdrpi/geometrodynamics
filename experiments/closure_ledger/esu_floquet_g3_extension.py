"""Prospective G3 extension (docs/esu_floquet_g3_extension_prereg.md, addendum 9a8ef99)."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves import esu_floquet as fl
from . import esu_floquet_probe as probe

ADDENDUM = '9a8ef99'
ARCHIVE = probe.ROOT/'experiments/closure_ledger/runs/20260926_esu_floquet/esu_floquet.json'
STEPS = (2**10, 2**11, 2**12)
FLOOR = 1e-9


def run(data):
    if not probe.replay(data, degrees=[]):
        raise ValueError('archived maps do not replay')
    out = dict(addendum=ADDENDUM, archive_sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
               steps=list(STEPS), floor=FLOOR, sectors={}, verdicts_ext={})
    for X in fl.SECTORS:
        rows = data['sectors'][X]['rows']
        mats = {s: probe.rk4_all(X, probe.DEGREES, s) for s in STEPS}
        per, ok, nonvacuous = [], True, False
        for i, r in enumerate(rows):
            e = [abs(float(np.trace(mats[s][i]))-r['trace']) for s in STEPS]
            ratios = []
            for a, b in zip(e, e[1:]):
                if a > FLOOR:
                    ratios.append(a/b if b > 0 else float('inf'))
            passed = all(8 <= q <= 32 for q in ratios)
            agree = r['integrators']['rk4_17_diff'] < 1e-7 and r['integrators']['secondary_trace_diff'] < 1e-7
            ok &= passed and agree
            nonvacuous |= bool(ratios) and r['n'] >= 20
            per.append(dict(n=r['n'], errors=e, ratios=ratios, passed=bool(passed and agree)))
        sec = dict(rows=per, passed=bool(ok and nonvacuous), nonvacuous=bool(nonvacuous))
        if sec['passed']:
            orig = data['sectors'][X]
            def stab(subset):
                h = [q['n'] for q in subset if q['stability'] == 'HYPERBOLIC']
                m = [q['n'] for q in subset if q['stability'] == 'MARGINAL']
                return f'HYPERBOLIC_AT{h}' if h else (f'MARGINAL_AT{m}' if m else 'ELLIPTIC_2_TO_80')
            odd = [q for q in rows if (q['n'] % 2 == 0 if X in 'TS' else q['n'] % 2 == 1)]
            out['verdicts_ext'][f'{X}_STABILITY_EXT'] = stab(rows)
            out['verdicts_ext'][f'{X}_REFOCUSING_EXT'] = probe.classify_refocusing(rows)[0]
            out['verdicts_ext'][f'{X}_STABILITY_ODD_EXT'] = stab(odd)
            out['verdicts_ext'][f'{X}_REFOCUSING_ODD_EXT'] = probe.classify_refocusing(odd)[0]
            if X == 'V':
                w = orig['wkb']
                out['verdicts_ext']['V_WKB_PREDICTION_EXT'] = 'PASS' if w['relative_error'] <= .1 else 'FAIL'
        else:
            out['verdicts_ext'][f'{X}_EXT'] = 'UNRESOLVED'
        out['sectors'][X] = sec
        print(X, 'passed' if sec['passed'] else 'FAILED', flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    args.output.write_text(json.dumps(dict(verdicts_ext={}, error='incomplete'))+'\n')
    try:
        result = run(json.loads(ARCHIVE.read_text()))
    except Exception as err:
        args.output.write_text(json.dumps(dict(verdicts_ext={}, error=str(err)))+'\n')
        raise
    args.output.write_text(json.dumps(result, indent=1, allow_nan=False)+'\n')
    print(json.dumps(result['verdicts_ext'], indent=1))


if __name__ == '__main__':
    main()
