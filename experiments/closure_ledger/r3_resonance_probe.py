"""R3 refocusing-resonance sign gate (docs/r3_refocusing_resonance_prereg.md, freeze 2e984ac,
correction 6b55c5b). Runs the registered ladder, scores gates N1-N5 and the decision gate."""
import argparse
import hashlib
import json
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from geometrodynamics.waves import esu_floquet as fl
from geometrodynamics.waves import r3_resonance as r3

FREEZE, CORRECTION = '2e984ac', '6b55c5b'
ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = ROOT/'experiments/closure_ledger/runs/20260929_r3_resonance'
SOURCES = ('geometrodynamics/waves/r3_resonance.py', 'experiments/closure_ledger/r3_resonance_probe.py',
           'geometrodynamics/waves/nonlinear_supported_tt.py')
LADDER = (.01, .02, .04, .08, .16)
POLS = (1, -1)
DECISION_EPS = (.01, .02, .04)


def sources():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def one(job):
    eps, pol, name = job
    try:
        out = r3.track(eps, pol, r3.TOLERANCES[name])
        out['error'] = None
    except Exception as err:   # recorded; the run fails its gates
        out = dict(error=f'{type(err).__name__}: {err}')
    return dict(eps=eps, pol=pol, integrator=name, **out)


def measure(workers=4):
    jobs = [(e, p, n) for n in r3.TOLERANCES for e in LADDER for p in POLS]
    with Pool(workers) as pool:
        runs = pool.map(one, jobs)
    return dict(runs=runs, trace_T2=float(np.trace(fl.monodromy('T', 2))))


def score(raw):
    runs = {(r['eps'], r['pol'], r['integrator']): r for r in raw['runs']}
    a = float(np.arccos(raw['trace_T2']/2)/(2*np.pi))
    ref = runs.get((.01, 1, 'primary'))
    failures, rows = [], []
    rho_ref = r3.birkhoff(ref['increments']) if ref and not ref['error'] else None
    rho0 = None if rho_ref is None else min((1+a, 2-a), key=lambda v: abs(v-rho_ref))
    if rho0 is None:
        failures.append('reference run failed')
    for eps in LADDER:
        for pol in POLS:
            p, s = runs[(eps, pol, 'primary')], runs[(eps, pol, 'secondary')]
            row = dict(eps=eps, pol=pol)
            if p['error'] or s['error']:
                failures.append(f'run failed eps={eps} pol={pol}: {p["error"] or s["error"]}')
                rows.append(row)
                continue
            rho, half, rho_sec = (r3.birkhoff(p['increments']), r3.birkhoff(p['increments'][:r3.K//2]),
                                  r3.birkhoff(s['increments']))
            err = max(abs(rho-half), abs(rho-rho_sec))
            row.update(rho=rho, rho_half=half, rho_secondary=rho_sec, err=err,
                       max_residual=max(max(p['residuals']), max(s['residuals'])),
                       max_A_dev=max(max(p['A_dev']), max(s['A_dev'])),
                       radius_ratio=(min(x[0] for x in p['radii']+s['radii'])/eps,
                                     max(x[1] for x in p['radii']+s['radii'])/eps),
                       widenings=max(w['widenings'] for w in p['windows']+s['windows']))
            if rho0 is not None:
                row['delta'] = rho-rho0
                row['c'] = (rho-rho0)/eps**2
                row['resolved'] = bool(abs(rho-rho0) >= 100*err)
            if row['max_residual'] > 1e-8:
                failures.append(f'N1 eps={eps} pol={pol}')
            if row['max_A_dev'] > .1 or not (.2 <= row['radius_ratio'][0] and row['radius_ratio'][1] <= 5):
                failures.append(f'N2 eps={eps} pol={pol}')
            if err > 1e-5:
                failures.append(f'N3 eps={eps} pol={pol} err={err:.2e}')
            rows.append(row)
    out = dict(rho0=rho0, frac_candidates=(1+a, 2-a), rows=rows, failures=failures)
    if rho0 is None:
        out['verdict'] = 'UNRESOLVED'
        return out
    get = {(r['eps'], r['pol']): r for r in rows}
    ref_row = get[(.01, 1)]
    if 'rho' in ref_row and abs(ref_row['rho']-rho0) > 1e-3:
        failures.append('N4 linear control')
    slopes = {}
    for pol in POLS:
        sl = []
        for e1, e2 in zip(DECISION_EPS, DECISION_EPS[1:]):
            r1, r2 = get[(e1, pol)], get[(e2, pol)]
            if r1.get('resolved') and r2.get('resolved'):
                sl.append(float(np.log(abs(r2['delta'])/abs(r1['delta']))/np.log(e2/e1)))
        slopes[pol] = sl
        if not sl or any(not 1.7 <= v <= 2.3 for v in sl):
            failures.append(f'N5 pol={pol} slopes={sl}')
    out['slopes'] = {str(k): v for k, v in slopes.items()}
    signs = {int(np.sign(get[(e, p)]['c'])) for e in DECISION_EPS for p in POLS if get[(e, p)].get('resolved')}
    out['c_signs'] = sorted(signs)
    crossings = [(r['eps'], r['pol']) for r in rows if 'rho' in r and (r['rho']-1.5)*(rho0-1.5) < 0]
    out['descriptive_crossings_of_3_2'] = crossings
    if failures or len(signs) != 1:
        out['verdict'] = 'UNRESOLVED'
        return out
    sign_c = signs.pop()
    if sign_c == np.sign(1.5-rho0):
        out['verdict'] = 'PASS'
        out['eps_star_estimate'] = float(np.sqrt((1.5-rho0)/ref_row['c']))
    else:
        out['verdict'] = 'FAIL'
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=RUN_DIR/'r3_resonance.json')
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    src = sources()
    args.output.write_text(json.dumps(dict(freeze=FREEZE, verdict='UNRESOLVED', error='incomplete'))+'\n')
    raw = measure(args.workers)
    if sources() != src:
        raise RuntimeError('sources changed during the run')
    result = score(raw)
    args.output.write_text(json.dumps(dict(freeze=FREEZE, correction=CORRECTION, sources=src,
                                           raw=raw, result=result), indent=1, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=1))
    for r in result['rows']:
        print({k: (round(v, 10) if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == '__main__':
    main()
