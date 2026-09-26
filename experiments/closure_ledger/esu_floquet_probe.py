"""Frozen Floquet spectrum and antipodal refocusing probe (docs/esu_floquet_refocusing_prereg.md).

Raw measurements (maps, integrator differences, residual scans, control numbers)
are separated from every derived quantity (observables, gates, labels, fits,
verdicts). `score` rebuilds the complete derived structure from raw evidence and a
validated G1 record; `replay` compares that reconstruction, and optionally a fresh
measurement, against an archive. No stored flag or label is trusted.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from geometrodynamics.waves import esu_floquet as fl
from . import esu_floquet_symbolic as symbolic

FREEZE = '4e65c3e'
ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('geometrodynamics/waves/esu_floquet.py', 'experiments/closure_ledger/esu_floquet_probe.py',
           'docs/esu_floquet_refocusing_prereg.md')
RUN_DIR = ROOT/'experiments/closure_ledger/runs/20260926_esu_floquet'
G1_ARCHIVE = RUN_DIR/'g1_symbolic.json'
DEGREES = list(range(2, 81))
RK4_STEPS = (2**16, 2**17)
PRIOR_T2_HALF_TRACE = -0.0963065402      # #294, one pi/2 period at phase zero
VERDICT_KEYS = ('T_STABILITY', 'V_STABILITY', 'S_STABILITY', 'T_REFOCUSING', 'V_REFOCUSING',
                'S_REFOCUSING', 'T_WKB_PREDICTION', 'V_WKB_PREDICTION')


def digest(p):
    return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def rk4_all(sector, degrees, steps, t1=np.pi, **kw):
    """Vectorised fixed-step RK4 fundamental matrices for all degrees at once."""
    d = fl.DIM.get(sector, 2)
    N = len(degrees)
    ncol = np.repeat(np.asarray(degrees, dtype=float), d)
    Y = np.tile(np.eye(d), (1, N))       # columns: degree-major blocks of identity
    f = fl.RHS[sector]
    h = t1/steps
    t = 0.
    for _ in range(steps):
        k1 = f(t, Y, ncol, **kw); k2 = f(t+h/2, Y+h*k1/2, ncol, **kw)
        k3 = f(t+h/2, Y+h*k2/2, ncol, **kw); k4 = f(t+h, Y+h*k3, ncol, **kw)
        Y = Y+h*(k1+2*k2+2*k3+k4)/6
        t += h
    return [Y[:, i*d:(i+1)*d] for i in range(N)]


def residual_scan(sector, n, rng, samples=60):
    if sector == 'T':
        return 0.
    d = fl.DIM[sector]
    y0 = rng.uniform(-1, 1, d)
    f = (lambda t, y: fl.rhs_scalar(t, y, n)) if sector == 'S' else (lambda t, y: fl.rhs_vector(t, y, n))
    sol = solve_ivp(f, (0, np.pi), y0, method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
    check = fl.scalar_trace_residual if sector == 'S' else fl.vector_constraint_residual
    return float(max(check(t, sol.sol(t), n) for t in np.linspace(0, np.pi, samples)))


def fit_power(ns, deficits):
    ns, deficits = np.asarray(ns, float), np.asarray(deficits, float)
    keep = deficits > 0
    if keep.sum() < 5:
        return None
    x, y = np.log(ns[keep]), np.log(deficits[keep])
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    dof = max(len(x)-2, 1)
    s2 = float(np.sum((y-A@coef)**2)/dof)
    cov = s2*np.linalg.inv(A.T@A)
    return dict(p=float(-coef[0]), p_se=float(np.sqrt(cov[0, 0])), logC=float(coef[1]), points=int(keep.sum()))


def classify_refocusing(rows):
    if max(r['defect'] for r in rows) < 1e-8:
        return 'EXACT', None
    sel = [r for r in rows if 20 <= r['n'] <= 80]
    fit = fit_power([r['n'] for r in sel], [1-r['fidelity'] for r in sel])
    high = [r['fidelity'] for r in rows if 40 <= r['n'] <= 80]
    if fit and fit['p']-2*fit['p_se'] > 0 and min(high) >= .99:
        return 'ASYMPTOTIC', fit
    tailv = [r['fidelity'] for r in rows if 60 <= r['n'] <= 80]
    tail = float(np.mean(tailv))
    # Operationalisation of "converges to a limit", fixed before any result was viewed:
    # the F_n spread over n in [60,80] must be below .01.
    converged = max(tailv)-min(tailv) < .01
    if fit and abs(fit['p']) < 2*fit['p_se'] and tail < .999 and converged:
        return 'PLATEAU', fit
    return 'ABSENT', fit


def odd_subset(X, rows):
    return [r for r in rows if (r['n'] % 2 == 0 if X in 'TS' else r['n'] % 2 == 1)]


def stability_label(subset):
    h = [r['n'] for r in subset if r['stability'] == 'HYPERBOLIC']
    m = [r['n'] for r in subset if r['stability'] == 'MARGINAL']
    if h:
        return f'HYPERBOLIC_AT{h}'
    if m:
        return f'MARGINAL_AT{m}'
    return 'ELLIPTIC_2_TO_80'


# ---------------------------------------------------------------- raw evidence
def measure_controls():
    C1 = [fl.observables('V', n, fl.monodromy('C', n))['defect'] for n in DEGREES]
    C2 = [fl.observables('V', n, fl.monodromy('V', n, coupled=False))['defect'] for n in DEGREES]
    C4 = []
    for n in DEGREES:
        theta = fl.observables('T', n, fl.monodromy('B', n))['theta']
        a = abs(np.pi*(np.sqrt(n*(n+2))-(n+1))) % (2*np.pi)
        C4.append(abs(theta-min(a, 2*np.pi-a)))
    return dict(C1_defects=[float(v) for v in C1], C2_defects=[float(v) for v in C2],
                C3_T2_half_trace=float(np.trace(fl.monodromy('T', 2, t1=np.pi/2))),
                C4_theta_errors=[float(v) for v in C4])


def measure_maps(X, degrees):
    return ([fl.monodromy(X, n).tolist() for n in degrees],
            [fl.monodromy(X, n, t1=np.pi/2).tolist() for n in degrees])


def measure_sector(X, rng):
    primary, half = measure_maps(X, DEGREES)
    secondary = [np.trace(fl.monodromy(X, n, rtol=1e-10, atol=1e-12)) for n in DEGREES]
    rk = {s: rk4_all(X, DEGREES, s) for s in RK4_STEPS}
    rows = []
    for i, n in enumerate(DEGREES):
        tp = np.trace(primary[i])
        rows.append(dict(n=n, matrix=primary[i], half_period_matrix=half[i],
                         secondary_trace_diff=float(abs(secondary[i]-tp)),
                         rk4_16_diff=float(abs(np.trace(rk[RK4_STEPS[0]][i])-tp)),
                         rk4_17_diff=float(abs(np.trace(rk[RK4_STEPS[1]][i])-tp)),
                         constraint_residual=residual_scan(X, n, rng)))
    phase = [float(abs(np.trace(fl.monodromy(X, n, t0=.3, t1=.3+np.pi))-np.trace(primary[n-2])))
             for n in range(2, 11)]
    return dict(rows=rows, phase_independence=phase)


def measure():
    rng = np.random.default_rng(2026092611)
    raw = dict(wkb=fl.wkb_masses(), controls=measure_controls(), sectors={})
    for X in fl.SECTORS:
        print('sector', X, flush=True)
        raw['sectors'][X] = measure_sector(X, rng)
    return raw


# ------------------------------------------------------------- derived result
def score(raw, g1_record):
    """Rebuild every observable, gate, label and verdict from raw evidence."""
    c = raw['controls']
    controls = dict(C1_max_defect=max(c['C1_defects']), C2_max_defect=max(c['C2_defects']),
                    C3_T2_half_trace=c['C3_T2_half_trace'], C3_prior=PRIOR_T2_HALF_TRACE,
                    C3_error=abs(c['C3_T2_half_trace']-PRIOR_T2_HALF_TRACE),
                    C4_max_theta_error=max(c['C4_theta_errors']))
    controls['pass'] = bool(len(c['C1_defects']) == len(c['C2_defects']) == len(c['C4_theta_errors']) == len(DEGREES)
                            and controls['C1_max_defect'] < 1e-10 and controls['C2_max_defect'] < 1e-10
                            and controls['C3_error'] < 1e-8 and controls['C4_max_theta_error'] < 1e-8)
    g1 = bool(symbolic.g1_valid(g1_record))
    wkb = raw['wkb']
    sectors = {}
    for X in fl.SECTORS:
        rs = raw['sectors'][X]
        if [r['n'] for r in rs['rows']] != DEGREES:
            raise ValueError('changed degree schedule')
        rows = []
        for r in rs['rows']:
            M, H = np.asarray(r['matrix'], float), np.asarray(r['half_period_matrix'], float)
            if not (np.isfinite(M).all() and np.isfinite(H).all()):
                raise ValueError('nonfinite map')
            o = fl.observables(X, r['n'], M)
            o['half_max_abs_multiplier'] = float(np.max(abs(np.linalg.eigvals(H))))
            o['stability'] = fl.stability(o['max_abs_multiplier'])
            e16, e17 = r['rk4_16_diff'], r['rk4_17_diff']
            o['integrators'] = dict(secondary_trace_diff=r['secondary_trace_diff'], rk4_16_diff=e16, rk4_17_diff=e17,
                                    rk4_ratio=float(e16/e17) if e16 > 1e-11 and e17 > 0 else None)
            o['constraint_residual'] = r['constraint_residual']
            rows.append(o)
        g2 = max(r['constraint_residual'] for r in rows)
        g3 = all(r['integrators']['rk4_17_diff'] < 1e-7 and r['integrators']['secondary_trace_diff'] < 1e-7
                 and (r['integrators']['rk4_ratio'] is None or 8 <= r['integrators']['rk4_ratio'] <= 32)
                 for r in rows)
        g4 = max(abs(r['det']-1) for r in rows)
        gates = dict(G1=g1, G2_max_residual=g2, G2=bool(g2 < 1e-8), G3=bool(g3), G4_max_det_error=g4,
                     G4=bool(g4 < 1e-8), phase_independence_max=float(max(rs['phase_independence'])),
                     controls=controls['pass'])
        gates['non_G3_pass'] = bool(g1 and gates['G2'] and gates['G4'] and controls['pass'])
        verdict_ok = gates['non_G3_pass'] and gates['G3']
        odd = odd_subset(X, rows)
        refoc, fit = classify_refocusing(rows)
        refoc_odd, fit_odd = classify_refocusing(odd)
        sec = dict(rows=rows, gates=gates, stability=stability_label(rows) if verdict_ok else 'UNRESOLVED',
                   stability_odd=stability_label(odd) if verdict_ok else 'UNRESOLVED',
                   refocusing=refoc if verdict_ok else 'UNRESOLVED', fit=fit,
                   refocusing_odd=refoc_odd if verdict_ok else 'UNRESOLVED', fit_odd=fit_odd)
        if X in 'TV':
            m2 = wkb['m_T2'] if X == 'T' else wkb['m_V2']
            pred = float(np.pi*abs(m2)/2)
            meas = float(np.mean([(r['n']+1)*r['theta'] for r in rows if 40 <= r['n'] <= 80 and r['theta'] is not None]))
            sec['wkb'] = dict(predicted=pred, measured=meas, relative_error=abs(meas-pred)/pred,
                              verdict=('PASS' if abs(meas-pred) <= .1*pred else 'FAIL') if verdict_ok else 'UNRESOLVED')
        sectors[X] = sec
    verdicts = dict(
        T_STABILITY=sectors['T']['stability'], V_STABILITY=sectors['V']['stability'],
        S_STABILITY=sectors['S']['stability'], T_REFOCUSING=sectors['T']['refocusing'],
        V_REFOCUSING=sectors['V']['refocusing'], S_REFOCUSING=sectors['S']['refocusing'],
        T_WKB_PREDICTION=sectors['T']['wkb']['verdict'], V_WKB_PREDICTION=sectors['V']['wkb']['verdict'])
    odd_verdicts = {f'{X}_{k}': sectors[X][k+'_odd'] for X in fl.SECTORS for k in ('stability', 'refocusing')}
    return dict(controls=controls, sectors=sectors, verdicts=verdicts, verdicts_odd_sector=odd_verdicts)


def provenance(g1_record):
    return dict(freeze=FREEZE, sources={p: digest(p) for p in SOURCES}, degrees=DEGREES,
                g1_sha256=hashlib.sha256(json.dumps(g1_record, sort_keys=True).encode()).hexdigest())


def load_g1(path=G1_ARCHIVE):
    return json.loads(Path(path).read_text())


def run(g1_record=None):
    g1_record = load_g1() if g1_record is None else g1_record
    raw = measure()
    return dict(**provenance(g1_record), raw=raw, result=score(raw, g1_record))


# ---------------------------------------------------------------------- replay
def close(a, b, tol):
    """Structural equality: labels/booleans/None exactly, numbers within scaled tol."""
    if isinstance(b, dict):
        return isinstance(a, dict) and a.keys() == b.keys() and all(close(a[k], v, tol) for k, v in b.items())
    if isinstance(b, (list, tuple)):
        return isinstance(a, (list, tuple)) and len(a) == len(b) and all(close(x, y, tol) for x, y in zip(a, b))
    if isinstance(b, (bool, str)) or b is None:
        return type(a) is type(b) and a == b
    if isinstance(a, bool) or not isinstance(a, (int, float)):
        return False
    return bool(np.isfinite(a) and abs(a-b) <= tol*max(1., abs(b)))


def replay(data, g1_record=None, degrees=None, full=False, tol=1e-9):
    """Evidence check.

    Always: provenance, freeze and source hashes; a validated G1 record whose hash
    matches the archive; and complete reconstruction of the derived structure (all
    rows, gates, labels, fits, WKB and verdicts, including odd-sector results) from
    raw evidence, compared exactly for labels and within tol for numbers.
    Maps: primary and half-period maps are recomputed for `degrees` (all by default).
    full=True: remeasure all raw evidence and compare it as well.
    """
    try:
        g1_record = load_g1() if g1_record is None else g1_record
        ok = close(data.get('freeze'), FREEZE, 0) and data.get('degrees') == DEGREES
        ok &= all(data['sources'].get(p) == digest(p) for p in SOURCES)
        ok &= data.get('g1_sha256') == provenance(g1_record)['g1_sha256']
        ok &= symbolic.g1_valid(g1_record)
        ok &= close(data['result'], score(data['raw'], g1_record), tol)
        if full:
            ok &= close(data['raw'], measure(), tol)
        else:
            for X in fl.SECTORS:
                rows = data['raw']['sectors'][X]['rows']
                for n in (DEGREES if degrees is None else degrees):
                    M, H = measure_maps(X, [n])
                    ok &= close(rows[n-2]['matrix'], M[0], tol) and close(rows[n-2]['half_period_matrix'], H[0], tol)
            if degrees is None:
                ok &= close(data['raw']['controls'], measure_controls(), tol)
        return bool(ok)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError):
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--replay', type=Path)
    ap.add_argument('--full', action='store_true', help='with --replay: remeasure all raw evidence')
    args = ap.parse_args()
    if args.replay:
        ok = replay(json.loads(args.replay.read_text()), full=args.full)
        print('replay evidence:', ok)
        raise SystemExit(0 if ok else 1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(result=dict(verdicts=dict.fromkeys(VERDICT_KEYS, 'UNRESOLVED')), error='incomplete'))+'\n')
    try:
        result = run()
    except Exception as err:
        args.output.write_text(json.dumps(dict(result=dict(verdicts=dict.fromkeys(VERDICT_KEYS, 'UNRESOLVED')), error=str(err)))+'\n')
        raise
    args.output.write_text(json.dumps(result, indent=1, allow_nan=False)+'\n')
    r = result['result']
    print(json.dumps(dict(verdicts=r['verdicts'], odd=r['verdicts_odd_sector'], controls=r['controls']), indent=1))


if __name__ == '__main__':
    main()
