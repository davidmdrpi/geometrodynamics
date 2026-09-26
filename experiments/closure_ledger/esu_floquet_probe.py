"""Frozen Floquet spectrum and antipodal refocusing probe (docs/esu_floquet_refocusing_prereg.md)."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from geometrodynamics.waves import esu_floquet as fl

FREEZE = '4e65c3e'
ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('geometrodynamics/waves/esu_floquet.py', 'experiments/closure_ledger/esu_floquet_probe.py',
           'docs/esu_floquet_refocusing_prereg.md')
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
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    dof = max(len(x)-2, 1)
    s2 = float(np.sum((y-A@coef)**2)/dof)
    cov = s2*np.linalg.inv(A.T@A)
    return dict(p=float(-coef[0]), p_se=float(np.sqrt(cov[0, 0])), logC=float(coef[1]), points=int(keep.sum()))


def classify_refocusing(rows):
    ns = [r['n'] for r in rows]
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


def run():
    rng = np.random.default_rng(2026092611)
    wkb = fl.wkb_masses()
    out = dict(freeze=FREEZE, sources={p: digest(p) for p in SOURCES}, degrees=DEGREES, wkb=wkb)
    # Controls
    C1 = [fl.observables('V', n, fl.monodromy('C', n)) for n in DEGREES]
    C2 = [fl.observables('V', n, fl.monodromy('V', n, coupled=False)) for n in DEGREES]
    t2half = float(np.trace(fl.monodromy('T', 2, t1=np.pi/2)))
    C4 = []
    for n in DEGREES:
        o = fl.observables('T', n, fl.monodromy('B', n))
        C4.append(dict(n=n, theta=o['theta'], analytic=abs(np.pi*(np.sqrt(n*(n+2))-(n+1))) % (2*np.pi)))
    c4_err = max(abs(r['theta']-min(r['analytic'], 2*np.pi-r['analytic'])) for r in C4)
    controls = dict(C1_max_defect=max(r['defect'] for r in C1), C2_max_defect=max(r['defect'] for r in C2),
                    C3_T2_half_trace=t2half, C3_prior=PRIOR_T2_HALF_TRACE,
                    C3_error=abs(t2half-PRIOR_T2_HALF_TRACE), C4_max_theta_error=float(c4_err))
    controls['pass'] = bool(controls['C1_max_defect'] < 1e-10 and controls['C2_max_defect'] < 1e-10
                            and controls['C3_error'] < 1e-8 and c4_err < 1e-8)
    out['controls'] = controls
    sectors = {}
    for X in fl.SECTORS:
        print('sector', X, flush=True)
        primary = [fl.monodromy(X, n) for n in DEGREES]
        secondary = [fl.monodromy(X, n, rtol=1e-10, atol=1e-12) for n in DEGREES]
        half = [fl.monodromy(X, n, t1=np.pi/2) for n in DEGREES]
        rk = {s: rk4_all(X, DEGREES, s) for s in RK4_STEPS}
        rows = []
        for i, n in enumerate(DEGREES):
            o = fl.observables(X, n, primary[i])
            o['matrix'] = primary[i].tolist()
            o['half_period_matrix'] = half[i].tolist()
            o['half_max_abs_multiplier'] = float(np.max(abs(np.linalg.eigvals(half[i]))))
            o['stability'] = fl.stability(o['max_abs_multiplier'])
            tp = np.trace(primary[i])
            e16, e17 = abs(np.trace(rk[RK4_STEPS[0]][i])-tp), abs(np.trace(rk[RK4_STEPS[1]][i])-tp)
            o['integrators'] = dict(secondary_trace_diff=float(abs(np.trace(secondary[i])-tp)),
                                    rk4_16_diff=float(e16), rk4_17_diff=float(e17),
                                    rk4_ratio=float(e16/e17) if e16 > 1e-11 and e17 > 0 else None)
            o['constraint_residual'] = residual_scan(X, n, rng)
            rows.append(o)
        phase = []
        for n in range(2, 11):
            Mp = fl.monodromy(X, n, t0=.3, t1=.3+np.pi)
            phase.append(abs(float(np.trace(Mp))-rows[n-2]['trace']))
        g2 = max(r['constraint_residual'] for r in rows)
        g3 = all(r['integrators']['rk4_17_diff'] < 1e-7 and r['integrators']['secondary_trace_diff'] < 1e-7
                 and (r['integrators']['rk4_ratio'] is None or 8 <= r['integrators']['rk4_ratio'] <= 32)
                 for r in rows)
        g4 = max(abs(r['det']-1) for r in rows)
        gates = dict(G2_max_residual=g2, G2=bool(g2 < 1e-8), G3=bool(g3), G4_max_det_error=g4,
                     G4=bool(g4 < 1e-8), phase_independence_max=float(max(phase)))
        verdict_ok = gates['G2'] and gates['G3'] and gates['G4'] and controls['pass']
        hyper = [r['n'] for r in rows if r['stability'] == 'HYPERBOLIC']
        marg = [r['n'] for r in rows if r['stability'] == 'MARGINAL']
        def stab(subset):
            h = [r['n'] for r in subset if r['stability'] == 'HYPERBOLIC']
            m = [r['n'] for r in subset if r['stability'] == 'MARGINAL']
            if h:
                return f'HYPERBOLIC_AT{h}'
            if m:
                return f'MARGINAL_AT{m}'
            return 'ELLIPTIC_2_TO_80'
        odd = [r for r in rows if (r['n'] % 2 == 0 if X in 'TS' else r['n'] % 2 == 1)]
        refoc, fit = classify_refocusing(rows)
        refoc_odd, fit_odd = classify_refocusing(odd)
        sec = dict(rows=rows, gates=gates, stability=stab(rows) if verdict_ok else 'UNRESOLVED',
                   stability_odd=stab(odd) if verdict_ok else 'UNRESOLVED',
                   refocusing=refoc if verdict_ok else 'UNRESOLVED', fit=fit,
                   refocusing_odd=refoc_odd if verdict_ok else 'UNRESOLVED', fit_odd=fit_odd)
        if X in 'TV':
            m2 = wkb['m_T2'] if X == 'T' else wkb['m_V2']
            pred = np.pi*abs(m2)/2
            meas = float(np.mean([(r['n']+1)*r['theta'] for r in rows if 40 <= r['n'] <= 80 and r['theta'] is not None]))
            sec['wkb'] = dict(predicted=pred, measured=meas, relative_error=abs(meas-pred)/pred,
                              verdict=('PASS' if abs(meas-pred) <= .1*pred else 'FAIL') if verdict_ok else 'UNRESOLVED')
        sectors[X] = sec
        print(X, sec['stability'], sec['refocusing'], gates, flush=True)
    out['sectors'] = sectors
    out['verdicts'] = dict(
        T_STABILITY=sectors['T']['stability'], V_STABILITY=sectors['V']['stability'],
        S_STABILITY=sectors['S']['stability'], T_REFOCUSING=sectors['T']['refocusing'],
        V_REFOCUSING=sectors['V']['refocusing'], S_REFOCUSING=sectors['S']['refocusing'],
        T_WKB_PREDICTION=sectors['T']['wkb']['verdict'], V_WKB_PREDICTION=sectors['V']['wkb']['verdict'])
    out['verdicts_odd_sector'] = {f'{X}_{k}': sectors[X][k+'_odd'] for X in fl.SECTORS for k in ('stability', 'refocusing')}
    return out


def replay(data, degrees=None, tol=1e-9):
    """Recompute stored maps and derived observables; any mismatch clears evidence."""
    ok = data.get('freeze') == FREEZE and data.get('degrees') == DEGREES
    ok &= all(data['sources'].get(p) == digest(p) for p in SOURCES)
    for X in fl.SECTORS:
        rows = data['sectors'][X]['rows']
        if [r['n'] for r in rows] != DEGREES:
            return False
        for r in rows:
            M = np.asarray(r['matrix'], float)
            if not np.isfinite(M).all():
                return False
            o = fl.observables(X, r['n'], M)
            for key in ('trace', 'det', 'fidelity', 'defect', 'max_abs_multiplier'):
                ok &= abs(o[key]-r[key]) <= tol*max(1, abs(r[key]))
            ok &= fl.stability(o['max_abs_multiplier']) == r['stability']
        for n in (degrees if degrees is not None else DEGREES):
            M = fl.monodromy(X, n)
            ok &= bool(np.max(abs(M-np.asarray(rows[n-2]['matrix']))) <= tol*max(1., np.max(abs(M))))
        refoc, _ = classify_refocusing(rows)
        if data['sectors'][X]['refocusing'] != 'UNRESOLVED':
            ok &= refoc == data['sectors'][X]['refocusing']
    return bool(ok)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--replay', type=Path)
    args = ap.parse_args()
    if args.replay:
        ok = replay(json.loads(args.replay.read_text()))
        print('replay evidence:', ok)
        raise SystemExit(0 if ok else 1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(verdicts=dict.fromkeys(VERDICT_KEYS, 'UNRESOLVED'), error='incomplete'))+'\n')
    try:
        result = run()
    except Exception as err:
        args.output.write_text(json.dumps(dict(verdicts=dict.fromkeys(VERDICT_KEYS, 'UNRESOLVED'), error=str(err)))+'\n')
        raise
    args.output.write_text(json.dumps(result, indent=1, allow_nan=False)+'\n')
    print(json.dumps(dict(verdicts=result['verdicts'], odd=result['verdicts_odd_sector'], controls=result['controls']), indent=1))


if __name__ == '__main__':
    main()
