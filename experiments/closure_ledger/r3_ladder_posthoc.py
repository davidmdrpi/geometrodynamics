"""Post-hoc diagnostics of the ladder (written after the registered labels were fixed; descriptive only).

sigma = -ln(Lambda)/Q_h is the per-harmonic decay exponent: for an analytic family, the order-Q resonant
coefficient behaves like exp(-Q sigma), sigma being the width of the analyticity strip in the angle. The
uniform-radius model of the specification is the special case sigma = ln(R/a). Also reported: the true
error of the double-precision machinery (|lambda_double - lambda_hp|), and agreement of the 2/5 rung with
the archived r3_breaking scan.
Output: runs/20261010_r3_ladder/posthoc.json.
"""
import json
import numpy as np
from experiments.closure_ledger import r3_ladder_probe as lp

OLD = lp.ROOT/'experiments/closure_ledger/runs/20261005_r3_breaking/scan_main.json'


def main():
    res = lp._read('result.json')['result']
    rows = []
    for p, q in lp.RUNGS:
        r = res['rungs'][lp.tag(p, q)]
        h = np.array(r['harmonics'])
        Q, a, L = r['Q_h'], r['a'], r['Lambda']
        sigma = -np.log(L)/Q
        rows.append(dict(rung=f'{p}/{q}', a=a, Q_h=Q, Lambda=L, sigma=float(sigma), R=float(a*np.exp(sigma)),
                         harmonic_Q=float(h[Q]), largest_other=float(np.delete(h[1:], Q-1).max()),
                         double_error=r['double_error']))
    a = np.array([x['a'] for x in rows])
    s = np.array([x['sigma'] for x in rows])
    order = np.argsort(a)
    fit = np.polyfit(a, s, 1)
    hp25 = [x['lam_float'] for x in lp._read('hp_2_5.json')['rows']]
    old = [x['lam'] for x in json.loads(OLD.read_text())['points']]
    out = dict(rungs=rows, sigma_monotone_decreasing_in_a=bool(np.all(np.diff(s[order]) < 0)),
               sigma_linear_fit=dict(slope=float(fit[0]), intercept=float(fit[1]),
                                     rms=float(np.sqrt(np.mean((np.polyval(fit, a)-s)**2)))),
               two_fifths_vs_archive=dict(max_diff=float(np.max(np.abs(np.array(hp25)-old))),
                                          correlation=float(np.corrcoef(hp25, old)[0, 1])))
    (lp.RUN_DIR/'posthoc.json').write_text(json.dumps(dict(sources=lp.sources(), **out), indent=1))
    for x in rows:
        print(f"{x['rung']:>5} a={x['a']:.4f} Q={x['Q_h']:2d} Lambda={x['Lambda']:.2e} sigma={x['sigma']:.3f} R={x['R']:.2f}"
              f" hQ/other={x['harmonic_Q']/x['largest_other']:.1f} double_err={x['double_error']:.1e}")
    print(json.dumps({k: v for k, v in out.items() if k != 'rungs'}, indent=1))


if __name__ == '__main__':
    main()
