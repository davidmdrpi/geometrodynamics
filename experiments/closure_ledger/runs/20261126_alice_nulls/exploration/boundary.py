# IN-SAMPLE: the reachable-depletion boundary at a1 = (0.5, 2) and a2 = (0.35, 2).
import importlib.util, json, sys, numpy as np
s = importlib.util.spec_from_file_location('v', 'experiments/closure_ledger/alice_vcurve_probe.py'); V = importlib.util.module_from_spec(s); s.loader.exec_module(V)
out = []
for a in (0.5, 0.35):
    spec = V.BASE.solve(a, 2, 2000)
    for kind in ('sin2', 'sin4'):
        for tau in (2.0, 2.5, 3.21, 4.0, 5.0, 6.0, 7.05, 9.0):
            g = V.Sweep(a, 2, 2000, 150, spec, tau, kind)
            strengths = 0.05 * 8.0 / tau * 1.25 ** np.arange(50)
            eta = np.array([g.absorbed_plus(x, g.mean) for x in strengths])
            i = int(np.argmax(eta))
            row = dict(a=a, profile=kind, tau=tau, eta_max=float(eta[i]), peak_rate_at_max=float(strengths[i]),
                       peak_rate_over_wbar=float(strengths[i] / g.mean), wbar_tau=float(g.mean * tau))
            # the peak rate needed for eta = 0.5 and 0.9 on the weak branch, where reachable
            for target in (0.5, 0.9):
                below = np.flatnonzero(eta[:i + 1] >= target)
                row[f"rate_over_wbar_for_{target}"] = None if not len(below) else float(strengths[below[0]] / g.mean)
            out.append(row)
            print(a, kind, tau, 'eta_max %.3f' % row['eta_max'], 'gpeak/wbar at max %.3f' % row['peak_rate_over_wbar'],
                  'needed for .5: %s, .9: %s' % (row['rate_over_wbar_for_0.5'], row['rate_over_wbar_for_0.9']), flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
