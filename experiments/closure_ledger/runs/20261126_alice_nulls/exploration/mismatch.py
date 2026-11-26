# IN-SAMPLE: actual minus perfect-selectivity doublet at a1 = (0.5, 2), exact nulls, eta = 0.5 and 0.9;
# and full-150 against doublet-only at the same points.
import importlib.util, json, sys, numpy as np
s = importlib.util.spec_from_file_location('d', '/tmp/claude-0/-home-user-geometrodynamics/1c9d2eb5-541f-57e7-9180-5d1a6778cfd2/scratchpad/doublet_dissect.py')
ns = {}; exec(open(s.origin).read().split("out = []")[0], ns)
V, doublet_only, variant, defect = ns['V'], ns['doublet_only'], ns['variant'], ns['defect']
out = []
spec = V.BASE.solve(0.5, 2, 2000); loop, sec, lam, Vm = spec
for k in (5, 8, 13, 20):
    g = V.Sweep(0.5, 2, 2000, 150, spec, 8.0, 'sin2'); tau = k * np.pi / g.mean
    g = V.Sweep(0.5, 2, 2000, 150, spec, tau, 'sin2'); d = doublet_only(g)
    Vk = Vm[:, [g.iA, g.iB]]; d._H = np.pi * 2 * (Vk.T @ (np.repeat(loop.mass, 2)[:, None] * V.BASE.rotate_pairs(Vk)))
    p = variant(d, perfect=True)
    row = dict(k=k, tau=float(tau))
    for eta in (0.5, 0.9):
        row[f"full_{eta}"] = defect(g, eta); row[f"doublet_{eta}"] = defect(d, eta); row[f"perfect_{eta}"] = defect(p, eta)
    out.append(row)
    print(k, 'tau %.3f' % tau, ' '.join('%s %.3e' % (kk, v) for kk, v in row.items() if kk not in ('k', 'tau')), flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
