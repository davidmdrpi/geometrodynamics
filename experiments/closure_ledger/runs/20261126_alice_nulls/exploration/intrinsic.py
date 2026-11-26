# IN-SAMPLE: the intrinsic null residual exactly at tau_k, doublet only, perfect selectivity, vs eta and tau.
import copy, importlib.util, json, sys, numpy as np
s = importlib.util.spec_from_file_location('d', __import__('pathlib').Path(__file__).with_name('doublet_dissect.py'))
src = open(s.origin).read().split("out = []")[0]
ns = {}; exec(src, ns)
V, doublet_only, variant, defect = ns['V'], ns['doublet_only'], ns['variant'], ns['defect']
from geometrodynamics.bulk import alice_handle_ports as ah
out = []
for a in (0.3, 0.5):
    spec = V.BASE.solve(a, 2, 2000)
    loop, sec, lam, Vm = spec
    wbar = None
    for k in (5, 6, 7, 8, 10, 13, 16, 20):
        g = V.Sweep(a, 2, 2000, 150, spec, 8.0, 'sin2')
        wbar = g.mean
        tau = k * np.pi / wbar
        g = V.Sweep(a, 2, 2000, 150, spec, tau, 'sin2')
        d = doublet_only(g)
        Vk = Vm[:, [g.iA, g.iB]]
        d._H = np.pi * 2 * (Vk.T @ (np.repeat(loop.mass, 2)[:, None] * V.BASE.rotate_pairs(Vk)))
        p = variant(d, perfect=True)
        row = dict(a=a, k=k, tau=float(tau), actual_doublet_eta05=defect(d, 0.5))
        for eta in (0.1, 0.25, 0.5, 0.9):
            row[f"perfect_eta{eta}"] = defect(p, eta)
        out.append(row)
        print(a, k, 'tau %.3f' % tau, 'actual %.2e' % row['actual_doublet_eta05'],
              ' '.join('eta%.2f %.2e' % (e, row[f'perfect_eta{e}']) for e in (0.1, 0.25, 0.5, 0.9)), flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
