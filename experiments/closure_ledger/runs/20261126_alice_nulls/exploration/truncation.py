# IN-SAMPLE truncation check at a1 = (0.5, 2), the null region of the V-curve run.
import copy, importlib.util, json, sys, numpy as np
sys.path.insert(0, '.')
s = importlib.util.spec_from_file_location('v', 'experiments/closure_ledger/alice_vcurve_probe.py'); V = importlib.util.module_from_spec(s); s.loader.exec_module(V)

def restrict(geo, idx):
    """The same geometry with the modal basis cut to idx (doublet first)."""
    idx = np.asarray(idx)
    assert idx[0] == geo.iA or idx[1] == geo.iA
    g = copy.copy(geo)
    g.keep = len(idx)
    g.omega_all = geo.omega_all[idx]
    g.W1 = geo.W1[np.ix_(idx, idx)]
    g.W0 = geo.W0[np.ix_(idx, idx)]
    g.iA, g.iB = int(np.flatnonzero(idx == geo.iA)[0]), int(np.flatnonzero(idx == geo.iB)[0])
    return g

def defect(g, eta=0.5):
    g0 = g.calibrate(eta, g.mean)
    return g.run(g0, g.mean)["max_abs_defect"]

spec = V.BASE.solve(0.5, 2, 2000)
even = [4, 5, 8, 9, 12, 13, 16, 17, 20, 21]
subsets = {
    "doublet": [0, 1],
    "+l4": [0, 1, 4, 5],
    "+l4,l6": [0, 1, 4, 5, 8, 9],
    "+l4..l8": [0, 1, 4, 5, 8, 9, 12, 13],
    "+even<=l12": [0, 1] + even,
    "full-150": list(range(150)),
    "full minus l4": [k for k in range(150) if k not in (4, 5)],
    "full minus l6": [k for k in range(150) if k not in (8, 9)],
}
out = {}
for tau in (6.4, 7.0, 7.7, 8.3, 9.0, 9.6, 10.3):
    base = V.Sweep(0.5, 2, 2000, 150, spec, tau, 'sin2')
    row = {name: defect(restrict(base, idx)) for name, idx in subsets.items()}
    out[str(tau)] = row
    print(tau, ' '.join('%s=%.2e' % (k, v) for k, v in row.items()), flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
