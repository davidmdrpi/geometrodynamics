# IN-SAMPLE: dissect the null residual inside the doublet-only model at a1 = (0.5, 2) and (0.3, 2).
import copy, importlib.util, json, sys, numpy as np
s = importlib.util.spec_from_file_location('v', 'experiments/closure_ledger/alice_vcurve_probe.py'); V = importlib.util.module_from_spec(s); s.loader.exec_module(V)

def doublet_only(geo):
    g = copy.copy(geo); idx = np.array([geo.iA, geo.iB])
    g.keep = 2; g.omega_all = geo.omega_all[idx]
    g.W1 = geo.W1[np.ix_(idx, idx)].copy(); g.W0 = geo.W0[np.ix_(idx, idx)].copy(); g.iA, g.iB = 0, 1
    return g

def variant(g, split_scale=1.0, perfect=False):
    h = copy.copy(g)
    if split_scale != 1.0:
        half = 0.5 * g.split * split_scale * np.sign(g.wB - g.wA)
        h.omega_all = np.array([g.mean - half, g.mean + half]); h.wA, h.wB = h.omega_all
        h.split = abs(h.wB - h.wA)
        h.E_d = np.diag([h.wA ** 2 / 2, h.wB ** 2 / 2, 0.5, 0.5]); h.Bmat = np.diag([np.sqrt(2) / h.wA, np.sqrt(2) / h.wB, np.sqrt(2), np.sqrt(2)])
        # channel projectors: same construction as Geometry with the new frequencies
        from scipy import linalg
        L_sym = np.zeros((4, 4)); H = g._H
        L_sym[2:, :2] = H / 2; L_sym[:2, 2:] = H.T / 2
        mu, xi = linalg.eigh(L_sym, h.E_d); plus = xi[:, mu > 0]
        h.P_plus = (h.E_d @ plus) @ (h.E_d @ plus).T
        pc = np.linalg.solve(h.Bmat, plus); h.plus_c = pc / np.sqrt(np.sum(pc ** 2, axis=0))
    if perfect:
        w = 0.5 * (g.W1[0, 0] + g.W1[1, 1])
        h.W1 = np.diag([w, w]); h.W0 = np.array([[0.0, np.sign(g.W0[0, 1]) * w], [np.sign(g.W0[1, 0]) * w, 0.0]])
    return h

def defect(h, eta=0.5):
    return h.run(h.calibrate(eta, h.mean), h.mean)["max_abs_defect"]

out = []
for a in (0.5, 0.3):
    spec = V.BASE.solve(a, 2, 2000)
    for tau in ((6.4, 7.7, 9.0, 10.3, 8.3) if a == 0.5 else (6.4, 7.7, 9.0, 10.2, 8.3)):
        geo = V.Sweep(a, 2, 2000, 150, spec, tau, 'sin2')
        d = doublet_only(geo)
        k = 2 * geo.keep  # recover the L-form block used for the channels
        G = geo  # H = pi m G_dd from Geometry: rebuild from P_plus is awkward, so recompute directly
        loop, sec, lam, Vm = spec
        from geometrodynamics.bulk import alice_handle_ports as ah
        Vk = Vm[:, [geo.iA, geo.iB]]
        Gm = Vk.T @ (np.repeat(loop.mass, 2)[:, None] * V.BASE.rotate_pairs(Vk))
        d._H = np.pi * 2 * Gm
        W = d.W1; W0 = d.W0
        mism = dict(wAA=float(W[0, 0]), wBB=float(W[1, 1]), w0AB=float(W0[0, 1]), w0BA=float(W0[1, 0]))
        row = dict(a=a, tau=tau, **mism,
                   doublet=defect(d), perfect_selectivity=defect(variant(d, perfect=True)),
                   split_x0p1=defect(variant(d, split_scale=0.1)), split_x0p1_perfect=defect(variant(d, split_scale=0.1, perfect=True)))
        out.append(row)
        print(a, tau, 'mismatch wAA/wBB-1=%.2e |w0|/wAA-1=%.2e' % (W[0,0]/W[1,1]-1, abs(W0[0,1])/W[0,0]-1),
              'doublet %.2e perfect %.2e split/10 %.2e split/10+perfect %.2e' % (row['doublet'], row['perfect_selectivity'], row['split_x0p1'], row['split_x0p1_perfect']), flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
