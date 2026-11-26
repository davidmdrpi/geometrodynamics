"""Zero-parameter check of the closed form for the intrinsic null term (proposed after the CONFIRMED run).

  C_int(eta, k) = (pi^2/4) (-ln(1-eta)) (1 - sqrt(1-eta)) / (1 - k^-2),   defect = C_int / (w tau_k)^3

Evaluated in the doublet-only model with identical strip profiles (no mismatch)
and the splitting cut 100-fold and 10,000-fold (exchange removed), exactly at tau_k = k pi / w.
"""
import importlib.util, json, sys
from pathlib import Path
import numpy as np

here = Path(__file__).resolve()
ns = {}
exec(open(here.with_name("doublet_dissect.py")).read().split("out = []")[0], ns)
V, doublet_only, variant, defect = ns["V"], ns["doublet_only"], ns["variant"], ns["defect"]


def closed_form(eta, k):
    return np.pi ** 2 / 4 * (-np.log(1 - eta)) * (1 - np.sqrt(1 - eta)) / (1 - k ** -2.0)


rows = []
for a in (0.3, 0.45, 0.55):
    spec = V.BASE.solve(a, 2, 2000)
    loop, sec, lam, Vm = spec
    for k in (5, 8, 13, 20):
        g = V.Sweep(a, 2, 2000, 150, spec, 8.0, "sin2")
        tau = k * np.pi / g.mean
        g = V.Sweep(a, 2, 2000, 150, spec, tau, "sin2")
        d = doublet_only(g)
        Vk = Vm[:, [g.iA, g.iB]]
        d._H = np.pi * 2 * (Vk.T @ (np.repeat(loop.mass, 2)[:, None] * V.BASE.rotate_pairs(Vk)))
        for scale in (1e-2, 1e-4):
            clean = variant(d, split_scale=scale, perfect=True)
            for eta in (0.1, 0.25, 0.5, 0.75, 0.9):
                measured = defect(clean, eta) * (clean.mean * tau) ** 3
                rows.append(dict(a=a, k=k, eta=eta, tau=float(tau), split_scale=scale, measured=float(measured),
                                 closed_form=float(closed_form(eta, k)), ratio=float(measured / closed_form(eta, k))))
                print(a, k, eta, scale, "measured %.4f closed form %.4f ratio %.4f" % (measured, closed_form(eta, k), measured / closed_form(eta, k)), flush=True)
json.dump(rows, open(here.with_name("closed_form_check.json"), "w"), indent=1)
