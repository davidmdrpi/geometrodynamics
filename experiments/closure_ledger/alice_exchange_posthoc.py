"""Post hoc, after the frozen verdict of the exchange test: where the defect floor comes from.

Not part of freeze 0fb3b6e and not a gate. Varies the detector strip width and
the window length at N = 2000, K = 150, eta = 0.5, and compares the A and B
mode amplitudes inside the strip.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("alice_exchange_probe", Path(__file__).with_name("alice_exchange_probe.py"))
P = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(P)
CELLS = [(0.3, 1), (0.3, 2)]


def _set_tau(tau):
    P.PRED.TAU = P.TAU = tau


def diagnose(a, m, n=2000, keep=150):
    spec = P.solve(a, m, n)
    loop, _, _, V = spec
    g = P.Geometry(a, m, n, keep, spec)
    band = np.flatnonzero(loop.on_band)
    theta = a + loop.ell[band]
    rA = g.wA * V[:, g.iA].reshape(-1, 2)[band, 0]
    rB = g.wB * V[:, g.iB].reshape(-1, 2)[band, 1]
    strip = np.abs(theta - np.pi / 2) < P.HALF
    mismatch = np.abs(rB[strip]) / np.abs(rA[strip]) - 1
    out = dict(a=a, m=m, one_minus_overlap=1 - g.band_overlap,
               strip_amplitude_mismatch=[float(mismatch.min()), float(mismatch.max())], widths={}, windows={})
    original = P.HALF
    for half in (0.15, 0.35, 0.6, 1.0):
        P.HALF = half
        gw = P.Geometry(a, m, n, keep, spec)
        out["widths"][str(half)] = gw.run(gw.calibrate(0.5, gw.mean), gw.mean)["max_abs_defect"]
    P.HALF = original
    for tau in (4.0, 8.0, 16.0, 32.0):
        _set_tau(tau)
        gt = P.Geometry(a, m, n, keep, spec)
        r = gt.run(gt.calibrate(0.5, gt.mean), gt.mean)
        out["windows"][str(tau)] = dict(defect=r["max_abs_defect"], leakage=r["leakage_max"],
                                        split_times_tau=float(gt.split * tau))
    _set_tau(8.0)
    return out


if __name__ == "__main__":
    rows = [diagnose(a, m) for a, m in CELLS]
    path = ROOT / "experiments/closure_ledger/runs/20261125_alice_exchange/posthoc.json"
    path.write_text(json.dumps(dict(note="post hoc, after the frozen verdict; not a gate", cells=rows), indent=1) + "\n")
    print(json.dumps(rows, indent=1))
