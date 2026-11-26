"""Post-result check: how the defect degrades when the medium's rotation misses the doublet frequency.

The dashpot force is -g W (xdot + s*wbar*C x); s = 1 is the tuned selective detector of #370.
Expected wrong-channel/right-channel drain-rate ratio r = ((1 - s)/(1 + s))^2, so the -m channel
loses 1 - (1 - eta)^r, which bounds the defect at leading order. Doublet only, identical strip
profiles, at an exact sin^2 null, so switching, mismatch and most of the intrinsic term are absent.
"""
import importlib.util, json
from pathlib import Path
import numpy as np

here = Path(__file__).resolve()
_s = importlib.util.spec_from_file_location("run", here.parents[3] / "alice_nulls_run.py")
R = importlib.util.module_from_spec(_s); _s.loader.exec_module(R)
V = R.V

rows = []
for a, k in ((0.3, 13), (0.45, 13)):
    spec = V.BASE.solve(a, 2, 2000)
    g = V.Sweep(a, 2, 2000, 150, spec, 8.0, "sin2")
    tau = k * np.pi / g.mean
    d = R.identical_profiles(R.doublet_only(V.Sweep(a, 2, 2000, 150, spec, tau, "sin2")))
    for eta in (0.5, 0.9):
        for s in (0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 2.0):
            rot = s * d.mean
            g0 = d.calibrate(eta, rot)
            res = d.run(g0, rot)
            r = ((1 - s) / (1 + s)) ** 2
            rows.append(dict(a=a, k=k, eta=eta, s=s, defect=res["max_abs_defect"], ratio_r=r,
                             wrong_channel_loss=float(1 - (1 - eta) ** r), g0_over_wbar=g0 / d.mean))
            print(a, k, eta, "s %.2f defect %.3e  1-(1-eta)^r %.3e  g0/wbar %.3f" % (s, res["max_abs_defect"], 1 - (1 - eta) ** r, g0 / d.mean), flush=True)
json.dump(rows, open(here.with_name("rotation_detuning.json"), "w"), indent=1)
