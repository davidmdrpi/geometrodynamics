"""Four-term predictor for the exchange defect, including the residual at the sin^2 spectral nulls.
Calls no solver.

  E = A(eta)|F(2w)| + c_x(eta) E_2L + C(eta)/(w tau)^3 + K(eta) Delta/w

The first two terms, and their frozen coefficients, are those of
alice_vcurve_predictions.py (freeze 291ab5f). Two terms are added; both live
inside the 4-dimensional doublet:

* intrinsic: second order in the drain rate, falling as (w tau)^-3. It
  remains with the strip profiles made identical.
* mismatch: the detector strip sees the A and B mode profiles differ by
  eps ~ 3.9 Delta/w, so it absorbs a little from the -m channel. It is
  roughly independent of tau.

The new coefficients C(eta) and K(eta) were fitted on in-sample cells only:
(0.3, 2) for C, and (0.5, 2) for K. They are frozen with docs/alice_nulls_prereg.md.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

_spec = importlib.util.spec_from_file_location("alice_vcurve_predictions", Path(__file__).with_name("alice_vcurve_predictions.py"))
VC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(VC)
NECK = VC.NECK

C = {0.5: 0.51, 0.9: 3.95}
K = {0.5: 0.26, 0.9: 0.42}

M = 2
CELLS = (0.45, 0.55)                 # held out: never run with a detector
ETAS = (0.5, 0.9)
NULL_KS = (5, 8, 13, 20)             # tau_k = k pi / w, sin^2 nulls of F(2w)
BOUNDARY_TAUS = (2.5, 3.21, 4.0)     # reachable-depletion map
BOUNDARY_IN_SAMPLE = {2.5: 0.797, 3.21: 0.877, 4.0: 0.943}   # eta_max, sin^2, at (0.5, 2) and (0.35, 2)
DOMAIN_RATE = 0.5                    # validity: the peak drain rate must not exceed 0.5 w


def terms(a, m, tau, kind, eta, omega=None):
    omega = np.sqrt(m * (m + 1)) if omega is None else omega
    split = NECK.prediction(a, m)["predicted_delta"]
    sw = VC.A[eta] * VC.spectrum(kind, tau, 2 * omega)
    ex = VC.C_X[eta] * VC.two_level(split, tau, kind, eta)
    intrinsic = C[eta] / (omega * tau) ** 3
    mismatch = K[eta] * split / omega
    return dict(a=a, m=m, tau=float(tau), profile=kind, eta=eta, switching=sw, exchange=ex,
                intrinsic=intrinsic, mismatch=mismatch, total=sw + ex + intrinsic + mismatch)


def predictions():
    omega = np.sqrt(M * (M + 1))
    rows = [terms(a, M, k * np.pi / omega, "sin2", eta) for a in CELLS for eta in ETAS for k in NULL_KS]
    return dict(nulls=rows, boundary=BOUNDARY_IN_SAMPLE, domain_rate=DOMAIN_RATE)


if __name__ == "__main__":
    print(json.dumps(predictions(), indent=1))
