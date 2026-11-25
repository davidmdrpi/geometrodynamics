"""Frozen predictions of the exact-neck splitting law (docs/klein_exact_neck_prereg.md).

The law E(a, m) replaces the WKB step through the neck in the Herring law H of
docs/klein_splitting_law_prereg.md by the exact solution of the neck equation.
It keeps H's band state, frequency and normalisation. This module computes no
splitting with the finite-element solver; the held-out measurements are run
only after the freeze is published.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gammaln, lpmv

_spec = importlib.util.spec_from_file_location(
    "klein_splitting_law_predictions", Path(__file__).with_name("klein_splitting_law_predictions.py"))
HERRING = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(HERRING)

# Every cell whose splitting has already been measured: #355's archive (six)
# and #362's archive (six gated, five diagnostic).
IN_SAMPLE = ([(0.3, m) for m in (1, 2, 3)] + [(0.5, m) for m in (1, 2, 3)]
             + [(0.2, 1), (0.2, 2), (0.4, 1), (0.4, 2), (0.4, 3), (0.5, 4)]
             + [(0.6, m) for m in range(1, 6)])
GATED = ([(0.25, m) for m in (1, 2)] + [(0.35, m) for m in (1, 2, 3)]
         + [(0.45, m) for m in (1, 2, 3, 4)] + [(0.55, m) for m in (1, 2, 3, 4, 5)])
DIAGNOSTIC = [(a, m) for a in (0.65, 0.7) for m in range(1, 6)]
RTOL = 1e-13


def bound(a):
    """G1 tolerance on |r - 1|: twice the in-sample envelope, plus a resolution floor."""
    return 0.2 * a ** 4 + 1e-3


def neck_state(a, m, omega, rtol=RTOL):
    """The neck solution D that vanishes at the far rim, normalised to 1 at the near rim.

    In the conformal coordinate v of the neck (s = b sinh v, rho = b cosh v,
    v in [-L, L] with L = asinh(cot a)) the sector equation is exactly
    D'' = (m^2 - omega^2 b^2 cosh^2 v) D. Returns D(0) D'(0) and the log
    derivative D'(-L)/D(-L), with D(-L) = 1.
    """
    b, half = np.sin(a) ** 2, np.arcsinh(1 / np.tan(a))
    rhs = lambda v, y: [y[1], (m * m - omega * omega * b * b * np.cosh(v) ** 2) * y[0]]
    opts = dict(method="DOP853", rtol=rtol, atol=1e-300, first_step=1e-4)
    mid = solve_ivp(rhs, (half, 0.0), [0.0, -1.0], **opts).y[:, -1]
    rim = solve_ivp(rhs, (half, -half), [0.0, -1.0], **opts).y[:, -1]
    return float(mid[0] * mid[1] / rim[0] ** 2), float(rim[1] / rim[0])


def prediction(a, m, n=0, rtol=RTOL):
    """|omega_+ - omega_-| of the n-th Klein pair, from the band state P_{m+n}^m."""
    ell = m + n
    omega = np.sqrt(ell * (ell + 1))
    psi_rim = lpmv(m, ell, np.cos(a))
    norm = 2 * np.exp(gammaln(ell + m + 1) - gammaln(ell - m + 1)) / (2 * ell + 1)
    flux, log_derivative = neck_state(a, m, omega, rtol)
    delta = 4 * psi_rim ** 2 * abs(flux) / (omega * norm)
    out = dict(a=a, m=m, n=n, omega=float(omega), predicted_delta=float(delta),
               predicted_split=float(delta / omega), neck_log_derivative=log_derivative,
               band_log_derivative=float(-m * np.cos(a)) if n == 0 else None)
    if n == 0:
        out["herring_delta"] = HERRING.prediction(a, m)["predicted_delta"]
        out["bound"] = bound(a)
    return out


if __name__ == "__main__":
    print(json.dumps({name: [prediction(a, m) for a, m in points] for name, points in
                      (("in_sample", IN_SAMPLE), ("gated", GATED),
                       ("diagnostic", DIAGNOSTIC))}, indent=1))
