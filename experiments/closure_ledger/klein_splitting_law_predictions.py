"""Frozen predictions of the Klein splitting law (docs/klein_splitting_law_prereg.md).

Evaluates the closed-form prediction H(a, m) for the lowest Klein pair of the
Alice-handle loop of PR #355. It computes no splitting with the solver; the
held-out measurements are run only after the freeze is published.
"""
import json

import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln

IN_SAMPLE = [(0.3, 1), (0.3, 2), (0.3, 3), (0.5, 1), (0.5, 2), (0.5, 3)]
GATED = [(0.2, 1), (0.2, 2), (0.4, 1), (0.4, 2), (0.4, 3), (0.5, 4)]
DIAGNOSTIC = [(0.6, m) for m in range(1, 6)]


def neck_action(a, m, omega):
    """W_neck = integral over the neck of sqrt(m^2/rho^2 - omega^2) ds."""
    b, half = np.sin(a) ** 2, np.sin(a) * np.cos(a)
    integrand = lambda s: np.sqrt(m * m / (s * s + b * b) - omega * omega)
    return quad(integrand, -half, half, limit=400, epsabs=0, epsrel=1e-13)[0]


def prediction(a, m):
    """|omega_+ - omega_-| of the lowest Klein pair, n = 0."""
    omega = np.sqrt(m * (m + 1))
    kappa_rim = np.sqrt(m * m - omega * omega * np.sin(a) ** 2)
    norm = np.exp(0.5 * np.log(np.pi) + gammaln(m + 1) - gammaln(m + 1.5))
    w = neck_action(a, m, omega)
    delta = 4 * kappa_rim * np.sin(a) ** (2 * m) * np.exp(-w) / (omega * norm)
    return dict(a=a, m=m, omega=float(omega), neck_action=float(w),
                kappa_rim=float(kappa_rim), predicted_delta=float(delta),
                predicted_split=float(delta / omega))


if __name__ == "__main__":
    print(json.dumps({name: [prediction(a, m) for a, m in points] for name, points in
                      (("in_sample", IN_SAMPLE), ("gated", GATED),
                       ("diagnostic", DIAGNOSTIC))}, indent=1))
