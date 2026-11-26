"""Two-part predictor for the exchange defect versus window length and envelope. Calls no solver.

  E_pred(tau, profile, eta) = A(eta) |F(2 w; tau, profile)| + c_x * E_two_level(Delta, tau, profile, eta)

* F is the envelope's Fourier transform at twice the doublet frequency,
  normalised by its area: the counter-rotating coupling between the +m
  channel and the conjugate of the -m channel.
* E_two_level is the rotating-frame two-level defect of
  alice_exchange_predictions.py with the same window and depletion.
* w = sqrt(m(m+1)) and Delta comes from the exact-neck law.
* A(eta) and c_x were fitted on the in-sample cells (0.3, 1) and (0.3, 2) and
  are frozen with docs/alice_vcurve_prereg.md.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

_spec = importlib.util.spec_from_file_location(
    "klein_exact_neck_predictions", Path(__file__).with_name("klein_exact_neck_predictions.py"))
NECK = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(NECK)


def profile(kind, t, tau):
    x = np.asarray(t, dtype=float) / tau
    inside = (x > 0) & (x < 1)
    if kind == "sin2":
        return np.where(inside, np.sin(np.pi * x) ** 2, 0.0)
    if kind == "sin4":
        return np.where(inside, np.sin(np.pi * x) ** 4, 0.0)
    if kind == "bump":
        y = 2 * x - 1
        with np.errstate(divide="ignore", over="ignore"):
            val = np.exp(1 - 1 / np.where(inside, 1 - y * y, 1.0))
        return np.where(inside, val, 0.0)
    raise ValueError(kind)


def spectrum(kind, tau, omega):
    """|∫ f(t) e^{i omega t} dt| / ∫ f(t) dt over the window."""
    f = lambda t: float(profile(kind, t, tau))
    re = quad(lambda t: f(t) * np.cos(omega * t), 0, tau, limit=1000)[0]
    im = quad(lambda t: f(t) * np.sin(omega * t), 0, tau, limit=1000)[0]
    return float(np.hypot(re, im) / quad(f, 0, tau, limit=1000)[0])


def _evolve(state, g0, split, tau, kind, t0, t1):
    def rhs(t, y):
        z, w = y[:len(y) // 2], y[len(y) // 2:]
        g = g0 * float(profile(kind, t, tau))
        return np.concatenate([0.5j * split * w - 0.5 * g * z, 0.5j * split * z])
    y0 = np.concatenate([state[0], state[1]]).astype(complex)
    y = solve_ivp(rhs, (t0, t1), y0, method="DOP853", rtol=1e-12, atol=1e-14).y[:, -1]
    k = len(y) // 2
    return np.stack([y[:k], y[k:]])


def _free(state, split, t):
    c, s = np.cos(split * t / 2), np.sin(split * t / 2)
    return np.stack([c * state[0] + 1j * s * state[1], 1j * s * state[0] + c * state[1]])


def two_level(split, tau, kind, eta):
    plus = np.array([[1.0], [0.0]])
    absorbed = lambda g: 1 - np.sum(np.abs(_evolve(plus, g, split, tau, kind, 0.0, tau)) ** 2)
    lo, g = 0.0, 0.05 * 8.0 / tau
    while absorbed(g) < eta:
        lo, g = g, 1.5 * g
    g0 = brentq(lambda x: absorbed(x) - eta, lo, g, xtol=1e-14)
    basis = np.array([[1, 1j, 0, 0], [0, 0, 1, 1j]], dtype=complex)
    u = _free(basis, split, tau / 2)
    v = _free(_evolve(basis, g0, split, tau, kind, 0.0, tau), split, tau / 2 + np.pi / split - tau)
    form = lambda vec: np.stack([vec.real, vec.imag]).reshape(2, 4).T @ np.stack([vec.real, vec.imag]).reshape(2, 4)
    Q = form(u[0]) + form(v[0]) - np.eye(4)
    return float(np.max(np.abs(np.linalg.eigvalsh(0.5 * (Q + Q.T)))))


# Frozen coefficients, fitted on the in-sample cells (0.3, 1) and (0.3, 2).
A = {0.5: 0.33, 0.9: 0.86}
C_X = {0.5: 1.3, 0.9: 1.1}

# Held-out design. Both cells have m = 2, and neither has been run with a detector.
A1, A2, M = 0.5, 0.35, 2
ETAS = (0.5, 0.9)
OMEGA = np.sqrt(M * (M + 1))
# V-curve grid on sin^2 sidelobe peaks, 2*OMEGA*tau = 2*pi*(k + 1/2), roughly doubling.
KS = (2.5, 5.5, 11.5, 24.5, 49.5, 99.5, 199.5)
GRID = tuple(float(k * np.pi / OMEGA) for k in KS)
# Spectral-null scan at A1: sin^2 nulls at tau_k = k*pi/OMEGA for k = 5..8.
NULL_KS = (5, 6, 7, 8)
NULL_SCAN = tuple(float(t) for t in np.round(np.arange(6.0, 11.0001, 0.1), 10))
# Envelope control at A2, in the switching-dominated regime.
ENVELOPE_TAUS = (8.0, 16.0)
PROFILES = ("sin2", "sin4", "bump")


def predict(a, m, tau, kind, eta):
    omega = np.sqrt(m * (m + 1))
    split = NECK.prediction(a, m)["predicted_delta"]
    sw = A[eta] * spectrum(kind, tau, 2 * omega)
    ex = C_X[eta] * two_level(split, tau, kind, eta)
    return dict(a=a, m=m, tau=tau, profile=kind, eta=eta, switching=sw, exchange=ex,
                total=sw + ex, split_times_tau=split * tau)


def predictions():
    vcurve = [predict(a, M, tau, "sin2", eta) for a in (A1, A2) for eta in ETAS for tau in GRID]
    envelope = [predict(A2, M, tau, kind, 0.5) for tau in ENVELOPE_TAUS for kind in PROFILES]
    nulls = [float(k * np.pi / OMEGA) for k in NULL_KS]
    return dict(vcurve=vcurve, envelope=envelope, null_taus=nulls)


if __name__ == "__main__":
    print(json.dumps(predictions(), indent=1))
