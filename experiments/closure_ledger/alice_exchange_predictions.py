"""Frozen two-level predictions for docs/alice_exchange_prereg.md. Calls no solver.

The lowest Klein doublet of #355's handle, written in its circulating basis
(z for +m, w for -m, energy-normalized, rotating frame at the mean frequency),
obeys dz/dt = i(D/2)w - (g(t)/2)z, dw/dt = i(D/2)z when a detector absorbs
from the +m channel only, at energy rate g(t) = g0 sin^2(pi t/tau) on 0<t<tau.
D is the doublet splitting, taken from the exact-neck law of
klein_exact_neck_predictions.py. With no detector, the evolution over pi/D
exchanges z and w, so the +m intensity read at the window centre (u) and at
pi/D later (v) satisfy u + v = 1 for every state. This module computes the
largest |u + v - 1| over all doublet states when the detector is on.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

_spec = importlib.util.spec_from_file_location(
    "klein_exact_neck_predictions", Path(__file__).with_name("klein_exact_neck_predictions.py"))
NECK = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(NECK)

CELLS = [(0.5, 1), (0.4, 1), (0.3, 1), (0.3, 2)]
TAU = 8.0
ETAS = (0.1, 0.5, 0.9)
STRIP_HALF_WIDTH = 0.35


def window(t):
    return np.where((t > 0) & (t < TAU), np.sin(np.pi * t / TAU) ** 2, 0.0)


def _evolve(state, g0, split, t0, t1):
    """state: complex (z, w) columns, shape (2, k). Rotating-frame two-level flow."""
    def rhs(t, y):
        z, w = y[:len(y) // 2], y[len(y) // 2:]
        g = g0 * window(t)
        return np.concatenate([0.5j * split * w - 0.5 * g * z, 0.5j * split * z])
    y0 = np.concatenate([state[0], state[1]]).astype(complex)
    sol = solve_ivp(rhs, (t0, t1), y0, method="DOP853", rtol=1e-12, atol=1e-14)
    y = sol.y[:, -1]
    k = len(y) // 2
    return np.stack([y[:k], y[k:]])


def absorbed_pure_plus(g0, split):
    end = _evolve(np.array([[1.0], [0.0]]), g0, split, 0.0, TAU)
    return float(1 - np.sum(np.abs(end) ** 2))


def calibrate(eta, split):
    return brentq(lambda g: absorbed_pure_plus(g, split) - eta, 1e-8, 50.0, xtol=1e-14)


def _free(state, split, t):
    c, s = np.cos(split * t / 2), np.sin(split * t / 2)
    z, w = state
    return np.stack([c * z + 1j * s * w, 1j * s * z + c * w])


def defect(eta, split):
    """Largest |u + v - 1| over normalized doublet states, as an eigenvalue problem."""
    g0 = calibrate(eta, split)
    # Real basis of C^2 = R^4: Re z, Im z, Re w, Im w.
    basis = np.array([[1, 1j, 0, 0], [0, 0, 1, 1j]], dtype=complex)
    centre = TAU / 2
    u_states = _free(basis, split, centre)
    v_states = _free(_evolve(basis, g0, split, 0.0, TAU), split, centre + np.pi / split - TAU)
    def form(vec):
        a = np.stack([vec.real, vec.imag], axis=0).reshape(2, 4)  # Re/Im of first component
        return a.T @ a
    Q = form(u_states[0]) + form(v_states[0]) - np.eye(4)
    vals = np.linalg.eigvalsh(0.5 * (Q + Q.T))
    return dict(g0=g0, max_abs_defect=float(np.max(np.abs(vals))),
                min_defect=float(vals.min()), max_defect=float(vals.max()))


def predictions():
    out = []
    for a, m in CELLS:
        split = NECK.prediction(a, m)["predicted_delta"]
        row = dict(a=a, m=m, predicted_split=split, split_times_tau=split * TAU,
                   exchange_time=float(np.pi / split), selective={})
        for eta in ETAS:
            row["selective"][str(eta)] = defect(eta, split)
        out.append(row)
    return out


if __name__ == "__main__":
    print(json.dumps(predictions(), indent=1))
