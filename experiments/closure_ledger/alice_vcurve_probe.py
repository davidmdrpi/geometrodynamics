"""Window-length and envelope sweeps of the Alice-handle exchange under detector depletion.

Reuses the geometry, detector, readings and defect of alice_exchange_probe.py
(#370) with two parameters made explicit: the window length tau and the
switching profile. The defect is the largest |u + v - 1| over all doublet states.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("alice_exchange_probe", Path(__file__).with_name("alice_exchange_probe.py"))
BASE = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BASE)


def profile(kind, t, tau):
    """Switching envelopes on 0 < t < tau, each with peak 1."""
    x = np.asarray(t, dtype=float) / tau
    inside = (x > 0) & (x < 1)
    if kind == "sin2":                       # C^1: the second derivative jumps at the edges
        return np.where(inside, np.sin(np.pi * x) ** 2, 0.0)
    if kind == "sin4":                       # C^3
        return np.where(inside, np.sin(np.pi * x) ** 4, 0.0)
    if kind == "bump":                       # C^infinity, compact support
        y = 2 * x - 1
        with np.errstate(divide="ignore", over="ignore"):
            val = np.exp(1 - 1 / np.where(inside, 1 - y * y, 1.0))
        return np.where(inside, val, 0.0)
    raise ValueError(kind)


class Sweep(BASE.Geometry):
    """BASE.Geometry with the window length and profile as parameters."""

    def __init__(self, a, m, n, keep, spectrum, tau, kind):
        super().__init__(a, m, n, keep, spectrum)
        self.tau, self.kind = float(tau), kind

    def window(self, g0, rotation):
        k, om2, tau = self.keep, self.omega_all ** 2, self.tau

        def rhs(t, y):
            Y = y[:8 * k].reshape(2 * k, 4)
            q, p = Y[:k], Y[k:]
            g = g0 * float(profile(self.kind, t, tau))
            force = -g * (self.W1 @ p + rotation * (self.W0 @ q))
            dD = -0.5 * (p.T @ force + force.T @ p)
            return np.concatenate([np.vstack([p, -om2[:, None] * q + force]).ravel(), dD.ravel()])
        y0 = np.concatenate([self.basis().ravel(), np.zeros(16)])
        sol = solve_ivp(rhs, (0.0, tau), y0, method="DOP853", rtol=1e-11, atol=1e-13)
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1]
        return y[:8 * k].reshape(2 * k, 4), y[8 * k:].reshape(4, 4)

    def calibrate(self, eta, rotation):
        """Smallest g0 reaching depletion eta; the strength scales as 1/tau at fixed eta."""
        f = lambda g: self.absorbed_plus(g, rotation) - eta
        lo, g = 0.0, 0.05 * 8.0 / self.tau
        while f(g) < 0:
            lo, g = g, 1.5 * g
            if g > 1e4:
                raise RuntimeError("depletion not reached")
        return brentq(f, lo, g, xtol=1e-14, rtol=1e-11)

    def run(self, g0, rotation, split_for_timing=None):
        split = self.split if split_for_timing is None else split_for_timing
        Y, D = self.window(g0, rotation)
        balance = float(np.max(np.abs(self.energy_form(Y) + D - np.eye(4))))
        centre = self.tau / 2
        S_u = self.doublet(self.free(self.basis(), centre))
        U = S_u.T @ self.P_plus @ S_u
        Y2 = self.free(Y, centre + np.pi / split - self.tau)
        S_v = self.doublet(Y2)
        Vf = S_v.T @ self.P_plus @ S_v
        sym = lambda X: 0.5 * (X + X.T)
        defect = np.linalg.eigvalsh(sym(U + Vf - np.eye(4)))
        leak = np.linalg.eigvalsh(sym(self.energy_form(Y2) - S_v.T @ self.E_d @ S_v))
        return dict(max_abs_defect=float(np.max(np.abs(defect))), leakage_max=float(leak.max()),
                    balance=balance, U=U.tolist(), V=Vf.tolist())


def point(spectrum, a, m, n, keep, tau, kind, eta):
    geo = Sweep(a, m, n, keep, spectrum, tau, kind)
    g0 = geo.calibrate(eta, geo.mean)
    res = geo.run(g0, geo.mean)
    res.update(a=a, m=m, N=n, K=keep, tau=tau, profile=kind, eta=eta, g0=g0,
               split=float(geo.split), split_times_tau=float(geo.split * tau), mean=float(geo.mean))
    return res
