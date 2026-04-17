"""
Spectral analysis classes for the geometrodynamic QCD model.

MöbiusSpectrum      — half-integer mode spectrum from Möbius topology
AlphaSStringAnsatz  — string-tension effective coupling diagnostic
ThroatBranchCrosswalk — 5D throat ↔ 1D branch mass matching
"""

from __future__ import annotations

from math import pi, sqrt, isfinite, nan
from typing import Dict, List

import numpy as np
from scipy.optimize import brentq as _brentq
from scipy.linalg import eig as _scipy_eig

from geometrodynamics.qcd.constants import HBAR_C, R_MID, R_OUTER


class MobiusSpectrum:
    """Half-integer mode spectrum of a Möbius flux tube."""

    def __init__(self, L=1.0, v=1.0, N=150):
        self.L = L
        self.v = v
        self.N = N

    def analytic_modes(self, n_max=5):
        return [
            {
                "n": n,
                "omega_mob": (n + 0.5) * pi * self.v / self.L,
                "omega_ord": (n + 1) * pi * self.v / self.L,
                "ratio": (n + 0.5) / (n + 1),
            }
            for n in range(n_max)
        ]

    def numerical_modes(self, n_max=5) -> np.ndarray:
        N = self.N
        ds = self.L / (N - 1)
        v2 = self.v ** 2
        Nin = N - 1
        d = np.ones(Nin) * (2 * v2 / ds ** 2)
        o = np.ones(Nin - 1) * (-v2 / ds ** 2)
        H = np.diag(d) + np.diag(o, 1) + np.diag(o, -1)
        H[-1, -1] = v2 / ds ** 2
        ev = np.linalg.eigvalsh(H)
        pos = np.sort(ev[ev > 0])[:n_max]
        return np.sqrt(pos)

    def report(self, n_max=4) -> Dict:
        an = self.analytic_modes(n_max)
        num = self.numerical_modes(n_max)
        checks = [
            {
                "n": k,
                "omega_num": float(num[k]),
                "omega_mob": an[k]["omega_mob"],
                "ratio": an[k]["ratio"],
                "pass": abs(num[k] - an[k]["omega_mob"]) / (an[k]["omega_mob"] + 1e-12) < 0.05,
            }
            for k in range(min(n_max, len(num)))
        ]
        return {
            "checks": checks,
            "all_pass": all(c["pass"] for c in checks),
            "zero_point_GeV": HBAR_C * (0.5 * pi * self.v / self.L),
        }


class ThroatBranchCrosswalk:
    """5D throat ↔ 1D branch mass matching with Q_maxwell calibration."""

    def __init__(self, R_mid=R_MID, delta=0.26, N_cheb=55):
        self.R_mid = R_mid
        self.R_outer = R_mid + delta
        self.N_cheb = N_cheb
        self._Q_maxwell = None

    def _rstar(self, r):
        rs = self.R_mid
        return r + rs / 2 * np.log(abs((r - rs) / (r + rs) + 1e-15))

    def _r_of_rstar(self, rstar):
        rs = self.R_mid
        def f(r):
            if r <= rs:
                return -1e30
            return r + rs / 2 * np.log(abs((r - rs) / (r + rs) + 1e-30)) - rstar
        try:
            return float(_brentq(f, rs + 1e-8, max(abs(rstar) + rs + 5, 2 * rs), xtol=1e-10))
        except Exception:
            return rs + 1e-6

    def _V(self, r, l):
        rs = self.R_mid
        f = 1 - (rs / r) ** 2
        return f * (l * (l + 2) / r ** 2 + 3 * rs ** 2 / r ** 4)

    def _cheb(self, N):
        x = np.cos(pi * np.arange(N + 1) / N)
        c = np.ones(N + 1)
        c[0] = c[N] = 2.0
        c *= (-1) ** np.arange(N + 1)
        X = np.tile(x, (N + 1, 1))
        dX = X - X.T
        D = (c[:, None] / c[None, :]) / (dX + np.eye(N + 1))
        D -= np.diag(np.sum(D, axis=1))
        return x, D

    def _slope(self, u, rg, n_fit=8):
        dr = rg[1:n_fit] - self.R_mid
        return float(np.dot(dr, u[1:n_fit]) / np.dot(dr, dr))

    def _solve(self, l):
        N = self.N_cheb
        rs_min = self._rstar(self.R_mid + 5e-4)
        rs_max = self._rstar(self.R_outer - 5e-4)
        _, D = self._cheb(N)
        D2 = D @ D
        Lscl = (rs_max - rs_min) / 2.0
        rsg = rs_min + Lscl * (1 - np.cos(pi * np.arange(N + 1) / N))
        rg = np.array([self._r_of_rstar(rs) for rs in rsg])
        Vg = self._V(rg, l)
        H = -(1 / Lscl ** 2) * D2 + np.diag(Vg)
        H_int = H[1:N, 1:N]
        ev, evec = _scipy_eig(H_int)
        ev = np.real(ev)
        pos = np.where(ev > 0)[0]
        idx = np.argsort(ev[pos])[0]
        u = np.zeros(N + 1)
        u[1:N] = np.real(evec[:, pos[idx]])
        if abs(u.min()) > u.max():
            u = -u
        u /= abs(u).max() + 1e-12
        return float(np.sqrt(ev[pos[idx]])), u, rg

    def solve_modes(self, l_values=(1, 3, 5, 7)) -> Dict[int, Dict]:
        results = {}
        ref_sl = None
        for l in l_values:
            om, u, rg = self._solve(l)
            sl = self._slope(u, rg)
            if l == 1:
                ref_sl = abs(sl)
            aq = sl / (ref_sl if ref_sl else 1.0)
            results[l] = {
                "l": l, "omega": float(om),
                "M_geom": float(HBAR_C * om),
                "slope": float(sl), "alpha_q": float(aq),
            }
        return results

    def calibrate_Q_maxwell(self, L_eq: float, v: float = 1.0) -> Dict:
        om_l1, u_l1, rg_l1 = self._solve(1)
        M_branch_n1 = HBAR_C * pi * v / L_eq
        M_throat_bare = HBAR_C * om_l1
        Q_calibrated = M_branch_n1 / (M_throat_bare + 1e-30)
        slope_l1 = abs(self._slope(u_l1, rg_l1, n_fit=10))
        Q_v22 = self.R_mid ** 2 * slope_l1 / (4.0 * pi)
        self._Q_maxwell = Q_calibrated
        return {
            "Q_calibrated": Q_calibrated,
            "Q_v22": Q_v22,
            "ratio_cal_v22": Q_calibrated / (Q_v22 + 1e-30),
            "M_branch_n1": M_branch_n1,
            "M_throat_bare_l1": M_throat_bare,
            "omega_l1_bare": float(om_l1),
        }
