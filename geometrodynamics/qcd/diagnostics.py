"""
Diagnostic classes for the QCD solver.

LatticeStringTension    — Schwinger-rate σ extraction from bridge nucleation
HybridModeShift         — loop-coupling shift of the hybrid mode spectrum
BridgeCouplingCalibrator — pre-fit utility for bridge coupling g
"""

from __future__ import annotations

from math import pi, sqrt, isfinite, nan
from typing import Dict, List, Optional

import numpy as np

from geometrodynamics.qcd.constants import (
    HBAR_C,
    SIGMA_QCD,
    SIGMA_FM,
    A_COULOMB,
    L_BREAK_LAT,
    L_STRING_ONSET_FM,
    R_HAD_FM,
    BRIDGE_ALPHA,
    BRIDGE_BETA,
    BRIDGE_GAMMA,
    BRIDGE_THRESHOLD,
    M_Q_SCHWINGER_GEV,
)
from geometrodynamics.qcd.bridge import (
    cornell_static_energy,
    cornell_equilibrium_amplitude,
)
from geometrodynamics.qcd.topology import make_meson_tube
from geometrodynamics.qcd.solver import HadronicNetworkSolver


class LatticeStringTension:
    """Lattice-style σ extraction from Schwinger bridge nucleation rate."""

    def __init__(self, L_vals=None, v=1.0, N=60, dt=0.004, t_max=80.0, bridge_g=2.0):
        self.L_vals = L_vals if L_vals is not None else list(np.linspace(1.0, 4.0, 16))
        self.v = v
        self.N = N
        self.dt = dt
        self.t_max = t_max
        self.bridge_g = bridge_g
        B_val = (pi - A_COULOMB) * HBAR_C
        self.m_q = sqrt(B_val * SIGMA_FM)

    def bridge_g_of_L(self, L: float) -> float:
        exponent = (pi * M_Q_SCHWINGER_GEV ** 2) / (SIGMA_FM * max(L, 0.1))
        tunneling_prob = np.exp(-exponent)
        V_ref = max(1e-12, cornell_static_energy(L_BREAK_LAT))
        V_c = max(0.0, cornell_static_energy(L))
        excess = max(0.0, (V_c - V_ref) / V_ref)
        cornell_boost = 1.0 + 0.55 * excess
        return self.bridge_g * 5.0 * tunneling_prob * cornell_boost

    def _t_nucleate(self, L: float) -> float:
        gL = self.bridge_g_of_L(L)
        net = make_meson_tube(L, v=self.v, N=self.N, dt=self.dt, bridge_g=gL)
        if net.branches[0].bridge is not None:
            net.branches[0].bridge.cornell_drive_scale = 0.08
        s = np.linspace(0, L, self.N)
        net.initialize_fields(psi0={0: 0.50 * np.sin(pi * s / L)})
        slv = HadronicNetworkSolver(net, antipodal_coupling=0.04)
        for _ in range(int(self.t_max / self.dt)):
            slv.step()
            if slv.nucleation_events:
                return slv.nucleation_events[0]["time"]
        return float("nan")

    def scan(self) -> Dict:
        results = []
        for L in self.L_vals:
            ampL = cornell_equilibrium_amplitude(L)
            V_c = cornell_static_energy(L)
            gL = self.bridge_g_of_L(L)
            t_nuc = self._t_nucleate(L)
            results.append({
                "L": L, "t_nuc": t_nuc,
                "Gamma": 1.0 / t_nuc if (isfinite(t_nuc) and t_nuc > 0) else nan,
                "amp_eq": ampL, "V_cornell": V_c, "bridge_g": gL,
            })
        valid = [
            (r["L"], r["Gamma"]) for r in results
            if isfinite(r["Gamma"]) and r["Gamma"] > 0 and r["L"] >= L_STRING_ONSET_FM
        ]
        sigma_eff = fit_slope = fit_r2 = nan
        if len(valid) >= 3:
            L_v = np.array([x[0] for x in valid])
            G_v = np.array([x[1] for x in valid])
            lnG = np.log(G_v + 1e-30)
            x = 1.0 / L_v
            try:
                p = np.polyfit(x, lnG, 1)
                fit_slope = float(p[0])
                pred = np.polyval(p, x)
                ss_res = float(np.sum((lnG - pred) ** 2))
                ss_tot = float(np.sum((lnG - np.mean(lnG)) ** 2))
                fit_r2 = 1.0 - ss_res / (ss_tot + 1e-30)
                sigma_eff = pi * self.m_q ** 2 / (abs(fit_slope) * HBAR_C + 1e-30)
            except Exception:
                pass
        return {
            "scan": results, "fit_slope": fit_slope,
            "sigma_eff": sigma_eff, "sigma_qcd": SIGMA_QCD,
            "ratio": sigma_eff / SIGMA_QCD if isfinite(sigma_eff) else nan,
            "fit_r2": fit_r2, "m_q_GeV": self.m_q,
        }


class HybridModeShift:
    """Hybrid mode-shift diagnostic using direct matrix eigensolver."""

    def __init__(self, meson_length=1.0, ring_circumference=1.5, v=1.0, N=120):
        self.mL = meson_length
        self.rC = ring_circumference
        self.v = v
        self.N = N
        self.La = meson_length / 2.0

    def _build_K(self, L_arm: float, kappa_junc: float = 0.0) -> np.ndarray:
        N = self.N
        ds = L_arm / (N - 1)
        v2 = self.v ** 2
        Nin = N - 1
        d = np.ones(Nin) * (2 * v2 / ds ** 2)
        o = np.ones(Nin - 1) * (-v2 / ds ** 2)
        K = np.diag(d) + np.diag(o, 1) + np.diag(o, -1)
        K[-1, -1] = v2 / ds ** 2 + kappa_junc
        return K

    def _omega_n(self, kappa_junc: float, n: int = 0) -> float:
        K = self._build_K(self.La, kappa_junc)
        ev = np.linalg.eigvalsh(K)
        pos = np.sort(ev[ev > 0])
        return float(np.sqrt(pos[n])) if len(pos) > n else nan

    def run(self) -> Dict:
        L_loop = self.rC / (2.0 * pi * self.v)
        omega_loop = 2.0 * pi * self.v / self.rC
        kappa_default = L_loop * omega_loop ** 2
        omega_free = self._omega_n(kappa_junc=0.0, n=0)
        omega_coupled = self._omega_n(kappa_junc=kappa_default, n=0)

        if isfinite(omega_free) and isfinite(omega_coupled) and omega_free > 1e-6:
            delta_omega = omega_coupled - omega_free
            delta_alpha = delta_omega / (omega_free + 1e-30)
        else:
            delta_omega = delta_alpha = nan

        return {
            "omega_free": omega_free,
            "omega_coupled": omega_coupled,
            "delta_omega": delta_omega,
            "delta_alpha_s": delta_alpha,
            "kappa_used": kappa_default,
        }


class BridgeCouplingCalibrator:
    """Pre-fit utility: find bridge coupling g that nucleates at t_target."""

    def __init__(self, L_ref=1.8, v=1.0, N=60, dt=0.005, t_max=100.0):
        self.L = L_ref
        self.v = v
        self.N = N
        self.dt = dt
        self.t_max = t_max
        self.t_target = L_BREAK_LAT / v

    def _t_nuc(self, g):
        psi = 0.8 * np.exp(
            -0.5 * ((np.linspace(0, self.L, self.N) - self.L / 2) / (self.L / 8)) ** 2
        )
        psi_old = psi.copy()
        ds = self.L / (self.N - 1)
        v2 = self.v ** 2
        eta = etadot = t = 0.0
        for _ in range(int(self.t_max / self.dt)):
            lap = np.zeros(self.N)
            lap[1:-1] = (psi[:-2] - 2 * psi[1:-1] + psi[2:]) / ds ** 2
            p_new = np.clip(
                2 * psi - psi_old + self.dt ** 2 * (v2 * lap - 0.02 * psi ** 3),
                -50, 50,
            )
            p_new[0] = p_new[-1] = 0.0
            psi_old = psi.copy()
            psi = p_new
            acc = (
                -BRIDGE_ALPHA * eta
                - BRIDGE_BETA * eta ** 3
                - BRIDGE_GAMMA * etadot
                + g * float(psi[self.N // 2] ** 2)
            )
            etadot += self.dt * acc
            eta += self.dt * etadot
            t += self.dt
            if abs(eta) > BRIDGE_THRESHOLD:
                return t
        return float("nan")

    def calibrate(self, g_min=0.1, g_max=3.0, n_pts=10) -> Dict:
        g_vals = np.linspace(g_min, g_max, n_pts)
        t_vals = [self._t_nuc(float(g)) for g in g_vals]
        t_arr = np.array(t_vals)
        g_star = float("nan")
        for k in range(len(t_arr) - 1):
            if (
                isfinite(t_arr[k])
                and isfinite(t_arr[k + 1])
                and (t_arr[k] - self.t_target) * (t_arr[k + 1] - self.t_target) < 0
            ):
                frac = (self.t_target - t_arr[k]) / (t_arr[k + 1] - t_arr[k])
                g_star = float(g_vals[k] + frac * (g_vals[k + 1] - g_vals[k]))
                break
        if not isfinite(g_star):
            fm = np.isfinite(t_arr)
            if fm.any():
                g_star = float(g_vals[fm][np.argmin(np.abs(t_arr[fm] - self.t_target))])
        t_star = self._t_nuc(g_star) if isfinite(g_star) else nan
        return {
            "g_star": g_star,
            "t_nuc_star": t_star,
            "t_target": self.t_target,
            "ratio": t_star / self.t_target if isfinite(t_star) else nan,
        }
