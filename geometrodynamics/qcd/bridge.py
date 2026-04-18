"""
Bridge field (double-well order parameter for string breaking)
and Cornell potential utilities.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import Callable

import numpy as np

from geometrodynamics.qcd.constants import (
    BRIDGE_ALPHA,
    BRIDGE_BETA,
    BRIDGE_GAMMA,
    BRIDGE_THRESHOLD,
    SIGMA_FM,
    A_COULOMB,
    HBAR_C,
    L_BREAK_LAT,
    R_HAD_FM,
)


@dataclass
class BridgeField:
    """Double-well order parameter η for string breaking.

    V(η) = ½α η² + ¼β η⁴,  α < 0 ⟹ minima at η = ±√(−α/β).
    Bridge "breaks" when |η| > BRIDGE_THRESHOLD.
    """

    eta: float = 0.0
    etadot: float = 0.0
    alpha: float = BRIDGE_ALPHA
    beta: float = BRIDGE_BETA
    gamma: float = BRIDGE_GAMMA
    g: float = 0.30
    broken: bool = False
    cornell_drive_scale: float = 0.0

    def dV(self) -> float:
        return self.alpha * self.eta + self.beta * self.eta ** 3

    def step(self, psi2: float, dt: float, v_cornell: float = 0.0) -> None:
        if self.broken:
            return
        acc = (
            -self.dV()
            - self.gamma * self.etadot
            + self.g * (psi2 + self.cornell_drive_scale * v_cornell)
        )
        self.etadot += dt * acc
        self.eta += dt * self.etadot

    @property
    def above_threshold(self) -> bool:
        return abs(self.eta) > BRIDGE_THRESHOLD

    @property
    def minima(self) -> float:
        return sqrt(-self.alpha / self.beta) if self.alpha < 0 else 0.0

    def energy(self) -> float:
        """E_η = ½η̇² + ½αη² + ¼βη⁴."""
        return (
            0.5 * self.etadot ** 2
            + 0.5 * self.alpha * self.eta ** 2
            + 0.25 * self.beta * self.eta ** 4
        )


# ── Cornell potential ────────────────────────────────────────────────────────

def cornell_static_energy(
    L: float,
    sigma_fm: float = SIGMA_FM,
    A_c: float = A_COULOMB,
    hbar_c: float = HBAR_C,
) -> float:
    """Cornell static energy: V(L) = σL − A·ℏc/L."""
    return sigma_fm * L - A_c * hbar_c / (L + 1e-12)


def cornell_equilibrium_amplitude(
    L: float,
    sigma_fm: float = SIGMA_FM,
    A_c: float = A_COULOMB,
    hbar_c: float = HBAR_C,
    r_had: float = R_HAD_FM,
    l_break: float = L_BREAK_LAT,
) -> float:
    """Cornell amplitude map with hadronic anchor and no hard ceiling."""
    Vc = max(0.0, cornell_static_energy(L, sigma_fm, A_c, hbar_c))
    Vb = max(1e-12, cornell_static_energy(l_break, sigma_fm, A_c, hbar_c))
    amp_had = 0.26
    amp_break = 0.98
    tail_gain = 0.42
    tail_pow = 0.55
    amp_cap = 3.0
    if L <= r_had:
        amp = amp_had
    elif L <= l_break:
        x = (L - r_had) / (l_break - r_had + 1e-12)
        s = x * x * (3.0 - 2.0 * x)
        amp = amp_had + (amp_break - amp_had) * s
    else:
        y = max(Vc / Vb, 1.0)
        amp = amp_break + tail_gain * (y ** tail_pow - 1.0)
    return float(np.clip(amp, 0.0, amp_cap))


def make_cornell_branch_potential(
    length: float,
    sigma_fm: float = SIGMA_FM,
    A_c: float = A_COULOMB,
    hbar_c: float = HBAR_C,
) -> Callable[[float, float], float]:
    """Local Cornell energy density on a branch."""
    eps = 0.03 * length + 1e-6

    def U(s: float, L: float) -> float:
        U_lin = sigma_fm / max(L, 1e-12)
        dL = max(s, eps)
        dR = max(L - s, eps)
        U_c = -0.5 * A_c * hbar_c * (1.0 / dL + 1.0 / dR) / max(L, 1e-12)
        return U_lin + U_c

    return U
