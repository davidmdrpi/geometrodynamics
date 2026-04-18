"""
Core dataclasses for the particle/wormhole-mouth model.

Each particle is a point on S³ with a wormhole mouth carrying
dynamical throat modes.  Gravitational waves propagate as expanding
shells on S³, triggering transactions when they reach the antipodal point.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np

from geometrodynamics.constants import GAMMA_MODE


@dataclass
class ThroatMode:
    """Single radial eigenmode of the wormhole throat."""

    l: int
    n: int
    omega: float
    alpha_q: float
    alpha_spin: float = 0.5
    gamma: float = GAMMA_MODE
    a: float = 0.0
    adot: float = 0.0
    phase: float = 0.0

    def step(self, drive: float, dt: float) -> None:
        acc = drive - 2.0 * self.gamma * self.adot - (self.omega ** 2) * self.a
        self.adot += dt * acc
        self.a += dt * self.adot
        self.phase = float(np.angle(self.a + 1j * self.adot / max(self.omega, 1e-9)))

    def complex_amplitude(self) -> complex:
        return complex(self.a * np.exp(1j * self.phase))


@dataclass
class MouthState:
    """State of one wormhole mouth, collecting its active throat modes."""

    orientation_sign: int
    modes: Dict[Tuple[int, int], ThroatMode] = field(default_factory=dict)
    _Q_maxwell: float = 1.0  # set externally after Maxwell solve

    def q_geom(self) -> float:
        """Geometric charge from d*F = J at the throat shell."""
        return self.orientation_sign * self._Q_maxwell * sum(
            m.alpha_q * m.a for m in self.modes.values()
        )

    def i_geom(self) -> float:
        """Time-derivative of q_geom — the geometrically derived current."""
        return self.orientation_sign * self._Q_maxwell * sum(
            m.alpha_q * m.adot for m in self.modes.values()
        )

    def mean_phase(self) -> float:
        if not self.modes:
            return 0.0
        ph = sum(np.exp(1j * m.phase) for m in self.modes.values()) / max(
            len(self.modes), 1
        )
        return float(np.angle(ph))


@dataclass
class Particle4:
    """A particle located at a point on S³ with a wormhole mouth."""

    pid: int
    p4: np.ndarray
    mouth: MouthState
    vel4: np.ndarray = field(default_factory=lambda: np.zeros(4))


@dataclass
class GravWave:
    """Expanding gravitational-wave shell on S³."""

    gid: int
    src_pid: int
    p0: np.ndarray
    t_emit: float
    radius: float = 0.0
    done: bool = False
    hit_set: set = field(default_factory=set)
    current_hits: List[Tuple[int, float, float]] = field(default_factory=list)
    triggered_pairs: set = field(default_factory=set)

    def step(self, t: float, c_gw: float) -> None:
        self.radius = min(np.pi, c_gw * (t - self.t_emit))
        self.current_hits = []
        if self.radius >= np.pi:
            self.done = True
