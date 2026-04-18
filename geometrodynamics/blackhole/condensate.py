"""
Coherent wormhole-throat condensate — the black-hole microscopic model.

Central hypothesis: a black hole of mass M is a macroscopic coherent state
of N non-orientable wormhole throats, each carrying the same mode spectrum
as a microscopic particle-antiparticle pair (ThroatMode).

The condensate locks N throat oscillators into phase, producing:
  - a macroscopic trapped surface (horizon) whose area A = 16πM²,
  - net charge Q from the orientation-sum of non-orientable throats,
  - entropy S = A/(4 l_P²) from the number of independent throat states,
  - a discrete resonance spectrum from the collective mode structure.

In geometric units (G = c = 1), with ℏ = l_P².
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from geometrodynamics.constants import GAMMA_MODE
from geometrodynamics.transaction.particles import ThroatMode


# ── Physical scales (geometric units G = c = 1) ─────────────────────────────
L_PLANCK_SQ = 1.0       # ℏ in geometric units; set to 1 for dimensionless work
A_THROAT_MIN = 4.0 * L_PLANCK_SQ / np.log(2)  # ≈ 5.77 l_P²  per throat


@dataclass
class ThroatState:
    """State of a single wormhole throat inside the condensate.

    Each throat is one microscopic wormhole — the same object that appears
    as a particle-antiparticle pair in the transaction module.

    Attributes
    ----------
    orientation : int
        +1 or -1.  The non-orientable topology assigns each throat an
        orientation sign; electric charge comes from the sum over the
        condensate.
    modes : dict
        (l, n) -> ThroatMode  mapping for active radial eigenmodes.
    phase : float
        Overall coherent phase of this throat relative to the condensate.
    """

    orientation: int = +1
    modes: Dict[Tuple[int, int], ThroatMode] = field(default_factory=dict)
    phase: float = 0.0

    def mode_energy(self) -> float:
        """Total energy stored in the throat's active modes."""
        return sum(
            0.5 * (m.adot ** 2 + m.omega ** 2 * m.a ** 2)
            for m in self.modes.values()
        )


@dataclass
class CoherentCondensate:
    """A black hole as a coherent condensate of N wormhole throats.

    Parameters
    ----------
    mass : float
        Total gravitational mass M (geometric units).
    throats : list of ThroatState
        The microscopic throat degrees of freedom.
    l_core : float or None
        Core length scale replacing the singularity.
        If None, computed from throat count: l_core ~ r_throat / N^{1/3}.
    """

    mass: float
    throats: List[ThroatState] = field(default_factory=list)
    l_core: Optional[float] = None

    # ── Derived geometric quantities ─────────────────────────────────────

    @property
    def N(self) -> int:
        """Number of throats in the condensate."""
        return len(self.throats)

    @property
    def r_horizon(self) -> float:
        """Schwarzschild horizon radius: r_H = 2M."""
        return 2.0 * self.mass

    @property
    def horizon_area(self) -> float:
        """Horizon area: A = 4π r_H² = 16π M²."""
        return 4.0 * np.pi * self.r_horizon ** 2

    @property
    def core_scale(self) -> float:
        """Core length scale replacing the singularity.

        Derived from the horizon area and throat count:

            l = 2M / √N

        With N = 4πM²/ln(2), this gives l = √(ln2/π) ≈ 0.47 l_P,
        a Planck-scale constant independent of black-hole mass.  This
        is the characteristic throat separation at the core, set by
        the area per throat a_min = A/(4πN) → l² ~ a_min.
        """
        if self.l_core is not None:
            return self.l_core
        if self.N > 0:
            return 2.0 * self.mass / np.sqrt(self.N)
        return 1.0

    # ── Charge from throat orientations ──────────────────────────────────

    @property
    def net_charge(self) -> int:
        """Net charge = sum of throat orientation signs.

        This realises Wheeler's 'charge without charge': each non-
        orientable throat carries ±1 unit of topological charge, and the
        macroscopic black-hole charge is their algebraic sum.
        """
        return sum(t.orientation for t in self.throats)

    @property
    def charge_fraction(self) -> float:
        """Fractional charge: Q/N.  A neutral BH has f → 0 for large N."""
        if self.N == 0:
            return 0.0
        return self.net_charge / self.N

    # ── Entropy from throat state counting ───────────────────────────────

    def entropy_from_area(self) -> float:
        """Bekenstein-Hawking entropy from the horizon area.

        S_BH = A / (4 l_P²).
        """
        return self.horizon_area / (4.0 * L_PLANCK_SQ)

    def entropy_from_throats(self, states_per_throat: int = 2) -> float:
        """Entropy from counting independent throat microstates.

        Each throat has `states_per_throat` configurations (default 2
        from non-orientable orientation ±1).  The entropy is
        S = N · ln(states_per_throat).

        The Bekenstein-Hawking formula is recovered when
        N = A / a_throat  and  a_throat · ln(k) = 4 l_P².
        """
        if self.N == 0:
            return 0.0
        return self.N * np.log(states_per_throat)

    def throat_count_from_area(self, states_per_throat: int = 2) -> int:
        """Number of throats predicted by matching S_BH = S_throat.

        N = A / (4 l_P² / ln(k))  =  S_BH / ln(k).
        """
        s_bh = self.entropy_from_area()
        return int(round(s_bh / np.log(states_per_throat)))

    # ── Coherent phase and mode spectrum ─────────────────────────────────

    @property
    def mean_phase(self) -> float:
        """Mean phase of the condensate (coherent order parameter)."""
        if self.N == 0:
            return 0.0
        phasor = sum(np.exp(1j * t.phase) for t in self.throats)
        return float(np.angle(phasor / self.N))

    @property
    def coherence(self) -> float:
        """Coherence fraction: |⟨e^{iφ}⟩|.  1 = fully coherent, 0 = thermal."""
        if self.N == 0:
            return 0.0
        phasor = sum(np.exp(1j * t.phase) for t in self.throats)
        return float(abs(phasor) / self.N)

    def collective_mode_spectrum(self) -> Dict[Tuple[int, int], float]:
        """Collective-mode frequencies from phase-locked throats.

        In the coherent condensate, each (l, n) mode gains a √N
        enhancement (Dicke superradiance analogue), shifting the
        collective frequency:  Ω_{l,n} = ω_{l,n} · √(1 + (N-1)·c²)
        where c is the coherence.
        """
        if self.N == 0:
            return {}
        # Collect all mode keys
        all_keys: set = set()
        for t in self.throats:
            all_keys.update(t.modes.keys())

        coh = self.coherence
        enhancement = np.sqrt(1.0 + (self.N - 1) * coh ** 2)

        spectrum = {}
        for key in sorted(all_keys):
            # Average bare frequency for this mode across throats
            omegas = [t.modes[key].omega for t in self.throats if key in t.modes]
            if omegas:
                omega_bar = float(np.mean(omegas))
                spectrum[key] = omega_bar * enhancement
        return spectrum

    def total_mode_energy(self) -> float:
        """Sum of mode energies over all throats."""
        return sum(t.mode_energy() for t in self.throats)

    # ── Temperature from surface gravity ─────────────────────────────────

    @property
    def hawking_temperature(self) -> float:
        """Hawking temperature from surface gravity.

        T_H = ℏ κ / (2π)  where κ = 1/(4M) for Schwarzschild.
        In our units (ℏ = l_P² = 1): T = 1/(8πM).

        The geometric interpretation: T is the lowest collective
        resonance frequency of the throat condensate at the horizon,
        not a particle-pair creation rate.
        """
        if self.mass <= 0:
            return float("inf")
        return 1.0 / (8.0 * np.pi * self.mass)


# ── Constructors ─────────────────────────────────────────────────────────────

def build_schwarzschild_condensate(
    mass: float,
    *,
    seed_modes: Optional[List[Tuple[int, int, float]]] = None,
    states_per_throat: int = 2,
    coherent: bool = True,
    rng: Optional[np.random.Generator] = None,
) -> CoherentCondensate:
    """Construct a Schwarzschild (Q=0) black-hole condensate.

    Populates exactly the number of throats required to reproduce
    S_BH = A/(4 l_P²) with `states_per_throat` states each.

    Parameters
    ----------
    mass : float
        Gravitational mass M.
    seed_modes : list of (l, n, omega) or None
        Bare mode spectrum to seed into every throat.
        Default: first three Tangherlini modes at l=1.
    states_per_throat : int
        Number of microstates per throat (default 2 for ±orientation).
    coherent : bool
        If True, all throat phases are aligned (zero temperature limit).
        If False, phases are randomised (thermal state).
    rng : Generator or None
        Random number generator for stochastic initialisation.

    Returns
    -------
    CoherentCondensate
    """
    if rng is None:
        rng = np.random.default_rng(42)

    # Number of throats from entropy matching
    bh = CoherentCondensate(mass=mass)
    N = bh.throat_count_from_area(states_per_throat)
    N = max(N, 1)

    # Default seed modes
    if seed_modes is None:
        seed_modes = [(1, 0, 3.14), (1, 1, 6.28), (3, 0, 5.50)]

    # Build N throats with random orientations summing to ~0 (neutral)
    orientations = rng.choice([-1, +1], size=N)
    # Ensure exact neutrality for even N
    if N > 1:
        excess = orientations.sum()
        if excess != 0:
            sign = int(np.sign(excess))
            flips = np.where(orientations == sign)[0]
            n_flip = min(abs(excess) // 2, len(flips))
            orientations[flips[:n_flip]] = -sign

    throats = []
    for i in range(N):
        modes = {}
        for l, n, omega in seed_modes:
            # Small random perturbation around coherent amplitude
            a0 = 0.01 * (1.0 + 0.05 * rng.standard_normal())
            tm = ThroatMode(
                l=l, n=n, omega=omega,
                alpha_q=1.0, gamma=GAMMA_MODE,
                a=a0, adot=0.0,
            )
            modes[(l, n)] = tm

        phase = 0.0 if coherent else float(rng.uniform(0, 2 * np.pi))
        throats.append(ThroatState(
            orientation=int(orientations[i]),
            modes=modes,
            phase=phase,
        ))

    return CoherentCondensate(mass=mass, throats=throats)


def build_charged_condensate(
    mass: float,
    charge: int,
    *,
    seed_modes: Optional[List[Tuple[int, int, float]]] = None,
    states_per_throat: int = 2,
    rng: Optional[np.random.Generator] = None,
) -> CoherentCondensate:
    """Construct a Reissner-Nordström-like charged condensate.

    Same as Schwarzschild condensate but with a net orientation bias
    producing charge Q.

    Parameters
    ----------
    mass : float
        Gravitational mass M.
    charge : int
        Net topological charge (sum of throat orientations).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    bh = CoherentCondensate(mass=mass)
    N = bh.throat_count_from_area(states_per_throat)
    N = max(N, max(abs(charge), 1))

    if seed_modes is None:
        seed_modes = [(1, 0, 3.14), (1, 1, 6.28), (3, 0, 5.50)]

    # Distribute orientations to achieve net charge Q
    n_plus = (N + charge) // 2
    n_minus = N - n_plus
    orientations = np.array([+1] * n_plus + [-1] * n_minus)
    rng.shuffle(orientations)

    throats = []
    for i in range(N):
        modes = {}
        for l, n, omega in seed_modes:
            a0 = 0.01 * (1.0 + 0.05 * rng.standard_normal())
            tm = ThroatMode(
                l=l, n=n, omega=omega,
                alpha_q=1.0, gamma=GAMMA_MODE,
                a=a0, adot=0.0,
            )
            modes[(l, n)] = tm

        throats.append(ThroatState(
            orientation=int(orientations[i]),
            modes=modes,
            phase=0.0,
        ))

    return CoherentCondensate(mass=mass, throats=throats)
