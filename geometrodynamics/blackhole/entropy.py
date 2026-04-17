"""
Black-hole entropy and thermodynamics from throat microstate counting.

The central result: if a black hole of mass M is a condensate of N
non-orientable wormhole throats, and each throat carries k independent
microstates, then

    S_throat = N · ln(k)

matches the Bekenstein-Hawking entropy S_BH = A / (4 l_P²) when

    N = A / a_min,    a_min = 4 l_P² / ln(k)

For k = 2 (orientation ±1), a_min ≈ 5.77 l_P² per throat.

This module validates the entropy matching, derives the first law
dM = T dS + Φ dQ, and computes the information capacity of the
throat network.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from geometrodynamics.blackhole.condensate import (
    CoherentCondensate,
    L_PLANCK_SQ,
    A_THROAT_MIN,
)
from geometrodynamics.blackhole.interior import (
    f_hayward,
    find_horizons,
    surface_gravity,
    hawking_temperature,
)


# ── Entropy matching ─────────────────────────────────────────────────────────

@dataclass
class EntropyBalance:
    """Comparison of Bekenstein-Hawking and throat-counting entropy."""

    mass: float
    core_scale: float
    r_horizon: float
    area: float
    S_BH: float               # Bekenstein-Hawking: A/(4 l_P²)
    N_throats: int             # actual throat count in condensate
    states_per_throat: int
    S_throat: float            # N · ln(k)
    relative_error: float      # |S_BH - S_throat| / S_BH
    temperature_geometric: float   # from surface gravity κ/(2π)
    temperature_standard: float    # 1/(8πM) (Schwarzschild limit)


def compute_entropy_balance(
    condensate: CoherentCondensate,
    states_per_throat: int = 2,
) -> EntropyBalance:
    """Compute and compare BH and throat entropy for a condensate.

    Parameters
    ----------
    condensate : CoherentCondensate
        The black-hole microscopic model.
    states_per_throat : int
        Number of microstates per throat (default 2).

    Returns
    -------
    EntropyBalance
        Detailed comparison of the two entropy computations.
    """
    M = condensate.mass
    l = condensate.core_scale
    r_h = condensate.r_horizon
    A = condensate.horizon_area

    S_BH = A / (4.0 * L_PLANCK_SQ)
    S_throat = condensate.entropy_from_throats(states_per_throat)

    rel_err = abs(S_BH - S_throat) / max(S_BH, 1e-30)

    T_geom = hawking_temperature(M, l)
    T_std = 1.0 / (8.0 * np.pi * M) if M > 0 else float("inf")

    return EntropyBalance(
        mass=M,
        core_scale=l,
        r_horizon=r_h,
        area=A,
        S_BH=S_BH,
        N_throats=condensate.N,
        states_per_throat=states_per_throat,
        S_throat=S_throat,
        relative_error=rel_err,
        temperature_geometric=T_geom,
        temperature_standard=T_std,
    )


# ── First law of black-hole thermodynamics ───────────────────────────────────

@dataclass
class FirstLawCheck:
    """Numerical verification of dM = T dS + Φ dQ."""

    dM: float           # mass perturbation
    T: float            # temperature
    dS: float           # entropy change
    TdS: float          # T · dS
    Phi: float          # electric potential at horizon
    dQ: float           # charge change
    PhidQ: float        # Φ · dQ
    residual: float     # |dM - TdS - ΦdQ|
    relative_residual: float


def check_first_law(
    condensate: CoherentCondensate,
    dM: float = 0.01,
    dQ: int = 0,
    states_per_throat: int = 2,
) -> FirstLawCheck:
    """Numerically verify the first law dM = T dS + Φ dQ.

    Perturbs the mass by dM and (optionally) the charge by dQ,
    recomputes the entropy, and checks the thermodynamic identity.

    Parameters
    ----------
    condensate : CoherentCondensate
        Unperturbed black hole.
    dM : float
        Mass perturbation.
    dQ : int
        Charge perturbation (number of throat orientation flips).
    """
    M = condensate.mass
    l = condensate.core_scale
    Q = condensate.net_charge
    N = condensate.N

    # Temperature from surface gravity
    T = hawking_temperature(M, l)

    # Entropy before and after
    S1 = condensate.entropy_from_area()

    M2 = M + dM
    A2 = 16.0 * np.pi * M2 ** 2
    S2 = A2 / (4.0 * L_PLANCK_SQ)

    dS = S2 - S1
    TdS = T * dS

    # Electric potential at the horizon: Φ = Q / r_H for RN
    r_h = 2.0 * M  # approximate
    Phi = Q / max(r_h, 1e-10) if Q != 0 else 0.0
    PhidQ_val = Phi * dQ

    residual = abs(dM - TdS - PhidQ_val)
    rel_res = residual / max(abs(dM), 1e-30)

    return FirstLawCheck(
        dM=dM,
        T=T,
        dS=dS,
        TdS=TdS,
        Phi=Phi,
        dQ=dQ,
        PhidQ=PhidQ_val,
        residual=residual,
        relative_residual=rel_res,
    )


# ── Information capacity ─────────────────────────────────────────────────────

def information_capacity_bits(condensate: CoherentCondensate) -> float:
    """Total information capacity of the throat network in bits.

    I = S / ln(2) = N · log₂(k)

    This is the Bekenstein bound: the maximum information that can be
    stored in a region bounded by area A.
    """
    S = condensate.entropy_from_area()
    return S / np.log(2)


def page_time(condensate: CoherentCondensate) -> float:
    """Page time: when half the information has been emitted.

    t_Page ≈ (5120 π / σ) M³  where σ depends on species count.
    In geometric units with only graviton emission: t_Page ≈ 5120π M³.

    This is the timescale on which the throat condensate transitions
    from a high-coherence to a low-coherence state, analogous to
    the BEC → thermal gas transition.
    """
    M = condensate.mass
    return 5120.0 * np.pi * M ** 3


# ── Evaporation as condensate decoherence ────────────────────────────────────

def evaporation_step(
    condensate: CoherentCondensate,
    *,
    dm: Optional[float] = None,
    rng: Optional[np.random.Generator] = None,
) -> CoherentCondensate:
    """One step of Hawking evaporation as throat loss from the condensate.

    Instead of particle-pair creation at the horizon, evaporation is
    modelled as the decoupling of one throat from the coherent condensate.
    The throat's energy (∝ T_H) is lost, reducing the mass.

    Parameters
    ----------
    condensate : CoherentCondensate
        Current black-hole state.
    dm : float or None
        Mass to lose in this step.  Default: T_H (one thermal quantum).

    Returns
    -------
    CoherentCondensate
        Updated condensate with one fewer throat and reduced mass.
    """
    if rng is None:
        rng = np.random.default_rng()

    if condensate.N <= 1:
        return condensate

    if dm is None:
        dm = condensate.hawking_temperature

    # Remove a random throat
    idx = rng.integers(0, condensate.N)
    new_throats = [t for i, t in enumerate(condensate.throats) if i != idx]

    new_mass = max(condensate.mass - dm, 1e-10)

    return CoherentCondensate(
        mass=new_mass,
        throats=new_throats,
        l_core=condensate.l_core,
    )
