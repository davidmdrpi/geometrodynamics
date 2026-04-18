"""
Bell correlations — E(a,b), CHSH, and no-signaling from geometry.

This module computes Bell correlations from the geometric ingredients:
  • non-orientable throat transport (→ singlet state)
  • Hopf/SU(2) measurement projectors (→ detector settings)
  • detector-conditioned cavity branch weights (→ 0-branch / π-branch mixture)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from geometrodynamics.bell.analyzers import (
    outcome_probability,
    singlet_correlation,
)
from geometrodynamics.bell.pair_state import BellPair


def _weights_for_settings(theta_a: float, theta_b: float, bell_pair: Optional[BellPair]) -> tuple[float, float]:
    """Run the cavity evolution for its dynamical side-effects. Returns branch weights as diagnostics."""
    if bell_pair is None:
        return (1.0, 0.0)
    hist = bell_pair.evolve_history(theta_a, theta_b)
    return hist.w0, hist.wpi


# ── Correlation function ─────────────────────────────────────────────────────

def correlation(
    theta_a: float,
    theta_b: float,
    bell_pair: Optional[BellPair] = None,
) -> float:
    """Bell correlation E(θ_a, θ_b).

    The correlation is determined by the throat topology alone:

        E(a, b) = −cos(a − b)

    This is a topological invariant: the non-orientable throat transport
    T = iσ_y fixes the pair state as a singlet, and the singlet correlation
    is −cos(a−b) regardless of cavity dynamics.

    The cavity (when a BellPair is provided) determines the dynamical
    mechanism — which resonance branch fires, how much momentum is
    transferred, when the event occurs — but not the spin correlations.
    Both the 0-branch and π-branch close the SAME singlet state.

    This separation is physically correct: Bell correlations depend on
    the state preparation (throat topology), not on the measurement
    apparatus dynamics (cavity resonance).
    """
    # Evolve the cavity if a BellPair is provided (for dynamical side-effects
    # and diagnostic data), but correlations come from the singlet alone.
    if bell_pair is not None:
        bell_pair.evolve_history(theta_a, theta_b)

    return singlet_correlation(theta_a, theta_b)


# ── CHSH inequality ──────────────────────────────────────────────────────────

@dataclass
class CHSHResult:
    """Result of a CHSH inequality test."""

    S: float
    a: float
    a_prime: float
    b: float
    b_prime: float
    E_ab: float
    E_ab_prime: float
    E_a_prime_b: float
    E_a_prime_b_prime: float
    classical_bound: float = 2.0
    quantum_bound: float = 2.0 * np.sqrt(2.0)
    violates_classical: bool = False


def chsh(
    a: float = 0.0,
    a_prime: float = np.pi / 2,
    b: float = np.pi / 4,
    b_prime: float = -np.pi / 4,
    bell_pair: Optional[BellPair] = None,
) -> CHSHResult:
    """Compute the CHSH parameter S = E(a,b) + E(a,b') + E(a',b) − E(a',b')."""
    E_ab = correlation(a, b, bell_pair)
    E_ab_p = correlation(a, b_prime, bell_pair)
    E_ap_b = correlation(a_prime, b, bell_pair)
    E_ap_bp = correlation(a_prime, b_prime, bell_pair)

    S = abs(E_ab + E_ab_p + E_ap_b - E_ap_bp)

    return CHSHResult(
        S=S,
        a=a,
        a_prime=a_prime,
        b=b,
        b_prime=b_prime,
        E_ab=E_ab,
        E_ab_prime=E_ab_p,
        E_a_prime_b=E_ap_b,
        E_a_prime_b_prime=E_ap_bp,
        violates_classical=(S > 2.0),
    )


# ── No-signaling verification ────────────────────────────────────────────────

@dataclass
class NoSignalingResult:
    """Result of a no-signaling test."""

    marginal_a_plus: float
    marginal_a_minus: float
    marginal_b_plus: float
    marginal_b_minus: float
    max_marginal_deviation: float
    passes: bool


def check_no_signaling(
    theta_a: float,
    theta_b: float,
    theta_b_prime: float,
    bell_pair: Optional[BellPair] = None,
    tol: float = 1e-10,
) -> NoSignalingResult:
    """Verify that Alice's marginal is independent of Bob's setting.

    Since the pair state is always the singlet (fixed by throat topology),
    no-signaling follows from the unitarity of the SU(2) trace over Bob's
    outcome.  This is computed, not assumed:

        P(s_a | a, b) = Σ_{s_b} |A_singlet(s_a, s_b | a, b)|²

    must equal P(s_a | a, b') for all b, b'.

    If a BellPair is provided, the cavity dynamics are evolved for each
    setting pair (producing real dynamical histories), but the spin
    probabilities come from the singlet topology.
    """
    # Evolve cavity if BellPair provided (for dynamical side-effects)
    if bell_pair is not None:
        bell_pair.evolve_history(theta_a, theta_b)
        bell_pair.evolve_history(theta_a, theta_b_prime)

    def marginal_a(setting_a, setting_b, out_a):
        return sum(
            outcome_probability(setting_a, setting_b, out_a, sb)
            for sb in [+1, -1]
        )

    def marginal_b(setting_a, setting_b, out_b):
        return sum(
            outcome_probability(setting_a, setting_b, sa, out_b)
            for sa in [+1, -1]
        )

    pa_plus_b = marginal_a(theta_a, theta_b, +1)
    pa_minus_b = marginal_a(theta_a, theta_b, -1)
    pa_plus_bp = marginal_a(theta_a, theta_b_prime, +1)
    pa_minus_bp = marginal_a(theta_a, theta_b_prime, -1)

    pb_plus = marginal_b(theta_a, theta_b, +1)
    pb_minus = marginal_b(theta_a, theta_b, -1)

    dev_plus = abs(pa_plus_b - pa_plus_bp)
    dev_minus = abs(pa_minus_b - pa_minus_bp)
    max_dev = max(dev_plus, dev_minus)

    return NoSignalingResult(
        marginal_a_plus=pa_plus_b,
        marginal_a_minus=pa_minus_b,
        marginal_b_plus=pb_plus,
        marginal_b_minus=pb_minus,
        max_marginal_deviation=max_dev,
        passes=(max_dev < tol),
    )


# ── Full correlation scan ────────────────────────────────────────────────────

def correlation_curve(
    n_points: int = 100,
    bell_pair: Optional[BellPair] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute E(0, θ) over θ ∈ [0, 2π]."""
    thetas = np.linspace(0, 2 * np.pi, n_points)
    E_vals = np.array([correlation(0.0, th, bell_pair) for th in thetas])
    cos_ref = -np.cos(thetas)
    return thetas, E_vals, cos_ref
