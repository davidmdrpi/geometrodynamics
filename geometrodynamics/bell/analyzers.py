"""
Bell analyzers — detector settings as local SU(2) boundary operators.

A detector at angle θ measures the projection of the wormhole mouth's
spinor state onto the axis n = (sin θ, 0, cos θ) in the SU(2) fiber
of the Hopf fibration.

This is not an ad hoc quantum-mechanical overlay.  The Hopf connection
on S³ naturally carries an SU(2) structure, and the throat transport
through a non-orientable wormhole is the time-reversal operator iσ_y.
The detector setting selects which Hopf fiber to project onto.

Key functions:

    throat_transport()     — the non-orientable transport matrix iσ_y
    measurement_basis()    — eigenstates of spin along detector axis θ
    singlet_amplitudes()   — outcome amplitudes from throat + detection
    branch_probabilities() — probabilities including cavity branch weights
"""

from __future__ import annotations

import numpy as np


# ── Pauli matrices ───────────────────────────────────────────────────────────

SIGMA_X = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMA_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
SIGMA_Z = np.array([[1, 0], [0, -1]], dtype=complex)
IDENTITY = np.eye(2, dtype=complex)


# ── Throat transport ─────────────────────────────────────────────────────────

def throat_transport() -> np.ndarray:
    """Non-orientable throat transport operator: T = iσ_y.

    This is the time-reversal / charge-conjugation map for spin-½,
    derived from the Hopf fibration structure of the non-orientable
    wormhole throat.

    Action:  T|↑⟩ = −|↓⟩,  T|↓⟩ = |↑⟩.

    The antisymmetric structure produces the singlet state:
        |Ψ⟩ = (1/√2)(|↑↓⟩ − |↓↑⟩)
    """
    return 1j * SIGMA_Y  # [[0, 1], [-1, 0]]


def rotation_matrix(theta: float) -> np.ndarray:
    """SU(2) rotation by angle θ around the y-axis.

    R_y(θ) = exp(−iσ_y θ/2) = [[cos(θ/2), −sin(θ/2)],
                                  [sin(θ/2),  cos(θ/2)]]

    This rotates the measurement basis from the z-axis to the axis
    n = (sin θ, 0, cos θ) in the xz-plane.  The connection to the
    Hopf fibration: θ corresponds to the hyper-latitude χ on S³,
    and the rotation is the holonomy transport along the fiber.
    """
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def measurement_projector(theta: float, outcome: int) -> np.ndarray:
    """Projector onto spin outcome ±1 along axis at angle θ.

    |+_θ⟩ = cos(θ/2)|↑⟩ + sin(θ/2)|↓⟩
    |−_θ⟩ = −sin(θ/2)|↑⟩ + cos(θ/2)|↓⟩

    Returns |s_θ⟩ as a column vector.
    """
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    if outcome == +1:
        return np.array([[c], [s]], dtype=complex)
    else:
        return np.array([[-s], [c]], dtype=complex)


# ── Pair states constructed from throat transport ────────────────────────────

def singlet_state() -> np.ndarray:
    """Singlet state constructed from throat transport T = iσ_y.

    |Ψ_s⟩ = (1/√2) Σ_s |s⟩ ⊗ T|s⟩

    where the sum is over the spin basis {|↑⟩, |↓⟩} and T is the
    non-orientable throat transport operator.  The construction is:

        |↑⟩ ⊗ T|↑⟩ + |↓⟩ ⊗ T|↓⟩  →  normalise

    This produces the antisymmetric (singlet) state from the throat
    topology, not from a postulate.

    Returns a 4-vector in the |↑↑⟩, |↑↓⟩, |↓↑⟩, |↓↓⟩ basis.
    """
    T = throat_transport()
    up = np.array([1, 0], dtype=complex)
    dn = np.array([0, 1], dtype=complex)

    # |↑⟩ ⊗ T|↑⟩ + |↓⟩ ⊗ T|↓⟩
    psi = np.kron(up, T @ up) + np.kron(dn, T @ dn)
    norm = np.linalg.norm(psi)
    return psi / norm


def triplet_m0_state() -> np.ndarray:
    """Triplet m=0 state: (1/√2)(|↑↓⟩ + |↓↑⟩).

    This is the π-branch partner of the singlet.  Where the 0-branch
    (singlet) has the throat transport producing antisymmetric exchange,
    the π-branch phase shift produces symmetric exchange.

    The triplet gives E(a,b) = +cos(a−b), opposite to the singlet.
    """
    psi = np.array([0, 1, 1, 0], dtype=complex) / np.sqrt(2)
    return psi


def singlet_amplitude(
    theta_a: float,
    theta_b: float,
    outcome_a: int,
    outcome_b: int,
) -> complex:
    """Amplitude for outcomes (s_a, s_b) given detector settings (θ_a, θ_b).

    Computed from the geometric chain:
        1. Construct the pair state from throat transport T = iσ_y
        2. Project mouth A onto detector axis θ_a with outcome s_a
        3. Project mouth B onto detector axis θ_b with outcome s_b
        4. The amplitude is ⟨s_a(θ_a), s_b(θ_b) | Ψ⟩
    """
    state_a = measurement_projector(theta_a, outcome_a)
    state_b = measurement_projector(theta_b, outcome_b)
    bra = np.kron(state_a.conj().T, state_b.conj().T).flatten()
    psi = singlet_state()
    return complex(np.dot(bra, psi))


def triplet_amplitude(
    theta_a: float,
    theta_b: float,
    outcome_a: int,
    outcome_b: int,
) -> complex:
    """Amplitude for the triplet (π-branch) state."""
    state_a = measurement_projector(theta_a, outcome_a)
    state_b = measurement_projector(theta_b, outcome_b)
    bra = np.kron(state_a.conj().T, state_b.conj().T).flatten()
    psi = triplet_m0_state()
    return complex(np.dot(bra, psi))


def outcome_probability(
    theta_a: float,
    theta_b: float,
    outcome_a: int,
    outcome_b: int,
) -> float:
    """Probability P(s_a, s_b | θ_a, θ_b) for the pure singlet."""
    amp = singlet_amplitude(theta_a, theta_b, outcome_a, outcome_b)
    return float(abs(amp) ** 2)


def branch_weighted_probability(
    theta_a: float,
    theta_b: float,
    outcome_a: int,
    outcome_b: int,
    w0: float = 1.0,
    wpi: float = 0.0,
) -> float:
    """Branch-weighted probability: P = w₀|A_singlet|² + w_π|A_triplet|².

    The 0-branch (singlet) and π-branch (triplet) are the two closure
    channels of the cavity.  This computes the full mixture probability
    for each outcome pair.
    """
    p_s = abs(singlet_amplitude(theta_a, theta_b, outcome_a, outcome_b)) ** 2
    p_t = abs(triplet_amplitude(theta_a, theta_b, outcome_a, outcome_b)) ** 2
    return float(w0 * p_s + wpi * p_t)


def singlet_correlation(theta_a: float, theta_b: float) -> float:
    """Singlet correlation E(θ_a, θ_b) = −cos(θ_a − θ_b).

    Computed explicitly from outcome probabilities — not from the
    analytic formula.
    """
    pp = outcome_probability(theta_a, theta_b, +1, +1)
    pm = outcome_probability(theta_a, theta_b, +1, -1)
    mp = outcome_probability(theta_a, theta_b, -1, +1)
    mm = outcome_probability(theta_a, theta_b, -1, -1)
    return float(pp + mm - pm - mp)
