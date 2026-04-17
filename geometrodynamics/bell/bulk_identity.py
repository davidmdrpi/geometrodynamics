"""
Kinematic Bell from shared bulk identity.

This module proves the Bell correlation law from pure topology alone,
with no cavity stepping, no history evolution, and no time-domain
handshake.  It provides a second, independent path to the same E(a,b)
that the dynamic cavity/history modules produce.

The argument is:

1. Alice and Bob are two boundary mouths of ONE continuous 5D
   topological object (a non-orientable wormhole throat).

2. The throat's non-orientable transport is T = iσ_y, derived from
   the unique orientation-reversing, Hopf-preserving isometry of S³.

3. A detector at angle θ imposes a boundary constraint on its mouth:
   the spin state is projected onto the axis n(θ) = (sinθ, 0, cosθ).

4. Because both mouths belong to the same object, the boundary
   constraints at A and B jointly determine the allowed state of
   the single continuous throat.

5. The correlation is a geometric projection law:
       E(a,b) = −cos(a − b)
   following from T and the SU(2) projection, with nothing else.

This is the "higher-dimensional locality" reading: the correlation
is not nonlocal action — it is two boundary samples of one object
that extends through the bulk.

The module is deliberately minimal: one class, one correlation
function, one set of tests.  No imports from cavity, history, or
transaction.  Only embedding topology and Hopf geometry.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from geometrodynamics.embedding.transport import (
    derive_throat_transport,
    derive_singlet_from_transport,
)


@dataclass(frozen=True)
class BulkConnectedPair:
    """One continuous non-orientable throat with two boundary mouths.

    This is NOT two particles.  It is one 5D topological defect
    observed from two boundary frames on S³.

    The pair state |Ψ⟩ is fixed by the bulk topology:
        |Ψ⟩ = N × Σ_s |s⟩_A ⊗ T|s⟩_B

    where T = iσ_y is the non-orientable transport through the throat.
    No dynamics, no cavity, no history — just topology.

    Attributes
    ----------
    pair_state : ndarray
        The 4-vector |Ψ⟩ in the |↑↑⟩, |↑↓⟩, |↓↑⟩, |↓↓⟩ basis.
    transport : ndarray
        The 2×2 throat transport matrix T.
    """

    pair_state: np.ndarray
    transport: np.ndarray

    @property
    def is_antisymmetric(self) -> bool:
        """Check that the pair state is antisymmetric under exchange."""
        psi = self.pair_state
        swapped = np.array([psi[0], psi[2], psi[1], psi[3]])
        return bool(np.allclose(swapped, -psi))

    @property
    def is_normalised(self) -> bool:
        return bool(abs(np.dot(self.pair_state.conj(), self.pair_state) - 1.0) < 1e-12)

    @property
    def is_single_object(self) -> bool:
        """The pair state is entangled (not separable).

        A separable state |ψ_A⟩⊗|ψ_B⟩ has Schmidt rank 1.
        The singlet has Schmidt rank 2.  This confirms that the
        two mouths cannot be described independently — they are
        one object.
        """
        # Reshape into 2×2 matrix and check rank
        rho = self.pair_state.reshape(2, 2)
        sv = np.linalg.svd(rho, compute_uv=False)
        # Schmidt rank > 1 means entangled (non-separable)
        return int(np.sum(sv > 1e-10)) > 1


def make_bulk_pair() -> BulkConnectedPair:
    """Create a bulk-connected pair from the derived throat transport.

    The pair state is constructed from T = iσ_y via:
        |Ψ⟩ = N × Σ_s |s⟩ ⊗ T|s⟩

    No cavity.  No history.  Pure topology.
    """
    T = derive_throat_transport()
    psi = derive_singlet_from_transport()
    return BulkConnectedPair(pair_state=psi, transport=T)


# ── Boundary constraints (detector settings) ────────────────────────────────

def boundary_projector(theta: float, outcome: int) -> np.ndarray:
    """Detector setting θ as a boundary constraint on one mouth.

    The detector selects the spin component along n(θ) = (sinθ, 0, cosθ).
    This is a boundary condition on the bulk throat, not a measurement
    on an independent particle.

    |+_θ⟩ = cos(θ/2)|↑⟩ + sin(θ/2)|↓⟩
    |−_θ⟩ = −sin(θ/2)|↑⟩ + cos(θ/2)|↓⟩
    """
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    if outcome == +1:
        return np.array([c, s], dtype=complex)
    else:
        return np.array([-s, c], dtype=complex)


# ── Geometric projection law ────────────────────────────────────────────────

def bulk_amplitude(
    pair: BulkConnectedPair,
    theta_a: float,
    theta_b: float,
    outcome_a: int,
    outcome_b: int,
) -> complex:
    """Amplitude from geometric projection of the bulk state.

    Both detectors constrain the SAME continuous throat.
    The amplitude is the overlap of the joint boundary constraint
    with the bulk pair state:

        A(s_a, s_b | a, b) = ⟨s_a(a), s_b(b) | Ψ⟩

    This is purely kinematic — no time evolution.
    """
    bra_a = boundary_projector(theta_a, outcome_a).conj()
    bra_b = boundary_projector(theta_b, outcome_b).conj()
    bra = np.kron(bra_a, bra_b)
    return complex(np.dot(bra, pair.pair_state))


def bulk_probability(
    pair: BulkConnectedPair,
    theta_a: float,
    theta_b: float,
    outcome_a: int,
    outcome_b: int,
) -> float:
    """Probability from the geometric projection law.

    P(s_a, s_b | a, b) = |A(s_a, s_b | a, b)|²
    """
    return float(abs(bulk_amplitude(pair, theta_a, theta_b, outcome_a, outcome_b)) ** 2)


def bulk_correlation(
    pair: BulkConnectedPair,
    theta_a: float,
    theta_b: float,
) -> float:
    """Bell correlation from pure bulk topology.

    E(a, b) = Σ_{s_a, s_b} s_a s_b P(s_a, s_b | a, b)

    No cavity.  No history.  No time stepping.
    Just: one object, two boundary constraints, geometric projection.
    """
    E = 0.0
    for sa in [+1, -1]:
        for sb in [+1, -1]:
            p = bulk_probability(pair, theta_a, theta_b, sa, sb)
            E += sa * sb * p
    return float(E)


def bulk_marginal(
    pair: BulkConnectedPair,
    theta_a: float,
    theta_b: float,
    outcome_a: int,
) -> float:
    """Alice's marginal from the bulk projection.

    P(s_a | a, b) = Σ_{s_b} P(s_a, s_b | a, b)

    No-signaling means this is independent of theta_b.
    """
    return sum(
        bulk_probability(pair, theta_a, theta_b, outcome_a, sb)
        for sb in [+1, -1]
    )


def bulk_chsh(
    pair: BulkConnectedPair,
    a: float = 0.0,
    a_prime: float = np.pi / 2,
    b: float = np.pi / 4,
    b_prime: float = -np.pi / 4,
) -> float:
    """CHSH parameter from pure bulk topology.

    S = |E(a,b) + E(a,b') + E(a',b) − E(a',b')|
    """
    E_ab = bulk_correlation(pair, a, b)
    E_abp = bulk_correlation(pair, a, b_prime)
    E_apb = bulk_correlation(pair, a_prime, b)
    E_apbp = bulk_correlation(pair, a_prime, b_prime)
    return float(abs(E_ab + E_abp + E_apb - E_apbp))
