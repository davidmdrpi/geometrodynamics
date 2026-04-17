"""
Tests for the kinematic shared-bulk Bell interpretation.

Proves Bell correlations from pure topology alone — no cavity,
no history evolution, no time stepping.  Every test uses only
BulkConnectedPair + geometric projection.

What this test class proves:

  1. The pair is ONE object (Schmidt rank > 1, non-separable)
  2. Detector settings act as boundary constraints on the same throat
  3. The non-orientable transport T = iσ_y enforces inversion
  4. E(a,b) = −cos(a−b) from geometric projection alone
  5. Equal settings → perfect anticorrelation
  6. CHSH = 2√2
  7. No-signaling: each marginal = ½, independent of remote setting
  8. Probabilities sum to 1 for all settings
  9. Detector-label symmetry: E(a,b) = E(b,a)
 10. The bulk result exactly matches the dynamic cavity/history result

These are the SAME correlations that the dynamic Bell/cavity/history
modules produce, but derived here without any time-domain simulation.
"""

import numpy as np
import pytest

from geometrodynamics.bell.bulk_identity import (
    BulkConnectedPair,
    make_bulk_pair,
    boundary_projector,
    bulk_amplitude,
    bulk_probability,
    bulk_correlation,
    bulk_marginal,
    bulk_chsh,
)


@pytest.fixture
def pair():
    return make_bulk_pair()


# ════════════════════════════════════════════════════════════════════════════
# 1.  THE PAIR IS ONE OBJECT
# ════════════════════════════════════════════════════════════════════════════

class TestSingleObject:
    """The pair is one continuous bulk defect, not two particles."""

    def test_pair_is_normalised(self, pair):
        assert pair.is_normalised

    def test_pair_is_antisymmetric(self, pair):
        assert pair.is_antisymmetric

    def test_pair_is_non_separable(self, pair):
        """Schmidt rank > 1 confirms the two mouths are one object."""
        assert pair.is_single_object

    def test_transport_squares_to_minus_identity(self, pair):
        """T² = −I: the double cover / 4π periodicity of spin-½."""
        assert np.allclose(pair.transport @ pair.transport, -np.eye(2))


# ════════════════════════════════════════════════════════════════════════════
# 2.  BOUNDARY CONSTRAINTS
# ════════════════════════════════════════════════════════════════════════════

class TestBoundaryConstraints:
    """Detector settings act as boundary conditions on the throat."""

    def test_projectors_are_orthonormal(self):
        """At any angle, |+⟩ and |−⟩ are orthonormal."""
        for theta in [0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, np.pi, 2.5]:
            plus = boundary_projector(theta, +1)
            minus = boundary_projector(theta, -1)
            assert abs(np.dot(plus.conj(), plus) - 1.0) < 1e-12
            assert abs(np.dot(minus.conj(), minus) - 1.0) < 1e-12
            assert abs(np.dot(plus.conj(), minus)) < 1e-12

    def test_projectors_span_space(self):
        """|+⟩⟨+| + |−⟩⟨−| = I at any angle."""
        for theta in [0, np.pi/4, np.pi/2, 1.0, np.pi]:
            plus = boundary_projector(theta, +1)
            minus = boundary_projector(theta, -1)
            identity = np.outer(plus, plus.conj()) + np.outer(minus, minus.conj())
            assert np.allclose(identity, np.eye(2))

    def test_z_axis_recovers_standard_basis(self):
        """θ = 0 gives |↑⟩ and |↓⟩."""
        plus = boundary_projector(0, +1)
        minus = boundary_projector(0, -1)
        assert np.allclose(plus, [1, 0])
        assert np.allclose(minus, [0, 1])


# ════════════════════════════════════════════════════════════════════════════
# 3.  PROBABILITIES
# ════════════════════════════════════════════════════════════════════════════

class TestProbabilities:
    """Outcome probabilities from the bulk projection law."""

    def test_probabilities_sum_to_one(self, pair):
        """Σ P(s_a, s_b | a, b) = 1 for all settings."""
        rng = np.random.default_rng(42)
        for _ in range(20):
            a, b = rng.uniform(0, 2 * np.pi, 2)
            total = sum(
                bulk_probability(pair, a, b, sa, sb)
                for sa in [+1, -1] for sb in [+1, -1]
            )
            assert abs(total - 1.0) < 1e-10

    def test_probabilities_non_negative(self, pair):
        """P ≥ 0 for all outcomes and settings."""
        rng = np.random.default_rng(42)
        for _ in range(20):
            a, b = rng.uniform(0, 2 * np.pi, 2)
            for sa in [+1, -1]:
                for sb in [+1, -1]:
                    p = bulk_probability(pair, a, b, sa, sb)
                    assert p >= -1e-15


# ════════════════════════════════════════════════════════════════════════════
# 4.  PERFECT ANTICORRELATION
# ════════════════════════════════════════════════════════════════════════════

class TestAnticorrelation:
    """Equal settings → perfect anticorrelation."""

    def test_equal_settings(self, pair):
        """E(θ, θ) = −1 for any θ."""
        for theta in [0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, np.pi, 2.0]:
            E = bulk_correlation(pair, theta, theta)
            assert abs(E - (-1.0)) < 1e-10, (
                f"E({theta:.3f}, {theta:.3f}) = {E:.6f}"
            )

    def test_opposite_settings(self, pair):
        """E(θ, θ+π) = +1."""
        for theta in [0, np.pi/4, np.pi/2, 1.0]:
            E = bulk_correlation(pair, theta, theta + np.pi)
            assert abs(E - 1.0) < 1e-10


# ════════════════════════════════════════════════════════════════════════════
# 5.  CORRELATION = −cos(a − b)
# ════════════════════════════════════════════════════════════════════════════

class TestCorrelationLaw:
    """The Bell curve from pure geometric projection."""

    def test_minus_cosine(self, pair):
        """E(a, b) = −cos(a − b) for random settings."""
        rng = np.random.default_rng(42)
        for _ in range(50):
            a, b = rng.uniform(0, 2 * np.pi, 2)
            E = bulk_correlation(pair, a, b)
            expected = -np.cos(a - b)
            assert abs(E - expected) < 1e-10, (
                f"E({a:.3f}, {b:.3f}) = {E:.6f}, expected {expected:.6f}"
            )

    def test_detector_symmetry(self, pair):
        """E(a, b) = E(b, a) — swapping boundary labels preserves the law."""
        rng = np.random.default_rng(42)
        for _ in range(20):
            a, b = rng.uniform(0, 2 * np.pi, 2)
            assert abs(bulk_correlation(pair, a, b) - bulk_correlation(pair, b, a)) < 1e-10


# ════════════════════════════════════════════════════════════════════════════
# 6.  CHSH VIOLATION
# ════════════════════════════════════════════════════════════════════════════

class TestCHSH:
    """CHSH from pure bulk topology."""

    def test_tsirelson_optimal(self, pair):
        """S = 2√2 at Tsirelson-optimal settings."""
        S = bulk_chsh(pair)
        assert abs(S - 2 * np.sqrt(2)) < 1e-10

    def test_exceeds_classical_bound(self, pair):
        """S > 2 violates the classical bound."""
        S = bulk_chsh(pair)
        assert S > 2.0

    def test_never_exceeds_tsirelson(self, pair):
        """S ≤ 2√2 at all settings."""
        rng = np.random.default_rng(42)
        for _ in range(50):
            a, ap, b, bp = rng.uniform(0, 2 * np.pi, 4)
            S = bulk_chsh(pair, a, ap, b, bp)
            assert S <= 2 * np.sqrt(2) + 1e-10


# ════════════════════════════════════════════════════════════════════════════
# 7.  NO-SIGNALING
# ════════════════════════════════════════════════════════════════════════════

class TestNoSignaling:
    """Each marginal = ½, independent of the remote setting."""

    def test_marginals_are_half(self, pair):
        """P(s_a | a, b) = ½ for all a, b, s_a."""
        rng = np.random.default_rng(42)
        for _ in range(20):
            a, b = rng.uniform(0, 2 * np.pi, 2)
            for sa in [+1, -1]:
                m = bulk_marginal(pair, a, b, sa)
                assert abs(m - 0.5) < 1e-10, (
                    f"P({sa:+d} | {a:.3f}, {b:.3f}) = {m:.6f}"
                )

    def test_alice_independent_of_bob(self, pair):
        """P(s_a | a, b) = P(s_a | a, b') for all b, b'."""
        for a in [0, np.pi/4, np.pi/3, 1.0]:
            for b, bp in [(0, np.pi/4), (np.pi/6, np.pi/2), (1.0, 2.5)]:
                for sa in [+1, -1]:
                    m1 = bulk_marginal(pair, a, b, sa)
                    m2 = bulk_marginal(pair, a, bp, sa)
                    assert abs(m1 - m2) < 1e-10


# ════════════════════════════════════════════════════════════════════════════
# 8.  PATH EQUIVALENCE
# ════════════════════════════════════════════════════════════════════════════

class TestPathEquivalence:
    """The bulk result matches the dynamic Bell/history result exactly."""

    def test_matches_bell_analyzers(self, pair):
        """Bulk E(a,b) = Bell analyzer E(a,b) for random settings."""
        from geometrodynamics.bell.analyzers import singlet_correlation
        rng = np.random.default_rng(42)
        for _ in range(20):
            a, b = rng.uniform(0, 2 * np.pi, 2)
            E_bulk = bulk_correlation(pair, a, b)
            E_bell = singlet_correlation(a, b)
            assert abs(E_bulk - E_bell) < 1e-10

    def test_matches_history_branches(self, pair):
        """Bulk E(a,b) = history branch E(a,b) for several settings."""
        from geometrodynamics.history.closure import enumerate_bell_branches
        for a, b in [(0, 0), (0, np.pi/4), (0, np.pi/2), (np.pi/3, np.pi/6)]:
            E_bulk = bulk_correlation(pair, a, b)
            E_hist = sum(
                x.outcome_a * x.outcome_b * x.probability
                for x in enumerate_bell_branches(a, b)
            )
            assert abs(E_bulk - E_hist) < 1e-10, (
                f"({a:.3f}, {b:.3f}): bulk={E_bulk:.6f} hist={E_hist:.6f}"
            )

    def test_chsh_matches_bell_module(self, pair):
        """Bulk CHSH = Bell module CHSH = 2√2."""
        from geometrodynamics.bell.correlations import chsh
        S_bulk = bulk_chsh(pair)
        S_bell = chsh().S
        assert abs(S_bulk - S_bell) < 1e-10
