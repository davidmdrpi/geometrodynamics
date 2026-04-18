"""
Tests for the Bell subpackage and supporting modules.

Make-or-break assertions:
  1. Equal settings → perfect anticorrelation
  2. CHSH exceeds 2 (classical bound violated)
  3. Local marginals independent of remote setting (no-signaling)
  4. Swapping detector order does not change statistics
  5. Correlation curve matches −cos(a−b)
  6. Throat transport has correct sign structure
  7. Cavity persists energy between steps
  8. Embedding topology: conjugate sectors carry opposite orientation
"""

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geometrodynamics.bell.analyzers import (
    throat_transport,
    singlet_state,
    singlet_amplitude,
    outcome_probability,
    singlet_correlation,
    measurement_projector,
)
from geometrodynamics.bell.pair_state import BellPair, make_bell_pair
from geometrodynamics.bell.correlations import (
    correlation,
    chsh,
    check_no_signaling,
    correlation_curve,
)
from geometrodynamics.embedding.topology import (
    ThroatDefect,
    ConjugatePair,
    make_singlet_pair,
)
from geometrodynamics.transaction.cavity import (
    CavityMode,
    AntipodalCavity,
    make_cavity,
)


# ════════════════════════════════════════════════════════════════════════════
# 1.  THROAT TRANSPORT
# ════════════════════════════════════════════════════════════════════════════

class TestThroatTransport:
    """Non-orientable throat transport = iσ_y."""

    def test_transport_matrix(self):
        """T = iσ_y: T|↑⟩ = −|↓⟩, T|↓⟩ = |↑⟩ (non-orientable sign flip)."""
        T = throat_transport()
        up = np.array([1, 0], dtype=complex)
        down = np.array([0, 1], dtype=complex)
        assert np.allclose(T @ up, -down), f"T|↑⟩ = {T @ up}, expected −|↓⟩"
        assert np.allclose(T @ down, up), f"T|↓⟩ = {T @ down}, expected |↑⟩"

    def test_transport_antiunitary(self):
        """T² = −I (double cover property of spin-½)."""
        T = throat_transport()
        assert np.allclose(T @ T, -np.eye(2)), "T² should be −I"


# ════════════════════════════════════════════════════════════════════════════
# 2.  SINGLET STATE
# ════════════════════════════════════════════════════════════════════════════

class TestSingletState:
    """Singlet state constructed from throat transport."""

    def test_singlet_normalisation(self):
        """⟨Ψ|Ψ⟩ = 1."""
        psi = singlet_state()
        assert abs(np.dot(psi.conj(), psi) - 1.0) < 1e-12

    def test_singlet_antisymmetric(self):
        """Singlet is antisymmetric under particle exchange."""
        psi = singlet_state()
        psi_swapped = np.array([psi[0], psi[2], psi[1], psi[3]])
        assert np.allclose(psi_swapped, -psi), "Singlet should be antisymmetric"

    def test_singlet_built_from_transport(self):
        """Singlet state is Σ |s⟩⊗T|s⟩ / norm, NOT hardcoded."""
        T = throat_transport()
        up = np.array([1, 0], dtype=complex)
        dn = np.array([0, 1], dtype=complex)
        psi_manual = np.kron(up, T @ up) + np.kron(dn, T @ dn)
        psi_manual = psi_manual / np.linalg.norm(psi_manual)
        psi_func = singlet_state()
        # Equal up to global phase
        overlap = abs(np.dot(psi_manual.conj(), psi_func))
        assert abs(overlap - 1.0) < 1e-12, (
            f"singlet_state() not built from T: overlap = {overlap}"
        )

    def test_probabilities_sum_to_one(self):
        """Σ P(s_a, s_b) = 1 for any settings."""
        for a, b in [(0, 0), (0, np.pi/4), (np.pi/3, np.pi/6), (1.0, 2.0)]:
            total = sum(
                outcome_probability(a, b, sa, sb)
                for sa in [+1, -1] for sb in [+1, -1]
            )
            assert abs(total - 1.0) < 1e-10, (
                f"Probabilities sum to {total} at a={a:.2f}, b={b:.2f}"
            )


# ════════════════════════════════════════════════════════════════════════════
# 3.  PERFECT ANTICORRELATION
# ════════════════════════════════════════════════════════════════════════════

class TestAnticorrelation:
    """Equal settings → perfect anticorrelation."""

    def test_equal_settings_anticorrelation(self):
        """E(θ, θ) = −1 for any θ."""
        for theta in [0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, np.pi]:
            E = singlet_correlation(theta, theta)
            assert abs(E - (-1.0)) < 1e-10, (
                f"E({theta:.3f}, {theta:.3f}) = {E:.6f}, expected −1"
            )

    def test_opposite_settings_correlation(self):
        """E(θ, θ+π) = +1."""
        for theta in [0, np.pi/4, np.pi/2, 1.0]:
            E = singlet_correlation(theta, theta + np.pi)
            assert abs(E - 1.0) < 1e-10, (
                f"E({theta:.3f}, {theta+np.pi:.3f}) = {E:.6f}, expected +1"
            )

    def test_correlation_is_minus_cosine(self):
        """E(a, b) = −cos(a − b)."""
        rng = np.random.default_rng(42)
        for _ in range(20):
            a = rng.uniform(0, 2 * np.pi)
            b = rng.uniform(0, 2 * np.pi)
            E = singlet_correlation(a, b)
            expected = -np.cos(a - b)
            assert abs(E - expected) < 1e-10, (
                f"E({a:.3f}, {b:.3f}) = {E:.6f}, expected {expected:.6f}"
            )


# ════════════════════════════════════════════════════════════════════════════
# 4.  CHSH VIOLATION
# ════════════════════════════════════════════════════════════════════════════

class TestCHSH:
    """CHSH S exceeds the classical bound of 2."""

    def test_chsh_tsirelson(self):
        """Tsirelson-optimal settings give S = 2√2."""
        result = chsh()
        expected = 2.0 * np.sqrt(2.0)
        assert abs(result.S - expected) < 1e-10, (
            f"S = {result.S:.6f}, expected {expected:.6f}"
        )
        assert result.violates_classical

    def test_chsh_exceeds_two(self):
        """S > 2 at the standard settings (violates classical bound)."""
        result = chsh()
        assert result.S > 2.0, f"S = {result.S:.6f} ≤ 2"

    def test_chsh_does_not_exceed_tsirelson(self):
        """S ≤ 2√2 at all settings (Tsirelson bound)."""
        rng = np.random.default_rng(42)
        for _ in range(50):
            a = rng.uniform(0, 2 * np.pi)
            ap = rng.uniform(0, 2 * np.pi)
            b = rng.uniform(0, 2 * np.pi)
            bp = rng.uniform(0, 2 * np.pi)
            result = chsh(a, ap, b, bp)
            assert result.S <= 2 * np.sqrt(2) + 1e-10, (
                f"S = {result.S:.6f} exceeds Tsirelson bound"
            )

    def test_chsh_with_bell_pair(self):
        """CHSH with BellPair gives S = 2√2 (topology determines correlations)."""
        bp = make_bell_pair()
        result = chsh(bell_pair=bp)
        expected = 2.0 * np.sqrt(2.0)
        assert abs(result.S - expected) < 1e-10, (
            f"S = {result.S:.6f}, expected {expected:.6f}"
        )
        assert result.violates_classical


# ════════════════════════════════════════════════════════════════════════════
# 5.  NO-SIGNALING
# ════════════════════════════════════════════════════════════════════════════

class TestNoSignaling:
    """Local marginals independent of remote setting."""

    def test_alice_marginal_independent_of_bob(self):
        """P(s_a | a, b) = P(s_a | a, b') for all b, b'."""
        for a in [0, np.pi/4, np.pi/3, np.pi/2, 1.0]:
            for b, bp in [(0, np.pi/4), (np.pi/6, np.pi/2), (1.0, 2.5)]:
                result = check_no_signaling(a, b, bp)
                assert result.passes, (
                    f"No-signaling violated: a={a:.2f}, b={b:.2f}, b'={bp:.2f}, "
                    f"deviation={result.max_marginal_deviation:.2e}"
                )

    def test_marginals_are_half(self):
        """Each marginal = ½ (maximum ignorance without remote info)."""
        for a in [0, np.pi/3, np.pi/2, 1.0]:
            result = check_no_signaling(a, 0, np.pi)
            assert abs(result.marginal_a_plus - 0.5) < 1e-10, (
                f"P(+|a={a:.2f}) = {result.marginal_a_plus}"
            )
            assert abs(result.marginal_b_plus - 0.5) < 1e-10

    def test_no_signaling_with_bell_pair(self):
        """No-signaling holds for branch-weighted BellPair (computed, not hardcoded)."""
        bp = make_bell_pair()
        for a in [0, np.pi/4, np.pi/3, 1.0]:
            for b, bprime in [(0, np.pi/4), (np.pi/6, np.pi/2), (1.0, 2.5)]:
                result = check_no_signaling(a, b, bprime, bell_pair=bp)
                assert result.passes, (
                    f"Branch-weighted no-signaling violated: a={a:.2f}, "
                    f"b={b:.2f}, b'={bprime:.2f}, "
                    f"deviation={result.max_marginal_deviation:.2e}"
                )
                assert abs(result.marginal_a_plus - 0.5) < 1e-10, (
                    f"Branch-weighted marginal P(+|a={a:.2f}) = "
                    f"{result.marginal_a_plus}"
                )

    def test_swapping_detectors(self):
        """E(a, b) = E(b, a) — swapping detector labels preserves statistics."""
        rng = np.random.default_rng(42)
        for _ in range(20):
            a = rng.uniform(0, 2 * np.pi)
            b = rng.uniform(0, 2 * np.pi)
            E_ab = singlet_correlation(a, b)
            E_ba = singlet_correlation(b, a)
            assert abs(E_ab - E_ba) < 1e-10, (
                f"E(a,b)={E_ab:.6f} ≠ E(b,a)={E_ba:.6f}"
            )


# ════════════════════════════════════════════════════════════════════════════
# 6.  EMBEDDING TOPOLOGY
# ════════════════════════════════════════════════════════════════════════════

class TestEmbeddingTopology:
    """Non-orientable throat topology."""

    def test_conjugate_opposite_orientation(self):
        """Conjugate pair mouths have opposite orientation."""
        pair = make_singlet_pair()
        assert pair.mouth_a.orientation == -pair.mouth_b.orientation

    def test_one_pass_flips(self):
        """One pass through throat flips orientation and wrap parity."""
        d = ThroatDefect(orientation=+1, wrap_parity=+1)
        d2 = d.one_pass_transport()
        assert d2.orientation == -1
        assert d2.wrap_parity == -1

    def test_two_pass_identity(self):
        """Two passes through throat = identity."""
        d = ThroatDefect(orientation=+1, wrap_parity=+1)
        d3 = d.two_pass_transport()
        assert d3.orientation == d.orientation
        assert d3.wrap_parity == d.wrap_parity

    def test_species_preserved(self):
        """Transport preserves species identity."""
        d = ThroatDefect(species_id="electron")
        assert d.one_pass_transport().species_id == "electron"

    def test_antipodal_quality(self):
        """Exact antipodal pair has quality ≈ 1."""
        pair = make_singlet_pair()
        assert pair.antipodal_quality > 0.99


# ════════════════════════════════════════════════════════════════════════════
# 7.  CAVITY PERSISTENCE
# ════════════════════════════════════════════════════════════════════════════

class TestCavityPersistence:
    """Cavity energy persists between steps."""

    def test_detector_settings_feed_branch_weights(self):
        """Detector settings change the cavity-derived branch weights."""
        bp = make_bell_pair()
        res_eq = bp.evolve_history(0.0, 0.0)
        res_sep = bp.evolve_history(0.0, np.pi / 2)
        assert abs(res_eq.wpi - res_sep.wpi) > 1e-4, (
            f"Detector settings did not change branch weights: "
            f"wpi_eq={res_eq.wpi}, wpi_sep={res_sep.wpi}"
        )

    def test_detector_history_generates_packets(self):
        """Detector-conditioned evolution should leave a packet/history trace."""
        bp = make_bell_pair()
        res = bp.evolve_history(0.0, np.pi / 4)
        assert (res.n_packets_0 + res.n_packets_pi) > 0, "Expected fired cavity packets"
        assert res.energy_final > 0.0, "Expected non-zero residual cavity energy"

    def test_cavity_stores_energy(self):
        """Driving a cavity mode deposits persistent energy."""
        cav = make_cavity((0, 1))
        for _ in range(100):
            cav.step({0: 0.1, 1: 0.0}, {0: 0.0, 1: 0.0}, dt=0.01)
        assert cav.energy() > 0, "Cavity should store energy"

    def test_cavity_energy_decays_slowly(self):
        """With no driving, energy decays at rate 2γ."""
        cav = make_cavity((0, 1))
        # Drive to build up amplitude
        for _ in range(200):
            cav.step({0: 0.1}, {}, dt=0.01)
        E_start = cav.energy()
        # Free ring-down
        for _ in range(100):
            cav.step({0: 0.0}, {}, dt=0.01)
        E_end = cav.energy()
        assert E_end < E_start, "Energy should decay"
        assert E_end > 0.1 * E_start, "Decay should be slow (τ_ring >> 1s)"

    def test_closure_check_structure(self):
        """Closure check returns (mismatch, branch, is_closed)."""
        mode = CavityMode(n=0, omega=1.0, b=0.1, bdot=0.0)
        mm, br, cl = mode.closure_check(tau_semi=np.pi, phi_spin=0, phi_throat=0)
        assert isinstance(mm, float)
        assert br in (0, 1)
        assert isinstance(cl, bool)


# ════════════════════════════════════════════════════════════════════════════
# 8.  CORRELATION CURVE
# ════════════════════════════════════════════════════════════════════════════

class TestCorrelationCurve:
    """Full correlation curve E(0, θ) = −cos(θ)."""

    def test_correlation_curve_matches_cosine(self):
        """E(0, θ) matches −cos(θ) over [0, 2π]."""
        thetas, E_vals, cos_ref = correlation_curve(n_points=50)
        assert np.allclose(E_vals, cos_ref, atol=1e-10), (
            f"Max deviation: {np.max(np.abs(E_vals - cos_ref)):.2e}"
        )
