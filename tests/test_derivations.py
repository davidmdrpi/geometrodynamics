"""
Tests for the derivation modules: transport, hopf phases, history/closure.

Validates the derivation chain:
  1. T = iσ_y derived from Hopf fibration (7 algebraic properties)
  2. Singlet state constructed from derived T matches Bell analyzer
  3. Bell phases from Hopf holonomy (baseline + setting-dependence)
  4. History closure reproduces singlet correlations (E, CHSH, no-signaling)
  5. Transaction histories conserve charges
"""

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geometrodynamics.embedding.transport import (
    derive_throat_transport,
    verify_transport_properties,
    verify_hopf_preservation,
    derive_singlet_from_transport,
)
from geometrodynamics.bell.analyzers import singlet_state
from geometrodynamics.bell.hopf_phases import (
    derived_phase_spin,
    geodesic_spin_phase,
    detector_holonomy_phase,
    compare_phase_formulas,
)
from geometrodynamics.history.closure import (
    make_bell_history,
    make_transaction_history,
    enumerate_bell_branches,
)


# ════════════════════════════════════════════════════════════════════════════
# 1.  THROAT TRANSPORT DERIVATION
# ════════════════════════════════════════════════════════════════════════════

class TestTransportDerivation:
    """T = iσ_y derived from Hopf fibration geometry."""

    def test_all_algebraic_properties(self):
        """T satisfies all 7 required algebraic properties."""
        T = derive_throat_transport()
        props = verify_transport_properties(T)
        for name, (err, passed) in props.items():
            assert passed, f"{name}: error = {err}"

    def test_hopf_fibration_preservation(self):
        """The orientation reversal σ preserves the Hopf bundle."""
        err = verify_hopf_preservation(n_samples=500)
        assert err < 1e-12, f"Hopf preservation error: {err}"

    def test_derived_singlet_matches_bell(self):
        """Singlet from Hopf → T → |Ψ⟩ matches Bell analyzer singlet."""
        psi_derived = derive_singlet_from_transport()
        psi_bell = singlet_state()
        overlap = abs(np.dot(psi_derived.conj(), psi_bell))
        assert abs(overlap - 1.0) < 1e-12, (
            f"Derived singlet overlap with Bell singlet: {overlap}"
        )

    def test_derived_singlet_is_antisymmetric(self):
        """Derived singlet is antisymmetric under particle exchange."""
        psi = derive_singlet_from_transport()
        psi_swap = np.array([psi[0], psi[2], psi[1], psi[3]])
        assert np.allclose(psi_swap, -psi)


# ════════════════════════════════════════════════════════════════════════════
# 2.  HOPF PHASE DERIVATION
# ════════════════════════════════════════════════════════════════════════════

class TestHopfPhases:
    """Bell phases derived from Hopf holonomy."""

    def test_antipodal_baseline(self):
        """Equal settings + antipodal pair → phase = π/2."""
        phi = derived_phase_spin(0.0, 0.0, theta_transport=np.pi)
        assert abs(phi - np.pi / 2) < 1e-10, f"baseline = {phi}"

    def test_geodesic_component(self):
        """Geodesic spin phase = α_spin × θ = π/2 for spin-½ antipodal."""
        assert abs(geodesic_spin_phase(np.pi, 0.5) - np.pi / 2) < 1e-12

    def test_detector_holonomy_zero_at_equal_settings(self):
        """Equal detector settings → zero holonomy difference."""
        for theta in [0, np.pi / 4, np.pi / 2, np.pi]:
            dh = detector_holonomy_phase(theta, theta)
            assert abs(dh) < 1e-12, f"θ={theta}: Δφ_det = {dh}"

    def test_detector_holonomy_nonzero_at_different_settings(self):
        """Different detector settings → nonzero holonomy difference."""
        dh = detector_holonomy_phase(0.0, np.pi / 2)
        assert abs(dh) > 0.1, f"Δφ_det = {dh}, should be nonzero"

    def test_settings_change_derived_phase(self):
        """Different settings produce different derived phases."""
        phi_eq = derived_phase_spin(0.0, 0.0)
        phi_sep = derived_phase_spin(0.0, np.pi / 2)
        assert abs(phi_eq - phi_sep) > 0.1, (
            f"eq={phi_eq}, sep={phi_sep}, should differ"
        )


# ════════════════════════════════════════════════════════════════════════════
# 3.  HISTORY CLOSURE — BELL
# ════════════════════════════════════════════════════════════════════════════

class TestHistoryBell:
    """Bell correlations from the history closure framework."""

    def test_bell_history_structure(self):
        """Bell history has 2 events, 1 worldline, conserves charge."""
        h = make_bell_history(0, 0, +1, -1)
        cl = h.check_closure()
        assert cl.n_events == 2
        assert cl.n_worldlines == 1
        assert cl.conservation_error < 1e-10

    def test_bell_history_closes(self):
        """Bell history satisfies phase closure."""
        cl = make_bell_history(0, 0, +1, -1).check_closure()
        assert cl.is_closed, f"Not closed: mismatch={cl.phase_mismatch}"

    def test_four_branches_normalised(self):
        """Four outcome branches sum to probability 1."""
        bs = enumerate_bell_branches(0, 0)
        assert len(bs) == 4
        total = sum(b.probability for b in bs)
        assert abs(total - 1.0) < 1e-10, f"Total = {total}"

    def test_equal_settings_anticorrelation(self):
        """E(0, 0) = −1 from history branches."""
        bs = enumerate_bell_branches(0, 0)
        E = sum(b.outcome_a * b.outcome_b * b.probability for b in bs)
        assert abs(E - (-1.0)) < 1e-10, f"E(0,0) = {E}"

    def test_correlation_is_minus_cosine(self):
        """E(a, b) = −cos(a−b) from history branches."""
        for a, b in [(0, np.pi / 4), (0, np.pi / 2), (np.pi / 3, np.pi / 6)]:
            bs = enumerate_bell_branches(a, b)
            E = sum(x.outcome_a * x.outcome_b * x.probability for x in bs)
            expected = -np.cos(a - b)
            assert abs(E - expected) < 1e-10, (
                f"E({a:.3f},{b:.3f}) = {E:.6f}, expected {expected:.6f}"
            )

    def test_chsh_from_history(self):
        """CHSH = 2√2 from history branch enumeration."""
        def hE(a, b):
            return sum(
                x.outcome_a * x.outcome_b * x.probability
                for x in enumerate_bell_branches(a, b)
            )

        S = abs(
            hE(0, np.pi / 4)
            + hE(0, -np.pi / 4)
            + hE(np.pi / 2, np.pi / 4)
            - hE(np.pi / 2, -np.pi / 4)
        )
        assert abs(S - 2 * np.sqrt(2)) < 1e-10, f"S = {S}"

    def test_no_signaling_from_history(self):
        """Marginals = ½ independent of remote setting."""
        for b1, b2 in [(0, np.pi / 4), (np.pi / 6, np.pi / 2)]:
            p1 = sum(
                x.probability
                for x in enumerate_bell_branches(0, b1)
                if x.outcome_a == +1
            )
            p2 = sum(
                x.probability
                for x in enumerate_bell_branches(0, b2)
                if x.outcome_a == +1
            )
            assert abs(p1 - 0.5) < 1e-10, f"P(+|b={b1}) = {p1}"
            assert abs(p1 - p2) < 1e-10, (
                f"P(+|b={b1}) = {p1} ≠ P(+|b={b2}) = {p2}"
            )


# ════════════════════════════════════════════════════════════════════════════
# 4.  HISTORY CLOSURE — TRANSACTIONS
# ════════════════════════════════════════════════════════════════════════════

class TestHistoryTransaction:
    """Transaction histories satisfy closure and conservation."""

    def test_transaction_structure(self):
        """Transaction has 2 events, 2 worldlines."""
        p4_a = np.array([0, 0, 0, 1.0])
        p4_b = np.array([0, 0, 0, -1.0])
        h = make_transaction_history(p4_a, p4_b, 1.0, -1.0, 0, 1)
        cl = h.check_closure()
        assert cl.n_events == 2
        assert cl.n_worldlines == 2

    def test_transaction_conservation(self):
        """Transaction conserves charge."""
        p4_a = np.array([0, 0, 0, 1.0])
        p4_b = np.array([0, 0, 0, -1.0])
        cl = make_transaction_history(
            p4_a, p4_b, 1.0, -1.0, 0, 1
        ).check_closure()
        assert cl.conservation_error < 1e-10
