"""
Tests for ``geometrodynamics.qcd.quark_spectrum``.

Covers the structural and regression criteria from v3 spec §7:

    1. Hermiticity of the Hamiltonian.
    2. Correct shape and basis ordering.
    3. Lepton-limit reduction (§7 criterion 4).
    4. Color independence (§7 criterion 5, structural).
    5. Partition-class sign convention.
    6. Basis-to-species mapping consistency with u_q(k) = k - 2.

The numerical-mass tests (criteria 2 and 3 of §7) are omitted here
because they require ``LOCKED_QUARK_PARAMS`` to be populated by the
calibration scripts.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.qcd import quark_spectrum as qs


# ════════════════════════════════════════════════════════════════════════
# STRUCTURAL TESTS
# ════════════════════════════════════════════════════════════════════════

class TestHamiltonianStructure:
    """Shape, Hermiticity, and basis-ordering invariants."""

    def test_default_params_shape(self):
        H = qs.build_quark_hamiltonian()
        assert H.shape == (6, 6)

    def test_hermiticity_default(self):
        H = qs.build_quark_hamiltonian()
        np.testing.assert_allclose(H, H.conj().T, atol=1e-14)

    def test_hermiticity_with_partition_mixing(self):
        params = qs.QuarkParams(partition_mixing=0.3, phase=0.01)
        H = qs.build_quark_hamiltonian(params)
        np.testing.assert_allclose(H, H.conj().T, atol=1e-14)

    def test_basis_states_length(self):
        assert len(qs.BASIS_STATES) == 6

    def test_basis_states_unique(self):
        assert len(set(qs.BASIS_STATES)) == 6

    def test_pass_counts_odd_only(self):
        for k, _ in qs.BASIS_STATES:
            assert k % 2 == 1, f"pass count {k} is not odd"

    def test_partition_classes_complete(self):
        classes = {p for _, p in qs.BASIS_STATES}
        assert classes == {"+", "-"}

    def test_basis_to_species_bijective(self):
        assert set(qs.BASIS_TO_SPECIES.keys()) == set(qs.BASIS_STATES)
        assert set(qs.BASIS_TO_SPECIES.values()) == set(qs.QUARK_SPECIES)


# ════════════════════════════════════════════════════════════════════════
# PARAMETER AND u_q-FORM TESTS
# ════════════════════════════════════════════════════════════════════════

class TestMinimalAnsatz:
    """§4: u_q(k) = k - 2 minimal topological ansatz."""

    def test_u_q_at_pass_counts(self):
        assert qs._u_q(1, "k_minus_2") == -1.0
        assert qs._u_q(3, "k_minus_2") == +1.0
        assert qs._u_q(5, "k_minus_2") == +3.0

    def test_u_q_zero_form(self):
        for k in qs.PASS_COUNTS:
            assert qs._u_q(k, "zero") == 0.0

    def test_u_q_unknown_form_raises(self):
        with pytest.raises(ValueError, match="Unknown u_q form"):
            qs._u_q(1, "not_a_real_form")

    def test_species_ordering_follows_minimal_ansatz(self):
        """
        With u_q(k) = k - 2 and γ_q > 0, within each pass count the
        lighter partition class should match the canonical
        BASIS_TO_SPECIES assignments (u<d, s<c, b<t).
        """
        params = qs.QuarkParams(gamma_q=1.0, partition_mixing=0.0)
        H = qs.build_quark_hamiltonian(params)
        diag = np.real(np.diag(H))

        # Indices into BASIS_STATES
        idx = {state: i for i, state in enumerate(qs.BASIS_STATES)}

        # At k=1: (1,+)=u should be lighter than (1,-)=d
        assert diag[idx[(1, "+")]] < diag[idx[(1, "-")]]
        # At k=3: (3,-)=s should be lighter than (3,+)=c
        assert diag[idx[(3, "-")]] < diag[idx[(3, "+")]]
        # At k=5: (5,-)=b should be lighter than (5,+)=t
        assert diag[idx[(5, "-")]] < diag[idx[(5, "+")]]


# ════════════════════════════════════════════════════════════════════════
# LEPTON-LIMIT REGRESSION (v3 §7 criterion 4)
# ════════════════════════════════════════════════════════════════════════

class TestLeptonLimit:
    """Single-topology-continuity: quark → lepton reduction."""

    def test_block_diagonal_structure(self):
        """
        With γ_q = 0 and partition_mixing = 0, the Hamiltonian must
        block-diagonalize into two decoupled 3×3 blocks, one per
        partition class.
        """
        params = qs.QuarkParams(
            action_base=2.0 * math.pi,
            gamma_q=0.0,
            partition_mixing=0.0,
        )
        H = qs.build_quark_hamiltonian(params)

        # Basis ordering is (1,+),(1,-),(3,+),(3,-),(5,+),(5,-).
        # The two 3x3 blocks are the "+"-partition states at indices
        # {0, 2, 4} and the "-"-partition states at indices {1, 3, 5}.
        # In this ordering, H is NOT block-diagonal as a contiguous
        # matrix, but it should be a permuted block-diagonal.
        plus_indices = [0, 2, 4]
        minus_indices = [1, 3, 5]

        # All cross-partition entries should be zero
        for i in plus_indices:
            for j in minus_indices:
                assert abs(H[i, j]) < 1e-14, (
                    f"nonzero cross-partition H[{i},{j}] = {H[i, j]} "
                    f"with γ_q=0 and partition_mixing=0"
                )

    def test_eigenvalue_pairing(self):
        """
        The two 3×3 blocks are identical (same diagonal structure, same
        off-diagonal structure), so their eigenvalues must pair up.
        """
        params = qs.QuarkParams(
            action_base=2.0 * math.pi,
            gamma_q=0.0,
            partition_mixing=0.0,
        )
        H = qs.build_quark_hamiltonian(params)
        eigvals = np.linalg.eigvalsh(H)
        for i in range(3):
            assert abs(eigvals[2 * i + 1] - eigvals[2 * i]) < 1e-12

    @pytest.mark.skipif(
        not pytest.importorskip(
            "geometrodynamics.tangherlini",
            reason="lepton_spectrum not importable; run in installed package.",
        ),
        reason="Requires geometrodynamics.tangherlini.solved_lepton_masses_mev",
    )
    def test_lepton_limit_check_ratios(self):
        """
        Run the public ``quark_lepton_limit_check`` against the live
        lepton module and assert ratio equality to the documented tol.

        NOTE: This test asserts RATIO equality, which requires only that
        the residual continuous knobs (phase, transport, pinhole,
        resistance, β) hold the same values in both modules. If the
        lepton module uses values different from the defaults in
        ``QuarkParams``, this test will fail with a message indicating
        which knobs must be synchronized.
        """
        result = qs.quark_lepton_limit_check(tol=1e-6)
        assert result["passed"], (
            f"Lepton-limit check failed. Details: {result}"
        )


# ════════════════════════════════════════════════════════════════════════
# COLOR-INDEPENDENCE STRUCTURAL CHECK (v3 §7 criterion 5)
# ════════════════════════════════════════════════════════════════════════

class TestColorIndependence:

    def test_color_check_reports_pass(self):
        result = qs.color_independence_check()
        assert result["passed"] is True
        assert result["color_multiplicity"] == 3

    def test_hamiltonian_signature_has_no_color_argument(self):
        """
        Color is a degeneracy label, not a Hamiltonian parameter.
        ``build_quark_hamiltonian`` must accept only ``params``.
        """
        import inspect
        sig = inspect.signature(qs.build_quark_hamiltonian)
        param_names = set(sig.parameters.keys())
        # Strictly: only ``params``. Stricter than needed, but documents intent.
        assert "color" not in param_names
        assert "c1" not in param_names
        assert "su3" not in param_names


# ════════════════════════════════════════════════════════════════════════
# CALIBRATION-STATE GUARDS
# ════════════════════════════════════════════════════════════════════════

class TestCalibrationState:
    """
    ``solved_quark_masses_mev`` must raise a helpful error until the
    calibration scripts have populated LOCKED_QUARK_PARAMS.
    """

    def test_uncalibrated_raises(self):
        # Ensure a clean uncalibrated state
        original = qs.LOCKED_QUARK_PARAMS
        try:
            qs.LOCKED_QUARK_PARAMS = None
            with pytest.raises(NotImplementedError, match="not yet calibrated"):
                qs.solved_quark_masses_mev()
        finally:
            qs.LOCKED_QUARK_PARAMS = original
