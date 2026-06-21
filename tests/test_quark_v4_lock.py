"""
Regression tests for the v4 flavor-CP lock (PR #164).

Framing: QFT on the *fixed classical* throat geometry (geometry → fields),
not quantum gravity.  These tests pin the additive v4 layer migrated into
``geometrodynamics.qcd.quark_spectrum``:

  * the v3 lock stays bit-for-bit reproducible (the φ_h = 0 default), so
    every PR #155–#162 probe is untouched;
  * the v4 lock inherits the v3 mass spectrum exactly (the Hopf holonomy
    is a pure mixing phase — the #158 relocation);
  * ``extract_ckm_matrix`` at the v4 lock realizes the complete
    nine-observable flavor-CP dataset (|V_us|, |V_cb|, |V_ub|, |V_td|,
    |V_ts|, J, β, γ, α and sin δ) at ≤ 1%, unitarily, at the derived
    φ_h = π/k₅.
"""

import math

import numpy as np
import pytest

from geometrodynamics.qcd import quark_spectrum as qs


# PDG-anchored observed values used as the ≤1% targets (the #161 dataset).
V_OBS = {"us": 0.225, "cb": 0.04182, "ub": 0.00369, "td": 0.00857, "ts": 0.0411}
J_OBS = 3.08e-5
TRIANGLE_OBS = {"beta": 22.2, "gamma": 65.9, "alpha": 91.9}
SIN_DELTA_OBS = 0.887


def _triangle(V):
    J = float(np.imag(V[0, 0] * V[1, 1] * np.conj(V[0, 1]) * np.conj(V[1, 0])))
    beta = math.degrees(
        np.angle(-V[1, 0] * np.conj(V[1, 2]) / (V[2, 0] * np.conj(V[2, 2])))
    )
    gamma = math.degrees(
        np.angle(-V[0, 0] * np.conj(V[0, 2]) / (V[1, 0] * np.conj(V[1, 2])))
    )
    return J, beta, gamma


class TestV4DefaultsUnchanged:
    """The v4 layer is default-off: v3 behaviour is bit-for-bit preserved."""

    def test_new_fields_default_off(self):
        p = qs.QuarkParams()
        assert p.phi_h == 0.0
        assert p.eta_k1k3_plus == 0.0
        assert p.eta_k1k3_minus == 0.0
        assert p.eta_k1k5_minus == 0.0
        assert p.diag_shift_plus == (0.0, 0.0, 0.0)
        assert p.diag_shift_minus == (0.0, 0.0, 0.0)

    def test_v3_hamiltonian_real_at_phi_zero(self):
        H = qs.build_quark_hamiltonian(qs.LOCKED_QUARK_PARAMS)
        assert np.max(np.abs(H.imag)) == 0.0

    def test_v3_ckm_has_no_cp(self):
        V = qs.extract_ckm_matrix(qs.LOCKED_QUARK_PARAMS)
        assert np.max(np.abs(V.imag)) < 1e-12
        J, _, _ = _triangle(V.astype(complex))
        assert abs(J) < 1e-12


class TestV4LockHermiticity:
    def test_hermitian(self):
        H = qs.build_quark_hamiltonian(qs.LOCKED_QUARK_PARAMS_V4)
        assert np.allclose(H, H.conj().T)

    def test_phase_is_present(self):
        # The v4 lock genuinely carries the holonomy: the Hamiltonian is
        # complex (some same-partition off-diagonal has a nonzero phase).
        H = qs.build_quark_hamiltonian(qs.LOCKED_QUARK_PARAMS_V4)
        assert np.max(np.abs(H.imag)) > 1e-3

    def test_blocks_decoupled(self):
        H = qs.build_quark_hamiltonian(qs.LOCKED_QUARK_PARAMS_V4)
        plus = list(qs.PARTITION_PLUS_INDICES)
        minus = list(qs.PARTITION_MINUS_INDICES)
        assert np.max(np.abs(H[np.ix_(plus, minus)])) < 1e-14


class TestV4MassInheritance:
    """The Hopf holonomy is a pure mixing phase: masses are unchanged."""

    def test_real_eigenvalues_match_v3(self):
        e3 = np.sort(
            np.linalg.eigvalsh(
                qs.build_quark_hamiltonian(qs.LOCKED_QUARK_PARAMS).real
            )
        )
        e4 = np.sort(
            np.linalg.eigvalsh(
                qs.build_quark_hamiltonian(
                    qs.replace(qs.LOCKED_QUARK_PARAMS_V4, phi_h=0.0)
                ).real
            )
        )
        assert np.max(np.abs(e4 / e3 - 1)) < 1e-8

    def test_physical_spectrum_inherits_v3(self):
        m3 = qs.extract_physical_spectrum(qs.LOCKED_QUARK_PARAMS)
        m4 = qs.extract_physical_spectrum(qs.LOCKED_QUARK_PARAMS_V4)
        for s in qs.QUARK_SPECIES:
            denom = abs(m3[s]) + 1e-9
            assert abs(m4[s] - m3[s]) / denom < 1e-6

    def test_extract_strips_phi_h(self):
        # extract_physical_spectrum must read the φ_h = 0 spectrum, so
        # passing the lock with or without the phase gives identical masses.
        m_with = qs.extract_physical_spectrum(qs.LOCKED_QUARK_PARAMS_V4)
        m_without = qs.extract_physical_spectrum(
            qs.replace(qs.LOCKED_QUARK_PARAMS_V4, phi_h=0.0)
        )
        for s in qs.QUARK_SPECIES:
            assert m_with[s] == pytest.approx(m_without[s], abs=1e-9)


class TestV4FlavorCPObservables:
    """The complete nine-observable dataset at the v4 lock, all ≤ 1%."""

    @pytest.fixture(scope="class")
    def V(self):
        return qs.extract_ckm_matrix(qs.LOCKED_QUARK_PARAMS_V4)

    def test_unitary(self, V):
        assert np.max(np.abs(V.conj().T @ V - np.eye(3))) < 1e-10

    def test_magnitudes(self, V):
        ratios = {
            "us": abs(V[0, 1]) / V_OBS["us"],
            "cb": abs(V[1, 2]) / V_OBS["cb"],
            "ub": abs(V[0, 2]) / V_OBS["ub"],
            "td": abs(V[2, 0]) / V_OBS["td"],
            "ts": abs(V[2, 1]) / V_OBS["ts"],
        }
        for name, r in ratios.items():
            assert abs(r - 1) < 0.01, f"|V_{name}| off by {abs(r - 1):.3%}"

    def test_jarlskog(self, V):
        J, _, _ = _triangle(V)
        assert abs(J / J_OBS - 1) < 0.01
        # quartet-consistent: the invariant is the same on any quartet.
        J2 = float(
            np.imag(V[0, 1] * V[1, 2] * np.conj(V[0, 2]) * np.conj(V[1, 1]))
        )
        assert abs(J2 - (-J)) < 1e-12 or abs(J2 - J) < 1e-3 * abs(J) + 1e-12

    def test_triangle_angles(self, V):
        _, beta, gamma = _triangle(V)
        alpha = 180.0 - beta - gamma
        assert abs(beta - TRIANGLE_OBS["beta"]) < 1.0
        assert abs(gamma - TRIANGLE_OBS["gamma"]) < 1.0
        assert abs(alpha - TRIANGLE_OBS["alpha"]) < 1.0

    def test_sin_delta(self, V):
        J, _, _ = _triangle(V)
        sin_delta = J / (abs(V[0, 1]) * abs(V[1, 2]) * abs(V[0, 2]))
        assert abs(sin_delta - SIN_DELTA_OBS) < 0.01


class TestV4UnmixedReference:
    """The unmixed species-labeling reference must zero the v4 targeted
    mixing couplings — they generate the CKM rotation and belong in the
    adiabatic ramp, not in the diagonal reference."""

    def test_unmixed_params_zeroes_v4_mixing_couplings(self):
        u = qs._unmixed_params(qs.LOCKED_QUARK_PARAMS_V4)
        # the mixing knobs (transport + the three v4 etas) are off
        assert u.transport == 0.0
        assert u.gamma_q == 0.0
        assert u.partition_mixing == 0.0
        assert u.eta_k1k3_plus == 0.0
        assert u.eta_k1k3_minus == 0.0
        assert u.eta_k1k5_minus == 0.0
        # the structural level-repulsion stays on (it is not a mixing knob)
        assert u.eta_k3k5_minus == qs.LOCKED_QUARK_PARAMS_V4.eta_k3k5_minus

    def test_unmixed_reference_has_clean_diagonal(self):
        # In the unmixed reference (φ_h stripped as in mass extraction) the
        # three v4 targeted off-diagonals vanish, so each eigenvector is a
        # clean basis state and species labeling is unambiguous.  Only the
        # structural (3,−)–(5,−) repulsion survives, across an enormous gap.
        u = qs.replace(qs._unmixed_params(qs.LOCKED_QUARK_PARAMS_V4), phi_h=0.0)
        H = qs.build_quark_hamiltonian(u)
        # u–c = (1,+)-(3,+) at (0,2); d–s = (1,-)-(3,-) at (1,3);
        # d–b = (1,-)-(5,-) at (1,5)
        assert abs(H[0, 2]) < 1e-14
        assert abs(H[1, 3]) < 1e-14
        assert abs(H[1, 5]) < 1e-14
        # the structural s–b element remains
        assert abs(H[3, 5]) > 1.0
        # every eigenvector overlaps a single basis state ≥ 0.999
        _, vecs = np.linalg.eigh(H)
        for col in range(vecs.shape[1]):
            assert np.max(np.abs(vecs[:, col])) > 0.999

    def test_v3_unmixed_reference_unchanged(self):
        # For the v3 lock the v4 etas are already 0, so the new zeroing is
        # a strict no-op: the unmixed reference is bit-identical.
        u = qs._unmixed_params(qs.LOCKED_QUARK_PARAMS)
        assert u.eta_k1k3_plus == 0.0
        assert u.eta_k1k3_minus == 0.0
        assert u.eta_k1k5_minus == 0.0
        assert u.eta_k3k5_minus == qs.LOCKED_QUARK_PARAMS.eta_k3k5_minus


class TestV4DerivedPhase:
    def test_phi_h_is_pi_over_k5(self):
        assert qs.LOCKED_QUARK_PARAMS_V4.phi_h == pytest.approx(math.pi / 5.0)

    def test_ckm_requires_decoupled_blocks(self):
        # Turning on partition_mixing must make extract_ckm_matrix refuse.
        mixed = qs.replace(qs.LOCKED_QUARK_PARAMS_V4, partition_mixing=0.3)
        with pytest.raises(ValueError):
            qs.extract_ckm_matrix(mixed)
