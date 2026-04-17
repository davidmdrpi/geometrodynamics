"""Tests for Hopf fibration geometry."""

import numpy as np
import pytest

from geometrodynamics.hopf import (
    hopf_connection,
    hopf_curvature,
    hopf_holonomy,
    hopf_circle,
    compute_c1,
    compute_spinor_monodromy,
)


class TestHopfConnection:
    def test_equatorial_vanishing(self):
        """A(χ=π/2) = 0: zero self-energy at equator."""
        assert abs(hopf_connection(np.pi / 2)) < 1e-15

    def test_polar_values(self):
        """A(0) = 0.5, A(π) = -0.5."""
        assert abs(hopf_connection(0.0) - 0.5) < 1e-15
        assert abs(hopf_connection(np.pi) + 0.5) < 1e-15

    def test_curvature_maximum_at_equator(self):
        """|F| is maximum at χ = π/2."""
        chi = np.linspace(0, np.pi, 1000)
        F = hopf_curvature(chi)
        assert np.argmax(F) == pytest.approx(500, abs=5)


class TestChernNumber:
    def test_c1_equals_one(self):
        """First Chern number |c₁| = 1 to high precision."""
        result = compute_c1()
        assert result["err_abs"] < 1e-9

    def test_c1_orientation(self):
        """c₁ in dχ∧dφ orientation is −1."""
        result = compute_c1()
        assert abs(result["c1_chiphi"] - (-1.0)) < 1e-9


class TestSpinorMonodromy:
    def test_sign_flip_at_2pi(self):
        """⟨ψ₀|U(2π)|ψ₀⟩ = −1 (spin-½ sign flip)."""
        result = compute_spinor_monodromy()
        assert result["signflip_err"] < 1e-10

    def test_return_at_4pi(self):
        """⟨ψ₀|U(4π)|ψ₀⟩ = +1 (full return)."""
        result = compute_spinor_monodromy()
        assert result["return_err"] < 1e-10

    def test_u1_phase_at_pole(self):
        """Hopf holonomy at χ=0 gives e^{iπ} = −1."""
        result = compute_spinor_monodromy()
        assert abs(result["u1_phase_2pi"] + 1.0) < 1e-10


class TestHopfCircle:
    def test_circle_shape(self):
        """Hopf circle returns three coordinate arrays of correct length."""
        x, y, z = hopf_circle(np.pi / 4, 0.0, N=64)
        assert len(x) == 64
        assert len(y) == 64
        assert len(z) == 64
