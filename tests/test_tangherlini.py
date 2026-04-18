"""Tests for 5D Tangherlini eigenmodes and Maxwell solver."""

import numpy as np
import pytest

from geometrodynamics.tangherlini import (
    solve_radial_modes,
    solve_maxwell_from_eigenmode,
)
from geometrodynamics.tangherlini.alpha_q import derive_alpha_q, throat_du_dr


@pytest.fixture
def computed_modes():
    """Pre-compute modes for l=1,3,5."""
    modes = {}
    for l in [1, 3, 5]:
        oms, fns, rg = solve_radial_modes(l)
        modes[l] = {"omega": oms, "funcs": fns, "r_grid": rg}
    return modes


class TestRadialModes:
    def test_positive_eigenfrequencies(self, computed_modes):
        for l, data in computed_modes.items():
            assert all(w > 0 for w in data["omega"])

    def test_mode_ordering(self, computed_modes):
        """ω₀ < ω₁ < ω₂ ..."""
        for l, data in computed_modes.items():
            oms = data["omega"]
            for i in range(len(oms) - 1):
                assert oms[i] < oms[i + 1]

    def test_dirichlet_bc(self, computed_modes):
        """u(R_MID) ≈ 0 at the throat (first grid point)."""
        for l, data in computed_modes.items():
            fn = data["funcs"][0]
            assert abs(fn["u_half"][0]) < 1e-6


class TestAlphaQ:
    def test_alpha_q_reference_mode(self, computed_modes):
        """α_q(1, 0) = ±1 by construction."""
        table = derive_alpha_q(computed_modes)
        assert abs(abs(table[(1, 0)]) - 1.0) < 1e-9

    def test_alpha_q_higher_modes(self, computed_modes):
        """Higher-l modes have well-defined ratios."""
        table = derive_alpha_q(computed_modes)
        for key, val in table.items():
            assert np.isfinite(val)


class TestMaxwellSolver:
    def test_coulomb_validation(self, computed_modes):
        """Solved field matches Q/r to relative error < 1e-6."""
        result = solve_maxwell_from_eigenmode(computed_modes)
        assert result["rel_err"] < 1e-6

    def test_positive_charge(self, computed_modes):
        """Extracted Q is positive."""
        result = solve_maxwell_from_eigenmode(computed_modes)
        assert result["Q"] > 0

    def test_linsys_residual(self, computed_modes):
        """Linear-system residual ||L·A − rhs||/||rhs|| is small.

        This replaces the previous PDE residual test, which was misleading:
        a raw FD Laplacian applied to the solved (or even exact) field
        grows with N near the throat due to stencil cancellation, not
        solver inaccuracy.  The linear-system residual is the correct
        acceptance metric for the sparse BVP solve.
        """
        result = solve_maxwell_from_eigenmode(computed_modes)
        assert result["linsys_res"] < 1e-6

    def test_electric_field_matches_coulomb(self, computed_modes):
        """E_r ≈ Q/r² in the interior of the radial domain.

        The analytic Coulomb potential is A(r) = Q/r − Q/R_OUTER,
        giving E_r = −dA/dr = Q/r².  We check this away from both
        boundaries (throat and outer edge) where finite-difference
        stencil edge effects are largest.

        This is the primary regression test for the E_r field that
        the transaction subsystem consumes via advanced_confirm_amplitude.
        """
        result = solve_maxwell_from_eigenmode(computed_modes)
        r = result["r"]
        E_r = result["E_r"]
        Q = result["Q"]
        N = len(r)

        # Interior band: skip 10% at each end to avoid boundary stencils
        lo = N // 10
        hi = N - N // 10
        r_int = r[lo:hi]
        E_r_int = E_r[lo:hi]
        E_r_exact = Q / r_int ** 2

        # Relative error in interior band
        rel_err = np.abs(E_r_int - E_r_exact) / E_r_exact
        max_rel_err = float(rel_err.max())
        mean_rel_err = float(rel_err.mean())

        assert max_rel_err < 1e-4, (
            f"E_r vs Q/r² max relative error {max_rel_err:.2e} in interior "
            f"(expected < 1e-4)"
        )
        assert mean_rel_err < 1e-5, (
            f"E_r vs Q/r² mean relative error {mean_rel_err:.2e} in interior "
            f"(expected < 1e-5)"
        )

    def test_electric_field_positive(self, computed_modes):
        """E_r is positive throughout the domain (repulsive Coulomb)."""
        result = solve_maxwell_from_eigenmode(computed_modes)
        # Skip first point (throat boundary) where the one-sided stencil
        # can produce a small numerical artifact
        assert np.all(result["E_r"][1:] > 0), "E_r should be positive (Q > 0)"

    def test_electric_field_monotone_decreasing(self, computed_modes):
        """E_r = Q/r² is strictly decreasing for r > R_MID."""
        result = solve_maxwell_from_eigenmode(computed_modes)
        E_r = result["E_r"]
        # Interior only (skip boundary stencil artifacts at both ends)
        N = len(E_r)
        lo, hi = N // 10, N - N // 10
        diffs = np.diff(E_r[lo:hi])
        assert np.all(diffs < 0), "E_r should decrease monotonically in interior"
