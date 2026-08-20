"""Tests for the initial-data geometry round.

The load-bearing one is
`test_the_constraint_eigenvalue_is_the_esu_dipole_level`: the whole reason this
calculation needs a resolved neck is that ``λ = 1 + R̂/2`` lands exactly on
``(n+1)²`` at ``n = 1``.  If that coincidence were not exact the round would
have a different shape, so it is pinned rather than described.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.finite_mouth import FOUR_PI
from geometrodynamics.waves.initial_data import (
    CONSTRAINT_EIGENVALUE, areal_response_from_conformal_factor,
    constraint_operator_eigenvalue,
    measure_removing_the_balls_lifts_the_degeneracy,
    measure_the_constraint_operator_is_degenerate_on_the_sphere,
)
from geometrodynamics.waves.neck import (
    exterior_dtn, exterior_dtn_monopole, exterior_log_derivative,
)


# ── the constraint ─────────────────────────────────────────────────────────
def test_the_constraint_eigenvalue_is_the_esu_dipole_level():
    """``λ = 1 + R̂/2 = 4 = (n+1)²`` at ``n = 1`` — exactly, not nearly.

    This is the coincidence the round turns on.  The Hamiltonian constraint's
    operator is degenerate with the ambient's dipole level, so on the full
    sphere it has a kernel and the fixed-ambient Green function has a literal
    pole at it.
    """
    assert constraint_operator_eigenvalue(6.0) == pytest.approx(4.0, abs=1e-15)
    assert CONSTRAINT_EIGENVALUE == pytest.approx((1 + 1) ** 2, abs=1e-15)


def test_the_operator_has_a_kernel_at_k_one_and_is_indefinite():
    m = measure_the_constraint_operator_is_degenerate_on_the_sphere()
    assert m["kernel_modes"] == [1]
    assert m["the_operator_has_a_kernel"]
    assert m["the_operator_is_indefinite"]
    assert m["it_sits_on_the_esu_dipole_level"]
    for row in m["rows"]:
        k = row["k"]
        assert row["operator_eigenvalue"] == pytest.approx(3 - k * (k + 2))
        assert row["free_esu_lambda"] == pytest.approx((k + 1) ** 2)


def test_the_areal_response_is_four_u():
    """``g = ψ⁴ĝ`` means areas scale as ``ψ⁴``, so ``δA/A = 4u``.

    The sign convention is what the round reports, so it is pinned: ``u < 0``
    shrinks the mouth and moves the geometry *toward* a neck.
    """
    assert areal_response_from_conformal_factor(0.25) == pytest.approx(1.0)
    assert float(areal_response_from_conformal_factor(-0.1)) < 0.0
    got = areal_response_from_conformal_factor(np.array([0.0, 0.5, -0.5]))
    assert np.allclose(got, [0.0, 2.0, -2.0])


# ── the neck lifts the degeneracy ──────────────────────────────────────────
def test_removing_the_balls_gives_four_pi_a_cubed():
    """``N₀(a, 4) = 4πa³`` — the closed form for the lifting.

    ``4πa`` at ``λ = 0`` against ``4πa³`` here: the extra ``a²`` *is* the
    degeneracy, and it sends the stiffness to zero in the point limit.
    """
    m = measure_removing_the_balls_lifts_the_degeneracy()
    assert m["it_is_four_pi_a_cubed"]
    assert m["it_is_nonzero_everywhere"]
    assert m["it_vanishes_in_the_point_limit"]
    for row in m["rows"]:
        assert row["over_a_cubed"] == pytest.approx(FOUR_PI, rel=2e-3)


def test_the_stiffness_is_positive_and_far_smaller_than_at_zero():
    for a in (0.02, 0.05, 0.2):
        n4 = exterior_dtn_monopole(a, CONSTRAINT_EIGENVALUE)
        n0 = exterior_dtn_monopole(a, 0.0)
        assert n4 > 0.0
        assert n4 < n0
        # the ratio is a^2 up to the 4pi's, which is the degeneracy showing
        assert n4 / n0 == pytest.approx(a ** 2, rel=0.15)


def test_the_higher_multipoles_are_not_degenerate():
    """The softness is confined to the one channel the coincidence produces —
    the ``cos χ`` member of the ``k = 1`` multiplet about the mouth."""
    for a in (0.05, 0.2):
        values = [exterior_dtn(a, CONSTRAINT_EIGENVALUE, l)
                  for l in (0, 1, 2, 3)]
        assert values[0] < 0.2 * values[1]
        assert all(v > 0.5 for v in values[1:])
        assert all(b > c for c, b in zip(values, values[1:]))


def test_the_shooting_integrator_still_agrees_at_the_eigenvalue():
    """It agrees, but far less well than PR #262's ``1e-14`` — and that is
    informative rather than alarming.

    ``N₀ ~ 4πa³`` is a near-cancellation of two ``O(a)`` terms at the degenerate
    eigenvalue, so the shooting solve loses digits exactly where the operator
    goes soft.  ``2e-06`` is the price of the degeneracy, not a solver defect.
    """
    for a in (0.02, 0.05, 0.15):
        shot = -(FOUR_PI * math.sin(a) ** 2) * exterior_log_derivative(
            a, CONSTRAINT_EIGENVALUE, 0)
        closed = exterior_dtn_monopole(a, CONSTRAINT_EIGENVALUE)
        assert shot == pytest.approx(closed, rel=1e-4)


def test_the_closed_form_comes_from_the_cos_chi_mode():
    """``N₀ = −4π sin²a · ψ'/ψ`` with ``ψ = cos χ`` gives ``4π sin²a tan a``.

    Checking the *mechanism* rather than the fitted coefficient: the monopole
    member of the ``k = 1`` multiplet about the mouth's centre is ``cos χ``, and
    it reproduces the measured stiffness at every radius including where
    ``4πa³`` has started to drift.
    """
    for a in (0.01, 0.05, 0.2, 0.35):
        predicted = FOUR_PI * math.sin(a) ** 2 * math.tan(a)
        assert exterior_dtn_monopole(a, CONSTRAINT_EIGENVALUE) == pytest.approx(
            predicted, rel=1e-9)
