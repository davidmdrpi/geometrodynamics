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
    CONSTRAINT_EIGENVALUE, KERNEL_PROJECTOR,
    areal_response_from_conformal_factor,
    constraint_operator_eigenvalue, kernel_component, regularised_green,
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


def test_the_dipole_partners_are_degenerate_too():
    """The ``k = 1`` multiplet has **four** members, and they do not match.

    About a mouth it splits into one ``ℓ = 0`` member ``cos χ`` and three
    ``ℓ = 1`` members ``sin χ``.  Both are degenerate at ``λ = 4``, but

        ``N₀ = +4π sin²a·tan a ~ +4πa³`` ,   ``N₁ = −2π sin 2a ~ −4πa``

    — opposite signs and two orders apart in ``a``.  A first version of this
    test asserted the ``ℓ ≥ 1`` channels were positive and ``O(1)``, i.e. not
    degenerate at all.  That was `neck` solving ``v'' + [λ − ℓ(ℓ+2)/sin²χ]v = 0``
    with the ``S³`` harmonic eigenvalue where the angular Laplacian on the
    two-sphere of directions gives ``ℓ(ℓ+1)`` — so **the test was validating the
    wrong differential equation**, and the three dipole partners of the
    degeneracy were hidden by it.
    """
    for a in (0.05, 0.2):
        values = [exterior_dtn(a, CONSTRAINT_EIGENVALUE, l)
                  for l in (0, 1, 2, 3)]
        assert values[0] > 0.0                     # monopole: positive, tiny
        assert values[1] < 0.0                     # dipoles: NEGATIVE
        assert values[2] > 0.5 and values[3] > 0.5  # l >= 2: genuinely free
        assert abs(values[0]) < abs(values[1])      # and two orders softer


def test_the_dipole_closed_form_is_minus_two_pi_sin_two_a():
    """``ψ = sin χ`` gives ``N₁ = −4π sin²a·cot a = −2π sin 2a``, exactly.

    This closed form did not exist before the ``ℓ(ℓ+1)`` correction, and it is
    the check that would have caught the bug: with ``ℓ(ℓ+2)`` the shooting solve
    misses it outright rather than by a tolerance.
    """
    for a in (0.02, 0.05, 0.15, 0.35):
        predicted = -2.0 * math.pi * math.sin(2.0 * a)
        assert exterior_dtn(a, CONSTRAINT_EIGENVALUE, 1) == pytest.approx(
            predicted, rel=1e-5)


def test_the_radial_equation_carries_the_two_sphere_eigenvalue():
    """The bug itself, pinned by exhibiting solutions rather than re-deriving.

    At ``λ = 4`` the ``S³`` harmonics ``x^A`` split about a point into
    ``ψ = cos χ`` (``ℓ = 0``) and ``ψ = sin χ`` (``ℓ = 1``), so ``v = ψ sin χ``
    must solve ``v'' + [4 − L²/sin²χ]v = 0``.  With ``L² = ℓ(ℓ+1)`` both do;
    with ``L² = ℓ(ℓ+2)`` the dipole leaves a residual of exactly ``−1``.
    """
    h = 1e-6

    def residual(v, l2, chi):
        return ((v(chi + h) - 2 * v(chi) + v(chi - h)) / h ** 2
                + (4.0 - l2 / math.sin(chi) ** 2) * v(chi))

    for chi in (0.4, 1.0, 1.9, 2.7):
        assert abs(residual(lambda x: math.cos(x) * math.sin(x), 0.0,
                            chi)) < 1e-3
        assert abs(residual(lambda x: math.sin(x) ** 2, 2.0, chi)) < 1e-3
        assert abs(residual(lambda x: math.sin(x) ** 2, 3.0, chi)) \
            == pytest.approx(1.0, rel=1e-3)


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


# ── the regularised Green function ─────────────────────────────────────────
def test_the_regularised_green_function_inverts_on_the_complement():
    """``(∇² + 3)G_⊥ = −δ³ + (2/π²)cos χ`` — the identity, checked away from
    the origin.

    The residual must be **exactly** the normalised kernel projector: the
    ``k = 1`` harmonics are ``x^A`` with ``‖x^A‖² = π²/2`` and
    ``Σ_A x^Ay^A = cos χ``, so the projector kernel is ``2cos χ/π²``.  A
    constant ratio is what says the leftover is the kernel and nothing else.
    """
    h = 1e-5
    for chi in (0.3, 0.8, 1.5, 2.2, 2.9):
        second = (regularised_green(chi + h) - 2 * regularised_green(chi)
                  + regularised_green(chi - h)) / h ** 2
        first = (regularised_green(chi + h)
                 - regularised_green(chi - h)) / (2 * h)
        lhs = float(second + 2.0 / math.tan(chi) * first
                    + 3.0 * regularised_green(chi))
        assert lhs / math.cos(chi) == pytest.approx(KERNEL_PROJECTOR, rel=1e-4)


def test_the_green_function_has_the_right_short_distance_normalisation():
    """``G_⊥ · 4πχ → 1`` — otherwise it is not a Green function in 3D."""
    for chi in (1e-3, 1e-4, 1e-5):
        assert float(regularised_green(chi)) * 4.0 * math.pi * chi == \
            pytest.approx(1.0, rel=1e-3)


def test_the_green_function_is_finite_at_the_antipode():
    """``G_⊥(π) = 1/(4π²)`` exactly — this is the ``k = 1`` removal working.

    The unregularised ``G`` diverges there; removing the degenerate mode is
    what makes the antipode finite.
    """
    for eps in (1e-5, 1e-7, 1e-9):
        assert float(regularised_green(math.pi - eps)) == pytest.approx(
            1.0 / (4.0 * math.pi ** 2), rel=1e-6)


def test_the_kernel_projector_is_two_over_pi_squared():
    assert KERNEL_PROJECTOR == pytest.approx(2.0 / math.pi ** 2, rel=1e-15)
    # normalisation of the k=1 harmonics on the unit S^3: ||x^A||^2 = pi^2/2
    rng = np.random.default_rng(3)
    pts = rng.normal(size=(200000, 4))
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    mean_sq = float(np.mean(pts[:, 0] ** 2))
    volume = 2.0 * math.pi ** 2
    assert mean_sq * volume == pytest.approx(math.pi ** 2 / 2.0, rel=0.02)


def test_the_kernel_component_detects_a_dipole_source():
    """The obstruction has to be *measurable*, or the solvability condition is
    decoration."""
    rng = np.random.default_rng(11)
    pts = rng.normal(size=(40000, 4))
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    w = np.full(len(pts), 2.0 * math.pi ** 2 / len(pts))
    # a pure dipole source: overlap must be nonzero
    dipole = kernel_component(pts[:, 2], w, pts)
    assert abs(dipole[2]) > 0.5 * math.pi ** 2 / 2.0
    # an even (monopole-like) source: overlap must vanish
    even = kernel_component(np.ones(len(pts)), w, pts)
    assert np.max(np.abs(even)) < 0.05 * math.pi ** 2 / 2.0
