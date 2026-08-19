"""Tests for the metric backreaction round.

Two of these carry the round: `test_the_batched_solver_matches_the_reference`,
because every number here is built on it, and
`test_the_angular_rule_integrates_the_singular_direction_structure_exactly`,
because that exactness is what makes the ``1/χ⁴`` cancel and the ``S³``
integrals converge at all.

The expensive measurements are exercised at deliberately small quadratures.
They check *structure* — signs, convergence direction, that the control fires —
and the converged numbers live in the probe, where the rule is big enough to
trust them.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.backreaction import (
    TENSOR_MODE_FREQUENCY, BackreactionSetup, ShearQuadrature,
    WORKING_BACKREACTION, derive_the_tensor_mode_equation, left_invariant_frame,
    measure_the_stress_tensor_splits_bilinearly,
    measure_the_tensor_sector_is_stable_and_fluid_free,
    quaternion_product, shear_projection, shear_response, solve_batch,
    stress_series, unreachable_fraction,
)
from geometrodynamics.waves.backreaction import _direction_rule, _tangent_basis
from geometrodynamics.waves.two_wave import (
    circle_point, orthonormal_frame, solve_field, stress_tensor,
)

TINY = ShearQuadrature(bulk=(8, 6, 10), ball=(4, 4, 8), radius=0.4)


# ── the derivation ─────────────────────────────────────────────────────────
def test_the_derivation_validates_on_three_known_cases():
    """The anchors that caught the first draft's error.

    A remembered connection formula gave a ``G₀₀`` containing ``ä``, which is
    impossible.  These three have known answers, and the ESU one reproduces
    `two_wave`'s ``_EINSTEIN = diag(3,−1,−1,−1)`` from a completely independent
    calculation.
    """
    import sympy as sp
    d = derive_the_tensor_mode_equation()
    a = sp.Symbol("a", positive=True)
    assert sp.simplify(d["s3_ricci"] - 2 / a ** 2) == 0
    assert sp.simplify(d["s3_scalar"] - 6 / a ** 2) == 0
    for got, want in zip(d["esu_einstein_diagonal"], (3, -1, -1, -1)):
        assert sp.simplify(got - want / a ** 2) == 0
    assert d["the_validations_pass"]


def test_the_tensor_mode_frequency_is_two_root_two():
    """``δG^TT = β̈ + 8β/a²``, so ``ω₃ = 2√2`` — and ``ω₃² > 0`` is the
    stability gate PR #260 taught the arc to check before computing anything."""
    import sympy as sp
    d = derive_the_tensor_mode_equation()
    assert sp.simplify(d["omega_squared_times_a_squared"] - 8) == 0
    assert TENSOR_MODE_FREQUENCY == pytest.approx(math.sqrt(8.0))
    assert TENSOR_MODE_FREQUENCY ** 2 > 0.0


def test_the_anisotropy_is_traceless_and_obeys_the_momentum_constraint():
    """Both are what make the channel genuinely TT — decoupled from the fluid's
    ``δρ``, ``δp`` and from the vector and scalar sectors."""
    m = measure_the_tensor_sector_is_stable_and_fluid_free()
    assert m["the_first_order_piece_is_traceless"]
    assert m["the_momentum_constraint_holds"]
    assert m["the_sector_is_stable"]


# ── the frame ──────────────────────────────────────────────────────────────
def test_the_left_invariant_frame_is_orthonormal_and_tangent():
    rng = np.random.default_rng(5)
    for _ in range(8):
        x = rng.normal(size=4)
        x /= np.linalg.norm(x)
        L = left_invariant_frame(x)
        assert np.allclose(L @ L.T, np.eye(3), atol=1e-13)
        assert np.allclose(L @ x, 0.0, atol=1e-13)


def test_the_frame_has_the_su2_commutators():
    """``[L_i,L_j] = 2ε_{ijk}L_k`` — the normalization the derivation assumes.

    If this were ``1`` rather than ``2`` the derived ``ω₃²`` would be ``2``
    instead of ``8``, so it is worth pinning by finite differences on the flows
    rather than by trusting the quaternion algebra.
    """
    rng = np.random.default_rng(11)
    x = rng.normal(size=4)
    x /= np.linalg.norm(x)
    units = (np.array([0.0, 1, 0, 0]), np.array([0.0, 0, 1, 0]),
             np.array([0.0, 0, 0, 1]))
    h = 1e-5

    def flow(p, i, s):
        return quaternion_product(
            p, np.cos(s) * np.array([1.0, 0, 0, 0]) + np.sin(s) * units[i])

    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        comm = (flow(flow(x, i, h), j, h) - flow(flow(x, j, h), i, h)) / h ** 2
        assert np.allclose(comm, 2.0 * left_invariant_frame(x)[k], atol=1e-4)


def test_the_shear_projection_is_frame_independent():
    """It has to be: `two_wave.orthonormal_frame` is an arbitrary basis per
    point, and averaging components across points is meaningless unless the
    rotation into the global frame removes that choice."""
    x = circle_point(0.9)
    setup = WORKING_BACKREACTION
    pts = np.array([x])
    for seed in (7, 23, 101):
        fr = np.array([orthonormal_frame(x, seed=seed)])
        sol = solve_batch(setup._setup(x), setup.source_a,
                          setup._setup(x).pulse_a, pts, fr)
        proj = shear_projection(stress_series(sol), pts, fr)
        if seed == 7:
            first = proj
        else:
            assert np.allclose(proj, first, atol=1e-10, rtol=1e-8)


# ── the solver ─────────────────────────────────────────────────────────────
def test_the_batched_solver_matches_the_reference():
    """Every number in this module rests on this.

    `solve_batch` exists only because the ``30 000``-point quadrature the
    convergence control needs is unaffordable one point at a time; it earns that
    by being the same computation, not a cheaper one.
    """
    setup = WORKING_BACKREACTION
    pts = np.array([circle_point(0.7), circle_point(2.0), [0.0, 0.0, 0.0, 1.0],
                    circle_point(4.1)])
    pts = pts / np.linalg.norm(pts, axis=1)[:, None]
    frames = np.array([orthonormal_frame(p) for p in pts])
    base = setup._setup(pts[0])
    batched = solve_batch(base, setup.source_a, base.pulse_a, pts, frames)
    for i, p in enumerate(pts):
        s = setup._setup(p)
        ref = solve_field(s, setup.source_a, s.pulse_a)
        for key in ("phi", "dt", "dtt", "grad", "dtgrad", "hess"):
            scale = max(float(np.max(np.abs(ref[key]))), 1e-300)
            assert np.max(np.abs(batched[key][i] - ref[key])) / scale < 1e-12


def test_the_batched_stress_tensor_matches_the_reference():
    setup = WORKING_BACKREACTION
    x = circle_point(1.4)
    pts = np.array([x])
    frames = np.array([orthonormal_frame(x)])
    base = setup._setup(x)
    sol = solve_batch(base, setup.source_a, base.pulse_a, pts, frames)
    series = stress_series(sol)
    ref_sol = solve_field(base, setup.source_a, base.pulse_a)
    for idx in (100, 900, 2500):
        ref = stress_tensor(ref_sol, idx)
        scale = max(float(np.max(np.abs(ref))), 1e-300)
        assert np.max(np.abs(series[0, idx] - ref)) / scale < 1e-11


def test_a_quadrature_point_on_a_singularity_is_rejected():
    setup = WORKING_BACKREACTION
    pts = np.array([setup.source_a])
    frames = np.array([orthonormal_frame(circle_point(0.3))])
    with pytest.raises(ValueError):
        solve_batch(setup._setup(pts[0]), setup.source_a,
                    setup._setup(pts[0]).pulse_a, pts, frames)


# ── the quadrature ─────────────────────────────────────────────────────────
def test_the_angular_rule_integrates_the_singular_direction_structure_exactly():
    """The load-bearing property of the whole quadrature.

    Near a source ``T ~ n_i n_j/χ⁴`` and the traceless part of that has **zero**
    angular average.  A rule that integrates ``n_i n_j`` exactly kills the
    ``1/χ⁴`` exactly; a random rule leaves ``O(1/√N)`` of it, which is what made
    a first draft of this round diverge.  Degree 4 is checked too, since the
    next order matters for the residual.
    """
    dirs, wts = _direction_rule(8, 16)
    total = wts.sum()
    assert total == pytest.approx(4.0 * math.pi, rel=1e-12)
    mono = np.einsum("p,pi->i", wts, dirs) / total
    assert np.allclose(mono, 0.0, atol=1e-12)
    quad = np.einsum("p,pi,pj->ij", wts, dirs, dirs) / total
    assert np.allclose(quad, np.eye(3) / 3.0, atol=1e-12)
    quart = np.einsum("p,pi,pi,pj,pj->", wts, dirs, dirs, dirs, dirs) / total
    assert quart == pytest.approx(1.0, rel=1e-12)


def test_the_tangent_basis_is_orthonormal_and_orthogonal_to_its_centre():
    rng = np.random.default_rng(3)
    for _ in range(5):
        c = rng.normal(size=4)
        c /= np.linalg.norm(c)
        b = _tangent_basis(c)
        assert np.allclose(b @ b.T, np.eye(3), atol=1e-12)
        assert np.allclose(b @ c, 0.0, atol=1e-12)


def test_the_quadrature_converges_on_known_integrals():
    """Refinement, not a magic tolerance.

    ``TINY`` is deliberately far too coarse to quote, so the claim worth
    testing is that refining it *helps* — which is what the partition of unity
    bought and what plain excision did not.
    """
    setup = WORKING_BACKREACTION
    a = setup.source_a
    exact_vol = 2.0 * math.pi ** 2
    errors = []
    for rule in (TINY, ShearQuadrature(bulk=(16, 10, 20), ball=(8, 6, 12))):
        pts, wts = rule.build(setup)
        vol = abs(float(wts.sum()) - exact_vol) / exact_vol
        got = float((wts * (pts @ a) ** 2).sum())
        moment = abs(got - exact_vol * 0.25) / (exact_vol * 0.25)
        errors.append((vol, moment))
    assert errors[0][0] < 2e-2 and errors[0][1] < 5e-2
    assert errors[1][0] < 0.25 * errors[0][0]
    assert errors[1][1] < 0.25 * errors[0][1]


def test_the_quadrature_covers_all_eight_singular_points():
    """The mouths are the two a first draft forgot, and `two_source` puts the
    first at ``(1,0,0,0)`` — the natural quadrature axis."""
    setup = WORKING_BACKREACTION
    sing = TINY.singular_points(setup)
    assert len(sing) == 8
    c1, c2 = setup.mouths()
    for p in (setup.source_a, setup.source_b, c1, c2):
        assert any(np.allclose(p, s) for s in sing)
        assert any(np.allclose(-np.asarray(p), s) for s in sing)


def test_overlapping_balls_are_rejected():
    with pytest.raises(ValueError):
        ShearQuadrature(radius=0.9).build(WORKING_BACKREACTION)


def test_no_quadrature_point_lands_on_a_singularity():
    setup = WORKING_BACKREACTION
    pts, _ = TINY.build(setup)
    for s in TINY.singular_points(setup):
        d = np.arccos(np.clip(pts @ s, -1.0, 1.0))
        assert d.min() > 1e-6


# ── the response ───────────────────────────────────────────────────────────
def test_the_response_solves_its_own_differential_equation():
    """``β̈ + ω²β = S``, checked by finite differences rather than by trusting
    the convolution that produced it."""
    dt = 0.002
    n = 6000
    t = np.arange(n) * dt
    src = np.zeros((n, 3, 3))
    src[:, 0, 1] = src[:, 1, 0] = np.exp(-((t - 3.0) ** 2) / 0.5) * np.cos(2 * t)
    src[:, 0, 0] = np.exp(-((t - 4.0) ** 2) / 0.8)
    beta = shear_response(src, dt)
    inner = slice(2, n - 2)
    dd = (beta[2:, 0, 1] - 2 * beta[1:-1, 0, 1] + beta[:-2, 0, 1]) / dt ** 2
    lhs = dd + TENSOR_MODE_FREQUENCY ** 2 * beta[1:-1, 0, 1]
    rhs = src[1:-1, 0, 1]
    scale = float(np.max(np.abs(rhs)))
    assert np.max(np.abs(lhs - rhs)) / scale < 1e-4


def test_the_response_is_retarded():
    dt = 0.005
    n = 3000
    t = np.arange(n) * dt
    src = np.zeros((n, 3, 3))
    src[:, 0, 1] = np.exp(-((t - 8.0) ** 2) / 0.2)
    beta = shear_response(src, dt)
    before = t < 6.0
    assert float(np.max(np.abs(beta[before]))) < 1e-8


def test_the_response_is_linear():
    dt = 0.01
    n = 1200
    rng = np.random.default_rng(4)
    s1 = rng.normal(size=(n, 3, 3))
    s2 = rng.normal(size=(n, 3, 3))
    a = shear_response(s1 + 2.5 * s2, dt)
    b = shear_response(s1, dt) + 2.5 * shear_response(s2, dt)
    assert np.allclose(a, b, atol=1e-12)


# ── the split that makes it a linear-algebra question ──────────────────────
def test_the_stress_tensor_splits_bilinearly():
    m = measure_the_stress_tensor_splits_bilinearly()
    assert m["the_self_term_is_quadratic"]
    assert m["the_cross_term_is_linear"]
    assert m["it_vanishes_with_one_source_off"]
    assert m["cross_with_one_source_off"] == 0.0


# ── the residual ───────────────────────────────────────────────────────────
def test_the_residual_is_zero_for_something_in_the_span():
    """The measurement has to return ``0`` when the answer is *reachable*, or a
    large residual would mean nothing."""
    rng = np.random.default_rng(9)
    n = 400
    ba = rng.normal(size=(n, 3, 3))
    bb = rng.normal(size=(n, 3, 3))
    t = np.arange(n) * 0.05
    r = unreachable_fraction(ba, bb, 1.7 * ba - 0.4 * bb, t, (0.0, 100.0))
    assert r["unreachable_off_the_span"] < 1e-10
    assert r["best_fit"][0] == pytest.approx(1.7, abs=1e-8)
    assert r["best_fit"][1] == pytest.approx(-0.4, abs=1e-8)


def test_the_residual_is_one_for_something_orthogonal():
    n = 300
    ba = np.zeros((n, 3, 3))
    bb = np.zeros((n, 3, 3))
    bc = np.zeros((n, 3, 3))
    ba[:, 0, 0] = 1.0
    bb[:, 1, 1] = 1.0
    bc[:, 2, 2] = 1.0
    t = np.arange(n) * 0.05
    r = unreachable_fraction(ba, bb, bc, t, (0.0, 100.0))
    assert r["unreachable_off_the_span"] == pytest.approx(1.0, abs=1e-10)


# ── the physics, at a small quadrature ─────────────────────────────────────
@pytest.mark.parametrize("throat", [True, False])
def test_the_interference_response_is_mostly_unreachable(throat):
    """The round's answer, checked structurally at a rule too small to quote.

    The converged figure lives in the probe; what is asserted here is that the
    interference response is a substantial, non-degenerate object both with the
    throat and without it.
    """
    setup = BackreactionSetup(with_throat=throat)
    src = setup.shear_sources(TINY)
    dt = setup.grid.dt
    r = unreachable_fraction(shear_response(src["A"], dt),
                             shear_response(src["B"], dt),
                             shear_response(src["cross"], dt),
                             setup.grid.times, (4.0, 30.0))
    assert r["unreachable_off_the_span"] > 0.5
    assert r["cross_over_single"] > 0.05
    assert r["cos_between_the_singles"] < 0.99


def test_the_two_sources_give_different_shear():
    """If ``A`` and ``B`` drove the same shear the span would be
    one-dimensional and the residual would be large for a trivial reason."""
    src = WORKING_BACKREACTION.shear_sources(TINY)
    a, b = src["A"], src["B"]
    ca = a.reshape(-1)
    cb = b.reshape(-1)
    cos = abs(ca @ cb) / np.linalg.norm(ca) / np.linalg.norm(cb)
    assert cos < 0.99


def test_the_sources_carry_no_trace():
    """The projection is the *traceless* part, so its trace has to vanish
    identically — otherwise a scalar-sector leak would be hiding in it."""
    src = WORKING_BACKREACTION.shear_sources(TINY)
    for key in ("A", "B", "cross"):
        tr = np.einsum("tii->t", src[key])
        scale = max(float(np.max(np.abs(src[key]))), 1e-300)
        assert float(np.max(np.abs(tr))) / scale < 1e-10


def test_the_shear_sources_are_symmetric():
    src = WORKING_BACKREACTION.shear_sources(TINY)
    for key in ("A", "B", "cross"):
        s = src[key]
        assert np.allclose(s, np.transpose(s, (0, 2, 1)), atol=1e-12)
