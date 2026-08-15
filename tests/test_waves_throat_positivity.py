"""
Tests for the positive sector of the point-throat boundary family.

The load-bearing claims:

* the operator is non-negative **iff** ``A ⪰ Γ(0)`` — one inequality in four
  parameters, checked against an actual negative-``λ`` scan;
* the same argument **counts**: ``#{eigenvalues < λ*}`` is the number of
  negative eigenvalues of ``A − Γ(λ*)``, for any ``λ*`` below the free ground
  state;
* the region is a **forward light cone** ``x₀ ≥ |x|`` in ``A − Γ(0) = x₀I + x·σ``,
  whose null boundary carries a **zero mode** and whose apex carries two;
* the growth rate turns on like ``√ε`` past the boundary, with the coefficient
  fixed by the eigenvalue slope rather than fitted;
* and PR #256's wedge is the ``x₂ = x₃ = 0`` slice — correct there, and wrong in
  general, which is measured rather than asserted.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.throat_operator import (
    MouthPair,
    negative_lambda_modes,
    stability_thresholds,
)
from geometrodynamics.waves.throat_positivity import (
    allowed_charge_basis,
    apex,
    boundary_pair_from_unitary,
    finite_boundary_from_unitary,
    gram_derivative,
    is_maximal_self_adjoint,
    is_non_negative_pair,
    measure_the_criterion_extends_to_the_boundary_strata,
    measure_the_monotonicity_is_a_gram_matrix,
    reduced_boundary_form,
    boundary_point,
    cone_coordinates,
    cone_fraction,
    count_modes_below,
    hermitian_from_parameters,
    inertia_below,
    is_non_negative,
    measure_the_boundary_of_the_cone_is_a_zero_mode,
    measure_the_exchange_symmetric_wedge_is_a_slice,
    measure_the_growth_rate_turns_on_with_a_square_root,
    measure_the_inertia_theorem_counts_modes_below_any_threshold,
    measure_the_positive_sector_is_a_shifted_psd_cone,
    measure_where_the_apex_sits_as_the_mouths_separate,
    positivity_defect,
    threshold_matrix,
    threshold_scaling,
    zero_mode,
)

D = 1.3


# ── the threshold matrix ────────────────────────────────────────────────────
def test_the_threshold_matrix_is_gamma_at_zero():
    g = threshold_matrix(D)
    th = stability_thresholds(D)
    assert g[0, 0] == pytest.approx(th["g_at_zero"])
    assert g[0, 1] == pytest.approx(th["G_d_at_zero"])
    assert g[0, 1] == pytest.approx(g[1, 0])
    assert np.abs(g.imag).max() == 0.0


def test_the_threshold_matrix_eigenvalues_are_the_channel_thresholds():
    th = stability_thresholds(D)
    ev = np.linalg.eigvalsh(threshold_matrix(D))
    assert ev[0] == pytest.approx(th["antisymmetric_threshold"], abs=1e-14)
    assert ev[1] == pytest.approx(th["symmetric_threshold"], abs=1e-14)


# ── the criterion ───────────────────────────────────────────────────────────
@pytest.mark.parametrize("a1,a2,beta,stable", [
    (0.05, 0.05, 0.03 + 0j, True),
    (0.2, 0.2, 0.15 + 0j, True),
    (0.2, -0.13, 0.15 + 0.07j, False),
    (-0.4, 0.07, -0.09 + 0.31j, False),
    (0.0, 0.0, 0.0 + 0j, False),
])
def test_the_criterion_agrees_with_the_scan(a1, a2, beta, stable):
    A = hermitian_from_parameters(a1, a2, beta)
    p = MouthPair(D, a1, a2, beta)
    assert is_non_negative(A, D) is stable
    assert (not negative_lambda_modes(p, lambda_min=-40000.0,
                                      n_grid=6000)) is stable


def test_the_criterion_matches_a_random_sweep():
    rng = np.random.default_rng(3)
    for _ in range(60):
        a1, a2 = rng.normal(0, 0.15, 2)
        b = complex(*rng.normal(0, 0.15, 2))
        A = hermitian_from_parameters(float(a1), float(a2), b)
        p = MouthPair(D, float(a1), float(a2), b)
        scan_ok = not negative_lambda_modes(p, lambda_min=-40000.0,
                                            n_grid=6000)
        assert is_non_negative(A, D) is scan_ok


def test_the_criterion_is_the_light_cone():
    rng = np.random.default_rng(9)
    for _ in range(80):
        a1, a2 = rng.normal(0, 0.2, 2)
        b = complex(*rng.normal(0, 0.2, 2))
        A = hermitian_from_parameters(float(a1), float(a2), b)
        c = cone_coordinates(A, D)
        assert c["inside_the_cone"] is is_non_negative(A, D)
        assert c["lightlike_defect"] == pytest.approx(
            c["x0"] - c["norm_x"], abs=1e-15)


def test_the_cone_coordinates_are_the_pauli_components():
    A = hermitian_from_parameters(0.2, -0.13, 0.15 + 0.07j)
    c = cone_coordinates(A, D)
    g = threshold_matrix(D)
    assert c["x0"] == pytest.approx(0.5 * (0.2 - 0.13) - g[0, 0].real)
    assert c["x"][0] == pytest.approx(0.15 - g[0, 1].real)
    assert c["x"][1] == pytest.approx(-0.07)
    assert c["x"][2] == pytest.approx(0.5 * (0.2 + 0.13))


def test_the_cone_is_convex_and_closed_under_positive_scaling():
    """It is a cone with apex Γ(0), so this is a real structural check."""
    g = threshold_matrix(D)
    rng = np.random.default_rng(17)
    stable = []
    while len(stable) < 6:
        a1, a2 = rng.normal(0, 0.3, 2)
        b = complex(*rng.normal(0, 0.3, 2))
        A = hermitian_from_parameters(float(a1), float(a2), b)
        if is_non_negative(A, D):
            stable.append(A)
    for i in range(0, len(stable), 2):
        p, q = stable[i], stable[i + 1]
        mid = g + 0.5 * ((p - g) + (q - g))
        assert is_non_negative(mid, D)
        assert is_non_negative(g + 3.7 * (p - g), D)


# ── the inertia theorem ─────────────────────────────────────────────────────
@pytest.mark.parametrize("lstar", [-2.0, 0.0, 0.5, 0.9])
def test_the_inertia_theorem_counts_modes_below_the_threshold(lstar):
    rng = np.random.default_rng(int(abs(lstar) * 100) + 5)
    for _ in range(25):
        a1, a2 = rng.normal(0, 0.15, 2)
        b = complex(*rng.normal(0, 0.15, 2))
        A = hermitian_from_parameters(float(a1), float(a2), b)
        p = MouthPair(D, float(a1), float(a2), b)
        assert inertia_below(A, D, lstar) == count_modes_below(p, lstar)


def test_gamma_grows_with_lambda_below_threshold():
    """The monotonicity the whole argument rests on."""
    for lam in (-100.0, -9.0, -1.0, -0.05, -0.001):
        h = abs(lam) * 1e-5
        dg = ((threshold_matrix(D, lam + h) - threshold_matrix(D, lam - h))
              / (2.0 * h))
        assert (np.linalg.eigvalsh(dg) > 0).all()


def test_far_below_threshold_the_krein_matrix_is_positive():
    """Both eigenvalues run to +∞ as λ → −∞, which is the other half."""
    A = hermitian_from_parameters(0.2, -0.13, 0.15 + 0.07j)
    for lam in (-100.0, -1000.0, -10000.0):
        ev = np.linalg.eigvalsh(A - threshold_matrix(D, lam))
        assert (ev > 0).all()
    small = np.linalg.eigvalsh(A - threshold_matrix(D, -100.0))[0]
    big = np.linalg.eigvalsh(A - threshold_matrix(D, -10000.0))[0]
    assert big > small


def test_two_negative_eigenvalues_means_two_growing_modes():
    A = hermitian_from_parameters(-0.6, -0.5, 0.0 + 0j)
    assert positivity_defect(A, D)["n_negative"] == 2
    p = MouthPair(D, -0.6, -0.5, 0.0)
    assert len(negative_lambda_modes(p, lambda_min=-40000.0,
                                     n_grid=6000)) == 2


# ── the boundary is a zero mode ─────────────────────────────────────────────
@pytest.mark.parametrize("direction", [(1, 0, 0), (0, 1, 0), (0, 0, 1),
                                       (0.6, -0.5, 0.62)])
@pytest.mark.parametrize("x0", [0.02, 0.4])
def test_a_boundary_point_carries_a_zero_mode(direction, x0):
    A = boundary_point(D, direction, x0)
    zm = zero_mode(A, D)
    assert zm["is_a_zero_mode"]
    assert zm["degeneracy"] == 1
    p = MouthPair(D, float(A[0, 0].real), float(A[1, 1].real),
                  complex(A[0, 1]))
    assert abs(p.secular(0.0)) < 1e-12
    assert not [r for r in negative_lambda_modes(p) if r["lmbda"] < -1e-10]


def test_the_apex_carries_two_zero_modes():
    ap = apex(D)
    zm = zero_mode(ap["Gamma_0"], D)
    assert zm["degeneracy"] == 2
    assert abs(zm["eigenvalue"]) < 1e-15


def test_a_boundary_point_is_lightlike():
    for x0 in (0.05, 0.3):
        c = cone_coordinates(boundary_point(D, (0.3, 0.9, -0.2), x0), D)
        assert c["lightlike_defect"] == pytest.approx(0.0, abs=1e-15)
        assert c["x0"] == pytest.approx(x0, abs=1e-15)


# ── the threshold scaling ───────────────────────────────────────────────────
def test_the_growth_rate_rises_like_a_square_root():
    r = threshold_scaling(D)
    assert r["asymptotic_exponent"] == pytest.approx(0.5, abs=0.01)
    good = [x for x in r["rows"] if x["lmbda"] is not None]
    assert len(good) >= 4
    assert good[-1]["lambda_over_epsilon"] == pytest.approx(
        r["predicted_from_the_eigenvalue_slope"], rel=0.02)


def test_the_growth_rate_is_continuous_at_the_boundary():
    r = threshold_scaling(D)
    sigmas = [x["sigma"] for x in r["rows"] if x["sigma"] is not None]
    assert sigmas == sorted(sigmas, reverse=True)
    assert sigmas[-1] < 1e-2


# ── the wedge is a slice ────────────────────────────────────────────────────
def test_the_exchange_symmetric_wedge_is_the_slice():
    th = stability_thresholds(D)
    for a in np.linspace(-0.1, 0.15, 9):
        for b in np.linspace(-0.15, 0.15, 11):
            A = hermitian_from_parameters(float(a), float(a), float(b))
            wedge = ((a + b) >= th["symmetric_threshold"]
                     and (a - b) >= th["antisymmetric_threshold"])
            assert is_non_negative(A, D) is bool(wedge)


def test_an_imaginary_beta_can_leave_the_cone_without_moving_the_wedge():
    """The dimension the wedge could not see."""
    a, b = 0.05, 0.03
    assert is_non_negative(hermitian_from_parameters(a, a, b + 0j), D)
    assert not is_non_negative(
        hermitian_from_parameters(a, a, complex(b, 0.2)), D)


def test_unequal_mouths_can_leave_the_cone_without_moving_the_wedge():
    """And the other one."""
    a, b = 0.05, 0.03
    assert is_non_negative(hermitian_from_parameters(a, a, b), D)
    assert not is_non_negative(
        hermitian_from_parameters(a + 0.25, a - 0.25, b), D)


# ── the apex ────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("sep", [0.2, 0.8, 1.3, 2.0, 3.0])
def test_the_apex_trace_does_not_depend_on_the_mouth_separation(sep):
    assert apex(sep)["trace"] == pytest.approx(
        -1.0 / (2.0 * math.pi ** 2), abs=1e-15)


@pytest.mark.parametrize("sep", [0.2, 0.8, 1.3, 2.0, 3.0])
def test_away_from_the_antipode_a_zero_boundary_matrix_never_works(sep):
    ap = apex(sep)
    assert ap["indefinite"]
    assert not ap["zero_matrix_is_stable"]
    assert not is_non_negative(np.zeros((2, 2), dtype=complex), sep)


def test_the_symmetric_threshold_closes_as_the_mouths_go_antipodal():
    near = apex(1.0)["eigenvalues"][1]
    far = apex(3.05)["eigenvalues"][1]
    assert 0.0 < far < near / 100.0


# ── the antipodal endpoint, which is a different statement ──────────────────
def test_the_green_function_is_finite_at_antipodal_mouth_separation():
    """``G_d`` has a removable singularity at ``d = π``, not a pole."""
    g = threshold_matrix(math.pi)
    assert g[0, 1] == pytest.approx(1.0 / (4.0 * math.pi ** 2), abs=1e-15)
    assert g[0, 0] == pytest.approx(-1.0 / (4.0 * math.pi ** 2), abs=1e-15)


def test_the_antipodal_apex_is_negative_semidefinite_not_indefinite():
    ap = apex(math.pi)
    assert ap["eigenvalues"][0] == pytest.approx(-1.0 / (2.0 * math.pi ** 2),
                                                 abs=1e-15)
    assert ap["eigenvalues"][1] == pytest.approx(0.0, abs=1e-15)
    assert ap["negative_semidefinite"]
    assert not ap["indefinite"]


def test_at_the_antipode_a_zero_boundary_matrix_is_marginally_stable():
    """The correction: ``A = 0`` is on the cone's boundary there, not outside."""
    zero = np.zeros((2, 2), dtype=complex)
    assert is_non_negative(zero, math.pi)
    assert positivity_defect(zero, math.pi)["min_eigenvalue"] == (
        pytest.approx(0.0, abs=1e-15))
    zm = zero_mode(zero, math.pi)
    assert zm["is_a_zero_mode"]
    q = np.abs(zm["q"])
    assert q[0] == pytest.approx(q[1], abs=1e-9)      # the symmetric channel


def test_the_antipodal_limit_is_continuous_from_below():
    """``G₀(d)`` falls to ``1/(4π²)`` from above, and the symmetric threshold
    ``g₀ + G₀`` closes to zero — no jump at the endpoint."""
    gpi = threshold_matrix(math.pi)[0, 1]
    gaps = []
    for d in (3.0, 3.10, 3.14, 3.1415):
        gaps.append(threshold_matrix(d)[0, 1] - gpi)
    assert all(g > 0 for g in gaps)
    assert gaps == sorted(gaps, reverse=True)
    assert apex(math.pi - 1e-9)["eigenvalues"][1] == pytest.approx(
        0.0, abs=1e-12)


# ── the monotonicity is a Gram matrix ───────────────────────────────────────
@pytest.mark.parametrize("sep", [1.3, 2.6, math.pi])
@pytest.mark.parametrize("lam", [-9.0, -1.0, 0.0, 0.5])
def test_the_gram_sum_is_the_closed_form_derivative(sep, lam):
    h = max(abs(lam), 1.0) * 1e-6
    fd = (threshold_matrix(sep, lam + h) - threshold_matrix(sep, lam - h)) \
        / (2.0 * h)
    gm = gram_derivative(lam, sep)
    assert np.abs(fd - gm).max() < 1e-9
    assert (np.linalg.eigvalsh(gm) > 0).all()


def test_the_gram_matrix_is_positive_definite_for_distinct_mouths():
    """PSD for free; PD because the two resolvent vectors are independent."""
    for sep in (0.3, 1.3, 2.9, math.pi):
        for lam in (-50.0, -1.0, 0.0, 0.8):
            assert (np.linalg.eigvalsh(gram_derivative(lam, sep)) > 0).all()


# ── beyond the finite-A chart ───────────────────────────────────────────────
def _haar(rng):
    x = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, r = np.linalg.qr(x)
    return q @ np.diag(np.diag(r) / np.abs(np.diag(r)))


def test_every_unitary_gives_a_maximal_self_adjoint_pair():
    rng = np.random.default_rng(31)
    for _ in range(30):
        b, c = boundary_pair_from_unitary(_haar(rng))
        info = is_maximal_self_adjoint(b, c)
        assert info["maximal"] and info["self_adjoint"]


def test_the_chart_covers_exactly_the_unitaries_without_eigenvalue_one():
    rng = np.random.default_rng(5)
    v = _haar(rng)
    u_bad = v @ np.diag([1.0, np.exp(1j * 1.1)]) @ v.conjugate().T
    with pytest.raises(ValueError):
        finite_boundary_from_unitary(u_bad)
    u_ok = v @ np.diag([np.exp(1j * 2.2), np.exp(1j * 1.1)]) @ v.conjugate().T
    a = finite_boundary_from_unitary(u_ok)
    assert np.abs(a - a.conjugate().T).max() < 1e-10


def test_the_general_criterion_reduces_to_the_cone_on_the_chart():
    rng = np.random.default_rng(77)
    for _ in range(40):
        u = _haar(rng)
        b, c = boundary_pair_from_unitary(u)
        if abs(np.linalg.det(b)) < 1e-6:
            continue
        a = finite_boundary_from_unitary(u)
        gen = is_non_negative_pair(b, c, D)
        assert gen["k"] == 2
        assert gen["non_negative"] is is_non_negative(a, D)


def test_a_dirichlet_stratum_drops_one_direction():
    rng = np.random.default_rng(13)
    v = _haar(rng)
    u = v @ np.diag([1.0, np.exp(1j * 2.0)]) @ v.conjugate().T
    b, c = boundary_pair_from_unitary(u)
    info = allowed_charge_basis(b, c)
    assert info["k"] == 1 and info["n_dirichlet"] == 1
    red = reduced_boundary_form(b, c)
    assert red["hermitian_defect"] < 1e-9
    assert red["row_space_defect"] < 1e-9


def test_the_free_stratum_has_no_mouth_active_spectrum():
    b = np.zeros((2, 2), dtype=complex)
    c = np.eye(2, dtype=complex)
    assert is_maximal_self_adjoint(b, c)["maximal"]
    got = is_non_negative_pair(b, c, D)
    assert got["k"] == 0 and got["non_negative"]
    assert "free" in got["stratum"]


def test_the_dirichlet_criterion_matches_a_root_scan():
    from geometrodynamics.waves.throat_operator import gamma_at
    rng = np.random.default_rng(101)
    for _ in range(12):
        th = float(rng.uniform(0.2, 6.0))
        v = _haar(rng)
        u = v @ np.diag([1.0, np.exp(1j * th)]) @ v.conjugate().T
        b, c = boundary_pair_from_unitary(u)
        lams = -np.geomspace(1e-8, 4000.0, 4000)[::-1]
        vals = [float(np.linalg.det(c - b @ gamma_at(float(x), D)).real)
                for x in lams]
        roots = sum(1 for i in range(len(vals) - 1)
                    if vals[i] * vals[i + 1] < 0.0)
        assert is_non_negative_pair(b, c, D)["non_negative"] is (roots == 0)


# ── the box fraction ────────────────────────────────────────────────────────
def test_the_cone_fraction_is_reported_with_its_box():
    r = cone_fraction(D, half_width=0.2, n_draws=800)
    assert 0.0 < r["fraction"] < 1.0
    assert r["half_width"] == 0.2
    assert "box-dependent" in r["caveat"]


# ── the measurements ────────────────────────────────────────────────────────
def test_measure_the_positive_sector_is_a_shifted_psd_cone():
    r = measure_the_positive_sector_is_a_shifted_psd_cone()
    assert r["the_criterion_is_exact"]
    assert r["n_mismatches"] == 0
    assert r["both_verdicts_occur"]
    assert r["the_light_cone_form_agrees"]
    assert r["n_with_complex_beta"] == r["n_draws"]


def test_measure_the_inertia_theorem_counts_modes_below_any_threshold():
    r = measure_the_inertia_theorem_counts_modes_below_any_threshold()
    assert r["the_inertia_theorem_holds"]
    assert r["n_mismatches"] == 0
    assert r["n_tested"] == 160
    assert r["d_gamma_d_lambda_is_positive_definite"]


def test_measure_the_boundary_of_the_cone_is_a_zero_mode():
    r = measure_the_boundary_of_the_cone_is_a_zero_mode()
    assert r["every_boundary_point_has_a_zero_mode"]
    assert r["the_secular_function_vanishes_there"]
    assert r["the_marginal_mode_sits_at_lambda_zero"]
    assert r["no_boundary_point_is_strictly_unstable"]
    assert r["the_apex_carries_two_zero_modes"]


def test_measure_the_growth_rate_turns_on_with_a_square_root():
    r = measure_the_growth_rate_turns_on_with_a_square_root()
    assert r["exponent_is_one_half"]
    assert r["lambda_is_linear_in_epsilon"]
    assert r["the_coefficient_matches_the_eigenvalue_slope"]


def test_measure_the_exchange_symmetric_wedge_is_a_slice():
    r = measure_the_exchange_symmetric_wedge_is_a_slice()
    assert r["the_wedge_is_exactly_the_slice"]
    assert r["the_slice_really_is_x2_equals_x3_equals_zero"]
    assert r["the_slice_rule_is_not_enough_in_general"]
    assert r["n_the_wedge_rule_gets_wrong"] > 0


def test_measure_where_the_apex_sits_as_the_mouths_separate():
    r = measure_where_the_apex_sits_as_the_mouths_separate()
    assert r["trace_is_separation_independent"]
    assert r["trace_matches_minus_one_over_two_pi_squared"]
    assert r["the_apex_is_indefinite_away_from_the_antipode"]
    assert r["the_zero_matrix_is_unstable_away_from_the_antipode"]
    assert r["eigenvalues_are_the_channel_thresholds"]
    assert r["the_apex_is_negative_semidefinite_at_the_antipode"]
    assert r["the_antipodal_endpoint_is_marginal"]
    assert r["at_the_antipode_A_zero_sits_on_the_boundary"]
    assert r["the_marginal_channel_is_symmetric"] == "symmetric"


def test_measure_the_monotonicity_is_a_gram_matrix():
    r = measure_the_monotonicity_is_a_gram_matrix()
    assert r["the_gram_sum_is_the_closed_form_derivative"]
    assert r["positive_definite_everywhere"]
    assert r["including_at_the_antipode"]
    assert "Gram matrix" in r["the_identity"]


def test_measure_the_criterion_extends_to_the_boundary_strata():
    r = measure_the_criterion_extends_to_the_boundary_strata()
    assert r["the_general_form_agrees_with_the_cone_on_the_chart"]
    assert r["the_general_form_agrees_with_the_scan_on_the_strata"]
    assert r["every_stratum_has_one_dirichlet_direction"]
    assert r["the_reduction_is_legitimate"]
    assert r["the_free_stratum_is_non_negative"]
    assert r["chart_mismatches"] == 0 and r["stratum_mismatches"] == 0


def test_the_module_says_what_is_still_put_in():
    from geometrodynamics.waves import throat_positivity
    doc = throat_positivity.__doc__ or ""
    assert "point-supported" in doc
    assert "no backreaction" in doc.lower()
    assert "chosen, not derived" in doc
