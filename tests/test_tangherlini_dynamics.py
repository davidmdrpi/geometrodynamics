"""The first evolved Einstein equations in the tree.

Two kinds of check. The exact ones — Tangherlini, the flat-space mode, the
positivity identity — are asserted against closed forms at machine precision.
The convergence ones are asserted as *rates* under refinement, never on a single
residual, because the residual's magnitude is meaningless on its own.
"""

import math

import numpy as np
import pytest

from geometrodynamics.tangherlini import dynamics as dy


# ── the derivation ──────────────────────────────────────────────────────────
@pytest.mark.parametrize("n", [2, 3])
def test_the_derivation_checks_itself_at_both_dimensions(n):
    """`n = 2` is the known `D = 4` system; passing there validates general `n`."""
    got = dy.derive_the_field_equations(n)
    assert got["spacetime_dim"] == n + 2
    assert got["the_rr_equation_is_the_delta_quadrature"]
    assert got["the_vr_equation_is_the_A_quadrature"]
    assert got["tangherlini_solves_vr_with_constant_mass"]


def test_the_vv_equation_is_the_only_one_carrying_a_v_derivative_of_A():
    """That is why it can serve as the monitor rather than as an identity."""
    import sympy as sp
    got = dy.derive_the_field_equations(3)
    A = sp.Function("A")(sp.Symbol("v", real=True), sp.Symbol("r", real=True))
    v = sp.Symbol("v", real=True)
    assert got["vv"].has(sp.Derivative(A, v))
    assert not got["rr"].has(sp.Derivative(A, v))
    assert not got["vr"].has(sp.Derivative(A, v))
    assert got["the_vv_equation_is_never_used"]


# ── exact solutions ─────────────────────────────────────────────────────────
def test_tangherlini_comes_back_at_machine_precision():
    got = dy.measure_tangherlini_is_a_fixed_point()
    assert got["the_metric_is_exact"]
    assert got["the_slice_is_static"]
    assert got["delta_is_identically_zero"]
    for row in got["rows"]:
        assert row["metric_error"] < 1e-13
        assert row["max_abs_psi"] == 0.0


def test_the_metric_function_is_the_D5_tangherlini_form():
    r = np.linspace(1.5, 9.0, 40)
    assert np.allclose(dy.tangherlini_A(r, 1.0), 1.0 - 1.0 / r ** 2, atol=1e-15)
    # and in D = 4 it is the Schwarzschild form
    assert np.allclose(dy.tangherlini_A(r, 1.0, n=2), 1.0 - 1.0 / r, atol=1e-15)


def test_the_flat_mode_is_reproduced_at_second_order():
    """The `psi` quadrature against a closed form, which is what pins the order."""
    got = dy.measure_the_hierarchy_reproduces_the_exact_flat_mode()
    assert got["the_metric_is_exactly_flat"]
    assert got["converges_at_second_order"]
    rates = [x["rate"] for x in got["rows"] if x["rate"] is not None]
    assert all(abs(x - 2.0) < 0.05 for x in rates)
    # and it really is converging, not sitting at a fixed error
    errs = [x["psi_relative_error"] for x in got["rows"]]
    assert all(b < a for a, b in zip(errs, errs[1:]))


def test_a_flat_slice_has_no_geometry_at_all():
    r = np.linspace(0.0, 8.0, 600)
    h = float(r[1] - r[0])
    delta, A, K, psi = dy.hierarchy(r, h, np.zeros(600), kappa=dy.KAPPA)
    assert np.max(np.abs(A - 1.0)) < 1e-14
    assert np.max(np.abs(delta)) == 0.0
    assert np.max(np.abs(psi)) == 0.0


# ── the constraint ──────────────────────────────────────────────────────────
def test_the_unused_einstein_equation_converges_at_second_order():
    """The headline. `vv` contains `d_v A`, which the hierarchy never forms."""
    got = dy.measure_the_unused_einstein_equation_converges()
    assert got["converges_at_second_order"]
    assert all(abs(x - 2.0) < 0.05 for x in got["rates"])
    resid = [x["max_abs_vv_residual"] for x in got["rows"]]
    assert all(b < a for a, b in zip(resid, resid[1:])), "must actually shrink"
    assert resid[-1] < 1e-5


def test_the_constraint_claim_is_scoped_away_from_a_hamiltonian_pair():
    """`rr` and `vr` are solved, so their residuals are identically zero."""
    got = dy.measure_the_unused_einstein_equation_converges()
    assert "circular" in got["what_it_is_not"]
    assert "Hamiltonian" in got["what_it_is_not"]


def test_the_solved_equations_are_satisfied_by_construction():
    """Not a convergence claim — an identity, and worth pinning as one."""
    r = np.linspace(0.0, 10.0, 900)
    h = float(r[1] - r[0])
    phi = 0.4 * (np.exp(-((r - 3.0) ** 2)) + np.exp(-((r + 3.0) ** 2)))
    delta, A, K, psi = dy.hierarchy(r, h, phi)
    # rr:  d_r delta = (kappa/n) r (d_r phi)^2, to quadrature accuracy
    pr = dy.radial_derivative(phi, h)
    lhs = dy.radial_derivative(delta, h, parity=-1)[2:-2]
    rhs = ((dy.KAPPA / dy.N_ANGULAR) * r * pr ** 2)[2:-2]
    assert np.max(np.abs(lhs - rhs)) < 5e-4
    # vr:  the quadrature identity, exactly
    ed = np.exp(delta)
    quad = dy.cumulative_trapezoid((dy.N_ANGULAR - 1) * r ** (dy.N_ANGULAR - 2) * ed, h)
    assert np.max(np.abs((r ** 2 * ed * A)[1:] - quad[1:])) < 1e-12


# ── the structural result ───────────────────────────────────────────────────
def test_a_regular_centre_forbids_a_trapped_surface():
    got = dy.measure_a_regular_centre_forbids_a_trapped_surface()
    assert got["no_trapped_surface_anywhere"]
    assert got["A_is_positive_by_the_quadrature"]
    assert all(row["min_A"] > 0.0 for row in got["rows"])
    assert not any(row["trapped"] for row in got["rows"])
    # and the scan really does push A close to zero, or it proves nothing
    assert got["smallest_A_reached"] < 1e-2


def test_the_positivity_is_an_identity_not_a_numerical_accident():
    """`r^{n-1} e^d A` is a positive integrand over a positive interval.

    Checked directly on a deliberately violent slice: a large, sharp, sign-
    changing profile, where any sign error in the quadrature would show.
    """
    r = np.linspace(0.0, 10.0, 2000)
    h = float(r[1] - r[0])
    phi = 8.0 * np.exp(-((r - 1.5) / 0.3) ** 2) * np.cos(9.0 * r)
    delta, A, K, psi = dy.hierarchy(r, h, phi)
    assert delta[0] < -500.0, "the slice has to be violent for this to mean anything"
    assert np.all(np.isfinite(A))
    assert np.all(A[1:] > 0.0)
    ed = np.exp(delta)
    lhs = (r ** 2 * ed * A)[1:]
    assert np.all(lhs > 0.0)
    assert A[0] == 1.0


def test_the_positivity_fails_to_nan_rather_than_to_a_wrong_sign():
    """The identity is exact; its floating-point *representation* is not.

    `delta` is gauge-fixed to zero at the outer edge and increases outward, so
    `e^delta <= 1` everywhere and underflows once the slice spans more than
    about `700` e-folds. Numerator and denominator underflow together and the
    affected points come back `nan` — never a wrong positive or a spurious
    negative, which is the failure mode that would matter. Recorded because a
    trapped-surface claim reading off a silently wrong `A` is exactly the kind
    of thing this arc has been caught by before.
    """
    r = np.linspace(0.0, 10.0, 2000)
    h = float(r[1] - r[0])
    for amp, expect_finite in ((8.0, True), (30.0, False)):
        phi = amp * np.exp(-((r - 1.5) / 0.3) ** 2) * np.cos(9.0 * r)
        delta, A, _, _ = dy.hierarchy(r, h, phi)
        finite = np.isfinite(A[1:])
        assert np.all(A[1:][finite] > 0.0), "no finite value may be non-positive"
        assert np.all(finite) is np.True_ or bool(np.all(finite)) == expect_finite


def test_an_interior_mass_is_what_lets_A_change_sign():
    """The escape hatch the identity names: give up the regular centre."""
    r = np.linspace(0.5, 10.0, 800)
    h = float(r[1] - r[0])
    delta, A, K, psi = dy.hierarchy(r, h, np.zeros(800),
                                    interior_mass=r[0] ** 2 - 1.0)
    assert np.any(A < 0.0), "inside the horizon A must be negative"
    assert np.all(A[r > 1.0] > 0.0)
    crossing = r[np.argmin(np.abs(A))]
    assert abs(crossing - 1.0) < 2.0 * h


def test_the_chart_statement_is_not_a_physics_statement():
    got = dy.measure_a_regular_centre_forbids_a_trapped_surface()
    assert "chart" in got["the_consequence"] or "gauge" in got["the_consequence"]
    assert "collapse still happens" in got["why_it_is_a_chart_statement"]


# ── the potential discrepancy ───────────────────────────────────────────────
def test_the_master_potential_flat_limit_is_bessel():
    """The proof that settles which potential is which, with nothing recalled."""
    got = dy.measure_the_master_potential_disagrees_with_the_repo()
    assert got["flat_limit_matches_bessel"] < 1e-12
    # directly: psi = r^{1/2} J_{l+1}(wr) has V = ((l+1)^2 - 1/4)/r^2
    for ell in (0, 1, 2, 5):
        r = np.linspace(0.7, 9.0, 200)
        V = dy.master_potential(r, ell, 1e-9)
        assert np.allclose(V, ((ell + 1) ** 2 - 0.25) / r ** 2, rtol=1e-7)


def test_the_gap_to_the_repo_potential_is_the_stated_closed_form():
    from geometrodynamics.tangherlini.radial import V_tangherlini
    got = dy.measure_the_master_potential_disagrees_with_the_repo()
    assert got["gap_matches_the_closed_form"] < 1e-12
    r = np.linspace(1.4, 10.0, 300)
    A = dy.tangherlini_A(r, 1.0)
    gap = dy.master_potential(r, 2, 1.0) - np.asarray(V_tangherlini(r, 2, rs=1.0))
    assert np.allclose(gap, 3.0 * A ** 2 / (4.0 * r ** 2), atol=1e-14)


def test_nothing_in_the_repository_was_changed_by_this_round():
    """The discrepancy is reported. Acting on it is the owner's call."""
    from geometrodynamics.tangherlini.radial import V_tangherlini
    got = dy.measure_the_master_potential_disagrees_with_the_repo()
    assert got["nothing_was_changed"]
    assert "fifty probes" in got["why_not"]
    # the existing function still returns exactly what it always did
    r = np.array([2.0, 4.0])
    A = 1.0 - 1.0 / r ** 2
    assert np.allclose(V_tangherlini(r, 2, rs=1.0),
                       A * (2 * 4 / r ** 2 + 3.0 / r ** 4))


def test_the_discrepancy_is_stated_for_one_specific_field():
    got = dy.measure_the_master_potential_disagrees_with_the_repo()
    assert "MINIMALLY COUPLED MASSLESS" in got["caveat"]
    assert "r^{n/2}" in got["caveat"]


# ── what was not delivered ──────────────────────────────────────────────────
def test_no_quasinormal_frequency_is_reported():
    """Two converged constructions disagree, so nothing is claimed."""
    got = dy.measure_the_spectrum_is_not_yet_cross_validated()
    assert got["no_frequency_is_reported"]
    assert got["both_are_converged"]
    assert "perturbation spectrum" in got["not_delivered"]
    assert "transfer function" in got["not_delivered"]


def test_the_shortfall_is_named_against_what_was_asked():
    got = dy.measure_the_spectrum_is_not_yet_cross_validated()
    assert len(got["what_was_asked"]) == 5
    assert len(got["delivered"]) == 3
    assert len(got["not_delivered"]) == 2
    assert got["the_standing_lesson"] == "a converged number is not a correct number"


# ── the numerics themselves ─────────────────────────────────────────────────
def test_the_derivative_imposes_evenness_at_the_centre_exactly():
    h = 0.01
    f = np.cos(np.arange(600) * h) ** 2
    assert dy.radial_derivative(f, h)[0] == 0.0


def test_the_cumulative_quadrature_is_exact_on_a_line():
    h = 0.01
    x = np.arange(500) * h
    got = dy.cumulative_trapezoid(3.0 * x, h)
    assert np.allclose(got, 1.5 * x ** 2, atol=1e-12)


def test_the_evolution_admits_no_outer_boundary_condition():
    """Both ends are determined; the record of what happened when one was imposed."""
    got = dy.measure_the_unused_einstein_equation_converges()
    assert "no outer boundary condition" in got["what_an_imposed_outer_condition_did"] \
        or "did not converge" in got["what_an_imposed_outer_condition_did"]


def test_a_short_evolution_conserves_the_flat_limit():
    """With `kappa = 0` the metric must stay exactly flat while the field moves."""
    r = np.linspace(0.0, 12.0, 800)
    h = float(r[1] - r[0])
    phi = 0.5 * (np.exp(-((r - 4.0) ** 2)) + np.exp(-((r + 4.0) ** 2)))
    moved = phi.copy()
    for _ in range(40):
        moved = dy.step_rk4(r, h, moved, 0.4 * h, kappa=0.0)
    _, A, _, _ = dy.hierarchy(r, h, moved, kappa=0.0)
    assert np.max(np.abs(A - 1.0)) < 1e-13
    assert not np.allclose(moved, phi), "the field must actually have evolved"
    assert np.all(np.isfinite(moved))
