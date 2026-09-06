"""Round 9 — positive-counting reachability, pre-registered in
``docs/positive_counting_prereg.md`` (`a987342`) before the open questions."""

import math

import numpy as np

from geometrodynamics.bulk.positive_counting import (
    reduction_residual, W, correlation_from_phi, chsh_from_phi, quantum_phi,
    monomial_W_coefficients, G_in_shifted_basis, symmetry_conditions,
    solve_minimal_degree, no_global_polynomial, threshold_phi_reaches_four,
    marginals_are_automatic, closed_form_W, sector_sum_identity,
    downstream_numbers_unchanged,
    dependency_ledger, verdict, D_MIN, D_MAX)


# ── the oracle ──────────────────────────────────────────────────────────────

def test_the_reduction_holds_as_an_equality_of_laws():
    r = reduction_residual()
    assert r["reduction_holds"], r["worst_sorted_residual"]
    assert r["t_pair_sums_to_two"]


def test_the_known_weights_reproduce_rounds_five_and_six():
    """Phi = |D| is round 5's positive branch; Phi = D is round 6's oriented sum."""
    assert abs(correlation_from_phi(np.abs, 1.0) + 0.3984966504) < 1e-9
    assert abs(correlation_from_phi(lambda d: d, 1.0) + math.cos(1.0)) < 1e-9


def test_marginals_are_a_property_of_the_class_not_of_a_member():
    m = marginals_are_automatic()
    assert m["marginals_forced_half"], m["flip_residual"]


# ── Q1: positivity bounds nothing ───────────────────────────────────────────

def test_a_nonnegative_threshold_weight_reaches_the_algebraic_maximum():
    q = threshold_phi_reaches_four()
    assert q["reaches_algebraic_maximum"]
    for row in q["rows"]:
        assert abs(row["CHSH"] - 4.0) < 1e-9, row
        assert abs(row["E_pi_over_4"] + 1.0) < 1e-9
        assert abs(row["E_3pi_over_4"] - 1.0) < 1e-9


# ── Q2: the quantum law is attained ─────────────────────────────────────────

def test_the_explicit_weight_is_nonnegative_on_the_physical_range():
    grid = np.linspace(D_MIN, D_MAX, 200001)
    assert quantum_phi(grid).min() > -1e-9
    # and it is negative beyond the physical range, which is why the scope matters
    assert quantum_phi(np.array([6.0]))[0] < 0.0


def test_the_explicit_weight_reproduces_minus_cos_to_machine_precision():
    """Pre-registered rule 1: only an explicit Phi verified at independent
    angles discharges Q2."""
    for gamma in (0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
        got = correlation_from_phi(quantum_phi, gamma)
        assert abs(got + math.cos(gamma)) < 1e-12, (gamma, got)


def test_the_explicit_weight_gives_tsirelson():
    assert abs(chsh_from_phi(quantum_phi) - 2.0 * math.sqrt(2.0)) < 1e-12


# ── the polynomial structure ────────────────────────────────────────────────

def test_monomial_transform_matches_direct_integration():
    for n in range(0, 6):
        coeffs = monomial_W_coefficients(n)
        for t in (0.3, 0.9, 1.7):
            closed = sum(c * t ** k for k, c in coeffs.items()) * 2 * math.pi
            direct = W(t, lambda d, n=n: np.asarray(d) ** n)
            assert abs(closed - direct) / max(1.0, abs(closed)) < 1e-6, (n, t)


def test_G1_is_constant_so_the_signed_weight_is_the_obvious_solution():
    s = symmetry_conditions(6)
    assert s["G1_is_constant"]
    assert s["leading_coefficients_are_one"]


def test_degree_three_is_the_minimal_nonnegative_solution():
    m = solve_minimal_degree()
    assert abs(m["G2_odd_u"] - 1.0) < 1e-12
    assert abs(m["G3_odd_u"] - 5.0) < 1e-12
    assert m["combined_is_even"]
    assert m["degree_2_forces_a2_zero"]
    assert m["phi_nonnegative_on_range"]


def test_no_polynomial_solution_is_globally_nonnegative():
    g = no_global_polynomial()
    assert g["every_even_top_degree_is_forced_zero"]


# ── rule 5 and the report ───────────────────────────────────────────────────

def test_nothing_in_rounds_five_to_eight_moved():
    d = downstream_numbers_unchanged()
    assert d["nothing_downstream_moved"], d["checks"]
    assert d["worst_delta"] < 1e-6


def test_the_verdicts_are_the_pre_registered_labels():
    v = verdict()
    assert v["Q1_positivity_bound"] == "NO_BOUND_ALGEBRAIC_MAXIMUM_REACHED"
    assert v["Q2_quantum_law"] == "QUANTUM_LAW_ATTAINABLE_BY_POSITIVE_COUNTING"
    assert v["Q2_worst_error_vs_minus_cos"] < 1e-12


def test_the_ledger_keeps_the_class_restriction_and_the_open_item():
    led = {r["input"]: r["status"] for r in dependency_ledger()}
    assert led["class restriction: W_s = int Phi(D_s) dsigma, Phi >= 0"] == "chosen"
    assert led["which Phi the physics selects"] == "open"


def test_round_six_framing_is_withdrawn_in_the_write_up():
    import pathlib
    text = pathlib.Path("docs/closure_measurement_dependence.md").read_text()
    assert "withdrawn by round 9" in text


# ── corrections from the round-9 review ─────────────────────────────────────

def test_a_failed_witness_reports_non_resolution_not_a_no_go():
    """Pre-registered rule 2: a no-go needs a dual certificate. Exercise the
    failure path by monkeypatching the witness to something that misses."""
    import geometrodynamics.bulk.positive_counting as pc
    original = pc.quantum_phi
    try:
        pc.quantum_phi = staticmethod(lambda d: np.abs(np.asarray(d, float)))
        v = pc.verdict()
        assert v["Q2_quantum_law"] == "UNRESOLVED_NUMERICALLY"
        assert v["Q2_witness"] is None
        assert v["no_go_requires_a_dual_certificate"]
    finally:
        pc.quantum_phi = original


def test_the_closed_form_proves_the_symmetry_outright():
    """W(t) = (2pi/5) t (5 + 2t - t^2); the second factor is invariant under
    t -> 2 - t, which proves E = -cos gamma without angle-by-angle checking."""
    for t in (0.2, 0.7, 1.0, 1.3, 1.9):
        assert abs(closed_form_W(t) - W(t, quantum_phi)) < 1e-9
    for t in (0.3, 0.8, 1.4):
        assert abs((5 + 2*t - t*t) - (5 + 2*(2-t) - (2-t)**2)) < 1e-12
    # the radial form frozen in the preregistration is STRICTLY stronger:
    # it would require W(1)/2pi = 1, and the cubic gives 1.2
    assert abs(closed_form_W(1.0) / (2 * math.pi) - 1.2) < 1e-12


def test_the_step_correlation_is_a_limiting_statement():
    """A finite bump does not give E = -sgn(cos gamma) at every angle."""
    v = 1.0 + math.sqrt(2.0)
    phi = lambda d: (np.abs(np.asarray(d, float) - v) < 0.2).astype(float)
    off = correlation_from_phi(phi, math.acos(0.1))
    assert abs(off + 1.0) > 0.4, off            # nowhere near the step value
    # but the standard angles, which is all CHSH = 4 needs, are exact
    assert abs(correlation_from_phi(phi, math.pi/4) + 1.0) < 1e-12
    assert abs(correlation_from_phi(phi, 3*math.pi/4) - 1.0) < 1e-12


def test_attaining_the_singlet_law_does_not_close_the_causality_gate():
    s = sector_sum_identity()
    assert s["identity_holds"], s["identity_residual"]
    assert s["orthogonal_sum_is_constant"]
    assert s["hazard_survives"]
