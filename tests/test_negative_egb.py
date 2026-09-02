"""Does negative-coupling EGB actually work? Checked against
`docs/negative_egb_prereg.md`, frozen before the module.

Five kinds. The *exterior* ones pin `R_kk`, `H_kk` and the `χ`-independence
that makes the outside constraint a single number. The *search* one hunts for a
surviving interval of `α_GB` rather than asserting there is none — a surviving
interval would reopen the branch. The *mechanism* ones pin that the two bounds
meet by continuity of one bracket, not by coincidence. The *cost* ones pin the
vacuum-form exterior — and, correcting this round's first draft, that the throat
matter there is **not** exotic. And the *principal-symbol* ones carry the actual
closure: the TT kinetic coefficient, derived on this product background rather
than recalled from the maximally symmetric formula.
"""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import finite_mouth as fm
from geometrodynamics.bulk import gauss_bonnet as gb
from geometrodynamics.bulk import negative_egb as ne


# ── the exterior ────────────────────────────────────────────────────────────

def test_the_exterior_null_ricci_is_positive_unlike_the_throat():
    """`+3/R²` outside against `−3f″/f` inside: the sign difference is the
    whole mechanism."""
    for R in (1.0, 1.3, 2.5):
        assert ne.exterior_ricci_null(R) == pytest.approx(3.0 / R ** 2, rel=1e-14)
        assert ne.exterior_ricci_null(R) > 0.0
    g = fm.geometry()
    assert float(gb.ricci_null_contraction(np.array([0.0]), g.neck_radius)[0]) < 0.0


def test_the_exterior_lanczos_null_is_twelve_over_R_to_the_fourth():
    for R in (1.0, 1.3, 2.5):
        assert ne.exterior_lanczos_null(R) == pytest.approx(12.0 / R ** 4, rel=1e-14)


def test_the_exterior_ratio_is_chi_independent():
    """`H_kk/R_kk = 4μ/f⁴ = 4/R²` for every `χ`, computed from the exterior's
    own `μ` and `f` rather than asserted. That constancy is what makes the
    outside constraint one number instead of a profile."""
    for R in (1.0, 2.5):
        chi = np.linspace(0.05, math.pi - 0.05, 101)
        ratio = ne.exterior_ratio(chi, R)
        assert np.ptp(ratio) < 1e-12
        assert float(ratio[0]) == pytest.approx(4.0 / R ** 2, rel=1e-12)


def test_the_exterior_null_stress_and_its_threshold():
    for R in (1.0, 2.5):
        for coupling in (0.0, -0.1, -0.25 * R ** 2, -R ** 2):
            assert ne.exterior_matter_null_stress(coupling, R) == pytest.approx(
                3.0 * (R ** 2 + 4.0 * coupling) / R ** 4, rel=1e-13)
        assert ne.exterior_threshold(R) == pytest.approx(-0.25 * R ** 2, rel=1e-14)
        assert ne.exterior_matter_null_stress(
            ne.exterior_threshold(R), R) == pytest.approx(0.0, abs=1e-13)


def test_the_two_bounds_are_the_same_number_pointing_opposite_ways():
    result = ne.measure_the_exterior_constrains_alpha_oppositely()
    assert result["the_two_bounds_coincide"]
    assert result["exterior_needs_alpha_at_least"] == pytest.approx(
        result["throat_needs_alpha_at_most"], rel=1e-14)


# ── the search ──────────────────────────────────────────────────────────────

def test_no_open_set_of_couplings_satisfies_both_regions():
    """The falsifier: a surviving interval would reopen the branch."""
    result = ne.measure_no_coupling_satisfies_both()
    assert result["the_surviving_set_is_measure_zero"]
    assert result["surviving_width"] < 1e-6
    assert result["survivors_are_at_the_critical_value"]
    assert result["samples"] >= 1000


def test_the_scan_actually_covers_both_signs_of_the_coupling():
    """A scan that never went positive would beg the question."""
    result = ne.measure_no_coupling_satisfies_both()
    low, high = result["scanned_range"]
    assert low < result["critical_coupling"] < 0.0 < high
    assert result["throat_satisfied_count"] > 0
    assert result["exterior_satisfied_count"] > 0


def test_the_throat_and_exterior_admissible_sets_are_complementary():
    """Directly: at any coupling one side or the other fails, except one."""
    R, a = 1.0, 0.3
    g = fm.geometry(R, a)
    s = np.linspace(-g.half_length, g.half_length, 201)
    for coupling in (-1.0, -0.5, -0.3, -0.2, -0.1, 0.0, 0.5):
        throat = bool(np.all(gb.matter_null_stress(s, g.neck_radius, coupling)
                             >= -1e-9))
        exterior = ne.exterior_matter_null_stress(coupling, R) >= -1e-9
        assert not (throat and exterior), \
            f"alpha={coupling} would reopen the branch"


# ── the mechanism ───────────────────────────────────────────────────────────

def test_the_ratio_is_continuous_across_the_seam():
    """`μ/f⁴ = 1/R²` from both sides — the same Misner–Sharp continuity as P2."""
    result = ne.measure_the_bracket_is_continuous_at_the_seam()
    assert result["ratio_is_continuous"]
    assert result["ratio_jump"] < 1e-12
    assert result["ratio_inside_at_the_seam"] == pytest.approx(
        result["common_value"], rel=1e-12)


def test_the_shared_bracket_vanishes_at_the_seam_at_the_critical_coupling():
    """It must be `≤ 0` from inside and `≥ 0` from outside while being
    continuous, so it is exactly `0`."""
    result = ne.measure_the_bracket_is_continuous_at_the_seam()
    assert result["bracket_vanishes_at_the_seam"]
    assert abs(result["bracket_inside"]) < 1e-12
    assert abs(result["bracket_outside"]) < 1e-12


def test_the_bracket_has_the_required_signs_on_each_side():
    """Below the critical coupling the throat is satisfied and the exterior is
    not; above, the reverse."""
    R, a = 1.0, 0.3
    g = fm.geometry(R, a)
    inside = float(gb.gauss_bonnet_ratio(np.array([g.half_length]),
                                         g.neck_radius)[0])
    outside = float(ne.exterior_ratio(np.array([a]), R)[0])
    critical = ne.exterior_threshold(R)
    assert float(ne.bracket(1.2 * critical, inside)) < 0.0     # throat happy
    assert float(ne.bracket(1.2 * critical, outside)) < 0.0    # exterior not
    assert float(ne.bracket(0.8 * critical, inside)) > 0.0     # throat not
    assert float(ne.bracket(0.8 * critical, outside)) > 0.0    # exterior happy


# ── the cost ────────────────────────────────────────────────────────────────

def test_gauss_bonnet_cannot_move_the_exterior_pressure():
    """`H^i_j = 0` on a maximally symmetric slice, so the whole correction
    lands in `ρ`."""
    for R in (1.0, 2.5):
        base = ne.exterior_pressure(0.0, R)
        for coupling in (-2.0, -0.25 * R ** 2, 0.0, 3.0):
            assert ne.exterior_pressure(coupling, R) == pytest.approx(
                base, rel=1e-15)
            assert base == pytest.approx(-3.0 / R ** 2, rel=1e-14)
    assert ne.measure_the_critical_exterior_is_empty()[
        "pressure_is_coupling_independent"]


def test_the_critical_exterior_is_exactly_vacuum_energy():
    """`w = −1` exactly, at half the Einstein-gravity density."""
    result = ne.measure_the_critical_exterior_is_empty()
    assert result["it_is_exactly_vacuum_energy"]
    assert result["critical_equation_of_state"] == pytest.approx(-1.0, abs=1e-12)
    assert result["density_is_halved"]
    for R in (1.0, 2.5):
        critical = ne.exterior_threshold(R)
        assert ne.exterior_density(critical, R) == pytest.approx(
            3.0 / R ** 2, rel=1e-13)


# ── the throat matter at criticality is NOT exotic (post-review) ────────────

def test_the_throat_matter_satisfies_nec_wec_and_dec_at_criticality():
    """Corrects this round's "exotic matter merely relocated". With
    `A = b²/f⁴`, `q = R²b²/f⁴ ≥ 1`: `ρ = 3Aq`, `p_s = −3A`, `p_Ω = A`."""
    result = ne.measure_the_throat_matter_is_not_exotic()
    assert result["nec_holds"] and result["wec_holds"] and result["dec_holds"]
    assert result["min_nec_radial"] >= -1e-9
    assert result["min_nec_angular"] > 0.0
    assert result["saturated_at_the_mouth"]
    assert result["q_at_the_mouth"] == pytest.approx(1.0, abs=1e-9)
    assert result["q_at_the_neck"] == pytest.approx(
        result["q_neck_closed_form"], rel=1e-9)


def test_the_module_records_the_withdrawn_relocation_claim():
    note = ne.measure_the_throat_matter_is_not_exotic()["what_this_corrects"]
    assert "merely relocated" in note
    assert "It is not" in note


# ── the classical linearised operator, which is what closes the branch ──────
#
# This is a classical PDE question throughout: the principal symbol of the
# linearised field equations and whether the Cauchy problem stays an evolution
# problem. Nothing here quantises the metric.

def test_the_kinetic_coefficient_is_one_plus_four_alpha_over_R_squared():
    """Derived by linearising on this product background, not recalled: the
    textbook `1 + 2α(D−3)(D−4)K` is for a maximally symmetric *spacetime*."""
    result = ne.measure_the_principal_symbol_degenerates()
    assert result["law_holds"]
    for row in result["rows"]:
        assert row["temporal_coefficient"] == pytest.approx(
            row["predicted_kinetic"], abs=1e-5)


def test_the_principal_symbol_is_a_measured_quadratic_form():
    """The symbol's shape is sampled off-axis rather than assumed. Without this
    the `C_t ω² + C_s κ²` reading would be an ansatz, not a measurement."""
    result = ne.measure_the_principal_symbol_degenerates()
    assert result["symbol_is_quadratic_and_isotropic"]
    assert result["worst_direction_error"] < 1e-4
    mixed = [r for r in result["direction_rows"] if "mixed" in r["direction"]]
    assert mixed, "a direction with both time and space components is required"
    for row in mixed:
        assert row["measured"] == pytest.approx(row["predicted"], abs=1e-5)


def test_the_leading_time_derivative_vanishes_at_the_critical_coupling():
    """The closure. Same coupling the NEC pins."""
    result = ne.measure_the_principal_symbol_degenerates()
    assert result["kinetic_vanishes_at_criticality"]
    assert abs(result["rows"][-1]["temporal_coefficient"]) < 1e-5
    for R in (1.0, 2.5):
        assert ne.tensor_kinetic_coefficient(
            ne.exterior_threshold(R), R) == pytest.approx(0.0, abs=1e-14)


def test_the_symbol_drops_from_degree_two_to_degree_zero_in_omega():
    """The precise classical statement: the Cauchy problem stops being an
    evolution problem in this sector."""
    result = ne.measure_the_principal_symbol_degenerates()
    assert result["degree_in_omega_drops_at_criticality"]
    assert result["rows"][0]["degree_in_omega"] == 2
    assert result["rows"][-1]["degree_in_omega"] == 0
    assert result["rows"][-1]["hyperbolic"] is False


def test_the_spatial_coefficient_is_untouched_so_it_is_a_degeneration():
    """`κ²` survives while `ω²` dies, and the lower-order term stays finite —
    the principal part degenerates rather than the whole equation vanishing."""
    result = ne.measure_the_principal_symbol_degenerates()
    assert result["spatial_coefficient_is_coupling_independent"]
    assert result["lower_order_term_stays_finite"]
    assert abs(result["rows"][-1]["spatial_coefficient"]) > 0.1


def test_the_open_interval_is_hyperbolic_but_acausal_not_ill_posed():
    """The distinction an earlier draft ran together. On `−R²/4 < α < 0` the
    operator is still hyperbolic; what fails is that the characteristic cone
    lies outside the metric null cone."""
    result = ne.measure_the_principal_symbol_degenerates()
    assert result["hyperbolic_on_the_open_interval"]
    assert result["cone_leaves_the_null_cone_before_criticality"]
    for row in result["rows"][1:-1]:
        assert row["hyperbolic"], "still well posed away from the endpoint"
        assert row["speed_squared"] > 1.0
    for fraction in (0.2, 0.5, 0.8):
        coupling = fraction * ne.exterior_threshold(1.0)
        assert ne.characteristic_speed_squared(coupling, 1.0) > 1.0
    assert ne.characteristic_speed_squared(0.0, 1.0) == pytest.approx(
        1.0, rel=1e-14)
    assert math.isinf(
        ne.characteristic_speed_squared(ne.exterior_threshold(1.0), 1.0))


def test_the_closed_form_symbol_agrees_with_the_measured_coefficients():
    """`principal_symbol` must reproduce the fitted numbers, and its zero must
    be the characteristic speed."""
    R = 1.0
    for fraction in (0.0, 0.5, 0.8):
        coupling = fraction * ne.exterior_threshold(R)
        speed = ne.characteristic_speed_squared(coupling, R)
        # A characteristic covector: omega^2 = c^2 kappa^2 makes P vanish.
        assert ne.principal_symbol(
            coupling, math.sqrt(speed), 1.0, R) == pytest.approx(0.0, abs=1e-12)
    critical = ne.exterior_threshold(R)
    # At criticality no finite omega solves P = 0 with kappa != 0.
    assert ne.principal_symbol(critical, 1e6, 1.0, R) == pytest.approx(0.5)


def test_the_module_retracts_the_quantum_wording():
    """BAM does not quantise the metric; the earlier 'graviton' framing is
    withdrawn and the record of the withdrawal is kept."""
    result = ne.measure_the_principal_symbol_degenerates()
    note = result["this_is_classical"]
    assert "No quantisation" in note
    assert "retracted" in note
    assert "calculation is unchanged" in note
    assert "graviton" not in ne.measure_the_negative_egb_ledger()["closed_by"]


def test_the_module_says_why_the_coefficient_had_to_be_derived():
    note = ne.measure_the_principal_symbol_degenerates()[
        "why_it_had_to_be_derived"]
    assert "maximally symmetric" in note
    assert "product" in note


def test_pushing_the_coupling_further_makes_the_exterior_exotic():
    """The exotic matter is relocated into the universe, not removed."""
    R = 1.0
    critical = ne.exterior_threshold(R)
    assert ne.exterior_matter_null_stress(1.5 * critical, R) < 0.0
    rows = {round(r["coupling"], 9): r
            for r in ne.measure_the_critical_exterior_is_empty()["rows"]}
    beyond = rows[round(1.5 * critical, 9)]
    assert beyond["sum"] < 0.0
    assert beyond["equation_of_state"] < -1.0


# ── the verdict ─────────────────────────────────────────────────────────────

def test_the_branch_is_closed_by_the_field_equations_not_the_matter():
    ledger = ne.measure_the_negative_egb_ledger()
    assert ledger["branch_is_closed"]
    assert "principal symbol" in ledger["closed_by"]
    assert "CLASSICAL" in ledger["closed_by"]
    assert "STRUCTURAL grounds" in ledger["verdict"]
    verdicts = " ".join(e["verdict"] for e in ledger["entries"]).lower()
    assert "coupling constant" in verdicts
    assert "nec, wec and dec all hold" in verdicts, \
        "the withdrawn exotic-matter claim must be recorded"
    assert "overreached" in verdicts, \
        "the withdrawn empty-universe claim must be recorded"


def test_the_ledger_names_what_it_does_not_test():
    """A source action, the scalar and vector sectors, dilatonic EGB, `f(R)`,
    and a different exterior."""
    note = ne.measure_the_negative_egb_ledger()["what_remains_untested"]
    assert "SOURCE" in note, "no source action was exhibited"
    assert "scalar and vector perturbation sectors" in note
    assert "dilatonic" in note
    assert "f(R)" in note
    assert "DIFFERENT exterior" in note


def test_the_ledger_keeps_five_branches_and_derives_its_numbers():
    ledger = ne.measure_the_negative_egb_ledger()
    assert len(ledger["the_remaining_branches"]) == 5
    scan = ne.measure_no_coupling_satisfies_both()
    entries = {e["claim"]: e for e in ledger["entries"]}
    row = next(e for k, e in entries.items() if "open set of couplings" in k)
    assert f"{scan['critical_coupling']:.6f}" in row["evidence"]


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import negative_egb_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
