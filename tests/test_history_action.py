"""Round 7 — the history-action question, pre-registered in
``docs/history_action_prereg.md`` (`a33a901`) before the module existed."""

import math

import numpy as np
import pytest

from geometrodynamics.bulk.history_action import (
    holonomy_trace, morse_bott_oracle, class_function_degeneracy,
    additive_functionals_have_no_critical_points, saddle_branch_ratio,
    morse_bott_component_masses, no_off_closure_critical_points,
    amplitude_dependence, excision_estimate, discrete_symmetry_extension,
    sector_symmetry_group, sector_orbits, fibre_action_is_weight_blind,
    detector_response_homogeneity, quadratic_readouts_disagree,
    local_square_mean_is_closed_form,
    radial_action_compatibility, source_observable_signalling,
    dependency_ledger, verdicts)


# ── A: the oracle frozen in the pre-registration ────────────────────────────

def test_the_frozen_oracle_reproduces_O1_to_O4():
    o = morse_bott_oracle()
    assert o["O1_trace_identity"] < 1e-12
    assert o["O1_closed_form"] < 1e-12
    assert o["O2_grad_on_closure"] < 1e-4          # finite-difference limited
    assert o["O3_transverse_curvature_rel"] < 1e-3
    assert o["O3_index_is_sign_D"]
    assert o["O4_saddle_magnitude_is_coarea_rel"] < 1e-3
    assert o["passes"]


def test_the_saddle_magnitude_is_the_positive_coarea_density():
    """O4 stated on its own: this is why S_H is the candidate at all."""
    assert morse_bott_oracle()["O4_saddle_magnitude_is_coarea_rel"] < 1e-3


def test_critical_set_agreement_is_not_evidence_of_selection():
    """Pre-registered non-evidence control: every F(cos theta) shares Crit."""
    d = class_function_degeneracy()
    assert d["all_share_the_critical_set"]
    assert len(d["rows"]) >= 5


def test_theta_is_additive_and_S_H_is_not():
    a = additive_functionals_have_no_critical_points()
    assert a["theta_is_additive_residual"] < 1e-12
    assert a["S_H_additivity_defect"] > 1e-3
    assert a["theta_additive"] and a["S_H_not_additive"]


def test_the_additive_functional_has_no_critical_point_on_closure():
    a = additive_functionals_have_no_critical_points()
    assert a["theta_has_no_critical_point"]
    assert a["min_|grad theta|_on_closure"] > 1e-6


@pytest.mark.parametrize("kappa,expect_real", [
    (math.pi / 4, True), (3 * math.pi / 4, True),
    (math.pi / 2, False), (1.0, False), (0.3, False)])
def test_F2_phase_factor_is_real_iff_4kappa_over_pi_is_an_odd_integer(kappa, expect_real):
    r = saddle_branch_ratio(kappa)
    assert r["phase_factor_is_real"] is expect_real
    assert r["odd_multiple_of_pi_over_4"] is expect_real


def test_the_component_masses_are_generically_unequal():
    """Morse-Bott stationary phase integrates over each component, so the
    ratio is M_pi/M_0, not 1. An earlier version claimed unit magnitude by
    mistaking the local phase prefactor for the component amplitude."""
    m = morse_bott_component_masses([0, 0, 1],
                                    [math.sin(1.0), 0, math.cos(1.0)])
    assert m["masses_are_unequal"]
    assert abs(m["mass_ratio"] - 1.0) > 0.5
    r = saddle_branch_ratio(math.pi / 4)
    assert abs(r["magnitude"] - r["mass_ratio"]) < 1e-12
    assert abs(r["magnitude"] - 1.0) > 0.5


def test_the_masses_are_conditional_on_the_round5_haar_measure():
    """The Morse-Bott coefficients use a = 1; a different preparation density
    moves them, so "nothing tuned" applies to the Hessian, not the measure."""
    a = amplitude_dependence()
    assert a["in_plane_amplitude_moves_the_aggregation"]
    assert a["masses_are_conditional_on_the_measure"]
    # an amplitude varying only along u x v is invisible on the closure circle
    assert a["normal_amplitude_is_invisible"]


def test_the_punctures_carry_no_mass_so_excision_is_safe():
    """S_H is singular at x = -u, -v; the excised mass is O(eps^2), so the
    Morse-Bott coefficients are genuine rather than formal."""
    e = excision_estimate()
    assert e["mass_is_O_eps_squared"]
    assert e["excision_is_safe"]
    assert e["worst_relative_error"] < 0.01


def test_no_identified_discrete_operation_mixes_like_and_unlike():
    """Detector exchange, history reversal and the Pin deck all preserve
    s_A s_B - an enumeration, explicitly not a classification."""
    d = discrete_symmetry_extension()
    assert d["all_preserve_like_unlike"]
    assert d["no_identified_operation_mixes_like_and_unlike"]
    assert d["is_a_classification_of_all_symmetries"] is False
    assert all(r["max_weight_change"] < 1e-9 for r in d["rows"])


def test_the_masses_reproduce_both_candidate_aggregations_exactly():
    """(M_0 - M_pi)|uxv| is the oriented sum and (M_0 + M_pi)|uxv| the positive
    count: stationary phase supplies both magnitudes, leaving only the phase."""
    rng = np.random.default_rng(2)
    for _ in range(4):
        u, v = (w / np.linalg.norm(w) for w in rng.normal(size=(2, 3)))
        m = morse_bott_component_masses(u, v)
        assert m["oriented_identity_residual"] < 1e-9
        assert m["positive_count_identity_residual"] < 1e-9


def test_the_undetermined_part_is_the_relative_phase_only():
    for kappa in (0.1, math.pi / 4, 1.0, 3 * math.pi / 4, 5.0):
        r = saddle_branch_ratio(kappa)
        assert abs(abs(r["phase_factor"]) - 1.0) < 1e-12   # phase only
        assert abs(r["arg_phase_factor"]
                   - math.atan2(math.sin(2 * kappa - math.pi / 2),
                                math.cos(2 * kappa - math.pi / 2))) < 1e-12


def test_F3_the_branch_sign_alternates_with_kappa_so_stationarity_cannot_choose():
    assert saddle_branch_ratio(math.pi / 4)["selects"].startswith("positive")
    assert saddle_branch_ratio(3 * math.pi / 4)["selects"].startswith("oriented")
    # the one value with any claim to naturalness is neither branch
    assert saddle_branch_ratio(1.0)["selects"].startswith("neither")


def test_there_are_no_critical_points_off_the_closure_set():
    """The half of Crit(S_H) = Gamma that checking on the closure set omits."""
    o = no_off_closure_critical_points(trials=6)
    assert o["min_grad_theta_off_closure"] > 1e-3
    assert o["required_x_p_equals_minus_sec_half_gamma"] < 1e-12
    assert o["every_required_x_p_exceeds_the_sphere"]
    assert o["no_off_closure_critical_points"]


def test_S_H_closed_form_matches_the_trace_definition():
    rng = np.random.default_rng(0)
    for _ in range(50):
        x, u, v = (rng.normal(size=3) for _ in range(3))
        x, u, v = (w / np.linalg.norm(w) for w in (x, u, v))
        from geometrodynamics.bulk.closure_measurement import solid_angle
        om = float(solid_angle(x, u, v)[0])
        assert abs(float(holonomy_trace(x, u, v)[0]) + math.cos(0.5 * om)) < 1e-14


# ── B: the sector coefficients ──────────────────────────────────────────────

def test_the_fixed_setting_group_has_two_sector_orbits_at_generic_angles():
    g = sector_symmetry_group([1, 0, 0], [math.cos(1.0), math.sin(1.0), 0])
    assert g["group_order"] == 4
    assert g["n_orbits"] == 2
    assert not g["r_forced"]
    assert sorted(map(sorted, g["orbits"])) == [[(-1, -1), (1, 1)],
                                                [(-1, 1), (1, -1)]]


def test_equal_sector_measure_is_forced_only_at_a_right_angle():
    o = sector_orbits()
    assert o["forced_only_at_right_angle"]
    assert not o["forced_at_any_chsh_angle"]


def test_orthogonal_settings_are_the_one_transitive_case():
    g = sector_symmetry_group([1, 0, 0], [0, 1, 0])
    assert g["transitive_on_sectors"] and g["r_forced"]
    assert g["group_order"] == 8


def test_a_larger_group_along_the_fibre_cannot_change_the_weights():
    assert fibre_action_is_weight_blind()["fibre_blind"]


# ── C: the readout ──────────────────────────────────────────────────────────

def test_all_six_audited_quantities_are_degree_two_homogeneous():
    """Six specific quantities, not every observable in the repository."""
    c = detector_response_homogeneity()
    assert len(c["rows"]) == 6
    assert c["all_quadratic"] and not c["any_linear"]
    for row in c["rows"]:
        assert abs(row["measured_degree"] - 2.0) < 1e-6, row


def test_two_ordinary_quadratic_readouts_disagree():
    """Degree-2 homogeneity does not name a readout, so C cannot be reported
    as a derived quadratic law."""
    q = quadratic_readouts_disagree()
    assert q["the_two_quadratics_disagree"]
    assert abs(q["square_of_integral"]["S_max"] - 8.0 * math.sqrt(2.0) / 3.0) < 1e-4
    assert abs(q["integral_of_square"]["S_max"]
               - 12.0 * math.sqrt(2.0) / 5.0) < 1e-4
    assert abs(q["closed_form_integral_of_square"]
               - 12.0 * math.sqrt(2.0) / 5.0) < 1e-12
    assert abs(q["linear"]["S_max"] - 2.0 * math.sqrt(2.0)) < 1e-8
    assert q["both_exceed_tsirelson"]


def test_both_quadratic_readouts_keep_exact_half_marginals():
    q = quadratic_readouts_disagree()
    assert q["square_of_integral"]["marginal_dev"] < 1e-12
    assert q["integral_of_square"]["marginal_dev"] < 1e-12


def test_the_local_square_mean_closed_form():
    assert local_square_mean_is_closed_form()["closed_form_holds"]


# ── D: compatibility with the existing radial action ────────────────────────

def test_the_radial_action_is_a_genuine_one_form_integral():
    d = radial_action_compatibility()
    assert d["both_legs_nonzero"], "the split must lie inside the allowed region"
    assert d["radial_action_additive_defect"] < 1e-6
    assert d["radial_action_is_a_one_form_integral"]


def test_the_module_agrees_with_the_repository_radial_action():
    d = radial_action_compatibility()
    assert d["agrees_with_repository_one_way"] < 1e-9
    assert d["closed_orbit_is_twice_one_way"] < 1e-9


# ── E: the causality gate ───────────────────────────────────────────────────

def test_the_conditioned_source_density_is_exactly_antipodally_even():
    e = source_observable_signalling(n=20000)
    assert e["density_is_antipodally_even_residual"] < 1e-15


def test_odd_means_and_signs_are_blind_but_full_distributions_are_not():
    e = source_observable_signalling(n=20000)
    assert e["odd_means_and_signs_are_blind"]
    assert not e["odd_full_distributions_are_blind"]
    assert e["odd_projection_variance_spread"] > 1e-3
    assert e["some_even_functions_separate_the_ensembles"]
    assert e["setting_information_present_at_source"]


def test_not_every_even_observable_separates():
    """The claim "every even observable signals" is false: constants and x.x
    are even and blind."""
    e = source_observable_signalling(n=20000)
    assert e["not_every_even_observable_separates"]
    assert e["blind_even_observable_spread"] < 1e-12


def test_no_operational_channel_is_claimed():
    """Neither the map from x to field configurations nor a source-local
    readout compatible with the two-boundary problem is constructed."""
    e = source_observable_signalling(n=20000)
    assert e["bam_couplings_shown_even_in_x"] is False
    assert e["operational_readout_constructed"] is False


def test_non_coplanar_settings_give_mutually_singular_source_measures():
    e = source_observable_signalling(n=20000)
    assert e["non_coplanar_supports_mutually_singular"]
    assert abs(float(e["non_coplanar_total_variation"]) - 1.0) < 1e-9


# ── the report ──────────────────────────────────────────────────────────────

def test_the_five_verdicts_are_the_pre_registered_labels():
    v = verdicts()
    assert v["A_action"] == "HOLONOMY_TRACE_IS_A_STATIONARY_FUNCTIONAL_NOT_A_DERIVED_ACTION"
    assert v["B_sectors"] == "NO_IDENTIFIED_SYMMETRY_FORCES_EQUAL_SECTOR_MEASURE"
    assert v["C_readout"] == "NO_BAM_DETECTOR_COUPLING_CURRENTLY_DEFINES_THE_READOUT"
    assert v["D_compatibility"] == "HISTORY_ACTION_INDEPENDENTLY_POSTULATED"
    assert v["E_causality"] == "SETTING_INFORMATION_IS_PRESENT_AT_SOURCE_READOUT_DYNAMICS_OPEN"


def test_the_headline_is_not_printed():
    assert "CLASSICAL_BAM_DERIVES_THE_SINGLET_PROBABILITY_LAW" not in str(verdicts())


def test_the_ledger_records_kappa_and_the_sector_prior_as_open():
    led = {row["input"]: row["status"] for row in dependency_ledger()}
    assert led["kappa, the normalisation in e^{i kappa S_H}"] == "open"
    assert led["equal sector prior r = 1"] == "open"
    assert led["current-to-frequency readout"] == "open"
    assert led["a map from the source variable x to field configurations"] == "open"
    assert led["antipodal scalar BC, eta, quotient-vs-cover"] == "not used"
