"""Round 7 — the history-action question, pre-registered in
``docs/history_action_prereg.md`` (`a33a901`) before the module existed."""

import math

import numpy as np
import pytest

from geometrodynamics.bulk.history_action import (
    holonomy_trace, morse_bott_oracle, class_function_degeneracy,
    additive_functionals_have_no_critical_points, saddle_branch_ratio,
    sector_symmetry_group, sector_orbits, fibre_action_is_weight_blind,
    detector_response_homogeneity, quadratic_readout_law,
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
def test_F2_ratio_is_real_iff_4kappa_over_pi_is_an_odd_integer(kappa, expect_real):
    r = saddle_branch_ratio(kappa)
    assert r["ratio_is_real"] is expect_real
    assert r["odd_multiple_of_pi_over_4"] is expect_real


def test_F1_the_saddle_magnitudes_never_separate_the_branches():
    for kappa in (0.1, math.pi / 4, 1.0, 3 * math.pi / 4, 5.0):
        assert abs(saddle_branch_ratio(kappa)["magnitude"] - 1.0) < 1e-12


def test_F3_the_branch_sign_alternates_with_kappa_so_stationarity_cannot_choose():
    assert saddle_branch_ratio(math.pi / 4)["selects"].startswith("positive")
    assert saddle_branch_ratio(3 * math.pi / 4)["selects"].startswith("oriented")
    # the one value with any claim to naturalness is neither branch
    assert saddle_branch_ratio(1.0)["selects"].startswith("neither")


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

def test_every_bam_coupling_is_degree_two_homogeneous():
    c = detector_response_homogeneity()
    assert c["all_quadratic"] and not c["any_linear"]
    for row in c["rows"]:
        assert abs(row["measured_degree"] - 2.0) < 1e-6, row


def test_the_quadratic_readout_is_superquantum():
    q = quadratic_readout_law()
    assert q["quadratic_exceeds_tsirelson"]
    assert abs(q["S_max_quadratic"] - 8.0 * math.sqrt(2.0) / 3.0) < 1e-4
    assert abs(q["S_max_linear"] - 2.0 * math.sqrt(2.0)) < 1e-8


def test_the_quadratic_readout_still_has_exact_half_marginals():
    assert quadratic_readout_law()["quadratic_marginal_deviation"] < 1e-12


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


def test_odd_source_observables_are_blind_but_even_ones_signal():
    e = source_observable_signalling(n=20000)
    assert e["odd_observables_are_blind"]
    assert e["even_observables_signal"]
    assert e["source_readout_signals"]


def test_non_coplanar_settings_give_mutually_singular_source_measures():
    e = source_observable_signalling(n=20000)
    assert e["non_coplanar_supports_mutually_singular"]
    assert abs(float(e["non_coplanar_total_variation"]) - 1.0) < 1e-9


# ── the report ──────────────────────────────────────────────────────────────

def test_the_five_verdicts_are_the_pre_registered_labels():
    v = verdicts()
    assert v["A_action"] == "HOLONOMY_TRACE_IS_A_STATIONARY_FUNCTIONAL_NOT_A_DERIVED_ACTION"
    assert v["B_sectors"] == "LIKE_UNLIKE_SECTOR_RATIO_REMAINS_FREE"
    assert v["C_readout"] == "CLASSICAL_DETECTOR_RESPONDS_QUADRATICALLY"
    assert v["D_compatibility"] == "HISTORY_ACTION_INDEPENDENTLY_POSTULATED"
    assert v["E_causality"] == "SOURCE_READOUT_SIGNALS_FUTURE_SETTINGS"


def test_the_headline_is_not_printed():
    assert "CLASSICAL_BAM_DERIVES_THE_SINGLET_PROBABILITY_LAW" not in str(verdicts())


def test_the_ledger_records_kappa_and_the_sector_prior_as_open():
    led = {row["input"]: row["status"] for row in dependency_ledger()}
    assert led["kappa, the normalisation in e^{i kappa S_H}"] == "open"
    assert led["equal sector prior r = 1"] == "open"
    assert led["linear current-to-frequency readout"] == "open"
    assert led["antipodal scalar BC, eta, quotient-vs-cover"] == "not used"
