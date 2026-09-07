"""Round 10 joint-closure tests. Freeze ``b78157a``."""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import joint_closure as jc
from experiments.closure_ledger import joint_closure_probe as probe


def _pair(gamma, sA=1, sB=1):
    return jc.Triangle((0.0, 0.0, 1.0), (math.sin(gamma), 0.0, math.cos(gamma)), sA, sB)


def test_triangle_rejects_collinear_non_unit_and_bad_signs():
    with pytest.raises(ValueError):
        jc.Triangle((0.0, 0.0, 1.0), (0.0, 0.0, 1.0), 1, 1)
    with pytest.raises(ValueError):
        jc.Triangle((0.0, 0.0, 2.0), (1.0, 0.0, 0.0), 1, 1)
    with pytest.raises(ValueError):
        jc.Triangle((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 0, 1)


@pytest.mark.parametrize("gamma", [0.35, 0.7, 1.0, math.pi / 2, 2.4, 2.7])
def test_closure_circle_has_vanishing_numerator_and_the_stated_D(gamma):
    for sA, sB in jc.SECTOR_SIGNS:
        tri = _pair(gamma, sA, sB)
        assert abs(float(np.linalg.norm(tri.s)) - math.sqrt(2 * tri.t)) < 1e-14
        for psi in np.linspace(0.0, 2 * math.pi, 41):
            N, D, _ = tri.invariants(tri.circle_point(psi))
            assert abs(N) < 1e-14
            assert abs(D - float(tri.D_of_psi(psi))) < 1e-13


def test_phase_gradient_is_correct_on_the_negative_D_branch():
    """Regression: ``atan2`` jumps by ``2 pi`` across ``D < 0``, and the whole
    closure circle sits on that cut. Without folding differences modulo ``pi``
    the measured gradient diverges as the step shrinks."""
    tri = _pair(1.0)
    negative = [p for p in np.linspace(0.05, 2 * math.pi - 0.05, 200)
                if float(tri.D_of_psi(p)) < -0.3]
    assert negative, "expected a negative-D arc at this angle"
    for psi in negative[::20]:
        x = tri.circle_point(psi)
        closed = float(np.linalg.norm(tri.q)) / abs(tri.invariants(x)[1])
        measured = float(np.linalg.norm(jc.phase_gradient(tri, x, 1e-5)))
        assert abs(measured - closed) / closed < 1e-7


def test_numerator_gradient_is_the_constant_round_eight_control():
    tri = _pair(1.3, -1, 1)
    for psi in np.linspace(0.0, 2 * math.pi, 17):
        g = float(np.linalg.norm(jc.numerator_gradient(tri, tri.circle_point(psi), 1e-5)))
        assert abs(g - float(np.linalg.norm(tri.q))) < 1e-8


@pytest.mark.parametrize("t", [0.2, 0.9, 1.0, 1.7, 1.999, 2.0, 2.5])
def test_W1_closed_matches_split_quadrature(t):
    assert abs(jc.W1_closed(t) - jc.W1_quadrature(t)) < 1e-11 * max(1.0, t)


def test_W1_closed_form_is_continuous_across_t_equals_two():
    below, above = jc.W1_closed(2.0 - 1e-9), jc.W1_closed(2.0)
    assert abs(below - above) < 1e-6
    assert abs(jc.W1_closed(1.0) - (math.pi + 4.0)) < 1e-13


@pytest.mark.parametrize("gamma", [0.35, 1.0, 2.4])
def test_punctures_are_minus_u_and_minus_w_with_slope_equal_to_q(gamma):
    tri = _pair(gamma)
    info = jc.puncture_geometry(tri)
    assert info["has_punctures"] == 1.0
    assert info["slope_minus_q"] < 1e-13
    assert max(info["puncture_is_minus_u"], info["puncture_is_minus_w"]) < 1e-13


def test_excised_coarea_mass_follows_the_two_eta_squared_law():
    result = jc.excision_masses(_pair(1.0))
    assert result["relative_error_vs_2eta2"] < 1e-3
    for row in result["rows"]:
        assert row["fraction"] < 1e-3


def test_joint_excision_bound_uses_inclusion_exclusion():
    t1, t2 = _pair(1.0), _pair(1.3, 1, -1)
    for row in jc.joint_excision_bound(t1, t2):
        assert row["joint_excluded_fraction"] >= max(row["factor1"], row["factor2"])
        assert row["joint_excluded_fraction"] <= row["factor1"] + row["factor2"]


def test_joint_gram_is_block_diagonal_and_second_order_accurate():
    t1, t2 = _pair(1.0), _pair(1.3, 1, -1)
    result = jc.gram_and_metric_checks(t1, t2)
    # Block-diagonality is structural: theta_i depends only on x_i. Recorded
    # as a regression oracle, not as independent evidence.
    assert result["max_offdiagonal_gram"] == 0.0
    errs = [result["per_step_relative_error"][s] for s in (1e-4, 5e-5, 2.5e-5)]
    assert result["finest_relative_error"] < 1e-6
    for coarse, fine in zip(errs[:-1], errs[1:]):
        assert 3.0 < coarse / fine < 5.0     # central differences are O(h^2)


def test_independent_analytic_route_reproduces_the_closed_jacobian():
    result = jc.analytic_gram_route(_pair(1.0), _pair(2.4, -1, 1))
    assert result["single_factor_relative_error"] < 1e-10
    assert result["joint_jacobian_relative_error"] < 1e-10


def test_density_is_covariant_under_a_common_frame_and_copy_exchange():
    a = (0.0, 0.0, 1.0)
    b1 = (math.sin(1.0), 0.0, math.cos(1.0))
    b2 = (math.sin(1.3), 0.0, math.cos(1.3))
    result = jc.covariance_and_exchange(a, b1, a, b2)
    assert max(result.values()) < 1e-10


def test_finite_windows_converge_quadratically_to_the_coarea_limit():
    a = (0.0, 0.0, 1.0)
    b1 = (math.sin(1.0), 0.0, math.cos(1.0))
    b2 = (math.sin(1.3), 0.0, math.cos(1.3))
    result = jc.window_convergence(a, b1, a, b2)
    assert result["final_window_discrepancy"] < 2e-3
    assert result["monotone"]
    symmetric = [r["max_discrepancy_from_coarea"] for r in result["rows"]
                 if r["epsilon"][0] == r["epsilon"][1]]
    for coarse, fine in zip(symmetric[:-1], symmetric[1:]):
        assert 3.0 < coarse / fine < 5.0


def test_sector_probabilities_agree_between_grids_and_the_closed_form():
    a = (0.0, 0.0, 1.0)
    b1 = (math.sin(0.35), 0.0, math.cos(0.35))
    b2 = (math.sin(0.7), 0.0, math.cos(0.7))
    assert jc.sector_probability_grids(a, b1, a, b2)["final_grid_error"] < 1e-6


def test_every_level_set_control_contains_distinct_joint_histories():
    """A level set is only informative if it holds more than one history."""
    controls = jc.level_set_controls()
    for key in ("reflection", "signed_product", "absolute_product"):
        assert controls[key]["chord_separation"] > 1e-3, key
    # the sign control deliberately holds factor 2 fixed; the product-space
    # separation must still be nonzero
    assert controls["absolute_product"]["chord_separation"] > 0.5


def test_reflection_preserves_the_pair_statistic_and_the_reference_density():
    c = jc.level_set_controls()["reflection"]
    assert c["statistic_match"] < 1e-12
    assert c["reference_density_gap"] < 1e-10
    assert c["cubic_gap"] < 1e-10          # the cubic is a function of the pair


def test_equal_products_separate_the_cubic_but_not_the_reference_density():
    c = jc.level_set_controls()["signed_product"]
    assert c["statistic_match"] < 1e-12
    assert c["reference_density_gap"] < 1e-10
    assert abs(c["cubic_gap"] - abs(jc.cubic_weight(1.0) * jc.cubic_weight(2.0)
                                    - jc.cubic_weight(math.sqrt(2.0)) ** 2)) < 1e-12
    assert c["cubic_gap"] > 1e-6


def test_opposite_sign_same_absolute_product_separates_the_cubic():
    c = jc.level_set_controls()["absolute_product"]
    assert c["statistic_match"] < 1e-12
    assert c["reference_density_gap"] < 1e-10
    # Phi(d) - Phi(-d) = -2 d^3 / 5, so the gap is Phi(1) * 2 (1/4)^3 / 5.
    assert abs(c["cubic_gap"] - jc.cubic_weight(1.0) * 2 * 0.25 ** 3 / 5) < 1e-14
    assert c["cubic_gap"] > 1e-6


def test_generic_closure_rule_accepts_a_union_of_non_closed_subsystems():
    r = jc.generic_closure_rule_is_rank_one()
    assert r["rank_one_cancellation_demonstrated"]
    assert r["union_is_closed"] and not any(r["subsystem_is_closed"])
    assert abs(r["union_total_phase"]) < 1e-12
    # and at the module default the phase gate cannot reject anything at all
    assert r["default_sigma_cannot_reject"]


def test_based_loop_additivity_requires_a_common_base_point():
    r = jc.based_loop_composition_scope()
    assert r["same_base_additivity_residual"] < 1e-12
    assert r["distinct_base_noncommutativity"] > 1e-3
    assert not r["theorem_applies_to_disconnected_pairs"]


def test_no_inspected_module_supplies_a_joint_weight_rule():
    audit = jc.repository_joint_rule_audit()
    assert not audit["any_module_supplies_joint_weight_rule"]
    assert len(audit["entries"]) == 5
    assert sum(e["applies_to_disconnected_pairs"] for e in audit["entries"]) == 1


def test_verdict_reports_unresolved_when_a_required_check_fails():
    audit = jc.repository_joint_rule_audit()
    controls = jc.level_set_controls()
    bad = jc.verdict({"anything": False}, audit, controls)
    assert set(bad.values()) == {"UNRESOLVED"}
    good = jc.verdict({"anything": True}, audit, controls)
    assert good["reference_composition"] == "INDEPENDENT_PHASE_PRODUCT_VERIFIED"
    assert good["reduction_of_specified_rules"] == "PRODUCT_STATISTIC_SUFFICIENT"
    assert good["additional_physical_rule"] == "JOINT_WEIGHT_RULE_UNSPECIFIED"
    assert good["consequence_for_selection"] == "NO_SELECTION_FROM_INDEPENDENCE"
    assert good["P3_hypotheses_established"].startswith("NO")


def test_structural_regressions_are_labelled_and_not_counted_as_evidence():
    report = probe.run_probe()
    assert report["checks_passed"]
    assert all(report["checks"].values())
    assert all("structural" in k for k in report["structural_regressions"])
    assert all(report["structural_regressions"].values())
    assert len(report["checks"]) == 15
    assert not any("structural" in k for k in report["checks"])
    assert report["weight_selection"] == "NOT_DERIVED"
    assert report["operational_source_readout"] == "NOT_DERIVED"
    assert jc.PUBLIC_PREREG in report["public_preregistration"]
