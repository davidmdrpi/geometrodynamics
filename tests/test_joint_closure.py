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


@pytest.mark.parametrize("gamma,sA,sB", [(0.35, 1, -1), (0.7, 1, -1), (1.0, 1, 1),
                                          (2.7, 1, 1), (math.pi / 2, 1, 1)])
def test_excision_uses_the_registered_domain_and_the_two_term_law(gamma, sA, sB):
    """N23 regression: the freeze excises ``|D|/|q| < eta``, not
    ``|psi - psi0| < eta``. Those agree only to leading order, so the excised
    mass is ``2 eta^2`` asymptotically and NOT exactly."""
    tri = _pair(gamma, sA, sB)
    result = jc.excision_masses(tri)
    assert result["max_relative_residual_vs_two_term"] < 1e-3
    assert result["residual_scales_as_eta6"]
    assert result["quartic_improvement_factor"] > 10.0
    for row in result["rows"]:
        assert row["fraction"] < 1e-3
        # the arc really is the |D|/|q| level set, not a fixed psi window
        for edge in row["arc"]:
            assert abs(abs(float(tri.D_of_psi(edge)))
                       / float(np.linalg.norm(tri.q)) - row["eta"]) < 1e-12


def test_the_excision_domains_differ_at_third_order():
    """The two domains are not interchangeable at the frozen widths."""
    tri = _pair(0.7, 1, -1)
    psi0 = math.acos(-tri.t / math.sqrt(2 * tri.t))
    lo, hi = jc.puncture_arc_bounds(tri, psi0, 0.02)
    assert abs((hi - lo) - 0.04) > 1e-6          # not the naive 2 eta arc
    measured = jc.excision_masses(tri, (0.02,))["rows"][0]
    assert measured["deviation_from_leading"] > 1e-3   # would fail a 2 eta^2 gate


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


def test_on_closure_holonomies_are_central_and_respect_the_window_bound():
    """N24 regression: the first version reported generic ``SU(2)``
    non-commutativity as the obstruction. On the closure locus
    ``theta_i in pi Z``, so the reduced holonomies are ``+-1`` and commute
    exactly; the generic figure is an off-closure control only."""
    r = jc.based_loop_composition_scope()
    assert r["same_base_additivity_residual"] < 1e-12
    assert r["on_closure_commutator"] < 1e-20
    assert r["on_closure_holonomy_is_central"] < 1e-12
    assert r["windows_respect_bound"]
    for row in r["window_rows"]:
        e1, e2 = row["epsilon"]
        assert row["max_commutator"] <= 2 * math.sin(e1) * math.sin(e2)
    assert r["off_closure_control_is_not_evidence"]
    assert not r["common_frame_transport_derived"]
    assert r["supplies_a_closure_condition_not_a_weight"]
    assert "composition_on_closure_is_rank_one" not in r   # withdrawn by N27


def test_composed_holonomy_condition_is_rank_two_not_rank_one():
    """N27 regression: centrality gives commutation ON closure, not its
    continuation away from it. ``d vec(G1 G2) = c1 c2 (x1 dtheta1 +
    x2 dtheta2)`` has rank two for non-parallel axes, so opposite small
    phases do not cancel and the scalar-sum rank-one defect does not
    transfer to the holonomy product."""
    r = jc.composed_holonomy_rank()
    assert r["rank_one_claim_withdrawn"]
    assert r["min_rank"] == 2
    assert r["max_analytic_residual"] < 1e-9
    assert r["min_opposite_phase_vec_norm"] > 1e-3
    for row in r["rows"]:
        assert row["rank"] == 2
        assert abs(row["axis_dot"]) < 1.0 - 1e-6      # non-parallel axes
        assert min(row["singular_values"]) > 1e-6     # genuinely rank two


def test_reviewed_configuration_reproduces_the_reported_singular_values():
    """The reviewed case: both triangles at ``t = 1`` with ``D = 2``, whose
    closure axes are orthogonal and whose normal Jacobian is ``0.5, 0.5``."""
    reviewed = jc.composed_holonomy_rank()["reviewed_configuration"]
    assert abs(reviewed["axis_dot"]) < 1e-12
    assert reviewed["rank"] == 2
    for value in reviewed["normal_singular_values"]:
        assert abs(value - 0.5) < 1e-9


def test_off_closure_control_is_a_real_commutator():
    """N29 regression: the previous version drew a fresh axis in each of the
    four slots, computing ``G1 G2 - G3 G4`` rather than ``[G1, G2]``."""
    import inspect
    source = inspect.getsource(jc.based_loop_composition_scope)
    tail = source[source.index("off_closure = 0.0"):]
    assert tail.count("_unit(rng.normal(size=3))") == 2
    assert "H1" in tail and "H2" in tail
    value = jc.based_loop_composition_scope()["off_closure_control_commutator"]
    assert 0.0 < value <= 2.0        # |[G1,G2]| = 2 |s1 s2| |x1 x x2| <= 2


def test_no_inspected_module_supplies_a_joint_weight_rule():
    audit = jc.repository_joint_rule_audit()
    assert not audit["any_module_supplies_joint_weight_rule"]
    assert len(audit["entries"]) == 5
    assert sum(e["applies_to_disconnected_pairs"] for e in audit["entries"]) == 2


def test_window_slice_finds_a_gap_narrower_than_any_sample_grid():
    """N28 regression: a 257-point scan only finds gaps containing a sample
    point, and missed this one entirely, returning ``(2.0, 1)``."""
    tri = _pair(0.1, 1, -1)
    psi, eps = 2.613587734905, 0.099365561552
    measure, components = jc.window_slice(tri, psi, eps)
    assert components == 2
    assert abs(measure - 1.997544152) < 1e-8
    edges = jc.window_boundary_z(tri, psi, eps)
    assert len(edges) == 2
    assert abs(edges[0] - 0.501339412) < 1e-8
    assert abs(edges[1] - 0.502567336) < 1e-8
    # the boundary really solves |z||q| = |D| tan(eps)
    qn = float(np.linalg.norm(tri.q))
    A = math.sqrt(2 * tri.t) * math.cos(psi)
    for z in edges:
        assert abs(qn * z - abs(tri.t + math.sqrt(1 - z * z) * A) * math.tan(eps)) < 1e-12


def test_frozen_grid_guarantees_a_single_window_interval():
    """``|q|/t > tan(epsilon)`` forces one interval on the positive half."""
    worst = min(float(np.linalg.norm(_pair(g, sA, sB).q)) / _pair(g, sA, sB).t
                for g in (0.35, 0.7, 1.0, 1.3, math.pi / 2, 2.4, 2.7)
                for sA, sB in jc.SECTOR_SIGNS)
    assert worst > math.tan(0.08)


def test_window_slice_handles_a_disconnected_accepted_set():
    """N25 regression: at ``gamma = 0.1``, sector ``(+,-)``, ``psi = pi``,
    ``epsilon = 0.1`` the accepted set is a neighbourhood of ``z = 0`` plus one
    of ``z = 1``. The first version returned the whole interval."""
    tri = _pair(0.1, 1, -1)
    measure, components = jc.window_slice(tri, math.pi, 0.1)
    assert components == 2
    assert measure < 1.0                       # not the full |z| <= 1 interval
    Aq = float(np.linalg.norm(tri.q))
    root = math.sqrt(2 * tri.t) * math.cos(math.pi)
    g = lambda z: Aq * z - abs(tri.t + math.sqrt(max(1 - z * z, 0.0)) * root) * math.tan(0.1)
    assert g(0.0) < 0 and g(0.2) > 0 and g(1.0) < 0   # genuinely disconnected
    assert not hasattr(jc, "_negative_set_measure")   # the scan is gone (N28)


def test_registered_windows_keep_the_accepted_set_connected():
    a = (0.0, 0.0, 1.0)
    b1 = (math.sin(0.35), 0.0, math.cos(0.35))
    b2 = (math.sin(0.7), 0.0, math.cos(0.7))
    assert jc.window_convergence(a, b1, a, b2)["max_accepted_components"] == 1


def test_verdict_reports_unresolved_when_a_required_check_fails():
    audit = jc.repository_joint_rule_audit()
    controls = jc.level_set_controls()
    bad = jc.verdict({"anything": False}, audit, controls)
    good_keys = set(jc.verdict({"anything": True}, audit, controls))
    # N26: the failure branch must carry every field the renderer reads
    assert set(bad) == good_keys
    assert bad["reference_composition"] == "UNRESOLVED"
    assert bad["additional_physical_rule"] == "UNRESOLVED"
    assert bad["consequence_for_selection"] == "UNRESOLVED"
    assert bad["failed_checks"] == ["anything"]
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
    assert len(report["checks"]) == 17
    assert not any("structural" in k for k in report["checks"])
    assert report["weight_selection"] == "NOT_DERIVED"
    assert report["operational_source_readout"] == "NOT_DERIVED"
    assert jc.PUBLIC_PREREG in report["public_preregistration"]


def test_failing_cli_path_renders_archives_and_exits_nonzero(tmp_path, monkeypatch):
    """N26 regression: a failed mandatory check used to raise ``KeyError`` in
    the renderer, so no report was written and a stale passing archive
    survived. The failing path must render, archive and exit nonzero."""
    real = probe.run_probe

    def failing(progress=lambda s: None):
        report = real(progress=progress)
        report["checks"]["Q1 injected failure"] = False
        report["checks_passed"] = False
        report["verdict"] = jc.verdict(report["checks"], report["repository_audit"],
                                       report["level_set_controls"])
        return report

    monkeypatch.setattr(probe, "run_probe", failing)
    out = tmp_path / "run"
    out.mkdir()
    (out / "probe.md").write_text("STALE PASSING ARCHIVE")
    assert probe.main(["--output-dir", str(out)]) == 1
    written = (out / "probe.md").read_text()
    assert "STALE" not in written
    assert "UNRESOLVED" in written
    assert "Q1 injected failure" in written
    import json
    archived = json.loads((out / "probe.json").read_text())
    assert archived["verdict"]["consequence_for_selection"] == "UNRESOLVED"
    assert archived["checks_passed"] is False
