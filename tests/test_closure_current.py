"""The closure-current fork round, against `docs/closure_current_prereg.md`."""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import closure_current as cc


def test_R1_pin_loop_reduces_to_the_triangle_with_the_partner_sign():
    r = cc.pin_loop_reduction(n=200)
    assert r["q5_equals_minus_q0_Ginv"] < 1e-12
    assert r["holonomy_fibre_independent"] < 1e-12
    assert r["frame_holonomy_is_Omega_of_triangle_x_u_minus_v"] < 1e-12
    assert r["lift_is_cos_half_Omega_plus_sin_half_Omega_x"] < 1e-12
    assert r["reduces_to_triangle"] and r["extra_from_J_squared"] == -1


def test_C1_branch_holonomy_is_the_sign_of_D_on_the_closure_set():
    r = cc.branch_holonomy_is_sign_D()
    assert r["exp_i_half_Omega_equals_sign_D"] < 1e-12
    assert r["loop_lift_equals_sign_D_on_circle"] < 1e-10
    assert r["signed_density_is_holonomy_times_coarea"]
    assert r["pi_branch_fraction_by_sector"][(1, 1)] == pytest.approx(1.0 / (2 * math.pi), abs=2e-3)


def test_R2_holonomy_weighted_current_is_the_quantum_joint_law_without_projectors():
    for g in (0.3, 1.0, 2.0):
        h = cc.holonomy_weighted_law(g)
        assert h["max_deviation"] < 1e-9 and h["E_is_cos"] and h["no_projectors_used"]
        assert h["circle_mean_of_x"] < 1e-12


def test_R3_sector_prior_moves_E_but_not_the_marginals():
    r = cc.sector_prior_control()
    assert r["marginals_stay_half"] and r["E_moves"]
    Es = [row["E"] for row in r["rows"]]
    assert Es == pytest.approx([0.0751, 0.3985, 0.6460], abs=2e-3)   # closed form; prereg values corrected
    assert r["symmetry_fixing_ratio"] is None


def test_R4_stationarity_is_unimplemented_and_its_set_is_not_the_closure_set():
    a = cc.stationarity_audit()
    assert a["named_in_module_docstring"] and not a["implemented"]
    s = a["stationary_set"]
    assert s["no_stationary_points"] and s["min_grad_on_grid"] > 0.1     # Lexell: no critical points
    assert not s["coincides_with_closure_set"]
    assert s["lexell"]["level_set_is_a_circle (plane residual)"] < 5e-3
    assert max(s["lexell"]["level_set_plane_contains_minus_u_minus_v (plane distances)"]) < 1e-2
    assert max(s["lexell"]["D_at_minus_u_minus_v"]) < 1e-12
    assert not a["stationarity_decides_the_fork"]


def test_R5_the_two_candidates_and_the_verdict_rule():
    c = cc.pin_label_versus_weight()["candidates"]
    assert c["POSITIVE_COAREA"]["S_max"] == pytest.approx(2.1423, abs=2e-3)
    assert c["HOLONOMY_WEIGHTED_COAREA"]["S_max"] == pytest.approx(2 * math.sqrt(2), abs=1e-3)
    assert cc.verdict(True, True, True, False, False) == "FORK_UNRESOLVED_BY_CURRENT_STRUCTURES"
    assert cc.verdict(True, True, True, True, False) == "HOLONOMY_WEIGHTED_COAREA"
    assert cc.verdict(False, True, True, False, False) == "NEITHER"
