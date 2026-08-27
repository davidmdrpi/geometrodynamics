"""Is the v4 CKM realization a prediction, or a locally flexible fit?

Four kinds of check. The *observable* ones pin that the four flavor
coordinates are genuinely independent and rephasing-invariant — the whole round
is worthless if the Jacobian is reading LAPACK's phase convention. The
*numerical* ones pin that both maps differentiate. The *structural* ones pin
the rank, which is the round's entire conclusion and is chart-independent. And
the *scoping* ones guard the negative verdicts against softening.
"""

import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.qcd import flavor_identifiability as fi
from geometrodynamics.qcd.quark_spectrum import (LOCKED_QUARK_PARAMS_V4,
                                                 extract_ckm_matrix,
                                                 extract_physical_spectrum)


def _read(rel: str) -> str:
    return (Path(__file__).resolve().parents[1] / rel).read_text()


# ── the observables ─────────────────────────────────────────────────────────

def test_there_are_exactly_four_flavor_coordinates():
    """A unitary 3×3 modulo rephasing carries four real parameters. No more."""
    assert len(fi.FLAVOR_OBSERVABLES) == fi.PHYSICAL_FLAVOR_DIMENSION == 4
    assert len(fi.MASS_OBSERVABLES) == 4


def test_the_ckm_is_unitary_so_the_ceiling_is_real():
    V = extract_ckm_matrix(LOCKED_QUARK_PARAMS_V4)
    assert np.max(np.abs(V.conj().T @ V - np.eye(3))) < 1e-12


def test_the_flavor_observables_are_rephasing_invariant():
    """Both eigenbases carry arbitrary per-column phases from `eigh`.

    If the extraction read them, every Jacobian column would be a function of
    LAPACK's convention rather than of the physics.
    """
    V = extract_ckm_matrix(LOCKED_QUARK_PARAMS_V4)
    j0 = fi.jarlskog()
    rng = np.random.default_rng(3)
    for _ in range(5):
        left = np.diag(np.exp(1j * rng.uniform(0, 2 * math.pi, 3)))
        right = np.diag(np.exp(1j * rng.uniform(0, 2 * math.pi, 3)))
        W = left @ V @ right
        assert np.max(np.abs(np.abs(W) - np.abs(V))) < 1e-12
        j_w = float(np.imag(W[0, 1] * W[1, 2]
                            * np.conj(W[0, 2]) * np.conj(W[1, 1])))
        assert abs(j_w - j0) < 1e-12


def test_the_flavor_observables_land_near_the_measured_ckm():
    """Not a fit check — a sanity check that we are reading the right object."""
    theta12, theta23, theta13, delta = fi.flavor_observables()
    assert theta12 == pytest.approx(0.2265, abs=0.01)
    assert theta23 == pytest.approx(0.0424, abs=0.005)
    assert theta13 == pytest.approx(0.00365, abs=0.001)
    assert delta == pytest.approx(1.144, abs=0.05)
    assert fi.jarlskog() == pytest.approx(3.08e-5, rel=0.1)


def test_the_derived_quantities_are_not_extra_observable_rows():
    obs = fi.measure_the_observables_are_four_and_rephasing_invariant()
    assert "NOT independent observable rows" in obs["validation_outputs_not_rows"]
    assert obs["physical_flavor_dimension"] == 4


# ── the maps differentiate ──────────────────────────────────────────────────

def test_both_jacobians_converge_under_step_refinement():
    result = fi.measure_both_jacobians_converge()
    assert result["both_converged"]
    assert result["ranks_stable"]
    assert result["mass_relative_spread"] < 1e-5
    assert result["flavor_relative_spread"] < 1e-5


def test_the_two_jacobians_share_one_parameter_chart():
    JM = fi.jacobian("mass")
    JF = fi.jacobian("flavor")
    assert JM.shape == JF.shape == (4, len(fi.COORD_NAMES))
    assert len(fi.COORD_NAMES) == 18


def test_the_v4_lock_inherits_the_v3_masses():
    """The holonomy is a pure phase; `extract_physical_spectrum` strips it."""
    from geometrodynamics.qcd.quark_spectrum import LOCKED_QUARK_PARAMS
    v3 = extract_physical_spectrum(LOCKED_QUARK_PARAMS)
    v4 = extract_physical_spectrum(LOCKED_QUARK_PARAMS_V4)
    for s in fi.MASS_OBSERVABLES:
        assert v4[s] / v4["d"] == pytest.approx(v3[s] / v3["d"], rel=1e-6)


# ── the structural result ───────────────────────────────────────────────────

def test_the_kernel_really_is_a_kernel():
    N = fi.mass_preserving_tangent_space()
    JM = fi.jacobian("mass")
    assert N.shape[1] == 14
    assert np.max(np.abs(JM @ N)) < 1e-10


def test_rank_K_is_four_the_full_physical_flavor_space():
    """The round's entire conclusion. Rank is chart-independent."""
    result = fi.measure_the_mass_preserving_flavor_rank()
    assert result["rank_K"] == fi.PHYSICAL_FLAVOR_DIMENSION == 4
    assert result["spans_the_whole_physical_flavor_space"]
    assert result["rank_stable_across_the_grid"]
    assert len(result["grid"]) == 9
    for row in result["grid"]:
        assert row["rank_K"] == 4


def test_the_rank_is_not_propped_up_by_a_near_zero_singular_value():
    """A rank claim needs its conditioning stated. 379× is a real spread."""
    result = fi.measure_the_mass_preserving_flavor_rank()
    sv = result["singular_values"]
    assert len(sv) == 4
    assert min(sv) > 1e-3
    assert result["singular_value_spread"] < 1e4


def test_an_arbitrary_flavor_displacement_is_reachable_at_fixed_masses():
    """Direct construction, not just a rank count."""
    result = fi.measure_the_mass_preserving_flavor_rank()
    assert len(result["direct_construction"]) == 3
    for trial in result["direct_construction"]:
        assert trial["flavor_miss"] < 1e-10
        assert trial["mass_disturbance"] < 1e-10


def test_there_is_no_left_null_vector_so_no_predicted_relation():
    """rank 4 ⟹ no `w` with `wᵀK = 0` ⟹ no invariant `wᵀδy_F = 0`."""
    result = fi.measure_the_mass_preserving_flavor_rank()
    assert result["left_null_vectors"] == 0
    K, _ = fi.flavor_response_on_mass_preserving_directions()
    left_null = fi.null_space(K.T, rcond=1e-8)
    assert left_null.shape[1] == 0


# ── the phi_h A/B test ──────────────────────────────────────────────────────

def test_fixing_the_derived_phi_h_does_not_lower_the_flavor_rank():
    """The A/B the library's own framing demands. Both cases give rank 4."""
    ab = fi.measure_the_phi_h_ab_test()
    assert ab["per_case"]["phi_h_released"]["rank_K"] == 4
    assert ab["per_case"]["phi_h_fixed"]["rank_K"] == 4
    assert not ab["fixing_phi_h_lowers_the_rank"]
    assert ab["per_case"]["phi_h_fixed"]["n_coordinates"] == 17
    assert ab["per_case"]["phi_h_released"]["n_coordinates"] == 18


def test_phi_h_is_the_most_efficient_cp_handle_even_though_not_identifying():
    """Efficiency and identifiability are different, and both get reported."""
    ab = fi.measure_the_phi_h_ab_test()
    assert ab["leading_singular_value_ratio"] > 4.0
    assert "efficiency is not identifiability" in ab["verdict"]


def test_phi_h_is_at_its_derived_value_in_the_lock():
    assert LOCKED_QUARK_PARAMS_V4.phi_h == pytest.approx(math.pi / 5.0, abs=1e-15)


# ── the census ──────────────────────────────────────────────────────────────

def test_the_trace_of_each_diagonal_shift_triple_is_an_exact_ckm_gauge():
    """It moves masses and no flavor observable — so it is a mass-only knob."""
    census = fi.measure_the_calibration_dof_census()
    assert census["the_trace_direction_is_flavor_blind"]
    for name, v in census["diagonal_shift_traces"].items():
        assert v["is_an_exact_ckm_gauge"]
        assert v["flavor_response_norm"] < 1e-8
        assert v["mass_response_norm"] > 1.0


def test_both_realised_diagonal_shift_triples_are_traceless():
    census = fi.measure_the_calibration_dof_census()
    assert census["both_realised_triples_are_traceless"]
    for v in census["diagonal_shift_traces"].values():
        assert abs(v["realised_trace"]) < 1e-9


def test_symbol_count_is_not_calibration_dof_count():
    """Three symbols per triple, measured flavor rank two."""
    census = fi.measure_the_calibration_dof_census()
    rows = {r["group"]: r for r in census["rows"]}
    for name in ("diag_shift_plus", "diag_shift_minus"):
        assert rows[name]["symbols"] == 3
        assert rows[name]["measured_flavor_rank"] == 2
    assert rows["v4 targeted etas"]["measured_flavor_rank"] == 3
    assert rows["phi_h"]["measured_flavor_rank"] == 1


def test_the_v4_additions_saturate_the_flavor_space_with_phi_h_fixed():
    census = fi.measure_the_calibration_dof_census()
    rows = {r["group"]: r for r in census["rows"]}
    fixed = rows["v4 additions, phi_h fixed"]
    assert fixed["symbols"] == 9
    assert fixed["measured_flavor_rank"] == fi.PHYSICAL_FLAVOR_DIMENSION


# ── scoping: the negative verdicts must not soften ──────────────────────────

def test_the_counting_claim_is_refuted_on_the_ceiling_alone():
    """'+5 independent observables' exceeds what a unitary CKM can supply."""
    cnt = fi.measure_the_counting_claim()
    assert cnt["claimed_new_independent_observables"] == 5
    assert cnt["the_observable_claim_exceeds_the_ceiling"]
    assert cnt["physical_flavor_dimension"] == 4
    assert cnt["measured_net_surplus_phi_h_fixed"] <= 0
    assert cnt["claimed_net_surplus"] == 2


def test_the_ledger_withdraws_the_prediction_claim():
    entries = {e["claim"]: e for e in fi.measure_the_flavor_ledger()["entries"]}
    prediction = next(e for k, e in entries.items() if "is a prediction" in k)
    assert prediction["verdict"] == "WITHDRAWN"
    relation = next(e for k, e in entries.items()
                    if "first-order flavor relation" in k)
    assert relation["verdict"] == "NONE EXISTS"
    counting = next(e for k, e in entries.items() if "+3 parameters" in k)
    assert counting["verdict"] == "REFUTED"


def test_the_ledger_states_what_the_result_does_not_say():
    """The <=1% agreement is real; the rank says it is not evidence FOR the model."""
    led = fi.measure_the_flavor_ledger()
    assert "not evidence FOR the Hamiltonian" in led["what_this_does_not_say"]
    assert "local" in led["scope"].lower()
    assert "first-order" in led["scope"]


def test_the_status_docs_carry_the_flavor_outcome():
    readme = _read("README.md")
    assert "flavor_identifiability" in readme
    assert "rank K = 4" in readme
    doc = _read("docs/flavor_identifiability.md")
    assert "not a prediction" in doc


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import flavor_identifiability_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
