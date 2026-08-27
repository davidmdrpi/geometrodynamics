"""The quark response Jacobian and the local geometry of its fit manifold.

Four kinds of check. The *algebraic* ones assert the null structure against the
closed forms it was derived from — an identity shift, an even function, a
unitary conjugation — at machine precision. The *numerical* ones pin the
Jacobian's convergence, because the adiabatic species relabelling could have
made it noise. The *structural* ones pin rank and identifiability. And the
*scoping* ones guard the three conclusions that must not soften back into
"the geometry roughly agrees".
"""

import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.qcd import response_jacobian as rj
from geometrodynamics.qcd.quark_spectrum import (BASIS_STATES, _SIGMA,
                                                 LOCKED_QUARK_PARAMS,
                                                 build_quark_hamiltonian)


def _read(rel: str) -> str:
    return (Path(__file__).resolve().parents[1] / rel).read_text()


# ── algebraic: the null structure, from its closed forms ────────────────────

def test_action_base_shifts_the_hamiltonian_by_the_identity():
    """`H(a) = H(0) + a·I` — the reason it is an exact gauge, not a soft direction."""
    n = len(BASIS_STATES)
    H0 = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, action_base=0.0))
    for a in (1.0, 7.3, -2.5):
        Ha = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                             action_base=a))
        assert np.max(np.abs((Ha - H0) - a * np.eye(n))) < 1e-11


def test_action_base_is_invariant_before_the_anchor_not_because_of_it():
    """The zero-point subtraction kills it, upstream of the d-anchor.

    This matters: if the d-anchor were doing the work, `action_base` would
    reappear in any absolute-scale observable. It does not.
    """
    H0 = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, action_base=0.0))
    Ha = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, action_base=7.3))
    w0, wa = np.linalg.eigvalsh(H0), np.linalg.eigvalsh(Ha)
    assert np.max(np.abs((wa - wa.min()) - (w0 - w0.min()))) < 1e-11
    assert "action_base" in rj.EXACT_GAUGE_KNOBS
    assert "action_base" not in rj.LOG_KNOBS


def test_action_base_is_flat_at_finite_displacement_not_merely_stationary():
    """A 3× change must move nothing. Stationarity would not survive this."""
    base = rj.anchored_log_masses()
    for factor in (1.05, 1.5, 3.0):
        moved = rj.anchored_log_masses(
            replace(LOCKED_QUARK_PARAMS,
                    action_base=LOCKED_QUARK_PARAMS.action_base * factor))
        assert np.max(np.abs(moved - base)) < 1e-10


def test_phase_is_even_at_the_lock_and_only_because_mixing_is_off():
    """The Z₂ is a property of the lock, not of the model."""
    for phi in (0.1, 0.37):
        Hp = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, phase=phi))
        Hm = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, phase=-phi))
        assert np.max(np.abs(Hp - Hm)) == 0.0
    # switch the different-partition element on and the evenness must break
    Hp = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, phase=0.37,
                                         partition_mixing=0.05))
    Hm = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS, phase=-0.37,
                                         partition_mixing=0.05))
    assert np.max(np.abs(Hp - Hm)) > 1e-3


def test_partition_mixing_sign_flip_is_a_unitary_conjugation():
    """`H(−p) = D H(p) D†`, `D = diag(σ)`. Isospectral, hence even in p."""
    D = np.diag([_SIGMA[p] for (_, p) in BASIS_STATES]).astype(float)
    assert np.array_equal(np.diag(D), [1, -1, 1, -1, 1, -1])
    for p in (0.05, 0.11):
        Hp = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                             partition_mixing=p))
        Hm = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                             partition_mixing=-p))
        assert np.max(np.abs(Hm - D @ Hp @ D)) == 0.0
        assert np.max(np.abs(np.linalg.eigvalsh(Hp)
                             - np.linalg.eigvalsh(Hm))) == 0.0


def test_the_quadratic_chart_is_linear_in_q_and_the_response_scales_as_x_squared():
    """A 10× step in `x` gives 100× the response — what licenses `q = x²`."""
    for knob in rj.QUADRATIC_KNOBS:
        base = rj.anchored_log_masses()
        small = np.max(np.abs(rj.anchored_log_masses(
            replace(LOCKED_QUARK_PARAMS, **{knob: 1e-3})) - base))
        large = np.max(np.abs(rj.anchored_log_masses(
            replace(LOCKED_QUARK_PARAMS, **{knob: 1e-2})) - base))
        assert large / small == pytest.approx(100.0, rel=0.02)
        # and d ln m / dq is then constant
        cols = [rj.quadratic_jacobian_column(knob, q)
                for q in (1e-6, 1e-7, 1e-8)]
        norms = [float(np.linalg.norm(c)) for c in cols]
        assert (max(norms) - min(norms)) / max(norms) < 1e-3


def test_the_quadratic_chart_refuses_to_run_away_from_the_symmetric_point():
    with pytest.raises(ValueError):
        rj.quadratic_jacobian_column(
            "phase", 1e-8, params=replace(LOCKED_QUARK_PARAMS, phase=0.2))


def test_the_three_nulls_are_reported_as_three_different_categories():
    n = rj.measure_the_null_structure_is_three_different_objects()
    assert "EXACT INVARIANCE" in n["action_base"]["status"]
    assert n["action_base"]["removed_upstream_of_the_d_anchor"]
    assert "QUADRATIC" in n["phase"]["status"]
    assert n["phase"]["the_Z2_is_a_property_of_the_lock_not_the_model"]
    assert "QUADRATIC" in n["partition_mixing"]["status"]
    assert n["partition_mixing"]["stronger_than_even_function"]
    assert n["action_base"]["status"] != n["phase"]["status"]


# ── numerical: the Jacobian is a derivative, not noise ──────────────────────

def test_the_jacobian_converges_under_step_refinement():
    """Adiabatic relabelling inside the stencil would show as step dependence."""
    result = rj.measure_the_jacobian_converges()
    assert result["all_converged_below_1e-3_relative"]
    for row in result["rows"]:
        assert row["relative_spread_over_1e-3_to_1e-6"] < 1e-3


def test_the_baseline_residual_matches_the_locked_ladder_error():
    """`r` must reproduce the v3 lock's published 1.6% max relative error."""
    r = rj.baseline_log_residual()
    assert np.max(np.abs(np.expm1(r))) == pytest.approx(0.0161, abs=5e-4)


# ── structural: rank, identifiability, and what repairs the masses ──────────

def test_the_rank_is_capped_by_the_number_of_scored_masses():
    """Four observables cannot constrain more than four combinations."""
    svd = rj.measure_the_singular_system_and_effective_rank()
    assert svd["rank"] == len(rj.SCORED_SPECIES) == 4
    assert svd["rank_is_capped_by_the_observable_count"]
    assert svd["n_first_order_knobs"] == 8
    assert svd["n_invisible_first_order_directions"] == 4


def test_it_is_not_a_sloppy_model_the_degeneracy_is_dimensional():
    """Guards against importing the wrong diagnosis: cond ~23, not 10³–10⁶."""
    svd = rj.measure_the_singular_system_and_effective_rank()
    assert svd["not_a_sloppy_model"]
    assert 10.0 < svd["condition_number"] < 100.0


def test_the_quadratic_directions_add_no_new_observable_direction():
    """And the projection residual being ~0 is automatic, not a discovery."""
    svd = rj.measure_the_singular_system_and_effective_rank()
    equiv = svd["quadratic_directions_are_not_new_observable_directions"]
    assert set(equiv) == set(rj.QUADRATIC_KNOBS)
    for v in equiv.values():
        assert v["projection_residual"] < 1e-10
    note = svd["the_projection_residual_is_automatic"]
    assert "column space is all of R^4" in note and "ANY 4-vector" in note


def test_the_minimum_norm_repair_reproduces_the_residual_exactly():
    J, _ = rj.response_jacobian()
    delta, c, S = rj.minimum_norm_correction()
    r = rj.baseline_log_residual()
    assert np.linalg.norm(J @ delta - r) < 1e-12
    assert np.linalg.norm(delta) == pytest.approx(0.0229, abs=5e-4)


def test_the_repair_is_carried_by_the_weakest_singular_direction():
    """A compensator, not a coherent physical combination."""
    rep = rj.measure_which_directions_repair_the_masses()
    assert rep["the_repair_rides_the_weakest_direction"]
    assert rep["dominant_direction"] == 4
    assert rep["dominant_share"] > 0.9
    shares = [row["share_of_correction_norm_squared"] for row in rep["rows"]]
    assert sum(shares) == pytest.approx(1.0)
    assert shares[0] < 0.01, "the stiffest direction contributes almost nothing"


def test_the_stiffest_direction_is_uplift_asymmetry_acting_on_the_bottom_quark():
    """ε = 1 − 1/k₅² is the most load-bearing number in the ladder.

    PR #272 classified it EXACTLY INVARIANT — true, it reads no radial
    operator. Being safe from the operator correction is not the same as being
    unimportant, and its elasticity had never been measured.
    """
    J, knobs = rj.response_jacobian()
    i, j = np.unravel_index(int(np.argmax(np.abs(J))), J.shape)
    assert knobs[j] == "uplift_asymmetry"
    assert rj.SCORED_SPECIES[i] == "b"
    assert abs(J[i, j]) > 15.0


# ── scoping: the conclusions that must not soften ───────────────────────────

def test_the_corrected_geometry_moves_away_from_what_the_masses_want():
    """The decisive test, and the round's headline. cos Θ must be NEGATIVE."""
    prj = rj.measure_the_geometric_displacement_against_what_the_masses_want()
    legacy = prj["per_operator"]["legacy"]
    corrected = prj["per_operator"]["scalar_correct"]
    assert legacy["cos_theta_observable_space"] > 0.0
    assert corrected["cos_theta_observable_space"] < 0.0
    assert prj["the_correction_flips_the_sign_of_the_overlap"]
    # and the same sign flip in parameter space
    assert legacy["cos_theta_parameter_space"] > 0.0
    assert corrected["cos_theta_parameter_space"] < 0.0
    # most of each displacement moves nothing at all
    for side in (legacy, corrected):
        assert side["fraction_of_displacement_invisible_to_the_masses"] > 0.5


def test_leave_one_species_out_explodes_while_the_fitted_three_do_not():
    """The holdout that separates local flexibility from structure."""
    hld = rj.measure_leave_one_species_out()
    assert hld["it_explodes_under_holdout"]
    assert hld["worst_held_out_error_percent"] > 5.0
    for row in hld["rows"]:
        assert row["fitted_species_max_error_percent_exact"] < 0.1
        assert abs(row["held_out_error_percent_exact"]) > 1.0
        # the linear prediction must track the exact re-run
        assert row["held_out_error_percent_linear"] == pytest.approx(
            row["held_out_error_percent_exact"], abs=1.0)
    worst = max(hld["rows"], key=lambda r: abs(r["held_out_error_percent_exact"]))
    assert worst["held_out"] == "b"


def test_the_ledger_withdraws_the_independent_constraint_claim():
    entries = {e["claim"]: e for e in rj.measure_the_jacobian_ledger()["entries"]}
    independent = next(e for k, e in entries.items()
                       if "independently constrained" in k)
    assert independent["verdict"] == "WITHDRAWN"
    refuted = next(e for k, e in entries.items()
                   if "per-knob improvements moved the ladder" in k)
    assert refuted["verdict"] == "REFUTED"
    conjecture = next(e for k, e in entries.items()
                      if "scalar relation" in k)
    assert conjecture["verdict"] == "REFUTED"


def test_the_previous_rounds_conjecture_is_recorded_as_wrong():
    """PR #272 guessed transport↔resistance. The nearest pair is gamma_q↔transport."""
    conjecture = next(e for e in rj.measure_the_jacobian_ledger()["entries"]
                      if "scalar relation" in e["claim"])
    assert "gamma_q/transport" in conjecture["evidence"]
    doc = _read("docs/quark_response_jacobian.md")
    assert "was wrong" in doc


def test_the_status_docs_carry_the_jacobian_outcome():
    readme = _read("README.md")
    assert "quark_response_jacobian" in readme
    for fragment in ("rank J = 4", "−0.616"):
        assert fragment in readme, f"missing {fragment}"


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import quark_response_jacobian_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
