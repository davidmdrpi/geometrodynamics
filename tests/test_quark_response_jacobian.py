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


def test_the_phase_symmetry_is_antiunitary_and_model_wide_not_lock_only():
    """`H(−φ, p) = H(φ, p)*` for ARBITRARY p, so the SPECTRUM is always even.

    The first draft asserted only `max|H(+φ) − H(−φ)| > 0` at nonzero mixing
    and concluded the Z₂ was a property of the lock. That mistook matrix
    inequality for spectral asymmetry. This test pins the real statement: the
    matrix does change, and the eigenvalues and anchored masses do not.
    """
    for p in (0.0, 0.05, 0.3):
        for phi in (0.1, 0.37):
            Hp = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                                 phase=phi,
                                                 partition_mixing=p))
            Hm = build_quark_hamiltonian(replace(LOCKED_QUARK_PARAMS,
                                                 phase=-phi,
                                                 partition_mixing=p))
            # the antiunitary identity holds at every mixing
            assert np.max(np.abs(Hm - np.conj(Hp))) == 0.0
            # hence isospectral, hence the masses are even in phase
            assert np.max(np.abs(np.linalg.eigvalsh(Hp)
                                 - np.linalg.eigvalsh(Hm))) == 0.0
            if p != 0.0:
                # the MATRIX does differ -- that is what misled the first draft
                assert np.max(np.abs(Hp - Hm)) > 1e-3
                assert np.max(np.abs(
                    rj.anchored_log_masses(replace(LOCKED_QUARK_PARAMS,
                                                   phase=phi,
                                                   partition_mixing=p))
                    - rj.anchored_log_masses(replace(LOCKED_QUARK_PARAMS,
                                                     phase=-phi,
                                                     partition_mixing=p)))) == 0.0


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
    assert "ANTIUNITARY" in n["phase"]["status"]
    assert n["phase"]["the_symmetry_is_model_wide_not_lock_only"]
    assert "UNITARY-CONJUGATION" in n["partition_mixing"]["status"]
    assert n["phase"]["status"] != n["partition_mixing"]["status"]
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
    assert svd["n_first_order_null_directions"] == 4
    assert svd["gram_exact_zeros"] == 4


def test_the_pathology_is_rank_deficiency_not_ill_conditioning():
    """Guards the wording: no long hierarchy among the NONZERO singular values.

    "Condition number 22.6" is only the ratio over the four nonzero values;
    `JᵀJ` has four exact zeros. The claim must be about rank, not conditioning.
    """
    svd = rj.measure_the_singular_system_and_effective_rank()
    assert svd["no_long_hierarchy_among_nonzero_singular_values"]
    assert 10.0 < svd["nonzero_singular_value_ratio"] < 100.0
    assert svd["gram_exact_zeros"] == 4
    assert "rank, not" in svd["the_dominant_pathology_is_structural_non_identifiability"]


def test_which_knobs_drift_is_itself_structure():
    """"Any knob would drift" is false: the null content varies hugely."""
    svd = rj.measure_the_singular_system_and_effective_rank()
    share = svd["identifiable_share_by_knob"]
    assert share["uplift_asymmetry"] == pytest.approx(1.0, abs=1e-3)
    assert share["eta_k3k5_minus"] < 0.01
    assert svd["most_identifiable_knob"] == "uplift_asymmetry"
    assert svd["least_identifiable_knob"] == "eta_k3k5_minus"


def test_every_right_space_claim_declares_its_metric():
    """The min-norm geometry is chart-dependent and must say so."""
    svd = rj.measure_the_singular_system_and_effective_rank()
    assert "log coordinates" in svd["metric_scope"]
    rep = rj.measure_which_directions_repair_the_masses()
    assert "rescaling" in rep["metric_caveat"]


def test_the_quadratic_directions_add_no_new_observable_direction():
    """And the projection residual being ~0 is automatic, not a discovery."""
    svd = rj.measure_the_singular_system_and_effective_rank()
    equiv = svd["quadratic_directions_add_no_observable_direction"]
    assert set(equiv) == set(rj.QUADRATIC_KNOBS)
    for v in equiv.values():
        assert v["projection_residual"] < 1e-10
    note = svd["the_projection_residual_is_automatic"]
    assert "column space is all of R^4" in note and "ANY 4-vector" in note
    # and it must NOT be used to justify dropping the quadratic coordinates
    assert "does NOT justify" in note


def test_the_minimum_norm_repair_reproduces_the_residual_exactly():
    J, _ = rj.response_jacobian()
    delta, c, S = rj.minimum_norm_correction()
    r = rj.baseline_log_residual()
    assert np.linalg.norm(J @ delta - r) < 1e-12
    assert np.linalg.norm(delta) == pytest.approx(0.0229, abs=5e-4)


def test_the_exact_nonlinear_repair_is_computed_not_quoted():
    """The 0.018% must come from a re-run, not from prose.

    #271/#272's central lesson is that an untested scientific number survives a
    green suite. The first draft returned this one as a literal.
    """
    rep = rj.measure_which_directions_repair_the_masses()
    assert rep["exact_max_rel_err_percent_after_repair"] < 0.05
    assert rep["locked_max_rel_err_percent"] == pytest.approx(1.61, abs=0.05)
    # recompute independently and require agreement
    delta, _, _ = rj.minimum_norm_correction()
    moved = replace(LOCKED_QUARK_PARAMS,
                    **{k: getattr(LOCKED_QUARK_PARAMS, k) * math.exp(v)
                       for k, v in zip(rj.LOG_KNOBS, delta)})
    observed = rj.anchored_log_masses() + rj.baseline_log_residual()
    exact = 100.0 * np.max(np.abs(np.expm1(
        rj.anchored_log_masses(moved) - observed)))
    assert rep["exact_max_rel_err_percent_after_repair"] == pytest.approx(
        exact, rel=1e-9)


def test_the_largest_knob_move_is_under_two_percent_not_one():
    """`beta` moves +1.78%, so "~1% or less" was wrong."""
    rep = rj.measure_which_directions_repair_the_masses()
    assert 1.0 < rep["largest_single_knob_move_percent"] < 2.0


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

def test_the_operator_induced_move_points_toward_the_residual():
    """The corrected headline: the move is `Δg`, not `g_corrected`.

    The first draft compared `J·g_corrected` to `r` and read `cos = −0.616` as
    "the correction moves away from the data". The move the correction makes is
    `Δg = g_corrected − g_legacy`, against the residual the legacy triple
    leaves, and that cosine is strongly POSITIVE.
    """
    prj = rj.measure_the_geometric_displacement_against_what_the_masses_want()
    move = prj["the_operator_induced_move"]
    assert move["cos_move_vs_residual_after_legacy"] == pytest.approx(
        0.873, abs=0.01)
    assert move["cos_move_vs_residual_after_legacy"] > 0.8
    # the candidate-from-the-lock cosines still have their old signs ...
    legacy = prj["per_operator"]["legacy"]
    corrected = prj["per_operator"]["scalar_correct"]
    assert legacy["cos_to_r_observable_space"] > 0.0
    assert corrected["cos_to_r_observable_space"] < 0.0
    # ... and the withdrawal must say so explicitly
    assert "Withdrawn" in prj["withdrawn"] or "withdrawn" in prj["withdrawn"]
    assert "displacement FROM" in prj["withdrawn"]
    for side in (legacy, corrected):
        assert side["fraction_of_displacement_in_the_null_space"] > 0.5


def test_the_two_objectives_disagree_and_both_are_reported():
    """L2 improves while max-rel-error worsens. Neither may be privileged."""
    prj = rj.measure_the_geometric_displacement_against_what_the_masses_want()
    obj = prj["the_two_objectives_disagree"]
    assert obj["l2_improves"] and obj["max_rel_err_worsens"]
    assert (obj["l2_log_residual"]["scalar_correct"]
            < obj["l2_log_residual"]["legacy"])
    assert (obj["max_rel_err_percent_linear"]["scalar_correct"]
            > obj["max_rel_err_percent_linear"]["legacy"])
    assert "objective-dependent" in prj["verdict"]


def test_the_holdout_measures_the_regulariser_not_the_hamiltonian():
    """Every held-out species is reachable from ker(J_keep) — so it predicts nothing.

    The first draft read the pseudoinverse's −10.4% miss as the Hamiltonian
    failing a holdout. With a 5-dimensional kernel and the held-out row outside
    the span of the other three, a kernel shift fits the withheld mass exactly.
    """
    hld = rj.measure_the_minimum_norm_regulariser_does_not_predict_a_held_out_mass()
    assert hld["every_held_out_species_is_reachable"]
    assert hld["a_kernel_shift_fits_the_held_out_mass_exactly"]
    for row in hld["rows"]:
        assert row["kernel_dimension"] == 5
        assert row["held_out_row_reachable_from_kernel"] > 1e-9
        shifted = row["after_a_kernel_shift"]
        assert abs(shifted["held_out_error_percent"]) < 1e-9
        assert shifted["fitted_three_max_error_percent"] < 1e-9
        assert shifted["norm_ratio_to_min_norm"] < 10.0
    assert "not of the Hamiltonian" in hld["verdict"]


def test_the_ledger_withdraws_the_independent_constraint_claim():
    entries = {e["claim"]: e for e in rj.measure_the_jacobian_ledger()["entries"]}
    independent = next(e for k, e in entries.items()
                       if "can identify the current quark" in k)
    assert independent["verdict"] == "WITHDRAWN"
    assert "rank J = 4" in independent["evidence"]
    conjecture = next(e for k, e in entries.items()
                      if "scalar relation" in k)
    assert conjecture["verdict"] == "REFUTED"


def test_the_ledger_records_this_rounds_own_withdrawals():
    """Three first-draft claims were wrong. The ledger has to say so."""
    entries = {e["claim"]: e for e in rj.measure_the_jacobian_ledger()["entries"]}
    own = [e for e in entries.values() if "OWN" in e["verdict"]]
    assert len(own) >= 3, "the round's own errors must be in its ledger"
    phase = next(e for k, e in entries.items() if "property of the lock" in k)
    assert "WITHDRAWN" in phase["verdict"]
    away = next(e for k, e in entries.items() if "moves AWAY" in k)
    assert "WITHDRAWN" in away["verdict"]
    holdout = next(e for k, e in entries.items() if "leave-one-species-out" in k)
    assert "WITHDRAWN" in holdout["verdict"]


def test_the_successor_claim_does_not_assume_ckm_is_column_free():
    """"Adds rows without adding columns" is not automatic — v4 has new knobs."""
    settle = rj.measure_the_jacobian_ledger()["what_would_settle_it"]
    assert "phi_h" in settle
    assert "NOT" in settle or "not" in settle


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
    assert "rank J = 4" in readme
    assert "+0.873" in readme, "the corrected headline cosine must be in README"


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import quark_response_jacobian_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
