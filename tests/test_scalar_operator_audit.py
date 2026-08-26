"""The radial scalar operator, corrected — and what the correction cost.

Three kinds of check. The *algebraic* ones (the gap, its ℓ-independence, the
flat limit) are asserted at machine precision against closed forms. The *impact*
ones are asserted as directions and magnitudes, not as frozen digits, because
their job is to characterise a change rather than to lock a number. And the
*scoping* ones pin the claims that must not silently drift back.
"""

import math

import numpy as np
import pytest

from geometrodynamics.constants import R_MID, R_OUTER
from geometrodynamics.tangherlini import operator_audit as oa
from geometrodynamics.tangherlini.radial import (V_scalar_tangherlini,
                                                 V_tangherlini,
                                                 V_tangherlini_legacy)


# ── the operator itself ─────────────────────────────────────────────────────
def test_the_canonical_name_now_carries_the_scalar_operator():
    r = np.linspace(1.02, 6.0, 200)
    for ell in (0, 1, 2, 5):
        assert np.array_equal(np.asarray(V_tangherlini(r, ell, 1.0)),
                              np.asarray(V_scalar_tangherlini(r, ell, 1.0)))


def test_the_legacy_operator_is_frozen_exactly_as_it_was():
    """It must keep returning what it always returned, digit for digit."""
    r = np.array([1.5, 2.0, 4.0])
    A = 1.0 - 1.0 / r ** 2
    for ell in (0, 1, 2, 3):
        assert np.allclose(V_tangherlini_legacy(r, ell, 1.0),
                           A * (ell * (ell + 2) / r ** 2 + 3.0 / r ** 4),
                           rtol=0, atol=0)


def test_the_gap_is_the_stated_closed_form_and_carries_no_ell():
    got = oa.measure_the_two_operators_and_their_exact_gap()
    assert got["gap_matches_the_closed_form"] < 1e-12
    assert got["the_gap_carries_no_ell"] < 1e-12
    r = np.linspace(1.05, 8.0, 300)
    A = 1.0 - 1.0 / r ** 2
    for ell in (0, 3, 7):
        gap = (np.asarray(V_scalar_tangherlini(r, ell, 1.0))
               - np.asarray(V_tangherlini_legacy(r, ell, 1.0)))
        assert np.allclose(gap, 3.0 * A ** 2 / (4.0 * r ** 2), atol=1e-14)


def test_the_flat_limit_is_bessel_which_is_what_settles_it():
    """`ψ = r^{1/2}J_{ℓ+1}(ωr)` ⟹ `V = ((ℓ+1)² − ¼)/r²`, with no `r_h` anywhere."""
    got = oa.measure_the_two_operators_and_their_exact_gap()
    assert got["flat_limit_matches_bessel"] < 1e-12
    r = np.linspace(0.6, 9.0, 200)
    for ell in (0, 1, 2, 4, 6):
        assert np.allclose(V_scalar_tangherlini(r, ell, 1e-9),
                           ((ell + 1) ** 2 - 0.25) / r ** 2, rtol=1e-7)
    # and the legacy operator does NOT reproduce it — the discriminator bites
    assert not np.allclose(V_tangherlini_legacy(r, 2, 1e-9),
                           ((2 + 1) ** 2 - 0.25) / r ** 2, rtol=1e-3)


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_the_general_n_form_matches_its_own_definition(n):
    """`V = A[ℓ(ℓ+n−1)/r² + n(n−2)A/(4r²) + nA'/(2r)]`, spelled out."""
    r = np.linspace(1.3, 7.0, 120)
    rh, ell = 1.0, 2
    A = 1.0 - (rh / r) ** (n - 1)
    Ap = (n - 1) * rh ** (n - 1) / r ** n
    want = A * (ell * (ell + n - 1) / r ** 2
                + (n * (n - 2) / 4.0) * A / r ** 2 + (n / 2.0) * Ap / r)
    assert np.allclose(V_scalar_tangherlini(r, ell, rh, n=n), want, atol=1e-15)


def test_it_agrees_bitwise_with_the_dynamics_derivation():
    """Two modules, one operator — the #270 derivation and this one."""
    from geometrodynamics.tangherlini.dynamics import master_potential
    r = np.linspace(1.02, 6.0, 400)
    for ell in (0, 1, 2, 3, 5):
        assert np.array_equal(np.asarray(V_scalar_tangherlini(r, ell, 1.0)),
                              np.asarray(master_potential(r, ell, 1.0)))


# ── what survives exactly ───────────────────────────────────────────────────
def test_the_cross_ell_operator_is_algebraically_unchanged():
    """`ΔV` has no `ℓ`, so `V_{ℓ₂} − V_{ℓ₁}` cannot move. The load-bearing one."""
    got = oa.measure_what_survives_exactly()
    assert got["the_cross_ell_operator_is_unchanged"] < 1e-12
    r = np.linspace(1.02, 5.0, 300)
    for l1 in range(0, 4):
        for l2 in range(l1 + 1, 6):
            d_old = (np.asarray(V_tangherlini_legacy(r, l2, 1.0))
                     - np.asarray(V_tangherlini_legacy(r, l1, 1.0)))
            d_new = (np.asarray(V_scalar_tangherlini(r, l2, 1.0))
                     - np.asarray(V_scalar_tangherlini(r, l1, 1.0)))
            assert np.allclose(d_new, d_old, atol=1e-13)


def test_the_matrix_elements_of_that_exact_operator_still_drift():
    """Structure invariant, numbers shifted — the distinction the audit turns on."""
    got = oa.measure_what_survives_exactly()
    assert got["structure_invariant_numbers_shifted"]
    assert got["largest_element_drift_percent"] > 0.0, \
        "if the elements did not move, the eigenfunctions would not have"
    assert got["largest_element_drift_percent"] < 50.0


# ── the impact ──────────────────────────────────────────────────────────────
def test_the_eigenvalues_move_only_at_the_tenth_of_a_percent_level():
    got = oa.measure_the_eigenvalue_shifts()
    assert got["all_shifts_below_a_fifth_of_a_percent"]
    assert got["eigenvectors_barely_move"]
    assert got["omega_1_0_correct"] > got["omega_1_0_legacy"], "the shift is up"


def test_the_eigenvalue_sensitivity_falls_monotonically_with_ell():
    """An `ℓ`-independent shift matters least where the centrifugal term wins."""
    got = oa.measure_the_eigenvalue_shifts()
    assert got["sensitivity_falls_with_ell"]
    shifts = [abs(r["ground_shift_percent"]) for r in got["rows"]]
    assert shifts[0] > shifts[-1] * 3.0, "the fall should be substantial"


def test_the_barrier_sums_move_more_than_the_eigenvalues():
    """A barrier reads the potential; an eigenvalue averages it."""
    eig = oa.measure_the_eigenvalue_shifts()
    gam = oa.measure_the_gamma_sums_and_the_r_outer_fixed_point()
    worst_sum = max(abs(100.0 * (r["sum_correct"] - r["sum_legacy"])
                        / r["sum_legacy"]) for r in gam["rows"])
    assert worst_sum > eig["largest_ground_shift_percent"]


def test_the_canonical_gamma_claim_improves_and_nothing_was_tuned():
    """`Σ V_max[1..5]` vs the locked `22.5`: `−2.2%` → `−0.75%`."""
    got = oa.measure_the_gamma_sums_and_the_r_outer_fixed_point()
    assert got["the_canonical_readme_claim_improves"]
    assert abs(got["canonical_residual_before"]) > 2.0
    assert abs(got["canonical_residual_after"]) < 1.0
    assert got["canonical_residual_after"] < 0.0, "still short of 22.5, not over"


def test_the_ell_zero_closure_claim_breaks_and_the_channel_set_swaps():
    """The one INTERPRETATION CHANGED entry, asserted as a reversal."""
    got = oa.measure_the_gamma_sums_and_the_r_outer_fixed_point()
    assert got["the_ell_zero_closure_claim_breaks"]
    assert got["ell_zero_residual_before"] < 0.0, "used to undershoot"
    assert got["ell_zero_residual_after"] > 0.0, "now overshoots"
    assert got["the_closest_channel_set_swaps"]
    assert got["closest_to_gamma_before"] == "l = 0..5"
    assert got["closest_to_gamma_after"] == "l = 1..5"


def test_the_r_outer_fixed_point_moves_for_both_channel_sets():
    got = oa.measure_the_gamma_sums_and_the_r_outer_fixed_point()
    for row in got["rows"]:
        assert row["r_outer_legacy"] is not None
        assert row["r_outer_correct"] is not None
        assert row["r_outer_correct"] < row["r_outer_legacy"], \
            "a taller barrier reaches gamma at a smaller radius"


def test_the_actions_sit_between_the_eigenvalues_and_the_barriers():
    eig = oa.measure_the_eigenvalue_shifts()
    act = oa.measure_the_wkb_action_shift()
    assert act["each_action_uses_its_own_ground_frequency"]
    assert act["largest_drift_percent"] > eig["largest_ground_shift_percent"]


def test_the_throat_flux_ratios_mostly_cancel_the_shift():
    got = oa.measure_the_eigenvector_derived_quantities()
    assert got["the_reference_mode_is_exactly_one"]
    assert got["ratios_absorb_most_of_the_shift"]
    assert "not automatically safe" in got["caveat"]


# ── the ledger ──────────────────────────────────────────────────────────────
def test_every_load_bearing_claim_carries_one_of_three_verdicts():
    got = oa.measure_the_downstream_ledger()
    allowed = {"EXACTLY INVARIANT", "NUMERICALLY SHIFTED", "INTERPRETATION CHANGED"}
    assert all(e["verdict"] in allowed for e in got["entries"])
    assert sum(got["counts"].values()) == len(got["entries"])
    for v in allowed:
        assert got["counts"][v] >= 1, f"no entry landed in {v}"


def test_the_topological_results_are_named_as_not_re_run():
    """Proximity is not dependence, and the audit says so out loud."""
    got = oa.measure_the_downstream_ledger()
    text = " ".join(e["claim"] + e["evidence"] for e in got["entries"])
    for name in ("Hopf", "Pin-", "odd-k", "antipodal parity"):
        assert name in text
    assert "proximity is not dependence" in got["not_re_run_and_why"]


def test_the_gamma_narrative_is_withdrawn_not_replaced():
    """The corrected sum landing nearer 22.5 is an observation, not a derivation."""
    got = oa.measure_the_downstream_ledger()
    assert "withdrawn, not replaced" in got["what_is_still_open"]
    assert "not a derivation of why" in got["what_is_still_open"]


def test_the_1054_factor_is_flagged_against_its_own_tolerance():
    got = oa.measure_the_downstream_ledger()
    entry = next(e for e in got["entries"] if "1.054" in e["claim"])
    assert entry["verdict"] == "INTERPRETATION CHANGED", \
        "there is no unique gamma-locked R_OUTER left to evaluate it at"
    fixed = next(e for e in got["entries"] if "FIXED R_OUTER" in e["claim"])
    assert fixed["verdict"] == "NUMERICALLY SHIFTED"
    assert "withdrawn" in entry["evidence"]


# ── the audit machinery itself ──────────────────────────────────────────────
def test_the_audit_grid_matches_the_production_solver():
    """The audit must measure the same operator the repository runs."""
    from geometrodynamics.tangherlini.radial import solve_radial_modes
    for ell in (0, 2):
        prod, _, rg_prod = solve_radial_modes(ell, N=80, n_modes=2)
        mine, _, rg_mine = oa.eigen_solve(ell, V_scalar_tangherlini,
                                          points=80, n_modes=2)
        assert np.allclose(rg_prod, rg_mine, atol=1e-12)
        assert np.allclose(prod[:2], mine[:2], rtol=1e-10)


def test_both_operators_are_registered_for_the_sweep():
    assert set(oa.OPERATORS) == {"legacy", "scalar_correct"}
    r = np.linspace(1.3, 4.0, 50)
    assert not np.allclose(oa.OPERATORS["legacy"](r, 2, 1.0),
                           oa.OPERATORS["scalar_correct"](r, 2, 1.0))


def test_the_locked_gamma_is_the_repository_constant():
    from geometrodynamics.tangherlini import LEPTON_BASELINE_PINHOLE
    assert oa.LOCKED_GAMMA == LEPTON_BASELINE_PINHOLE == 22.5


# ── the one narrow downstream re-derivation ─────────────────────────────────
def test_the_locked_lepton_block_is_reproduced_before_anything_is_varied():
    """If the baseline row is not the locked spectrum, the harness is wrong."""
    got = oa.measure_which_geometry_preserves_the_lepton_ladder()
    base = next(r for r in got["rows"] if r["case"].startswith("baseline"))
    assert base["gamma"] == 22.5
    assert abs(base["mu_error_percent"]) < 0.5
    assert abs(base["tau_error_percent"]) < 0.5
    # and the legacy fixed-R row reproduces the documented "within 3.8%"
    leg = next(r for r in got["rows"] if r["case"] == "legacy R=1.26, gamma[0..5]")
    assert 3.0 < leg["mu_error_percent"] < 4.5


def test_B_and_C_are_bit_identical_because_r_outer_is_not_an_input():
    """The locked block discards `r_outer`; only `gamma` reaches it."""
    got = oa.measure_which_geometry_preserves_the_lepton_ladder()
    b = next(r for r in got["rows"] if r["case"].startswith("B "))
    c = next(r for r in got["rows"] if r["case"].startswith("C "))
    assert b["r_outer"] != c["r_outer"], "the two geometries really do differ"
    assert b["m_mu"] == c["m_mu"]
    assert b["m_tau"] == c["m_tau"]
    assert got["B_and_C_are_bit_identical"]
    assert got["so_the_channel_set_is_invisible_to_the_observables"]


def test_gamma_is_the_selector_and_fixing_r_outer_breaks_the_ladder():
    """The outcome inverts the anticipated one: R_OUTER is downstream of gamma."""
    got = oa.measure_which_geometry_preserves_the_lepton_ladder()
    assert got["gamma_is_the_selector_r_outer_is_downstream"]
    assert got["fixing_r_outer_breaks_the_ladder"]
    errs = got["corrected_fixed_R_mu_errors"]
    assert errs[0] > 10.0 and errs[1] < -10.0, "one overshoots, one undershoots"
    assert got["the_outcome_was_not_one_of_the_three_anticipated"]


def test_the_correction_weakens_the_geometry_supplies_gamma_story():
    """Legacy R=1.26 landed within 3.8%; corrected lands at 15-21%, either way."""
    got = oa.measure_which_geometry_preserves_the_lepton_ladder()
    assert got["the_correction_weakens_the_geometry_supplies_gamma_story"]
    assert abs(got["legacy_fixed_R_mu_error"]) < 5.0
    assert all(abs(e) > 10.0 for e in got["corrected_fixed_R_mu_errors"])


def test_a_subpercent_gamma_residual_is_not_a_small_residual():
    """`d ln m_mu / d ln gamma` is large and negative, so gamma carries weight."""
    got = oa.measure_which_geometry_preserves_the_lepton_ladder()
    assert got["local_d_ln_m_mu_over_d_ln_gamma_at_22p5"] < -10.0
    assert got["secant_elasticity_over_22p331_to_22p836"] < -10.0
    assert got["so_a_subpercent_residual_is_not_small"]


def test_the_channel_set_question_is_left_undecided_on_purpose():
    got = oa.measure_which_geometry_preserves_the_lepton_ladder()
    assert got["nothing_was_retuned"]
    assert "not decidable by the lepton observables" in got["what_this_does_not_settle"]
    assert "only ever saw the scalar" in got["what_this_does_not_settle"]


# ── the review's fixes, pinned so they cannot drift back ────────────────────
def test_the_transport_element_is_the_off_diagonal_historical_object():
    """`⟨u_ℓ₁|V_ℓ₂−V_ℓ₁|u_ℓ₂⟩`, not a diagonal expectation.

    An earlier version computed `dot(v_ℓ², ΔV)` — no `u_{ℓ+2}` in it at all —
    and quoted the result as the drift of the *cross-ℓ transport* elements.
    That is a different object. This pins the real one.
    """
    got = oa.measure_what_survives_exactly()
    assert "off-diagonal" in got["the_element_measured"]
    assert "u_l2" in got["the_element_measured"]
    pairs = [e["pair"] for e in got["matrix_elements"]]
    assert pairs == ["<u_1|V_3-V_1|u_3>", "<u_3|V_5-V_3|u_5>", "<u_1|V_5-V_1|u_5>"]
    assert "diagonal" in got["what_an_earlier_version_measured"]
    # the elements are genuinely off-diagonal: swapping the states changes nothing,
    # but dropping one of them would change everything
    for e in got["matrix_elements"]:
        assert abs(e["element_legacy"]) > 1e-3
        assert abs(e["drift_percent"]) < 1.0


def test_the_closed_orbit_action_carries_the_return_leg():
    """`S_full = 2∫√(ω²−V)dr*` — the ledger's `∮p dq`, not the one-way integral."""
    one = oa.one_way_wkb_action(1.1, V_scalar_tangherlini, 1)
    both = oa.closed_orbit_action(1.1, V_scalar_tangherlini, 1)
    assert both == pytest.approx(2.0 * one, rel=1e-12)
    assert one > 0.0
    # and the audit reports the doubled one
    got = oa.measure_the_wkb_action_shift()
    for row in got["rows"]:
        w, _, _ = oa.eigen_solve(row["ell"], V_scalar_tangherlini, n_modes=1)
        assert row["action_correct"] == pytest.approx(
            2.0 * oa.one_way_wkb_action(float(w[0]), V_scalar_tangherlini,
                                        row["ell"]), rel=1e-9)


def test_the_secant_and_the_local_derivative_are_reported_separately():
    """One is a finite difference over a range; the other is `d/d` at the lock."""
    got = oa.measure_which_geometry_preserves_the_lepton_ladder()
    secant = got["secant_elasticity_over_22p331_to_22p836"]
    local = got["local_d_ln_m_mu_over_d_ln_gamma_at_22p5"]
    assert secant < -10.0 and local < -10.0
    assert secant != local, "a secant over a finite range is not the derivative"
    assert got["the_headline_number_is_the_local_derivative"]
    assert "finite difference" in got["why_both_are_reported"]


def test_the_gamma_locked_1054_factor_is_withdrawn_not_re_quoted():
    """No unique γ-locked `R_OUTER` survives, so the factor cannot be re-evaluated."""
    got = oa.measure_the_downstream_ledger()
    locked = next(e for e in got["entries"] if "GAMMA-LOCKED" in e["claim"])
    assert locked["verdict"] == "INTERPRETATION CHANGED"
    assert "1.24614" in locked["evidence"] and "1.26788" in locked["evidence"]
    fixed = next(e for e in got["entries"] if "FIXED R_OUTER" in e["claim"])
    assert fixed["verdict"] == "NUMERICALLY SHIFTED"
    assert fixed["claim"] != locked["claim"], "two different quantities"


def test_the_status_docs_no_longer_assert_the_reopened_claims():
    """README and THESIS must not still call the reopened claims Verified."""
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[1]
    readme = (root / "README.md").read_text()
    for row in ("R_OUTER selected by cross-species fixed point",
                "Pinhole γ ≈ Σ V_max[1..5] on Chebyshev grid"):
        line = next(l for l in readme.splitlines() if row in l)
        assert "**Verified**" not in line, f"still asserted as Verified: {row}"
        assert "REOPENED" in line
    thesis = (root / "docs" / "THESIS.md").read_text()
    assert "REOPENED BY PR #271" in thesis
    # the present-tense claim must be gone; the past-tense record may stay
    assert "is the\n  unique physical selection" not in thesis
    assert "was the\n  unique physical selection" in thesis
    # and the reopened marker has to actually follow the claim it reopens
    head, tail = thesis.split("REOPENED BY PR #271", 1)
    assert "cross-species self-consistency loop" in head[-1200:]
    assert "bit-identical" in tail[:1200]
