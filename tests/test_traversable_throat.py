"""The whole-throat S-matrix of a supported traversable 5D throat.

Four kinds of check. The *exact* ones pin the closed forms the geometry
supplies — the potential, its throat value, its tail, and the stress tensor.
The *structural* ones pin what the solver must not assume: unitarity imposed
nowhere, reciprocity from the evenness of `V`, and no Fabry–Pérot split. The
*oracle* one pins the threshold law, which is deliberately **not** `T₀(0) = g`.
And the *verdict* ones guard what this module alone can say — the NEC price,
the dispersion that denies `closure_offset` a constant `τ_th`, and the absence
of any perfect-transmission point in the scanned band. Two of those are
deliberately **weaker** than an earlier draft: the raw `Δ_phase`/`Δ_group` gap
is branch dependent, and `V > 0` does not forbid `|T| = 1`. The branch-free
closure verdict is `tests/test_derived_network.py`, which needs `network.py`'s
own convention to state.
"""

import math

import numpy as np
import pytest

from geometrodynamics.tangherlini import traversable_throat as tt


# ── exact: the geometry ─────────────────────────────────────────────────────

def test_the_profile_solves_the_scalar_flat_condition():
    """`f = √(s²+a²)` satisfies `f′² = 1 − a²/f²`."""
    s = np.array([0.0, 0.5, 1.0, 3.0, 10.0])
    a = 1.0
    f = tt.throat_radius(s, a)
    derivative = s / f
    assert np.allclose(derivative ** 2, 1.0 - a * a / f ** 2, atol=1e-14)


def test_the_potential_at_the_throat_is_exact():
    """`V_ℓ(0) = [ℓ(ℓ+2) + 3/2]/a²`."""
    for a in (0.5, 1.0, 2.0):
        for ell in (0, 1, 2, 5):
            assert float(tt.throat_potential(0.0, ell, a)) == pytest.approx(
                (ell * (ell + 2) + 1.5) / (a * a), rel=1e-14)


def test_the_potential_tail_matches_pr_275s_centrifugal_form():
    """`V → [(ℓ+1)² − ¼]/s²` — the same tail, so the same Jost condition."""
    for ell in (0, 1, 2, 3):
        limit = (ell + 1) ** 2 - 0.25
        assert float(tt.throat_potential(1e5, ell, 1.0)) * 1e10 == \
            pytest.approx(limit, rel=1e-8)


def test_the_potential_is_a_single_even_positive_barrier():
    """One smooth symmetric maximum — no cavity, so no Fabry–Pérot split."""
    s = np.linspace(-30.0, 30.0, 3001)
    values = tt.throat_potential(s, 1, 1.0)
    assert np.allclose(values, values[::-1], atol=1e-13), "V must be even"
    assert np.all(values > 0.0), "V must be positive"
    assert np.argmax(values) == len(values) // 2, "the maximum is at s = 0"
    # strictly decreasing away from the throat: one barrier, not two
    half = values[len(values) // 2:]
    assert np.all(np.diff(half) < 0.0)


# ── exact: the GR price ─────────────────────────────────────────────────────

def test_the_density_vanishes_identically():
    stress = tt.stress_tensor(np.linspace(-5.0, 5.0, 101), 1.0)
    assert np.all(stress["density"] == 0.0)


def test_the_pressures_are_the_derived_closed_forms():
    """`p_s = −3a²/f⁴`, `p_Ω = +a²/f⁴` (times `1/8πG₅`)."""
    s = np.array([0.0, 0.7, 2.0, 6.0])
    a = 1.3
    quartic = (s * s + a * a) ** 2
    stress = tt.stress_tensor(s, a)
    assert np.allclose(stress["radial_pressure"], -3.0 * a * a / quartic)
    assert np.allclose(stress["angular_pressure"], a * a / quartic)


def test_the_radial_null_energy_condition_is_violated_everywhere():
    stress = tt.stress_tensor(np.linspace(-20.0, 20.0, 401), 1.0)
    assert np.all(stress["radial_nec"] < 0.0)


def test_the_null_energy_integral_is_exactly_minus_three_sixteenths_over_a():
    for a in (0.5, 1.0, 3.0):
        assert tt.null_energy_integral(a) == pytest.approx(
            -3.0 / (16.0 * a), rel=1e-15)


def test_the_exoticity_scales_inversely_with_the_throat_radius():
    """A wider throat is cheaper — worth stating, since it is the only knob."""
    narrow = abs(tt.null_energy_integral(0.5))
    wide = abs(tt.null_energy_integral(2.0))
    assert narrow == pytest.approx(4.0 * wide, rel=1e-12)


def test_the_price_is_not_claimed_to_be_the_cost_of_the_clock_offset():
    """Separate questions; a frozen metric can carry a large `Δ` for free."""
    note = tt.measure_the_null_energy_price()["what_this_is_NOT"]
    assert "NOT the energy cost of the clock offset" in note
    assert "moving-mouth dynamics" in note


# ── structural: the scattering solver ───────────────────────────────────────

def test_flux_is_conserved_across_the_whole_throat():
    result = tt.measure_the_scattering_is_symmetric_and_unitary()
    assert result["unitarity_holds"]
    assert result["worst_unitarity_residual"] < 1e-10


def test_unitarity_holds_at_frequencies_the_probe_does_not_sample():
    omega = np.array([0.07, 0.33, 0.8, 1.7, 4.0, 9.0])
    reflected, transmitted = tt.scattering_matrix(omega, 1)
    residual = np.abs(reflected) ** 2 + np.abs(transmitted) ** 2 - 1.0
    assert np.max(np.abs(residual)) < 1e-10


def test_the_solver_is_second_order_in_the_step():
    assert tt.measure_the_scattering_is_symmetric_and_unitary()[
        "second_order_in_the_step"]


def test_the_throat_blocks_dc_and_becomes_transparent():
    omega = np.array([0.02, 20.0])
    _, transmitted = tt.scattering_matrix(omega, 0)
    assert abs(transmitted[0]) < 1e-3, "threshold: nearly total reflection"
    assert abs(transmitted[1]) == pytest.approx(1.0, abs=1e-9)


def test_no_fabry_perot_factorisation_is_imposed():
    """The geometry gives one barrier; PR #216's split is a modelling choice."""
    note = tt.measure_the_geometry_is_derived_not_recalled()[
        "it_is_one_barrier_not_two"]
    assert "modelling choice" in note
    assert "does not impose it" in note


# ── the oracle: the threshold law ───────────────────────────────────────────

def test_the_static_conductance_is_the_exact_closed_form():
    """`I₃ = 2/a²`, `g = 2π²/I₃ = π²a²`."""
    for a in (0.5, 1.0, 2.0):
        static = tt.monopole_conductance(a)
        assert static["resistance"] == pytest.approx(2.0 / (a * a), rel=1e-14)
        assert static["conductance"] == pytest.approx(
            math.pi ** 2 * a * a, rel=1e-14)
        assert static["conductance"] == pytest.approx(
            2.0 * math.pi ** 2 / static["resistance"], rel=1e-14)


def test_the_threshold_magnitude_law_holds():
    """`|T₀| → (π/8)(aω)²`, the coefficient the static solution controls."""
    result = tt.measure_the_threshold_law()
    assert result["the_magnitude_law_holds"]
    assert result["magnitude_ratio"][0] == pytest.approx(1.0, abs=5e-3)
    # and it must converge as omega falls, not merely be close at one point
    ratios = result["magnitude_ratio"]
    assert abs(ratios[0] - 1.0) < abs(ratios[-1] - 1.0)


def test_the_threshold_phase_offset_is_a_pure_convention():
    """A constant `+π/2` — i.e. a factor of `i` — across the sampled band."""
    result = tt.measure_the_threshold_law()
    assert result["the_ratio_is_a_constant_phase"]
    assert result["the_convention_factor_is_i"]
    for value in result["phase_of_ratio_over_pi"]:
        assert value == pytest.approx(0.5, abs=5e-3)


def test_the_static_conductance_is_not_asserted_to_be_transmission_at_zero():
    """That would be a false regression, and the module says why."""
    note = tt.measure_the_threshold_law()["the_static_conductance_is_not_T_at_zero"]
    assert "false regression" in note
    assert "COEFFICIENT" in note
    # and the two really are different numbers
    static = tt.monopole_conductance(1.0)["conductance"]
    _, transmitted = tt.scattering_matrix(np.array([0.01]), 0)
    assert abs(abs(transmitted[0]) - static) > 1.0


# ── the verdict ─────────────────────────────────────────────────────────────

def test_phase_closure_and_group_closure_demand_different_offsets():
    """The raw dispersion data behind PR #216's single tuned `Δ`."""
    result = tt.measure_the_closure_offsets_disagree()
    assert result["the_two_offsets_disagree"]
    assert result["worst_disagreement"] > 1.0


def test_the_raw_gap_table_says_it_is_not_the_invariant_statement():
    """It is evaluated at `n = 0`; the branch-free test is in
    `derived_network.closure_residual`, and this table must point there."""
    note = tt.measure_the_closure_offsets_disagree()[
        "THIS_IS_NOT_THE_INVARIANT_STATEMENT"]
    assert "BRANCH DEPENDENT" in note
    assert "derived_network.closure_residual" in note


def test_there_is_no_frequency_independent_transit_time():
    """What survives branch choice: the throat is dispersive, so
    `closure_offset(d_A, d_B, tau_th)` has no constant `tau_th` to take."""
    rows = tt.measure_the_closure_offsets_disagree()["rows"]
    delays = [r["wigner_delay"] for r in rows]
    assert max(delays) > 5.0 * min(delays), "the Wigner delay is not constant"
    assert all(abs(r["second_derivative"]) > 1e-2 for r in rows)


def test_the_closure_disagreement_shrinks_as_the_throat_becomes_transparent():
    """It is a dispersion effect, so it must die where the throat stops
    dispersing — otherwise it would be a numerical artefact."""
    rows = tt.measure_the_closure_offsets_disagree()["rows"]
    gaps = [r["disagreement"] for r in rows]
    assert gaps[0] > gaps[-1]
    assert gaps == sorted(gaps, reverse=True), "monotone in frequency"


def test_the_round_replaces_the_tuned_closure_offset():
    note = tt.measure_the_closure_offsets_disagree()["what_this_replaces"]
    assert "SOLVES for the offset" in note
    assert "frequency dependent" in note


def test_no_perfect_transmission_point_is_found_but_none_is_forbidden():
    """`|Λ| = 1` needs `|T| = 1`. None was found — but positivity of `V` does
    not forbid one, and the module must not claim that it does."""
    result = tt.measure_whether_the_loop_can_close()
    assert result["transmission_is_below_one_at_finite_frequency"]
    verdict = result["so_lambda_cannot_equal_one_exactly"]
    assert "NOT a theorem" in verdict or "NOT A THEOREM" in verdict
    assert "CAN support perfect-transmission resonances" in verdict


def test_the_ultraviolet_slope_is_predicted_not_fitted():
    """First Born on `Ṽ₀(q) = (3π/8a)(3+a|q|)e^{−a|q|}` gives `−4a`."""
    note = tt.measure_whether_the_loop_can_close()[
        "the_slope_is_predicted_not_fitted"]
    assert "Born" in note
    assert "-4a" in note or "-4.0" in note


def test_the_approach_to_unit_transmission_is_exponential():
    """Which makes it a UV problem, not a benign resonance."""
    result = tt.measure_whether_the_loop_can_close()
    assert result["the_approach_to_unity_is_exponential"]
    assert result["log_slope"] < -2.0


def test_the_causal_classification_is_stated():
    classification = tt.measure_whether_the_loop_can_close()[
        "causal_classification"]
    assert "CTC" in classification["D_loop < 0"]
    assert "chronology" in classification["D_loop = 0"]


def test_the_ledger_derives_its_numbers_from_the_measurements():
    ledger = tt.measure_the_traversable_throat_ledger()
    closure = tt.measure_the_closure_offsets_disagree()
    entries = {e["claim"]: e for e in ledger["entries"]}
    row = next(e for k, e in entries.items() if "branch dependent" in
               e["verdict"].lower())
    assert f"{closure['worst_disagreement']:.4f}" in row["evidence"]
    assert "derived_network.closure_residual" in row["evidence"]


def test_the_ledger_records_that_an_advanced_leg_alone_is_not_enough():
    entries = {e["claim"]: e for e in
               tt.measure_the_traversable_throat_ledger()["entries"]}
    row = next(e for k, e in entries.items() if "advanced leg alone" in k)
    assert row["verdict"] == "NO -- PROVED, NOT COMPUTED"


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import traversable_throat_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
