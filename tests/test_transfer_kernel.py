"""The retarded outer→inner transfer kernel.

Four kinds of check. The *exact* ones pin the three closed forms `D = 5`
supplies — the rational barrier peak, the integrated potential, and hence the
kernel's high-frequency phase. The *gate* ones are the three conditions this
object had to pass before being read physically: causal support, flux
conservation, and a ringdown matching the external published value. The
*cross-validation* one pins the kernel against a method sharing no code with
it. And the *verdict* ones guard the conclusion — the transfer is not rigid,
and the late-time tail is not claimed.
"""

import math

import numpy as np
import pytest

from geometrodynamics.tangherlini import transfer_kernel as tk


# ── exact: the three closed forms ───────────────────────────────────────────

def test_the_barrier_peak_at_ell_one_is_exactly_one_hundred_over_eighty_one():
    """`V_max = 100/81` at `r² = 9/5`, the mode #270's codes disagreed on."""
    peak = tk.barrier_peak(1)
    assert peak["r_squared"] == pytest.approx(1.8, abs=1e-13)
    assert peak["V_max"] == pytest.approx(100.0 / 81.0, rel=1e-12)


def test_only_ell_one_has_a_rational_barrier_peak():
    """A theorem, not spot checks: `16m⁴+28m²+73` is square only at `m = 2`."""
    result = tk.measure_only_ell_one_has_a_rational_peak()
    assert result["only_m_equals_two"]
    assert result["squares_found_in_range"] == [2]
    assert result["bracketing_holds_for_m_at_least_four"]
    assert result["small_cases"]["2"]["value"] == 441      # = 21^2
    assert not result["small_cases"]["1"]["is_square"]
    assert not result["small_cases"]["3"]["is_square"]


def test_the_bracketing_step_of_the_theorem_is_tight():
    """`(4m²+3)² < 16m⁴+28m²+73 < (4m²+4)²` needs `4m² > 57`, i.e. `m ≥ 4` —
    and must genuinely fail at `m = 3`, or the proof would prove too much."""
    def form(m):
        return 16 * m ** 4 + 28 * m ** 2 + 73
    assert not form(3) < (4 * 9 + 4) ** 2, "the bound must fail at m = 3"
    for m in (4, 5, 10, 100):
        assert (4 * m * m + 3) ** 2 < form(m) < (4 * m * m + 4) ** 2


def test_the_barrier_peak_approaches_the_photon_sphere():
    """`r² → 2`, which PR #274 pinned exactly. A consistency bridge."""
    radii = [tk.barrier_peak(ell)["r_squared"] for ell in (0, 1, 2, 3, 10, 40)]
    assert radii == sorted(radii), "must increase toward the photon sphere"
    assert radii[-1] == pytest.approx(2.0, abs=0.02)
    assert all(r < 2.0 for r in radii), "and approach it from below"


def test_the_integrated_potential_is_ell_times_ell_plus_two_plus_three_halves():
    for ell in (0, 1, 2, 3, 7):
        assert tk.integrated_potential(ell) == pytest.approx(
            ell * (ell + 2) + 1.5, abs=1e-12)


def test_the_high_frequency_phase_constant_is_nine_quarters_at_ell_one():
    assert tk.high_frequency_phase_constant(1) == pytest.approx(2.25, abs=1e-15)


def test_the_quadrature_deficit_is_the_predicted_truncated_tail():
    """The domain truncation is accounted for, not assumed small."""
    result = tk.measure_the_exact_background_anchors()
    for row in result["integrated_potential"]:
        deficit = row["exact"] - row["quadrature_on_the_truncated_domain"]
        assert deficit == pytest.approx(
            row["missing_tail_beyond_the_outer_edge"], rel=0.05)


def test_the_anchors_measurement_agrees_with_the_direct_checks():
    result = tk.measure_the_exact_background_anchors()
    assert result["ell_1_peak_is_exactly_100_over_81"]
    assert result["ell_1_peak_radius_squared_is_exactly_9_over_5"]
    assert result["c_of_ell_1_is_exactly_nine_quarters"]


# ── gate 2: flux conservation ───────────────────────────────────────────────

def test_flux_is_conserved_to_near_machine_precision():
    """`|R|² + |T|² = 1`, imposed nowhere, so it measures the computation."""
    result = tk.measure_the_scattering_is_well_conditioned()
    assert result["unitarity_holds"]
    assert result["worst_unitarity_residual"] < 1e-10


def test_unitarity_holds_at_frequencies_the_probe_does_not_sample():
    omega = np.array([0.05, 0.45, 0.9, 1.2, 2.5, 6.0, 17.0])
    reflected, transmitted = tk.scattering_amplitudes(omega, 1)
    residual = np.abs(reflected) ** 2 + np.abs(transmitted) ** 2 - 1.0
    assert np.max(np.abs(residual)) < 1e-10


def test_the_barrier_blocks_dc_and_passes_high_frequencies():
    result = tk.measure_the_scattering_is_well_conditioned()
    assert result["transmission_goes_to_zero_at_dc"]
    assert result["transmission_goes_to_one_at_high_frequency"]


def test_the_transfer_matrix_is_second_order_in_the_step():
    assert tk.measure_the_scattering_is_well_conditioned()[
        "second_order_in_the_step"]


def test_the_round_says_why_this_is_not_the_shoot_that_failed():
    """Real `ω` is a different problem from `Im ω < 0`, not a repair of it."""
    note = tk.measure_the_scattering_is_well_conditioned()[
        "why_this_is_not_the_failed_shoot"]
    assert "REAL" in note
    assert "well-posed" in note


# ── gate 1: causality ───────────────────────────────────────────────────────

def test_the_kernel_vanishes_before_the_front():
    result = tk.measure_the_kernel_is_causal()
    assert result["the_kernel_is_causal"]
    assert result["worst_acausal_value"] < 1e-3
    assert result["noise_floor_far_from_the_front"] < 1e-5


def test_causality_holds_at_times_the_probe_does_not_sample():
    before = tk.transfer_kernel(np.array([-7.0, -33.0, -71.0, -150.0]), 1)
    assert np.max(np.abs(before)) < 1e-4


def test_the_exact_zero_is_used_as_an_error_bar():
    note = tk.measure_the_kernel_is_causal()["the_exact_zero_is_a_free_error_bar"]
    assert "noise floor" in note
    assert "no reference value" in note


# ── gate 3: the published ringdown ──────────────────────────────────────────

def test_the_kernel_carries_the_published_ringdown():
    result = tk.measure_the_kernel_reproduces_the_published_ringdown()
    assert result["the_ringdown_is_the_published_one"]
    assert result["every_fit_lands_on_the_fundamental_mode"]
    assert result["the_best_fit_is_sub_tenth_percent"]


def test_the_ringdown_is_scored_against_an_external_value():
    """Never against this repository's own solver."""
    result = tk.measure_the_kernel_reproduces_the_published_ringdown()
    assert "arXiv:2107.04815" in result["source"]
    assert result["published"][0] == pytest.approx(1.01601691149, rel=1e-11)


def test_the_ringdown_band_is_reported_rather_than_the_best_row():
    result = tk.measure_the_kernel_reproduces_the_published_ringdown()
    assert len(result["rows"]) >= 6
    assert result["real_part_band"]["max"] > result["real_part_band"]["min"]


# ── the independent cross-check ─────────────────────────────────────────────

def test_the_kernel_convolution_reproduces_the_time_domain_evolution():
    """Two methods sharing no code."""
    result = tk.measure_the_kernel_against_the_time_domain_evolution()
    assert result["the_two_methods_agree"]
    assert result["best"]["peak_relative_max_difference"] < 0.02
    assert result["best"]["peak_relative_rms_difference"] < 0.005


def test_pr_274s_launch_radius_is_recorded_as_inadequate_here_only():
    """It was fine for a quasinormal frequency and wrong for an amplitude."""
    note = tk.measure_the_kernel_against_the_time_domain_evolution()[
        "what_this_exposed"]
    assert "harmless for a quasinormal frequency" in note
    assert "fatal for a transmission ratio" in note


# ── the verdict ─────────────────────────────────────────────────────────────

def test_the_dc_sum_rule_holds():
    """`∫K_reg dt = −1` exactly, from `T(0) = 0`. A constraint, not a fit."""
    result = tk.measure_the_transfer_is_not_rigid()
    assert result["the_sum_rule_holds"]
    assert result["sum_rule_measured"] == pytest.approx(-1.0, abs=5e-3)


def test_the_transfer_is_not_rigid():
    """The deliverable."""
    result = tk.measure_the_transfer_is_not_rigid()
    assert result["the_kernel_is_not_rigid"]
    assert result["memory_absolute_mass"] > 1.5, "memory exceeds the delta"
    assert result["transmission_at_dc"] < 1e-4, "DC is blocked completely"


def test_the_rigidity_claim_states_its_scope():
    """It is a statement about the geometry, not about any BAM kernel."""
    note = tk.measure_the_transfer_is_not_rigid()["scope_of_that_claim"]
    assert "separate question" in note


def test_the_late_time_tail_is_not_claimed():
    """It sits below the causality noise floor, and no exponent is quoted."""
    entries = {e["claim"]: e for e in
               tk.measure_the_transfer_kernel_ledger()["entries"]}
    tail = next(e for k, e in entries.items() if "power-law tail" in k)
    assert tail["verdict"] == "NO"
    assert "No exponent is quoted" in tail["evidence"]


def test_the_ledger_records_what_the_causality_gate_caught():
    caught = tk.measure_what_the_causality_gate_caught()
    for key in ("missing_dc_cell", "gibbs_ringing_from_the_slow_tail"):
        assert "would_have_been_read_as" in caught[key]
    assert "free" in caught["the_general_point"]


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import transfer_kernel_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)


# ── the two review patches ──────────────────────────────────────────────────

def test_the_outer_jost_solutions_have_the_exact_wronskian():
    """`h₊h₋′ − h₋h₊′ = −2i` exactly, at every radius. Checked, not assumed."""
    for x in (0.5, 5.0, 50.0, 500.0):
        plus, minus, d_plus, d_minus = tk.outer_jost_solutions(1, np.array([x]))
        wronskian = (plus * d_minus - minus * d_plus)[0]
        assert wronskian == pytest.approx(-2j, abs=1e-10)


def test_the_jost_solutions_reduce_to_plane_waves_at_large_argument():
    """Normalised so the high-frequency convention is untouched."""
    for x in (200.0, 2000.0):
        plus, minus, _, _ = tk.outer_jost_solutions(1, np.array([x]))
        assert plus[0] / np.exp(1j * x) == pytest.approx(1.0, abs=2e-2)
        assert minus[0] / np.exp(-1j * x) == pytest.approx(1.0, abs=2e-2)


def test_the_subtracted_coefficient_is_the_exact_one_not_a_fit():
    """A fitted `c` would leave `−i(c_exact−c_fit)/ω`, still `1/ω`."""
    result = tk.measure_the_high_frequency_tail_is_the_exact_one()
    assert result["exact"] == pytest.approx(2.25, abs=1e-15)
    assert result["agrees_with_the_exact_value"]
    assert "still falls only as 1/w" in result[
        "why_the_exact_value_is_the_one_subtracted"]


def test_the_fitted_coefficient_is_independent_of_the_outer_edge():
    """Edge-independence is what distinguishes a bounded shortfall from a
    truncation drift — the plane-wave version drifted."""
    result = tk.measure_the_high_frequency_tail_is_the_exact_one()
    assert result["independent_of_the_outer_edge"]
    assert result["spread_across_outer_edges"] < 1e-4


def test_the_lowest_frequency_bin_is_inside_the_centrifugal_tail():
    """Which is why plane-wave matching is the wrong basis there."""
    result = tk.measure_the_low_frequency_outer_matching_is_converged()
    assert result["the_lowest_bin_is_inside_the_centrifugal_tail"]
    assert result["outer_turning_scale"][0] > 2.0 * tk.OUTER_EDGE


def test_the_low_frequency_spectrum_is_converged_in_the_outer_edge():
    result = tk.measure_the_low_frequency_outer_matching_is_converged()
    assert result["converged_in_the_outer_edge"]
    assert max(result["relative_spread_across_outer_edges"]) < 1e-3


def test_the_cross_check_residual_halves_with_the_launch_radius():
    """Tracking the exactly-known `L/r*_launch` is what identifies the residual
    as launch placement rather than a disagreement between methods."""
    result = tk.measure_the_kernel_against_the_time_domain_evolution()
    assert result["the_residual_halves_as_the_launch_radius_doubles"]
    assert all(1.7 < r < 2.3 for r in result["successive_ratios"])


def test_the_earlier_cross_check_number_is_recorded_as_flattered():
    """`0.92%` was two errors partly cancelling; `2.73%` is the honest one."""
    note = tk.measure_the_kernel_against_the_time_domain_evolution()[
        "an_earlier_number_was_flattered_by_cancellation"]
    assert "0.92%" in note and "2.73%" in note
    assert "larger number is the honest one" in note
