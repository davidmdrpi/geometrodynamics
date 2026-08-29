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
    """`ℓ = 0, 2, 3` carry `√13`, `√1621`, `√57`. Checked by irrationality of
    the discriminant, not by eyeballing decimals."""
    for ell in (0, 2, 3):
        limit = (ell + 1) ** 2 - 0.25
        discriminant = 4.0 * (2.25 - limit) ** 2 + 27.0 * limit
        root = math.sqrt(discriminant)
        assert abs(root - round(root)) > 1e-6, f"l={ell} unexpectedly rational"
    limit = 3.75
    discriminant = 4.0 * (2.25 - limit) ** 2 + 27.0 * limit
    assert math.sqrt(discriminant) == pytest.approx(10.5, abs=1e-12)


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


def test_the_agreement_improves_as_the_pulse_is_launched_further_out():
    """The residual tracks `V` at the launch point — which is what identifies
    it as a launch-placement effect rather than a solver disagreement."""
    result = tk.measure_the_kernel_against_the_time_domain_evolution()
    assert result["the_error_tracks_the_potential_at_launch"]
    rows = result["rows"]
    assert rows[0]["potential_at_launch"] > rows[-1]["potential_at_launch"]
    assert rows[0]["peak_relative_max_difference"] > 10 * \
        rows[-1]["peak_relative_max_difference"]


def test_pr_274s_launch_radius_is_recorded_as_inadequate_here_only():
    """It was fine for a quasinormal frequency and wrong for an amplitude."""
    note = tk.measure_the_kernel_against_the_time_domain_evolution()[
        "what_this_exposed"]
    assert "ringdown did not" in note


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
