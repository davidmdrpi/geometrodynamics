"""Tests for multi-pass knotted lepton spectrum surrogate."""

import numpy as np

from geometrodynamics.tangherlini.lepton_spectrum import (
    Crossing,
    calibrate_electron_predict_heavier,
    compute_knotted_lepton_spectrum,
    compute_tunneling_envelope,
    tune_transport_and_resistance,
)


def test_knotted_spectrum_returns_requested_depths():
    spec = compute_knotted_lepton_spectrum(depths=(1, 3), n_points=24)
    assert set(spec) == {1, 3}
    assert spec[1] > 0.0
    assert spec[3] > 0.0


def test_resistance_amplifies_depth_ladder():
    flat = compute_knotted_lepton_spectrum(
        depths=(1, 3, 5),
        n_points=24,
        resistance_model="none",
        resistance_scale=0.0,
    )
    amped = compute_knotted_lepton_spectrum(
        depths=(1, 3, 5),
        n_points=24,
        resistance_model="curvature",
        resistance_scale=0.9,
    )
    # Stronger re-entry resistance should raise high-depth separation.
    assert (amped[5] - amped[1]) > (flat[5] - flat[1])


def test_exponential_resistance_amplifies_more_than_curvature():
    curvature = compute_knotted_lepton_spectrum(
        depths=(1, 5),
        n_points=24,
        resistance_model="curvature",
        resistance_scale=0.5,
    )
    exponential = compute_knotted_lepton_spectrum(
        depths=(1, 5),
        n_points=24,
        resistance_model="exponential",
        resistance_scale=0.5,
    )
    assert (exponential[5] - exponential[1]) > (curvature[5] - curvature[1])


def test_projected_crossings_do_not_change_operator():
    base = compute_knotted_lepton_spectrum(depths=(3,), n_points=24)
    projected = compute_knotted_lepton_spectrum(
        depths=(3,),
        n_points=24,
        crossings=[Crossing(pass_a=1, pass_b=2, identified=False, coupling=10.0)],
    )
    assert abs(base[3] - projected[3]) < 1e-9


def test_hard_pinhole_delta_targets_k3_k5_only():
    base = compute_knotted_lepton_spectrum(depths=(1, 3, 5), n_points=24)
    pinhole = compute_knotted_lepton_spectrum(
        depths=(1, 3, 5),
        n_points=24,
        hard_pinhole_gamma=8.0,
        hard_pinhole_depths=(3, 5),
    )
    d1 = abs(base[1] - pinhole[1])
    d3 = abs(base[3] - pinhole[3])
    d5 = abs(base[5] - pinhole[5])
    assert d3 > 1e-6
    assert d5 > 1e-6
    assert d3 >= d1
    assert d5 >= d1


def test_winding_cost_hits_deeper_depth_more_strongly():
    base = compute_knotted_lepton_spectrum(depths=(3, 5), n_points=24, hard_pinhole_gamma=0.0)
    pin = compute_knotted_lepton_spectrum(depths=(3, 5), n_points=24, hard_pinhole_gamma=10.0)
    d3 = abs(pin[3] - base[3])
    d5 = abs(pin[5] - base[5])
    assert d3 > 1e-6
    assert d5 > 1e-6
    assert d5 > 0.5 * d3


def test_identified_crossings_modify_spectrum():
    base = compute_knotted_lepton_spectrum(depths=(3,), n_points=24)
    identified = compute_knotted_lepton_spectrum(
        depths=(3,),
        n_points=24,
        crossings=[Crossing(pass_a=1, pass_b=2, identified=True, coupling=10.0)],
    )
    assert abs(base[3] - identified[3]) > 1e-5


def test_electron_calibration_is_exact_by_construction():
    fit = calibrate_electron_predict_heavier(
        depths=(1, 3, 5),
        n_points=24,
    )
    assert abs(fit.predicted_mev[1] - 0.51099895) < 1e-12


def test_tuning_returns_grid_member_parameters():
    phase_grid = (0.0, 0.5)
    resistance_grid = (0.2, 0.6)
    fit = tune_transport_and_resistance(
        phase_grid=phase_grid,
        resistance_grid=resistance_grid,
        n_points=24,
    )
    assert fit.phase_per_pass in phase_grid
    assert fit.resistance_scale in resistance_grid


def test_depth_cost_ablation_modes_change_spectrum():
    both = compute_knotted_lepton_spectrum(depths=(1, 3, 5), n_points=24, depth_cost_mode="both")
    diag = compute_knotted_lepton_spectrum(depths=(1, 3, 5), n_points=24, depth_cost_mode="diag_only")
    tun = compute_knotted_lepton_spectrum(depths=(1, 3, 5), n_points=24, depth_cost_mode="tunnel_only")
    assert any(abs(both[k] - diag[k]) > 1e-10 for k in (1, 3, 5))
    assert all(np.isfinite(tun[k]) for k in (1, 3, 5))


def test_winding_mode_and_envelope_available():
    env_delta = compute_tunneling_envelope(depth=3, phase_per_pass=0.1, winding_mode="delta")
    env_max = compute_tunneling_envelope(depth=3, phase_per_pass=0.1, winding_mode="max")
    assert env_max.shape == (3, 3)
    assert env_max[0, 1] > 0.0
    assert env_delta[0, 2] > env_max[0, 2]


def test_k_uplift_targets_k5_sector():
    base = compute_knotted_lepton_spectrum(depths=(1, 3, 5), n_points=24, k_uplift_beta=0.0)
    uplift = compute_knotted_lepton_spectrum(depths=(1, 3, 5), n_points=24, k_uplift_beta=50.0)
    d1 = abs(base[1] - uplift[1])
    d3 = abs(base[3] - uplift[3])
    d5 = abs(base[5] - uplift[5])
    assert uplift[5] > base[5]
    assert d5 > d3
    assert d5 > d1
