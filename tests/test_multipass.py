"""
Tests for the multi-pass throat-knot lepton-mass hypothesis.

Verifies calibration invariants (anchor exact, RMS-log residual consistent
with reported rel errors), amplifier laws (zero-param vs one-param
behaviour), scan ordering, and the central physics claim of this module:
Möbius odd-depth (1, 3, 5) with a single-parameter power-law amplifier
reduces the bare ~99% lepton-mass residual to ≲ 15% on both μ and τ.
"""

from __future__ import annotations

import numpy as np
import pytest

from geometrodynamics.tangherlini import (
    AMPLIFIERS,
    DEFAULT_ASSIGNMENT_RADIAL,
    KNOT_DEPTHS,
    MultiPassFit,
    PDG_LEPTON_MASSES_MEV,
    amp_exponential,
    amp_linear,
    amp_power,
    evaluate_multipass,
    format_multipass_fit,
    format_multipass_scan,
    scan_multipass,
)


# Default N=100 is fast and accurate enough for all of these checks.
_N = 100


@pytest.fixture
def best_fit() -> MultiPassFit:
    return evaluate_multipass(
        DEFAULT_ASSIGNMENT_RADIAL, depth="odd", amplifier="power", N=_N,
    )


class TestAmplifiers:
    def test_linear_is_identity_in_k(self):
        for k in (1, 2, 3, 5, 7):
            assert amp_linear(k) == float(k)

    def test_power_reduces_to_unity_at_p0(self):
        for k in (1, 2, 5):
            assert amp_power(k, 0.0) == pytest.approx(1.0)

    def test_power_matches_k_at_p1(self):
        for k in (1, 2, 5):
            assert amp_power(k, 1.0) == pytest.approx(float(k))

    def test_exponential_at_beta0_is_one(self):
        for k in (1, 2, 5):
            assert amp_exponential(k, 0.0) == pytest.approx(1.0)

    def test_exponential_is_positive_and_monotone(self):
        prev = 0.0
        for k in (1, 2, 3, 4):
            v = amp_exponential(k, 0.5)
            assert v > prev
            prev = v

    def test_amplifier_registry_has_three_entries(self):
        assert set(AMPLIFIERS) == {"linear", "power", "exponential"}


class TestKnotDepths:
    def test_canonical_sequences_have_length_three(self):
        for name, seq in KNOT_DEPTHS.items():
            assert len(seq) == 3, f"{name} length {len(seq)}"

    def test_first_element_is_one_for_all(self):
        # Electron is always k=1 by convention (single-pass defect).
        for name, seq in KNOT_DEPTHS.items():
            assert seq[0] == 1, f"{name} starts at {seq[0]}"

    def test_sequences_are_non_decreasing(self):
        for name, seq in KNOT_DEPTHS.items():
            for a, b in zip(seq, seq[1:]):
                assert a <= b, f"{name} non-monotone"


class TestCalibrationInvariants:
    def test_anchor_exactly_matches_pdg(self, best_fit):
        e_idx = next(
            i for i, a in enumerate(best_fit.assignment) if a.name == "electron"
        )
        assert best_fit.masses_pred_mev[e_idx] == pytest.approx(
            PDG_LEPTON_MASSES_MEV["electron"], rel=1e-12
        )
        assert abs(best_fit.rel_errors[e_idx]) < 1e-12

    def test_rms_log_residual_matches_reported_errors(self, best_fit):
        # RMS log over non-anchor modes must equal the stored field.
        anchor_idx = next(
            i for i, a in enumerate(best_fit.assignment)
            if a.name == best_fit.calibration_lepton
        )
        others = [i for i in range(len(best_fit.assignment)) if i != anchor_idx]
        log_res = np.array(
            [
                np.log(best_fit.masses_pred_mev[i])
                - np.log(best_fit.masses_pdg_mev[i])
                for i in others
            ]
        )
        expected = float(np.sqrt(np.mean(log_res ** 2)))
        assert best_fit.rms_log_residual == pytest.approx(expected, rel=1e-9)

    def test_linear_amp_has_zero_fit_param(self):
        fit = evaluate_multipass(
            DEFAULT_ASSIGNMENT_RADIAL, depth="odd", amplifier="linear", N=_N,
        )
        assert fit.fit_param == 0.0


class TestPhysicsClaim:
    """The headline result: Möbius odd-depth × power amp brings the bare
    ~99% residual down to ≲ 15% on both μ and τ predictions."""

    def test_best_fit_brings_residuals_under_15pct(self, best_fit):
        # electron is fixed; check μ and τ
        by_name = {a.name: i for i, a in enumerate(best_fit.assignment)}
        assert abs(best_fit.rel_errors[by_name["muon"]]) < 0.15
        assert abs(best_fit.rel_errors[by_name["tau"]]) < 0.15

    def test_power_exponent_in_physical_range(self, best_fit):
        # Depth=odd, amp=power should land near p ~ 4.4 with these ω's.
        assert 3.5 < best_fit.fit_param < 5.5

    def test_odd_power_beats_single_power_substantially(self):
        odd = evaluate_multipass(
            DEFAULT_ASSIGNMENT_RADIAL, depth="odd", amplifier="power", N=_N,
        )
        single = evaluate_multipass(
            DEFAULT_ASSIGNMENT_RADIAL, depth="single", amplifier="power", N=_N,
        )
        # depth=single kills the amplifier → same as the bare model
        assert odd.rms_log_residual < 0.5 * single.rms_log_residual

    def test_exponential_cannot_fit_lepton_spacing_well(self):
        # PDG has τ/μ ≈ 17 while μ/e ≈ 207, so a pure exponential in k
        # cannot reproduce both ratios — rms_log should stay >> 0.
        fit = evaluate_multipass(
            DEFAULT_ASSIGNMENT_RADIAL,
            depth="linear", amplifier="exponential", N=_N,
        )
        assert fit.rms_log_residual > 0.3


class TestScan:
    def test_scan_returns_ordered_by_rms(self):
        results = scan_multipass(DEFAULT_ASSIGNMENT_RADIAL, N=_N)
        for a, b in zip(results, results[1:]):
            assert a.rms_log_residual <= b.rms_log_residual

    def test_scan_best_entry_matches_evaluate_multipass(self):
        results = scan_multipass(DEFAULT_ASSIGNMENT_RADIAL, N=_N)
        top = results[0]
        # Should be odd / power combo — the global best for radial.
        assert top.depth_name == "odd"
        assert top.amplifier_name == "power"

    def test_scan_enumerates_all_combinations(self):
        results = scan_multipass(DEFAULT_ASSIGNMENT_RADIAL, N=_N)
        assert len(results) == len(KNOT_DEPTHS) * len(AMPLIFIERS)


class TestAPIValidation:
    def test_mismatched_depth_length_raises(self):
        with pytest.raises(ValueError):
            evaluate_multipass(
                DEFAULT_ASSIGNMENT_RADIAL, depth=(1, 2), N=_N,
            )

    def test_unknown_amplifier_raises(self):
        with pytest.raises(ValueError):
            evaluate_multipass(
                DEFAULT_ASSIGNMENT_RADIAL, amplifier="nope", N=_N,
            )

    def test_unknown_calibrate_to_raises(self):
        with pytest.raises(ValueError):
            evaluate_multipass(
                DEFAULT_ASSIGNMENT_RADIAL, calibrate_to="quark", N=_N,
            )


class TestFormatting:
    def test_format_fit_contains_all_lepton_names(self, best_fit):
        text = format_multipass_fit(best_fit)
        for name in ("electron", "muon", "tau"):
            assert name in text
        assert "power" in text

    def test_format_scan_has_header_and_rows(self):
        results = scan_multipass(DEFAULT_ASSIGNMENT_RADIAL, N=_N)
        text = format_multipass_scan(results)
        assert "depth" in text
        assert "rms_log" in text
        # One row per (depth, amp) combination.
        assert text.count("\n") >= len(results)
