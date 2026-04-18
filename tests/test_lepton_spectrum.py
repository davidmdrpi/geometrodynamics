"""
Tests for the charged-lepton mass calibration harness.

These verify solver wiring, calibration invariants (anchor hits PDG
exactly, ω-ratios reproduced by mass ratios, throat length positive),
and the sign/ordering of predictions.  They do **not** assert that the
bare Tangherlini spectrum reproduces the measured lepton ratios — it
does not, and quantifying that gap is the point of the module.
"""

from __future__ import annotations

import numpy as np
import pytest

from geometrodynamics.tangherlini import (
    DEFAULT_ASSIGNMENT_ANGULAR,
    DEFAULT_ASSIGNMENT_RADIAL,
    HBAR_C_MEV_FM,
    PDG_LEPTON_MASSES_MEV,
    LeptonAssignment,
    LeptonSpectrumReport,
    compute_lepton_spectrum,
    format_report,
)


# Reuse a single solver run per assignment (fast) but recompute per test
# to keep tests independent — N=100 is accurate enough for the checks.
@pytest.fixture
def radial_report() -> LeptonSpectrumReport:
    return compute_lepton_spectrum(DEFAULT_ASSIGNMENT_RADIAL, N=100)


@pytest.fixture
def angular_report() -> LeptonSpectrumReport:
    return compute_lepton_spectrum(DEFAULT_ASSIGNMENT_ANGULAR, N=100)


class TestCalibrationInvariants:
    def test_anchor_matches_pdg_exactly(self, radial_report):
        e = next(m for m in radial_report.modes if m.name == "electron")
        assert abs(e.rel_error) < 1e-12
        assert e.mass_pred_mev == pytest.approx(PDG_LEPTON_MASSES_MEV["electron"])

    def test_throat_length_positive_and_finite(self, radial_report):
        R = radial_report.calibration_length_fm
        assert np.isfinite(R) and R > 0
        # Sanity: R should be O(Compton wavelength of the electron ~ 386 fm)
        assert 100.0 < R < 1000.0

    def test_mass_scale_consistent_with_hbar_c(self, radial_report):
        # For the anchor: m_pred = ω · (ℏc / R_throat)
        e = next(m for m in radial_report.modes if m.name == "electron")
        scale = HBAR_C_MEV_FM / radial_report.calibration_length_fm
        assert e.omega * scale == pytest.approx(e.mass_pred_mev, rel=1e-12)

    def test_calibrate_to_muon(self):
        rep = compute_lepton_spectrum(
            DEFAULT_ASSIGNMENT_RADIAL, calibrate_to="muon", N=100
        )
        assert rep.calibration_lepton == "muon"
        mu = next(m for m in rep.modes if m.name == "muon")
        assert abs(mu.rel_error) < 1e-12


class TestSolverHookup:
    def test_radial_assignment_yields_three_modes(self, radial_report):
        assert len(radial_report.modes) == 3
        names = [m.name for m in radial_report.modes]
        assert names == ["electron", "muon", "tau"]

    def test_angular_assignment_yields_three_modes(self, angular_report):
        assert len(angular_report.modes) == 3

    def test_eigenfrequencies_monotone_increasing(self, radial_report):
        oms = [m.omega for m in radial_report.modes]
        for a, b in zip(oms, oms[1:]):
            assert a < b

    def test_predicted_mass_ordering_matches_frequency_ordering(self, radial_report):
        # Predictions should preserve the ω ordering.
        ms = [m.mass_pred_mev for m in radial_report.modes]
        for a, b in zip(ms, ms[1:]):
            assert a < b

    def test_ratios_match_frequency_ratios(self, radial_report):
        # m_pred_i / m_pred_anchor = ω_i / ω_anchor exactly
        anchor = next(
            m for m in radial_report.modes
            if m.name == radial_report.calibration_lepton
        )
        for m in radial_report.modes:
            if m.name == anchor.name:
                continue
            key = f"{m.name}/{anchor.name}"
            assert radial_report.ratios_pred[key] == pytest.approx(
                m.omega / anchor.omega, rel=1e-12
            )


class TestPhysicsGap:
    """Document (rather than assert) the current model's mass-ratio gap.

    The bare Tangherlini ladder is near-linear in (l, n); the lepton
    ratios are near-geometric with huge steps.  We *assert* the gap
    exists and is in the expected regime — if a future change makes
    the bare model agree to better than 50%, one of these will fail
    and prompt a README update.
    """

    def test_radial_tau_muon_both_far_below_pdg(self, radial_report):
        # Anchor is the electron — muon and tau predictions should be
        # heavily under PDG with the bare radial ladder.
        for name in ("muon", "tau"):
            m = next(mm for mm in radial_report.modes if mm.name == name)
            assert m.rel_error < -0.5, (
                f"{name} prediction {m.mass_pred_mev:.3g} MeV unexpectedly close "
                f"to PDG {m.mass_pdg_mev:.3g} MeV"
            )

    def test_angular_gap_is_even_larger_than_radial(
        self, radial_report, angular_report
    ):
        # Angular ladder is flatter than the radial ladder, so its worst
        # relative error should be no smaller.
        assert angular_report.worst_mass_rel_err >= radial_report.worst_mass_rel_err


class TestAPIValidation:
    def test_unknown_name_raises(self):
        with pytest.raises(ValueError):
            compute_lepton_spectrum(
                (LeptonAssignment("quark", l=1, n=0),),
            )

    def test_calibrate_to_missing_name_raises(self):
        with pytest.raises(ValueError):
            compute_lepton_spectrum(
                DEFAULT_ASSIGNMENT_RADIAL, calibrate_to="tau_prime"
            )

    def test_format_report_mentions_all_leptons(self, radial_report):
        text = format_report(radial_report)
        for name in ("electron", "muon", "tau"):
            assert name in text
        assert "R_throat" in text
