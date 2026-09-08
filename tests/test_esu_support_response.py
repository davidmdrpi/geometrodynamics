"""Proper clock, propagated constraints, support dependence, and failure gates."""

import copy
import json
import math

import numpy as np
import pytest

from experiments.closure_ledger import esu_support_response_probe as probe
from geometrodynamics.waves import esu_support_response as esu


@pytest.fixture(scope="module")
def report():
    return probe.run_probe()


def test_metric_variation_retains_lapse_and_curvature(report):
    derivation = report["metric_derivation"]
    assert derivation["all_exact_zero"]
    assert derivation["Hamiltonian_residual"] == "0"
    assert derivation["momentum_residuals"] == ["0"]*3
    assert derivation["spatial_Einstein_residuals"] == ["0"]*9
    assert derivation["wave_operator_residual"] == "0"
    assert max(r["Richardson_scaled_error"] for r in report["metric_variation_checks"]) < 1e-7


def test_full_improved_stress_fixes_scalar_anisotropic_projection(report):
    assert {r["points"] for r in report["stress_checks"]} == {2048, 6912}
    assert max(r[k] for r in report["stress_checks"] for k in (
        "rho", "pressure", "momentum", "anisotropic_scalar_projection")) < 1e-9
    assert report["spatial_refinement_error"] < 1e-9


def test_proper_clock_is_not_replaced_by_coordinate_acceleration(report):
    r = report["exact_initial_certificate"]
    assert r["proper_total"] == "-7976/875"
    assert r["coordinate_total"] == "55096/875"
    assert r["proper_by_degree"][0] == "-40"  # The mean is retained.
    for row in report["initial_force_checks"]+report["radius_checks"]:
        assert row["proper_quadrature_coefficient"] == pytest.approx(-7976/875, abs=1e-9)
        assert row["coordinate_quadrature_coefficient"] == pytest.approx(55096/875, abs=1e-9)
        assert row["proper_conversion_coefficient"] == pytest.approx(-7976/875, abs=1e-9)


def test_constraints_propagate_in_the_independent_fluid_evolution(report):
    for row in report["support_controls"]:
        assert row["independent_fluid_scaled_error"] < 1e-8
        assert row["unused_Hamiltonian_residual"] < 1e-8
        assert row["unused_momentum_residual"] < 1e-8
        assert row["curvature_trace_residual"] < 1e-9
        assert row["ODE_refinement_scaled_error"] < 1e-8
        # Rigidity is initial data only; a later responding density was evolved.
        density = np.array(row["support_density_group_amplitudes"])
        np.testing.assert_array_equal(density[0], np.zeros(4))
        assert np.linalg.norm(density[-1]) > 1e-6


def test_same_metric_geodesic_hessian_recovers_both_clock_signs(report):
    for row in report["initial_force_checks"]+report["radius_checks"]:
        clocks = row["same_metric_two_clocks"]
        assert clocks["coordinate"]["Richardson_coefficient"] == pytest.approx(55096/875, abs=1e-7)
        assert clocks["proper"]["Richardson_coefficient"] == pytest.approx(-7976/875, abs=1e-7)
        for clock in clocks.values():
            assert clock["Richardson_scaled_error"] < 1e-7


def test_geodesic_clock_gate_detects_a_coordinate_clock_substitution(monkeypatch):
    model = esu.SupportResponse()
    spatial = esu.SpatialFields(model)
    actual = spatial.exact_initial_clock_accelerations

    def wrong_clock(epsilon):
        clocks = actual(epsilon)
        return {**clocks, "proper": clocks["coordinate"]}

    monkeypatch.setattr(spatial, "exact_initial_clock_accelerations", wrong_clock)
    row = probe.initial_force_checks(model, spatial)
    # The old first-order conversion still passes; the independent check
    # rejects the wrong derivative by the full separation of the coefficients.
    assert row["clock_conversion_scaled_error"] < 1e-9
    assert row["same_metric_two_clocks"]["proper"]["coefficient_error"] > 70.
    assert not probe.clock_conversion_passes([row])
    checks = dict.fromkeys(esu.REQUIRED_CHECKS, True)
    checks["clock_conversion"] = probe.clock_conversion_passes([row])
    verdict = esu.verdict(checks)
    assert verdict["cancellation"] == "UNRESOLVED"
    assert verdict["failed_checks"] == ["clock_conversion"]


@pytest.mark.parametrize("radius,cs2", [(.7, 0.), (1., .2), (1., 1/3), (2., 1.)])
def test_fluid_only_frequency_from_independent_linearized_equations(radius, cs2):
    model = esu.SupportResponse(radius=radius, sound_speed_squared=cs2)
    initial = model.initial(independent_fluid=True)
    for index, degree in enumerate(esu.DEGREES):
        # Difference the affine, driven RHS to cancel the scalar source.
        # This perturbation obeys the homogeneous Hamiltonian constraint.
        perturbation = np.zeros(16)
        perturbation[index] = 1.
        perturbation[8+index] = 2*model.L[index]/model.kappa
        actual = (model.fluid_rhs(0., initial+perturbation)-model.fluid_rhs(0., initial))[4+index]
        expected = (1-cs2*(degree*(degree+2)-3))/radius**2
        assert actual == pytest.approx(expected, abs=1e-12)
    assert (1+3*cs2)/radius**2 > 0  # Homogeneous scale instability is present.


def test_initial_preparation_is_EOS_independent_but_later_response_is_not(report):
    controls = report["support_controls"]
    assert len({r["sound_speed_squared"] for r in controls}) == 4
    for r in controls:
        assert r["initial_proper_coefficient"] == pytest.approx(-7976/875, abs=1e-10)
        assert r["TT_force_over_cubic_scale"][0] == 0.
    first = np.array(controls[0]["coordinate_scalar_force_over_cubic_scale"])
    last = np.array(controls[-1]["coordinate_scalar_force_over_cubic_scale"])
    assert np.linalg.norm(first-last) > 1.  # The later parameter dependence was not frozen away.


def test_cubic_scaling_and_induced_TT_control(report):
    assert report["force_scaling_error"] < 1e-8
    assert report["metric_scaling_error"] < 1e-8
    assert report["induced_TT_scaled_error"] < 1e-8
    for r in report["radius_checks"]:
        assert r["proper_coefficient_error"] < 1e-9


def test_continuous_metric_bounds_keep_the_mean_and_cover_samples(report):
    for r in report["support_controls"]:
        b = r["continuous_all_space_potential_bounds"]
        assert r["sampled_max_abs_psi"] <= b["psi"] < .05
        assert r["sampled_max_abs_alpha"] <= b["alpha"] < .05


def test_verdict_requires_every_gate_and_remains_scoped(report):
    assert report["checks_passed"], report["checks"]
    assert set(report["checks"]) == set(esu.REQUIRED_CHECKS)
    assert report["verdict"]["cancellation"] == "IDENTICAL_CANCELLATION_EXCLUDED_FOR_FROZEN_PREPARATION"
    assert report["verdict"]["BAM_support_selection"] == "NOT_DERIVED"
    assert report["verdict"]["complete_evolution"] == "NOT_ESTABLISHED"
    missing = dict(report["checks"])
    missing.pop("clock_conversion")
    result = esu.verdict(missing)
    assert all(result[k] == "UNRESOLVED" for k in esu.VERDICT_FIELDS)
    assert "missing: clock_conversion" in result["failed_checks"]


def test_failed_cli_overwrites_a_stale_pass(report, monkeypatch, tmp_path):
    failed = copy.deepcopy(report)
    failed["checks"].pop("clock_conversion")
    failed["checks_passed"] = not esu.failed_checks(failed["checks"])
    failed["verdict"] = esu.verdict(failed["checks"])
    monkeypatch.setattr(probe, "run_probe", lambda progress: failed)
    (tmp_path/"probe.md").write_text("STALE PASS")
    (tmp_path/"probe.json").write_text('{"checks_passed":true}')
    assert probe.main(["--output-dir", str(tmp_path)]) == 1
    saved = json.loads((tmp_path/"probe.json").read_text())
    assert saved["checks_passed"] is False
    assert all(saved["verdict"][key] == "UNRESOLVED" for key in esu.VERDICT_FIELDS)
    assert "missing: clock_conversion" in (tmp_path/"probe.md").read_text()
    assert "STALE" not in (tmp_path/"probe.md").read_text()
