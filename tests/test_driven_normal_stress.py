"""Field source, signed rotor, complete compatibility, and failure paths."""

import copy
import json

import numpy as np
import pytest

from experiments.closure_ledger import driven_normal_stress_probe as probe
from geometrodynamics.waves import driven_normal_stress as d


@pytest.fixture(scope="module")
def report():
    return probe.run_probe()


def test_exact_polynomials_certify_source_and_all_constraint_charges(report):
    c = report["exact_certificate"]
    assert c["all_exact_zero"]
    assert c["norms"] == ["1", "1"] and c["overlap"] == "0"
    assert c["gradient_f"] == c["gradient_g"] == [["3", "0", "0"], ["0", "3", "0"], ["0", "0", "9"]]
    assert c["symmetric_cross_gradient"] == [["0"]*3 for _ in range(3)]
    assert c["six_rotation_charges"] == ["0"]*6
    assert c["four_gradient_charges"] == c["Hamiltonian_dipoles"] == ["0"]*4
    assert c["tensor_equation_residuals"] == ["0"]*9


def test_full_improved_stress_matches_the_action_with_nonzero_momentum(report):
    assert {row["points"] for row in report["spatial_checks"]} == {2048, 6912}
    assert max(row["source_error"] for row in report["spatial_checks"]) < 1e-9
    assert all(row["compatible"] for row in report["spatial_checks"])
    candidate = d.DrivenRotor()
    q, p = candidate.scalar(0.)
    assert np.linalg.norm(q) == pytest.approx(candidate.amplitude)
    assert np.linalg.norm(p) == pytest.approx(4*candidate.amplitude)
    assert abs(q @ p) < 1e-12


def test_compatibility_rejects_both_omitted_charge_controls_and_partial_schema(report):
    controls = report["negative_constraint_controls"]
    left = controls["missed_rotations"]
    assert max(abs(x) for x in left["three_invariant_charges"]) < 1e-9
    assert max(abs(x) for x in left["charges"]["rotations"]) > 1e-3
    assert left["rejected"] and controls["missed_gradients"]["rejected"]
    good = report["spatial_checks"][0]["charges_over_s2"]
    assert d.compatible(good)
    assert not d.compatible({**good, "rotations": good["rotations"][:3]})
    assert not d.compatible({k: v for k, v in good.items() if k != "gradient_dipoles"})


def test_all_equations_hold_on_the_negative_amplitude_branch(report):
    assert report["preparation"]["tensor_A"] < 0
    assert report["exact_certificate"]["A_times_C_over_s2"] == "-3/2"
    assert max(report["balance"]["normalized_residuals"].values()) < 1e-9
    assert min(report["balance"]["necessity_control_residuals"].values()) > 1e-3


def test_unrestricted_field_trajectory_stays_uniaxial_and_director_advances(report):
    e = report["evolution"]
    for key in ("trajectory_error", "refinement_error", "max_cone_distance",
                "max_eigenvalue_distance", "director_projector_error",
                "initial_vs_half_period_axis_overlap", "scalar_free_evolution_error"):
        assert e[key] < 1e-8
    assert np.asarray(e["tensor_coefficients_over_abs_A"]).shape == (401, 5)


def test_free_scalar_and_tensor_do_not_form_a_closed_joint_history():
    candidate = d.DrivenRotor()
    q0, _ = candidate.scalar(0.)
    q1, _ = candidate.scalar(candidate.period)
    # The tensor returns, but even allowing a scalar sign reversal does not
    # close the joint history. The analytic frequency ratio is irrational.
    np.testing.assert_allclose(candidate.orbit(0.)["beta"], candidate.orbit(candidate.period)["beta"], atol=1e-15)
    assert min(np.linalg.norm(q1-q0), np.linalg.norm(q1+q0))/candidate.amplitude > .1


def test_scaling_and_constant_source_exception_are_retained(report):
    assert report["scaling_error"] < 1e-9
    assert max(max(r["normalized_residuals"].values()) for r in report["radius_controls"]) < 1e-9
    c = report["frequency_certificate"]
    assert c["all_exact"] and c["original_speed_excluded"]
    assert c["quadratic_Fourier_residual"] == "0"
    assert c["constant_source_exception_coefficient"] == "0"
    assert c["original_speed_oscillatory_coefficient"] == "184/25"


def test_verdict_is_preparation_and_order_scoped(report):
    assert report["checks_passed"] and set(report["checks"]) == set(d.REQUIRED_CHECKS)
    assert report["verdict"]["normal_balance"] == "LEADING_ORDER_FIELD_SUPPORTED_ROTOR"
    assert report["verdict"]["full_Einstein_matter_evolution"] == "NOT_ESTABLISHED"
    for key in ("preparation_selection", "triangle_map", "Phi_selection"):
        assert report["verdict"][key] == "NOT_DERIVED"
    assert all(d.verdict({"unrelated": True})[key] == "UNRESOLVED" for key in d.VERDICT_FIELDS)


@pytest.mark.parametrize("gate", d.REQUIRED_CHECKS)
@pytest.mark.parametrize("remove", [False, True])
def test_every_failed_or_missing_gate_overwrites_stale_success_and_exits(report, gate, remove, monkeypatch, tmp_path):
    failed = copy.deepcopy(report)
    if remove:
        failed["checks"].pop(gate)
    else:
        failed["checks"][gate] = False
    # Deliberately leave stale checks_passed and verdict fields in memory;
    # the CLI must recalculate both from the mandatory gate list.
    monkeypatch.setattr(probe, "run_probe", lambda progress: failed)
    (tmp_path/"probe.json").write_text('{"checks_passed": true}')
    (tmp_path/"probe.md").write_text("STALE PASS")
    assert probe.main(["--output-dir", str(tmp_path)]) == 1
    saved = json.loads((tmp_path/"probe.json").read_text())
    assert not saved["checks_passed"]
    assert saved["verdict"]["failed_checks"] == ["missing: "+gate if remove else gate]
    assert all(saved["verdict"][key] == "UNRESOLVED" for key in d.VERDICT_FIELDS)
    assert "STALE" not in (tmp_path/"probe.md").read_text()
