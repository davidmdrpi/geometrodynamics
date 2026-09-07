"""Independent field/action checks for the d8dc90d reciprocal TT freeze."""

import copy
import json
import math

import numpy as np
import pytest

from geometrodynamics.bulk.tt_triangle_rotor import TensorModel
from geometrodynamics.waves import reciprocal_scalar_tt as rt
from experiments.closure_ledger import reciprocal_scalar_tt_probe as probe


@pytest.fixture(scope="module")
def report():
    return probe.run_probe()


@pytest.mark.parametrize("degree", [1, 3, 5])
def test_complete_multiplet_and_casimir(degree):
    h = rt.harmonic_multiplet(degree)
    info = h.algebra_checks()
    assert info.pop("dimension") == info.pop("expected_dimension") == (degree+1)**2
    assert max(info.values()) < 1e-10
    # Check the Lie bracket as well as the Casimir; a wrong sign convention
    # can preserve the latter and spoil frame covariance.
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        assert np.linalg.norm(h.D[i] @ h.D[j]-h.D[j] @ h.D[i]-2*h.D[k]) < 1e-10


def test_polynomial_frame_matches_the_inherited_quaternion_geometry():
    from geometrodynamics.waves.backreaction import left_invariant_frame
    points, _ = rt.sphere_quadrature(3, 8)
    for x in points[::11]:
        got = np.einsum("iab,b->ia", rt.QUATERNION_DERIVATIVES, x)
        np.testing.assert_allclose(got, left_invariant_frame(x), atol=1e-14)


@pytest.mark.parametrize("radius,kappa", [(.7, .4), (1., 1.), (2., 3.)])
def test_source_and_action_normalization(radius, kappa):
    model = rt.ReciprocalModel(radius=radius, kappa=kappa)
    inherited = TensorModel(radius, kappa)
    assert model.C == inherited.normalization
    assert model.omega_tensor2 == inherited.omega2
    assert model.omega_scalar2 == 16/radius**2
    axis = np.array([1., -2., 3.])/math.sqrt(14)
    q = model.multiplet.coherent(axis)
    expected = 6*(np.outer(axis, axis)-np.eye(3)/3)/radius**2
    np.testing.assert_allclose(rt.tensor(model.source(q)), expected, atol=1e-12)
    # Kinetic normalization is checked by physical-space integration.
    points, weights = rt.sphere_quadrature(8, 16, radius)
    jets, _ = rt.scalar_jets(model, q, q, -model.omega_scalar2*q, points)
    assert abs(weights @ jets["phi"][:, 0]**2-1) < 1e-12
    assert abs(weights @ jets["dt"][:, 0]**2-1) < 1e-12


def test_degree_one_null_coupling_is_an_operator_identity():
    assert np.linalg.norm(rt.ReciprocalModel(1).F) == 0.


def test_improved_stress_and_pointwise_scalar_equation(report):
    assert max(r["source_scaled_error"] for r in report["stress_checks"]) < 1e-9
    assert report["pointwise_wave_error"] < 1e-9
    # The nonminimal term is present in the inherited routine; agreement is
    # not established by comparing two copies of q^T F q.
    assert {r["points"] for r in report["stress_checks"]} == {2048, 6912}


def test_static_metric_variation_includes_conformal_curvature(report):
    assert all(3.5 < r < 4.5 for r in report["action_expansion"]["halving_ratios"])
    model = rt.ReciprocalModel()
    b, P, q, p = model.unpack(rt.primary_data(model))
    expected = (p @ p-model.omega_scalar2*(q @ q))/2
    assert abs(model.static_scalar_lagrangian(np.zeros((3, 3)), q, p)-expected) < 1e-13


def test_reciprocal_forces_follow_one_hamiltonian(report):
    r = report["hamiltonian_derivatives"]
    assert r["hamilton_equations_scaled_error"] < 1e-7
    assert r["mixed_derivative_scaled_error"] < 1e-10
    assert r["one_way_mixed_derivative_defect"] > 1e-3


def test_energy_refinement_and_free_controls(report):
    h = report["history"]
    assert max(h["fine_relative_energy_drift"], h["coarse_relative_energy_drift"]) < 1e-8
    assert h["max_absolute_state_difference_between_tolerances"] < 1e-8
    assert h["max_beta_frobenius"] < .05
    assert max(report["free_controls"].values()) < 1e-9


def test_covariance_and_director_even_observable(report):
    assert max(report["covariance"].values()) < 1e-9
    assert report["initial_Q_identity_error"] < 1e-10
    n = np.array([1., 2., 3.])/math.sqrt(14)
    m = np.array([2., -1., 0.])/math.sqrt(5)
    for sign in (-1, 1):
        beta = .01*(np.outer(sign*n, sign*n)-np.eye(3)/3)
        assert abs(m @ beta @ m-.01*((m @ n)**2-1/3)) < 1e-14


def test_nonzero_omitted_constraint_has_a_field_certificate(report):
    r = report["constraint_certificate"]
    assert r["certificate_error"] < 1e-10
    assert r["predicted_contrast_44_amplitude2_over_volume"] > 0
    assert r["north_energy_density"] > r["transverse_energy_density"]
    assert max(c["quadrature_disagreement"] for c in report["omitted_constraints"]) < 1e-9
    assert report["verdict"]["einstein_constraints"] == "OMITTED_METRIC_AND_SUPPORT_RESPONSE_REQUIRED"
    assert report["verdict"]["triangle_history_map"] == "NOT_DERIVED"
    assert report["verdict"]["source_local_readout"] == "NOT_DERIVED"
    assert report["verdict"]["probability_selection"] == "NOT_DERIVED"


def test_primary_history_keeps_five_tensor_components(report):
    b = np.array(report["history"]["beta_components"])
    assert b.shape == (401, 5)
    # This checks the recorded full shape; no uniaxial projection is fed back.
    assert max(report["history"]["distance_to_uniaxial_cone"]) > 1e-5


def test_probe_gates_and_empty_verdict(report):
    assert report["checks_passed"]
    assert all(report["checks"].values())
    assert set(rt.verdict({})) == set(report["verdict"])
    assert all(rt.verdict({})[key] == "UNRESOLVED" for key in rt.VERDICT_FIELDS)


def test_failed_cli_overwrites_stale_passing_archive(report, monkeypatch, tmp_path):
    failed = copy.deepcopy(report)
    failed["checks"]["injected reciprocity failure"] = False
    failed["checks_passed"] = False
    failed["verdict"] = rt.verdict(failed["checks"])
    monkeypatch.setattr(probe, "run_probe", lambda progress: failed)
    (tmp_path/"probe.md").write_text("STALE PASSING REPORT")
    (tmp_path/"probe.json").write_text('{"checks_passed":true}')
    assert probe.main(["--output-dir", str(tmp_path)]) == 1
    md = (tmp_path/"probe.md").read_text()
    archived = json.loads((tmp_path/"probe.json").read_text())
    assert "STALE" not in md and "UNRESOLVED" in md
    assert "injected reciprocity failure" in md
    assert archived["checks_passed"] is False
    assert all(archived["verdict"][k] == "UNRESOLVED" for k in rt.VERDICT_FIELDS)
