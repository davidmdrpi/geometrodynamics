"""Independent constraints, sharp bounds, compatibility and failure behavior."""

import copy
import json
import math

import numpy as np
import pytest

from experiments.closure_ledger import scalar_tt_constraints_probe as probe
from geometrodynamics.waves import scalar_tt_constraints as sc


@pytest.fixture(scope="module")
def report():
    return probe.run_probe()


@pytest.mark.parametrize("a,kappa", [(.7, .4), (1., 1.), (2., 3.)])
def test_operator_sign_radius_and_sharp_bounds(a, kappa):
    degrees = np.array([0, 2, 4, 6])
    rho = np.array([.3, -1., 2., -.7])
    u = sc.solve_hamiltonian(rho, degrees, a, kappa)
    np.testing.assert_allclose((3-degrees*(degrees+2))*u/a**2, -kappa*rho/4, atol=1e-14)
    assert u[0] < 0  # Keeping the mean changes the sign of this denominator.
    assert np.linalg.norm(u) <= kappa*a*a*np.linalg.norm(rho)/12
    assert kappa*a*a*np.linalg.norm(rho[1:])/180 <= np.linalg.norm(u[1:])
    assert np.linalg.norm(u[1:]) <= kappa*a*a*np.linalg.norm(rho[1:])/20
    for degree, constant in ((0, 12), (2, 20), (6, 180)):
        actual = abs(sc.solve_hamiltonian([1.], [degree], a, kappa)[0])
        assert actual == pytest.approx(kappa*a*a/constant)


def test_sourced_dipole_is_rejected_even_when_small():
    for amplitude in (1., 1e-16):
        with pytest.raises(ValueError, match="dipole"):
            sc.solve_hamiltonian([amplitude], [1])
    np.testing.assert_array_equal(sc.solve_hamiltonian([0., 1.], [1, 2]), [0., .05])


def test_reconstruction_and_improved_stress(report):
    assert report["checks"]["full_even_reconstruction"]
    assert report["checks"]["independent_improved_stress"]
    assert report["checks"]["two_quadrature_rules"]
    # Parseval: the integrated energy of the free normalized mode is 8 s^2.
    expected_mean = -.2**2*16/(24*2*math.pi**2)
    for r in report["history"]:
        assert r["u_mean"] == pytest.approx(expected_mean, abs=1e-13)


def test_differentiated_constraints_and_kernel_controls(report):
    assert report["checks"]["Hamiltonian_and_momentum_residuals"]
    assert report["checks"]["dipole_and_CKV_compatibility"]
    assert report["checks"]["K_norm_identity"]
    r = report["compatibility_controls"]
    assert r["dipole_rejected"] and r["charged_momentum_rejected"]
    assert r["left_Killing_charge"] < -.01
    assert r["charge_error"] < 1e-10


def test_envelopes_and_same_order_scaling(report):
    assert report["checks"]["spectral_bounds_and_sharp_constants"]
    assert report["checks"]["continuous_envelopes"]
    assert report["checks"]["quadratic_amplitude_scaling"]
    assert report["checks"]["independent_forced_TT_solution"]
    env = report["envelopes"]
    # A rigorous metric size comparison; not a comparison of scalar forces.
    assert env["scalar_metric_inhomogeneous_all_time_lower"] > env["induced_TT_metric_all_time_upper"]
    assert report["coherent_square_certificate"]["powers_times_volume"] == ["1", "27/25", "1/5", "201/175"]
    assert report["coherent_square_certificate"]["exact"]
    assert report["square_power_error"] < 1e-10


def test_support_can_cancel_or_reinforce_as_constraint_data():
    rho = np.array([.3, -.1, .2, .4])
    degrees = np.array([0, 2, 4, 6])
    scalar_u = sc.solve_hamiltonian(rho, degrees)
    for coefficient in (-1., 0., 10.):
        total = sc.solve_hamiltonian(rho+coefficient*rho, degrees)
        np.testing.assert_allclose(total, (1+coefficient)*scalar_u, atol=1e-14)
    # This is freedom in constraint source data, not a constructed fluid.


def test_pointwise_longitudinal_trace_and_radius_scaling():
    g = sc.EvenHarmonicGrid(4, 10, radius=2.)
    rng = np.random.default_rng(sc.SEED)
    w = rng.normal(size=84)
    A = g.longitudinal_coefficients(w)
    np.testing.assert_allclose(np.trace(A, axis1=0, axis2=1), 0., atol=1e-12)
    div = g.divergence_coefficients(A)
    expected = np.array([(4/3)*d @ ((3/g.radius**2-g.lambdas)*w) for d in g.D])
    np.testing.assert_allclose(div, expected, atol=2e-10)


def test_source_shape_and_computed_TT_zeros(report):
    c = report["homogeneous_TT_certificate"]
    assert c["exact"]
    assert c["linear_spatial_curvature"] == "0"
    assert c["momentum_connection_contraction"] == ["0", "0", "0"]
    assert report["degree_one_anticommutator_error"] == 0.
    # No positive lower gate: a computed source-shape zero would be admissible.
    for key in ("initial", "at_time_2"):
        row = report["primary_source_shape"][key]
        assert 0 <= row["uniaxial_distance"] <= row["source_norm"]


def test_success_does_not_claim_a_complete_evolution(report):
    assert report["checks_passed"], report["checks"]
    assert report["verdict"]["complete_scalar_backreaction_bound"] == "NOT_ESTABLISHED_WITHOUT_SUPPORT_AND_EVOLUTION_CLOSURE"
    assert report["verdict"]["readout"] == "NOT_DERIVED"
    assert all(sc.verdict({})[k] == "UNRESOLVED" for k in sc.VERDICT_FIELDS)


def test_failed_cli_overwrites_old_success(report, monkeypatch, tmp_path):
    failed = copy.deepcopy(report)
    failed["checks"]["injected incompatible dipole accepted"] = False
    failed["checks_passed"] = False
    failed["verdict"] = sc.verdict(failed["checks"])
    monkeypatch.setattr(probe, "run_probe", lambda progress: failed)
    (tmp_path/"probe.md").write_text("STALE SUCCESS")
    (tmp_path/"probe.json").write_text('{"checks_passed": true}')
    assert probe.main(["--output-dir", str(tmp_path)]) == 1
    saved = json.loads((tmp_path/"probe.json").read_text())
    assert saved["checks_passed"] is False
    assert all(saved["verdict"][key] == "UNRESOLVED" for key in sc.VERDICT_FIELDS)
    assert "STALE" not in (tmp_path/"probe.md").read_text()
    assert "injected incompatible dipole accepted" in (tmp_path/"probe.md").read_text()
