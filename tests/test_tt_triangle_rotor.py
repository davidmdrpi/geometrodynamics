"""Independent action, unreduced dynamics and evidence-gate checks."""

import json
import math

import numpy as np
import pytest

from geometrodynamics.bulk import tt_triangle_rotor as rotor
from experiments.closure_ledger import tt_triangle_rotor_probe as probe


def test_adm_curvature_derivation_has_the_round_anchor_and_tensor_frequency():
    result = rotor.adm_quadratic_derivation()
    assert result["round_anchor"] and result["linear_vanishes"]
    assert result["kinetic_is_frobenius"]
    assert result["omega2_radius2"] == pytest.approx(8., abs=1e-13)
    assert rotor.TensorModel(radius=2.).omega2 == 2.


def test_action_normalization_tracks_volume_and_gravitational_coupling():
    model = rotor.TensorModel(radius=2., kappa=3.)
    assert model.normalization == pytest.approx(16 * math.pi ** 2 / 3)
    J = rotor.stf(np.diag([1., 2., 0.]))
    B = rotor.embedding(.02, [0., 0., 1.])
    S = rotor.normalized_source(J, model)
    difference = rotor.lagrangian(B, np.zeros((3, 3)), S, model) - rotor.lagrangian(B, np.zeros((3, 3)), model=model)
    assert difference == pytest.approx(np.sum(B * J), abs=1e-13)


def test_independent_chart_variation_recovers_both_restricted_equations():
    import sympy as sp
    A, theta, phi = q = sp.symbols("A theta phi", real=True)
    adot, tdot, pdot = vel = sp.symbols("adot tdot pdot", real=True)
    addot, tddot, pddot = acc = sp.symbols("addot tddot pddot", real=True)
    n = sp.Matrix([sp.sin(theta) * sp.cos(phi), sp.sin(theta) * sp.sin(phi), sp.cos(theta)])
    source = sp.Matrix([[sp.Rational(1, 50), sp.Rational(1, 100), 0],
                        [sp.Rational(1, 100), -sp.Rational(3, 100), sp.Rational(1, 200)],
                        [0, sp.Rational(1, 200), sp.Rational(1, 100)]])
    L = adot ** 2 / 3 + A ** 2 * (tdot ** 2 + sp.sin(theta) ** 2 * pdot ** 2) - 8 * A ** 2 / 3 + A * (n.T * source * n)[0]
    equations = []
    for qi, vi in zip(q, vel):
        momentum = sp.diff(L, vi)
        equations.append(sum(sp.diff(momentum, z) * dz for z, dz in zip(q + vel, vel + acc)) - sp.diff(L, qi))
    th, ph, td, pd = 1.1, .7, .2, -.1
    axis = np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)])
    et = np.array([math.cos(th) * math.cos(ph), math.cos(th) * math.sin(ph), -math.sin(th)])
    ep = np.array([-math.sin(ph), math.cos(ph), 0.])
    v = td * et + math.sin(th) * pd * ep
    aa, nd = rotor.restricted_acceleration(.02, .003, axis, v, np.array(source, dtype=float))
    ta = et @ nd + math.sin(th) * math.cos(th) * pd ** 2
    pa = (ep @ nd - 2 * math.cos(th) * td * pd) / math.sin(th)
    substitutions = dict(zip(q + vel + acc, (.02, th, ph, .003, td, pd, aa, ta, pa)))
    assert max(abs(float(e.subs(substitutions))) for e in equations) < 1e-12


def test_matrix_finite_differences_and_noncommuting_adm_kinetic_agree():
    result = probe.algebra_controls()
    assert result["metric_finite_difference"] < 1e-7
    assert result["adm_kinetic"] < 1e-10
    assert result["normal_formula"] < 1e-10
    assert result["normal_basis"] < 1e-10
    assert result["covariance"] < 1e-10


def test_biaxial_projector_is_orthogonal_to_every_tangent_velocity():
    n = np.array([.6, 0., .8])
    tangent = np.array([.8, 0., -.6])
    M = rotor.stf(np.array([[.1, .3, .4], [.3, .2, .1], [.4, .1, -.3]]))
    N = rotor.normal_projection(M, n)
    assert abs(np.sum(N * rotor.field_velocity(.03, .2, n, tangent))) < 1e-14
    assert abs(np.trace(N)) < 1e-14
    np.testing.assert_allclose(N @ n, 0., atol=1e-14)
    np.testing.assert_allclose(sum(rotor.tensor_components(M, n).values()), M, atol=1e-14)


@pytest.mark.parametrize("amplitude", [-.02, 0., .02])
def test_cone_minimization_includes_negative_amplitude_and_the_apex(amplitude):
    B = rotor.embedding(amplitude, [0., 0., 1.])
    result = rotor.nearest_uniaxial(B)
    assert result["distance"] < 1e-14
    assert result["amplitude"] == pytest.approx(amplitude, abs=1e-14)
    assert result["axis_defined"] == (amplitude != 0)


def test_cone_minimum_is_global_over_axes_not_a_fixed_axis_comparison():
    rng = np.random.default_rng(431)
    for _ in range(8):
        B = rotor.stf(rng.normal(size=(3, 3)))
        result = rotor.nearest_uniaxial(B)
        assert result["distance"] == pytest.approx(result["eigenvalue_distance"], abs=1e-13)
        axes = rng.normal(size=(1024, 3))
        axes /= np.linalg.norm(axes, axis=1)[:, None]
        amplitudes = 1.5 * np.einsum("ni,ij,nj->n", axes, B, axes)
        candidates = amplitudes[:, None, None] * (axes[:, :, None] * axes[:, None, :] - np.eye(3) / 3)
        assert result["distance"] <= np.linalg.norm(candidates - B, axis=(1, 2)).min() + 1e-13


def test_director_is_antipodally_identified_and_not_defined_at_zero_amplitude():
    n, v = np.array([0., 0., 1.]), np.array([.2, .1, 0.])
    np.testing.assert_array_equal(rotor.embedding(.01, n), rotor.embedding(.01, -n))
    np.testing.assert_array_equal(rotor.field_velocity(.01, .002, n, v), rotor.field_velocity(.01, .002, -n, -v))
    with pytest.raises(ValueError, match="apex"):
        rotor.restricted_acceleration(0., .01, n, v)
    assert np.linalg.matrix_rank(rotor.pullback_metric(0., 1.)) == 1


@pytest.mark.parametrize("speed", [0., .4, .8])
def test_unconstrained_tensor_ode_agrees_with_exact_flow_and_detects_departure(speed):
    n, v = np.array([0., 0., 1.]), np.array([speed, 0., 0.])
    B0, Bd0 = rotor.embedding(.01, n), rotor.field_velocity(.01, 0., n, v)
    times = np.linspace(0, .25, 31)
    B, Bd = rotor.exact_free_flow(B0, Bd0, times)
    reference, reference_dot = probe.ode_flow(B0, Bd0, times)
    np.testing.assert_allclose(B, reference, atol=1e-11, rtol=0)
    np.testing.assert_allclose(Bd, reference_dot, atol=1e-11, rtol=0)
    distance = rotor.nearest_uniaxial(B[-1])["distance"]
    if speed:
        assert distance > 1e-6
    else:
        assert distance < 1e-14


def test_departure_is_quadratic_in_time_and_linear_in_field_amplitude():
    primary, half = probe.free_motion(), probe.free_motion(amplitude=.005)
    assert all(3.8 < r < 4.2 for r in primary["halving_ratios"])
    assert primary["distances_over_t2"][-1] == pytest.approx(primary["predicted_coefficient"], rel=.01)
    np.testing.assert_allclose(half["distances"], .5 * np.array(primary["distances"]), atol=1e-14)


def test_fixed_director_with_free_amplitude_is_invariant_without_freezing_shape_by_force():
    n = np.array([.6, 0., .8])
    B0 = rotor.embedding(.01, n)
    Bd0 = rotor.field_velocity(.01, -.004, n, np.zeros(3))
    times = np.linspace(0, .2, 41)
    B, _ = rotor.exact_free_flow(B0, Bd0, times)
    assert max(rotor.nearest_uniaxial(b)["distance"] for b in B) < 1e-14
    # Freezing A as well would require a nonzero radial source.
    assert np.linalg.norm(rotor.TensorModel().omega2 * B0) > .01


def test_manufactured_rotation_is_a_forced_control_not_a_no_go_for_all_sources():
    result = probe.dynamics_controls()
    assert result["manufactured_ode_error"] < 1e-10
    assert result["manufactured_full_equation_residual"] < 1e-10
    assert result["manufactured_min_radial_force"] > .01
    assert result["manufactured_min_normal_force"] > 1e-4


def test_free_energy_and_rotation_charge_are_conserved_for_generic_biaxial_fields():
    rng = np.random.default_rng(456)
    B0, Bd0 = [rotor.stf(rng.normal(size=(3, 3))) for _ in range(2)]
    B, Bd = rotor.exact_free_flow(B0, Bd0, np.linspace(0, 2, 151))
    np.testing.assert_allclose(rotor.free_energy(B, Bd), rotor.free_energy(B0, Bd0), atol=1e-12)
    np.testing.assert_allclose(rotor.rotational_charge(B, Bd),
                               np.broadcast_to(rotor.rotational_charge(B0, Bd0), B.shape), atol=1e-12)


def test_repository_response_kernel_matches_the_normalized_source_ode():
    result = probe.inherited_response_control()
    assert result["source_variation_identity_residual"] < 1e-13
    assert result["errors_513_1025"][1] < result["errors_513_1025"][0] < 1e-8


@pytest.mark.parametrize("speed", [.2, .4, .8])
def test_post_review_exact_spectrum_and_distance_hold_over_a_full_period(speed):
    times = np.linspace(0, 2 * math.pi / math.sqrt(8), 301)
    for amplitude in (-.01, .01):
        B0 = rotor.embedding(amplitude, [0., 0., 1.])
        Bd0 = rotor.field_velocity(amplitude, 0., [0., 0., 1.], [speed, 0., 0.])
        B, _ = rotor.exact_free_flow(B0, Bd0, times)
        closed = rotor.closed_form_free_spectrum(times, amplitude, speed)
        np.testing.assert_allclose(closed["eigenvalues"], np.linalg.eigvalsh(B), atol=1e-14)
        numerical = [rotor.nearest_uniaxial(b)["distance"] for b in B]
        np.testing.assert_allclose(closed["distance"], numerical, atol=1e-14)


def test_post_review_periodic_returns_are_isolated_not_invariant_motion():
    omega = math.sqrt(8)
    times = np.arange(7) * math.pi / (2 * omega)
    distances = rotor.closed_form_free_spectrum(times)["distance"]
    assert np.max(distances[::2]) < 1e-15
    assert np.min(distances[1::2]) > 1e-4


def test_post_review_eigenline_advances_monotonically_despite_the_fitted_axis_switch():
    omega = math.sqrt(8)
    times = np.linspace(0, 2 * math.pi / omega, 301)
    line = rotor.continuous_eigenline(times)
    assert np.min(np.diff(line["angle"])) > 0
    assert np.min(line["angular_speed"]) > 0
    assert line["angle"][-1] == pytest.approx(math.pi, abs=1e-14)
    # In particular it passes 45 degrees, rather than reversing there.
    assert rotor.continuous_eigenline(np.array([.75 * math.pi / omega]))["angle"][0] > math.pi / 4
    h = 1e-6
    numerical_rate = (rotor.continuous_eigenline(times + h)["angle"]
                      - rotor.continuous_eigenline(times - h)["angle"]) / (2 * h)
    np.testing.assert_allclose(numerical_rate, line["angular_speed"], atol=1e-8)


def test_post_review_nearest_axis_has_two_distinct_optima_at_the_switch():
    result = probe.post_review_controls()
    assert all(result["checks"].values())
    assert not result["nearest_axis_defined_at_quarter_period"]
    assert result["nearest_amplitude_changes_sign"]
    assert max(result["two_optima_distance_errors"]) < 1e-14


def test_failed_or_uncertified_numerics_cannot_be_reported_as_an_obstruction():
    assert rotor.physical_verdict({"required": False}, True) == "UNRESOLVED"
    assert rotor.physical_verdict({"required": True}, False) == "UNRESOLVED"
    assert rotor.physical_verdict({}, True) == "UNRESOLVED"
    assert rotor.physical_verdict({"required": True}, True) == "FREE_ROTATING_UNIAXIAL_TT_FAMILY_NOT_INVARIANT"


@pytest.mark.parametrize("passed", [False, True])
def test_probe_exit_and_archive_keep_failed_checks(monkeypatch, tmp_path, passed):
    report = {"checks": {"required": passed}, "passed": passed}
    monkeypatch.setattr(probe, "run_probe", lambda progress=None: report)
    monkeypatch.setattr(probe, "render", lambda report: str(report))
    assert probe.main(["--output-dir", str(tmp_path)]) == (0 if passed else 1)
    assert json.loads((tmp_path / "probe.json").read_text()) == report
