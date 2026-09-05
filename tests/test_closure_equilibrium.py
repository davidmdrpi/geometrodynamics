"""Independent oracles for the preregistered classical equilibrium model."""

import json
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import erf

from geometrodynamics.bulk import closure_equilibrium as ce
from geometrodynamics.bulk.closure_current import minimal_rotation_lift
from geometrodynamics.bulk.closure_measurement import closed_form_correlation
from geometrodynamics.bulk.mouth_spin_frame import qinv, qmul, spin_frame


def _unit(v):
    return v / np.linalg.norm(v)


def test_restoring_energy_from_actual_quaternion_frame_transport():
    """P1: compose geodesic lifts; compare returned frames, without atan2."""
    rng = np.random.default_rng(281)
    for _ in range(80):
        q0 = _unit(rng.normal(size=4))
        x = spin_frame(q0)[0]
        u, w = [_unit(v) for v in rng.normal(size=(2, 3))]
        G = qmul(minimal_rotation_lift(w, x),
                 qmul(minimal_rotation_lift(u, w), minimal_rotation_lift(x, u)))
        initial = np.asarray(spin_frame(q0))
        returned = np.asarray(spin_frame(qmul(q0, qinv(G))))
        measured = float(np.sum((returned - initial) ** 2) / 16.0)
        assert ce.rotor_potential(x, u, w) == pytest.approx(measured, abs=1e-11)
        # J^2 changes the lift sign, never the classical frame or its energy.
        negative = np.asarray(spin_frame(qmul(q0, qinv(-G))))
        assert np.max(np.abs(negative - returned)) < 1e-12
        # Simultaneous spatial rotations and an arbitrary initial frame.
        R = np.asarray(spin_frame(_unit(rng.normal(size=4))))
        assert ce.rotor_potential(R @ x, R @ u, R @ w) == pytest.approx(measured, abs=1e-11)
        q1 = _unit(rng.normal(size=4))
        f1 = np.asarray(spin_frame(q1))
        f2 = np.asarray(spin_frame(qmul(q1, qinv(G))))
        assert np.sum((f2 - f1) ** 2) / 16.0 == pytest.approx(measured, abs=1e-11)


@pytest.mark.parametrize("c", [-0.8, -0.25, 0.4, 0.85])
@pytest.mark.parametrize("model", ce.MODELS)
def test_closed_form_masses_against_independent_piecewise_quadrature(c, model):
    """P2/P3/P4: split at the actual zeros, including the negative D arc."""
    root = math.acos(-math.sqrt((1.0 + c) / 2.0))
    edges = [0.0, root, 2 * math.pi - root, 2 * math.pi]
    value = sum(quad(lambda t: ce.limiting_density(t, c, model=model), lo, hi,
                     epsabs=1e-11)[0] for lo, hi in zip(edges, edges[1:]))
    assert value == pytest.approx(ce.limiting_mass(c, model=model), rel=1e-11)


@pytest.mark.parametrize("gamma", [math.pi / 4, 1.0, 3 * math.pi / 4])
def test_whole_sphere_integral_approaches_the_predicted_measure(gamma):
    """P2 thresholds frozen in 76ed50e; test the mass as well as ratios."""
    beta = 4096
    target = ce.joint_probabilities(gamma, None)
    actual = ce.joint_probabilities(gamma, beta)
    assert max(abs(actual[k] - target[k]) for k in target) < 2e-3
    for c in (-math.cos(gamma), math.cos(gamma)):
        scaled = math.sqrt(beta / (2 * math.pi)) * ce.canonical_partition(c, beta)
        oracle = ce.limiting_mass(c) / (4 * math.pi)
        assert scaled == pytest.approx(oracle, rel=0.01)
    # Independently refined coordinates: convergence is not the oracle.
    normal = ce.joint_probabilities(gamma, beta, n_normal=256)
    azimuth = ce.joint_probabilities(gamma, beta, n_azimuth=2048)
    assert max(abs(actual[k] - normal[k]) for k in target) < 1e-4
    assert max(abs(actual[k] - azimuth[k]) for k in target) < 1e-4


def test_low_temperature_law_is_the_existing_positive_measure_with_partner_sign():
    for gamma in (0.3, 0.8, 1.0, 1.8, 2.5):
        got = ce.correlation(ce.joint_probabilities(gamma, None))
        assert got == pytest.approx(-float(closed_form_correlation(gamma)), abs=1e-12)
    assert ce.standard_chsh() == pytest.approx(2.142283163168201, abs=1e-12)
    assert abs(ce.standard_chsh() - 2 * math.sqrt(2)) > 0.6


def test_partition_quadrature_against_unreduced_spherical_monte_carlo():
    """P2 independent construction; all sectors and ratio uncertainty measured."""
    mc = ce.monte_carlo_joint(1.0, 16, n_samples=120000)
    reference = ce.joint_probabilities(1.0, 16)
    for sector in reference:
        error = abs(mc["probabilities"][sector] - reference[sector])
        assert error < 6 * mc["standard_errors"][sector] + 1e-5
    assert abs(mc["correlation"] - ce.correlation(reference)) < (
        6 * mc["correlation_standard_error"] + 1e-5)


@pytest.mark.parametrize("beta", [0, 1, 16, 4096])
@pytest.mark.parametrize("c", [-0.8, 0.0, 0.8])
def test_normal_residual_partition_is_a_one_dimensional_gaussian(beta, c):
    a = beta * (1.0 - c * c) / 2.0
    exact = 1.0 if a == 0 else math.sqrt(math.pi) * erf(math.sqrt(a)) / (2 * math.sqrt(a))
    assert ce.canonical_partition(c, beta, model="normal") == pytest.approx(exact, abs=1e-12)
    assert ce.correlation(ce.joint_probabilities(1.0, beta, model="normal")) == pytest.approx(0, abs=1e-12)


def test_same_zero_locus_but_different_restoring_energies_select_different_ensembles():
    u, w = np.array([0.0, 0.0, 1.0]), np.array([math.sin(1.0), 0.0, math.cos(1.0)])
    # Sample the regular circle, excluding both undefined geodesic punctures.
    t = 2 * math.pi * (np.arange(100) + 0.37) / 100
    x = np.column_stack([np.sin(t), np.zeros_like(t), np.cos(t)])
    assert np.max(ce.rotor_potential(x, u, w, model="frame")) < 1e-20
    assert np.max(ce.rotor_potential(x, u, w, model="normal")) < 1e-20
    assert ce.correlation(ce.joint_probabilities(1.0, None, model="normal")) == 0
    assert ce.correlation(ce.joint_probabilities(1.0, None)) < -0.39


def test_reparametrization_with_covariant_stiffness_preserves_finite_temperature_measure():
    rng = np.random.default_rng(282)
    x = rng.normal(size=(1000, 3))
    x /= np.linalg.norm(x, axis=1)[:, None]
    u, w = [_unit(v) for v in rng.normal(size=(2, 3))]
    for epsilon in (-0.4, 0.25, 0.4):
        baseline = ce.rotor_potential(x, u, w, model="normal")
        covariant = ce.rotor_potential(x, u, w, model="normal_covariant", epsilon=epsilon)
        assert np.max(np.abs(baseline - covariant)) < 1e-12
        for c in (-0.6, 0.6):
            a = ce.canonical_partition(c, 64, model="normal")
            b = ce.canonical_partition(c, 64, model="normal_covariant", epsilon=epsilon)
            assert a == pytest.approx(b, abs=1e-12)
    # Leaving the stiffness fixed changes the response, with a closed-form target.
    altered = ce.joint_probabilities(1.0, None, model="normal_rescaled")
    assert ce.correlation(altered) == pytest.approx(-0.038650674985766, abs=1e-12)
    finite = ce.joint_probabilities(1.0, 4096, model="normal_rescaled")
    assert max(abs(finite[k] - altered[k]) for k in altered) < 2e-3


def test_transverse_coarea_covariance_from_finite_differenced_residuals():
    """The product sqrt(a)|grad F| is invariant, including variable g(x)."""
    u, w = np.array([0.0, 0.0, 1.0]), np.array([math.sin(1.0), 0.0, math.cos(1.0)])
    h = 1e-6
    for t in (0.3, 1.4, 2.1, 4.0):
        x = np.array([math.sin(t), 0.0, math.cos(t)])
        xp = math.sqrt(1 - h * h) * x + np.array([0, h, 0])
        xm = math.sqrt(1 - h * h) * x - np.array([0, h, 0])
        Np, _ = ce.triangle_data(xp, u, w)
        Nm, _ = ce.triangle_data(xm, u, w)
        gp, gm, g0 = [1 + 0.25 * q @ (u + w) for q in (xp, xm, x)]
        dF, dGF = (Np - Nm) / (2 * h), (gp * Np - gm * Nm) / (2 * h)
        assert ce.coarea_weight(abs(dF)) == pytest.approx(
            ce.coarea_weight(abs(dGF), 1 / g0 ** 2), rel=1e-10)


@pytest.mark.parametrize("gamma", [math.pi / 4, 1.0, 3 * math.pi / 4])
def test_higher_order_frame_stiffness_preserves_the_low_temperature_limit(gamma):
    p = ce.joint_probabilities(gamma, 4096, quartic=2)
    target = ce.joint_probabilities(gamma, None)
    assert max(abs(p[k] - target[k]) for k in p) < 2e-3
    for c in (-math.cos(gamma), math.cos(gamma)):
        scaled = math.sqrt(4096 / (2 * math.pi)) * ce.canonical_partition(c, 4096, quartic=2)
        assert scaled == pytest.approx(ce.limiting_mass(c) / (4 * math.pi), rel=0.01)
    # The models are not the same away from the low-temperature limit.
    a, b = ce.joint_probabilities(gamma, 16), ce.joint_probabilities(gamma, 16, quartic=2)
    assert max(abs(a[k] - b[k]) for k in a) > 1e-4


def test_geodesic_punctures_carry_no_limiting_mass():
    c = math.cos(1.0)
    root = math.acos(-math.sqrt((1 + c) / 2))
    for eps in (0.1, 0.05, 0.025):
        missing = sum(quad(lambda t: ce.limiting_density(t, c), center - eps,
                           center + eps, points=[center], epsabs=1e-12)[0]
                      for center in (root, 2 * math.pi - root))
        assert missing == pytest.approx(2 * eps * eps, rel=1e-3)


def test_sector_priors_remain_inputs_while_paired_symmetry_preserves_marginals():
    Es = []
    for ratio in (0.5, 1.0, 2.0):
        p = ce.joint_probabilities(1.0, 64, prior_ratio=ratio)
        assert all(v > 0 for v in p.values())
        assert sum(v for (sa, sb), v in p.items() if sa == 1) == pytest.approx(0.5, abs=1e-12)
        assert sum(v for (sa, sb), v in p.items() if sb == 1) == pytest.approx(0.5, abs=1e-12)
        Es.append(ce.correlation(p))
    assert Es[2] - Es[0] > 0.4


@pytest.mark.parametrize("gamma", [0, math.pi, -0.1, float("nan")])
def test_collinear_and_nonphysical_angles_are_not_silently_regularized(gamma):
    with pytest.raises(ValueError):
        ce.joint_probabilities(gamma, 16)


def test_undefined_geodesic_and_invalid_preparations_are_rejected():
    u, w = np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="undefined"):
        ce.triangle_data(-u, u, w)
    with pytest.raises(ValueError, match="prior"):
        ce.joint_probabilities(1.0, 16, prior_ratio=0)
    with pytest.raises(ValueError, match="positive"):
        ce.canonical_partition(0.4, 16, model="normal_rescaled", epsilon=0.5)
    with pytest.raises(ValueError, match="nonnegative"):
        ce.canonical_partition(0.4, -1)


@pytest.mark.parametrize("passed,expected_status", [(False, 1), (True, 0)])
def test_probe_exit_status_uses_actual_checks_and_preserves_failure_evidence(
        monkeypatch, tmp_path, passed, expected_status):
    from experiments.closure_ledger import closure_equilibrium_probe as probe

    # A stale success count must not hide a failed numerical criterion.
    report = {"checks": [{"name": "injected numerical criterion", "passed": passed}],
              "passed": 1, "total": 1}
    monkeypatch.setattr(probe, "run_probe", lambda: report)
    monkeypatch.setattr(probe, "render", lambda r: "injected probe report\n")
    assert probe.main(["--output-dir", str(tmp_path)]) == expected_status
    assert json.loads((tmp_path / "probe.json").read_text())["checks"][0]["passed"] is passed
    assert (tmp_path / "probe.md").read_text() == "injected probe report\n"
