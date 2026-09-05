"""Independent oracles for the preregistered causality/grounding controls."""

import json
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation

from geometrodynamics.bulk import source_readout as readout
from geometrodynamics.bulk.closure_measurement import source_density_on_circle
from geometrodynamics.bulk.closure_grounding import tube_frame_energy, inertia_controls, momentum_integrals
from geometrodynamics.bulk.closure_consistency import run_consistency
from experiments.closure_ledger import source_readout_probe as probe


@pytest.mark.parametrize("choice,variance", [(0, .1), (1, .4)])
def test_odd_readout_mean_does_not_determine_its_law(choice, variance):
    a, b = readout.FIXED_ANALYZER, readout.FUTURE_ANALYZERS[choice]
    x, w = readout.source_circle(a, b)
    r = readout.record_statistics(x, w)
    assert abs(r["mean"]) < 1e-14
    assert r["positive_probability"] == pytest.approx(.5, abs=1e-14)
    assert r["variance"] == pytest.approx(variance, abs=1e-14)
    # Orthogonal-setting closed form, independent of sector summation.
    expected = (1 + np.abs(x @ a) + np.abs(x @ b)) / (1 + 4 / math.pi)
    np.testing.assert_allclose(w * len(w), expected, atol=1e-8, rtol=0)
    np.testing.assert_allclose(w * len(w), source_density_on_circle(a, b, x), atol=1e-14)
    assert r["tail_probability"] == pytest.approx(readout.orthogonal_noiseless_tail()[choice], abs=1e-4)


@pytest.mark.parametrize("beta", [None, 64.])
def test_same_gaussian_kernel_retains_information_while_constant_kernel_erases_it(beta):
    records, constant = [], []
    for b in readout.FUTURE_ANALYZERS:
        x, w = (readout.source_circle(readout.FIXED_ANALYZER, b) if beta is None
                else readout.source_sphere(readout.FIXED_ANALYZER, b, beta))
        records.append(readout.record_statistics(x, w, noise=.15))
        constant.append(float(w @ np.full(len(w), .37)))
    assert records[1]["tail_probability"] - records[0]["tail_probability"] > .2
    assert records[1]["variance"] - records[0]["variance"] > .15
    assert constant == pytest.approx([.37, .37], abs=1e-14)


@pytest.mark.parametrize("choice", [0, 1])
def test_finite_temperature_variance_converges_in_both_coordinates(choice):
    a, b = readout.FIXED_ANALYZER, readout.FUTURE_ANALYZERS[choice]
    values = [readout.record_statistics(*readout.source_sphere(a, b, 64., nz, nphi))["variance"]
              for nz, nphi in ((128, 512), (256, 512), (128, 1024))]
    assert max(abs(values[i] - values[0]) for i in (1, 2)) < 1e-4


J = np.kron(np.eye(3), np.array([[0., 1.], [-1., 0.]]))


def _interaction(z):
    chi, _, phi, _, _, P = z
    x = np.array([np.sin(chi) * np.cos(phi), np.sin(chi) * np.sin(phi), np.cos(chi)])
    return P * (readout.READOUT_AXIS @ x)


@pytest.mark.parametrize("momentum", [0., -.3, .7])
def test_exact_kick_matches_independent_canonical_hamilton_ode(momentum):
    initial = np.array([1.1, .2, .7, -.4, .1, momentum])

    def rhs(t, state):
        gradient = np.empty(6)
        for i in range(6):
            z = state.astype(complex)
            z[i] += 1e-30j
            gradient[i] = _interaction(z).imag / 1e-30
        return J @ gradient

    result = solve_ivp(rhs, (0, .8), initial, rtol=1e-12, atol=1e-14)
    assert result.success
    np.testing.assert_allclose(result.y[:, -1], readout.pointer_kick(initial, .8), atol=1e-12, rtol=0)


def test_kick_is_symplectic_and_finite_momentum_recoils():
    initial = np.array([1.1, .2, .7, -.4, .1, .3])
    h = 1e-6
    jac = np.column_stack([(readout.pointer_kick(initial + h * e, .8)
                            - readout.pointer_kick(initial - h * e, .8)) / (2 * h)
                           for e in np.eye(6)])
    np.testing.assert_allclose(jac.T @ J @ jac, J, atol=3e-10, rtol=0)
    assert np.linalg.det(jac) == pytest.approx(1., abs=3e-10)
    zero = initial.copy()
    zero[5] = 0.
    np.testing.assert_array_equal(readout.pointer_kick(zero)[:4], zero[:4])
    assert abs(readout.pointer_kick(initial)[1] - initial[1]) > 1e-3
    twice = initial.copy()
    twice[5] *= 2
    np.testing.assert_allclose(readout.pointer_kick(twice)[[1, 3]] - twice[[1, 3]],
                               2 * (readout.pointer_kick(initial)[[1, 3]] - initial[[1, 3]]))


def test_finite_pulse_on_existing_periodic_hamiltonian_preserves_both_source_boundaries():
    r = probe.duffing_pointer_demo()
    assert r["baseline_closure_error"] < 1e-9
    assert r["source_endpoint_change"] < 1e-9
    assert r["source_trajectory_change"] < 1e-9
    assert r["energy_range"] < 1e-9
    assert r["pointer_momentum_change"] == 0
    assert r["pointer_shift"] > 1e-6
    # A finite P changes the source: no smooth-preparation no-recoil claim.
    recoil = probe.duffing_pointer_demo(.001)
    assert recoil["source_endpoint_change"] > 1e-6


@pytest.mark.parametrize("angle", [0., .2, 1., 2., math.pi])
def test_scalar_tube_extension_matches_rotation_energy(angle):
    R = np.array([[math.cos(angle), -math.sin(angle), 0.],
                  [math.sin(angle), math.cos(angle), 0.], [0., 0., 1.]])
    for n in (31, 63):
        r = tube_frame_energy(R, n_interior=n)
        expected = 4 * .7 / 1.3 * math.sin(angle / 2) ** 2
        assert r["finite_difference_energy"] == pytest.approx(expected, abs=1e-12)
        assert r["dtn_energy"] == pytest.approx(expected, abs=1e-10)
        assert r["predicted_energy"] == pytest.approx(r["K"] / 2 * math.sin(angle / 2) ** 2)


def test_biinvariance_does_not_select_the_exact_potential():
    R, S, L, T = Rotation.random(4, random_state=31).as_matrix()
    norm = np.linalg.norm(R - S) ** 2
    transformed = np.linalg.norm(L @ R @ T - L @ S @ T) ** 2
    assert transformed == pytest.approx(norm, abs=1e-13)
    assert transformed ** 2 == pytest.approx(norm ** 2, abs=1e-12)
    assert abs(norm ** 2 - norm) > .1


def test_haar_is_weaker_than_round_quadratic_rotor_premises():
    integrals = momentum_integrals()
    assert integrals["gaussian"] == pytest.approx(2 * math.pi, abs=1e-10)
    assert integrals["quartic"] == pytest.approx(math.pi ** 1.5, abs=1e-10)
    for x in ([0, 0, 1], [0, 0, -1], [1, 0, 0], [1, 2, 3]):
        m = inertia_controls(x)
        assert m["variable_volume"] == pytest.approx(m["expected_variable_volume"], abs=1e-14)
        assert m["anisotropic_volume"] == pytest.approx(1., abs=1e-14)
    assert inertia_controls([1, 0, 0])["anisotropy"] > .1
    assert inertia_controls([0, 0, 1])["variable_volume"] == pytest.approx(1.5)
    assert inertia_controls([0, 0, -1])["variable_volume"] == pytest.approx(.5)


def test_cross_round_values_agree_in_all_24_sectors():
    report = run_consistency()
    assert len(report["mass_rows"]) == 24
    for key, value in report["max_residuals"].items():
        tolerance = 1e-7 if key in ("positive_mass", "phase_gradient", "frame_hessian") else 1e-9
        assert value < tolerance, (key, value)


@pytest.mark.parametrize("passed", [False, True])
def test_probe_exit_code_and_saved_evidence_follow_numerical_criteria(monkeypatch, tmp_path, passed):
    report = {"passed": passed, "checks": {"numeric": passed}}
    monkeypatch.setattr(probe, "run_probe", lambda: report)
    monkeypatch.setattr(probe, "render", lambda r: str(r["checks"]))
    assert probe.main(["--output-dir", str(tmp_path)]) == (0 if passed else 1)
    assert json.loads((tmp_path / "probe.json").read_text()) == report
