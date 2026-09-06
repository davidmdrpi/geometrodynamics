"""Independent mechanics, geometry and preparation tests for pointer spread."""

from dataclasses import replace
import json
import math

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from geometrodynamics.bulk import pointer_spread as spread
from geometrodynamics.bulk.closure_current import minimal_rotation_lift
from geometrodynamics.bulk.mouth_spin_frame import qmul
from experiments.closure_ledger import pointer_spread_probe as probe


@pytest.mark.parametrize("P", [-.3, 0., .2])
def test_split_flow_matches_independent_hamilton_and_transport_ode(P):
    x = np.array([.3, .4, math.sqrt(.75)])
    p = np.array([.2, -.1, 0.])
    p -= (p @ x) * x
    reference = probe.ode_reference(x, p, P)
    errors = []
    for n in (32, 64, 128):
        result = spread.evolve(x[None], p[None], [P], steps=n)
        errors.append(max(float(np.max(abs(result[k][0] - reference[k]))) for k in
                          ("x_read", "p_read", "record_shift", "x_final", "p_final", "transport")))
        assert max(result["diagnostics"].values()) < 1e-10
    assert errors[-1] < 1e-3
    if P != 0:
        assert 3 < errors[0] / errors[1] < 5
        assert 3 < errors[1] / errors[2] < 5
    assert reference["work_residual"] < 1e-10
    assert reference["P_change"] == 0


def test_free_source_geodesic_and_transport_agree_with_existing_lift():
    x, p = spread.source_nodes(5)
    C = np.tile([1., 0., 0., 0.], (len(x), 1))
    end, _, transport = spread.geodesic_drift(x, p, C, 1.)
    for i in range(len(x)):
        np.testing.assert_allclose(transport[i], minimal_rotation_lift(x[i], end[i]), atol=1e-12)
    paths = spread.evolve(x, p, np.zeros(len(x)))
    free = spread.free_path(x, p)
    for key in ("x_final", "p_final", "transport"):
        np.testing.assert_allclose(paths[key], free[key], atol=1e-12)


def test_zero_width_reduces_to_previous_triangle_and_record():
    result = probe.geometry_controls()
    assert max(result.values()) < 1e-11


def test_augmented_holonomy_retains_path_transport_even_at_identical_endpoints():
    x = np.array([[.2, .3, math.sqrt(.87)]])
    identity = np.array([[1., 0., 0., 0.]])
    # A closed source excursion can carry holonomy while leaving x fixed.
    C = np.column_stack(([math.cos(.4)], math.sin(.4) * x))
    plain = {"x0": x, "x_final": x, "transport": identity}
    loop = {**plain, "transport": C}
    a, b = spread.FIXED_ANALYZER, spread.FUTURE_ANALYZERS[0]
    original = spread.history_weights(plain, a, b)
    changed = spread.history_weights(loop, a, b)
    assert np.max(abs(changed["energies"] - original["energies"])) > .01
    expected = []
    for sa in (-1, 1):
        for sb in (-1, 1):
            u, w = sa * a, -sb * b
            G = qmul(minimal_rotation_lift(w, x[0]), qmul(minimal_rotation_lift(u, w),
                qmul(minimal_rotation_lift(x[0], u), C[0])))
            expected.append(.5 * (1 - G[0] ** 2))
    np.testing.assert_allclose(changed["energies"][0], expected, atol=1e-14)
    assert changed["loop_axis_residual"] < 1e-14


def test_complete_instrument_rotation_covariance_gives_blind_axis_control():
    # Rotate every source state by the symmetry exchanging the future settings.
    # With m=a this also leaves the instrument invariant; records and masses match.
    R = Rotation.from_rotvec([0., 0., math.pi / 2]).as_matrix()
    x, p = spread.source_nodes(5)
    P = np.linspace(-.2, .2, len(x))
    model = replace(spread.PointerModel(), axis=(0., 0., 1.))
    first = spread.evolve(x, p, P, model)
    second = spread.evolve(x @ R.T, p @ R.T, P, model)
    np.testing.assert_allclose(first["record_shift"], second["record_shift"], atol=1e-12)
    w0 = spread.history_weights(first, spread.FIXED_ANALYZER, spread.FUTURE_ANALYZERS[0])
    w1 = spread.history_weights(second, spread.FIXED_ANALYZER, spread.FUTURE_ANALYZERS[1])
    np.testing.assert_allclose(w0["weight"], w1["weight"], atol=1e-11)
    assert np.var(first["record_shift"]) > .1


def test_fixed_preparation_and_joint_conditioning_are_distinct_normalized_laws():
    pn, prior = spread.pointer_nodes(.2, 12)
    weight = (1 + pn[:, None] ** 2) * np.array([[.2, .3, .5]])
    measures = spread.normalized_measures(weight, prior, np.array([.2, .3, .5]))
    np.testing.assert_allclose(measures["fixed"].sum(axis=1), prior, atol=1e-14)
    assert np.max(abs(measures["joint"].sum(axis=1) - prior)) > 1e-3
    for key in ("fixed", "joint", "frozen"):
        assert measures[key].sum() == pytest.approx(1., abs=1e-14)
        assert np.sum(measures[key] * .37) == pytest.approx(.37, abs=1e-14)
    assert prior @ pn == pytest.approx(0., abs=1e-14)
    assert prior @ (pn ** 2) == pytest.approx(.04, abs=1e-14)


def test_early_record_is_independent_of_the_time_of_later_free_evolution():
    x, p = spread.source_nodes(4)
    P = np.full(len(x), .2)
    first = spread.evolve(x, p, P)
    later = spread.evolve(x, p, P, replace(spread.PointerModel(), future_time=3.))
    np.testing.assert_array_equal(first["record_shift"], later["record_shift"])
    assert np.max(abs(first["x_final"] - later["x_final"])) > .01


def test_nonzero_spread_recoils_preserves_prepared_marginal_and_keeps_parity():
    r = spread.simulate_spread(.1, power=6, hermite=8, steps=32)
    assert r["diagnostics"]["reference_momentum_recoil_rms"] > .01
    assert r["diagnostics"]["reference_position_recoil_rms"] > .001
    assert r["diagnostics"]["record_antipodal_residual"] < 1e-10
    assert r["diagnostics"]["weight_parity_residual"] < 1e-10
    for row in r["choices"]:
        assert row["fixed"]["P_mean"] == pytest.approx(0., abs=1e-12)
        assert row["fixed"]["P_variance"] == pytest.approx(.01, abs=1e-12)
        assert row["fixed"]["positive_probability"] == pytest.approx(.5, abs=1e-12)
        assert row["fixed"]["mean"] == pytest.approx(0., abs=1e-12)


def test_finite_pulse_coordinate_flow_is_symplectic_at_finite_pointer_momentum():
    initial = np.array([1.1, .2, .7, -.1, .1, .3])

    def flow(z):
        chi, pc, phi, pp, Q, P = z
        x = np.array([np.sin(chi) * np.cos(phi), np.sin(chi) * np.sin(phi), np.cos(chi)])
        ec = np.array([np.cos(chi) * np.cos(phi), np.cos(chi) * np.sin(phi), -np.sin(chi)])
        ep = np.array([-np.sin(phi), np.cos(phi), 0.])
        p = pc * ec + pp / np.sin(chi) * ep
        r = spread.evolve(x[None], p[None], [P], steps=32)
        x, p = r["x_read"][0], r["p_read"][0]
        chi, phi = np.arccos(x[2]), np.arctan2(x[1], x[0])
        ec = np.array([np.cos(chi) * np.cos(phi), np.cos(chi) * np.sin(phi), -np.sin(chi)])
        ep = np.array([-np.sin(phi), np.cos(phi), 0.])
        return np.array([chi, p @ ec, phi, np.sin(chi) * (p @ ep), Q + r["record_shift"][0], P])

    h = 2e-6
    derivative = np.column_stack([(flow(initial + h * e) - flow(initial - h * e)) / (2 * h) for e in np.eye(6)])
    J = np.kron(np.eye(3), [[0., 1.], [-1., 0.]])
    np.testing.assert_allclose(derivative.T @ J @ derivative, J, atol=3e-7, rtol=0)


@pytest.mark.parametrize("passed", [False, True])
def test_probe_preserves_failed_criteria_in_exit_code_and_archive(monkeypatch, tmp_path, passed):
    report = {"checks": {"primary": passed}, "passed": passed}
    monkeypatch.setattr(probe, "run_probe", lambda progress=None: report)
    monkeypatch.setattr(probe, "render", lambda r: str(r))
    assert probe.main(["--output-dir", str(tmp_path)]) == (0 if passed else 1)
    assert json.loads((tmp_path / "probe.json").read_text()) == report


def test_antipodal_geodesic_puncture_is_not_given_a_weight():
    with pytest.raises(ValueError, match="puncture"):
        spread.minimal_lift([0., 0., 1.], [0., 0., -1.])
