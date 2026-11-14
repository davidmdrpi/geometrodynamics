"""The wave-mouth-network visualization is drawn from a solved Hamiltonian."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from geometrodynamics.viz import wave_mouth_network as w


@pytest.fixture(scope="module")
def main_run():
    return w.scenario()


@pytest.fixture(scope="module")
def small():
    return w.Network(sites=160, length=40.0)


def test_equations_of_motion_are_hamiltons(small):
    """Finite-difference gradient of H against the vector field."""
    m = small
    rng = np.random.default_rng(0)
    z = 0.2 * rng.standard_normal(m.dim)
    n, k = m.sites, m.mouths
    z[2 * n:2 * n + k] = [8.0, 20.0, 31.0]
    rhs, eps = m.rhs(z), 1e-6
    grad = np.array([(m.energy(z + eps * e) - m.energy(z - eps * e)) / (2 * eps)
                     for e in np.eye(m.dim)])
    assert np.abs(-grad[:n] / m.dx - rhs[n:2 * n]).max() < 1e-6
    assert np.abs(grad[n:2 * n] / m.dx - rhs[:n]).max() < 1e-6
    for pos, mom in ((2 * n, 2 * n + k), (2 * n + 2 * k, 2 * n + 3 * k)):
        assert np.abs(-grad[pos:pos + k] - rhs[mom:mom + k]).max() < 1e-6
        assert np.abs(grad[mom:mom + k] - rhs[pos:pos + k]).max() < 1e-6


def test_support_is_metastable_only_where_energy_is_stored():
    m = w.Network(stored=(4.0, 0.0, 0.0))
    one, zero = np.ones(3), np.zeros(3)
    assert np.allclose(m.support(one), [4.0, 0.0, 0.0])
    assert np.allclose(m.support(zero)[0], 0.0)
    assert np.abs(m.support_force(one)).max() < 1e-12          # equilibrium
    grid = np.linspace(0, 1, 201)
    well = np.array([m.support(np.full(3, q))[0] for q in grid])
    assert well.max() > 4.0 + 0.01                             # a barrier
    assert 0.5 < grid[well.argmax()] < 1.0
    with pytest.raises(ValueError):
        w.Network(barrier=0.1)


def test_ledger_closes(main_run):
    s = w.summary(main_run)
    assert s["energy_max_relative_error"] < 1e-4
    parts = sum(v.sum(axis=1) for v in main_run.ledger.values())
    assert np.allclose(parts, main_run.total)


def test_wave_moves_every_mouth(main_run):
    s = w.summary(main_run)
    assert all(abs(d) > 0.5 for d in s["displacements"])
    assert 0.1 < s["max_speed"] < 1.0                          # subluminal


def test_collapse_emits_a_larger_wave(main_run):
    s = w.summary(main_run)
    assert s["collapse_times"][0] is not None
    assert s["collapse_times"][1:] == [None, None]
    assert s["final_openings"][0] < 0.3
    assert s["wave_energy_final"] > 5 * s["wave_energy_initial"]
    # The gain is paid for by the stored support energy, not created.
    gain = s["wave_energy_final"] - s["wave_energy_initial"]
    assert gain < s["stored_released"] < 4.0 + 1e-9


def test_collapse_follows_the_arrival_of_the_wave(main_run):
    """Causal order: the packet starts 15 away and travels below unit speed."""
    assert w.summary(main_run)["collapse_times"][0] > 15.0


def test_no_wave_no_motion():
    s = w.summary(w.scenario(0.0, duration=20.0))
    assert s["displacements"] == [0.0, 0.0, 0.0]
    assert s["collapse_times"] == [None, None, None]
    assert s["energy_max_relative_error"] == 0.0


def test_below_threshold_moves_mouths_without_collapse():
    s = w.summary(w.scenario(w.SUBTHRESHOLD))
    assert s["collapse_times"] == [None, None, None]
    assert abs(s["displacements"][0]) > 0.5
    assert s["wave_energy_final"] < s["wave_energy_initial"]


def test_without_stored_energy_the_wave_cannot_grow():
    s = w.summary(w.scenario(stored=(0.0, 0.0, 0.0)))
    assert s["collapse_times"] == [None, None, None]
    assert s["wave_energy_max"] <= s["wave_energy_initial"] * (1 + 1e-6)


def test_weak_push_is_second_order_in_amplitude():
    d = [w.summary(w.scenario(a, duration=40.0))["displacements"][0]
         for a in (0.02, 0.04)]
    assert np.log2(d[1] / d[0]) == pytest.approx(2.0, abs=0.15)


def test_second_order_convergence_against_an_independent_integrator(small):
    m = small
    z0 = m.start(positions=(12.0, 22.0, 31.0), amplitude=0.3, centre=4.0)
    ref = solve_ivp(lambda t, z: m.rhs(z), (0, 12), z0, method="DOP853",
                    rtol=1e-10, atol=1e-12).y[:, -1]
    err = [np.linalg.norm(w.evolve(m, z0, 12.0, h=h, sample=int(round(12 / h))
                                   ).states[-1] - ref) for h in (0.01, 0.005)]
    assert err[1] < 1e-3 * np.linalg.norm(ref)
    assert err[0] / err[1] == pytest.approx(4.0, rel=0.1)


def test_render_writes_a_figure(tmp_path, small):
    hist = w.evolve(small, small.start(positions=(12.0, 22.0, 31.0),
                                       amplitude=0.3, centre=4.0), 6.0)
    out = tmp_path / "figure.png"
    summary = w.render(hist, out)
    assert out.stat().st_size > 10_000
    assert summary["energy_max_relative_error"] < 1e-4
