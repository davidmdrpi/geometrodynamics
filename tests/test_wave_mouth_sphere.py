"""The solved S² model behind the sphere-with-Alice-handles drawing."""
import numpy as np
import pytest

from geometrodynamics.viz import wave_mouth_sphere as v
from geometrodynamics.viz.klein_handles import euler_characteristic, is_orientable
from geometrodynamics.waves import wave_mouth_sphere as m


def rotate(x, axis, angle):
    x = np.asarray(x)
    return (x * np.cos(angle) + np.cross(axis, x) * np.sin(angle)
            + axis * (axis @ x) * (1 - np.cos(angle)))


@pytest.fixture(scope="module")
def main_run():
    return m.scenario()


def perturbed():
    model = m.WaveMouthSphere()
    s = model.start()
    model.drift(s, 2.9)
    s.q[:] = [[-0.3, 0.05], [0.02, -0.04]]
    s.pq[:] = [[0.1, -0.2], [0.05, 0.03]]
    s.angular[:] = [[0.3, -0.2, 0.0], [0.0, 0.1, 0.4]]
    return model, s


def test_kick_is_minus_the_gradient_of_the_potential():
    model, s = perturbed()
    before = s.copy()
    model.kick(s, 1e-6)
    force_q = (s.pq - before.pq) / 1e-6
    torque = (s.angular - before.angular) / 1e-6
    for h in range(2):
        for k in range(2):
            plus, minus = before.copy(), before.copy()
            plus.q[h, k] += 1e-6
            minus.q[h, k] -= 1e-6
            grad = (model.energy(plus) - model.energy(minus)) / 2e-6
            assert force_q[h, k] == pytest.approx(-grad, abs=1e-7)
        for i, axis in enumerate(np.eye(3)):
            plus, minus = before.copy(), before.copy()
            plus.x[h] = rotate(before.x[h], axis, 1e-6)
            minus.x[h] = rotate(before.x[h], axis, -1e-6)
            grad = (model.energy(plus) - model.energy(minus)) / 2e-6
            assert torque[h, i] == pytest.approx(-grad, abs=1e-7)


def test_relativistic_pair_drift_is_a_geodesic_at_the_right_speed():
    model = m.WaveMouthSphere()
    s = model.start(0.0)
    s.angular[0] = [0.0, 8.0, 0.0]                 # |L| > M: still subluminal
    model.drift(s, 0.5)
    speed = 8.0 / np.hypot(model.par.mass, 8.0)
    assert speed < 1
    np.testing.assert_allclose(s.x[0], [np.sin(speed * 0.5), 0, np.cos(speed * 0.5)],
                               atol=1e-14)


def test_support_is_metastable_only_where_loaded():
    model = m.WaveMouthSphere()
    q = np.zeros((2, 2))
    np.testing.assert_allclose(model.support(q), model.stored)
    collapsed = model.support(q - 1)
    assert collapsed[0, 0] == pytest.approx(0, abs=1e-15)
    assert (collapsed.ravel()[1:] > 0).all()        # unloaded: plain well
    np.testing.assert_allclose(model.support_force(q), 0, atol=1e-15)
    top = model.barrier_top()
    assert -0.5 < top < 0
    assert model.support(np.full((2, 2), top))[0, 0] > model.stored[0, 0]


def test_ledger_closes_and_angular_momentum_is_conserved(main_run):
    model, d = main_run
    s = v.summary(model, d)
    assert s["energy_max_relative_error"] < 2e-4
    assert s["angular_momentum_drift"] < 1e-11


def test_step_is_reversible():
    model, s = perturbed()
    before = s.copy()
    for _ in range(20):
        model.step(s, 0.01)
    for _ in range(20):
        model.step(s, -0.01)
    for key in ("a", "p", "x", "angular", "q", "pq"):
        np.testing.assert_allclose(getattr(s, key), getattr(before, key), atol=1e-12)


def test_wave_moves_both_pairs_and_tips_only_the_loaded_throat(main_run):
    model, d = main_run
    s = v.summary(model, d)
    assert all(x > 0.5 for x in s["pair_displacements"])
    assert 0.05 < s["max_speed"] < 1.0
    assert s["transition_times"][0] is not None
    assert s["transition_times"][1:] == [None, None, None]
    assert s["transition_times"][0] > 2.0          # after the wave has travelled
    assert s["field_energy_final"] > 3 * s["field_energy_initial"]
    gain = s["field_energy_final"] - s["field_energy_initial"]
    assert gain < s["stored_released"] < 2.0 + 1e-9
    assert s["field_energy_final_l0"] < 0.2 * s["field_energy_final"]
    assert s["min_pair_separation"] > 2 * v.COLLAR


def test_controls():
    model, d = m.scenario(0.0, duration=5.0)
    s = v.summary(model, d)
    assert s["pair_displacements"] == [0.0, 0.0]
    assert s["transition_times"] == [None] * 4
    assert s["energy_max_relative_error"] == 0.0
    for kw in (dict(amplitude=m.SUBTHRESHOLD),
               dict(stored=((0.0, 0.0), (0.0, 0.0)))):
        model, d = m.scenario(**kw)
        s = v.summary(model, d)
        assert s["transition_times"] == [None] * 4
        assert s["pair_displacements"][0] > 0.1    # the wave still pushes
        assert s["field_energy_final"] < 1.1 * s["field_energy_initial"]


def test_energy_error_is_second_order():
    _, coarse = m.scenario(duration=8.0, h=0.01, sample=10)
    _, fine = m.scenario(duration=8.0, h=0.005, sample=20)
    err = [np.abs(d["total"] - d["total"][0]).max() for d in (coarse, fine)]
    assert err[0] / err[1] == pytest.approx(4.0, rel=0.25)


def test_drawn_surface_is_the_two_klein_handle_sphere_at_the_solved_axes(main_run):
    model, d = main_run
    scene = v.SphereScene(model, n_sphere=700, n_ring=16)
    for j in (0, len(d["t"]) // 2, -1):
        surf, geom, drawn, field = scene.state(d["a"][j], d["x"][j], 0.2)
        assert euler_characteristic(surf.tris) == -2
        assert not is_orientable(surf.tris)
        for h in range(2):
            for sign, rim in zip((1, -1), surf.rims[h]):
                np.testing.assert_allclose(surf.points[rim].mean(axis=0),
                                           sign * surf.neck.ca * d["x"][j][h],
                                           atol=1e-12)
        # The colour is the solved field at the drawn sphere vertices.
        np.testing.assert_allclose(
            field, model.sample(d["a"][j], surf.points[geom.sphere]))
        assert np.abs(np.linalg.norm(drawn[geom.sphere], axis=1) - 1).max() <= 0.1 + 1e-12


def test_figure_and_animation_render(main_run, tmp_path):
    model, d = main_run
    scene = v.SphereScene(model, n_sphere=500, n_ring=12)
    out = tmp_path / "sphere.png"
    s = v.render(model, d, out, times=(0.0, 6.0), scene=scene)
    assert out.stat().st_size > 20_000 and len(s["frame_times"]) == 2
    gif = tmp_path / "sphere.gif"
    assert v.animate(model, d, gif, frames=3, scene=scene) == 3
    assert gif.stat().st_size > 5_000


# ---- outward energy flux through a fixed circle on the sphere ---------------

@pytest.fixture(scope="module")
def cap_flux():
    return m.cap_flux_ledger()


def test_cap_quadrature_is_exact_on_the_whole_sphere():
    model = m.WaveMouthSphere()
    s = model.start()
    model.drift(s, 1.7)
    whole = m.CapFlux(model, radius=np.pi - 1e-9)
    assert whole.energy(s.a, s.p) == pytest.approx(model.wave_energy(s.a, s.p),
                                                   rel=1e-12)


def test_free_field_obeys_the_continuity_equation_through_the_cap_edge():
    model = m.WaveMouthSphere()
    s = model.start()
    model.drift(s, 1.7)
    cap = m.CapFlux(model)
    plus, minus = s.copy(), s.copy()
    model.drift(plus, 1e-5)
    model.drift(minus, -1e-5)
    rate = (cap.energy(plus.a, plus.p) - cap.energy(minus.a, minus.p)) / 2e-5
    assert abs(rate) > 0.1
    assert rate == pytest.approx(-cap.flux(s.a, s.p), rel=1e-7)


def test_cap_flux_balances_the_energy_inside(cap_flux):
    f = cap_flux
    assert f["balance_max_residual"] < 5e-4
    assert f["net_outward_final"] == pytest.approx(
        f["inside_initial"] - f["inside_final"], abs=5e-4)
    assert f["loaded_mouth_margin"] > 0.44       # three footprint widths
    assert f["other_mouth_margin"] > 0.44


def test_packet_enters_the_cap_then_more_energy_leaves(cap_flux):
    f = cap_flux
    assert f["net_outward_minimum"] < -0.5
    assert f["time"][f["net_outward"].argmin()] < 5.0
    assert f["net_outward_final"] > 1.0
    # A closed sphere returns what it emits: the curve is not monotonic.
    assert f["net_flux_fall"] > 1.0
    assert np.diff(f["net_outward"]).min() < 0 < np.diff(f["net_outward"]).max()


def test_cap_absorbs_when_there_is_no_transition():
    for kw in (dict(amplitude=m.SUBTHRESHOLD),
               dict(stored=((0.0, 0.0), (0.0, 0.0)))):
        f = m.cap_flux_ledger(**kw)
        assert f["net_outward_final"] < 0
        assert f["balance_max_residual"] < 1e-4
    assert m.cap_flux_ledger(amplitude=0.0, duration=2.0)["net_outward_final"] == 0


def test_archived_sphere_flux():
    import json
    from pathlib import Path
    path = (Path(__file__).resolve().parents[1] / "docs" / "images"
            / "wave_mouth_sphere" / "flux.json")
    f = json.loads(path.read_text())
    full, half = f["wave_mouth_sphere"], f["half_step"]
    assert full["duration"] == m.CAP_WINDOW
    assert full["net_outward_final"] == pytest.approx(
        half["net_outward_final"], abs=5e-4)
    assert half["balance_max_residual"] < full["balance_max_residual"]


def test_sphere_flux_figure_renders(cap_flux, tmp_path):
    out = tmp_path / "flux.png"
    v.render_flux(cap_flux, out, transition=3.6)
    assert out.stat().st_size > 10_000


def test_gross_flows_are_split_pointwise_on_the_circle(cap_flux):
    """Review of 324c2e4: inflow on one arc must not cancel outflow on another."""
    f = cap_flux
    assert f["gross_outward"] - f["gross_inward"] == pytest.approx(
        f["net_outward_final"], abs=1e-12)
    assert f["gross_outward"] == pytest.approx(3.319, abs=2e-3)
    assert f["gross_inward"] == pytest.approx(2.152, abs=2e-3)
    # The sign-split of the integrated net flux is smaller: it cancels.
    assert f["net_flux_rise"] < f["gross_outward"] - 0.5
    assert f["net_flux_fall"] < f["gross_inward"] - 0.5


def test_simultaneous_inflow_and_outflow_do_not_cancel():
    """A free wave crossing the cap enters on one arc and leaves on another."""
    model = m.WaveMouthSphere(m.replace(m.PARAMETERS, coupling=0.0))
    cap = m.CapFlux(model, centre=(1.0, 0.0, 0.0), radius=0.6)
    s = model.start()
    model.drift(s, 1.2)                  # the front is sweeping across this cap
    density = cap.flux_density(s.a, s.p)
    out, back = np.clip(density, 0, None).sum(), np.clip(-density, 0, None).sum()
    assert out > 0.01 and back > 0.01
    assert abs(density.sum()) < out + back - 0.02
    assert density.sum() == pytest.approx(cap.flux(s.a, s.p))
    # The net flux is band-limited and so independent of the edge resolution.
    fine = m.CapFlux(model, centre=(1.0, 0.0, 0.0), radius=0.6, edge_points=416)
    assert fine.flux(s.a, s.p) == pytest.approx(cap.flux(s.a, s.p), rel=1e-11)


def test_footprint_tail_work_restores_second_order_balance():
    import json
    from pathlib import Path
    path = (Path(__file__).resolve().parents[1] / "docs" / "images"
            / "wave_mouth_sphere" / "flux.json")
    f = json.loads(path.read_text())
    res = [f[k]["balance_max_residual_with_tails"]
           for k in ("wave_mouth_sphere", "half_step", "quarter_step")]
    assert res[0] / res[1] == pytest.approx(4.0, rel=0.02)
    assert res[1] / res[2] == pytest.approx(4.0, rel=0.02)
    assert f["wave_mouth_sphere"]["tail_work_final"] == pytest.approx(5.11e-5, rel=0.02)
    # Gross flows do not depend on the edge resolution.
    for key in ("gross_outward", "gross_inward"):
        assert f["edge_416"][key] == pytest.approx(f["wave_mouth_sphere"][key],
                                                   abs=1e-5)


def test_readme_code_fences_are_paired():
    from pathlib import Path
    text = (Path(__file__).resolve().parents[1] / "README.md").read_text()
    fences = [line for line in text.splitlines() if line.startswith("```")]
    assert len(fences) % 2 == 0
    start = text.index("## A visualization drawn from a solved model")
    assert sum(l.startswith("```") for l in text[:start].splitlines()) % 2 == 0
