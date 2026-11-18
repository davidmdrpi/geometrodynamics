"""Alice handle ports: the separable solver, its ledgers, and the archive."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.bulk import alice_handle_ports as ah

PROBE = (Path(__file__).resolve().parents[1] / "experiments" / "closure_ledger"
         / "alice_handle_ports_probe.py")
ARCHIVE = (Path(__file__).resolve().parents[1] / "experiments" / "closure_ledger"
           / "runs" / "20261118_alice_handle_ports" / "probe.json")


@pytest.fixture(scope="module")
def loop():
    return ah.Loop(0.5, 300)


def test_loop_radius_is_continuous_at_both_rims(loop):
    rho = loop.rho_exact()
    assert rho[loop.rim_a] == pytest.approx(np.sin(loop.a), rel=1e-14)
    assert rho[loop.rim_b] == pytest.approx(np.sin(loop.a), rel=1e-14)
    assert loop.n_band + loop.n_neck == loop.n
    assert loop.h_band == pytest.approx(loop.h_neck, rel=0.02)


@pytest.mark.parametrize("gluing", [ah.TORUS, ah.Gluing("torus", 1.0), ah.KLEIN])
def test_energy_form_is_the_stiffness_quadratic_form(loop, gluing):
    sec = ah.Sector(loop, 2, gluing)
    v = np.random.default_rng(1).standard_normal((sec.K.shape[0], 3))
    np.testing.assert_allclose(sec.energy_form(v),
                               np.einsum("ik,ij,jk->k", v, sec.K, v), rtol=1e-12)
    np.testing.assert_allclose(sec.K, sec.K.T, atol=1e-12)


def test_gluing_maps():
    np.testing.assert_allclose(ah.KLEIN.matrix(3), np.diag([1.0, -1.0]))
    T = ah.Gluing("torus", 0.7).matrix(2)
    np.testing.assert_allclose(T @ T.T, np.eye(2), atol=1e-15)
    np.testing.assert_allclose(T @ ah.J2, ah.J2 @ T, atol=1e-15)      # keeps U(1)
    K = ah.KLEIN.matrix(1)
    np.testing.assert_allclose(K @ ah.J2, -ah.J2 @ K, atol=1e-15)     # reverses it
    for g in (ah.TORUS, ah.KLEIN):
        np.testing.assert_allclose(g.matrix(0), [[1.0]])


@pytest.mark.parametrize("gluing,factor", [(ah.TORUS, 0.0), (ah.Gluing("torus", 1.0), 0.0),
                                           (ah.KLEIN, -2.0)])
def test_semi_discrete_angular_momentum_ledger(loop, gluing, factor):
    """dL/dt = −π m Σ (K x)_i · J x_i equals −2 Φ_B for Klein and 0 for torus."""
    m = 2
    sec = ah.Sector(loop, m, gluing)
    x = np.random.default_rng(2).standard_normal(sec.K.shape[0])
    Kx = (sec.K @ x).reshape(loop.n, 2)
    X = x.reshape(loop.n, 2)
    rate = -np.pi * m * np.sum(Kx * (X @ ah.J2.T))
    i, j = loop.rim_b, loop.rim_b - 1                 # out of the neck at rim B
    phi_b = ah.lflux(sec, i, j, X[i], X[j])
    assert abs(phi_b) > 1e-3
    assert rate == pytest.approx(factor * phi_b, abs=1e-9 * abs(phi_b))
    # and the two regions balance through their edges
    band, neck = ah.band_nodes(loop), ah.neck_nodes(loop)
    rate_neck = -np.pi * m * np.sum(Kx[neck] * (X[neck] @ ah.J2.T))
    phi_a_in = -ah.lflux(sec, loop.n - 1, 0, X[loop.n - 1], X[0])
    assert rate_neck == pytest.approx(phi_a_in - phi_b, rel=1e-10)


def test_m0_cannot_tell_the_gluings_apart(loop):
    K = [ah.Sector(loop, 0, g).K for g in (ah.TORUS, ah.Gluing("torus", 1.0), ah.KLEIN)]
    assert all(np.array_equal(K[0], k) for k in K[1:])


def test_torus_conserves_and_klein_is_the_union(loop):
    m = 2
    t0, t1 = ah.Sector(loop, m, ah.TORUS), ah.Sector(loop, m, ah.Gluing("torus", 1.0))
    tp, kl = ah.Sector(loop, m, ah.Gluing("torus", np.pi / m)), ah.Sector(loop, m, ah.KLEIN)
    r = ah.p2_torus(t1)
    assert r["degeneracy"] < 1e-12 and r["L_drift"] < 1e-10
    u = ah.p3_union(kl, t0, tp)
    assert max(u["union"], u["cos_vs_c0"], u["sin_vs_cpi"]) < 1e-12
    assert u["counts"] == [loop.n, loop.n]


def test_klein_doublets_reverse_at_the_predicted_time(loop):
    rows = ah.p4_p5(ah.Sector(loop, 1, ah.KLEIN), pairs=3)
    for r in rows:
        assert r["split"] > 1e-3
        assert r["timing_relative"] <= r["timing_bound"]
        assert r["identity"] < 1e-9
        assert r["port_fraction"] > 0.9
        assert np.sign(r["L_end"]) == -np.sign(r["L0"])


def test_trig_sum_folds_and_merges():
    f = ah.TrigSum(np.array([2.0, -2.0, 0.5]), np.array([1 + 1j, 1 - 1j, 0.3]))
    assert len(f.Omega) == 2
    t = np.linspace(0, 5, 11)
    direct = np.real((1 + 1j) * np.exp(2j * t) + (1 - 1j) * np.exp(-2j * t)
                     + 0.3 * np.exp(0.5j * t))
    np.testing.assert_allclose(f(t), direct, atol=1e-14)
    # cancelling pair merges to zero amplitude
    g = ah.TrigSum(np.array([3.0, -3.0]), np.array([1j, 1j]))
    assert len(g.Omega) == 0


def test_certified_first_sign_change_on_a_beat():
    w1, w2 = 5.0, 5.0 + 1e-6
    # ω̄ cos Δt − (Δ/2) cos 2ω̄t, written as two exponentials each
    mean, delta = 0.5 * (w1 + w2), w2 - w1
    f = ah.TrigSum(np.array([delta, 2 * mean]), np.array([mean, -delta / 2]))
    t = ah.first_sign_change(f, 2 * np.pi / delta, cut=mean, fine=0.5 / mean)
    assert t == pytest.approx(np.pi / (2 * delta), rel=delta / mean)


def test_sign_search_finds_a_shallow_crossing_between_samples():
    """Review of #355: f = 1 − 1e-7 − cos(t − 0.005) dips below zero only on a
    9e-4-wide window starting at t = 0.0045527864; sampling alone missed it."""
    f = ah.TrigSum(np.array([0.0, 1.0]), np.array([1 - 1e-7, -np.exp(-0.005j)]))
    t = ah.first_sign_change(f, 1.0, cut=2.0, fine=0.1)
    assert t == pytest.approx(0.005 - np.arccos(1 - 1e-7), abs=1e-12)
    # a function that keeps its sign is certified to, not merely sampled
    g = ah.TrigSum(np.array([0.0, 1.0]), np.array([1 + 1e-9, -1.0]))
    assert ah.first_sign_change(g, 20.0, cut=2.0, fine=0.1) is None


def _probe():
    spec = importlib.util.spec_from_file_location("ahp_probe", PROBE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _size(fail=(), worst=None):
    keys = ("P1", "P2 degeneracy", "P2 conservation", "P3", "P4 split", "P4 timing",
            "P5 identity", "P5 non-vacuous", "P6 sign", "P6 tolerance")
    w = dict.fromkeys(keys, 1.0)
    w.update(worst or {})
    return {"gates": {k: k not in fail for k in keys}, "worst": w}


GOOD_P7 = {"geometry": True, "spectrum": True}


def test_verdict_needs_both_registered_resolutions():
    pr = _probe()
    one = {2000: _size(fail=("P4 split",))}
    assert pr.verdict(one, GOOD_P7)[0] == "NO_VERDICT"
    both = {2000: _size(fail=("P4 split",)), 4000: _size(fail=("P4 split",))}
    assert pr.verdict(both, GOOD_P7)[:2] == ("COUNTEREXAMPLE", ["P4 split"])
    once = {2000: _size(fail=("P4 split",), worst={"P4 split": 1e-7}),
            4000: _size(worst={"P4 split": 2e-6})}
    assert pr.verdict(once, GOOD_P7)[0] == "UNRESOLVED_NUMERICALLY"
    assert pr.verdict(both, {"geometry": None, "spectrum": None})[0] == "NO_VERDICT"
    assert pr.verdict(both, {"geometry": False, "spectrum": True})[0] == "INCOMPLETE"


def test_unresolved_needs_the_miss_to_shrink():
    pr = _probe()
    worse = {2000: _size(fail=("P5 identity",), worst={"P5 identity": 4e-8}),
             4000: _size(fail=("P5 identity",), worst={"P5 identity": 4e-5})}
    v, why, notes = pr.verdict(worse, GOOD_P7)
    assert v == "NO_VERDICT" and why == ["P5 identity"] and notes
    better = {2000: _size(fail=("P5 identity",), worst={"P5 identity": 4e-8}),
              4000: _size(fail=("P5 identity",), worst={"P5 identity": 2e-8})}
    assert pr.verdict(better, GOOD_P7)[0] == "UNRESOLVED_NUMERICALLY"
    # with the one allowed doubling the comparison is 4000 → 8000
    redo = {**worse, 8000: _size(worst={"P5 identity": 1e-9})}
    assert pr.verdict(redo, GOOD_P7)[0] == "UNRESOLVED_NUMERICALLY"
    assert pr.verdict({**redo, 16000: _size()}, GOOD_P7)[0] == "NO_VERDICT"
    # a structural failure at both is reported with the other misses noted
    both = {n: _size(fail=("P4 split", "P5 identity"),
                     worst={"P5 identity": w}) for n, w in ((2000, 4e-8), (4000, 2e-7))}
    v, why, notes = pr.verdict(both, GOOD_P7)
    assert v == "COUNTEREXAMPLE" and why == ["P4 split"]
    assert notes == ["also missed: P5 identity (not shrinking)"]


def test_packet_through_the_handle_arrives_reversed():
    r = ah.p6_packet(0.5, 400, 1, samples=81)
    assert r["antisymmetry"] < 1e-3 < r["symmetry"]
    assert r["peak_fraction"] > 1e-2
    assert r["identity"] < 1e-8


def test_mesh_is_the_modelled_surface():
    for klein in (False, True):
        surf, vals, zero = ah.mesh_spectrum(klein, 2000, 24, count=4)
        g = ah.mesh_gluing(surf)
        if klein:
            assert g["sum_spread"] < 1e-12
        else:
            assert g["difference_spread"] < 1e-12
            assert g["difference"] == pytest.approx(np.pi, abs=1e-12)
        assert abs(zero) < 1e-10


def test_archived_probe():
    rep = json.loads(ARCHIVE.read_text())
    assert rep["prereg"] == "36ba9c2"
    assert rep["verdict"] in ("CONFIRMED", "COUNTEREXAMPLE",
                              "UNRESOLVED_NUMERICALLY", "INCOMPLETE")
    assert rep["sizes"] == [2000, 4000]
    assert rep["P7"]["geometry"] is True
    for n in ("2000", "4000"):
        g = rep["per_size"][n]["gates"]
        for key in ("P1", "P2 degeneracy", "P2 conservation", "P3", "P4 timing",
                    "P5 non-vacuous", "P6 sign", "P6 tolerance"):
            assert g[key], (n, key)
