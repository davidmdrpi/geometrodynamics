"""
Tests for the S³ wormhole loop and its ledger.

The load-bearing claims are three, and none of them is the ledger closing —
that is an assumption, and it is tested as one:

* the antipodal staging is **geometry**: a geodesic sphere has area
  ``4π sin²χ``, so a shell refocuses exactly at ``χ = π``, and the same fact
  placed both the future mouth and the receiver;
* a **linear** wave on a time-displaced loop has exactly one self-consistent
  amplitude, for every ``κ ≠ 1`` — which is a fact about linear equations, so
  the quadratic control that breaks it is tested too;
* the drawn shell **is** the geodesic level set, and its drawn size is
  ``√(A/4π)`` — the picture measures rather than illustrates.

Also held down, because both were got wrong first:

* a receiver placed at a generic point is *arrived at*, not collapsed onto —
  ``dA/dχ`` has the wrong sign there — so "collapse" has to be a measured sign;
* a stereographic chart is unbounded at its own pole, and a shell sweeps all of
  ``S³``, so no pole is safe; the orthographic shadow is bounded everywhere.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.transaction.s3_geometry import antipode4, geo4, nrm4
from geometrodynamics.viz.wormhole_ledger import (
    ANTIPODE,
    Mouth,
    WormholeLoop,
    enclosed_volume,
    geo_point_at,
    geodesic_shell,
    measure_nonlinearity_is_what_would_break_it,
    measure_only_the_resonance_obstructs,
    measure_the_delay_does_not_enter_the_ledger,
    measure_the_drawn_shell_is_the_geodesic_level_set,
    measure_the_ledger_closes,
    measure_the_loop_has_one_self_consistent_solution,
    measure_the_receiver_is_struck_by_a_collapsing_shell,
    measure_the_shell_focuses_at_the_antipode,
    measure_two_local_events_one_conserved_wave,
    shadow,
    shell_amplitude,
    shell_area,
)


# ── the area law, which is what makes the staging geometry ──────────────────
def test_shell_area_vanishes_at_both_poles():
    assert shell_area(0.0) == 0.0
    assert abs(shell_area(math.pi)) < 1e-28


def test_shell_area_is_maximal_at_the_equator():
    assert shell_area(math.pi / 2) == pytest.approx(4.0 * math.pi)
    for chi in (0.3, 1.0, 2.0, 2.9):
        assert shell_area(chi) < shell_area(math.pi / 2)


def test_enclosed_volume_is_the_integral_of_the_area():
    h = 1e-6
    for chi in (0.4, 1.1, 2.0, 2.8):
        d = (enclosed_volume(chi + h) - enclosed_volume(chi - h)) / (2 * h)
        assert d == pytest.approx(shell_area(chi), rel=1e-8)


def test_total_volume_of_the_three_sphere():
    assert enclosed_volume(math.pi) == pytest.approx(2.0 * math.pi ** 2)


def test_amplitude_diverges_at_the_refocus():
    assert shell_amplitude(math.pi / 2) == pytest.approx(1.0)
    assert shell_amplitude(math.pi - 1e-8) > 1e7


# ── the loop's geometry ─────────────────────────────────────────────────────
def test_the_future_mouth_is_exactly_the_antipode():
    assert WormholeLoop().emitter_to_future_mouth == pytest.approx(ANTIPODE,
                                                                   abs=1e-12)


def test_the_default_receiver_is_the_past_mouths_antipode():
    """Not decoration: it is why the arriving shell is collapsing."""
    loop = WormholeLoop()
    assert loop.past_mouth_to_receiver == pytest.approx(ANTIPODE, abs=1e-12)


def test_a_displaced_receiver_is_still_accepted():
    loop = WormholeLoop()
    other = WormholeLoop(receiver=Mouth(
        geo_point_at(loop.past_mouth.position, loop.emitter.position, 1.2),
        0.0, "receiver"))
    assert other.past_mouth_to_receiver == pytest.approx(1.2, abs=1e-12)


def test_geo_point_at_lands_at_the_requested_distance():
    a = nrm4([0.2, 0.5, -0.3, 0.7])
    b = nrm4([-0.4, 0.1, 0.8, 0.2])
    for chi in (0.1, 1.0, 2.5, 3.0):
        assert geo4(a, geo_point_at(a, b, chi)) == pytest.approx(chi,
                                                                 abs=1e-12)


def test_geo_point_at_refuses_antipodal_endpoints():
    a = nrm4([0.0, 0.0, 0.0, 1.0])
    with pytest.raises(ValueError):
        geo_point_at(a, antipode4(a), 1.0)


def test_the_wave_path_only_ever_moves_forward():
    """Path length is monotone; it is the TIME the throat moves backwards."""
    knots = WormholeLoop(delay=-12.0).wave_path()
    s = [k["s"] for k in knots]
    assert s == sorted(s)
    assert s[-1] == pytest.approx(2.0 * ANTIPODE, abs=1e-12)
    assert knots[2]["t"] < knots[1]["t"]      # the throat jumps back in time
    assert knots[2]["s"] == pytest.approx(knots[1]["s"])


def test_arrival_times_report_the_ordering_rather_than_assume_it():
    assert WormholeLoop(delay=-12.0).arrival_times()[
        "receiver_before_emitter"] is True
    assert WormholeLoop(delay=-1.0).arrival_times()[
        "receiver_before_emitter"] is False


# ── the self-consistent amplitude ───────────────────────────────────────────
@pytest.mark.parametrize("kappa", [0.2, -0.5, 0.9, 0.3 + 0.4j])
def test_the_fixed_point_solves_the_loop(kappa):
    loop = WormholeLoop(kappa=kappa)
    a = loop.self_consistent_amplitude()
    assert abs(a - (loop.source + kappa * a)) < 1e-12


def test_the_fixed_point_matches_brute_iteration():
    loop = WormholeLoop(kappa=0.9)
    assert abs(loop.iterate_loop(4000)
               - loop.self_consistent_amplitude()) < 1e-12


def test_the_resonance_is_refused():
    with pytest.raises(ValueError):
        WormholeLoop(kappa=1.0).self_consistent_amplitude()


def test_the_fixed_point_exists_where_the_iteration_diverges():
    """|κ| > 1: the sum blows up and the solution is still there and unique."""
    loop = WormholeLoop(kappa=2.5)
    a = loop.self_consistent_amplitude()
    assert math.isfinite(abs(a))
    assert abs(a - (loop.source + 2.5 * a)) < 1e-12
    assert abs(loop.iterate_loop(200)) > 1e30


# ── the projection, which is what the picture rests on ──────────────────────
def test_the_drawn_shell_is_the_geodesic_level_set():
    c = nrm4([0.3, -0.5, 0.7, 0.4])
    for chi in (0.15, 1.0, math.pi / 2, 2.6):
        for p in geodesic_shell(c, chi, 8, 16):
            assert geo4(c, p) == pytest.approx(chi, abs=1e-12)
            assert np.linalg.norm(p) == pytest.approx(1.0, abs=1e-12)


def test_the_shadow_never_leaves_the_unit_ball():
    rng = np.random.default_rng(5)
    for _ in range(200):
        q = nrm4(rng.normal(size=4))
        p3, depth = shadow(q)
        assert np.linalg.norm(p3) <= 1.0 + 1e-12
        assert -1.0 - 1e-12 <= depth <= 1.0 + 1e-12
        assert np.linalg.norm(p3) ** 2 + depth ** 2 == pytest.approx(1.0)


def test_the_drawn_size_is_the_square_root_of_the_area():
    """One constant across χ — so the apparent size IS √(A/4π)."""
    c = nrm4([0.3, -0.5, 0.7, 0.4])
    ratios = []
    for chi in (0.2, 0.8, math.pi / 2, 2.4, 3.0):
        s = np.asarray([shadow(p)[0] for p in geodesic_shell(c, chi, 12, 26)])
        ratios.append(float((s.max(axis=0) - s.min(axis=0)).max())
                      / math.sin(chi))
    assert max(ratios) - min(ratios) < 1e-12


def test_no_stereographic_pole_is_safe_for_this_scene():
    """Whatever pole is chosen, the shell crosses it and the chart diverges."""
    e = WormholeLoop().emitter.position
    pole = nrm4([0.2, 0.7, -0.4, 0.5])
    chi_star = geo4(e, pole)
    assert 0.0 < chi_star < math.pi          # the shell does reach it
    prev = None
    for eps in (1e-2, 1e-3, 1e-4):
        p = geo_point_at(e, pole, chi_star + eps)
        r = float(np.linalg.norm(p - np.dot(p, pole) * pole)) / (
            1.0 - float(np.dot(p, pole)))
        assert r * eps == pytest.approx(2.0, abs=1e-3)   # diverges as 2/ε
        if prev is not None:
            assert r > 8.0 * prev                        # and does not settle
        prev = r
        assert np.linalg.norm(shadow(p)[0]) <= 1.0 + 1e-12


# ── the measurements ────────────────────────────────────────────────────────
def test_measure_the_shell_focuses_at_the_antipode():
    r = measure_the_shell_focuses_at_the_antipode(samples=60_000)
    assert r["the_area_law_holds"]
    assert r["derivative_of_volume_is_the_area"]
    assert r["total_volume_is_two_pi_squared"]
    assert r["area_vanishes_at_both_poles"]
    assert r["tested_against_the_sampling_error_not_a_fixed_percent"]


def test_measure_the_receiver_is_struck_by_a_collapsing_shell():
    r = measure_the_receiver_is_struck_by_a_collapsing_shell()
    assert r["default_arrival_is_the_antipode"]
    assert r["the_shell_collapses_only_at_the_refocus"]
    assert r["refocused_receiver"]["collapsing_all_the_way_in"]
    assert r["displaced_receiver"]["expanding_all_the_way_in"]
    assert r["the_focused_arrival_has_vanishing_area"]
    assert r["the_displaced_arrival_does_not"]


def test_measure_the_drawn_shell_is_the_geodesic_level_set():
    r = measure_the_drawn_shell_is_the_geodesic_level_set(n_random=20)
    assert r["the_drawn_points_are_the_level_set"]
    assert r["the_shadow_never_leaves_the_unit_ball"]
    assert r["the_drawn_size_is_sqrt_of_the_area"]
    assert r["the_stereographic_radius_diverges_as_two_over_epsilon"]
    assert r["it_does_not_converge_under_refinement"]
    assert r["while_the_shadow_stays_in_the_unit_ball"]


def test_measure_the_loop_has_one_self_consistent_solution():
    # the default 4000 rounds is not slack: the iteration's own truncation is
    # |κ|^rounds × the amplification, so at κ = 0.99 (amplification ×100) a
    # 2000-round sum is still 2e-7 from the fixed point by construction
    r = measure_the_loop_has_one_self_consistent_solution()
    assert r["the_fixed_point_is_what_the_loop_does"]
    assert r["worst_abs_error"] < 1e-9


def test_the_iteration_truncates_geometrically():
    """Why the round count matters — the residual is |κ|^n, not noise."""
    loop = WormholeLoop(kappa=0.9)
    exact = loop.self_consistent_amplitude()
    for n in (50, 100, 200):
        assert abs(loop.iterate_loop(n) - exact) == pytest.approx(
            abs(exact) * 0.9 ** n, rel=1e-9)


def test_measure_only_the_resonance_obstructs():
    r = measure_only_the_resonance_obstructs()
    assert r["failures"] == 0
    assert r["the_fixed_point_exists_everywhere_but_one_point"]
    assert r["the_resonance_is_refused"]


def test_measure_nonlinearity_is_what_would_break_it():
    """The control that stops uniqueness reading as a claim about wormholes."""
    r = measure_nonlinearity_is_what_would_break_it()
    assert r["two_solutions_below_threshold"]
    assert r["none_above_it"]
    assert r["so_uniqueness_is_a_linearity_result"]
    assert r["every_reported_root_solves_the_loop"]
    assert r["threshold_source"] == pytest.approx(1.0 / (4.0 * r["kappa"]))


def test_measure_the_ledger_closes():
    r = measure_the_ledger_closes()
    assert r["every_ledger_closes"]
    assert r["worst_residual"] < 1e-12
    # and it is labelled as the assumption it is
    assert "put in" in r["but_that_is_the_assumption_not_the_result"]


def test_measure_the_delay_does_not_enter_the_ledger():
    r = measure_the_delay_does_not_enter_the_ledger()
    assert r["the_delay_changes_no_conserved_quantity"]
    assert r["but_it_decides_the_ordering"]
    assert r["amplitude_spread"] < 1e-12


def test_measure_two_local_events_one_conserved_wave():
    r = measure_two_local_events_one_conserved_wave()
    assert r["it_is_exactly_the_antipode"]
    assert r["the_pair_conserves"]
    assert not r["emitter_alone_conserves"]
    assert not r["receiver_alone_conserves"]
    assert r["the_pair_closes_to"] < 1e-12
