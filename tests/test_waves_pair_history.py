"""
Tests for two closed histories sewn at one interaction.

The load-bearing claims are ranks, so the tests are about ranks:

* five equations in five unknowns — the interaction event is **isolated**, with
  the Jacobian at full rank, so isolation is a property of the system and not of
  the solver;
* remove one incoming wave and the rank drops to 4: the solutions become a
  **one-parameter family**. They do not vanish, which is the point;
* a single shared throat cannot carry the conjugate pair — infeasible when
  traversed oppositely, rank-deficient when traversed the same way.

And the non-circularity check is tested as hard as the result: with the throat
delays free rather than given, *every* event on both fronts closes, so the whole
thing depends on the throat being data.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.pair_history import (
    PairHistorySystem,
    Throat,
    feasible_delay_band,
    geodesic_distance,
    measure_a_shared_throat_cannot_carry_the_pair,
    measure_closure_is_a_geodesic_ellipsoid,
    measure_removing_a_wave_removes_the_selection,
    measure_the_conjugacy_is_carried_not_derived,
    measure_the_delays_must_be_given_not_solved_for,
    measure_the_event_is_selected_not_inserted,
    measure_the_threshold_is_a_separate_condition,
)


def _nrm(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def _throat(rng, label="T"):
    mp, mm = _nrm(rng.normal(size=4)), _nrm(rng.normal(size=4))
    lo, hi = feasible_delay_band(geodesic_distance(mp, mm))
    return Throat(mp, mm, -float(rng.uniform(lo + 0.05, hi - 0.05)), label)


def _system(seed=0, shared=None):
    rng = np.random.default_rng(seed)
    sa, sb = _nrm(rng.normal(size=4)), _nrm(rng.normal(size=4))
    ta = _throat(rng, "A")
    if shared == "same":
        tb = Throat(ta.m_plus, ta.m_minus, ta.delay, "B")
    elif shared == "opposite":
        tb = ta.reversed()
    else:
        tb = _throat(rng, "B")
    return PairHistorySystem(sa, sb, ta, tb, 0.0, 0.3)


# ── the closure locus ───────────────────────────────────────────────────────
def test_feasible_band_is_the_mouth_separation_and_its_complement():
    lo, hi = feasible_delay_band(1.2)
    assert lo == pytest.approx(1.2)
    assert hi == pytest.approx(2 * math.pi - 1.2)


def test_a_delay_outside_the_band_is_infeasible():
    a, b = _nrm([0, 0, 0, 1.0]), _nrm([1.0, 0, 0, 0])
    d = geodesic_distance(a, b)
    assert not Throat(a, b, -(0.5 * d)).is_feasible()
    assert not Throat(a, b, -(2 * math.pi - 0.5 * d)).is_feasible()
    assert Throat(a, b, -(d + 0.4)).is_feasible()


def test_closure_residual_vanishes_exactly_on_the_ellipsoid():
    rng = np.random.default_rng(4)
    th = _throat(rng)
    # a point on the geodesic between the mouths has summed distance = d
    tight = Throat(th.m_plus, th.m_minus, -th.mouth_separation)
    assert tight.closure_residual(th.m_plus) == pytest.approx(0.0, abs=1e-12)
    assert tight.closure_residual(th.m_minus) == pytest.approx(0.0, abs=1e-12)


def test_reversed_throat_flips_the_delay_and_the_mouths():
    rng = np.random.default_rng(6)
    th = _throat(rng)
    rv = th.reversed()
    assert rv.delay == pytest.approx(-th.delay)
    assert np.allclose(rv.m_plus, th.m_minus)
    assert rv.mouth_separation == pytest.approx(th.mouth_separation)


# ── the system and its rank ─────────────────────────────────────────────────
def test_the_system_is_five_equations_in_five_unknowns():
    sysm = _system(1)
    u = np.concatenate([_nrm([0.2, -0.5, 0.7, 0.4]), [1.3]])
    assert len(sysm.residuals(u)) == 5
    assert sysm.jacobian(u).shape == (5, 5)


def test_selected_events_are_isolated_at_full_rank():
    for seed in (1, 2, 3, 5):
        sysm = _system(seed)
        sols = sysm.solve(n_starts=220, seed=seed)
        for s in sols:
            assert s["rank"] == 5
            assert s["worst_residual"] < 1e-9
        if sols:
            return                      # at least one configuration must solve
    pytest.fail("no configuration admitted a pair-history")


def test_removing_a_wave_drops_the_rank_to_four():
    for seed in (1, 2, 3, 5, 8):
        sysm = _system(seed)
        one = sysm.solve(n_starts=220, seed=seed, with_b_wave=False,
                         dedupe=1e-4)
        if len(one) > 5:
            assert {s["rank"] for s in one} == {4}
            return
    pytest.fail("no configuration produced a one-wave family")


def test_the_event_must_follow_both_launches():
    sysm = _system(2)
    for s in sysm.solve(n_starts=200, seed=2):
        assert s["t"] > max(sysm.tau_a, sysm.tau_b)


def test_an_opposite_shared_throat_is_infeasible():
    sysm = _system(3, shared="opposite")
    # B's delay is positive, so its closure would need a negative path length
    assert sysm.throat_b.delay > 0
    assert not sysm.throat_b.is_feasible()


def test_a_same_sense_shared_throat_duplicates_a_constraint():
    sysm = _system(3, shared="same")
    u = np.concatenate([_nrm([0.1, 0.4, -0.6, 0.5]), [1.1]])
    r = sysm.residuals(u)
    assert r[2] == pytest.approx(r[3])          # the two closures coincide
    assert sysm.rank_at(u) <= 4


def test_a_same_sign_pair_is_refused():
    rng = np.random.default_rng(9)
    with pytest.raises(ValueError):
        PairHistorySystem(_nrm(rng.normal(size=4)), _nrm(rng.normal(size=4)),
                          _throat(rng), _throat(rng), orientations=(+1, +1))


def test_net_orientation_cancels():
    assert _system(4).net_orientation() == 0


# ── the invariant is a separate question ────────────────────────────────────
def test_the_invariant_at_a_selected_event_is_well_defined():
    for seed in (1, 2, 3, 5):
        sysm = _system(seed)
        sols = sysm.solve(n_starts=220, seed=seed)
        if not sols:
            continue
        inv = sysm.invariant_at(sols[0], energy=1.0)
        assert 0.0 <= inv["opening_angle"] <= math.pi
        assert inv["s"] == pytest.approx(
            2.0 * (1.0 - math.cos(inv["opening_angle"])), rel=1e-12)
        assert inv["above_threshold"] == bool(inv["s"] >= 4.0)
        return
    pytest.fail("no configuration admitted a pair-history")


# ── the measurements ────────────────────────────────────────────────────────
def test_measure_closure_is_a_geodesic_ellipsoid():
    r = measure_closure_is_a_geodesic_ellipsoid(samples=6000)
    assert r["sampling_never_goes_below_the_band"]
    assert r["sampling_never_goes_above_the_band"]
    assert r["an_infeasible_delay_is_rejected"]


def test_measure_the_event_is_selected_not_inserted():
    r = measure_the_event_is_selected_not_inserted(n_configs=4, n_starts=150)
    assert r["the_event_is_selected_not_inserted"]
    assert r["every_event_is_nondegenerate"]
    assert r["equations"] == r["unknowns"] == 5


def test_measure_removing_a_wave_removes_the_selection():
    """The falsification: isolation is lost, existence is not."""
    r = measure_removing_a_wave_removes_the_selection(n_configs=4,
                                                      n_starts=150)
    assert r["the_selection_requires_both_waves"]
    assert r["two_waves_give_isolated_events"]
    assert r["one_wave_gives_a_one_parameter_family"]
    assert r["nullity_with_one_wave"] == 1
    assert r["the_solutions_do_not_vanish_they_stop_being_isolated"]


def test_measure_a_shared_throat_cannot_carry_the_pair():
    r = measure_a_shared_throat_cannot_carry_the_pair(n_configs=3,
                                                      n_starts=140)
    assert r["opposite_traversal_is_infeasible"]
    assert r["same_traversal_loses_a_rank"]
    assert r["so_the_pair_needs_two_distinct_throats"]


def test_measure_the_delays_must_be_given_not_solved_for():
    """If Δ were solved for, every event would close and nothing is selected."""
    r = measure_the_delays_must_be_given_not_solved_for(n_starts=120)
    assert r["with_free_delays_almost_any_event_closes"]
    assert r["nullity_with_delays_free"] == 2
    assert r["fraction_closable_by_choosing_a_delay"] > 0.9


def test_measure_the_threshold_is_a_separate_condition():
    r = measure_the_threshold_is_a_separate_condition(n_configs=6,
                                                      n_starts=140)
    assert r["closure_selects_where_the_invariant_decides_whether"]
    if r["selected_events"]:
        assert r["none_clear_threshold_at_energy_equal_mass"]
        by_e = {row["energy_over_mass"]: row["fraction_clearing_threshold"]
                for row in r["rows"]}
        assert by_e[1.0] == 0.0
        assert by_e[3.0] >= by_e[1.5] >= by_e[1.0]


def test_measure_the_conjugacy_is_carried_not_derived():
    r = measure_the_conjugacy_is_carried_not_derived()
    assert r["the_labels_cancel"]
    assert r["a_same_sign_pair_is_refused"]
    assert "bookkeeping" in r["but_nothing_here_derives_charge"]
