"""
Tests for two closed histories sewn at one interaction.

The load-bearing claims are ranks, so the tests are about ranks — and about the
scope those ranks are stated in:

* five equations in five unknowns, so on a **fixed branch** the allowed events
  are discrete and locally isolated. Not that all roots were found, and not that
  the event is unique;
* the **branch scope is load-bearing**: `d` is the principal geodesic distance,
  and inside the principal delay band it is the only feasible branch, so every
  other test here is principal-branch by construction. Off it a mixed branch
  fixes the *difference* of distances — a hyperboloid, not an ellipsoid — and
  discreteness still survives per branch;
* removing a wave drops the rank to 4, but so does removing a **closure**, so
  that test is a dimensionality control and not a statement about photons;
* a shared throat fails both ways in this model — infeasible on **every** branch
  when traversed oppositely, rank-deficient on the same branch — with the
  same-traversal half scanned over branches rather than argued.

And the non-circularity check is tested as hard as the result: with the throat
delays free rather than given, *every* event on both fronts closes, so the whole
thing depends on the throat being data.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.pair_history import (
    PRINCIPAL,
    PairHistorySystem,
    Throat,
    all_branches,
    leg_length,
    winding_bound,
    measure_the_branch_set_is_finite_and_bounded_by_the_delay,
    measure_discrete_events_survive_branch_completion,
    measure_the_shared_throat_obstruction_survives_branch_completion,
    measure_the_delay_dependence_survives_branch_completion,
    measure_the_results_are_scoped_to_the_principal_branch,
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


def test_measure_the_events_are_discrete_and_locally_isolated():
    r = measure_the_event_is_selected_not_inserted(n_configs=4, n_starts=150)
    assert r["events_are_discrete_and_locally_isolated"]
    assert r["every_event_is_nondegenerate"]
    assert r["equations"] == r["unknowns"] == 5
    # and the claim is scoped, not universal
    assert r["branch_scope"].startswith("principal")
    assert "does not show all roots" in r["but_not_proved_exhaustive"]


def test_removing_a_wave_is_a_dimensionality_control_not_physics():
    """Deleting ANY one equation drops the rank — including a closure."""
    r = measure_removing_a_wave_removes_the_selection(n_configs=5,
                                                      n_starts=180)
    assert r["the_square_system_behaves_nondegenerately"]
    assert r["two_waves_give_isolated_events"]
    assert r["one_wave_gives_a_one_parameter_family"]
    assert r["deleting_a_closure_instead_drops_the_rank_the_same_way"]
    assert r["nullity_with_one_wave"] == 1
    assert "not a statement about waves" in \
        r["this_is_a_dimensionality_control_not_a_physics_result"] or \
        "nothing here singles out the wave" in \
        r["this_is_a_dimensionality_control_not_a_physics_result"]


def test_measure_a_shared_throat_cannot_carry_the_pair():
    r = measure_a_shared_throat_cannot_carry_the_pair(n_configs=3,
                                                      n_starts=140)
    assert r["opposite_traversal_is_infeasible"]
    assert r["same_traversal_loses_a_rank"]
    assert r["so_in_this_model_the_pair_needs_two_distinct_throats"]
    # the same-traversal half is scanned over branches, not argued
    assert r["no_branch_pair_rescues_a_shared_throat"]
    assert "not excluded" in r["the_same_traversal_result_is_scoped"]


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


# ── the branch scope, which every other claim depends on ────────────────────
def test_leg_length_covers_short_long_and_winding():
    d = 1.0
    assert leg_length(d, 0, 0) == pytest.approx(d)
    assert leg_length(d, 1, 0) == pytest.approx(2 * math.pi - d)
    assert leg_length(d, 0, 1) == pytest.approx(d + 2 * math.pi)
    assert leg_length(d, 1, 1) == pytest.approx(4 * math.pi - d)


def test_all_branches_enumerates_the_labels():
    assert PRINCIPAL in all_branches(0)
    assert len(all_branches(0)) == 4
    assert len(all_branches(1)) == 16


def test_inside_the_principal_band_only_one_branch_is_feasible():
    rng = np.random.default_rng(2)
    th = _throat(rng)                    # drawn inside [D, 2pi - D]
    assert th.feasible_branches(1) == [PRINCIPAL]


def test_off_the_principal_branch_the_locus_changes_kind():
    """A mixed branch fixes the DIFFERENCE of distances, not the sum."""
    rng = np.random.default_rng(2)
    mp, mm = _nrm(rng.normal(size=4)), _nrm(rng.normal(size=4))
    d = geodesic_distance(mp, mm)
    th = Throat(mp, mm, -(2 * math.pi))          # the difference window
    mixed = [b for b in th.feasible_branches(1) if b[0] != b[2]]
    assert mixed, "no mixed branch feasible at |delta| = 2pi"
    # a mixed branch really is a difference condition
    b = mixed[0]
    c = _nrm(rng.normal(size=4))
    d1, d2 = geodesic_distance(c, mp), geodesic_distance(mm, c)
    got = th.closure_residual(c, b)
    expect = leg_length(d1, b[0], b[1]) + leg_length(d2, b[2], b[3]) + th.delay
    assert got == pytest.approx(expect)


def test_measure_the_results_are_scoped_to_the_principal_branch():
    r = measure_the_results_are_scoped_to_the_principal_branch(n_configs=3,
                                                               n_starts=140)
    assert r["inside_the_band_only_the_principal_branch_is_feasible"]
    assert r["outside_it_more_branches_open"]
    assert r["off_branch_loci_are_difference_type"]
    assert "difference" in r["locus_kinds_off_the_principal_branch"]
    assert r["discreteness_survives_on_a_fixed_off_branch"]


# ── branch completeness, which is what finishes the round ───────────────────
def test_the_winding_is_bounded_by_the_delay():
    """l = (d or 2pi-d) + 2pi k >= 2pi k and l1+l2 = |delta|."""
    assert winding_bound(-0.5) == 0
    assert winding_bound(-2 * math.pi) == 1
    assert winding_bound(-6.9 * math.pi) == 3
    assert winding_bound(0.0) == 0


def test_no_feasible_branch_exceeds_the_bound():
    rng = np.random.default_rng(12)
    for _ in range(40):
        mp, mm = _nrm(rng.normal(size=4)), _nrm(rng.normal(size=4))
        delay = -float(rng.uniform(0.05, 5 * 2 * math.pi))
        th = Throat(mp, mm, delay, "t")
        bound = winding_bound(delay)
        for b in all_branches(9):
            if th.branch_is_feasible(b):
                assert b[1] + b[3] <= bound


def test_feasible_branches_is_complete_by_default():
    """No cutoff to choose: the default IS the whole feasible set."""
    rng = np.random.default_rng(13)
    mp, mm = _nrm(rng.normal(size=4)), _nrm(rng.normal(size=4))
    th = Throat(mp, mm, -3.0 * 2 * math.pi, "t")
    complete = th.feasible_branches()
    generous = [b for b in all_branches(9) if th.branch_is_feasible(b)]
    assert sorted(complete) == sorted(generous)


def test_branch_complete_solving_is_a_superset():
    rng = np.random.default_rng(14)
    for _ in range(6):
        sysm = _wide(rng)
        full = sysm.solve_branch_complete(n_starts=120, seed=1)
        if not full:
            continue
        assert all(e["rank"] == 5 for e in full)
        assert all("branch_a" in e and "branch_b" in e for e in full)
        return
    pytest.fail("no wide configuration admitted a branch-complete event")


def _wide(rng):
    from geometrodynamics.waves.pair_history import _wide_system
    return _wide_system(rng)


def test_measure_the_branch_set_is_finite_and_bounded_by_the_delay():
    r = measure_the_branch_set_is_finite_and_bounded_by_the_delay(n_random=60,
                                                                  brute_winding=8)
    assert r["no_feasible_branch_exceeds_the_bound"]
    assert r["the_branch_set_is_finite"]
    assert r["worst_excess_over_the_bound"] <= 0


def test_measure_discrete_events_survive_branch_completion():
    """The roadmap's gate: isolated events must persist, or stop here."""
    r = measure_discrete_events_survive_branch_completion(n_configs=5,
                                                          n_starts=120)
    assert r["so_discreteness_survives_branch_completion"]
    assert r["every_event_still_isolated"]
    assert r["events_at_full_rank_5"] == r["total_events"] > 0


def test_measure_the_shared_throat_obstruction_survives_branch_completion():
    r = measure_the_shared_throat_obstruction_survives_branch_completion(
        n_configs=5, n_starts=130)
    assert r["the_obstruction_survives_branch_completion"]
    assert r["pairs_that_restored_full_rank"] == 0
    assert r["distinct_branch_pairs_tested"] > 0


def test_measure_the_delay_dependence_survives_branch_completion():
    r = measure_the_delay_dependence_survives_branch_completion(n_points=200)
    assert r["so_the_delay_must_still_be_given"]
    assert r["fraction_closable_by_choosing_a_delay"] == 1.0
