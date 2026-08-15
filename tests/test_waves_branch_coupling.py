"""
Tests for branch-resolved field coupling — the throat solved for, not applied.

The load-bearing claims:

* the branch (image) series sums in **closed form**, and its poles are the
  conformal ESU eigenfrequencies with residues equal to the mode functions over
  ``2ω`` — so the image and mode representations are one function;
* ``1/(1 − L)`` is the sum over **every** throat traversal, and PR #254's
  post-processed answer is its ``n = 0`` term;
* the solved field therefore has arrivals at times the free-branch ledger does
  not contain, on a ``κⁿ`` ladder, signed by every Maslov factor in the word;
* the primitive ``K_ab`` is indexed by a **pair** of branches: the amplitude
  factorizes over that index (rank one) and PR #253's closure condition does
  not;
* and the expansion PR #254 implicitly truncated has an expansion parameter with
  no bound as the regulator is removed — ``κ_c ∝ γ``.
"""

from __future__ import annotations

import cmath
import math

import numpy as np
import pytest

from geometrodynamics.waves.branch_coupling import (
    TWO_PI,
    CoupledThroat,
    branch_pair_matrix,
    branch_transfer,
    coupled_arrivals,
    coupled_propagator,
    coupled_waveform,
    critical_coupling,
    esu_mode_weight,
    free_branch_propagator,
    leg_branches,
    measure_closure_is_broadband_coherence,
    measure_solving_the_throat_resums_every_traversal,
    measure_the_closed_form_transfer_is_the_branch_sum,
    measure_the_coupled_field_has_arrivals_the_free_branches_do_not,
    measure_the_expansion_fails_at_the_eigenfrequencies,
    measure_the_primitive_is_rank_one_for_one_throat_and_not_for_two,
    mouth_transfer,
    traversal_series,
)


# ── one leg ─────────────────────────────────────────────────────────────────
def test_leg_branches_are_the_ray_ledger():
    chi = 1.1
    got = leg_branches(chi, n_k=3)
    assert len(got) == 6
    assert [r["ell"] for r in got] == sorted(r["ell"] for r in got)
    for k in range(3):
        assert any(abs(r["ell"] - (chi + TWO_PI * k)) < 1e-12 for r in got)
        assert any(abs(r["ell"] - (TWO_PI * (k + 1) - chi)) < 1e-12
                   for r in got)


def test_the_short_way_is_always_plus_and_the_long_way_always_minus():
    """The whole reason the winding sum is two geometric series."""
    for r in leg_branches(0.8, n_k=5):
        assert r["sign"] == (-1) ** r["crossings"]
        assert r["sign"] == (-1 if r["long_way"] else 1)


@pytest.mark.parametrize("bad", [0.0, math.pi, -0.3, 4.0])
def test_a_leg_needs_a_nondegenerate_geodesic(bad):
    with pytest.raises(ValueError):
        leg_branches(bad)


# ── the transfer function ───────────────────────────────────────────────────
@pytest.mark.parametrize("chi", [0.5, 1.3, 2.2])
@pytest.mark.parametrize("omega", [0.4, 1.7, 3.3])
def test_the_closed_form_is_the_branch_series(chi, omega):
    a = branch_transfer(omega, chi, damping=0.06, n_k=400)
    b = mouth_transfer(omega, chi, damping=0.06)
    assert abs(a - b) < 1e-12


def test_the_transfer_is_singular_at_the_poles_of_the_sphere():
    for bad in (0.0, math.pi):
        with pytest.raises(ValueError):
            mouth_transfer(1.0, bad)


def test_the_transfer_blows_up_at_the_esu_eigenfrequencies():
    """``1 − e^{−2πu}`` vanishes at ``ω ∈ ℤ`` when the regulator goes."""
    mags = [abs(mouth_transfer(2.0, 1.3, g)) for g in (1e-2, 1e-3, 1e-4)]
    assert mags[1] > 8.0 * mags[0]
    assert mags[2] > 8.0 * mags[1]


def test_the_transfer_is_finite_at_zero_frequency():
    """``ω = 0`` would be a pole, except the numerator vanishes there too."""
    assert abs(mouth_transfer(0.0, 1.3, 0.02)) < 10.0


@pytest.mark.parametrize("m", [1, 2, 3, 4])
def test_the_residues_are_the_esu_mode_functions_over_two_omega(m):
    chi, g = 1.3, 1e-6
    u = complex(g, float(m))
    res = (mouth_transfer(float(m), chi, g)
           * (1.0 - cmath.exp(-u * TWO_PI)) / (TWO_PI * 1j))
    want = -esu_mode_weight(m, chi) / (2.0 * m)
    assert res.real == pytest.approx(want, rel=1e-4)
    assert abs(res.imag) < 1e-4 * abs(want)


def test_conjugation_symmetry_makes_the_field_real():
    a = mouth_transfer(2.3, 1.1, 0.03)
    b = mouth_transfer(-2.3, 1.1, 0.03)
    assert b == pytest.approx(a.conjugate(), abs=1e-12)


# ── the throat, solved ──────────────────────────────────────────────────────
def test_the_resolvent_is_the_geometric_series():
    th = CoupledThroat(1.3, 1.0, +1, 0.3)
    for w in (0.4, 1.1, 2.7, 5.3):
        L = th.loop_transfer(w, 0.05)
        assert abs(L) < 1.0
        walked = sum(L ** n for n in range(400))
        assert abs(th.resolvent(w, 0.05) - walked) < 1e-14


def test_a_decoupled_throat_is_the_identity():
    th = CoupledThroat(1.3, 1.0, +1, 0.0)
    assert th.resolvent(2.0, 0.05) == pytest.approx(1.0 + 0j, abs=1e-15)
    assert abs(free_branch_propagator(2.0, 1.2, 0.9, th, 0.05)) == 0.0


def test_the_orientation_flips_the_loop_gain_sign():
    a = CoupledThroat(1.3, 1.0, +1, 0.3).loop_transfer(2.0, 0.05)
    b = CoupledThroat(1.3, 1.0, -1, 0.3).loop_transfer(2.0, 0.05)
    assert b == pytest.approx(-a, abs=1e-15)


def test_the_solved_propagator_is_the_traversal_sum():
    th = CoupledThroat(1.3, 1.0, +1, 0.45)
    for w in (0.4, 1.1, 2.7, 5.3):
        assert abs(coupled_propagator(w, 1.2, 0.9, th, 0.05)
                   - traversal_series(w, 1.2, 0.9, th, 0.05, 400)) < 1e-14


def test_post_processing_is_the_first_term_and_misses_by_the_loop_gain():
    """``|1/(1−L) − 1| / |1/(1−L)| = |L|`` exactly."""
    th = CoupledThroat(1.3, 1.0, +1, 0.4)
    for w in (0.9, 2.2, 4.4):
        s = coupled_propagator(w, 1.2, 0.9, th, 0.05)
        f = free_branch_propagator(w, 1.2, 0.9, th, 0.05)
        assert abs(s - f) / abs(s) == pytest.approx(
            abs(th.loop_transfer(w, 0.05)), rel=1e-12)


# ── the primitive ───────────────────────────────────────────────────────────
def test_the_branch_pair_matrix_sums_to_the_one_traversal_propagator():
    th = CoupledThroat(1.3, 1.0, +1, 0.3)
    got = branch_pair_matrix(2.1, 1.2, 0.9, th, damping=0.12, n_k=300)
    want = free_branch_propagator(2.1, 1.2, 0.9, th, damping=0.12)
    assert got["pair_sum"] == pytest.approx(want, abs=1e-12)
    assert complex(got["K"].sum()) == pytest.approx(want, abs=1e-12)


def test_the_primitive_is_indexed_by_a_pair_of_branches():
    th = CoupledThroat(1.3, 1.0, +1, 0.3)
    got = branch_pair_matrix(2.1, 1.2, 0.9, th, n_k=4)
    assert got["K"].shape == (len(got["rows"]), len(got["cols"]))
    assert got["K"].shape == (8, 8)


def test_one_throat_makes_the_primitive_rank_one():
    th = CoupledThroat(1.3, 1.0, +1, 0.3)
    K = branch_pair_matrix(2.1, 1.2, 0.9, th, n_k=6)["K"]
    s = np.linalg.svd(K, compute_uv=False)
    assert s[1] / s[0] < 1e-12


def test_two_throats_make_it_rank_two():
    a = branch_pair_matrix(2.3, 1.2, 0.9,
                           CoupledThroat(1.3, 1.0, +1, 0.3), n_k=6)
    b = branch_pair_matrix(2.3, 0.7, 1.6,
                           CoupledThroat(2.0, -0.4, -1, 0.22), n_k=6)
    s = np.linalg.svd(a["K"] + b["K"], compute_uv=False)
    assert s[1] / s[0] > 1e-3
    assert s[2] / s[0] < 1e-12


def test_each_entry_carries_the_phase_of_its_branch_pair():
    """``K_ab ∝ e^{−iω(ℓ_a + Δ + ℓ_b)}`` — which is why closure is coherence."""
    th = CoupledThroat(1.3, -2.1, +1, 0.3)
    w0, w1 = 1.0, 3.0
    g0 = branch_pair_matrix(w0, 1.2, 0.9, th, damping=0.0, n_k=3)
    g1 = branch_pair_matrix(w1, 1.2, 0.9, th, damping=0.0, n_k=3)
    for i, a in enumerate(g0["rows"]):
        for j, b in enumerate(g0["cols"]):
            lag = a["ell"] + th.delay + b["ell"]
            ratio = g1["K"][i, j] / g0["K"][i, j]
            assert ratio == pytest.approx(cmath.exp(-1j * (w1 - w0) * lag),
                                          abs=1e-10)


def test_a_closed_pair_is_frequency_independent():
    """PR #253's ``ℓ_a + Δ + ℓ_b = 0``, restated."""
    chi_in, chi_out = 1.2, 0.9
    th = CoupledThroat(1.3, -(chi_in + chi_out), +1, 0.3)
    vals = [branch_pair_matrix(w, chi_in, chi_out, th, damping=0.0,
                               n_k=3)["K"][0, 0]
            for w in (0.5, 2.0, 7.0)]
    assert vals[1] == pytest.approx(vals[0], abs=1e-12)
    assert vals[2] == pytest.approx(vals[0], abs=1e-12)


# ── the time domain ─────────────────────────────────────────────────────────
def test_history_words_arrive_where_their_lengths_say():
    th = CoupledThroat(1.3, 1.0, +1, 0.5)
    for w in coupled_arrivals(1.2, 0.9, th, 12.0, n_traversal=2, n_k=3):
        assert 0.0 <= w["t"] <= 12.0
        assert w["traversals"] >= 1
        assert w["sign"] in (+1, -1)
        assert w["amplitude"] > 0.0


def test_the_amplitude_ladder_is_kappa_to_the_number_of_traversals():
    th = CoupledThroat(1.3, 1.0, +1, 0.5)
    words = coupled_arrivals(1.2, 0.9, th, 40.0, n_traversal=3, n_k=2,
                             damping=0.0, rel_floor=0.0)
    a1 = 1.0 / (4.0 * math.pi * math.sin(1.2))
    a2 = 1.0 / (4.0 * math.pi * math.sin(0.9))
    ad = 1.0 / (4.0 * math.pi * math.sin(1.3))
    for w in words:
        n = w["traversals"] - 1
        assert w["amplitude"] == pytest.approx(
            0.5 ** (n + 1) * a1 * a2 * ad ** n, rel=1e-12)


def test_only_the_solve_produces_multi_traversal_words():
    """The control has one traversal and nothing else, by construction."""
    th = CoupledThroat(1.3, 1.0, +1, 0.5)
    words = coupled_arrivals(1.2, 0.9, th, 12.0, n_traversal=3, n_k=3)
    assert any(w["traversals"] > 1 for w in words)
    assert any(w["traversals"] == 1 for w in words)


def test_the_solved_waveform_is_the_sum_over_history_words():
    """The strongest statement the module makes about its own bookkeeping."""
    th = CoupledThroat(1.3, 1.0, +1, 0.6)
    width = 0.05
    ts, phi = coupled_waveform(11.0, 1.2, 0.9, th, 0.02, width)
    words = coupled_arrivals(1.2, 0.9, th, 14.0, n_traversal=5, n_k=5,
                             damping=0.02, rel_floor=1e-10)
    rec = np.zeros_like(ts)
    for w in words:
        rec += (w["sign"] * w["amplitude"]
                * np.exp(-(ts - w["t"]) ** 2 / (2.0 * width ** 2))
                / (width * math.sqrt(TWO_PI)))
    assert np.abs(rec - phi).max() / np.abs(phi).max() < 1e-4


def test_the_control_waveform_has_nothing_at_the_echo_time():
    th = CoupledThroat(1.3, 1.0, +1, 0.6)
    ts, solved = coupled_waveform(11.0, 1.2, 0.9, th, 0.02, 0.05)
    _, free = coupled_waveform(11.0, 1.2, 0.9, th, 0.02, 0.05,
                               resolvent=False)
    m = np.abs(ts - 5.4) < 0.1          # ℓ_a + ℓ_c + ℓ_b + 2Δ
    assert np.abs(solved[m]).max() > 1e5 * np.abs(free[m]).max()


def test_the_waveform_is_real_and_causal_before_the_first_arrival():
    th = CoupledThroat(1.3, 1.0, +1, 0.5)
    ts, phi = coupled_waveform(11.0, 1.2, 0.9, th, 0.02, 0.05)
    early = phi[ts < 2.5]               # first arrival is at 1.2 + 1.0 + 0.9
    assert np.abs(early).max() < 1e-4 * np.abs(phi).max()


# ── the critical coupling ───────────────────────────────────────────────────
def test_the_critical_coupling_is_where_the_loop_gain_reaches_one():
    c = critical_coupling(1.3, 0.03)
    th = CoupledThroat(1.3, 1.0, +1, c["kappa_critical"])
    assert abs(th.loop_transfer(c["omega_of_the_peak"], 0.03)) == (
        pytest.approx(1.0, rel=1e-9))


def test_the_critical_coupling_falls_linearly_with_the_regulator():
    a = critical_coupling(1.3, 0.04)["kappa_critical"]
    b = critical_coupling(1.3, 0.02)["kappa_critical"]
    assert a / b == pytest.approx(2.0, rel=1e-3)


def test_the_gain_peaks_on_an_esu_eigenfrequency():
    w = critical_coupling(1.3, 0.02)["omega_of_the_peak"]
    assert abs(w - round(w)) < 1e-3
    assert round(w) >= 1


# ── the measurements ────────────────────────────────────────────────────────
def test_measure_the_closed_form_transfer_is_the_branch_sum():
    r = measure_the_closed_form_transfer_is_the_branch_sum(
        chis=(1.3,), omegas=(1.7, 3.3), n_k=300)
    assert r["the_series_is_the_closed_form"]
    assert r["the_poles_are_the_esu_spectrum"]
    assert r["the_residues_are_the_mode_functions"]


def test_measure_solving_the_throat_resums_every_traversal():
    r = measure_solving_the_throat_resums_every_traversal(
        kappas=(0.2, 0.45), omegas=(1.1, 2.7))
    assert r["the_resolvent_is_the_traversal_sum"]
    assert r["worst_relative_miss_of_post_processing"] > 0.0
    for row in r["rows"]:
        assert row["relative_miss_of_one_traversal"] == pytest.approx(
            row["loop_gain"], rel=1e-9)


def test_measure_the_coupled_field_has_arrivals_the_free_branches_do_not():
    r = measure_the_coupled_field_has_arrivals_the_free_branches_do_not()
    assert r["the_waveform_is_the_sum_over_history_words"]
    assert r["every_isolated_echo_stands_above_the_control"]
    assert r["every_echo_carries_its_maslov_word_sign"]
    assert r["isolated_echoes"]
    ladder = {row["traversals"]: row["over_first"] for row in
              r["amplitude_ladder"]}
    assert ladder[1] == pytest.approx(1.0)
    assert ladder[2] < 0.1 and ladder[3] < ladder[2]


def test_measure_closure_is_broadband_coherence():
    r = measure_closure_is_broadband_coherence()
    assert r["closed_pairs_are_broadband_coherent"]
    assert r["every_other_pair_dephases"]
    assert r["worst_closed_coherence"] == pytest.approx(1.0, abs=1e-9)


def test_the_closure_condition_does_not_factorize_but_the_amplitude_does():
    """Why the pair is the primitive."""
    r = measure_closure_is_broadband_coherence()
    assert r["n_closed"] == 3
    assert r["pairs_a_single_index_rule_would_select"] == 9
    assert r["the_condition_does_not_factorize"]
    assert r["the_amplitude_does_factorize"]


def test_measure_the_primitive_is_rank_one_for_one_throat_and_not_for_two():
    r = measure_the_primitive_is_rank_one_for_one_throat_and_not_for_two()
    assert r["one_throat_is_rank_one"]
    assert r["two_throats_are_rank_two"]
    assert r["rank_one_throat"] == 1 and r["rank_two_throats"] == 2
    assert r["cross_term_agrees"]
    assert r["the_cross_term_is_a_full_fringe"]
    assert "not the two-source invariant" in r["scope"]


def test_measure_the_expansion_fails_at_the_eigenfrequencies():
    r = measure_the_expansion_fails_at_the_eigenfrequencies()
    assert r["kappa_critical_scales_like_damping"]
    assert r["the_peak_sits_on_an_esu_eigenfrequency"]
    assert r["resonance_is_where_post_processing_is_worst"]
    assert r["mean_exponent"] == pytest.approx(1.0, abs=1e-2)


def test_the_module_says_what_it_still_puts_in():
    r = measure_the_primitive_is_rank_one_for_one_throat_and_not_for_two()
    assert "throat" in r["why_it_vanishes_without_a_second_throat"]
    from geometrodynamics.waves import branch_coupling
    doc = branch_coupling.__doc__ or ""
    assert "identification map" in doc
    assert "no backreaction" in doc.lower()
