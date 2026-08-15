"""
Tests for the solved field on the Einstein static universe.

The load-bearing claims:

* the **conformal** frequencies are exactly ``n+1``, which is why the Green
  function is periodic and its support is an image set — the branches are exact,
  not asymptotic;
* a truncated **mode** sum equals the closed-form **image** sum, so the solve
  really is the thing the ray ledger was describing;
* the peaks of the solved field land on PR #253's branch times, with amplitude
  ``1/(4π sin χ)`` — PR #251's shell law, derived here rather than imposed;
* and every arrival carries the **Maslov sign** ``(−1)^crossings``, which is the
  one thing a path-length ledger structurally cannot produce.

Also pinned: the sharp ledger belongs to the *conformally* coupled field. The
minimally coupled one has irrational frequencies and no image structure, and the
test asserts the contrast rather than leaving it implied.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.field_solve import (
    TWO_PI,
    branch_arrivals,
    esu_frequency,
    field_peaks,
    image_field,
    measure_minimal_coupling_has_no_branch_structure,
    measure_the_amplitude_is_the_shell_law,
    measure_the_field_support_is_the_branch_ledger,
    measure_the_phases_are_the_maslov_index,
    measure_the_spectral_solve_matches_the_image_sum,
    measure_the_throat_reproduces_the_closure_condition,
    spectral_field,
    through_throat_arrivals,
)


# ── the spectrum, which is why any of this is sharp ─────────────────────────
@pytest.mark.parametrize("n", [0, 1, 2, 5, 11])
def test_conformal_frequencies_are_the_integers(n):
    assert esu_frequency(n, conformal=True) == pytest.approx(n + 1.0, abs=1e-12)


@pytest.mark.parametrize("n", [1, 2, 5, 11])
def test_minimal_frequencies_are_irrational(n):
    om = esu_frequency(n, conformal=False)
    assert om == pytest.approx(math.sqrt(n * (n + 2)))
    assert abs(om - round(om)) > 1e-3          # never an integer


def test_the_zero_mode_is_massless_only_without_conformal_coupling():
    assert esu_frequency(0, conformal=False) == 0.0
    assert esu_frequency(0, conformal=True) == 1.0


# ── the solve against the image sum ─────────────────────────────────────────
@pytest.mark.parametrize("chi", [0.4, 1.1, 2.0, 2.7])
def test_mode_sum_equals_image_sum(chi):
    for t in np.linspace(0.05, 4.0 * math.pi, 25):
        assert spectral_field(chi, float(t)) == pytest.approx(
            image_field(chi, float(t)), abs=1e-10)


def test_the_field_is_singular_at_the_poles():
    for bad in (0.0, math.pi):
        with pytest.raises(ValueError):
            spectral_field(bad, 1.0)
        with pytest.raises(ValueError):
            image_field(bad, 1.0)


# ── the branch ledger ───────────────────────────────────────────────────────
def test_branch_arrivals_are_the_ray_ledger():
    chi = 1.1
    got = branch_arrivals(chi, 4.0 * math.pi)
    times = [r["t"] for r in got]
    assert times == sorted(times)
    for k in (0, 1):
        assert any(abs(t - (chi + TWO_PI * k)) < 1e-12 for t in times)
        assert any(abs(t - (TWO_PI * (k + 1) - chi)) < 1e-12 for t in times)


def test_the_maslov_sign_alternates_with_focal_crossings():
    for r in branch_arrivals(1.1, 6.0 * math.pi):
        assert r["sign"] == (-1) ** r["crossings"]
        # short way with k windings crosses 2k foci; the long way, 2k+1
        assert r["crossings"] == (2 * r["winding"] + (1 if r["long_way"]
                                                      else 0))


@pytest.mark.parametrize("chi", [0.7, 1.1, 1.9])
def test_solved_peaks_land_on_the_branch_times_with_the_right_sign(chi):
    t_max = 4.0 * math.pi
    peaks = field_peaks(chi, t_max, width=0.05)
    ledger = branch_arrivals(chi, t_max)
    assert len(peaks) == len(ledger)
    for p in peaks:
        best = min(ledger, key=lambda r: abs(r["t"] - p["t"]))
        assert abs(best["t"] - p["t"]) < 1e-2          # grid-limited
        assert (1 if p["value"] > 0 else -1) == best["sign"]


# ── the amplitude ───────────────────────────────────────────────────────────
def test_peak_amplitude_is_one_over_four_pi_sin_chi():
    width = 0.05
    expect = 1.0 / (4.0 * math.pi * width * math.sqrt(TWO_PI))
    for chi in (0.35, 1.1, math.pi / 2, 2.6):
        peak = spectral_field(chi, chi, width)
        assert peak * math.sin(chi) == pytest.approx(expect, rel=1e-9)


# ── the throat ──────────────────────────────────────────────────────────────
def test_through_throat_arrivals_are_leg_plus_delay_plus_leg():
    a, b, delay = 1.2, 0.9, -2.1
    got = through_throat_arrivals(a, b, delay, 4.0 * math.pi)
    assert got
    for r in got:
        assert r["t"] == pytest.approx(r["leg_in"] + delay + r["leg_out"])
        assert r["sign"] in (+1, -1)


def test_closure_puts_an_arrival_back_on_the_emission_event():
    """PR #253's ℓ₁ + Δ + ℓ₂ = 0, now as a statement about a field arrival."""
    a, b = 1.2, 0.9
    delay = -(a + b)                     # the closure condition
    got = through_throat_arrivals(a, b, delay, 4.0 * math.pi)
    assert any(abs(r["t"]) < 1e-9 for r in got)


def test_the_orientation_flips_the_closing_sign():
    a, b = 1.2, 0.9
    delay = -(a + b)
    plus = [r for r in through_throat_arrivals(a, b, delay, 1.0, +1)
            if abs(r["t"]) < 1e-9][0]
    minus = [r for r in through_throat_arrivals(a, b, delay, 1.0, -1)
             if abs(r["t"]) < 1e-9][0]
    assert plus["sign"] == -minus["sign"]


# ── the measurements ────────────────────────────────────────────────────────
def test_measure_the_spectral_solve_matches_the_image_sum():
    r = measure_the_spectral_solve_matches_the_image_sum(chis=(0.6, 1.6),
                                                         n_t=20)
    assert r["the_two_constructions_agree"]
    assert r["worst_abs_error"] < 1e-10


def test_measure_the_field_support_is_the_branch_ledger():
    r = measure_the_field_support_is_the_branch_ledger(chis=(1.1,),
                                                       t_max=2.5 * math.pi)
    assert r["so_the_field_reproduces_the_ray_ledger"]
    assert r["every_branch_has_a_peak_and_no_peak_is_spurious"]


def test_measure_the_phases_are_the_maslov_index():
    """The one thing the ray ledger could not have supplied."""
    r = measure_the_phases_are_the_maslov_index(chis=(1.1,),
                                                t_max=2.5 * math.pi)
    assert r["every_sign_is_the_maslov_factor"]
    assert r["rows"]


def test_measure_the_amplitude_is_the_shell_law():
    r = measure_the_amplitude_is_the_shell_law()
    assert r["the_product_is_constant"]
    assert r["matches_one_over_four_pi_sin_chi"]
    assert r["relative_spread"] < 1e-9


def test_measure_minimal_coupling_has_no_branch_structure():
    r = measure_minimal_coupling_has_no_branch_structure(n_grid=1200)
    assert r["conformal_is_sharp"]
    assert r["minimal_is_not"]
    assert r["the_ledger_belongs_to_the_conformal_field"]
    assert r["minimal"]["ratio"] > 100.0 * r["conformal"]["ratio"]


def test_measure_the_throat_reproduces_the_closure_condition():
    r = measure_the_throat_reproduces_the_closure_condition()
    assert r["closure_puts_an_arrival_on_the_emission_event"]
    assert r["every_closing_sign_is_eta_times_the_maslov_factors"]
    assert "identification map" in r["what_is_still_put_in"]
