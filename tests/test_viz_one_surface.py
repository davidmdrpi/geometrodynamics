"""One field on one surface — and where the plane-wave model stops applying.

The load-bearing checks here are of two kinds, and they are asserted
differently.  Identities (the trig collapse, the coincident cancellation, the
zonal bound) are asserted at machine precision because that is what an identity
means.  Statements about *which model the visualisation is in* — the zonal
optimum, the pulse plateau — are measurements, and are asserted at tolerances
that would fail if the finding reversed.
"""

import math

import numpy as np
import pytest

from geometrodynamics.viz import one_surface as osf
from geometrodynamics.viz.circle_slice import ANTIPODAL_TIME, TWO_PI


@pytest.fixture()
def surf():
    p = osf.OneSurfaceSlice(offset=math.pi, signs=osf.OPPOSED)
    p.slice_.reset()
    p.slice_.advance_to(ANTIPODAL_TIME)
    return p


# ── one field, one curve ────────────────────────────────────────────────────
def test_the_signs_must_be_plus_or_minus_one():
    with pytest.raises(ValueError):
        osf.OneSurfaceSlice(signs=(1.0, 0.0))


def test_there_is_exactly_one_curve(surf):
    """The whole point of the module: `radius` returns one array, not two."""
    r = surf.radius(gain=0.3)
    assert r.ndim == 1
    assert len(r) == len(surf.sigma)
    u = surf.field()
    mid = surf.slice_.bulk.r_mid
    assert np.max(np.abs(r - (mid + 0.3 * u))) < 1e-15


def test_the_field_is_the_sum_of_the_signed_contributions(surf):
    a, b = surf.contributions()
    assert np.max(np.abs(surf.field() - (a + b))) == 0.0


def test_the_bisector_is_a_node_of_the_opposed_field():
    """`u_A = u_B` on the equidistant axis, so the opposed field vanishes."""
    for f in (0.2, 0.5, 0.9):
        p = osf.OneSurfaceSlice(offset=f * math.pi, signs=osf.OPPOSED)
        p.slice_.reset()
        p.slice_.advance_to(0.5 * f * math.pi)
        assert abs(float(p.field(p.bisector)[0])) < 1e-13


# ── the repair costs nothing ────────────────────────────────────────────────
def test_v66s_separation_was_the_one_surface_field():
    """Same array, up to `R_mid`'s own rounding — so v66's numbers survive."""
    got = osf.measure_the_two_curve_picture_was_one_field_all_along()
    assert got["they_are_the_same_array"]
    assert got["residue_in_ulps"] <= 4.0
    assert got["the_residue_is_the_mid_radius_not_the_fields"]


def test_the_labels_swap_between_the_two_constructions():
    """v66 'like' is one-surface OPPOSED; v66 'opposed' is one-surface LIKE."""
    from geometrodynamics.viz.two_wave_slice import (BOTH_OUTWARD,
                                                     TwoWaveSlice)
    alpha = 0.5 * math.pi
    v66 = TwoWaveSlice(offset=alpha, signs=BOTH_OUTWARD)
    one = osf.OneSurfaceSlice(offset=alpha, signs=osf.OPPOSED)
    one.slice_ = v66.slice_
    v66.slice_.reset()
    v66.slice_.advance_to(ANTIPODAL_TIME)
    ulp = float(np.spacing(v66.slice_.bulk.r_mid))
    assert np.max(np.abs(v66.separation(gain=1.0, closed=True)
                         - one.field())) <= 4.0 * ulp


# ── coincidence cancels ─────────────────────────────────────────────────────
def test_coincident_opposed_foci_cancel_exactly():
    got = osf.measure_coincident_foci_cancel_exactly()
    assert got["the_opposed_field_is_identically_zero"]
    assert got["worst_absolute_field"] == 0.0
    assert got["no_gain_connects_it"]
    assert not math.isfinite(got["opposed_span_gain"])


def test_the_like_pair_is_unaffected_by_coincidence():
    """Only the opposed combination cancels; the like one reinforces."""
    got = osf.measure_coincident_foci_cancel_exactly()
    assert got["the_like_pair_is_unaffected"]
    assert math.isfinite(got["like_signed_cheapest_span_gain"])


# ── the monochromatic law ───────────────────────────────────────────────────
def test_the_trig_collapse_is_exact():
    """`u = -2A sin(m alpha/2) sin(m theta - omega t)`, checked on a grid."""
    th = -math.pi + (np.arange(4001) + 0.5) * TWO_PI / 4001
    for m in (1, 2, 3, 5, 8):
        for f in (0.13, 0.5, 0.87, 1.0):
            alpha = f * math.pi
            p = osf.MonochromaticPair(mode=m, offset=alpha)
            for t in (0.0, 0.7, 2.3):
                got = p.field(th, t=t)
                want = -2.0 * math.sin(0.5 * m * alpha) * np.sin(m * th - t)
                assert np.max(np.abs(got - want)) < 1e-12


def test_the_amplitude_factor_is_the_measured_amplitude():
    got = osf.measure_the_optimum_is_half_a_wavelength()
    assert got["the_closed_form_is_the_measured_amplitude"]
    assert got["worst_amplitude_error"] < 1e-6


def test_the_first_optimum_is_half_a_wavelength():
    got = osf.measure_the_optimum_is_half_a_wavelength()
    assert got["the_first_optimum_is_half_a_wavelength"]
    for r in got["rows"]:
        assert abs(r["first_optimum_over_pi"] - 1.0 / r["mode"]) < 1e-12
        assert abs(r["wavelength_over_pi"] - 2.0 / r["mode"]) < 1e-12


def test_the_required_amplitude_diverges_at_coincidence():
    gap = 0.52
    assert not math.isfinite(
        osf.MonochromaticPair(mode=1, offset=0.0).required_amplitude(gap))
    at_pi = osf.MonochromaticPair(mode=1, offset=math.pi)
    assert abs(at_pi.required_amplitude(gap) - gap / 4.0) < 1e-15
    for deg, want in ((90, 0.183848), (60, 0.260000), (30, 0.502281)):
        p = osf.MonochromaticPair(mode=1, offset=math.radians(deg))
        assert abs(p.required_amplitude(gap) - want) < 1e-6


def test_the_extrema_do_not_move_with_the_offset():
    """A plane wave has no centre, so the foci do not pull the extrema."""
    for m in (1, 2, 3):
        base = osf.MonochromaticPair(mode=m, offset=0.3).extrema()
        for f in (0.5, 0.9, 1.0):
            got = osf.MonochromaticPair(mode=m,
                                        offset=f * math.pi).extrema()
            assert np.max(np.abs(np.array(got) - np.array(base))) < 1e-12


# ── the parity ──────────────────────────────────────────────────────────────
def test_opposite_foci_are_maximal_for_odd_modes_and_cancel_for_even():
    got = osf.measure_the_antipode_is_parity_dependent()
    assert got["odd_modes_are_maximal_at_the_antipode"]
    assert got["even_modes_cancel_at_the_antipode"]


def test_the_parity_is_in_the_real_spectrum_not_the_reduction():
    """`Z_n(pi) = (-1)^n` on S^3 gives the same odd/even split, exactly."""
    got = osf.measure_the_antipode_is_parity_dependent()
    assert got["the_zonal_spectrum_agrees"]
    assert got["zonal_odd_doubles"]
    assert got["zonal_even_cancels"]


def test_the_zonal_harmonic_is_bounded_and_has_both_endpoints_right():
    """The endpoint guard is load-bearing: `sin chi` vanishes at `pi` too."""
    chi = np.linspace(1e-9, math.pi - 1e-9, 20001)
    for n in range(0, 10):
        z = osf.zonal_harmonic(n, chi)
        assert np.max(np.abs(z)) <= 1.0 + 1e-9, f"|Z_{n}| must not exceed 1"
        assert abs(float(osf.zonal_harmonic(n, np.array([0.0]))[0]) - 1.0) < 1e-12
        assert abs(float(osf.zonal_harmonic(n, np.array([math.pi]))[0])
                   - (-1.0) ** n) < 1e-12


def test_the_zonal_harmonic_solves_its_eigenvalue_problem():
    """`f'' + 2 cot(chi) f' = -n(n+2) f`, away from the nodes."""
    h = 1e-5
    chi = np.linspace(0.4, math.pi - 0.4, 41)
    for n in (1, 2, 3, 5, 8):
        f0 = osf.zonal_harmonic(n, chi)
        fp = osf.zonal_harmonic(n, chi + h)
        fm = osf.zonal_harmonic(n, chi - h)
        # the S^3 zonal Laplacian is f'' + 2 cot(chi) f' -- the 2 matters,
        # and dropping it turns the eigenvalue into -(n^2 + 2n - 1)
        lap = ((fp - 2 * f0 + fm) / h**2
               + 2.0 * (fp - fm) / (2 * h) / np.tan(chi))
        ok = np.abs(f0) > 0.05
        assert np.max(np.abs(lap[ok] / f0[ok] + n * (n + 2))) < 2e-3


# ── where the two models part company ───────────────────────────────────────
def test_the_zonal_optimum_is_the_antipode_not_half_a_wavelength():
    """The measured disagreement with the plane-wave picture."""
    got = osf.measure_the_zonal_optimum_is_the_antipode()
    assert got["the_bound_is_two"]
    assert got["odd_orders_peak_at_the_antipode"]
    assert got["and_they_saturate_the_bound"]
    assert got["half_a_wavelength_does_not_reproduce_it"]


def test_even_zonal_orders_never_reach_the_bound():
    got = osf.measure_the_zonal_optimum_is_the_antipode()
    assert got["even_orders_cancel_at_the_antipode"]
    assert got["and_never_reach_the_bound"]
    assert got["best_even_strength"] < 1.4


def test_the_even_maximiser_location_is_not_unique():
    """So the strength is reported as a result and the location as a caveat.

    A 1501-point sweep and a 601-point sweep disagree about where `n = 6`
    peaks (`0.794π` against `0.205π`) while agreeing on the height to seven
    digits. Reporting the location without the degeneracy would have shipped a
    number that moves with the grid.
    """
    got = osf.measure_the_zonal_optimum_is_the_antipode()
    assert got["the_odd_maximiser_is_unique"]
    assert got["the_even_maximiser_is_not"]


def test_the_zonal_bound_of_two_is_attained_only_antipodally_and_odd():
    """`|Z_A - Z_B| <= 2` needs `Z_A = +1` and `Z_B = -1` at one point."""
    for n in (1, 3, 5):
        assert abs(osf.zonal_pair_strength(n, math.pi) - 2.0) < 1e-4
        for f in (0.2, 0.5, 0.8):
            assert osf.zonal_pair_strength(n, f * math.pi) < 2.0 - 1e-3
    for n in (2, 4, 6):
        assert osf.zonal_pair_strength(n, math.pi) < 1e-6


def test_a_pulse_saturates_where_a_mode_diverges():
    got = osf.measure_a_pulse_saturates_where_a_mode_diverges()
    assert got["the_pulse_threshold_saturates"]
    assert got["the_monochromatic_law_keeps_falling"]
    assert got["it_still_cancels_at_coincidence"]
    assert got["plateau_spread"] / got["plateau_threshold"] < 0.01


def test_the_pulse_cancellation_window_is_about_its_own_width():
    got = osf.measure_a_pulse_saturates_where_a_mode_diverges()
    assert 0.3 < got["cancellation_window_in_pulse_widths"] < 4.0


# ── the chord, and what it costs ────────────────────────────────────────────
def test_the_chord_is_the_law_of_cosines():
    got = osf.measure_the_chord_shrinks_with_frequency()
    assert got["the_closed_form_is_the_law_of_cosines"]
    assert got["worst_closed_form_error"] < 1e-12


def test_the_chord_shrinks_to_the_radial_gap():
    got = osf.measure_the_chord_shrinks_with_frequency()
    assert got["the_chord_shrinks_monotonically"]
    assert abs(got["chord_at_the_fundamental"] - 2.0) < 1e-9
    assert got["chord_at_the_highest_mode"] > got["limit_is_the_radial_gap"]
    far = osf.bulk_chord(0.74, 1.26, 1e-6)
    assert abs(far - 0.52) < 1e-9


def test_the_span_is_the_same_across_the_half_wavelength_family():
    got = osf.measure_the_chord_shrinks_with_frequency()
    assert got["the_span_is_the_same_at_every_mode"]
    for r in got["rows"]:
        assert abs(r["radial_span_over_A"] - 2.0) < 1e-6


def test_fixed_energy_is_not_fixed_amplitude():
    got = osf.measure_fixed_energy_reverses_the_preference()
    assert got["span_is_flat_at_fixed_amplitude"]
    assert got["span_falls_like_one_over_omega"]
    assert got["the_chord_and_the_span_both_shrink"]
    assert got["highest_mode_that_still_spans_at_fixed_energy"] is not None


# ── the surface reaching both boundaries ────────────────────────────────────
def test_the_span_gain_is_set_by_the_smaller_reach(surf):
    out, inn = surf.reach()
    assert abs(surf.span_gain() - 0.5 * surf.gap / min(out, inn)) < 1e-15


def test_at_the_span_gain_the_curve_touches_both_boundaries(surf):
    got = surf.connection(gain=surf.span_gain())
    assert got["connected"]
    assert got["reaches_r_outer"] and got["reaches_r_inner"]
    b = surf.slice_.bulk
    assert max(abs(got["radius_out"] - b.r_outer),
               abs(got["radius_in"] - b.r_inner)) < 1e-9


def test_a_one_signed_field_can_never_span(surf):
    """It must reach *both* boundaries; reaching one twice is not a spanning."""
    p = osf.OneSurfaceSlice(offset=0.0, signs=osf.OPPOSED)
    p.slice_.reset()
    p.slice_.advance_to(ANTIPODAL_TIME)
    assert not math.isfinite(p.span_gain())


# ── the decomposition: what sits where on the one surface ───────────────────
def test_the_decomposition_sums_to_the_field(surf):
    d = osf.decompose(surf)
    assert np.max(np.abs(d["contribution_a"] + d["contribution_b"]
                         - d["field"])) == 0.0
    assert len(d["sigma"]) == len(surf.sigma)


def test_reinforcing_means_the_contributions_share_a_sign(surf):
    d = osf.decompose(surf)
    prod = d["contribution_a"] * d["contribution_b"]
    assert np.all(prod[d["reinforcing"]] > 0)
    assert np.all(prod[~d["reinforcing"]] <= 0)


def test_the_contributions_barely_overlap_once_the_foci_are_apart():
    """The offset does not turn interference on -- it turns cancellation off."""
    got = osf.measure_how_one_surface_answers_two_fronts()
    assert got["the_contributions_barely_overlap"]
    lo, hi = got["amplification_when_apart"]
    assert 1.0 < lo and hi < 1.05


def test_the_field_peak_sits_on_one_contributions_peak():
    """Two near-independent dents, so the total peaks where a single one does."""
    got = osf.measure_how_one_surface_answers_two_fronts()
    assert got["the_field_peak_sits_on_a_contribution_peak"]


def test_only_coincidence_gives_a_total_overlap_and_it_cancels():
    got = osf.measure_how_one_surface_answers_two_fronts()
    assert got["coincidence_is_the_only_total_overlap"]
    assert got["rows"][0]["peak_field"] == 0.0
    assert got["rows"][0]["amplification"] == 0.0
