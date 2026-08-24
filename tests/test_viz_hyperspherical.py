"""What higher dimension does to the bulk/intersection picture.

Two kinds of check here. The measure and orientation facts are exact and are
asserted against closed forms or against a determinant; the Monte-Carlo ones
(random directions) are asserted only where the law predicts a scaling, never
on a single sampled number.
"""

import math

import numpy as np
import pytest

from geometrodynamics.viz import hyperspherical as hs
from geometrodynamics.viz import regularized_center as rc


# ── the peak, and the unit it hides ─────────────────────────────────────────
def test_the_ball_peaks_at_five_and_the_sphere_two_later():
    """The familiar claim is about the *ball*; the sphere peaks at `S⁶ ⊂ ℝ⁷`."""
    assert max(range(1, 60), key=hs.ball_volume) == 5
    assert max(range(1, 60), key=hs.sphere_area) == 7
    got = hs.measure_the_peak_is_the_ball_not_the_sphere_and_needs_a_unit()
    assert got["the_ball_peaks_at_five"]
    assert got["the_sphere_peaks_two_later"]
    assert got["the_peaking_sphere_is"] == "S^6 in R^7"


def test_the_peak_is_not_a_fact_about_dimension_alone():
    """Doubling `R` moves the ball's peak from `d = 5` to `d = 24`.

    `V_d` and `A_{d−1}` carry different units at different `d`, so comparing
    across `d` picks a length scale. Any "spheres peak at N dimensions" claim
    is a claim about one radius.
    """
    assert max(range(1, 80), key=lambda d: hs.ball_volume(d, 0.5)) == 1
    assert max(range(1, 80), key=lambda d: hs.ball_volume(d, 1.0)) == 5
    assert max(range(1, 80), key=lambda d: hs.ball_volume(d, 2.0)) == 24
    got = hs.measure_the_peak_is_the_ball_not_the_sphere_and_needs_a_unit()
    assert got["the_peak_moves_with_the_unit"]
    peaks = {s["ball_peak_dim"] for s in got["radius_scan"]}
    assert len(peaks) == len(got["radius_scan"]), "every radius its own peak"


@pytest.mark.parametrize("d,expect", [(1, 2.0), (2, math.pi), (3, 4.0 * math.pi / 3.0)])
def test_the_ball_volume_matches_the_elementary_cases(d, expect):
    assert abs(hs.ball_volume(d) - expect) < 1e-12


@pytest.mark.parametrize("d,expect", [(2, 2.0 * math.pi), (3, 4.0 * math.pi),
                                      (4, 2.0 * math.pi ** 2)])
def test_the_sphere_area_matches_the_elementary_cases(d, expect):
    """`|S¹| = 2π`, `|S²| = 4π`, `|S³| = 2π²`."""
    assert abs(hs.sphere_area(d) - expect) < 1e-12


def test_the_ball_and_sphere_are_related_by_the_radial_integral():
    """`V_d = ∫₀¹ A_{d−1}(r) dr`, an independent tie between the two."""
    from scipy.integrate import quad
    for d in (2, 3, 5, 7):
        got = quad(lambda r: hs.sphere_area(d, r), 0.0, 1.0, limit=200)[0]
        assert abs(got / hs.ball_volume(d) - 1.0) < 1e-10


# ── concentration at the equator ────────────────────────────────────────────
def test_the_shell_measure_concentrates_with_width_one_over_root_n():
    got = hs.measure_the_shell_measure_concentrates_at_the_equator()
    assert got["the_mean_is_always_the_equator"]
    assert got["std_times_sqrt_n_tends_to_one"]
    assert got["the_band_narrows_as_one_over_sqrt_n"]
    tail = [r for r in got["rows"] if r["sphere_dim"] >= 200]
    for row in tail:
        assert abs(row["std_times_sqrt_n"] - 1.0) < 1e-4


def test_the_gaussian_approximation_to_the_shell_measure():
    """`sin^{n−1}(π/2+δ) ≈ e^{−(n−1)δ²/2}` — the reason for the `1/√n`."""
    for n in (50, 200, 1000):
        for delta in (0.2 / math.sqrt(n), 0.5 / math.sqrt(n)):
            exact = hs.shell_measure(math.pi / 2 + delta, n)
            approx = math.exp(-(n - 1) * delta ** 2 / 2.0)
            assert abs(exact / approx - 1.0) < 1e-2


def test_the_shell_measure_is_symmetric_about_the_equator():
    for n in (3, 7, 40):
        for delta in (0.1, 0.4, 1.0):
            a = hs.shell_measure(math.pi / 2 - delta, n)
            b = hs.shell_measure(math.pi / 2 + delta, n)
            assert abs(a - b) < 1e-14 * max(a, 1e-300)


# ── the antipode is non-generic ─────────────────────────────────────────────
def test_random_directions_are_nearly_orthogonal():
    got = hs.measure_the_antipode_is_a_vanishing_measure_relation(
        samples=60000, seed=11)
    assert got["std_scales_as_one_over_sqrt_n"]
    assert got["the_mean_angle_is_a_right_angle"]
    for row in got["rows"]:
        assert abs(row["std_dot_times_sqrt_n"] - 1.0) < 0.03


def test_the_near_antipodal_fraction_collapses_with_dimension():
    """Asserted as a *trend*, not on one sampled number."""
    got = hs.measure_the_antipode_is_a_vanishing_measure_relation(
        samples=60000, seed=5)
    fracs = [r["fraction_near_antipodal"] for r in got["rows"]]
    assert fracs[0] > 0.0, "S^2 does have near-antipodal pairs"
    assert all(b <= a for a, b in zip(fracs, fracs[1:])), "monotone"
    assert fracs[-1] == 0.0, "and none at all by n = 500"
    assert got["near_antipodal_fraction_collapses"]


def test_the_measurement_says_what_the_non_genericity_does_not_prove():
    got = hs.measure_the_antipode_is_a_vanishing_measure_relation(
        samples=20000, seed=2)
    caveat = got["what_is_not_claimed"].lower()
    assert "correct" in caveat


# ── the wavefront ───────────────────────────────────────────────────────────
def test_the_front_vanishes_at_both_poles_and_peaks_at_the_equator():
    """The antipode is machine-zero *relative to the peak*, not absolutely.

    At `n = 2` the exponent is 1, so `sin(π)` is not squared and the value is
    `7.7e-16` — small because `sin(π)` is, not because the law drives it to
    zero. The meaningful statement is the ratio to the equator.
    """
    for n in (2, 3, 5):
        peak = hs.front_measure(math.pi / 2, n)
        assert hs.front_measure(0.0, n) == 0.0
        assert hs.front_measure(math.pi, n) < 1e-14 * peak
        for chi in (0.2, 0.9, 2.2, 2.9):
            assert hs.front_measure(chi, n) < peak


def test_the_front_on_s3_is_the_sin_squared_law():
    """`S³` is the repo's case: front `∝ sin²χ`."""
    for chi in (0.3, 1.0, 2.4):
        a = hs.front_measure(chi, 3)
        b = hs.sphere_area(3) * math.sin(chi) ** 2
        assert abs(a - b) < 1e-12


def test_the_middle_dominates_more_as_dimension_rises():
    got = hs.measure_the_front_is_a_point_then_a_shell_then_a_point()
    assert got["the_peak_is_always_the_equator"]
    assert got["the_poles_vanish"]
    assert got["the_middle_dominates_more_with_dimension"]
    rows = {r["sphere_dim"]: r for r in got["rows"]}
    assert rows[2]["measure_within_a_tenth_of_the_equator"] < 0.11
    assert rows[7]["measure_within_a_tenth_of_the_equator"] > 0.19


@pytest.mark.parametrize("radius", [0.5, 1.0, 3.0])
def test_the_front_scales_with_the_radius_as_the_power(radius):
    for n in (2, 3, 4):
        a = hs.front_measure(1.0, n, radius)
        b = hs.front_measure(1.0, n, 1.0) * radius ** (n - 1)
        assert abs(a / b - 1.0) < 1e-12


# ── the collapse that changes the picture ───────────────────────────────────
@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_the_patch_collapses_as_the_nth_power(n):
    assert abs(hs.patch_collapse(1e-3, 1.0, n) - 10.0 ** (-3 * n)) < 1e-18


def test_the_two_dimensional_drawing_understates_the_collapse():
    """The headline correction: same angular footprint, `f₀ⁿ` physical overlap."""
    got = hs.measure_the_patch_collapses_as_f_to_the_n()
    assert got["the_angular_footprint_is_untouched"]
    assert got["the_drawn_2d_length_understates_it"]
    assert got["the_collapse_steepens_with_dimension"]
    assert abs(got["understatement_factor_at_n_three"] - 1e6) < 1.0
    assert "ribbons" in got["what_the_picture_is_of"]
    assert "angular question" in got["the_yes_no_criterion_is_unaffected"]


def test_the_collapse_exponent_is_the_same_weight_regularized_center_uses():
    """`f^n` here is the `f^q` the previous round's review forced.

    `monopole_resistance` weights by the volume element `f^q`; the overlap
    region is weighted the same way. One fact, two consequences — so they must
    agree on the exponent.
    """
    f0, F = 1e-3, 0.5
    for q in (1, 2, 3, 4):
        # the resistance integrand carries f^{-q}; the patch carries f^{+q}
        near = rc.monopole_resistance(F, f0, q)
        assert near > 0.0
        assert abs(math.log(hs.patch_collapse(f0, F, q))
                   / math.log(f0 / F) - q) < 1e-12


def test_the_angular_criterion_is_untouched_by_the_collapse():
    """Whether two fronts meet never involved `f₀`, and still does not."""
    got = rc.measure_the_bearing_replaces_collision_with_overlap()
    assert got["the_verdict_does_not_depend_on_the_neck"]
    for row in got["rows"]:
        reach = 0.5 * (row["angular_width_a"] + row["angular_width_b"])
        assert row["they_meet"] == (row["angular_separation"] < reach)


# ── orientability is a parity effect ────────────────────────────────────────
@pytest.mark.parametrize("n,expect", [(1, True), (2, False), (3, True),
                                      (4, False), (5, True), (6, False)])
def test_projective_space_is_orientable_iff_n_is_odd(n, expect):
    assert hs.projective_space_is_orientable(n) is expect
    assert (round(np.linalg.det(-np.eye(n + 1))) > 0) is expect


def test_the_two_quotients_this_repo_uses_are_always_opposite():
    """Spatial `ℝP^d` and exchange `ℝP^{d−1}` are one apart, hence opposite.

    At `d = 3` the spatial quotient is orientable and the exchange space is
    not — and it is the exchange space that carries the Pin⁻ structure. Raising
    the spatial dimension swaps them.
    """
    got = hs.measure_the_quotient_flips_orientability_with_parity()
    assert got["orientable_iff_odd"]
    assert got["the_two_quotients_are_always_opposite"]
    assert got["raising_the_dimension_swaps_them"]
    here = {r["spatial_dim"]: r for r in got["rows"]}[3]
    assert here["spatial_orientable"] is True
    assert here["exchange_orientable"] is False
    assert here["exchange_quotient"] == "RP^2"
    nxt = {r["spatial_dim"]: r for r in got["rows"]}[4]
    assert nxt["spatial_orientable"] is False
    assert nxt["exchange_orientable"] is True


def test_two_dimensional_space_fails_differently():
    """`ℝP¹ ≃ S¹` has `π₁ = ℤ` — the braid group, not `ℤ₂`."""
    got = hs.measure_the_quotient_flips_orientability_with_parity()
    two = {r["spatial_dim"]: r for r in got["rows"]}[2]
    assert "braid" in two["exchange_pi1"]
    assert all("Z_2" == r["exchange_pi1"] for r in got["rows"]
               if r["spatial_dim"] >= 3)


def test_the_parity_claim_is_scoped_to_what_it_would_imply():
    got = hs.measure_the_quotient_flips_orientability_with_parity()
    assert got["it_is_parity_not_size"]
    caveat = got["what_is_not_claimed"].lower()
    assert "should move dimension" in caveat
    assert "re-derived" in caveat


# ── S³ is exceptional ───────────────────────────────────────────────────────
def test_s3_carries_a_global_orthonormal_frame():
    """Parallelizable: `q·i, q·j, q·k` are tangent, orthonormal, nowhere zero."""
    got = hs.measure_s3_is_not_a_generic_sphere(samples=4000, seed=1)
    assert got["the_frame_is_global_and_nowhere_zero"]
    assert got["frame_orthonormality_error"] < 1e-12
    assert got["smallest_frame_vector_norm"] > 0.999


def test_the_frame_is_tangent_to_the_sphere_at_every_point():
    rng = np.random.default_rng(4)
    q = rng.normal(size=(3000, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    frame = hs.quaternion_frame(q)
    assert frame.shape == (3000, 3, 4)
    radial = np.einsum("nik,nk->ni", frame, q)
    assert np.abs(radial).max() < 1e-13, "orthogonal to the radial direction"


def test_the_hopf_map_lands_on_the_sphere_with_circle_fibres():
    got = hs.measure_s3_is_not_a_generic_sphere(samples=4000, seed=6)
    assert got["the_fibre_is_a_circle_not_a_point"]
    assert got["hopf_image_on_the_sphere_error"] < 1e-12
    assert got["hopf_fibre_error"] < 1e-12
    assert got["median_distance_moved_along_a_fibre"] > 0.1


def test_the_hopf_image_covers_the_whole_two_sphere():
    """A weaker check than surjectivity, but it rules out a degenerate image."""
    rng = np.random.default_rng(9)
    q = rng.normal(size=(40000, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    h = hs.hopf_map(q)
    assert h[:, 2].min() < -0.99 and h[:, 2].max() > 0.99
    assert np.abs(np.linalg.norm(h, axis=1) - 1.0).max() < 1e-12


def test_the_module_says_nothing_extrapolates_smoothly():
    got = hs.measure_s3_is_not_a_generic_sphere(samples=2000, seed=8)
    assert "S^1, S^3 and S^7" in got["so_nothing_here_extrapolates_smoothly"]
    assert "hairy-ball" in got["s2_has_no_such_frame"]
    assert any("parallelizable" in s for s in got["identifications"])


# ── the interpretation, marked as one ───────────────────────────────────────
def test_the_bearing_reading_is_marked_as_an_interpretation():
    """It is a reading of PR #268, not a new result, and says so."""
    got = hs.measure_the_bearing_is_the_blown_up_direction_space()
    assert got["every_bearing_measure_is_positive"]
    assert "interpretation" in " ".join(got.keys())
    text = got["this_is_an_interpretation_not_a_result"]
    assert "stands on its own geometry" in text
    assert "nothing in" in text


def test_the_bearing_measure_is_the_direction_space_at_finite_size():
    got = hs.measure_the_bearing_is_the_blown_up_direction_space()
    for row in got["rows"]:
        expect = hs.sphere_area(row["angular_dim"] + 1) * \
            row["neck"] ** row["angular_dim"]
        assert abs(row["bearing_measure"] / expect - 1.0) < 1e-12
        assert row["bearing_measure"] > 0.0
    # and it collapses to nothing only in the limit
    tiny = hs.sphere_area(4) * 1e-30 ** 3
    assert tiny > 0.0


def test_the_peak_scan_reports_ties_instead_of_picking_a_side():
    """At `R = ½` the sphere measure is exactly tied between `d = 2` and `3`.

    `2πR = 4πR²` at `R = ½`, so a bare `max` picks a side by floating-point
    accident — and changing how the value is computed flips it, which is how
    this was noticed. The scan reports the tie.
    """
    got = hs.measure_the_peak_is_the_ball_not_the_sphere_and_needs_a_unit()
    half = next(s for s in got["radius_scan"] if s["radius"] == 0.5)
    assert half["sphere_tied_dims"] == [2, 3]
    assert abs(hs.sphere_area(2, 0.5) - hs.sphere_area(3, 0.5)) < 1e-15
    assert abs(hs.sphere_area(2, 0.5) - math.pi) < 1e-15
    # and the unit-radius peaks, which are the ones quoted, are not ties
    assert got["the_unit_peaks_are_unique"]
    unit = next(s for s in got["radius_scan"] if s["radius"] == 1.0)
    assert unit["ball_tied_dims"] == [5] and unit["sphere_tied_dims"] == [7]


def test_no_reported_peak_is_the_search_ceiling():
    """A peak equal to the top of the scanned range is a range, not a peak."""
    got = hs.measure_the_peak_is_the_ball_not_the_sphere_and_needs_a_unit()
    assert got["no_peak_is_at_the_search_ceiling"]
    big = next(s for s in got["radius_scan"] if s["radius"] == 4.0)
    assert big["ball_peak_dim"] == 100, "interior, and it moved from 79 when " \
        "the range was widened -- which is why the guard exists"


# ── which n belongs to which object ─────────────────────────────────────────
def test_n_is_the_dimension_of_the_objects_own_transverse_sphere():
    """The scope pin: `n` is a fact about the object, not a modelling choice.

    The millionfold figure is the `S³` **bearing**'s. PR #265's spatial throat
    has an `S²` cross-section, so its own understatement against the 2-D
    drawing is a thousand — and without this stated, the same `f^n` law
    migrates between objects that do not share an `n`.
    """
    got = hs.measure_which_n_is_physical_for_which_object()
    by_dim = {r["angular_dim"]: r for r in got["rows"]}
    assert set(by_dim) == {1, 2, 3, 4}
    for n, row in by_dim.items():
        assert row["cross_section"] == f"S^{n}"
        assert row["patch_measure"] == pytest.approx(
            got["squeeze"] ** n, rel=1e-12)
        assert row["understatement_vs_the_drawing"] == pytest.approx(
            got["squeeze"] / got["squeeze"] ** n, rel=1e-12)
    assert "#265" in by_dim[2]["object"]
    assert "4-spatial" in by_dim[3]["object"]
    assert by_dim[2]["understatement_vs_the_drawing"] == pytest.approx(
        1e3, rel=1e-9)
    assert by_dim[3]["understatement_vs_the_drawing"] == pytest.approx(
        1e6, rel=1e-9)


def test_the_265_throat_is_an_f_squared_object_measured_not_asserted():
    """`physical_throat`'s neck area is `4π f₀²`, which is what `n = 2` means."""
    got = hs.measure_which_n_is_physical_for_which_object()
    assert got["physical_throat_is_n_equals_two"]
    assert got["physical_throat_neck_area"] == pytest.approx(
        got["four_pi_f0_squared"], rel=1e-12)
    from geometrodynamics.waves.physical_throat import VacuumThroat
    throat = VacuumThroat(mouth_radius=0.05)
    f0 = throat.neck_radius()
    assert got["physical_throat_neck_area"] == pytest.approx(
        4.0 * math.pi * f0 ** 2, rel=1e-12)
    assert got["physical_throat_neck_area"] == pytest.approx(
        throat.neck_area(), rel=1e-15)


def test_the_million_is_not_attributed_to_the_throat():
    got = hs.measure_which_n_is_physical_for_which_object()
    assert got["the_million_is_not_the_throats"]
    assert got["the_throats_own_understatement"] == pytest.approx(1e3, rel=1e-9)
    assert "bearing" in got["the_millionfold_figure_belongs_to"]
    assert "migrate" in got["why_it_matters"]


def test_the_doc_and_the_renderer_do_not_quote_the_million_unattributed():
    """Every place the millionfold figure appears must name its object.

    Six overstatements reached six files in PR #268 before being caught; this
    keeps the count at zero for this one by checking the text, not the number.
    """
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[1]
    for rel in ("docs/hyperspherical.md", "README.md",
                "scripts/geometrodynamics_v70_hyperspherical.py",
                "geometrodynamics/viz/hyperspherical.py"):
        text = (root / rel).read_text()
        for line_no, line in enumerate(text.splitlines(), 1):
            low = line.lower()
            if "million" not in low:
                continue
            window = "\n".join(text.splitlines()[
                max(0, line_no - 6):line_no + 6]).lower()
            assert "bearing" in window or "s^3" in window or "s³" in window, \
                f"{rel}:{line_no} quotes a millionfold figure without naming " \
                f"the object it belongs to"


# ── the bearing as a routing manifold ───────────────────────────────────────
def test_direction_capacity_is_the_reciprocal_cap_fraction():
    from scipy.integrate import quad
    for n in (1, 2, 3, 5):
        eps = 0.3
        cap = quad(lambda c: math.sin(c) ** (n - 1), 0.0, eps)[0]
        whole = quad(lambda c: math.sin(c) ** (n - 1), 0.0, math.pi)[0]
        assert hs.direction_capacity(eps, n) == pytest.approx(
            whole / cap, rel=1e-10)


def test_direction_capacity_grows_like_one_over_sin_eps_to_the_n():
    """Exponential in the dimension at fixed resolution."""
    eps = math.radians(20.0)
    caps = [hs.direction_capacity(eps, n) for n in range(1, 9)]
    ratios = [b / a for a, b in zip(caps, caps[1:])]
    assert all(r > 1.0 for r in ratios)
    # the asymptotic ratio is 1/sin(eps) = 2.924
    assert ratios[-1] == pytest.approx(1.0 / math.sin(eps), rel=0.10)
    # and it is still approaching from above, not settled
    assert ratios[-1] < ratios[0]


def test_direction_capacity_rejects_a_non_angle():
    for bad in (0.0, -0.2, math.pi, 4.0):
        with pytest.raises(ValueError):
            hs.direction_capacity(bad, 3)


def test_capacity_is_dimensionless_while_the_measure_collapses():
    """The two are independent, and the scan is a real one, not a tautology.

    An earlier draft of this check compared a value to itself. It varies the
    neck over six decades and confirms the proper measure moves while the
    capacity returns the identical float.
    """
    got = hs.measure_the_bearing_is_a_routing_manifold_not_a_hub()
    scan = got["neck_scan"]
    caps = {r["capacity_at_20_deg"] for r in scan}
    assert len(caps) == 1, "capacity must not depend on the neck at all"
    measures = [r["proper_measure"] for r in scan]
    assert measures[0] / measures[-1] > 1e17
    for row in scan:
        assert row["proper_measure"] == pytest.approx(
            hs.sphere_area(4) * row["neck"] ** 3, rel=1e-12)
    assert got["capacity_does_not_involve_the_neck"]
    assert got["proper_measure_spans_decades_meanwhile"]


def test_the_routing_reading_is_a_capacity_claim_not_a_crowding_one():
    got = hs.measure_the_bearing_is_a_routing_manifold_not_a_hub()
    assert got["capacity_at_20_deg_on_s3"] == pytest.approx(113.529, rel=1e-4)
    assert got["capacity_at_20_deg_on_s20"] > 1e10
    assert got["capacity_grows_exponentially"]
    assert got["a_thousand_directions_stay_near_orthogonal"]
    assert "not a hub" in got["the_reading"]


def test_high_dimensional_direction_families_stay_near_orthogonal():
    """A compressed direction space is not a crowded one."""
    got = hs.measure_the_bearing_is_a_routing_manifold_not_a_hub()
    fams = {f["ambient_dim"]: f for f in got["near_orthogonal_families"]}
    assert fams[1000]["max_pairwise_cos_of_1024"] < 0.2
    # and it is monotone in the dimension, which is the 1/sqrt(n) law of section 2
    dims = sorted(fams)
    worst = [fams[d]["max_pairwise_cos_of_1024"] for d in dims]
    assert all(b <= a for a, b in zip(worst, worst[1:]))


# ── what the limit separates ────────────────────────────────────────────────
def test_the_limit_keeps_the_angles_and_takes_only_the_measure():
    """Angular incidence survives; the proper interaction region does not."""
    got = hs.measure_the_limit_separates_three_things()
    angles = {r["overlap_angle"] for r in got["rows"]}
    assert len(angles) == 1, "the angular overlap must not depend on f0"
    assert all(r["they_meet"] for r in got["rows"])
    for r in got["rows"]:
        assert r["proper_overlap_n_2"] == pytest.approx(
            (r["overlap_angle"] * r["neck"]) ** 2, rel=1e-12)
        assert r["proper_overlap_n_3"] == pytest.approx(
            (r["overlap_angle"] * r["neck"]) ** 3, rel=1e-12)
    first, last = got["rows"][0], got["rows"][-1]
    assert first["proper_overlap_n_3"] / last["proper_overlap_n_3"] > 1e20
    assert got["angular_incidence_survives"]
    assert got["the_overlap_verdict_survives"]
    assert got["the_proper_measure_vanishes"]


def test_the_limit_merges_sizes_not_labels():
    """The reading the measurements support, and the one they replace."""
    got = hs.measure_the_limit_separates_three_things()
    assert len(got["the_three_things"]) == 3
    joined = " ".join(got["the_three_things"]).lower()
    assert "survives" in joined and "f0^n" in joined
    assert "zero proper measure" in got["so_the_origin_is"].lower()
    assert "sizes, not their labels" in got["what_the_first_reading_got_wrong"]
