"""
Tests for the gluing rule at the seam — the scaling was a free choice.

v46 folded by translation, ``r → r − gap``. That carries a radial *offset*
across unchanged, and the two boundary circles have different circumferences,
so the emerging feature is squashed by ``R_outer/R_inner``. The conformal rule
``r → r·(R_inner/R_outer)`` scales the offset with the boundary it crosses and
returns a faithful copy.

Three things are pinned down here:

* what the choice **changes** — the emerging feature's shape, the sheet ladder
  (and whether it survives the origin at all), and what a winding number would
  look like if there were one;
* what it **does not change** — the winding number itself, which stays
  identically zero on every rule, because a single-valued height has degree
  zero whichever coordinate the seam translates in;
* and that the two rules **agree to first order**, so earlier results were not
  wrong where the wave was small.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.viz.circle_slice import (
    RETURN_TIME,
    BulkAnnulus,
    CircleSlice,
    measure_the_emerging_feature_keeps_its_shape,
    measure_the_rule_does_not_rescue_the_winding,
    measure_the_translate_rule_runs_off_the_origin,
    measure_winding_is_a_magnification,
    winding_curve,
)
from geometrodynamics.viz.warped_sphere import NestedShells


# ── the rules themselves ────────────────────────────────────────────────────
def test_an_unknown_mode_is_rejected():
    with pytest.raises(ValueError):
        BulkAnnulus(mode="whatever")
    with pytest.raises(ValueError):
        CircleSlice(radial_law="whatever")


def test_the_conformal_seam_scales_and_the_translate_seam_shifts():
    sh = NestedShells()
    t, c = BulkAnnulus(sh, "translate"), BulkAnnulus(sh, "conformal")
    assert t.scale_per_crossing == pytest.approx(1.0)
    assert c.scale_per_crossing == pytest.approx(sh.r_inner / sh.r_outer)
    assert t.period == pytest.approx(sh.r_outer - sh.r_inner)
    assert c.period == pytest.approx(math.log(sh.r_outer / sh.r_inner))


def test_both_rules_land_inside_the_annulus():
    sh = NestedShells()
    for mode, r in (("translate", np.linspace(-3.0, 6.0, 4000)),
                    ("conformal", np.geomspace(1e-3, 60.0, 4000))):
        b = BulkAnnulus(sh, mode)
        drawn, _ = b.wrap(r)
        assert np.all(drawn >= b.r_inner - 1e-12)
        assert np.all(drawn < b.r_outer + 1e-12)


def test_the_conformal_seam_refuses_a_nonpositive_radius():
    """It is a translation in ``ln r``, so there is nothing to translate."""
    b = BulkAnnulus(mode="conformal")
    with pytest.raises(ValueError):
        b.wrap(np.array([-0.1]))
    with pytest.raises(ValueError):
        b.wrap(np.array([0.0]))


# ── what the choice changes ─────────────────────────────────────────────────
def test_the_emerging_feature_is_squashed_or_faithful():
    r = measure_the_emerging_feature_keeps_its_shape()
    assert r["translate_squashes_the_feature"] is True
    assert r["conformal_preserves_the_shape"] is True
    assert r["translate_distortion"] == pytest.approx(r["circumference_ratio"],
                                                      rel=1e-9)
    assert r["conformal_distortion"] == pytest.approx(1.0, abs=1e-9)


def test_the_translate_sheets_walk_through_the_origin():
    r = measure_the_translate_rule_runs_off_the_origin()
    assert r["translate_reaches_negative_radius"] is True
    assert r["conformal_stays_positive"] is True
    assert min(r["translate"]["inward"]) < 0.0
    assert min(r["conformal"]["inward"]) > 0.0
    # geometric sheets shrink by a constant factor, so they never arrive
    inner = r["conformal"]["inward"]
    ratios = [b / a for a, b in zip(inner, inner[1:])]
    assert max(ratios) - min(ratios) < 1e-9


def test_the_multiplicative_radial_law_is_always_positive():
    s = CircleSlice(radial_law="multiplicative", n_sigma=361)
    s.reset()
    s.advance_to(2.98)
    for gain in (0.5, 5.0, 50.0):
        assert np.all(s.radius(gain=gain) > 0.0)
    add = CircleSlice(radial_law="additive", n_sigma=361)
    add.reset()
    add.advance_to(2.98)
    assert np.any(add.radius(gain=50.0) <= 0.0)          # ...and this one is not


def test_the_two_radial_laws_agree_to_first_order():
    mul = CircleSlice(radial_law="multiplicative", n_sigma=361)
    add = CircleSlice(radial_law="additive", n_sigma=361)
    for s in (mul, add):
        s.reset()
        s.advance_to(1.4)
    prev = None
    for gain in (1e-2, 1e-3):
        diff = float(np.max(np.abs(mul.radius(gain=gain)
                                   - add.radius(gain=gain))))
        if prev is not None:
            assert diff < 0.02 * prev                    # second order in ε
        prev = diff


# ── what winding would look like ────────────────────────────────────────────
def test_a_wound_curve_climbs_exactly_one_sheet():
    for mode in ("translate", "conformal"):
        b = BulkAnnulus(mode=mode)
        _, r = winding_curve(b, turns=1, r_start=b.r_inner * 1.01)
        _, sheet = b.wrap(r)
        assert int(sheet[-1] - sheet[0]) == 1


def test_conformal_winding_is_a_start_independent_magnification():
    r = measure_winding_is_a_magnification()
    assert r["conformal_turns_winding_into_a_scale"] is True
    assert r["translate_gives_no_scale_at_all"] is True
    assert r["conformal_spread"] < 1e-9
    assert r["translate_spread"] > 0.1
    sh = NestedShells()
    assert r["conformal_magnification"] == pytest.approx(sh.r_outer / sh.r_inner,
                                                         rel=1e-9)


def test_two_turns_magnify_by_the_square():
    r = measure_winding_is_a_magnification(turns=2)
    sh = NestedShells()
    assert r["conformal_magnification"] == pytest.approx(
        (sh.r_outer / sh.r_inner) ** 2, rel=1e-9)
    assert r["a_wound_curve_climbs_the_sheets"] is True


def test_the_wound_curve_is_not_a_graph():
    """Which is the whole reason the wave cannot be it.

    Its unwrapped radius at ``σ = +π`` differs from ``σ = −π``, so no
    single-valued ``f(σ)`` on the circle can produce it.
    """
    b = BulkAnnulus(mode="conformal")
    _, r = winding_curve(b, turns=1, r_start=1.0)
    assert abs(float(r[-1] - r[0])) > 0.5
    assert float(r[-1] / r[0]) == pytest.approx(b.r_outer / b.r_inner, rel=1e-9)


# ── what the choice does not change ─────────────────────────────────────────
def test_no_gluing_rescues_the_winding_number():
    """The check the scaling objection deserves — and the answer is still zero."""
    r = measure_the_rule_does_not_rescue_the_winding(
        gains=(1.0, 6.0, 40.0), frames=60)
    assert r["every_rule_gives_zero"] is True
    assert r["the_negative_result_survives_the_choice"] is True
    assert r["worst_signed_or_winding"] == 0
    assert max(row["unsigned"] for row in r["rows"]) > 10   # genuinely driven
    for row in r["rows"]:
        assert row["signed"] == 0
        assert row["winding"] == 0


def test_the_conformal_slice_still_crosses_in_pairs():
    s = CircleSlice(bulk=BulkAnnulus(mode="conformal"),
                    radial_law="multiplicative", n_sigma=721)
    s.reset()
    seen = 0
    for i in range(40):
        s.advance_to((i + 1) * RETURN_TIME / 40)
        for gain in (0.6, 2.0):
            c = s.seam_crossings(gain=gain)
            assert c["signed"] == 0
            assert c["unsigned"] % 2 == 0
            assert c["signed"] == s.winding_number(gain=gain)
            seen = max(seen, c["unsigned"])
    assert seen >= 2


def test_the_unrolled_strip_uses_the_seams_own_coordinate():
    for mode, law in (("translate", "additive"), ("conformal", "multiplicative")):
        s = CircleSlice(bulk=BulkAnnulus(mode=mode), radial_law=law,
                        n_sigma=361)
        s.reset()
        s.advance_to(1.7)
        sigma, rho = s.unrolled(gain=1.0)
        assert np.allclose(sigma, s.sigma)
        assert np.all((rho >= 0.0) & (rho < 1.0 + 1e-12))
        drawn, _ = s.bulk.wrap(s.radius(gain=1.0))
        assert np.allclose(np.linalg.norm(s.points(gain=1.0), axis=1), drawn)
