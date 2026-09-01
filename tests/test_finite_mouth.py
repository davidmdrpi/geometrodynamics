"""The finite-mouth scalar-flat handle in a round `S⁴_R` bulk.

These tests check the module against `docs/finite_mouth_prereg.md`, which was
committed **before** the module existed. That ordering is the point: the repo
audit found 45 probe runs with 45 passes because checks were authored after the
answer was known.

Four kinds. The *forced* ones pin that nothing in the geometry is choosable.
The *junction* ones pin `[h] = [K] = 0` computed from the two sides separately.
The *theorem* ones attack the neck NEC price with hostile lapses instead of
confirming it with a friendly one. And the *regression* ones pin the closed-form
admittance against a solver that shares no algebra with it.
"""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import finite_mouth as fm


# ── P1: nothing is choosable ────────────────────────────────────────────────

def test_the_profile_is_scalar_flat_identically():
    """`⁴R = 0` is what forces `f = √(s²+b²)`, so it must vanish everywhere."""
    for a in (0.05, 0.3, 0.8, 1.2):
        g = fm.geometry(1.0, a)
        s = np.linspace(-g.half_length, g.half_length, 401)
        assert np.max(np.abs(fm.spatial_ricci_scalar(s, g.neck_radius))) < 1e-9


def test_the_matching_constants_are_the_pre_registered_ones():
    """`b = R sin²a`, `S = R sin a cos a`, `L = R sin 2a`."""
    for R in (1.0, 2.5):
        for a in (0.05, 0.3, 0.8, 1.2):
            g = fm.geometry(R, a)
            assert g.neck_radius == pytest.approx(R * math.sin(a) ** 2, rel=1e-14)
            assert g.half_length == pytest.approx(
                R * math.sin(a) * math.cos(a), rel=1e-14)
            assert g.proper_length == pytest.approx(R * math.sin(2 * a), rel=1e-14)
            assert g.mouth_radius == pytest.approx(R * math.sin(a), rel=1e-14)


def test_the_throat_meets_the_exterior_in_radius_and_slope():
    """Darmois: both conditions, checked against the exterior's own values."""
    for R, a in ((1.0, 0.3), (2.5, 0.8)):
        g = fm.geometry(R, a)
        assert float(fm.throat_radius(g.half_length, g.neck_radius)) == \
            pytest.approx(R * math.sin(a), abs=1e-14)
        assert float(fm.throat_radius_derivative(
            g.half_length, g.neck_radius)) == pytest.approx(math.cos(a), abs=1e-14)


def test_no_second_matching_pair_exists():
    """The falsifier for P1: two conditions on two constants has one solution,
    so perturbing `b` while holding the areal radius must break the slope."""
    R, a = 1.0, 0.3
    g = fm.geometry(R, a)
    for fraction in (0.05, -0.05, 0.2):
        b2 = g.neck_radius * (1.0 + fraction)
        s2 = math.sqrt(max(g.mouth_radius ** 2 - b2 ** 2, 0.0))
        assert abs(s2 / g.mouth_radius - math.cos(a)) > 1e-3


def test_the_mouth_angle_is_restricted_to_the_physical_range():
    with pytest.raises(ValueError):
        fm.geometry(1.0, 0.0)
    with pytest.raises(ValueError):
        fm.geometry(1.0, math.pi)
    with pytest.raises(ValueError):
        fm.ThroatGeometry(-1.0, 0.3)


def test_the_rapidity_satisfies_tanh_X_equals_cos_a():
    """`X = arcosh(1/sin a)` is what makes the admittance close in `cos a`."""
    for a in (0.05, 0.3, 0.8, 1.2):
        assert math.tanh(fm.geometry(1.0, a).rapidity) == pytest.approx(
            math.cos(a), abs=1e-14)


# ── P2/P3: the seam ─────────────────────────────────────────────────────────

def test_the_misner_sharp_mass_is_constant_inside_and_matches_outside():
    for R, a in ((1.0, 0.3), (1.0, 0.8), (2.5, 0.3)):
        g = fm.geometry(R, a)
        s = np.linspace(-g.half_length, g.half_length, 201)
        inside = fm.misner_sharp(s, g.neck_radius)
        assert np.ptp(inside) < 1e-12, "mu must be constant in the throat"
        assert float(inside[0]) == pytest.approx(g.neck_radius ** 2, rel=1e-13)
        assert float(inside[-1]) == pytest.approx(
            float(fm.misner_sharp_exterior(a, R)), rel=1e-13)


def test_the_seam_carries_no_israel_layer():
    """`[h] = [K] = 0`, each side computed from its own geometry."""
    for R, a in ((1.0, 0.3), (1.0, 0.8), (2.5, 0.3), (1.0, 0.05)):
        j = fm.junction_jumps(R, a)
        assert j["induced_metric_jump"] < 1e-12
        assert j["extrinsic_curvature_jump"] < 1e-12
        assert j["surface_stress_vanishes"]


def test_the_normal_pressure_agrees_across_the_seam():
    """Gauss–Codazzi requires it when no shell is present: `−3/(8πG₅R²)`."""
    for R, a in ((1.0, 0.3), (1.0, 0.8), (2.5, 0.3)):
        j = fm.junction_jumps(R, a)
        assert j["normal_pressure_jump"] < 1e-12
        assert j["normal_pressure_inside"] == pytest.approx(-3.0 / R ** 2, rel=1e-13)


def test_the_geometry_is_c1_and_not_c2_and_says_so():
    note = fm.junction_jumps()["second_derivative_jumps"]
    assert "C^1 and not C^2" in note
    assert "not a delta function" in note


# ── P4: the theorem, attacked ───────────────────────────────────────────────

def test_the_density_vanishes_for_every_lapse():
    """Gauss–Codazzi on a time-symmetric slice with `⁴R = 0`. This is what
    makes the neck result independent of the lapse."""
    g = fm.geometry()
    s = np.linspace(-g.half_length, g.half_length, 101)
    for lapse in (fm.lapse_ultrastatic, fm.lapse_vacuum,
                  lambda x, b: 1.0 + 0.7 * x + 2.0 * x ** 2):
        assert np.all(fm.stress_tensor(s, g.neck_radius,
                                       lapse=lapse)["density"] == 0.0)


@pytest.mark.parametrize("lapse", [
    lambda s, b: np.ones_like(s),
    lambda s, b: 1.0 + 0.7 * s,                       # asymmetric: N'(0) != 0
    lambda s, b: 1.0 - 2.0 * s ** 2 + 5.0 * s ** 3,
    lambda s, b: np.exp(4.0 * s),
    lambda s, b: 2.0 + np.cos(9.0 * s),
    lambda s, b: 0.05 + 8.0 * s ** 2,                 # nearly null at the neck
])
def test_no_smooth_traversable_lapse_evades_the_neck_price(lapse):
    """P4's falsifier. The lapse enters `p_s` only through `3f′N′/(fN)`, and
    `f′(0) = 0` is what *makes* `s=0` a neck — so no `N′` can help, and no
    reflection symmetry is needed."""
    g = fm.geometry()
    value = float(fm.stress_tensor(np.array([0.0]), g.neck_radius,
                                   lapse=lapse)["radial_nec"][0])
    assert value < 0.0, "smooth + traversable must violate the radial NEC"
    assert value == pytest.approx(fm.null_energy_at_neck(g.neck_radius), rel=1e-9)


def test_the_neck_price_is_minus_three_over_b_squared():
    for a in (0.1, 0.3, 0.8):
        g = fm.geometry(1.0, a)
        assert fm.null_energy_at_neck(g.neck_radius) == pytest.approx(
            -3.0 / (1.0 * math.sin(a) ** 2) ** 2, rel=1e-13)


def test_the_only_escape_is_a_vanishing_lapse_at_the_neck():
    """`N(0) = 0` is the Tangherlini branch — vacuum, and non-traversable."""
    g = fm.geometry()
    assert float(fm.lapse_vacuum(np.array([0.0]), g.neck_radius)[0]) == 0.0
    assert float(fm.lapse_ultrastatic(np.array([0.0]), g.neck_radius)[0]) == 1.0


def test_the_finite_anec_is_the_closed_form_and_diverges_for_small_mouths():
    """Finite, unlike PR #276's two-infinite-end version, and singular as
    `a → 0` — so finite mouths are not merely numerical regulators."""
    for R, a in ((1.0, 0.3), (2.5, 0.8)):
        expected = -3.0 * (1.0 / math.tan(a)
                           + (0.5 * math.pi - a) / math.sin(a) ** 2) / R
        assert fm.null_energy_integral(R, a) == pytest.approx(expected, rel=1e-13)
    scaled = [fm.null_energy_integral(1.0, a) * a ** 2
              for a in (0.03, 0.01, 0.003)]
    for value in scaled:
        assert value == pytest.approx(-1.5 * math.pi, rel=5e-3)


# ── P5: the admittance regression ───────────────────────────────────────────

def test_the_static_admittance_matches_an_independent_solver():
    """The BVP solve never uses the `sinh`/`cosh` reduction."""
    for ell in (0, 1, 2, 3, 5):
        closed = fm.static_admittance(ell)
        numeric = fm.solve_admittance(ell, steps=4000)
        rel = np.max(np.abs(numeric - closed)) / np.max(np.abs(closed))
        assert rel < 1e-5, f"l={ell} disagreement {rel:.2e}"


def test_the_independent_solver_converges_at_second_order():
    """Guards against the shooting basis that was tried first, whose error
    *grew* under refinement from growing-mode contamination."""
    closed = fm.static_admittance(2)
    errors = []
    for steps in (1000, 2000, 4000):
        errors.append(float(np.max(np.abs(
            fm.solve_admittance(2, steps=steps) - closed))))
    assert errors[0] / errors[1] > 3.5
    assert errors[1] / errors[2] > 3.5


def test_the_admittance_is_symmetric_and_its_monopole_row_sums_vanish():
    """A constant is in the kernel at `ℓ = 0`: no static monopole shunt."""
    for ell in (0, 1, 2, 3):
        y = fm.static_admittance(ell)
        assert y[0, 1] == pytest.approx(y[1, 0], rel=1e-14), "reciprocity"
        assert y[0, 0] == pytest.approx(y[1, 1], rel=1e-14), "mouth symmetry"
    assert np.max(np.abs(fm.static_admittance(0).sum(axis=1))) < 1e-12


def test_the_monopole_conductance_has_both_closed_forms():
    """`G = π²R²sin⁴a/cos a = 2π²/I₃` with `I₃ = 2cos a/(R²sin⁴a)`."""
    for R, a in ((1.0, 0.3), (1.0, 0.8), (2.5, 0.3)):
        g_value = fm.monopole_conductance(R, a)
        assert g_value == pytest.approx(
            math.pi ** 2 * R ** 2 * math.sin(a) ** 4 / math.cos(a), rel=1e-14)
        assert fm.static_resistance(R, a) == pytest.approx(
            2.0 * math.cos(a) / (R ** 2 * math.sin(a) ** 4), rel=1e-14)
        assert g_value == pytest.approx(
            2.0 * math.pi ** 2 / fm.static_resistance(R, a), rel=1e-13)
        y0 = fm.static_admittance(0, R, a)
        assert y0[0, 0] == pytest.approx(g_value, rel=1e-12)
        assert y0[0, 1] == pytest.approx(-g_value, rel=1e-12)


def test_the_static_resistance_is_the_profile_integral():
    """`I₃ = ∫ds/f³` computed by quadrature, not from the closed form."""
    for R, a in ((1.0, 0.3), (2.5, 0.8)):
        g = fm.geometry(R, a)
        s = np.linspace(-g.half_length, g.half_length, 200001)
        quad = np.trapezoid(1.0 / fm.throat_radius(s, g.neck_radius) ** 3, s)
        assert quad == pytest.approx(fm.static_resistance(R, a), rel=1e-8)


# ── P6: one spatial geometry, two lapses ────────────────────────────────────

def test_the_vacuum_lapse_reproduces_the_tangherlini_metric_function():
    """`N_vac² = 1 − b²/r²` with `r² = s²+b²`, on the same spatial profile."""
    g = fm.geometry()
    s = np.linspace(-g.half_length, g.half_length, 401)
    r_sq = s * s + g.neck_radius ** 2
    assert np.allclose(fm.lapse_vacuum(s, g.neck_radius) ** 2,
                       1.0 - g.neck_radius ** 2 / r_sq, atol=1e-14)


def test_the_vacuum_branch_is_stress_free_away_from_the_horizon():
    g = fm.geometry()
    s = np.linspace(0.3 * g.half_length, g.half_length, 60)
    stress = fm.stress_tensor(s, g.neck_radius, lapse=fm.lapse_vacuum)
    assert np.max(np.abs(stress["radial_nec"])) < 1e-6


def test_both_branches_share_the_identical_spatial_profile():
    """The physical fork is `N(0)`, not the geometry."""
    g = fm.geometry()
    s = np.linspace(-g.half_length, g.half_length, 101)
    profile = fm.throat_radius(s, g.neck_radius)
    assert np.max(np.abs(fm.spatial_ricci_scalar(s, g.neck_radius))) < 1e-9
    assert np.all(profile > 0.0)
    assert float(profile[len(s) // 2]) == pytest.approx(g.neck_radius, rel=1e-12)


# ── the probe declares its checks ───────────────────────────────────────────

def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import finite_mouth_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)


# ── post-review corrections ─────────────────────────────────────────────────

def test_the_vacuum_branch_does_not_match_the_ultrastatic_exterior():
    """P6 corrected. The shared object is the *spatial* profile, not a single
    global spacetime: the Tangherlini lapse carries `K^t_t = −tan(a)/R ≠ 0` at
    the seam against `0` outside, so it needs its own exterior."""
    for R, a in ((1.0, 0.3), (1.0, 0.8), (2.5, 0.3)):
        j = fm.junction_jumps(R, a)
        assert not j["vacuum_matches_the_ultrastatic_exterior"]
        assert j["vacuum_timelike_curvature_inside"] == pytest.approx(
            -math.tan(a) / R, rel=1e-12)
        assert j["vacuum_lapse_at_seam"] == pytest.approx(math.cos(a), rel=1e-12)
        assert j["vacuum_lapse_jump"] > 1e-3
        # the ultrastatic branch, by contrast, still matches
        assert j["surface_stress_vanishes"]


def test_the_module_records_that_the_fork_is_not_merely_N_at_zero():
    note = fm.junction_jumps()["the_fork_is_not_merely_N_at_zero"]
    assert "SPATIAL profile only" in note
    assert "Israel layer" in note


def test_the_scalar_flat_condition_is_named_as_an_independent_input():
    """"One assumption, everything else forced" hid two more: `⁴R = 0` in the
    handle does not follow from the `S⁴_R` completion, and `a` is free."""
    doc = fm.__doc__
    assert "Three** inputs, not one" in doc or "**Three** inputs" in doc
    assert "does *not* follow from (1)" in doc
    assert "free dimensionless shape parameter" in doc
