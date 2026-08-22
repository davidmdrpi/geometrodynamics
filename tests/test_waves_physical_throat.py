"""Which throat is physical, and the sign it gives.

The geometry checks here are against *exact* cases — the round sphere, the
closed-form length and resistance, the identity ``(f²u')' = 0`` — because the
whole point of this round is that the throat stopped being a free parameter.
A parameter you have fixed by derivation has to be checked against something
other than itself.
"""

import math

import numpy as np
import pytest

from geometrodynamics.waves import areal, physical_throat as pt


# ── the curvature, from the metric ──────────────────────────────────────────
def test_the_round_sphere_has_scalar_curvature_six():
    for s in (0.2, 0.7, 1.4, 2.5, 3.0):
        got = pt.profile_scalar_curvature(math.sin(s), math.cos(s), -math.sin(s))
        assert abs(got - 6.0) < 1e-12


def test_the_vacuum_profile_has_zero_scalar_curvature():
    """``f'² = 1 − f₀/f`` implies ``f'' = f₀/2f²``, and then ``R̂ = 0`` exactly."""
    for f0 in (1e-4, 1e-2, 0.3):
        for f in (1.5 * f0, 4.0 * f0, 0.5, 2.0):
            if f <= f0:
                continue
            got = pt.profile_scalar_curvature(
                f, math.sqrt(1.0 - f0 / f), f0 / (2.0 * f ** 2))
            assert abs(got) < 1e-8 * max(1.0, 1.0 / f ** 2)


def test_the_curvature_measurement_reports_both_exact_cases():
    got = pt.measure_the_curvature_formula_is_derived_not_quoted()
    assert got["both_are_exact"]
    assert got["worst_sphere_error"] < 1e-12


# ── the gluing forces everything ────────────────────────────────────────────
@pytest.mark.parametrize("a", [0.02, 0.05, 0.10, 0.20])
def test_smooth_gluing_forces_the_neck_radius(a):
    """``cos²a = 1 − f₀/sin a`` has exactly one root, and it is ``sin³a``."""
    t = pt.VacuumThroat(mouth_radius=a)
    f0 = t.neck_radius()
    assert abs(f0 - math.sin(a) ** 3) < 1e-15
    assert abs((1.0 - f0 / math.sin(a)) - math.cos(a) ** 2) < 1e-15
    assert 0.0 < f0 < math.sin(a)


@pytest.mark.parametrize("a", [0.02, 0.05, 0.10, 0.20])
def test_the_length_and_resistance_closed_forms_are_the_integrals(a):
    t = pt.VacuumThroat(mouth_radius=a)
    assert abs(t.length() / t.length_by_quadrature() - 1.0) < 1e-9
    assert abs(t.resistance() / t.resistance_by_quadrature() - 1.0) < 1e-9


@pytest.mark.parametrize("a", [0.02, 0.05, 0.10, 0.20])
def test_the_conductance_is_exactly_a_quarter_of_the_exteriors(a):
    """``4π/I = π sin³a/cos a`` against ``N₀ = 4π sin²a tan a``.

    An exact ratio with no free parameter in it, at every mouth radius.
    """
    t = pt.VacuumThroat(mouth_radius=a)
    n0 = 4.0 * math.pi * math.sin(a) ** 2 * math.tan(a)
    assert abs(t.conductance() / n0 - 0.25) < 1e-12


def test_the_throat_is_short_and_its_length_is_not_a_parameter():
    got = pt.measure_the_gluing_forces_the_neck_radius()
    assert got["length_is_about_twice_the_mouth_radius"]
    assert got["the_conductance_is_a_quarter_of_the_exteriors"]
    assert got["worst_length_drift"] < 1e-9
    assert got["worst_resistance_drift"] < 1e-9


# ── what the other throats would need ───────────────────────────────────────
def test_a_product_tube_of_area_four_pi_needs_a_third_the_ambient_density():
    assert abs(pt.product_tube_density_ratio(4 * math.pi) - 1.0 / 3.0) < 1e-12


def test_the_matched_tube_needs_matter_a_hundred_times_denser():
    ratio = pt.product_tube_density_ratio(4 * math.pi * math.sin(0.05) ** 2)
    assert ratio > 100.0


def test_only_one_product_area_carries_the_ambient_fluid():
    assert abs(pt.product_tube_density_ratio(4 * math.pi / 3.0) - 1.0) < 1e-12
    got = pt.measure_the_product_tubes_need_anomalous_matter()
    assert got["neither_used_area_is_the_ambient_fluid"]
    assert got["vacuum_throat_density_ratio"] == 0.0


# ── no cavity ───────────────────────────────────────────────────────────────
def test_the_vacuum_throats_monopole_admittance_is_the_closed_form():
    """Deliberately the *Riccati* solve, not `admittance`.

    Since the closed form became production, `admittance(0)` and
    `monopole_admittance_closed_form` are the same expression, and comparing
    them would check nothing.  The solve is fully resolved at ``ℓ = 0``, so it
    is the independent number here.
    """
    for a in (0.05, 0.10):
        t = pt.VacuumThroat(mouth_radius=a)
        num = t.admittance_riccati(0)
        closed = t.monopole_admittance_closed_form()
        assert np.max(np.abs(num - closed)) < 1e-10 * np.max(np.abs(closed))


def test_the_monopole_admittance_is_singular_because_flux_is_conserved():
    """``(f²u')' = 0`` makes ``Y·(1,1) = 0`` an identity, not a cancellation.

    With the closed form as production the *matrix* identities are now exact
    to the last bit — which is the reason ``ℓ = 0`` is special-cased as
    ``G·[[−1,1],[1,−1]]`` rather than evaluated through the general
    ``k coth(2kX) − tanh X``.

    The determinant is the exception, and deliberately not asserted exact:
    `np.linalg.det` is an LU factorisation with pivoting and carries its own
    rounding, so it returns ``1.4e-21`` on a matrix whose rows are bitwise
    negatives of each other.  That is the routine's error, not the throat's,
    and the relative tolerance is where it belongs.
    """
    for a in (0.02, 0.05, 0.10, 0.20):
        t = pt.VacuumThroat(mouth_radius=a)
        y = t.admittance(0)
        scale = float(np.max(np.abs(y)))
        assert np.max(np.abs(y @ np.ones(2))) == 0.0
        assert t.shunt() == 0.0
        assert y[0, 0] == -y[0, 1] and y[1, 1] == -y[1, 0]
        assert abs(float(np.linalg.det(y))) < 1e-14 * scale ** 2


def test_a_tube_with_matter_in_it_shunts_and_a_vacuum_one_does_not():
    assert abs(pt.VacuumThroat().shunt()) < 1e-12
    assert abs(areal.TubeModel().shunt()) > 1.0


def test_the_vacuum_throat_has_no_resonances_to_flip_a_sign_across():
    got = pt.measure_the_vacuum_throat_has_no_cavity()
    assert got["vacuum_throat_has_resonances"] is False
    assert got["it_matches_the_closed_form"]
    assert got["the_symmetric_channel_is_exactly_dead"]
    assert got["tube_scalar_curvature"] == 0.0
    assert len(got["wide_tube_has_resonances_at"]) == 3


def test_the_admittance_is_symmetric_in_both_channels():
    t = pt.VacuumThroat(mouth_radius=0.05)
    for ell in (0, 1):
        y = t.admittance(ell)
        assert abs(y[0, 1] - y[1, 0]) < 1e-14 * max(1.0, abs(y[0, 0]))


def test_the_dipole_channel_is_evanescent_so_the_mouths_decouple():
    """The neck attenuates ``ℓ = 1`` by ``1e+09``; the through-coupling is tiny."""
    t = pt.VacuumThroat(mouth_radius=0.05)
    y = t.admittance(1)
    assert abs(y[0, 1]) < 1e-6 * abs(y[0, 0])
    n1 = -2.0 * math.pi * math.sin(2 * 0.05)
    assert abs(y[0, 0] / n1 - 1.0) < 0.01, "and it nearly matches the exterior"


# ── the mechanism ───────────────────────────────────────────────────────────
def test_the_shunt_decides_the_sign_and_the_conductance_does_not():
    got = pt.measure_the_shunt_decides_the_sign()
    assert got["conductance_never_changes_the_sign"]
    assert got["the_shunt_does"]
    assert got["the_shunt_is_the_tubes_matter"]
    assert got["vacuum_shunt"] == 0.0
    assert got["wide_shunt"] > 1.0
    corners = {c["label"]: c["sign"] for c in got["corners"]}
    assert corners["vacuum throat"] == ["opens", "opens"]
    assert corners["wide conductance, no shunt"] == ["opens", "opens"]
    assert corners["vacuum conductance, wide shunt"] == ["closes", "closes"]
    assert corners["wide conductance and wide shunt"] == ["closes", "closes"]


def test_the_wide_tube_itself_reproduces_the_previous_round():
    """The comparison must contain PR #264 exactly, or it is not one.

    Note this is a *separate* entry from the 2x2 corners: those hold the
    l=1 admittance fixed at the vacuum throat's so that only the l=0
    decomposition varies, which is the controlled experiment. Swapping the
    l=1 channel too is what reproduces #264 to machine precision.
    """
    got = pt.measure_the_shunt_decides_the_sign()
    corner = got["the_actual_wide_tube"]
    m = areal.INTERFERENCE_MOMENTS[1]
    ref = np.asarray(areal.solve_matching(
        areal.MOUTHS, m.radius, areal.WORKING_TUBE, m.as_source(),
        m.signed_obstruction())["areal_response"], float)
    assert np.max(np.abs(np.array(corner) / ref - 1.0)) < 1e-9


# ── the answer ──────────────────────────────────────────────────────────────
def test_the_physical_throat_opens_both_mouths():
    got = pt.measure_the_signed_response_on_the_physical_throat()
    assert got["every_variant_agrees_in_sign"]
    assert got["sign"] == ["opens", "opens"]
    assert all(v > 0 for v in got["areal_response"])
    assert got["worst_residual"] < 1e-10
    assert got["quadrature_spread_at_fixed_radius"] < 0.05


def test_it_reverses_the_wide_tube_rather_than_agreeing_with_it():
    physical = pt.measure_the_signed_response_on_the_physical_throat()
    wide = areal.measure_the_signed_areal_response()
    assert physical["sign"] == ["opens", "opens"]
    assert wide["sign"] == ["closes", "closes"]
    assert (np.sign(physical["areal_response"][0])
            != np.sign(wide["areal_response"][0]))


def test_the_response_is_reported_with_its_cost_not_without_it():
    """The magnitude is large and grows as a^-3; that must be visible."""
    got = pt.measure_the_signed_response_on_the_physical_throat()
    assert "a^-3" in got["what_it_costs"]
    coarse = [r for r in got["rows"] if r["radius"] == 0.10
              and r["gluing"] == "transported"][-1]
    fine = [r for r in got["rows"] if r["radius"] == 0.05
            and r["gluing"] == "transported"][-1]
    grew = fine["areal_response"][0] / coarse["areal_response"][0]
    assert grew > 8.0, "halving the mouth radius must magnify it steeply"
    assert got["worst_condition_number"] > 1e6


def test_the_sign_survives_the_whole_vacuum_family_not_just_the_glued_point():
    """The gluing condition removes the last freedom, so the answer must not
    depend on hitting it exactly.  Four orders in the neck radius."""
    m = areal.INTERFERENCE_MOMENTS[1]
    a = m.radius
    basis = areal.basis_channels(areal.MOUTHS, a)
    glued = math.sin(a) ** 3

    class _Free(pt.VacuumThroat):
        pass

    signs = []
    for frac in (0.02, 0.2, 1.0, 5.0, 50.0, 300.0):
        f0 = glued * frac
        if f0 >= math.sin(a) * 0.99:
            continue
        throat = pt.VacuumThroat(mouth_radius=a)
        object.__setattr__(throat, "mouth_radius", a)
        y0 = _family_admittance(a, f0, 0)
        y1 = _family_admittance(a, f0, 1)

        class _Fixed:
            def admittance(self, ell, _y0=y0, _y1=y1):
                return _y0 if int(ell) == 0 else _y1

        v = np.asarray(areal.solve_matching(
            areal.MOUTHS, a, _Fixed(), m.as_source(), m.signed_obstruction(),
            basis=basis)["areal_response"], float)
        signs.append(np.sign(v))
    assert len(signs) >= 5
    assert all(np.array_equal(s, signs[0]) for s in signs)
    assert signs[0][0] > 0


def _family_admittance(a, f0, ell):
    """The vacuum family with the neck radius left free (gluing not imposed)."""
    from scipy.integrate import solve_ivp
    m, c = math.sin(a), float(ell * (ell + 1))
    tmax = math.sqrt(m - f0)

    def run(parity):
        if parity == "sym":
            rhs = lambda t, y: [2 * math.sqrt(f0 + t * t)
                                * (c - y[0] ** 2 / (f0 + t * t) ** 2)]
            s = solve_ivp(rhs, (0, tmax), [0.0], rtol=1e-12, atol=1e-16,
                          method="DOP853")
            return -4 * math.pi * float(s.y[0, -1])
        rhs = lambda t, y: [2 * math.sqrt(f0 + t * t)
                            * (1.0 / (f0 + t * t) ** 2 - c * y[0] ** 2)]
        s = solve_ivp(rhs, (0, tmax), [0.0], rtol=1e-12, atol=1e-18,
                      method="DOP853")
        return -4 * math.pi / float(s.y[0, -1])

    p, q = run("sym"), run("anti")
    return np.array([[0.5 * (p + q), 0.5 * (p - q)],
                     [0.5 * (p - q), 0.5 * (p + q)]])


# ── what the forced throat turns out to be ──────────────────────────────────
@pytest.mark.parametrize("a", [0.02, 0.05, 0.10, 0.20])
def test_the_profile_is_the_time_symmetric_schwarzschild_slice(a):
    """`ds² = dr²/(1−2M/r) + r²dΩ²` in proper distance is `f'² = 1 − 2M/f`."""
    t = pt.VacuumThroat(mouth_radius=a)
    m = t.mass()
    for frac in (1.2, 3.0, 50.0):
        f = min(t.neck_radius() * frac, t.mouth_f())
        schwarzschild = 1.0 - 2.0 * m / f
        forced = 1.0 - t.neck_radius() / f
        assert abs(schwarzschild - forced) < 1e-15


@pytest.mark.parametrize("a", [0.02, 0.05, 0.10, 0.20])
def test_the_neck_radius_is_twice_the_mass(a):
    t = pt.VacuumThroat(mouth_radius=a)
    assert abs(2.0 * t.mass() - t.neck_radius()) < 1e-18
    assert abs(t.mass() - math.sin(a) ** 3 / 2.0) < 1e-18


@pytest.mark.parametrize("a", [0.02, 0.05, 0.10, 0.20])
def test_three_independent_masses_agree(a):
    """Schwarzschild parameter, irreducible mass, Hawking mass."""
    t = pt.VacuumThroat(mouth_radius=a)
    ms = [t.mass(), t.irreducible_mass(), t.hawking_mass_in_the_tube()]
    assert (max(ms) - min(ms)) < 1e-12 * ms[0]


@pytest.mark.parametrize("a", [0.02, 0.05, 0.10, 0.20])
def test_the_neck_area_is_sixteen_pi_m_squared(a):
    t = pt.VacuumThroat(mouth_radius=a)
    assert abs(t.neck_area() - 16.0 * math.pi * t.mass() ** 2) < 1e-24


def test_the_hawking_mass_is_constant_along_the_vacuum_throat():
    """`(f/2)(1 − f'²) = f₀/2` at every areal radius, not just at the ends."""
    t = pt.VacuumThroat(mouth_radius=0.05)
    f0 = t.neck_radius()
    for frac in (1.01, 1.5, 8.0, 100.0, 390.0):
        f = f0 * frac
        if f > t.mouth_f():
            continue
        got = pt.hawking_mass(f, math.sqrt(1.0 - f0 / f))
        # `1 - f'^2` recovers `f0/f` by cancellation between numbers near 1,
        # so the achievable relative accuracy is ~eps*f/f0 -- 4e-14 at the
        # mouth. The exact statement `m_H = f0/2` needs no cancellation; this
        # check deliberately goes through the definition instead.
        assert abs(got - t.mass()) < 1e-12 * t.mass()


@pytest.mark.parametrize("a", [0.02, 0.05, 0.10, 0.20])
def test_the_gluing_condition_is_hawking_mass_continuity(a):
    """The ambient's own Hawking mass at radius `a` is `sin³a/2` — the same
    number the gluing forces the neck to carry."""
    t = pt.VacuumThroat(mouth_radius=a)
    ambient = pt.hawking_mass(math.sin(a), math.cos(a))
    assert abs(ambient - t.hawking_mass_in_the_tube()) < 1e-12 * ambient
    assert abs(ambient - t.mass()) < 1e-15


def test_the_hawking_mass_formula_against_its_definition():
    """`m_H = √(A/16π)(1 − (1/16π)∮H²dA)` with `H = 2f'/f`, `A = 4πf²`."""
    for f, fp in ((0.3, 0.7), (1.0, 0.2), (0.05, 0.99)):
        area = 4.0 * math.pi * f ** 2
        integral = (2.0 * fp / f) ** 2 * area
        expected = (math.sqrt(area / (16.0 * math.pi))
                    * (1.0 - integral / (16.0 * math.pi)))
        assert abs(pt.hawking_mass(f, fp) - expected) < 1e-14


def test_the_mass_law_inverts_and_has_the_expected_small_mouth_limit():
    got = pt.measure_the_throat_is_an_einstein_rosen_bridge()
    for row in got["rows"]:
        a = row["mouth_radius"]
        assert abs(row["mouth_from_the_mass"] - a) < 1e-12
        assert abs(row["mass"] / row["small_mouth_law"] - 1.0) < 0.05


def test_the_bridge_measurement_states_what_it_is_not():
    """The claim is strong, so the things it does not say are asserted."""
    got = pt.measure_the_throat_is_an_einstein_rosen_bridge()
    assert got["it_is_an_einstein_rosen_bridge"]
    assert got["the_mass_has_no_free_parameter"]
    assert got["the_neck_is_a_marginal_sphere"]
    assert got["schwarzschild_slope_error"] < 1e-15
    assert got["three_masses_agree"] < 1e-12
    assert got["the_gluing_is_hawking_mass_continuity"] < 1e-12
    caveats = " ".join(got["what_it_is_not"]).lower()
    for word in ("adm", "dimensionless", "handle", "apparent horizon",
                 "traversab"):
        assert word in caveats


def test_the_marginal_sphere_claim_is_the_identity_and_not_the_horizon():
    """`H = 0` with `K_ij = 0` gives `θ_± = 0`.  That is a statement about one
    surface in one slice.

    *An apparent horizon is the **outermost** MOTS*, which needs a global
    condition on the slice that nothing here evaluates; and *non-traversability
    is a property of the Lorentzian development*, which needs a lapse that
    nothing here chooses.  The measurement has to disclaim both itself rather
    than leave it to a reader — the earlier prose asserted the stronger
    conclusion directly, which is what this pins.
    """
    got = pt.measure_the_throat_is_an_einstein_rosen_bridge()
    identity = got["the_marginal_sphere_identity"].lower()
    assert "h = 0" in identity and "k_ij = 0" in identity
    assert "theta_+ = theta_- = 0" in identity
    assert "mots" in identity

    caveats = " ".join(got["what_it_is_not"]).lower()
    assert "not an apparent horizon" in caveats
    assert "outermost" in caveats
    assert "not shown to be non-traversable" in caveats
    assert "no lapse" in caveats
    # and the one condition under which the stronger claim does follow
    assert "schwarzschild" in caveats and "added assumption" in caveats
    # the key that used to assert the horizon outright is gone, not renamed
    assert "the_neck_is_an_apparent_horizon" not in got


def test_the_marginal_sphere_predicate_keeps_its_old_name_working():
    t = pt.VacuumThroat(mouth_radius=0.05)
    assert t.neck_is_a_marginal_sphere()
    assert t.neck_is_a_minimal_surface(), "the old name still resolves"
    assert abs(t.mouth_f() ** 3 - t.neck_radius()) < 1e-18


# ════════════════════════════════════════════════════════════════════════════
# THE CLOSED-FORM TWO-PORT  (PR #267)
#
# Three kinds of check, deliberately separated:
#
#   1. agreement with the Riccati solve *in the regime the solve resolves* --
#      that is what says the closed form is the same object, not a new one;
#   2. exact structural identities, which hold to the last bit and do not
#      depend on any solver at all;
#   3. the resolution boundary itself -- pinned as a boundary, not as an
#      error factor, because the factor is floating-point dependent and the
#      boundary is not.
# ════════════════════════════════════════════════════════════════════════════

# ── 1. the two routes agree where the solve is resolved ─────────────────────
@pytest.mark.parametrize("a", [0.05, 0.10, 0.20, 0.30])
@pytest.mark.parametrize("ell", [0, 1])
def test_the_closed_form_reproduces_the_riccati_solve_where_it_resolves(a, ell):
    """``ℓ = 0`` and ``ℓ = 1`` are above the solver's floor at every radius here.

    The diagonal is the sharp comparison -- it is ``O(1)`` and carries no
    cancellation, so it must agree to ``1e-10`` relative.  The cross term is
    checked more loosely because the *solve* loses digits to the ``½(s − t)``
    subtraction long before the closed form does; ``ℓ = 1`` at ``a = 0.05`` is
    already down to four figures.
    """
    t = pt.VacuumThroat(mouth_radius=a)
    cf, ric = t.admittance_closed_form(ell), t.admittance_riccati(ell)
    assert abs(cf[0, 0] / ric[0, 0] - 1.0) < 1e-10
    assert abs(cf[0, 1] / ric[0, 1] - 1.0) < 1e-3
    assert cf[0, 1] * ric[0, 1] > 0.0, "and they agree in sign"


def test_the_monopole_channel_agrees_to_ten_digits_by_both_routes():
    """The one channel where nothing cancels in either route."""
    for a in (0.02, 0.05, 0.10, 0.20, 0.50):
        t = pt.VacuumThroat(mouth_radius=a)
        cf, ric = t.admittance_closed_form(0), t.admittance_riccati(0)
        assert np.max(np.abs(cf / ric - 1.0)) < 1e-10


def test_production_admittance_is_the_closed_form_not_the_solve():
    for a in (0.05, 0.20):
        t = pt.VacuumThroat(mouth_radius=a)
        for ell in (0, 1, 2):
            assert np.array_equal(t.admittance(ell),
                                  t.admittance_closed_form(ell))


def test_the_riccati_solve_is_retained_as_a_validator():
    """It is demoted, not deleted -- the check above needs something to check
    *against*, and an implementation with no second opinion is not verified."""
    assert hasattr(pt.VacuumThroat, "admittance_riccati")
    assert hasattr(pt.VacuumThroat, "admittance_closed_form")
    y = pt.VacuumThroat(mouth_radius=0.05).admittance_riccati(0)
    assert np.all(np.isfinite(y))


# ── 2. exact structural identities ──────────────────────────────────────────
@pytest.mark.parametrize("a", [0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.0])
def test_the_zero_shunt_identity_is_exact_in_the_closed_form(a):
    """``Y₀·(1,1)ᵀ = 0`` to the last bit, because ``ℓ = 0`` is written as
    ``G·[[−1,1],[1,−1]]`` rather than assembled from ``D₀`` and ``C₀``."""
    y = pt.VacuumThroat(mouth_radius=a).admittance_closed_form(0)
    assert np.array_equal(y @ np.ones(2), np.zeros(2))
    assert y[0, 0] == -y[0, 1]
    assert y[0, 1] == y[1, 0]


@pytest.mark.parametrize("a", [0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.0])
def test_the_monopole_cross_term_is_four_pi_over_the_resistance(a):
    """``C₀ = 4π/I`` with ``I = ∫ds/f² = 4cos a/sin³a`` -- the DC conductance
    of the tube, arrived at from the wave equation rather than assumed."""
    t = pt.VacuumThroat(mouth_radius=a)
    assert t.admittance_closed_form(0)[0, 1] == 4.0 * math.pi / t.resistance()
    assert abs(t.admittance_closed_form(0)[0, 1]
               - math.pi * math.sin(a) ** 3 / math.cos(a)) < 1e-18


@pytest.mark.parametrize("a", [0.005, 0.01, 0.02, 0.05, 0.10, 0.50])
def test_the_general_formula_reproduces_the_monopole_special_case(a):
    """The special case is a *stability* choice, not a different answer.

    The general ``q``-form's cross term is exact — it is the same product.
    Its *diagonal* is not: ``k coth(2kX) − tanh X`` is a difference of two
    quantities that both tend to ``1`` as ``a → 0``, so the answer is built
    out of the cancelling tail.  It still agrees to ``1e-11``, which is why
    this is a stability choice rather than a correction, but the error grows
    as the mouth closes — ``3e-14`` at ``a = 0.02``, ``6e-12`` at
    ``a = 0.005`` — and that is the near-cancellation the special case exists
    to avoid.  The next test shows what it buys.
    """
    k, q = 1.0, math.tan(0.5 * a) ** 2
    w = 1.0 - q * q
    d = -2.0 * math.pi * math.sin(a) * (k * (1.0 + q * q) / w - math.cos(a))
    c = 4.0 * math.pi * math.sin(a) * k * q / w
    y = pt.VacuumThroat(mouth_radius=a).admittance_closed_form(0)
    assert abs(c / y[0, 1] - 1.0) < 1e-15
    assert abs(d / y[0, 0] - 1.0) < 1e-11


def test_the_special_case_is_what_makes_the_zero_shunt_identity_exact():
    """Side by side: the general form's shunt residue against the special
    case's, which is exactly zero at every radius."""
    worst_general = 0.0
    for a in (0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50):
        q = math.tan(0.5 * a) ** 2
        w = 1.0 - q * q
        d = -2.0 * math.pi * math.sin(a) * ((1.0 + q * q) / w - math.cos(a))
        c = 4.0 * math.pi * math.sin(a) * q / w
        worst_general = max(worst_general, abs(d + c) / abs(c))
        assert pt.VacuumThroat(mouth_radius=a).shunt() == 0.0
    assert worst_general > 1e-13, "the general form really does leave a residue"
    assert worst_general < 1e-10, "and it is small -- a stability choice, not a fix"


@pytest.mark.parametrize("ell", [1, 2, 3, 4])
@pytest.mark.parametrize("a", [0.01, 0.02, 0.05])
def test_the_small_mouth_asymptotic_is_the_leading_term(ell, a):
    """``C_ℓ → 4π(2ℓ+1) sin a · tan^{4ℓ+2}(a/2)``.

    The correction is ``O(q²) = O(a^{8ℓ+4})``, so at these radii the leading
    term *is* the answer to machine precision.
    """
    exact = pt.VacuumThroat(mouth_radius=a).admittance_closed_form(ell)[0, 1]
    lead = (4.0 * math.pi * (2 * ell + 1) * math.sin(a)
            * math.tan(0.5 * a) ** (4 * ell + 2))
    assert abs(exact / lead - 1.0) < 1e-12


@pytest.mark.parametrize("a", [0.05, 0.20, 0.50])
@pytest.mark.parametrize("ell", [0, 1, 2, 3])
def test_the_closed_form_two_port_is_symmetric_and_the_cross_term_positive(a, ell):
    y = pt.VacuumThroat(mouth_radius=a).admittance_closed_form(ell)
    assert y[0, 1] == y[1, 0]
    assert y[0, 0] == y[1, 1]
    assert y[0, 1] > 0.0 and y[0, 0] < 0.0


def test_the_exponent_hierarchy_is_four_powers_per_unit_of_angular_momentum():
    got = pt.measure_the_mouth_to_mouth_hierarchy()
    assert got["the_exponent_law_holds"]
    assert got["the_asymptotic_is_the_leading_term"]
    assert [f["expected"] for f in got["fits"]] == [3, 7, 11, 15]
    assert 8.0e5 < got["monopole_beats_dipole_by"] < 9.0e5
    row = got["rows"][0]
    assert row["mouth_radius"] == 0.05
    assert 5e-10 < row["ell_one_transmission"] < 1e-9


def test_the_hierarchy_states_the_narrow_claim_and_not_the_wide_one():
    """The proved statement is about one operator on one slice.  It is *not*
    that BAM cannot transmit orientation through the throat, and the
    measurement has to say so itself rather than leaving it to a reader."""
    got = pt.measure_the_mouth_to_mouth_hierarchy()
    narrow = got["the_narrow_statement"].lower()
    assert "static scalar laplacian" in narrow
    assert "l=1" in narrow and "monopole" in narrow
    assert "a = 0.05" in narrow
    disclaimer = got["what_is_not_claimed"].lower()
    assert "orientation" in disclaimer
    assert "lapse" in disclaimer
    assert "small rather than zero" in disclaimer


# ── 3. the resolution boundary ──────────────────────────────────────────────
def test_the_closed_form_still_produces_a_dipole_scale_answer_at_ell_two():
    """``C₂ ≈ 3.00105e-16`` at ``a = 0.05`` -- and, crucially, **not** built
    as a difference of the two eigenchannels.

    ``3e-16`` is roughly ``eps`` times the diagonal ``1.2565``, which is
    exactly why no subtraction of two ``O(1)`` numbers can produce it.  The
    closed form reaches it as a product, ``4π sin a · 5q/(1−q²)``, and every
    factor there is representable.
    """
    t = pt.VacuumThroat(mouth_radius=0.05)
    c2 = float(t.admittance_closed_form(2)[0, 1])
    assert abs(c2 / 3.001054624604e-16 - 1.0) < 1e-9
    # a product, not a difference: the same number from the raw factors
    q = math.tan(0.025) ** 10
    by_hand = 4.0 * math.pi * math.sin(0.05) * 5.0 * q / (1.0 - q * q)
    assert c2 == by_hand
    # and it is far below what a subtraction on this diagonal could carry
    assert c2 < 1e-15 * abs(float(t.admittance_closed_form(2)[0, 0]))


def test_the_riccati_solve_stops_resolving_the_cross_term_at_ell_two():
    got = pt.measure_where_the_riccati_solve_stops_resolving()
    assert got["riccati_is_trustworthy_through_ell"] == 1
    assert got["the_cross_term_fails_at_ell_two"]
    assert got["the_diagonal_was_never_affected"]
    rows = {r["ell"]: r for r in got["rows"]}
    assert rows[0]["cross_relative_error"] < 1e-12
    assert rows[1]["cross_relative_error"] < 1e-3
    assert rows[2]["cross_relative_error"] > 1.0
    assert not rows[2]["signs_agree"]


def test_the_boundary_is_pinned_and_the_error_factor_is_not():
    """A test asserting *39x* would be a test of one build of SciPy.

    What is asserted is the pair of statements that survive a change of
    solver: the sign is wrong, and the magnitude is above the true value by
    more than an order of magnitude, at a channel whose honest size is below
    the solver's floor.
    """
    got = pt.measure_where_the_riccati_solve_stops_resolving()
    row = {r["ell"]: r for r in got["rows"]}[2]
    assert row["cross_closed_form"] > 0.0
    assert row["cross_riccati"] < 0.0
    assert abs(row["cross_riccati"]) > 10.0 * row["cross_closed_form"]
    assert row["cross_closed_form"] < got["the_floor_is_about"]
    assert "not reproducible" in got["what_is_not_claimed"]


def test_the_downstream_matching_is_unchanged_by_the_repair():
    """`areal.solve_matching` consumes only ``ℓ = 0`` and ``ℓ = 1``, both of
    which the solve did resolve -- so the repair must *not* move the #265
    answer.  If it did, the old number was being carried by the broken term."""
    m = areal.INTERFERENCE_MOMENTS[1]
    a = m.radius
    basis = areal.basis_channels(areal.MOUTHS, a)
    throat = pt.VacuumThroat(mouth_radius=a)

    class _Riccati:
        def admittance(self, ell):
            return throat.admittance_riccati(int(ell))

    new = np.asarray(areal.solve_matching(
        areal.MOUTHS, a, throat, m.as_source(), m.signed_obstruction(),
        basis=basis)["areal_response"], float)
    old = np.asarray(areal.solve_matching(
        areal.MOUTHS, a, _Riccati(), m.as_source(), m.signed_obstruction(),
        basis=basis)["areal_response"], float)
    assert np.max(np.abs(new / old - 1.0)) < 1e-8
    assert np.all(new > 0.0), "both mouths still open"
