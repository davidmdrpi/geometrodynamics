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
    for a in (0.05, 0.10):
        t = pt.VacuumThroat(mouth_radius=a)
        num = t.admittance(0)
        closed = t.monopole_admittance_closed_form()
        assert np.max(np.abs(num - closed)) < 1e-10 * np.max(np.abs(closed))


def test_the_monopole_admittance_is_singular_because_flux_is_conserved():
    """``(f²u')' = 0`` makes ``Y·(1,1) = 0`` an identity, not a cancellation."""
    for a in (0.02, 0.05, 0.10, 0.20):
        t = pt.VacuumThroat(mouth_radius=a)
        y = t.admittance(0)
        scale = float(np.max(np.abs(y)))
        # the claim is that det vanishes *relative to the matrix*, so the
        # tolerance has to scale with it -- an absolute one is a different,
        # stronger claim that happens to hold only at small a
        assert abs(float(np.linalg.det(y))) < 1e-14 * scale ** 2
        assert np.max(np.abs(y @ np.ones(2))) < 1e-14 * scale
        assert abs(t.shunt()) < 1e-14 * scale


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
    """The claim is strong, so the four things it does not say are asserted."""
    got = pt.measure_the_throat_is_an_einstein_rosen_bridge()
    assert got["it_is_an_einstein_rosen_bridge"]
    assert got["the_mass_has_no_free_parameter"]
    assert got["the_neck_is_an_apparent_horizon"]
    assert got["schwarzschild_slope_error"] < 1e-15
    assert got["three_masses_agree"] < 1e-12
    assert got["the_gluing_is_hawking_mass_continuity"] < 1e-12
    caveats = " ".join(got["what_it_is_not"]).lower()
    for word in ("adm", "dimensionless", "handle", "apparent horizon"):
        assert word in caveats
