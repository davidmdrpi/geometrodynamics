"""The mouth matching and the signed areal response.

Every check that can be run in both angular sectors is run in both.  That is
the direct lesson of PR #263's review: `neck`'s ``ℓ(ℓ+2)`` error survived a
full suite because every closed-form check in it was ``ℓ = 0``, so the suite
was structurally blind to the term that was wrong.
"""

import math

import numpy as np
import pytest

from geometrodynamics.waves import areal
from geometrodynamics.waves.initial_data import (CONSTRAINT_EIGENVALUE,
                                                 KERNEL_PROJECTOR,
                                                 regularised_green)


# ── the Green function and its tail ─────────────────────────────────────────
def test_the_derivative_is_the_closed_form_of_the_green_function():
    x = np.linspace(0.2, 2.9, 21)
    h = 1e-6
    numeric = (regularised_green(x + h) - regularised_green(x - h)) / (2 * h)
    assert np.max(np.abs(areal.regularised_green_derivative(x) - numeric)) < 1e-7


def test_the_green_functions_tail_is_the_kernel_projector():
    got = areal.measure_the_kernel_projector_is_the_green_functions_own_tail()
    assert got["it_is_the_projector"]
    assert got["residual_matches_the_tail"] < 1e-5
    assert abs(got["kernel_norm_squared"] - math.pi ** 2 / 2) < 1e-9
    assert abs(got["one_over_the_norm"] - KERNEL_PROJECTOR) < 1e-12
    assert abs(got["unit_source_normalisation"] - 1.0) < 1e-5


def test_the_projector_kernel_is_the_cosine_of_the_geodesic_angle():
    """``Σ_A x^A y^A = cos χ(x,y)`` — the identity the projector rests on."""
    rng = np.random.default_rng(11)
    x = rng.standard_normal((40, 4)); x /= np.linalg.norm(x, axis=1)[:, None]
    y = rng.standard_normal((40, 4)); y /= np.linalg.norm(y, axis=1)[:, None]
    chi = np.arccos(np.clip(np.einsum("pa,pa->p", x, y), -1, 1))
    assert np.max(np.abs(np.einsum("pa,pa->p", x, y) - np.cos(chi))) < 1e-12


# ── the homogeneous pairs, in both sectors ──────────────────────────────────
@pytest.mark.parametrize("ell", [0, 1])
def test_both_homogeneous_solutions_solve_the_radial_equation(ell):
    f1, f2 = areal._homogeneous_pair(ell)
    cent = 0.0 if ell == 0 else 2.0
    for f in (f1, f2):
        for x in (0.4, 1.0, 1.9, 2.6):
            h = 1e-4
            lhs = ((f(x + h) - 2 * f(x) + f(x - h)) / h ** 2
                   + 2 / math.tan(x) * (f(x + h) - f(x - h)) / (2 * h)
                   + (3.0 - cent / math.sin(x) ** 2) * f(x))
            assert abs(lhs) < 1e-5 * max(1.0, abs(f(x)))


@pytest.mark.parametrize("ell", [0, 1])
def test_the_regular_member_is_a_k_equals_one_harmonic(ell):
    """Both sectors are degenerate, not just ``ℓ = 0``.

    ``cos χ`` and ``sin χ`` are both ``k = 1`` harmonics, and ``4 = (n+1)²``
    at ``n = 1``.  The three ``ℓ = 1`` partners are exactly as degenerate as
    the ``ℓ = 0`` one, and were hidden while the centrifugal term was wrong.
    """
    f1, _ = areal._homogeneous_pair(ell)
    x = np.linspace(0.05, math.pi - 0.05, 400)
    expected = np.cos(x) if ell == 0 else np.sin(x)
    assert np.max(np.abs([f1(float(t)) for t in x] - expected)) < 1e-12
    assert CONSTRAINT_EIGENVALUE == 4.0


@pytest.mark.parametrize("ell", [0, 1])
def test_the_wronskian_combination_is_one(ell):
    f1, f2 = areal._homogeneous_pair(ell)
    h = 1e-6
    for x in (0.5, 1.3, 2.4):
        d1 = (f1(x + h) - f1(x - h)) / (2 * h)
        d2 = (f2(x + h) - f2(x - h)) / (2 * h)
        assert abs(math.sin(x) ** 2 * (f1(x) * d2 - f2(x) * d1) - 1.0) < 1e-6


# ── frames, spheres and channels ────────────────────────────────────────────
def test_the_tangent_frame_is_orthonormal_and_perpendicular_to_the_mouth():
    for c in areal.MOUTHS + ((0.3, -0.4, 0.5, 0.7),):
        c = np.asarray(c, float); c /= np.linalg.norm(c)
        t = areal.tangent_frame(c)
        assert np.max(np.abs(t @ t.T - np.eye(3))) < 1e-13
        assert np.max(np.abs(t @ c)) < 1e-13


def test_the_tangent_frame_is_deterministic_and_not_svd_derived():
    doc = areal.tangent_frame.__doc__ or ""
    assert "Householder" in doc
    import inspect
    body = inspect.getsource(areal.tangent_frame).replace(doc, "")
    assert "svd" not in body.lower(), "the code, not the docstring, must be SVD-free"
    c = (0.3, -0.4, 0.5, 0.7)
    assert np.array_equal(areal.tangent_frame(c), areal.tangent_frame(c))


def test_the_direction_rule_weights_sum_to_one_and_kill_the_dipole():
    dirs, w = areal.direction_rule()
    assert abs(float(w.sum()) - 1.0) < 1e-13
    assert np.max(np.abs(np.einsum("p,pi->i", w, dirs))) < 1e-13
    assert np.max(np.abs(3 * np.einsum("p,pi,pj->ij", w, dirs, dirs)
                         - np.eye(3))) < 1e-12


def test_the_mouth_sphere_sits_at_the_stated_geodesic_radius():
    dirs, _ = areal.direction_rule(6, 8)
    for c in areal.MOUTHS:
        c = np.asarray(c, float)
        pts = areal.mouth_sphere(c, 0.05, dirs)
        assert np.max(np.abs(np.linalg.norm(pts, axis=1) - 1.0)) < 1e-13
        chi = np.arccos(np.clip(pts @ c, -1, 1))
        assert np.max(np.abs(chi - 0.05)) < 1e-12


def test_the_channel_split_separates_the_two_sectors_exactly():
    dirs, w = areal.direction_rule()
    m = np.array([0.6, -0.8, 0.0])
    v0, v1 = areal.channel_split(np.full(len(dirs), 2.5), dirs, w)
    assert abs(v0 - 2.5) < 1e-13 and np.max(np.abs(v1)) < 1e-12
    v0, v1 = areal.channel_split(dirs @ m, dirs, w)
    assert abs(v0) < 1e-13 and np.max(np.abs(v1 - m)) < 1e-12


# ── the gluing ──────────────────────────────────────────────────────────────
def test_parallel_transport_is_a_rotation_carrying_one_mouth_to_the_other():
    c1, c2 = [np.asarray(c, float) for c in areal.MOUTHS]
    p = areal.parallel_transport(c1, c2)
    assert np.max(np.abs(p @ p.T - np.eye(4))) < 1e-12
    assert abs(np.linalg.det(p) - 1.0) < 1e-12
    assert np.max(np.abs(p @ c1 - c2)) < 1e-12


def test_the_gluing_matrix_is_orthogonal_in_both_variants():
    for reflect in (False, True):
        r = areal.gluing_matrix(*areal.MOUTHS, reflect=reflect)
        assert np.max(np.abs(r @ r.T - np.eye(3))) < 1e-12


# ── the assembly against exact solves, in both sectors ──────────────────────
@pytest.mark.slow
def test_the_matching_reproduces_exact_radial_solves_in_both_sectors():
    got = areal.measure_the_matching_reproduces_an_exact_radial_solve()
    assert got["both_sectors_agree"], got["worst_overall"]
    assert got["worst_l0"] < 1e-8
    assert got["worst_l1"] < 1e-8
    assert len(got["rows"]) == 6
    assert {r["sector"] for r in got["rows"]} == {"l=0", "l=1"}


def test_the_stencil_is_fourth_order_because_two_point_was_not_enough():
    """The two-point rule left a systematic ``1e-06`` against the reference.

    Its relative truncation error on a ``1/χ`` field is exactly ``step²``, so
    the discrepancy was the assembly's, not the reference's — which is why the
    reference-side tolerances refused to move it.
    """
    import inspect
    body = inspect.getsource(areal.basis_channels)
    assert "8.0" in body and "12.0" in body
    assert "(-2, -1, 1, 2)" in body


# ── the matching system ─────────────────────────────────────────────────────
def test_the_system_is_square_and_solved_cleanly():
    m = areal.INTERFERENCE_MOMENTS[1]
    got = areal.solve_matching(areal.MOUTHS, m.radius, areal.WORKING_TUBE,
                               m.as_source(), m.signed_obstruction())
    assert got["residual"] < 1e-10
    assert got["condition_number"] < 1e8
    assert np.asarray(got["coefficients"]).shape == (12,)


def test_the_dipole_coefficients_are_tangent_to_their_own_mouths():
    m = areal.INTERFERENCE_MOMENTS[1]
    got = areal.solve_matching(areal.MOUTHS, m.radius, areal.WORKING_TUBE,
                               m.as_source(), m.signed_obstruction())
    for d, c in zip(np.asarray(got["dipoles"]), areal.MOUTHS):
        assert abs(float(d @ np.asarray(c, float))) < 1e-12


def test_the_solvability_condition_holds_on_the_solution():
    m = areal.INTERFERENCE_MOMENTS[1]
    got = areal.solve_matching(areal.MOUTHS, m.radius, areal.WORKING_TUBE,
                               m.as_source(), m.signed_obstruction())
    total = sum(a * np.asarray(c, float)
                for a, c in zip(np.asarray(got["monopoles"]), areal.MOUTHS))
    total = total + np.asarray(got["dipoles"]).sum(axis=0)
    assert np.max(np.abs(total - m.signed_obstruction())) < 1e-12


def test_the_measured_moments_are_perpendicular_to_their_own_mouths():
    """A structural check on the source quadrature, free and worth taking.

    ``u1_j`` is a gradient at ``c_j``, so it has no component along ``c_j``.
    """
    for m in areal.INTERFERENCE_MOMENTS:
        for g, c in zip(m.gradient, areal.MOUTHS):
            assert abs(float(np.asarray(g, float) @ np.asarray(c, float))) < 1e-16


# ── what the answer is made of ──────────────────────────────────────────────
def test_the_dipole_layers_are_required_for_the_problem_to_have_a_solution():
    got = areal.measure_the_dipole_layers_are_required_not_optional()
    assert got["monopole_span_dimension"] == 2
    assert got["with_dipoles_dimension"] == 4
    assert got["monopoles_alone_cannot_close_it"]
    assert got["the_dipoles_close_it"]
    assert got["monopole_only_shortfall"] > 0.6


def test_and_then_the_off_plane_obstruction_moves_the_answer_not_at_all():
    got = areal.measure_the_dipole_layers_are_required_not_optional()
    assert got["and_then_they_do_not_move_the_answer"]
    assert max(abs(v) for v in got["response_to_the_off_plane_part"]) < 1e-15
    whole = np.asarray(got["response_to_the_whole"])
    plane = np.asarray(got["response_to_the_in_plane_part"])
    assert np.max(np.abs(plane / whole - 1.0)) < 1e-10


def test_the_obstruction_carries_the_answer_and_the_local_data_do_not():
    got = areal.measure_the_obstruction_carries_the_answer()
    assert got["the_obstruction_is_the_answer"]
    assert got["obstruction_alone_reproduces_it"] < 5e-3
    assert got["without_the_obstruction"] < 1e-2
    assert got["worst_drift_from_the_dipole_moments"] < 5e-3
    assert got["linear_in_the_obstruction"] < 1e-10


# ── the headline ────────────────────────────────────────────────────────────
def test_the_areal_response_is_negative_at_both_mouths():
    got = areal.measure_the_signed_areal_response()
    assert got["every_variant_agrees_in_sign"]
    assert got["sign"] == ["closes", "closes"]
    assert all(v < 0 for v in got["areal_response"])
    assert got["quadrature_spread_at_fixed_radius"] < 0.05
    assert got["worst_residual"] < 1e-10


def test_the_gluing_cannot_reach_the_areal_response():
    """A full ``2π`` twist of the transverse frames, not just a reflection."""
    m = areal.INTERFERENCE_MOMENTS[1]
    basis = areal.basis_channels(areal.MOUTHS, m.radius)
    plain = areal.gluing_matrix
    base = None
    norms = []
    try:
        for th in (0.0, 0.7, 1.9, 3.4, 5.0):
            ct, st = math.cos(th), math.sin(th)
            twist = np.array([[1, 0, 0], [0, ct, -st], [0, st, ct]], float)
            areal.gluing_matrix = (
                lambda a, b, reflect=False, _t=twist: plain(a, b, reflect) @ _t)
            got = areal.solve_matching(areal.MOUTHS, m.radius,
                                       areal.WORKING_TUBE, m.as_source(),
                                       m.signed_obstruction(), basis=basis)
            got_v = np.asarray(got["areal_response"], float)
            norms.append(float(np.linalg.norm(got["dipoles"][0])))
            if base is None:
                base = got_v
            else:
                assert np.max(np.abs(got_v / base - 1.0)) < 1e-10
    finally:
        areal.gluing_matrix = plain
    assert max(norms) / min(norms) - 1.0 > 1e-3, "the twist must actually bite"


def test_the_conformal_factor_and_the_areal_response_differ_by_four():
    m = areal.INTERFERENCE_MOMENTS[1]
    got = areal.solve_matching(areal.MOUTHS, m.radius, areal.WORKING_TUBE,
                               m.as_source(), m.signed_obstruction())
    assert np.max(np.abs(np.asarray(got["areal_response"])
                         - 4.0 * np.asarray(got["conformal_factor"]))) < 1e-15


def test_the_answer_is_linear_in_the_coupling():
    got = areal.measure_the_signed_areal_response(coupling=2.0)
    once = areal.measure_the_signed_areal_response(coupling=1.0)
    assert np.max(np.abs(np.asarray(got["areal_response"])
                         / np.asarray(once["areal_response"]) - 2.0)) < 1e-9


# ── the throat is a cavity ──────────────────────────────────────────────────
def test_the_tube_channels_share_a_rate_and_differ_in_kind():
    t = areal.WORKING_TUBE
    assert abs(t.wavenumber() - 1.0) < 1e-15
    assert abs(float(t.monopole_rate()) - t.dipole_rate()) < 1e-15
    _, cc, ss = t.monopole_transfer()
    assert abs(cc - math.cos(0.9)) < 1e-13 and abs(ss - math.sin(0.9)) < 1e-13
    _, ch, sh = t.dipole_transfer()
    assert abs(ch - math.cosh(0.9)) < 1e-13 and abs(sh - math.sinh(0.9)) < 1e-13


def test_the_tube_area_sets_the_wavenumber_through_the_scalar_curvature():
    for area in (math.pi, 4 * math.pi, 16 * math.pi):
        t = areal.TubeModel(area=area)
        r = math.sqrt(area / (4 * math.pi))
        assert abs(t.wavenumber() - 1.0 / r) < 1e-13


@pytest.mark.slow
def test_the_sign_flips_at_the_throats_own_standing_waves():
    got = areal.measure_the_throat_is_a_resonant_cavity()
    assert got["flips_land_on_the_poles"]
    assert got["the_working_throat_is_off_resonance"]
    assert abs(got["working_phase"] - 0.9) < 1e-12
    for d in got["distance_to_the_closed_form"][:2]:
        assert d <= 2 * got["grid_spacing"]


# ── the headline is a statement about a wide tube ───────────────────────────
def test_the_matched_tube_wavenumber_is_one_over_sin_a():
    """`A = 4 pi sin^2 a` makes the tube as narrow as its own mouths."""
    for a in (0.05, 0.10, 0.2):
        t = areal.TubeModel(area=4 * math.pi * math.sin(a) ** 2)
        assert abs(t.wavenumber() - 1.0 / math.sin(a)) < 1e-12


def test_the_sign_does_not_survive_a_matched_tube():
    """The result is conditional on the throat, and this is the counterexample."""
    got = areal.measure_the_sign_does_not_survive_a_matched_tube()
    assert got["the_wide_throat_always_closes"]
    assert got["the_matched_throat_does_not"]
    assert got["the_matched_mouths_can_disagree"]
    assert got["the_sign_is_a_property_of_the_throat"]
    assert got["worst_residual"] < 1e-9
    # and it must be a real answer, not a singular system read out anyway
    assert got["worst_condition_number"] < 1e10
    phases = {r["radius"]: r["phase_over_pi"]
              for r in got["rows"] if r["model"] == "matched"}
    assert abs(phases[0.05] - 5.732) < 5e-3
    assert abs(phases[0.10] - 2.870) < 5e-3


def test_the_evanescent_channel_is_written_without_forming_a_growing_exponential():
    """A `cosh`/`sinh` transfer matrix costs a condition number of `e^{2 kL}`.

    At the matched area and `a = 0.05` that is `4.4e+15` — the system goes
    singular to double precision, and no answer can be read out of it. Carrying
    the tube's two end amplitudes as unknowns instead keeps every coefficient
    bounded by one.
    """
    t = areal.TubeModel(area=4 * math.pi * math.sin(0.05) ** 2, length=0.9)
    k, x = t.dipole_attenuation()
    assert 0.0 < x <= 1.0
    assert abs(x - math.exp(-k * 0.9)) < 1e-300 + 1e-15 * x
    m = areal.INTERFERENCE_MOMENTS[1]
    got = areal.solve_matching(areal.MOUTHS, 0.05, t, m.as_source(),
                               m.signed_obstruction())
    assert got["condition_number"] < 1e10, "the stable form must stay solvable"
    assert got["residual"] < 1e-9
    assert np.asarray(got["tube_amplitudes"]).shape == (6,)


def test_the_stable_form_agrees_with_the_transfer_matrix_where_both_work():
    """Where `cosh kL` is not yet enormous, the two formulations must agree.

    `dipole_transfer` is kept precisely so this comparison exists: the
    reformulation was made for conditioning, so it has to be shown to be a
    change of form and not of content.
    """
    got = areal.measure_the_matching_reproduces_an_exact_radial_solve()
    assert got["both_sectors_agree"]
    assert got["worst_overall"] < 1e-8
