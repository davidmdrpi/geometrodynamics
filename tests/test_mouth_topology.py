"""The finite-mouth topology round, against `docs/finite_mouth_topology_prereg.md`.

Committed after the pre-registration and checked against its numbers. Six
groups, each of which can fail: topology (H1), the involution (H2), the
bundle lift (H5), the analytic DtN oracle (B2), convergence, and the two
controls (topology and antipodal).
"""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import finite_mouth as fm
from geometrodynamics.bulk import mouth_topology as mt

TARGETS = {  # from the pre-registration, R = 1, a = 0.3
    (0, "Neumann"): 0.000000000, (0, "Dirichlet"): 0.157587622,
    (1, "Neumann"): 1.797266559, (1, "Dirichlet"): 1.804461992,
    (2, "Neumann"): 3.524607516, (2, "Dirichlet"): 3.524854051,
    (3, "Neumann"): 5.248595411, (3, "Dirichlet"): 5.248602920,
}


# ── H1: topology of the two-mouth handle ────────────────────────────────────

def test_antipodal_gluing_is_orientable_in_the_bulk():
    """T2: `det(-I_4) = +1`, so the antipodally glued handle has `w_1 = 0`."""
    c = mt.classify_gluing(mt.antipodal_gluing())
    assert c.det_bulk == +1 and c.bulk_orientable
    assert mt.mapping_torus_w1(mt.antipodal_gluing()) == +1


def test_antipodal_gluing_twists_the_brane_and_its_normal():
    """T3: the antipodal class is `(det m_3, eps) = (-, -)`."""
    c = mt.classify_gluing(mt.antipodal_gluing())
    assert c.det_brane == -1 and c.normal_sign == -1
    assert not c.brane_orientable and not c.normal_bundle_trivial


def test_the_four_classes_are_realised_and_w1_is_their_product():
    """T3: `w_1(bulk) = det m_3 * eps` on every class."""
    seen = set()
    for brane_flip in (False, True):
        for normal_flip in (False, True):
            m = np.eye(4)
            if brane_flip:
                m[0, 0] = -1.0
            if normal_flip:
                m[3, 3] = -1.0
            c = mt.classify_gluing(m, normal_index=3)
            assert c.det_bulk == c.det_brane * c.normal_sign
            seen.add((c.det_brane, c.normal_sign))
    assert len(seen) == 4


def test_harmonic_action_of_the_antipode_is_the_parity_computed():
    """T4: `D_l(-I) = (-1)^l I` comes out of a least-squares fit, with the
    right dimension `(l+1)^2`."""
    for ell in range(5):
        rep = mt.harmonic_representation(mt.antipodal_gluing(), ell)
        assert rep["dimension"] == (ell + 1) ** 2
        assert rep["fit_residual"] < 1e-9
        assert rep["is_scalar"] and rep["scalar"] == pytest.approx((-1) ** ell)


def test_a_reflection_does_not_act_as_a_scalar_on_harmonics():
    """The reflection class is not a parity grading: it mixes within `l`."""
    rep = mt.harmonic_representation(mt.reflection_gluing(3), 2)
    assert not rep["is_scalar"]
    assert rep["fit_residual"] < 1e-9


# ── H2: the involution ──────────────────────────────────────────────────────

def test_the_antipodal_involution_is_free_isometric_and_orientation_reversing():
    """I1: free, an isometry for both lapses, `det d iota = -1` on the bulk
    and `+1` on the brane slice."""
    inv = mt.antipodal_involution()
    assert inv.is_involution() and inv.is_free
    assert inv.is_isometry(0.09, fm.lapse_ultrastatic)
    assert inv.is_isometry(0.09, fm.lapse_vacuum)
    assert inv.tangent_determinant() == -1
    assert inv.brane_determinant() == +1


def test_only_minus_identity_is_free_on_the_four_sphere():
    """I2: of the six involution classes of `O(5)`, exactly one is free."""
    rows = mt.free_involutions_of_o5()
    free = [r for r in rows if r["free_on_S4"]]
    assert len(free) == 1 and free[0]["minus_ones"] == 5


def test_pin_types_of_the_quotient_pieces():
    """I3: `RP^2` is Pin^- only, `RP^3` is spin, `RP^4` is Pin^+ only."""
    assert mt.pin_structures_rp(2) == {"n": 2, "orientable": False, "spin": False,
                                       "pin_plus": False, "pin_minus": True}
    assert mt.pin_structures_rp(3)["spin"]
    p4 = mt.pin_structures_rp(4)
    assert p4["pin_plus"] and not p4["pin_minus"] and not p4["orientable"]


def test_the_static_operator_is_parity_symmetric_for_both_lapses():
    """L1: the neck sector labels do not depend on the lapse."""
    for lapse in (fm.lapse_ultrastatic, fm.lapse_vacuum):
        assert mt.static_operator_commutes_with_parity(0.09, lapse) < 1e-8


# ── H3: the sector labels and the oracle ────────────────────────────────────

def test_the_neck_sectors_reproduce_pr129_for_eta_plus_and_swap_for_eta_minus():
    """B1."""
    assert [mt.neck_sector(l) for l in range(6)] == \
        ["Neumann", "Dirichlet"] * 3
    assert [mt.neck_sector(l, eta=-1) for l in range(6)] == \
        ["Dirichlet", "Neumann"] * 3


def test_the_half_tube_solve_hits_the_pre_registered_targets():
    """B2: independent tridiagonal solve versus the closed form, `< 1e-5`."""
    for (ell, cond), target in TARGETS.items():
        oracle = mt.half_tube_admittance_oracle(ell, cond)
        assert oracle == pytest.approx(target, abs=2e-9)
        num = mt.half_tube_admittance(ell, cond, steps=2000)
        assert num == pytest.approx(oracle, abs=1e-5 * max(1.0, abs(oracle)))


def test_the_oracle_is_the_symmetric_and_antisymmetric_sector_of_pr277():
    """The half-tube values are the `(1, ±1)` eigenvalues of the two-mouth
    `Y_l(0)`, so the quotient sector is a restriction, not a new solve."""
    for ell in range(4):
        Y = fm.static_admittance(ell)
        d, c = Y[0, 0], Y[0, 1]
        assert mt.half_tube_admittance_oracle(ell, "Neumann") == pytest.approx(d + c, rel=1e-12)
        assert mt.half_tube_admittance_oracle(ell, "Dirichlet") == pytest.approx(d - c, rel=1e-12)


def test_the_neumann_monopole_is_exactly_zero_and_dirichlet_is_twice_G():
    assert mt.half_tube_admittance_oracle(0, "Neumann") == pytest.approx(0.0, abs=1e-14)
    assert mt.half_tube_admittance_oracle(0, "Dirichlet") == \
        pytest.approx(2.0 * fm.monopole_conductance(), rel=1e-12)


def test_second_order_convergence_of_the_half_tube_solve():
    """Refinement ratios `>= 3.5` at `l = 2` Neumann and `l = 1` Dirichlet."""
    for cond, ell in (("Neumann", 2), ("Dirichlet", 1)):
        oracle = mt.half_tube_admittance_oracle(ell, cond)
        errs = [abs(mt.half_tube_admittance(ell, cond, steps=n) - oracle)
                for n in (500, 1000, 2000, 4000)]
        for e0, e1 in zip(errs, errs[1:]):
            assert e0 / e1 >= 3.5


# ── controls ────────────────────────────────────────────────────────────────

def test_topology_control_a_reflection_gluing_flips_w1():
    assert mt.mapping_torus_w1(mt.identity_gluing()) == +1
    assert mt.mapping_torus_w1(mt.reflection_gluing(3)) == -1
    assert mt.mapping_torus_w1(mt.antipodal_gluing()) == +1


def test_antipodal_control_the_grading_disappears_and_the_action_is_not_free():
    """B3: with `Omega -> Omega` the involution fixes the neck and every `l`
    gets the same condition."""
    inv = mt.neck_reflection_involution()
    assert not inv.is_free
    assert [mt.neck_sector(l, involution=inv) for l in range(6)] == ["Neumann"] * 6
    assert [mt.neck_sector(l, eta=-1, involution=inv) for l in range(6)] == ["Dirichlet"] * 6


def test_geometry_control_labels_invariant_admittances_vary():
    """B4."""
    values = []
    for a in (0.05, 0.3, 0.8, 1.2):
        assert [mt.neck_sector(l) for l in range(4)] == ["Neumann", "Dirichlet"] * 2
        values.append(mt.half_tube_admittance_oracle(1, "Dirichlet", 1.0, a))
    assert len(set(round(v, 6) for v in values)) == 4


# ── H5: the lift ────────────────────────────────────────────────────────────

def test_J_is_left_multiplication_by_minus_j_with_det_plus_one():
    """S2/S3: built from `transport.py`'s own map; equals `L_{-j}`; `det = +1`."""
    c = mt.complex_structure_commutation()
    assert c["J_equals_left_mult_by_minus_j"]
    assert c["det_J"] == pytest.approx(1.0)
    assert c["J_squared_is_minus_identity"]
    assert c["hopf_base_antipode"]


def test_J_is_antilinear_for_the_hopf_structure_and_linear_for_the_spinor_one():
    """S2: anticommutes with `L_i`, commutes with `R_i`, and in the `R_i`
    basis it is the complex-linear matrix `[[0, 1], [-1, 0]]`."""
    c = mt.complex_structure_commutation()
    assert c["anticommutes_with_L_i"] and not c["commutes_with_L_i"]
    assert c["commutes_with_R_i"]
    assert np.allclose(c["linear_matrix_in_R_i_basis"], [[0, 1], [-1, 0]])


def test_pin_lifts_of_the_neck_reflection_are_real_chirality_exchanging_and_not_unique():
    """S1/S4: four lifts (two per Pin type), each a real matrix anticommuting
    with the volume element; squares `+1` in Pin^+ and `-1` in Pin^-."""
    out = mt.pin_lifts_of_reflection()
    assert out["count"] == 4 and out["lifts_within_a_type_are_distinct"]
    for row in out["lifts"]:
        assert row["anticommutes_with_volume"] and row["is_real_matrix"]
        assert row["square"] == (+1 if row["pin_type"] == "Pin+" else -1)


def test_spin_lifts_of_the_antipode_are_the_chirality_signs():
    out = mt.spin_lifts_of_antipode()
    assert out["anticommutes_with_generators"]
    assert all(r["square_is_identity"] for r in out["lifts"])
    assert out["eigenvalues"] == [-1.0, 1.0]


def test_J_is_the_lift_of_a_rotation_not_of_iota_or_the_antipode():
    """S3, corrected (R2): `J` is `L_{-j}`, an `SO(4)` rotation with the
    Spin(4) lift `(-j, 1)`; it is antilinear only for `L_i`. The lifts of
    `iota` are odd (chirality-exchanging) and the lifts of `-I_4` are the
    volume element, so `J` is a lift of neither."""
    c = mt.complex_structure_commutation()
    assert c["det_J"] == pytest.approx(1.0)      # a rotation
    assert c["anticommutes_with_L_i"] and c["commutes_with_R_i"]
    for row in mt.pin_lifts_of_reflection()["lifts"]:
        assert row["anticommutes_with_volume"]   # odd: not a chirality-preserving map
    assert mt.spin_lifts_of_antipode()["eigenvalues"] == [-1.0, 1.0]


# ── the mouth Pin holonomy (docs/mouth_pin_holonomy_prereg.md) ──────────────

def test_P1_deck_generator_reverses_both_normals_and_reflects_the_tangent():
    h = mt.neck_holonomy()
    assert np.allclose(h["holonomy"], np.diag([-1.0, -1.0, 1.0, -1.0]), atol=1e-9)
    assert h["both_normals_reversed"] and h["tangent_is_reflection"]
    assert h["det"] == pytest.approx(-1.0, abs=1e-9)


def test_P2_ambient_pin_plus_induces_intrinsic_pin_minus_through_lambda_plus_lambda():
    w = mt.restricted_stiefel_whitney()
    assert w["w_nu"] == (1, 0, 1) and w["w_TN"] == (1, 1, 1)
    assert w["w_TM_restricted"] == (1, 1, 0)
    assert w["ambient_pin_plus_compatible"]
    assert w["intrinsic_pin_minus"] and not w["intrinsic_pin_plus"]


def test_P3_twisted_tangent_generators_are_quaternionic():
    t = mt.twisted_tangent_generators()
    assert t["normal_volume_squared"] == -1
    assert t["twisted_t1_squared"] == -1 and t["twisted_t2_squared"] == -1
    assert t["twisted_anticommute"] and t["generates_quaternions"]


def test_P4_holonomy_lift_squares_to_minus_one_in_pin_plus_only():
    lifts = mt.holonomy_lifts()
    assert lifts["square_in_pin_plus"] == -1
    assert lifts["square_in_pin_minus"] == +1
    for row in lifts["lifts"]:
        assert row["odd_anticommutes_with_volume"] and row["equals_twisted_t2"]
        assert row["normal_part_squared"] == -1


def test_P5_closed_loop_spin_holonomy_is_minus_one_from_the_unwrapped_angle():
    out = mt.closed_loop_spin_holonomy()
    assert out["equator_rotation_angle"] == pytest.approx(2.0 * math.pi, abs=1e-6)
    assert out["max_deviation_from_classical_law"] < 1e-6
    assert out["spin_holonomy_is_minus_one"]


def test_P6_mouth_holonomy_matches_J_only_at_the_intrinsic_level_up_to_U1_and_sign():
    c = mt.compare_mouth_holonomy_with_transport()
    assert not c["vector_level"]["same_O4_class"]
    assert not c["ambient_spinor_level"]["conjugate_in_Pin4"]
    intrinsic = c["intrinsic_level"]
    assert intrinsic["L_minus_j_conjugate_to_L_j_by_spin2"]
    assert intrinsic["conjugating_angle"] == pytest.approx(0.5 * math.pi, abs=1e-3)
    assert intrinsic["also_equal_up_to_sign"]
    assert intrinsic["fibre_reversing_component_all_order_four"]
    assert c["outcome"] == "B"


def test_P7_the_K_is_anticommutation_with_the_spin2_generator():
    m = mt.intrinsic_pin2_module()
    assert m["holonomy_squared_minus_identity"]
    assert m["spin2_generator_is_L_k"]
    assert m["holonomy_anticommutes_with_spin2_generator"]
    assert m["holonomy_commutes_with_right_i"]
    assert m["pin2_is_normaliser_of_spin2"]


def test_P8_the_sign_is_a_pin_structure_and_is_not_selected():
    counts = mt.pin_structure_counts()
    assert counts["RP2_pin_minus_structures"] == 2
    assert counts["RP4_connected_sum_pin_plus_structures"] == 4
    assert not counts["sign_of_holonomy_fixed_by_geometry"]
