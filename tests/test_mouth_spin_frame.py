"""The mouth spin-frame round, against `docs/mouth_spin_frame_prereg.md` (T1-T6)."""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import mouth_spin_frame as sf


def test_T1_the_frame_map_is_the_spin_double_cover_of_the_frame_bundle():
    c = sf.frame_map_checks()
    assert c["orthonormal_error"] < 1e-12 and c["tangent_error"] < 1e-12
    assert c["orientation_error"] < 1e-12
    assert c["q_and_minus_q_coincide"] < 1e-12
    assert c["double_cover"]


def test_T2_fibre_angle_phi_rotates_the_frame_by_2phi():
    r = sf.fibre_rotates_frame_twice()
    assert r["base_point_fixed"] < 1e-12
    assert r["frame_angle_minus_two_phi"] < 1e-10
    assert r["angle_is_2phi"]


def test_T2_levi_civita_form_is_twice_the_hopf_connection():
    r = sf.levi_civita_versus_hopf_connection()
    assert r["worst_residual_of_omega_plus_2A"] < 1e-8
    assert r["omega_is_minus_twice_A"]


def test_T2_chern_number_one_against_euler_number_two():
    r = sf.chern_versus_euler()
    assert r["c1_hopf"] == pytest.approx(1.0, abs=1e-6)
    assert r["euler_TS2"] == pytest.approx(2.0, abs=1e-6)
    assert r["ratio_is_two"]


def test_T3_exactly_the_units_perpendicular_to_the_fibre_generator_cover_the_antipode():
    d = sf.deck_lifts_of_antipode()
    assert d["exactly_the_perp_units_cover"]
    for row in d["rows"]:
        assert row["L_u_squared_is_minus_one"]
        if row["perp_to_i"]:
            assert row["covers_antipode"] and row["image_frame_oriented"]
    assert d["frame_image_is_dA_times_one_reflection"]
    assert d["fibre_rotation_conjugates_within_perp_circle"]


def test_T3_T4_two_sheeted_involution_and_two_pin_minus_structures():
    t = sf.two_sheeted_involution()
    for eps in ("epsilon=+1", "epsilon=-1"):
        assert t[eps]["A_tilde_squared_is_identity"]
        assert t[eps]["single_sheet_L_u_squared_is_minus_one"]
    assert t["number_of_pin_minus_structures"] == 2
    assert not t["preferred_by_geometry"]


def test_T5_abk_invariants_are_plus_and_minus_one_mod_8():
    s = sf.pin_minus_structures_rp2()
    assert [r["ABK_mod8"] for r in s] == [1, 7]
    for r in s:
        assert r["is_quadratic"] and r["modulus_is_sqrt2"] and r["H1_action_shifts_q_by_2"]


def test_T6_opposite_sectors_bound_equal_sectors_and_single_mouths_do_not():
    r = sf.pin_bordism_pairing_rule()
    assert r["opposite_sectors_bound"]
    assert r["equal_sectors_do_not"]
    assert r["single_mouth_cannot_bound"]


def test_T7_no_quantum_objects_are_imported():
    """The module must not import the Bell, history or transaction packages,
    nor call any singlet, Born-rule or CHSH function. Checked on the import
    lines and the call names, not on the prose (which names the exclusions)."""
    import ast
    import inspect
    tree = ast.parse(inspect.getsource(sf))
    imported = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            imported.append(node.module or "")
        elif isinstance(node, ast.Import):
            imported.extend(a.name for a in node.names)
    for name in imported:
        assert not any(b in name for b in ("bell", "history", "transaction"))
    calls = {n.func.attr if isinstance(n.func, ast.Attribute) else getattr(n.func, "id", "")
             for n in ast.walk(tree) if isinstance(n, ast.Call)}
    for banned in ("singlet", "born", "chsh", "bulk_correlation", "bulk_chsh"):
        assert not any(banned in c.lower() for c in calls)
