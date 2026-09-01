"""Does negative-coupling EGB actually work? Checked against
`docs/negative_egb_prereg.md`, frozen before the module.

Four kinds. The *exterior* ones pin `R_kk`, `H_kk` and the `χ`-independence
that makes the outside constraint a single number. The *search* one hunts for a
surviving interval of `α_GB` rather than asserting there is none — a surviving
interval would reopen the branch. The *mechanism* ones pin that the two bounds
meet by continuity of one bracket, not by coincidence. And the *cost* ones pin
what the single surviving coupling does to the observed universe.
"""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import finite_mouth as fm
from geometrodynamics.bulk import gauss_bonnet as gb
from geometrodynamics.bulk import negative_egb as ne


# ── the exterior ────────────────────────────────────────────────────────────

def test_the_exterior_null_ricci_is_positive_unlike_the_throat():
    """`+3/R²` outside against `−3f″/f` inside: the sign difference is the
    whole mechanism."""
    for R in (1.0, 1.3, 2.5):
        assert ne.exterior_ricci_null(R) == pytest.approx(3.0 / R ** 2, rel=1e-14)
        assert ne.exterior_ricci_null(R) > 0.0
    g = fm.geometry()
    assert float(gb.ricci_null_contraction(np.array([0.0]), g.neck_radius)[0]) < 0.0


def test_the_exterior_lanczos_null_is_twelve_over_R_to_the_fourth():
    for R in (1.0, 1.3, 2.5):
        assert ne.exterior_lanczos_null(R) == pytest.approx(12.0 / R ** 4, rel=1e-14)


def test_the_exterior_ratio_is_chi_independent():
    """`H_kk/R_kk = 4μ/f⁴ = 4/R²` for every `χ`, computed from the exterior's
    own `μ` and `f` rather than asserted. That constancy is what makes the
    outside constraint one number instead of a profile."""
    for R in (1.0, 2.5):
        chi = np.linspace(0.05, math.pi - 0.05, 101)
        ratio = ne.exterior_ratio(chi, R)
        assert np.ptp(ratio) < 1e-12
        assert float(ratio[0]) == pytest.approx(4.0 / R ** 2, rel=1e-12)


def test_the_exterior_null_stress_and_its_threshold():
    for R in (1.0, 2.5):
        for coupling in (0.0, -0.1, -0.25 * R ** 2, -R ** 2):
            assert ne.exterior_matter_null_stress(coupling, R) == pytest.approx(
                3.0 * (R ** 2 + 4.0 * coupling) / R ** 4, rel=1e-13)
        assert ne.exterior_threshold(R) == pytest.approx(-0.25 * R ** 2, rel=1e-14)
        assert ne.exterior_matter_null_stress(
            ne.exterior_threshold(R), R) == pytest.approx(0.0, abs=1e-13)


def test_the_two_bounds_are_the_same_number_pointing_opposite_ways():
    result = ne.measure_the_exterior_constrains_alpha_oppositely()
    assert result["the_two_bounds_coincide"]
    assert result["exterior_needs_alpha_at_least"] == pytest.approx(
        result["throat_needs_alpha_at_most"], rel=1e-14)


# ── the search ──────────────────────────────────────────────────────────────

def test_no_open_set_of_couplings_satisfies_both_regions():
    """The falsifier: a surviving interval would reopen the branch."""
    result = ne.measure_no_coupling_satisfies_both()
    assert result["the_surviving_set_is_measure_zero"]
    assert result["surviving_width"] < 1e-6
    assert result["survivors_are_at_the_critical_value"]
    assert result["samples"] >= 1000


def test_the_scan_actually_covers_both_signs_of_the_coupling():
    """A scan that never went positive would beg the question."""
    result = ne.measure_no_coupling_satisfies_both()
    low, high = result["scanned_range"]
    assert low < result["critical_coupling"] < 0.0 < high
    assert result["throat_satisfied_count"] > 0
    assert result["exterior_satisfied_count"] > 0


def test_the_throat_and_exterior_admissible_sets_are_complementary():
    """Directly: at any coupling one side or the other fails, except one."""
    R, a = 1.0, 0.3
    g = fm.geometry(R, a)
    s = np.linspace(-g.half_length, g.half_length, 201)
    for coupling in (-1.0, -0.5, -0.3, -0.2, -0.1, 0.0, 0.5):
        throat = bool(np.all(gb.matter_null_stress(s, g.neck_radius, coupling)
                             >= -1e-9))
        exterior = ne.exterior_matter_null_stress(coupling, R) >= -1e-9
        assert not (throat and exterior), \
            f"alpha={coupling} would reopen the branch"


# ── the mechanism ───────────────────────────────────────────────────────────

def test_the_ratio_is_continuous_across_the_seam():
    """`μ/f⁴ = 1/R²` from both sides — the same Misner–Sharp continuity as P2."""
    result = ne.measure_the_bracket_is_continuous_at_the_seam()
    assert result["ratio_is_continuous"]
    assert result["ratio_jump"] < 1e-12
    assert result["ratio_inside_at_the_seam"] == pytest.approx(
        result["common_value"], rel=1e-12)


def test_the_shared_bracket_vanishes_at_the_seam_at_the_critical_coupling():
    """It must be `≤ 0` from inside and `≥ 0` from outside while being
    continuous, so it is exactly `0`."""
    result = ne.measure_the_bracket_is_continuous_at_the_seam()
    assert result["bracket_vanishes_at_the_seam"]
    assert abs(result["bracket_inside"]) < 1e-12
    assert abs(result["bracket_outside"]) < 1e-12


def test_the_bracket_has_the_required_signs_on_each_side():
    """Below the critical coupling the throat is satisfied and the exterior is
    not; above, the reverse."""
    R, a = 1.0, 0.3
    g = fm.geometry(R, a)
    inside = float(gb.gauss_bonnet_ratio(np.array([g.half_length]),
                                         g.neck_radius)[0])
    outside = float(ne.exterior_ratio(np.array([a]), R)[0])
    critical = ne.exterior_threshold(R)
    assert float(ne.bracket(1.2 * critical, inside)) < 0.0     # throat happy
    assert float(ne.bracket(1.2 * critical, outside)) < 0.0    # exterior not
    assert float(ne.bracket(0.8 * critical, inside)) > 0.0     # throat not
    assert float(ne.bracket(0.8 * critical, outside)) > 0.0    # exterior happy


# ── the cost ────────────────────────────────────────────────────────────────

def test_gauss_bonnet_cannot_move_the_exterior_pressure():
    """`H^i_j = 0` on a maximally symmetric slice, so the whole correction
    lands in `ρ`."""
    for R in (1.0, 2.5):
        base = ne.exterior_pressure(0.0, R)
        for coupling in (-2.0, -0.25 * R ** 2, 0.0, 3.0):
            assert ne.exterior_pressure(coupling, R) == pytest.approx(
                base, rel=1e-15)
            assert base == pytest.approx(-3.0 / R ** 2, rel=1e-14)
    assert ne.measure_the_critical_exterior_is_empty()[
        "pressure_is_coupling_independent"]


def test_the_critical_exterior_is_exactly_vacuum_energy():
    """`w = −1` exactly, at half the Einstein-gravity density."""
    result = ne.measure_the_critical_exterior_is_empty()
    assert result["it_is_exactly_vacuum_energy"]
    assert result["critical_equation_of_state"] == pytest.approx(-1.0, abs=1e-12)
    assert result["density_is_halved"]
    for R in (1.0, 2.5):
        critical = ne.exterior_threshold(R)
        assert ne.exterior_density(critical, R) == pytest.approx(
            3.0 / R ** 2, rel=1e-13)


def test_pushing_the_coupling_further_makes_the_exterior_exotic():
    """The exotic matter is relocated into the universe, not removed."""
    R = 1.0
    critical = ne.exterior_threshold(R)
    assert ne.exterior_matter_null_stress(1.5 * critical, R) < 0.0
    rows = {round(r["coupling"], 9): r
            for r in ne.measure_the_critical_exterior_is_empty()["rows"]}
    beyond = rows[round(1.5 * critical, 9)]
    assert beyond["sum"] < 0.0
    assert beyond["equation_of_state"] < -1.0


# ── the verdict ─────────────────────────────────────────────────────────────

def test_the_branch_is_reported_closed_on_physical_grounds():
    ledger = ne.measure_the_negative_egb_ledger()
    assert ledger["branch_is_closed"]
    assert "PHYSICAL grounds" in ledger["verdict"]
    assert "coupling constant" in " ".join(
        e["verdict"] for e in ledger["entries"]).lower()


def test_the_ledger_names_what_it_does_not_test():
    """Stability, the graviton kinetic term, dilatonic EGB, `f(R)`, and a
    different exterior."""
    note = ne.measure_the_negative_egb_ledger()["what_remains_untested"]
    assert "Stability" in note and "graviton kinetic term" in note
    assert "dilatonic" in note
    assert "f(R)" in note
    assert "DIFFERENT exterior" in note


def test_the_ledger_keeps_five_branches_and_derives_its_numbers():
    ledger = ne.measure_the_negative_egb_ledger()
    assert len(ledger["the_remaining_branches"]) == 5
    scan = ne.measure_no_coupling_satisfies_both()
    entries = {e["claim"]: e for e in ledger["entries"]}
    row = next(e for k, e in entries.items() if "open interval" in k)
    assert f"{scan['critical_coupling']:.6f}" in row["evidence"]


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import negative_egb_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
