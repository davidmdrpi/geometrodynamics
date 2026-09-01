"""Does Gauss–Bonnet reopen the finite-mouth throat?

Checked against `docs/gauss_bonnet_prereg.md`, frozen before the module.

Four kinds. The *validation* ones pin the Lanczos implementation where the
answer is already known — `H_ab ≡ 0` in `D = 4` — because everything else
depends on it. The *sign* ones pin the result that decides the branch, against
an independently differentiated Riemann tensor. The *coupling* ones pin what
the matter NEC would demand and that no positive value works. And the *scope*
ones pin what this round does **not** refute.
"""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import finite_mouth as fm
from geometrodynamics.bulk import gauss_bonnet as gb


# ── validation: D = 4 is topological ────────────────────────────────────────

def test_the_lanczos_tensor_vanishes_in_four_dimensions():
    """The Gauss–Bonnet term is a total derivative in `D = 4`."""
    result = gb.measure_the_lanczos_tensor_is_correct()
    assert result["topological_in_four_dimensions"]
    for row in result["rows"]:
        assert abs(row["lanczos_kk"]) < 1e-4


def test_the_four_dimensional_zero_is_not_a_vacuum_artefact():
    """Two of the three test metrics are non-vacuum, so `R_ab = 0` cannot be
    what is producing the zero."""
    rows = gb.measure_the_lanczos_tensor_is_correct()["rows"]
    non_vacuum = [r for r in rows if "non-vacuum" in r["metric"]]
    assert len(non_vacuum) >= 2
    for row in non_vacuum:
        assert abs(row["lanczos_kk"]) < 1e-4


def test_the_null_vectors_used_are_actually_null():
    result = gb.measure_the_lanczos_tensor_is_correct()
    assert result["worst_null_residual"] < 1e-12


# ── the sign that decides the branch ────────────────────────────────────────

def test_gauss_bonnet_has_the_same_sign_as_einstein_at_the_neck():
    """The decisive result, and the opposite of what the branch was for."""
    result = gb.measure_gauss_bonnet_reinforces_einstein()
    assert result["reinforces"]
    assert result["ratio_at_neck"] > 0.0
    assert result["ratio_at_neck"] == pytest.approx(
        result["expected_ratio_at_neck"], rel=1e-12)
    assert result["ricci_at_neck"] < 0.0
    assert result["lanczos_at_neck"] < 0.0


def test_the_ratio_is_four_over_the_neck_radius_squared():
    for a in (0.1, 0.3, 0.8):
        g = fm.geometry(1.0, a)
        ratio = float(gb.gauss_bonnet_ratio(np.array([0.0]), g.neck_radius)[0])
        assert ratio == pytest.approx(4.0 / g.neck_radius ** 2, rel=1e-12)


def test_the_ratio_is_the_misner_sharp_parameter_over_f_to_the_fourth():
    """`H_kk/R_kk = 4μ/f⁴` links this round to the seam continuity of P2."""
    g = fm.geometry()
    s = np.linspace(-g.half_length, g.half_length, 101)
    f = fm.throat_radius(s, g.neck_radius)
    expected = 4.0 * fm.misner_sharp(s, g.neck_radius) / f ** 4
    assert np.allclose(gb.gauss_bonnet_ratio(s, g.neck_radius), expected,
                       rtol=1e-12)


def test_the_reinforcement_holds_along_the_whole_throat():
    """Not only at the neck: `μ = b² > 0` everywhere here."""
    result = gb.measure_gauss_bonnet_reinforces_einstein()
    assert result["ratio_is_everywhere_positive"]
    assert result["misner_sharp_is_positive"]


def test_the_closed_forms_match_an_independently_differentiated_riemann():
    """The numerical Lanczos shares no algebra with the closed forms."""
    result = gb.measure_gauss_bonnet_reinforces_einstein()
    assert result["worst_relative_error"] < 1e-2
    for row in result["independent_checks"]:
        assert row["ricci_relative_error"] < 1e-2
        assert row["lanczos_relative_error"] < 1e-2


@pytest.mark.parametrize("lapse", [
    lambda s, b: np.ones_like(s),
    lambda s, b: 1.0 + 0.7 * s,                    # asymmetric: N'(0) != 0
    lambda s, b: 2.0 + np.cos(9.0 * s),
    lambda s, b: np.exp(3.0 * s),
])
def test_the_lapse_drops_out_of_both_contractions_at_the_neck(lapse):
    """`N′` multiplies `f′`, and `f′(0) = 0` is what makes `s = 0` a neck —
    the same structure that made P4 lapse-independent."""
    g = fm.geometry()
    ricci = float(gb.ricci_null_contraction(
        np.array([0.0]), g.neck_radius, lapse=lapse)[0])
    lanczos = float(gb.lanczos_null_contraction_closed_form(
        np.array([0.0]), g.neck_radius, lapse=lapse)[0])
    assert ricci == pytest.approx(-3.0 / g.neck_radius ** 2, rel=1e-6)
    assert lanczos == pytest.approx(-12.0 / g.neck_radius ** 4, rel=1e-6)


def test_the_ricci_contraction_agrees_with_the_source_audit():
    """`R_kk` here must be the same object `source_audit` used for the
    flare-out theorem, or the two rounds are describing different geometries."""
    from geometrodynamics.bulk import source_audit as sa
    g = fm.geometry()
    assert float(gb.ricci_null_contraction(np.array([0.0]), g.neck_radius)[0]) \
        == pytest.approx(sa.required_null_stress(), rel=1e-9)


# ── the coupling ────────────────────────────────────────────────────────────

def test_no_positive_coupling_satisfies_the_matter_nec():
    result = gb.measure_the_required_coupling()
    assert result["no_positive_coupling_works"]
    assert result["heterotic_makes_it_worse"]


def test_the_threshold_is_minus_a_quarter_of_the_neck_radius_squared():
    for a in (0.1, 0.3, 0.8):
        g = fm.geometry(1.0, a)
        assert gb.coupling_threshold(g.neck_radius) == pytest.approx(
            -0.25 * g.neck_radius ** 2, rel=1e-14)
        assert gb.coupling_threshold(g.neck_radius) == pytest.approx(
            -0.25 * math.sin(a) ** 4, rel=1e-13)


def test_the_frozen_numerical_threshold_reproduces():
    """`α_GB ≤ −0.001906728` at `R = 1, a = 0.3`, written down before the code."""
    g = fm.geometry(1.0, 0.3)
    assert gb.coupling_threshold(g.neck_radius) == pytest.approx(
        -0.001906728, abs=1e-9)


def test_the_matter_stress_crosses_zero_exactly_at_the_threshold():
    g = fm.geometry()
    threshold = gb.coupling_threshold(g.neck_radius)
    assert float(gb.matter_null_stress(np.array([0.0]), g.neck_radius,
                                       threshold)[0]) == pytest.approx(0.0, abs=1e-9)
    assert float(gb.matter_null_stress(np.array([0.0]), g.neck_radius,
                                       0.5 * threshold)[0]) < 0.0
    assert float(gb.matter_null_stress(np.array([0.0]), g.neck_radius,
                                       2.0 * threshold)[0]) > 0.0


def test_the_heterotic_sign_is_recorded_as_positive_and_deepens_it():
    """`α_GB = α′/8 > 0`, so the best-motivated value fails hardest."""
    assert gb.HETEROTIC_SIGN == +1
    g = fm.geometry()
    heterotic = float(gb.matter_null_stress(
        np.array([0.0]), g.neck_radius, +0.25 * g.neck_radius ** 2)[0])
    einstein = float(gb.matter_null_stress(np.array([0.0]), g.neck_radius, 0.0)[0])
    assert heterotic < einstein < 0.0


# ── what a negative coupling actually costs ─────────────────────────────────

def test_lovelock_terminates_at_gauss_bonnet_in_five_dimensions():
    """Why the withdrawn G6 was wrong: there is no further tower in `D = 5`.
    The `k`-th density antisymmetrises `2k` indices, so it vanishes for
    `2k > D` and is topological at `2k = D`."""
    assert gb.lovelock_status(2, 4) == "topological"
    assert gb.lovelock_status(2, 5) == "dynamical"
    assert gb.lovelock_status(3, 5) == "identically zero"
    assert gb.lovelock_status(3, 6) == "topological"
    assert gb.lovelock_status(3, 7) == "dynamical"
    assert gb.measure_the_negative_coupling_requirement()[
        "tower_terminates_at_gauss_bonnet"]


def test_the_withdrawn_claim_is_recorded_as_withdrawn():
    note = gb.measure_the_negative_coupling_requirement()["what_was_withdrawn"]
    assert "IDENTICALLY ZERO" in note
    assert "complete classical theory" in note


def test_a_neck_only_coupling_fails_the_nec_elsewhere_on_the_throat():
    """The neck is the *easiest* point: the pointwise requirement
    `alpha <= -f^4/(4 mu)` strengthens outward as `f` grows."""
    result = gb.measure_the_negative_coupling_requirement()
    assert not result["neck_only_satisfies_nec_globally"]
    assert result["neck_only_min_over_throat"] < -1.0


def test_the_global_threshold_is_minus_a_quarter_of_the_bulk_radius_squared():
    """`f_m^4/(4b^2) = R^2/4` exactly, independent of the mouth angle."""
    for R in (1.0, 2.5):
        assert gb.global_coupling_threshold(R) == pytest.approx(
            -0.25 * R ** 2, rel=1e-14)
        for a in (0.1, 0.3, 0.8):
            g = fm.geometry(R, a)
            assert -g.mouth_radius ** 4 / (4.0 * g.neck_radius ** 2) == \
                pytest.approx(-0.25 * R ** 2, rel=1e-12), "a must cancel"
    result = gb.measure_the_negative_coupling_requirement()
    assert result["global_threshold_is_minus_quarter_R_squared"]
    assert result["global_satisfies_nec"]


def test_the_global_requirement_is_a_cosmological_gauss_bonnet_length():
    """`sqrt|alpha| >= R/2` — half the closed universe, not a short scale."""
    result = gb.measure_the_negative_coupling_requirement()
    assert result["length_over_bulk_radius"] == pytest.approx(0.5, rel=1e-12)
    assert result["ratio_global_to_neck"] == pytest.approx(
        result["expected_ratio"], rel=1e-9)


# ── scope ───────────────────────────────────────────────────────────────────

def test_the_branch_is_reported_narrowed_and_not_closed():
    """The corrected verdict. A negative coupling *does* work, so claiming
    closure would be the same overreach the audit was written about."""
    ledger = gb.measure_the_gauss_bonnet_ledger()
    assert ledger["branch_is_narrowed_not_closed"]
    assert "NARROWED, not closed" in ledger["verdict"]
    assert "not closed" in ledger["what_this_narrows"]


def test_the_ledger_keeps_negative_coupling_egb_as_an_open_branch():
    ledger = gb.measure_the_gauss_bonnet_ledger()
    assert len(ledger["what_remains"]) == 5
    branches = " ".join(ledger["what_remains"].keys()).lower()
    assert "negative-coupling egb" in branches, \
        "the branch this round narrowed must stay listed as open"
    detail = " ".join(ledger["what_remains"].values())
    assert "not refuted" in detail
    assert "stability" in detail, "global existence/stability must be flagged open"
    untested = ledger["not_refuted_here"]
    assert "dilatonic" in untested.lower()
    assert "f(R)" in untested
    assert "terminates at Gauss-Bonnet" in untested, \
        "the Lovelock tower is not a separate premise in D = 5"


def test_the_ledger_derives_its_numbers_from_the_measurements():
    ledger = gb.measure_the_gauss_bonnet_ledger()
    sign = gb.measure_gauss_bonnet_reinforces_einstein()
    entries = {e["claim"]: e for e in ledger["entries"]}
    row = next(e for k, e in entries.items() if "supply the negative" in k)
    assert f"{sign['ratio_at_neck']:.4f}" in row["evidence"]
    requirement = gb.measure_the_negative_coupling_requirement()
    row = next(e for k, e in entries.items() if "works at the neck suffices" in k)
    assert f"{requirement['neck_only_min_over_throat']:.1f}" in row["evidence"]


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import gauss_bonnet_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)


# ── the spatial Lanczos block (added after review of PR #278) ────────────────

def test_the_spatial_block_vanishes_for_every_ultrastatic_product():
    """`H^i_j = 0` is the 4D Euler tensor, not a consequence of symmetry.

    Three docstrings in this programme credited the vanishing to the `S^4_R`
    slice being maximally symmetric. It holds for a throat slice, which is not,
    and for a slice with no symmetry at all.
    """
    result = gb.measure_the_spatial_block_vanishes()
    assert result["vanishes_for_every_ultrastatic_case"]
    assert result["worst_ultrastatic_residual"] < 1e-4
    for row in result["rows"]:
        if row["ultrastatic"]:
            assert row["refinement_ratio"] > 8.0, (
                f"{row['metric']}: a discretisation zero must shrink like "
                f"step^2, got ratio {row['refinement_ratio']:.1f}")


def test_the_nonconstant_lapse_control_does_not_vanish():
    """Without this the check could not distinguish a result from a bug."""
    result = gb.measure_the_spatial_block_vanishes()
    assert result["control_does_not_vanish"]
    assert result["control_residual"] > 1e-3
    control = next(r for r in result["rows"] if not r["ultrastatic"])
    assert control["refinement_ratio"] < 1.5, \
        "a structural nonzero must NOT shrink under refinement"


def test_the_reason_given_is_the_euler_tensor_not_maximal_symmetry():
    result = gb.measure_the_spatial_block_vanishes()
    why = result["why"].lower()
    assert "euler" in why and "topological" in why
    assert "maximal symmetry is not used" in why


def test_the_throat_lanczos_time_component_matches_the_closed_form():
    """`H_tt = -12 b^4/f^8` along the scalar-flat handle, checked at two
    points including the neck."""
    neck_radius = 0.7
    for s in (0.0, 0.9):
        point = np.array([0.0, s, 1.0, 0.8, 0.0])
        mixed, _ = gb.lanczos_mixed(gb._throat_metric(neck_radius), point,
                                    step=5e-4)
        f_squared = s * s + neck_radius ** 2
        expected = -12.0 * neck_radius ** 4 / f_squared ** 4
        # g_tt = -1, so H^t_t = -H_tt.
        assert mixed[0, 0] == pytest.approx(-expected, rel=1e-4), \
            f"at s = {s}"


def test_the_negative_egb_module_no_longer_credits_maximal_symmetry():
    from geometrodynamics.bulk import negative_egb as egb
    empty = egb.measure_the_critical_exterior_is_empty()
    why = empty["why_the_pressure_cannot_move"]
    assert "ultrastatic" in why.lower()
    assert "maximally symmetric" not in why.lower()
