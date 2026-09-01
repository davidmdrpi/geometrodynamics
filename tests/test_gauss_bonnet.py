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


# ── the regime of validity ──────────────────────────────────────────────────

def test_at_the_threshold_the_correction_equals_the_leading_term():
    result = gb.measure_the_expansion_breaks_down()
    assert result["it_equals_the_leading_term"]
    assert result["relative_size_at_threshold"] == pytest.approx(-1.0, abs=1e-12)


def test_the_required_gauss_bonnet_length_is_half_the_neck_radius():
    result = gb.measure_the_expansion_breaks_down()
    assert result["length_ratio"] == pytest.approx(0.5, rel=1e-12)


# ── scope ───────────────────────────────────────────────────────────────────

def test_the_branch_is_reported_closed():
    ledger = gb.measure_the_gauss_bonnet_ledger()
    assert ledger["branch_is_closed"]
    assert "does NOT reopen" in ledger["verdict"]


def test_the_ledger_leaves_four_branches_and_names_what_is_untested():
    ledger = gb.measure_the_gauss_bonnet_ledger()
    assert len(ledger["what_remains"]) == 4
    untested = ledger["not_refuted_here"]
    assert "Dilatonic" in untested
    assert "Lovelock" in untested
    assert "f(R)" in untested


def test_the_ledger_derives_its_numbers_from_the_measurements():
    ledger = gb.measure_the_gauss_bonnet_ledger()
    sign = gb.measure_gauss_bonnet_reinforces_einstein()
    entries = {e["claim"]: e for e in ledger["entries"]}
    row = next(e for k, e in entries.items() if "supply the negative" in k)
    assert f"{sign['ratio_at_neck']:.4f}" in row["evidence"]


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import gauss_bonnet_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
