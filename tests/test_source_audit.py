"""Can any classical field already in BAM supply the neck's `T_kk < 0`?

Checked against `docs/source_audit_prereg.md`, frozen before the module.

Four kinds. The *theorem* ones pin that the neck price follows from
Raychaudhuri alone, so it is not an artefact of `N = 1` or of symmetry — and
that this is a **known** result, recovered rather than discovered. The
*candidate* ones build each null stress from its actual tensor. The *loophole*
one pins the single coupling that could have flipped the sign, and why it does
not in any dimension. And the *dynamic* one attacks the Raychaudhuri no-go by
looking for a turning point rather than confirming its absence.
"""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import finite_mouth as fm
from geometrodynamics.bulk import source_audit as sa


# ── the theorem ─────────────────────────────────────────────────────────────

def test_the_neck_price_follows_from_raychaudhuri_alone():
    """No reference to `p_s`, `N = 1`, or reflection symmetry."""
    result = sa.measure_the_flare_out_requirement()
    assert result["identity_holds"]
    assert result["raychaudhuri_residual"] < 1e-9
    assert result["theta_at_neck"] == pytest.approx(0.0, abs=1e-12)
    assert result["dtheta_at_neck"] == pytest.approx(
        result["expected_dtheta"], rel=1e-9)
    assert result["flare_out_is_positive"]


def test_the_two_derivations_of_the_neck_price_agree():
    """The component route (`ρ+p_s`) and the Raychaudhuri route must land on
    the same number, having shared no algebra."""
    g = fm.geometry()
    component = fm.null_energy_at_neck(g.neck_radius)
    geometric = sa.measure_the_flare_out_requirement()["ricci_kk_at_neck"]
    assert geometric == pytest.approx(component, rel=1e-9)
    assert sa.required_null_stress() == pytest.approx(component, rel=1e-12)


def test_the_screen_dimension_is_three_in_five_dimensions():
    """The `θ²/3` coefficient is `1/(D−2)`; getting it wrong would move the
    neck price."""
    assert sa.measure_the_flare_out_requirement()["screen_dimension"] == 3


def test_the_theorem_is_attributed_and_not_claimed_as_new():
    """The repo has ~1 external anchor for 324 claims; miscrediting a known
    theorem would waste a second one."""
    note = sa.measure_the_flare_out_requirement()["attribution"]
    assert "NOT a new theorem" in note
    assert "Morris-Thorne" in note and "Hochberg-Visser" in note
    assert "VALIDATES" in note


# ── the candidates ──────────────────────────────────────────────────────────

def test_every_existing_bam_field_gives_non_negative_null_stress():
    result = sa.measure_every_candidate_stress()
    assert result["all_non_negative"]
    assert result["candidates_with_negative_null_stress"] == []
    assert len(result["rows"]) >= 9


def test_the_minimal_scalar_null_stress_is_a_square():
    values = sa.null_stress_minimal_scalar(np.array([-3.0, 0.0, 2.5]))
    assert np.all(values >= 0.0)
    assert values[1] == 0.0


def test_the_potential_drops_out_of_the_null_contraction():
    """Including a negative, symmetry-breaking region: `V` is pure `g_ab`."""
    metric = np.diag([-1.0, 1.0, 1.0, 1.0, 1.0])
    null_vector = np.array([1.0, 1.0, 0.0, 0.0, 0.0])
    assert abs(null_vector @ metric @ null_vector) < 1e-15
    covector = np.array([0.4, -1.1, 0.7, 0.0, 0.2])
    kinetic = float(covector @ np.linalg.inv(metric) @ covector)
    reference = None
    for potential in (-9.0, -1.0, 0.0, 4.0):
        stress = (np.outer(covector, covector)
                  - 0.5 * metric * (kinetic + 2.0 * potential))
        contracted = float(null_vector @ stress @ null_vector)
        assert contracted == pytest.approx(
            float(covector @ null_vector) ** 2, abs=1e-12)
        if reference is None:
            reference = contracted
        assert contracted == pytest.approx(reference, abs=1e-12)


def test_the_throat_order_field_cannot_pay_its_own_bill():
    """C2 is the field introduced to *represent* the throat."""
    values = sa.null_stress_complex_order_field(
        np.array([1.0 + 2.0j, -0.5j, 0.0]), stiffness=1.7)
    assert np.all(values >= 0.0)
    note = sa.measure_every_candidate_stress()[
        "the_order_field_cannot_pay_its_own_bill"]
    assert "cannot support the object it was introduced to represent" in note


def test_the_maxwell_null_stress_is_non_negative_and_orthogonal():
    """`V_a = F_ab k^b` satisfies `V·k = 0`, so `V` is spacelike or null."""
    rng = np.random.default_rng(7)
    metric = np.diag([-1.0, 1.0, 2.0, 2.0, 2.0])
    null_vector = np.array([math.sqrt(1.0), 1.0, 0.0, 0.0, 0.0])
    assert abs(null_vector @ metric @ null_vector) < 1e-15
    for _ in range(50):
        raw = rng.normal(size=(5, 5))
        strength = raw - raw.T
        assert abs(float((strength @ null_vector) @ null_vector)) < 1e-12
        assert sa.null_stress_maxwell(strength, metric, null_vector) >= -1e-12


def test_the_perfect_fluid_is_non_negative_exactly_when_the_nec_holds():
    assert np.all(sa.null_stress_perfect_fluid(
        np.array([1.0, 2.0]), np.array([0.0, -0.5]), np.array([1.0, 1.0])) >= 0)
    violating = sa.null_stress_perfect_fluid(1.0, -2.0, 1.0)
    assert violating < 0.0, "the fluid row is conditional on the NEC, not free"


def test_the_projected_weyl_row_says_it_is_the_wrong_equation():
    rows = {r["id"]: r for r in sa.measure_every_candidate_stress()["rows"]}
    assert "5D Einstein equation" in rows["C9"]["reason"]
    assert "T^(5)_ab = 0" in rows["C9"]["reason"]


def test_the_spinor_sector_is_recorded_as_transport_not_a_source():
    rows = {r["id"]: r for r in sa.measure_every_candidate_stress()["rows"]}
    assert "no stress" in rows["C8"]["reason"]
    assert "transport object" in rows["C8"]["reason"]


# ── the loophole ────────────────────────────────────────────────────────────

def test_the_conformal_coupling_never_flips_the_sign_in_any_dimension():
    """`1 − 2ξ_c = D/(2(D−1)) > 0` for every `D ≥ 2`."""
    result = sa.measure_the_nonminimal_loophole()
    assert result["conformal_never_flips"]
    for row in result["dimensions"]:
        assert row["one_minus_two_xi"] > 0.0
        assert row["closed_form_matches"], "the closed form must be exact"


def test_the_five_dimensional_conformal_coupling_is_three_sixteenths():
    assert sa.conformal_coupling(5) == pytest.approx(3.0 / 16.0, rel=1e-15)
    assert sa.CONFORMAL_COUPLING_5D == pytest.approx(3.0 / 16.0, rel=1e-15)
    assert 1.0 - 2.0 * sa.conformal_coupling(5) == pytest.approx(0.625, rel=1e-15)
    assert sa.one_minus_two_xi_conformal(5) == pytest.approx(0.625, rel=1e-15)


def test_a_sign_flip_needs_a_coupling_beyond_one_half():
    assert sa.sign_flip_threshold() == 0.5
    grad = np.array([1.5])
    assert sa.null_stress_nonminimal_at_node(grad, 0.0)[0] > 0.0
    assert sa.null_stress_nonminimal_at_node(grad, 3.0 / 16.0)[0] > 0.0
    assert sa.null_stress_nonminimal_at_node(grad, 0.5)[0] == pytest.approx(0.0)
    assert sa.null_stress_nonminimal_at_node(grad, 0.8)[0] < 0.0


def test_the_node_is_why_the_prefactor_is_one():
    note = sa.measure_the_nonminimal_loophole()["why_the_node_matters"]
    assert "prefactor" in note and "exactly 1" in note
    assert "defect core AT the node" in note


# ── the dynamic escape, attacked ────────────────────────────────────────────

def test_no_non_negative_ricci_profile_produces_a_turning_point():
    """The falsifier: a single `θ: − → +` would mean a sign error upstream."""
    result = sa.measure_the_dynamic_escape_fails()
    assert not result["any_turning_point"]
    assert result["turning_cases"] == []
    assert all(not row["turned_positive"] for row in result["rows"])


def test_every_ray_focuses_within_its_analytic_bound():
    """`λ_caustic ≤ (D−2)/|θ₀| = 3/|θ₀|`, and non-negative `R_kk` shortens it."""
    result = sa.measure_the_dynamic_escape_fails()
    assert result["all_focused_within_the_analytic_bound"]
    for row in result["rows"]:
        assert row["measured_caustic"] is not None
        assert row["measured_caustic"] <= row["focusing_bound"] * 1.05


def test_the_vacuum_caustic_matches_the_closed_form():
    """In vacuum `dθ/dλ = −θ²/3` integrates to a caustic at exactly `3/|θ₀|`."""
    rows = {(r["profile"], r["theta_initial"]): r
            for r in sa.measure_the_dynamic_escape_fails()["rows"]}
    for theta0 in (-0.05, -0.5, -2.0):
        row = rows[("vacuum R_kk = 0", theta0)]
        assert row["measured_caustic"] == pytest.approx(
            3.0 / abs(theta0), rel=0.02)


def test_the_earlier_fixed_window_is_recorded_as_an_artefact():
    """A `λ ≤ 6` window missed the `θ₀ = −0.05` caustic at `λ = 60`; the
    module has to say that was the window, not the theorem."""
    note = sa.measure_the_dynamic_escape_fails()["why_the_window_matters"]
    assert "artefact of the window" in note
    assert "never depended on the window" in note


# ── the verdict ─────────────────────────────────────────────────────────────

def test_the_verdict_is_negative_and_says_what_it_is_about():
    ledger = sa.measure_the_source_audit_ledger()
    assert ledger["verdict_is_no"]
    assert "cannot support a traversable particle throat" in ledger["verdict"]
    assert "CURRENT field content" in ledger["what_this_forces"]


def test_the_ledger_names_five_open_branches_and_refutes_none():
    ledger = sa.measure_the_source_audit_ledger()
    branches = ledger["the_remaining_branches"]
    assert len(branches) == 5
    joined = " ".join(branches.values())
    assert "Gauss-Bonnet" in joined
    assert "ghost" in joined or "wrong-sign" in joined
    assert "Casimir" in joined
    assert "refutes none" in ledger["what_the_audit_does_not_do"]


def test_the_ledger_derives_its_numbers_from_the_measurements():
    ledger = sa.measure_the_source_audit_ledger()
    requirement = sa.measure_the_flare_out_requirement()
    assert ledger["required_null_stress"] == pytest.approx(
        requirement["required_null_stress"], rel=1e-12)
    entries = {e["claim"]: e for e in ledger["entries"]}
    row = next(e for k, e in entries.items() if "artefact of N = 1" in k)
    assert f"{requirement['raychaudhuri_residual']:.1e}" in row["evidence"]


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import source_audit_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
