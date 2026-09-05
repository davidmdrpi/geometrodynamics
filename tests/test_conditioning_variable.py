"""Round 8 — the conditioning-variable correction, pre-registered in
``docs/conditioning_variable_prereg.md`` (`39be3e3`) before the module."""

import math

import numpy as np

from geometrodynamics.bulk.conditioning_variable import (
    n_window_is_sector_blind, gradient_magnitudes_on_closure,
    two_conditioning_families, phase_axiom_citation,
    downstream_numbers_unchanged, velocity_field_non_uniqueness,
    kappa_is_route_local, dependency_ledger, verdict)


# ── the oracle: why the two windows differ ──────────────────────────────────

def test_an_N_window_cannot_distinguish_the_outcome_sectors():
    """O1/O4 — u x w = -s_A s_B (a x b), so |N| is sector-independent."""
    s = n_window_is_sector_blind()
    assert s["counts_all_equal"], s["counts"]
    assert s["cross_product_is_pm_q_residual"] < 1e-12
    assert s["E_is_zero"]


def test_grad_N_is_constant_on_closure_but_grad_theta_is_not():
    """O2/O3 — the coarea density is 1/|grad(conditioning variable)|."""
    g = gradient_magnitudes_on_closure()
    assert g["grad_N_is_constant_on_closure"]
    assert g["theta_density_varies"]
    for row in g["rows"]:
        assert abs(row["|grad N|"] - math.sin(1.0)) < 1e-12
        lo, hi = row["theta_density_range"]
        assert hi - lo > 1e-3


def test_same_support_two_limits():
    """Borel-Kolmogorov: the zero set does not determine the conditional."""
    t = two_conditioning_families(n=500000)
    assert t["same_support_different_limits"]
    assert t["N_window_is_zero_at_every_eps"]
    assert t["phase_window_tracks_the_repository_law"]
    assert abs(t["phase_window_limit_closed_form"] - 0.3984966504) < 1e-8


def test_the_phase_axiom_exists_and_is_cited_from_source():
    """Rule 1 — the justification must already be in the repository."""
    a = phase_axiom_citation()
    assert a["axiom_is_stated_on_phase"]
    assert "Phase closure" in a["text"]
    assert a["matches"][0][0] < 40


# ── pre-registered rule 3 ───────────────────────────────────────────────────

def test_no_round_5_to_7_number_moves():
    """This corrects the status of an input, not any computation."""
    d = downstream_numbers_unchanged()
    assert d["nothing_downstream_moved"], d["checks"]
    assert d["worst_delta"] < 1e-6
    assert len(d["checks"]) == 5


# ── C: the velocity field is not unique ─────────────────────────────────────

def test_a_divergence_free_addition_leaves_continuity_unchanged():
    v = velocity_field_non_uniqueness()
    assert v["K_is_divergence_free"]
    assert v["continuity_unchanged"]
    assert v["uniqueness_claim_at_line_73_is_false"]


def test_the_mean_velocity_check_cannot_exclude_it_either():
    """int K d^3x = 0 for compactly supported curl A."""
    assert velocity_field_non_uniqueness()["mean_velocity_check_cannot_exclude_it"]


# ── B and the report ────────────────────────────────────────────────────────

def test_kappa_is_route_local_not_a_fourth_universal_input():
    k = kappa_is_route_local()
    assert len(k["three_universal_underived_inputs"]) == 3
    assert "holonomy-trace route" in k["kappa"]


def test_the_ledger_marks_the_conditioning_variable_chosen():
    led = {row["input"]: row["status"] for row in dependency_ledger()}
    assert led["the conditioning VARIABLE (phase, not N)"] == "chosen"
    # the density IS derived once the variable is fixed
    assert led["coarea density |D|/(2|u x v|) GIVEN the phase variable"] == "derived"
    assert led["equal prior on the four outcome sectors"] == "chosen"


def test_the_verdict_is_the_pre_registered_label():
    v = verdict()
    assert v["A_conditioning_variable"] == (
        "CONDITIONING_VARIABLE_IS_A_CHOSEN_INPUT_JUSTIFIED_BY_THE_PHASE_AXIOM")
    assert v["downstream_numbers_unchanged"]


def test_the_module_no_longer_asserts_that_the_set_supplies_the_measure():
    """The module docstring must not state the withdrawn claim at all."""
    import pathlib
    text = pathlib.Path(
        "geometrodynamics/bulk/closure_measurement.py").read_text()
    assert "Haar conditioned on ``N = 0`` is the coarea measure" not in text
    # and it must say which variable it actually conditions on
    assert "Haar conditioned **on the phase**" in text
    assert "conditioning *variable* is an input" in text


def test_the_ledger_entry_moved_from_derived_to_chosen():
    """The load-bearing correction: the round-5 ledger line."""
    import pathlib
    text = pathlib.Path("docs/closure_measurement_dependence.md").read_text()
    assert "coarea conditioning [derived; window limit]" not in text
    assert "conditioning variable = phase, not N [chosen" in text
    assert "coarea density given that variable [derived" in text


def test_the_prose_files_carry_the_correction_not_just_the_old_claim():
    """README and the round-5 write-up may quote the withdrawn sentence, but
    each must also state that the conditioning variable is chosen."""
    import pathlib
    for path in ("README.md", "docs/closure_measurement_dependence.md"):
        text = pathlib.Path(path).read_text()
        assert "conditioning" in text and "variable" in text, path
        assert ("chosen" in text and "Borel" in text), path


def test_the_velocity_uniqueness_claim_is_withdrawn():
    import pathlib
    text = pathlib.Path("docs/born_rule_equivariance.md").read_text()
    assert "it is the unique velocity field whose current closes" not in text
    assert "does not say *unique*" in text
