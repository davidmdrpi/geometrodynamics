"""
Tests for the closure-phase ledger experiment.

These are intentionally minimal:

- Verify the experiment runs to completion without exceptions, even if
  every repo import fails (i.e., on a fresh clone without the
  geometrodynamics package installed).
- Verify Layer 1 closure: under T² + χ=0, all three lepton ledger
  rows have available_total_mod_2pi ≈ 0.
- Verify Layer 2 produces a non-empty blocker report.
- Verify JSON and markdown writers succeed.

The intent is that these tests pass without needing to install the
full geometrodynamics package — the experiment should fall back to
README-published values cleanly when imports fail. Once Claude Code
verifies the import paths against the real repo, additional integration
tests against actual `geometrodynamics.*` symbols can be added.
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import pytest

from experiments.closure_ledger import run_experiment
from experiments.closure_ledger.ledger import (
    _load_repo_constants,
    compute_lepton_ledger,
)
from experiments.closure_ledger.sk_bridge import (
    DEFAULT_SK_CANDIDATE,
    WIRED_CANDIDATES,
    WKB_CONVENTION,
    Mode,
    phi_convergence_table,
    phi_radial_for_mode,
    phi_radial_from_sk,
    s_k_membership,
)


TAU = 2.0 * math.pi


def test_repo_constants_load_or_fall_back():
    """Constants either import from the repo or fall back cleanly."""
    constants, errors = _load_repo_constants()
    assert constants.action_base > 0
    assert constants.beta_lepton > 0
    assert constants.beta_quark > 0
    assert constants.lepton_quanta == 100
    assert constants.quark_quanta == 466
    # Errors list is fine to be empty or populated; just must be a list.
    assert isinstance(errors, list)


def test_lepton_ledger_universality_under_t2():
    """Layer 1: available terms close mod 2π across e, μ, τ under T²."""
    rows, _constants, _errors = compute_lepton_ledger(
        chi=0.0, transport_power=2,
    )
    assert len(rows) == 3
    for row in rows:
        # Every row should report a partial-or-full closure status that
        # contains the substring "closes" — i.e. the available terms
        # close mod 2π.
        assert "closes" in row.closure_status, (
            f"{row.label} (k={row.k}) status: {row.closure_status}"
        )
        assert math.isclose(
            row.available_total_mod_2pi, 0.0, abs_tol=1e-9,
        ) or math.isclose(
            row.available_total_mod_2pi, TAU, abs_tol=1e-9,
        ), (
            f"{row.label} (k={row.k}) "
            f"available mod 2π = {row.available_total_mod_2pi}"
        )


def test_lepton_ledger_t1_diagnostic_also_universal():
    """T¹ diagnostic: should also be universal across leptons (at 3π/2)."""
    rows, _c, _e = compute_lepton_ledger(chi=0.0, transport_power=1)
    mods = [r.available_total_mod_2pi for r in rows]
    spread = max(mods) - min(mods)
    assert spread < 1e-9, f"T¹ mods = {mods}, spread = {spread}"


def test_layer2_blocker_present_and_structured():
    """Layer 2: blocker report has verdict, evidence, candidates."""
    result = run_experiment()
    blocker = result.sk_bridge_blocker
    assert "verdict" in blocker and blocker["verdict"]
    assert "evidence" in blocker and len(blocker["evidence"]) >= 1
    assert "candidates" in blocker and len(blocker["candidates"]) >= 1
    assert "next_steps" in blocker and len(blocker["next_steps"]) >= 1
    # Each candidate must have the expected fields.
    for cand in blocker["candidates"]:
        for key in (
            "name", "formula", "physical_picture",
            "advantages", "open_questions",
        ):
            assert key in cand and cand[key], (
                f"candidate {cand.get('name')} missing field {key}"
            )


def test_full_experiment_runs_and_serializes():
    """End-to-end: run experiment, write JSON + markdown, reload JSON."""
    result = run_experiment()
    assert result.overall_status
    with tempfile.TemporaryDirectory() as td:
        out = Path(td)
        result.write_json(out / "result.json")
        result.write_summary_markdown(out / "summary.md")
        assert (out / "result.json").exists()
        assert (out / "summary.md").exists()
        loaded = json.loads((out / "result.json").read_text())
        # Default run wires the S(k) bridge → Layer 2.
        assert loaded["experiment_name"] == "closure_ledger.layer2"
        assert loaded["sk_candidate"] == DEFAULT_SK_CANDIDATE
        assert "rows" in loaded and len(loaded["rows"]) == 3
        assert "sk_bridge_blocker" in loaded
        # Markdown is non-empty and contains expected headers.
        md = (out / "summary.md").read_text()
        assert "# Closure-phase ledger — run summary" in md
        assert "## Per-lepton ledger" in md
        assert "## Layer 2 blocker" in md
        assert "## Radial bulk channel — per-mode breakdown" in md


def test_layer1_only_run_reproduces_pre_bridge_universality():
    """Passing sk_candidate='none' reproduces the Layer-1 universality result."""
    result = run_experiment(sk_candidate="none")
    assert result.experiment_name == "closure_ledger.layer1"
    assert result.universality_check["universal"] is True
    assert math.isclose(
        result.universality_check["universal_value"], 0.0, abs_tol=1e-9,
    )


def test_sk_membership_candidate_a():
    """Candidate A: S(k) is the odd-l ground states up to l=k."""
    assert s_k_membership(1, "A_lowest_radial_per_l") == [Mode(l=1, n=0)]
    assert s_k_membership(3, "A_lowest_radial_per_l") == [
        Mode(l=1, n=0), Mode(l=3, n=0),
    ]
    assert s_k_membership(5, "A_lowest_radial_per_l") == [
        Mode(l=1, n=0), Mode(l=3, n=0), Mode(l=5, n=0),
    ]


def test_sk_membership_candidate_b1():
    """Candidate B1: single l=k angular mode at n=0."""
    assert s_k_membership(1, "B1_single_angular_mode") == [Mode(l=1, n=0)]
    assert s_k_membership(3, "B1_single_angular_mode") == [Mode(l=3, n=0)]
    assert s_k_membership(5, "B1_single_angular_mode") == [Mode(l=5, n=0)]


def test_sk_membership_candidate_b2():
    """Candidate B2: single l=1 radial excitation at n = (k-1)/2."""
    assert s_k_membership(1, "B2_single_radial_excitation") == [Mode(l=1, n=0)]
    assert s_k_membership(3, "B2_single_radial_excitation") == [Mode(l=1, n=1)]
    assert s_k_membership(5, "B2_single_radial_excitation") == [Mode(l=1, n=2)]


def test_sk_membership_candidate_c_is_not_implemented():
    """Candidate C is intentionally NotImplementedError until defined."""
    with pytest.raises(NotImplementedError):
        s_k_membership(1, "C_eigenvector_weighted")


def test_sk_membership_none_is_empty():
    """sk_candidate='none' yields an empty membership."""
    assert s_k_membership(1, "none") == []
    assert s_k_membership(5, "none") == []


def test_phi_radial_for_mode_l1_n0_is_finite_and_positive():
    """The Bohr-Sommerfeld radial action for the (l=1, n=0) mode resolves."""
    result = phi_radial_for_mode(l=1, n=0)
    assert result.status == "computed"
    assert result.phi is not None and math.isfinite(result.phi)
    assert result.phi > 0
    # Eigenfrequency for the lowest l=1 mode should be near 1 in geometric units.
    assert result.omega is not None and 0.5 < result.omega < 2.0


def test_phi_radial_from_sk_candidate_a_total_is_sum_of_modes():
    """Φ_radial(k) is the sum of Φ(l, n) over S(k)."""
    result = phi_radial_from_sk(5, "A_lowest_radial_per_l")
    assert result.status == "computed"
    assert result.total_phi is not None
    expected = sum(m.phi for m in result.modes)
    assert math.isclose(result.total_phi, expected, abs_tol=1e-12)
    assert len(result.modes) == 3   # l ∈ {1, 3, 5}


def test_radial_channel_wired_in_default_run():
    """Default run has no row blocked on radial_bulk_phase."""
    result = run_experiment()
    for row in result.rows:
        names = [t["name"] for t in row["terms"]]
        assert "radial_bulk_phase" in names
        radial = next(t for t in row["terms"] if t["name"] == "radial_bulk_phase")
        assert radial["status"] == "available"
        assert radial["value"] is not None
        assert "radial_bulk_phase" not in row["blocking_terms"]


def test_layer2_blocker_marks_implemented_candidate():
    """When candidate A is wired, blocker reports it as the implemented one."""
    result = run_experiment(sk_candidate="A_lowest_radial_per_l")
    blocker = result.sk_bridge_blocker
    assert blocker["implemented_candidate"] == "A_lowest_radial_per_l"
    impls = {
        c["name"]: c["implementation_status"] for c in blocker["candidates"]
    }
    assert impls["A_lowest_radial_per_l"] == "implemented"
    assert impls["B1_single_angular_mode"] == "open"
    assert impls["B2_single_radial_excitation"] == "open"
    assert impls["C_eigenvector_weighted"] == "open"


def test_phi_radial_uses_wkb_convention_label():
    """ModePhase carries the WKB convention label by default."""
    mp = phi_radial_for_mode(l=1, n=0)
    assert mp.convention == WKB_CONVENTION
    assert mp.maslov_correction == 0.0


def test_phi_grid_convergence_for_l1_n0():
    """Φ(l=1, n=0) is grid-stable across N: variation under ~1e-3."""
    table = phi_convergence_table(1, 0, Ns=(60, 80, 100, 120))
    phis = [row["phi"] for row in table]
    assert all(p is not None for p in phis)
    # All values within 1e-3 of the largest-grid value.
    ref = phis[-1]
    assert max(abs(p - ref) for p in phis) < 1e-3, (
        f"Φ(l=1, n=0) grid spread = {phis} (ref={ref})"
    )


def test_maslov_shift_is_additive():
    """A maslov_correction adds a constant to phi at fixed (l, n, N)."""
    base = phi_radial_for_mode(l=1, n=0)
    shifted = phi_radial_for_mode(l=1, n=0, maslov_correction=math.pi / 4.0)
    assert math.isclose(
        shifted.phi - base.phi, math.pi / 4.0, abs_tol=1e-12,
    )
    assert shifted.maslov_correction == math.pi / 4.0


def test_b1_b2_falsify_universality_under_wkb():
    """Both B1 and B2 break universal closure mod 2π under WKB convention."""
    for cand in ("B1_single_angular_mode", "B2_single_radial_excitation"):
        result = run_experiment(sk_candidate=cand)
        assert result.universality_check["universal"] is False, (
            f"{cand} unexpectedly universal: {result.universality_check}"
        )
        # Spread should be visibly larger than numeric tolerance.
        assert result.universality_check["spread"] > 1e-3


def test_run_comparison_covers_layer1_baseline_and_wired_candidates():
    """run_comparison runs Layer-1 + every wired candidate and labels each."""
    from experiments.closure_ledger.runner import run_comparison

    comparison = run_comparison()
    candidates_run = comparison["candidates_run"]
    assert "none" in candidates_run
    for cand in WIRED_CANDIDATES:
        assert cand in candidates_run

    by_cand = {row["candidate"]: row for row in comparison["status_table"]}
    assert by_cand["none"]["result"] == "PASS"
    for cand in WIRED_CANDIDATES:
        assert by_cand[cand]["result"] == "FAIL"


def test_quark_sector_quanta_gap_is_366():
    """The quark/lepton lock-quanta gap is the 366 the chat-pass derived."""
    result = run_experiment()
    q = result.quark_sector
    assert q["lepton_lock_quanta"] == 100
    assert q["quark_lock_quanta"] == 466
    assert q["lock_quanta_gap"] == 366


def test_p3_is_downgraded():
    """The 366-quanta prediction is recorded as downgraded, not asserted."""
    result = run_experiment()
    blocker = result.sk_bridge_blocker
    assert "downgraded_predictions" in blocker
    assert any(
        "366" in dp for dp in blocker["downgraded_predictions"]
    ), f"downgraded_predictions = {blocker['downgraded_predictions']}"
