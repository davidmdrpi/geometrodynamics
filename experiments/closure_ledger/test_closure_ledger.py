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
        assert loaded["experiment_name"] == "closure_ledger.layer1"
        assert "rows" in loaded and len(loaded["rows"]) == 3
        assert "sk_bridge_blocker" in loaded
        # Markdown is non-empty and contains expected headers.
        md = (out / "summary.md").read_text()
        assert "# Closure-phase ledger — run summary" in md
        assert "## Per-lepton ledger" in md
        assert "## Layer 2 blocker" in md


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
