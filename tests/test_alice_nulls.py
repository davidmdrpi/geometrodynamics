"""Residual at the spectral nulls of the Alice-handle exchange: frozen predictions and in-sample record."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

LEDGER = Path(__file__).resolve().parents[1] / "experiments" / "closure_ledger"
RUN = LEDGER / "runs" / "20261126_alice_nulls"
_spec = importlib.util.spec_from_file_location("alice_nulls_predictions", LEDGER / "alice_nulls_predictions.py")
PRED = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(PRED)


def test_predictions_module_calls_no_solver():
    assert "alice_handle_ports" not in (LEDGER / "alice_nulls_predictions.py").read_text()


def test_held_out_cells_are_new():
    tested = {0.3, 0.35, 0.5}
    assert not set(PRED.CELLS) & tested


def test_frozen_table_reproduces():
    frozen = json.loads((RUN / "frozen_predictions.json").read_text())
    for row in frozen["nulls"][::3]:
        again = PRED.terms(row["a"], row["m"], row["tau"], row["profile"], row["eta"])
        assert again["total"] == pytest.approx(row["total"], rel=1e-6)
        assert again["switching"] < 1e-12          # exact sin^2 nulls


def test_radial_truncation_record():
    data = json.loads((RUN / "exploration" / "radial_truncation_0.5_2.json").read_text())
    for row in data.values():
        assert max(row.values()) / min(row.values()) < 1.01


def test_in_sample_fit_quality():
    mis = json.loads((RUN / "exploration" / "mismatch_exact_nulls.json").read_text())
    intr = json.loads((RUN / "exploration" / "intrinsic_exact_nulls.json").read_text())
    ratios = [r[f"full_{e}"] / PRED.terms(0.5, 2, r["tau"], "sin2", e)["total"] for r in mis for e in PRED.ETAS]
    ratios += [r["actual_doublet_eta05"] / PRED.terms(0.3, 2, r["tau"], "sin2", 0.5)["total"] for r in intr if r["a"] == 0.3]
    assert 0.84 <= min(ratios) and max(ratios) <= 1.10


def test_boundary_is_aperture_independent_in_sample():
    rows = json.loads((RUN / "exploration" / "depletion_boundary.json").read_text())
    for tau, eta_max in PRED.BOUNDARY_IN_SAMPLE.items():
        found = [r["eta_max"] for r in rows if r["profile"] == "sin2" and abs(r["tau"] - tau) < 1e-9]
        assert len(found) == 2 and max(abs(np.array(found) - eta_max)) < 2e-3
