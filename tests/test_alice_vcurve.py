"""V-curve, aperture scaling and envelope control for the Alice-handle exchange."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

LEDGER = Path(__file__).resolve().parents[1] / "experiments" / "closure_ledger"
RUN = LEDGER / "runs" / "20261125_alice_vcurve"


def _module(name):
    spec = importlib.util.spec_from_file_location(name, LEDGER / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PRED = _module("alice_vcurve_predictions")


def test_predictions_module_calls_no_solver():
    assert "alice_handle_ports" not in (LEDGER / "alice_vcurve_predictions.py").read_text()


def test_grid_sits_on_sidelobe_peaks_and_nulls_are_sin2_zeros():
    for tau in PRED.GRID:
        assert PRED.spectrum("sin2", tau, 2 * PRED.OMEGA) > 1e-4 * (8.0 / tau) ** 3
    for tau in (k * np.pi / PRED.OMEGA for k in PRED.NULL_KS):
        assert PRED.spectrum("sin2", tau, 2 * PRED.OMEGA) < 1e-9


def test_frozen_prediction_table():
    frozen = json.loads((RUN / "frozen_predictions.json").read_text())
    for row in frozen["vcurve"][::3] + frozen["envelope"]:
        again = PRED.predict(row["a"], row["m"], row["tau"], row["profile"], row["eta"])
        assert again["total"] == pytest.approx(row["total"], rel=1e-6)


def test_in_sample_fit_quality_is_as_stated():
    points = json.loads((RUN / "exploration" / "fit_points.json").read_text())
    for eta, lo, hi in ((0.5, 0.77, 1.91), (0.9, 0.89, 1.48)):
        ratios = []
        for p in points:
            if p["eta"] != eta or p["profile"] != "sin2" or abs(np.log2(p["tau"]) % 1) > 1e-9:
                continue
            pred = PRED.A[eta] * p["F"] + PRED.C_X[eta] * p["T"]
            ratios.append(p["max_abs_defect"] / pred)
        assert lo - 0.03 <= min(ratios) and max(ratios) <= hi + 0.03
