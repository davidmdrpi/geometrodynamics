"""Alice-handle exchange under detector depletion: two-level predictions, probe and archive."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

LEDGER = Path(__file__).resolve().parents[1] / "experiments" / "closure_ledger"


def _module(name):
    spec = importlib.util.spec_from_file_location(name, LEDGER / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PRED = _module("alice_exchange_predictions")

# The two-level table of docs/alice_exchange_prereg.md, at eta = 0.5.
FROZEN_HALF = {(0.5, 1): 1.068e-3, (0.4, 1): 3.459e-4, (0.3, 1): 9.708e-5, (0.3, 2): 2.828e-7}


def test_predictions_module_calls_no_solver():
    assert "alice_handle_ports" not in (LEDGER / "alice_exchange_predictions.py").read_text()


@pytest.mark.parametrize("cell", sorted(FROZEN_HALF))
def test_two_level_table(cell):
    split = PRED.NECK.prediction(*cell)["predicted_delta"]
    assert PRED.defect(0.5, split)["max_abs_defect"] == pytest.approx(FROZEN_HALF[cell], rel=2e-3)


def test_free_exchange_is_exact_complementarity():
    split = 0.02
    basis = np.array([[1, 1j, 0, 0], [0, 0, 1, 1j]], dtype=complex)
    swapped = PRED._free(basis, split, np.pi / split)
    np.testing.assert_allclose(np.abs(swapped[0]), np.abs(basis[1]), atol=1e-12)


def test_calibration_hits_the_requested_depletion():
    split = PRED.NECK.prediction(0.4, 1)["predicted_delta"]
    for eta in PRED.ETAS:
        g0 = PRED.calibrate(eta, split)
        assert PRED.absorbed_pure_plus(g0, split) == pytest.approx(eta, abs=1e-9)


def test_defect_vanishes_with_the_exchange_during_the_window():
    values = [PRED.defect(0.5, s)["max_abs_defect"] for s in (1e-2, 1e-3, 1e-4)]
    assert values[0] > values[1] > values[2]
    assert values[2] < 1e-6
