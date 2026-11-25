"""Alice-handle exchange under detector depletion: two-level predictions, probe and archive."""
import importlib.util
import json
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
    # linear in the splitting: each decade lowers the defect about tenfold
    assert values[0] > values[1] > values[2]
    assert 8 < values[0] / values[1] < 12 and 8 < values[1] / values[2] < 12
    assert values[2] < 3e-6


ARCHIVE = LEDGER / "runs" / "20261125_alice_exchange" / "probe.json"
needs_archive = pytest.mark.skipif(not ARCHIVE.exists(), reason="exchange run not yet archived")


def _defect(result):
    """Largest |u + v - 1| recomputed from the stored reading forms."""
    U, V = np.array(result["U"]), np.array(result["V"])
    Q = U + V - np.eye(4)
    return float(np.max(np.abs(np.linalg.eigvalsh(0.5 * (Q + Q.T)))))


def _verdict(rep):
    """The frozen gates, recomputed from the stored forms."""
    cells = rep["cells"]
    by = {(c["a"], c["m"], c["N"], c["K"]): c for c in cells}
    fails = {}
    for n in (2000, 4000):
        rows = [c for c in cells if c["N"] == n and c["K"] == 150]
        g1 = all(_defect(c["selective"][str(e)]) <= 0.05 * e for c in rows for e in PRED.ETAS)
        small = by[(0.3, 2, n, 150)]
        g2 = all(_defect(small["selective"][str(e)]) <= 1e-4 for e in PRED.ETAS)
        g3 = all(_defect(c["control"]) >= 0.25 for c in rows)
        fails[n] = dict(G1=not g1, G2=not g2, G3=not g3)
    k_ok = n_ok = True
    for a, m in PRED.CELLS:
        for n in (2000, 4000):
            lo, hi = by[(a, m, n, 150)], by[(a, m, n, 300)]
            pairs = [(lo["selective"][str(e)], hi["selective"][str(e)]) for e in PRED.ETAS] + [(lo["control"], hi["control"])]
            k_ok &= all(abs(_defect(x) - _defect(y)) <= max(1e-5, 0.1 * _defect(y)) for x, y in pairs)
        for e in PRED.ETAS:
            x, y = _defect(by[(a, m, 2000, 150)]["selective"][str(e)]), _defect(by[(a, m, 4000, 150)]["selective"][str(e)])
            n_ok &= abs(x - y) <= max(1e-4, 0.2 * y)
    undetected = max(_defect(c["undetected"]) for c in cells)
    balance = max(max([c["control"]["balance"]] + [r["balance"] for r in c["selective"].values()]) for c in cells)
    g4 = undetected <= 1e-8 and balance <= 1e-8 and k_ok and n_ok
    if any(fails[2000][k] and fails[4000][k] for k in ("G1", "G2", "G3")):
        return "COUNTEREXAMPLE"
    if not any(fails[2000].values()) and not any(fails[4000].values()) and g4:
        return "CONFIRMED"
    return "UNRESOLVED_NUMERICALLY"


@needs_archive
def test_archive_verdict_follows_from_the_stored_forms():
    rep = json.loads(ARCHIVE.read_text())
    assert rep["freeze"] == "0fb3b6e"
    assert {(c["a"], c["m"]) for c in rep["cells"]} == set(PRED.CELLS)
    assert len(rep["cells"]) == len(PRED.CELLS) * 2 * 2
    for c in rep["cells"]:
        for res in list(c["selective"].values()) + [c["control"], c["undetected"]]:
            assert _defect(res) == pytest.approx(res["max_abs_defect"], rel=1e-9, abs=1e-15)
        for e in PRED.ETAS:
            assert c["selective"][str(e)]["eta_check"] == pytest.approx(e, abs=1e-8)
    assert rep["verdict"] == _verdict(rep)


@needs_archive
def test_archive_test_catches_a_corrupted_form():
    rep = json.loads(ARCHIVE.read_text())
    cell = next(c for c in rep["cells"] if (c["a"], c["m"], c["N"], c["K"]) == (0.3, 2, 2000, 150))
    V = np.array(cell["selective"]["0.5"]["V"])
    cell["selective"]["0.5"]["V"] = (0.9 * V).tolist()     # a 10% loss of the second reading
    with pytest.raises(AssertionError):
        assert _defect(cell["selective"]["0.5"]) == pytest.approx(cell["selective"]["0.5"]["max_abs_defect"], rel=1e-9)
