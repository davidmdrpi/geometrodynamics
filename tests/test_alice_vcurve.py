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


ARCHIVE = RUN / "probe.json"
needs_archive = pytest.mark.skipif(not ARCHIVE.exists(), reason="held-out run not archived")


def _recomputed_verdict(rep, frozen):
    """Every gate recomputed from the stored defects and the frozen predictions."""
    pred = {(r["a"], r["eta"], round(r["tau"], 6)): r["total"] for r in frozen["vcurve"]}
    vc = [r for r in rep["vcurve"] if not r.get("unreachable")]
    g1 = all(1 / 3 <= r["max_abs_defect"] / pred[(r["a"], r["eta"], round(r["tau"], 6))] <= 3 for r in vc)
    grid = [round(t, 6) for t in PRED.GRID]
    g2, minima = True, {}
    for a in (PRED.A1, PRED.A2):
        for eta in PRED.ETAS:
            rows = [r for r in vc if r["a"] == a and r["eta"] == eta]
            best = min(rows, key=lambda r: r["max_abs_defect"])
            minima[(a, eta)] = best["max_abs_defect"]
            i = grid.index(round(best["tau"], 6))
            j = grid.index(round(min((r for r in frozen["vcurve"] if r["a"] == a and r["eta"] == eta),
                                     key=lambda r: r["total"])["tau"], 6))
            g2 &= 0 < i < len(grid) - 1 and abs(i - j) <= 1
    g3 = all(minima[(PRED.A2, e)] <= minima[(PRED.A1, e)] / 3 and minima[(PRED.A2, e)] <= 1e-4 for e in PRED.ETAS)
    taus = np.array([r["tau"] for r in rep["null_scan"]])
    vals = np.array([r["max_abs_defect"] for r in rep["null_scan"]])
    g4 = True
    for tk in frozen["null_taus"]:
        win = np.abs(taus - tk) <= 0.64
        i = np.argmin(np.where(win, vals, np.inf))
        g4 &= abs(taus[i] - tk) <= 0.15 and vals[i] <= 0.4 * vals[win].max()
    env = {(r["tau"], r["profile"]): r["max_abs_defect"] for r in rep["envelope"]}
    penv = {(r["tau"], r["profile"]): r["total"] for r in frozen["envelope"]}
    g5 = all(1 / 2.5 <= (env[(t, k)] / env[(t, "sin2")]) / (penv[(t, k)] / penv[(t, "sin2")]) <= 2.5
             for t in PRED.ENVELOPE_TAUS for k in ("sin4", "bump"))
    g5 &= env[(16.0, "sin4")] <= 0.2 * env[(16.0, "sin2")] and env[(16.0, "bump")] > env[(16.0, "sin2")]
    failed = [name for name, ok in zip(("G1", "G2", "G3", "G4", "G5"), (g1, g2, g3, g4, g5)) if not ok]
    unreachable = [r for r in rep["vcurve"] if r.get("unreachable")]
    if failed:
        return "COUNTEREXAMPLE", failed
    return ("NO_VERDICT" if unreachable else "CONFIRMED"), []


@needs_archive
def test_archive_verdict_follows_from_stored_values():
    rep = json.loads(ARCHIVE.read_text())
    frozen = json.loads((RUN / "frozen_predictions.json").read_text())
    assert rep["freeze"] == "291ab5f"
    assert len(rep["vcurve"]) == 2 * 2 * len(PRED.GRID)
    assert len(rep["null_scan"]) == len(PRED.NULL_SCAN)
    verdict, failed = _recomputed_verdict(rep, frozen)
    assert verdict == rep["verdict"]
    if verdict == "COUNTEREXAMPLE":
        assert failed == rep["reasons"]


@needs_archive
def test_archive_test_catches_a_corrupted_value():
    rep = json.loads(ARCHIVE.read_text())
    frozen = json.loads((RUN / "frozen_predictions.json").read_text())
    for row in rep["vcurve"]:                         # lift the smaller aperture's curve tenfold
        if row["a"] == PRED.A2 and row["eta"] == 0.9 and not row.get("unreachable"):
            row["max_abs_defect"] *= 10
    _, failed = _recomputed_verdict(rep, frozen)
    assert "G3" in failed
