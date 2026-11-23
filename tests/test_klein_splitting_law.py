"""Klein splitting law: frozen predictions, small-aperture limit, and the archive."""
import importlib.util
import json
from math import gamma, pi, sqrt, tan
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "experiments" / "closure_ledger"
ARCHIVE = LEDGER / "runs" / "20261123_klein_splitting_law" / "probe.json"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, LEDGER / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PRED = _load("klein_splitting_law_predictions")

# The held-out table of docs/klein_splitting_law_prereg.md (freeze 56055ab).
FROZEN = {(0.2, 1): 8.4421e-4, (0.2, 2): 4.9980e-7, (0.4, 1): 1.3262e-2,
          (0.4, 2): 1.3728e-4, (0.4, 3): 1.1766e-6, (0.5, 4): 3.8280e-7}


@pytest.mark.parametrize("cell", sorted(FROZEN))
def test_predictions_match_the_frozen_table(cell):
    assert PRED.prediction(*cell)["predicted_delta"] == pytest.approx(FROZEN[cell], rel=1e-4)


@pytest.mark.parametrize("m", [1, 2, 3])
def test_small_aperture_limit_is_the_leading_law(m):
    # C(m) tan^{4m}(a/2), the leading law registered in PR #360.
    norm = sqrt(pi) * gamma(m + 1) / (2 * gamma(m + 1.5))
    ratios = []
    for a in (0.04, 0.02, 0.01):
        leading = 2 * 4 ** m / ((m + 1) * norm) * tan(a / 2) ** (4 * m)
        ratios.append(PRED.prediction(a, m)["predicted_split"] / leading)
    assert abs(ratios[-1] - 1) < 1e-4
    assert abs(ratios[0] - 1) > abs(ratios[1] - 1) > abs(ratios[2] - 1)


def test_archive_verdict_follows_from_the_frozen_gates():
    rep = json.loads(ARCHIVE.read_text())
    assert rep["freeze"] == "56055ab"
    assert rep["reproduction"]["passed"] and rep["reproduction"]["worst_relative"] < 1e-6
    for n in ("2000", "4000"):
        values = {(g["a"], g["m"]): g["sizes"][n]["m_r_minus_1"] for g in rep["gated"]}
        assert all(0 < v <= 0.25 for v in values.values())
        c02 = sum(v for (a, _), v in values.items() if a == 0.2) / 2
        c04 = sum(v for (a, _), v in values.items() if a == 0.4) / 3
        assert c02 < c04
        assert all(0.048 <= v <= 0.182 for (a, _), v in values.items() if a == 0.4)
    assert all(g["resolution"] <= 1e-3 for g in rep["gated"])
    assert rep["verdict"] == "CONFIRMED" and rep["reasons"] == []
    assert sorted((g["a"], g["m"]) for g in rep["gated"]) == sorted(FROZEN)
