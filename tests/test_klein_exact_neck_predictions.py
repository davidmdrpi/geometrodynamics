"""Exact-neck Klein splitting law: frozen predictions, limits, and the in-sample record."""
import importlib.util
import json
from math import gamma, pi, sin, sqrt, tan
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "experiments" / "closure_ledger"
MODULE = LEDGER / "klein_exact_neck_predictions.py"
ARCHIVE_355 = LEDGER / "runs" / "20261118_alice_handle_ports" / "probe.json"
ARCHIVE_362 = LEDGER / "runs" / "20261123_klein_splitting_law" / "probe.json"

_spec = importlib.util.spec_from_file_location("klein_exact_neck_predictions", MODULE)
E = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E)

# The held-out table of docs/klein_exact_neck_prereg.md.
FROZEN = {(0.25, 1): 2.1249e-3, (0.25, 2): 3.0666e-6, (0.35, 1): 8.3890e-3,
          (0.35, 2): 4.7691e-5, (0.35, 3): 2.3234e-7, (0.45, 1): 2.3761e-2,
          (0.45, 2): 3.8202e-4, (0.45, 3): 5.2604e-6, (0.45, 4): 6.7961e-8,
          (0.55, 1): 5.5167e-2, (0.55, 2): 2.0630e-3, (0.55, 3): 6.6050e-5,
          (0.55, 4): 1.9820e-6, (0.55, 5): 5.7432e-8}


def _measured():
    """Every splitting already measured, by cell, at N = 2000 and 4000."""
    out = {}
    a355 = json.loads(ARCHIVE_355.read_text())
    for n in ("2000", "4000"):
        for row in a355["per_size"][n]["measured"]["P4"]:
            pair = row["pairs"][0]
            out.setdefault((row["a"], row["m"]), {})[n] = abs(pair["omega_cos"] - pair["omega_sin"])
    a362 = json.loads(ARCHIVE_362.read_text())
    for g in a362["gated"] + a362["diagnostic"]:
        out[(g["a"], g["m"])] = {n: g["sizes"][n]["delta"] for n in ("2000", "4000")}
    return out


@pytest.mark.parametrize("cell", sorted(FROZEN))
def test_predictions_match_the_frozen_table(cell):
    assert E.prediction(*cell)["predicted_delta"] == pytest.approx(FROZEN[cell], rel=1e-4)


def test_cell_sets():
    assert sorted(E.GATED) == sorted(FROZEN)
    assert not set(E.GATED) & set(E.IN_SAMPLE)
    assert not set(E.DIAGNOSTIC) & set(E.IN_SAMPLE)
    assert set(E.IN_SAMPLE) == set(_measured())


def test_the_module_does_not_call_the_solver():
    assert "alice_handle_ports" not in MODULE.read_text()


def test_predictions_are_converged_in_the_integrator():
    for a, m in E.GATED + E.DIAGNOSTIC:
        coarse = E.prediction(a, m, rtol=1e-11)["predicted_delta"]
        assert coarse == pytest.approx(E.prediction(a, m)["predicted_delta"], rel=1e-9)


@pytest.mark.parametrize("m", [1, 2, 3])
def test_small_aperture_limit_is_the_leading_law(m):
    # C(m) tan^{4m}(a/2) of PR #360, approached at relative O(a^4).
    norm = sqrt(pi) * gamma(m + 1) / (2 * gamma(m + 1.5))
    errors = [abs(E.prediction(a, m)["predicted_split"]
                  / (2 * 4 ** m / ((m + 1) * norm) * tan(a / 2) ** (4 * m)) - 1)
              for a in (0.04, 0.02, 0.01)]
    assert errors[-1] < 1e-7
    assert errors[0] / errors[1] > 10 and errors[1] / errors[2] > 10


def test_ground_state_construction_matches_the_herring_band_state():
    # P_m^m(cos a)^2 / its Legendre norm equals sin^{2m}(a) / N_m.
    for a, m in E.GATED:
        p = E.prediction(a, m)
        flux, _ = E.neck_state(a, m, p["omega"])
        norm = sqrt(pi) * gamma(m + 1) / gamma(m + 1.5)
        assert p["predicted_delta"] == pytest.approx(
            4 * sin(a) ** (2 * m) * abs(flux) / (p["omega"] * norm), rel=1e-12)


def test_in_sample_record():
    """The numbers the pre-registration quotes from the in-sample cells."""
    measured = _measured()
    worst = {"362": 0.0, "0.6": 0.0}
    for (a, m), d in measured.items():
        p = E.prediction(a, m)
        for n in ("2000", "4000"):
            r, r_h = d[n] / p["predicted_delta"], d[n] / p["herring_delta"]
            assert abs(r - 1) <= 0.66 * p["bound"]
            assert abs(r - 1) <= 0.092 * abs(r_h - 1)
        r4 = d["4000"] / p["predicted_delta"] - 1
        assert abs(r4) <= 0.11 * a ** 4
        if a == 0.6:
            worst["0.6"] = max(worst["0.6"], abs(r4))
        elif (a, m) in [(0.2, 1), (0.2, 2), (0.4, 1), (0.4, 2), (0.4, 3), (0.5, 4)]:
            worst["362"] = max(worst["362"], abs(r4))
    # the review of #362 reported 0.234% and 1.02%
    assert worst["362"] == pytest.approx(2.332e-3, rel=1e-3)
    assert worst["0.6"] == pytest.approx(1.017e-2, rel=1e-3)


ARCHIVE = LEDGER / "runs" / "20261124_klein_exact_neck" / "probe.json"
needs_archive = pytest.mark.skipif(not ARCHIVE.exists(), reason="held-out run not yet archived")


def _module(name):
    spec = importlib.util.spec_from_file_location(name, LEDGER / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _recomputed(g, n):
    """r - 1, r_H - 1 and |r - 1|/|r_H - 1| of one gated cell, from its measured splitting."""
    p = E.prediction(g["a"], g["m"])
    delta = g["sizes"][n]["delta"]
    r, r_h = delta / p["predicted_delta"] - 1, delta / p["herring_delta"] - 1
    return r, r_h, abs(r) / abs(r_h)


def _check_stored_fields(rep):
    """Every derived field the probe stored must agree with its recomputation."""
    for g in rep["gated"]:
        for n in ("2000", "4000"):
            r, r_h, fraction = _recomputed(g, n)
            assert r == pytest.approx(g["sizes"][n]["r_minus_1"], rel=1e-9, abs=1e-12)
            assert r_h == pytest.approx(g["sizes"][n]["r_h_minus_1"], rel=1e-9, abs=1e-12)
            assert fraction == pytest.approx(g["sizes"][n]["neck_fraction"], rel=1e-9, abs=1e-12)
        resolution = abs(g["sizes"]["2000"]["delta"] / g["sizes"]["4000"]["delta"] - 1)
        assert resolution == pytest.approx(g["resolution"], rel=1e-9, abs=1e-15)


def _verdict_from_measurements(rep):
    """The frozen gates and verdict rules, recomputed from the measured splittings alone."""
    fails = {}
    for n in ("2000", "4000"):
        g1 = g2 = False
        for g in rep["gated"]:
            r, _, fraction = _recomputed(g, n)
            g1 |= not abs(r) <= 0.2 * g["a"] ** 4 + 1e-3
            g2 |= not fraction <= 0.2
        fails[n] = dict(G1=g1, G2=g2)
    g3 = any(not abs(g["sizes"]["2000"]["delta"] / g["sizes"]["4000"]["delta"] - 1) <= 1e-3
             for g in rep["gated"])
    both = [k for k in ("G1", "G2") if fails["2000"][k] and fails["4000"][k]]
    if both:
        return "COUNTEREXAMPLE"
    if not any(fails["2000"].values()) and not any(fails["4000"].values()) and not g3:
        return "CONFIRMED"
    if g3:
        return "UNRESOLVED_NUMERICALLY"
    return "UNRESOLVED_NUMERICALLY" if not any(fails["4000"].values()) else "NO_VERDICT"


def _corrupt_consistently(rep, cell, factor, sizes):
    """Scale one measured splitting and rewrite its stored fields to match, as the probe would."""
    g = next(x for x in rep["gated"] if (x["a"], x["m"]) == cell)
    for n in sizes:
        g["sizes"][n]["delta"] *= factor
        r, r_h, fraction = _recomputed(g, n)
        g["sizes"][n].update(r_minus_1=r, r_h_minus_1=r_h, neck_fraction=fraction)
    signed = g["sizes"]["2000"]["delta"] / g["sizes"]["4000"]["delta"] - 1
    g.update(resolution=abs(signed), resolution_signed=signed)
    return rep


@needs_archive
def test_archive_verdict_follows_from_the_frozen_gates():
    rep = json.loads(ARCHIVE.read_text())
    assert rep["freeze"] == "834ac70"
    assert sorted((g["a"], g["m"]) for g in rep["gated"]) == sorted(FROZEN)
    measured = _measured()
    assert len(rep["reproduction"]["cells"]) == 2 * len(E.IN_SAMPLE)
    for cell in rep["reproduction"]["cells"]:
        assert cell["archive_delta"] == measured[(cell["a"], cell["m"])][str(cell["N"])]
        assert abs(cell["delta"] / cell["archive_delta"] - 1) < 1e-6
    _check_stored_fields(rep)
    assert rep["verdict"] == _verdict_from_measurements(rep)


@needs_archive
def test_archive_test_rejects_an_inconsistent_splitting():
    # A splitting changed without its derived fields fails the consistency check.
    rep = json.loads(ARCHIVE.read_text())
    rep["gated"][3]["sizes"]["4000"]["delta"] *= 1.05
    with pytest.raises(AssertionError):
        _check_stored_fields(rep)


@needs_archive
@pytest.mark.parametrize("sizes, expected", [(("2000", "4000"), "COUNTEREXAMPLE"),
                                             (("4000",), "UNRESOLVED_NUMERICALLY")])
def test_gates_catch_a_consistent_corruption(sizes, expected):
    # A 5% error at (0.35, 2), with every stored field rewritten to match, passes the
    # consistency check; the recomputed gates must then change the verdict. At both
    # sizes G1 and G2 fail twice; at N = 4000 only, G3 fails first.
    rep = _corrupt_consistently(json.loads(ARCHIVE.read_text()), (0.35, 2), 1.05, sizes)
    _check_stored_fields(rep)
    assert _verdict_from_measurements(rep) == expected
    assert _module("klein_exact_neck_probe").verdict(rep["gated"])[1] == expected


CONTINUUM = LEDGER / "runs" / "20261124_klein_exact_neck" / "continuum_check.json"


@pytest.mark.skipif(not CONTINUUM.exists(), reason="continuum cross-check not archived")
def test_continuum_cross_check():
    """#360's continuum solver agrees with the finite-element splittings, and the frozen
    gates also pass on its values. Recomputed from the archived continuum splittings."""
    archive = {(g["a"], g["m"]): g for g in json.loads(ARCHIVE.read_text())["gated"]}
    check = json.loads(CONTINUUM.read_text())
    gated = [x for x in check["cells"] if x["gated"]]
    assert sorted((x["a"], x["m"]) for x in gated) == sorted(FROZEN)
    for x in gated:
        p = E.prediction(x["a"], x["m"])
        fem = archive[(x["a"], x["m"])]["sizes"]["4000"]["delta"]
        assert -2e-4 < fem / x["continuum_delta"] - 1 < 0
        r = x["continuum_delta"] / p["predicted_delta"] - 1
        r_h = x["continuum_delta"] / p["herring_delta"] - 1
        assert 0 < r <= 0.2 * x["a"] ** 4 + 1e-3
        assert abs(r) <= 0.2 * abs(r_h)
        if x["m"] >= 2:
            kink = abs(p["neck_log_derivative"] / p["band_log_derivative"] - 1)
            assert 0.9 < r / kink < 1.1


@pytest.mark.skipif(not CONTINUUM.exists(), reason="continuum cross-check not archived")
@pytest.mark.parametrize("cell", [(0.25, 2), (0.55, 5)])
def test_continuum_values_reproduce(cell):
    cc = _module("klein_exact_neck_continuum_check")
    row = next(x for x in json.loads(CONTINUUM.read_text())["cells"] if (x["a"], x["m"]) == cell)
    assert cc.continuum_delta(*cell)[0] == pytest.approx(row["continuum_delta"], rel=1e-9)
