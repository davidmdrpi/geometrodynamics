"""Parity-solvability tests. Freeze ``495f1f1``."""

import math

import numpy as np
import pytest

from geometrodynamics.waves import parity_solvability as ps
from experiments.closure_ledger import parity_solvability_probe as probe


def test_kernel_is_the_four_degree_one_modes():
    """``lambda_1 = 3/a^2`` is the only eigenvalue meeting ``Delta + 3/a^2``."""
    assert ps.lambda_n(1) == pytest.approx(3.0)
    for n in range(0, 8):
        if n != 1:
            assert abs(ps.lambda_n(n) - 3.0) > 1e-9
    assert ps.harmonic_multiplet(1).dimension == 4


@pytest.mark.parametrize("n,m", [(n, m) for n in range(1, 6) for m in range(1, 6)])
def test_selection_rule_pairwise(n, m):
    """P3: the degree-1 triple overlap is nonzero only for ``|n-m| = 1``."""
    norm = float(np.linalg.norm(ps.triple_overlap(n, m)))
    if abs(n - m) == 1:
        assert norm > 1e-3
    else:
        assert norm < 1e-12


def test_the_table_reports_both_arms():
    """A table of zeros alone would not establish a selection rule."""
    table = ps.selection_rule_table(5)
    assert table["selection_rule_holds"]
    assert table["min_adjacent_norm"] > 1e-3
    assert table["max_non_adjacent_norm"] < 1e-12
    assert any(r["adjacent"] for r in table["rows"])
    assert any(not r["adjacent"] for r in table["rows"])


@pytest.mark.parametrize("degrees", [[1, 4], [2, 5], [3, 6]])
def test_mixed_parity_without_adjacent_degrees_is_unobstructed(degrees):
    """P4: the decisive counterexamples. Mixed parity, yet solvable."""
    even, odd = ps.parity_split(degrees)
    assert even and odd, "these must genuinely mix parity"
    assert not ps.has_adjacent_pair(degrees)
    rng = np.random.default_rng(7)
    for _ in range(50):
        for amplitude in (1.0, 0.3):
            state = ps.random_state(degrees, rng, amplitude)
            assert np.max(np.abs(ps.dipole_from_overlap(state))) < 1e-12


@pytest.mark.parametrize("degrees", [[1, 2], [2, 3], [3, 4]])
def test_adjacent_mixed_parity_is_obstructed(degrees):
    rng = np.random.default_rng(8)
    worst = max(np.max(np.abs(ps.dipole_from_overlap(ps.random_state(degrees, rng))))
                for _ in range(50))
    assert worst > 1e-3


@pytest.mark.parametrize("degrees", [[1, 3], [2, 4], [1, 3, 5]])
def test_parity_pure_data_is_always_unobstructed(degrees):
    rng = np.random.default_rng(9)
    for _ in range(50):
        assert np.max(np.abs(ps.dipole_from_overlap(
            ps.random_state(degrees, rng)))) < 1e-12


@pytest.mark.parametrize("degrees", [[1, 2], [2, 3], [1, 2, 3]])
def test_reduction_matches_the_independent_improved_stress(degrees):
    """P2, and the check that found correction C1 (a factor of exactly two)."""
    rng = np.random.default_rng(10)
    for _ in range(3):
        agree = ps.route_agreement(ps.random_state(degrees, rng))
        assert not agree["below_floor"]
        assert agree["relative"] < 1e-10


def test_route_agreement_uses_a_scale_floor():
    """A bare ratio is meaningless when both routes return rounding noise."""
    rng = np.random.default_rng(11)
    agree = ps.route_agreement(ps.random_state([1, 4], rng))
    assert agree["below_floor"]
    assert agree["relative"] == 0.0
    assert agree["absolute"] < 1e-12


@pytest.mark.parametrize("degrees", [[1, 2], [2, 3]])
def test_obstruction_is_bilinear_through_the_independent_route(degrees):
    result = ps.bilinearity_in_amplitude(degrees, samples=6)
    assert result["route"] == "improved_stress"
    assert result["max_halving_error"] < 1e-8


@pytest.mark.parametrize("n", [1, 2, 3])
def test_no_nonzero_subspace_evades_an_adjacent_partner(n):
    result = ps.subspace_maximality(n, n + 1)
    assert result["kernel_is_trivial"]
    assert result["smallest_singular_value"] > 1e-9


@pytest.mark.parametrize("n", [1, 2])
def test_no_nonzero_graph_subspace_evades_an_adjacent_pair(n):
    result = ps.graph_subspace_search(n, n + 1)
    assert result["only_trivial_graph"]
    assert result["nullity"] == 0


def test_momentum_obstruction_is_not_a_parity_condition():
    """Parity-pure data kills the dipole yet still carries Killing charge."""
    report = ps.momentum_sector_report()
    assert report["parity_pure_kills_dipole"]
    assert report["parity_pure_can_carry_charge"]
    assert report["note"] == "NOT_PREDICTED_IN_ADVANCE"
    for row in report["rows"]:
        assert row["hamiltonian_dipole"] < 1e-12
        assert row["killing_charge_aligned"] > 1e-3


def test_verdict_has_a_stable_schema_and_names_failures():
    good = ps.verdict({"a": True}, "note")
    bad = ps.verdict({"a": True, "b": False}, "note")
    assert set(bad) == set(good)
    assert bad["antipodal_parity_status"] == "UNRESOLVED"
    assert bad["failed_checks"] == ["b"]
    assert good["antipodal_parity_status"] == "SUFFICIENT_NOT_NECESSARY"
    assert good["f6_consequence"].startswith("CONSTRAINT_SOLVABILITY_DOES_NOT")


def test_failing_cli_path_renders_archives_and_exits_nonzero(tmp_path, monkeypatch):
    real = probe.run_probe

    def failing(progress=lambda s: None):
        report = real(progress=progress)
        report["checks"]["injected failure"] = False
        report["checks_passed"] = False
        report["verdict"] = ps.verdict(report["checks"], "note")
        return report

    monkeypatch.setattr(probe, "run_probe", failing)
    out = tmp_path / "run"
    out.mkdir()
    (out / "probe.md").write_text("STALE PASSING ARCHIVE")
    assert probe.main(["--output-dir", str(out)]) == 1
    written = (out / "probe.md").read_text()
    assert "STALE" not in written and "UNRESOLVED" in written
    assert "injected failure" in written
    import json
    archived = json.loads((out / "probe.json").read_text())
    assert archived["verdict"]["f6_consequence"] == "UNRESOLVED"


def test_probe_records_the_implementation_correction():
    report = probe.run_probe()
    assert report["checks_passed"]
    assert any("C1" in c for c in report["implementation_corrections"])
    assert report["triangle_map"] == "NOT_DERIVED"
    assert ps.PUBLIC_PREREG.startswith("495f1f1")
