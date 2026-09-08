"""Parity-solvability tests. Freeze ``495f1f1``."""

import math

import numpy as np
import pytest

from geometrodynamics.waves import parity_solvability as ps
from experiments.closure_ledger import parity_solvability_probe as probe


@pytest.fixture(scope="module")
def report():
    """One probe run shared by the tests that inspect its output.

    Each run does 63 fine-grid quadratures and 600 coarse ones, so running it
    per test dominated both this file and the full suite.
    """
    return probe.run_probe()


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


def _all_required():
    return {name: True for name in ps.REQUIRED_CHECKS}


def test_verdict_has_a_stable_schema_and_names_failures():
    good = ps.verdict(_all_required(), "note")
    bad = dict(_all_required()); bad[ps.REQUIRED_CHECKS[0]] = False
    bad = ps.verdict(bad, "note")
    assert set(bad) == set(good)
    assert bad["antipodal_parity_status"] == "UNRESOLVED"
    assert bad["failed_checks"] == [ps.REQUIRED_CHECKS[0]]
    assert good["antipodal_parity_status"] == "SUFFICIENT_NOT_NECESSARY"
    assert good["f6_consequence"].startswith("CONSTRAINT_SOLVABILITY_DOES_NOT")


def test_a_missing_required_check_cannot_yield_an_affirmative_verdict():
    """C3 regression: the first version gated only on ``all(values)``, so an
    unrelated truthy dict returned the full affirmative result and a probe
    that silently dropped every check would still have passed."""
    v = ps.verdict({"unrelated": True}, "x")
    assert v["antipodal_parity_status"] == "UNRESOLVED"
    assert len(v["missing_checks"]) == len(ps.REQUIRED_CHECKS)
    for dropped in ps.REQUIRED_CHECKS:
        partial = {k: True for k in ps.REQUIRED_CHECKS if k != dropped}
        result = ps.verdict(partial, "x")
        assert result["f6_consequence"] == "UNRESOLVED", dropped
        assert result["missing_checks"] == [dropped]


def _synthetic_report(checks):
    """Minimal report with the shape the renderer reads.

    The CLI tests exercise rendering, archival and exit codes, so they must not
    pay for a full probe run.
    """
    return {
        "public_preregistration": ps.PUBLIC_PREREG, "seed": ps.SEED,
        "selection_rule": {"rows": [], "min_adjacent_norm": 3.46,
                           "max_non_adjacent_norm": 1e-14,
                           "selection_rule_holds": True},
        "degree_set_scan": [], "route_agreement": [], "halving": {},
        "maximality": [], "graph_subspaces": [],
        "momentum_sector": {"rows": [], "note": "x"},
        "complete_momentum_audit": {}, "small_projection_counterexample": {},
        "quadrature_exactness": {},
        "checks": checks, "structural_regressions": {},
        "checks_passed": all(checks.values())
                         and not ps.verdict(checks, "n")["missing_checks"],
        "verdict": ps.verdict(checks, "note"),
        "implementation_corrections": ["C1: ..."],
        "scope": "test", "triangle_map": "NOT_DERIVED", "readout": "NOT_DERIVED",
    }


def _drive_cli(tmp_path, monkeypatch, checks):
    monkeypatch.setattr(probe, "run_probe",
                        lambda progress=lambda s: None: _synthetic_report(checks))
    out = tmp_path / "run"
    out.mkdir()
    (out / "probe.md").write_text("STALE PASSING ARCHIVE")
    code = probe.main(["--output-dir", str(out)])
    return code, (out / "probe.md").read_text(), out


def test_dropping_a_check_from_the_probe_fails_the_cli(tmp_path, monkeypatch):
    """C3 regression, end to end: a dropped check must produce UNRESOLVED,
    an overwritten report and a nonzero exit, not a silent pass."""
    checks = {k: True for k in ps.REQUIRED_CHECKS if k != ps.REQUIRED_CHECKS[0]}
    code, text, _ = _drive_cli(tmp_path, monkeypatch, checks)
    assert code == 1
    assert "STALE" not in text and "UNRESOLVED" in text


def test_failing_cli_path_renders_archives_and_exits_nonzero(tmp_path, monkeypatch):
    checks = {k: True for k in ps.REQUIRED_CHECKS}
    checks[ps.REQUIRED_CHECKS[0]] = False
    code, text, out = _drive_cli(tmp_path, monkeypatch, checks)
    assert code == 1
    assert "STALE" not in text and "UNRESOLVED" in text
    assert ps.REQUIRED_CHECKS[0] in text
    import json
    archived = json.loads((out / "probe.json").read_text())
    assert archived["verdict"]["f6_consequence"] == "UNRESOLVED"


def test_passing_cli_path_exits_zero(tmp_path, monkeypatch):
    code, text, _ = _drive_cli(tmp_path, monkeypatch,
                               {k: True for k in ps.REQUIRED_CHECKS})
    assert code == 0
    assert "UNRESOLVED" not in text


def test_probe_records_the_implementation_correction(report):
    assert report["checks_passed"]
    assert any("C1" in c for c in report["implementation_corrections"])
    assert report["triangle_map"] == "NOT_DERIVED"
    assert ps.PUBLIC_PREREG.startswith("495f1f1")


def test_a_small_projection_subspace_evades_within_adjacent_degrees():
    """C4 regression. The first version claimed the parity answer fails "only
    globally". It does not: ``span{x_0, x_1 x_2}`` sits inside the adjacent
    pair ``V_1 + V_2``, mixes parity, and has an identically vanishing dipole
    because ``int x^A x_0 x_1 x_2 dV = 0`` for every ambient coordinate."""
    result = ps.small_projection_counterexample()
    assert result["mixed_parity"] and result["degrees_adjacent"]
    assert result["max_dipole_overlap_route"] < 1e-12
    assert result["max_dipole_stress_route"] < 1e-12
    assert result["generic_adjacent_dipole"] > 1e-3
    assert result["evades"]


def test_all_six_killing_generators_are_used():
    """C5 regression: only three of six were audited."""
    assert len(ps.killing_generators(3)) == 6
    flat = np.array([g.ravel() for g in ps.killing_generators(3)])
    assert np.linalg.matrix_rank(flat, tol=1e-9) == 6


def test_the_incomplete_momentum_audit_could_report_a_false_zero():
    audit = ps.complete_momentum_audit()
    assert audit["originally_reported_three"] < 1e-10
    assert audit["complete_six_on_the_same_data"] > 1e-3
    assert audit["incomplete_audit_would_report_zero"]


def test_gradient_conformal_charges_are_an_independent_sector():
    """Killing charges are diagonal in degree; the gradient charges follow the
    adjacency rule and pair field against momentum."""
    audit = ps.complete_momentum_audit()
    assert audit["split_degree_killing"] < 1e-12
    assert audit["split_degree_hamiltonian_dipole"] < 1e-12
    assert audit["split_degree_gradient_charge"] > 1e-3
    assert audit["gradient_sector_is_independent"]
    assert audit["gradient_same_degree"] < 1e-12
    assert audit["gradient_non_adjacent"] < 1e-12
    assert audit["gradient_obeys_adjacency"]
    assert audit["killing_is_diagonal_in_degree"]


def test_sphere_rules_are_exact_for_these_polynomial_integrands():
    """Licenses the coarse grid used by the high-sample bilinearity scan."""
    result = ps.quadrature_exactness()
    assert result["exact"]
    assert result["fine_vs_coarse_relative"] < 1e-12


def test_subspace_maximality_uses_the_corrected_normalization():
    """C2 regression: the stale factor of two from before C1."""
    import inspect
    assert "2.0 * coeff" not in inspect.getsource(ps.subspace_maximality)
    result = ps.subspace_maximality(1, 2)
    assert abs(result["smallest_singular_value"] - 6.928203) < 1e-5
    assert result["kernel_is_trivial"]


def test_probe_covers_every_pair_through_degree_six_with_momentum(report):
    """C6 regression: the reduction check sampled seven hand-picked sets with
    field data only, where the freeze demands all pairs and nonzero momentum."""
    rows = report["route_agreement"]
    pairs = {tuple(r["degrees"]) for r in rows}
    expected = {(n,) if n == m else (n, m)
                for n in range(1, ps.MAX_DEGREE + 1)
                for m in range(n, ps.MAX_DEGREE + 1)}
    assert pairs == expected
    assert {r["kind"] for r in rows} == {"field", "momentum", "both"}
    assert len(rows) == 3 * len(expected)
    for r in rows:
        assert (r["absolute"] < 1e-12 if r["below_floor"] else r["relative"] < 1e-10)
