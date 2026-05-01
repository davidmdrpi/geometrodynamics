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
from experiments.closure_ledger.sk_bridge import (
    D_OPERATOR_VARIANTS,
    DEFAULT_SK_CANDIDATE,
    LEPTON_DEPTHS,
    WIRED_CANDIDATES,
    WKB_CONVENTION,
    Mode,
    _locked_lepton_eigenvectors,
    _operator_phase_matrix,
    lepton_eigenvector_weights,
    phi_convergence_table,
    phi_radial_for_mode,
    phi_radial_from_sk,
    s_k_membership,
    s_k_weighted_modes,
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
        # Default run wires the S(k) bridge → Layer 2.
        assert loaded["experiment_name"] == "closure_ledger.layer2"
        assert loaded["sk_candidate"] == DEFAULT_SK_CANDIDATE
        assert "rows" in loaded and len(loaded["rows"]) == 3
        assert "sk_bridge_blocker" in loaded
        # Markdown is non-empty and contains expected headers.
        md = (out / "summary.md").read_text()
        assert "# Closure-phase ledger — run summary" in md
        assert "## Per-lepton ledger" in md
        assert "## Layer 2 blocker" in md
        assert "## Radial bulk channel — per-mode breakdown" in md


def test_layer1_only_run_reproduces_pre_bridge_universality():
    """Passing sk_candidate='none' reproduces the Layer-1 universality result."""
    result = run_experiment(sk_candidate="none")
    assert result.experiment_name == "closure_ledger.layer1"
    assert result.universality_check["universal"] is True
    assert math.isclose(
        result.universality_check["universal_value"], 0.0, abs_tol=1e-9,
    )


def test_sk_membership_candidate_a():
    """Candidate A: S(k) is the odd-l ground states up to l=k."""
    assert s_k_membership(1, "A_lowest_radial_per_l") == [Mode(l=1, n=0)]
    assert s_k_membership(3, "A_lowest_radial_per_l") == [
        Mode(l=1, n=0), Mode(l=3, n=0),
    ]
    assert s_k_membership(5, "A_lowest_radial_per_l") == [
        Mode(l=1, n=0), Mode(l=3, n=0), Mode(l=5, n=0),
    ]


def test_sk_membership_candidate_b1():
    """Candidate B1: single l=k angular mode at n=0."""
    assert s_k_membership(1, "B1_single_angular_mode") == [Mode(l=1, n=0)]
    assert s_k_membership(3, "B1_single_angular_mode") == [Mode(l=3, n=0)]
    assert s_k_membership(5, "B1_single_angular_mode") == [Mode(l=5, n=0)]


def test_sk_membership_candidate_b2():
    """Candidate B2: single l=1 radial excitation at n = (k-1)/2."""
    assert s_k_membership(1, "B2_single_radial_excitation") == [Mode(l=1, n=0)]
    assert s_k_membership(3, "B2_single_radial_excitation") == [Mode(l=1, n=1)]
    assert s_k_membership(5, "B2_single_radial_excitation") == [Mode(l=1, n=2)]


def test_sk_membership_candidate_c1_is_full_depth_basis():
    """C1: membership covers all B1 modes across the depth basis."""
    expected = [Mode(l=1, n=0), Mode(l=3, n=0), Mode(l=5, n=0)]
    assert s_k_membership(1, "C1_eigenvector_weighted_B1") == expected
    assert s_k_membership(3, "C1_eigenvector_weighted_B1") == expected
    assert s_k_membership(5, "C1_eigenvector_weighted_B1") == expected


def test_sk_membership_candidate_c2_is_full_depth_basis():
    """C2: membership covers all B2 modes across the depth basis."""
    expected = [Mode(l=1, n=0), Mode(l=1, n=1), Mode(l=1, n=2)]
    assert s_k_membership(1, "C2_eigenvector_weighted_B2") == expected
    assert s_k_membership(3, "C2_eigenvector_weighted_B2") == expected
    assert s_k_membership(5, "C2_eigenvector_weighted_B2") == expected


def test_sk_membership_candidate_c1_maslov_matches_c1():
    """C1_maslov_standard shares the C1 mode set; only the policy differs."""
    for k in LEPTON_DEPTHS:
        assert (
            s_k_membership(k, "C1_maslov_standard")
            == s_k_membership(k, "C1_eigenvector_weighted_B1")
        )


def test_sk_membership_candidate_b2_maslov_matches_b2():
    """B2_maslov_standard shares the B2 mode set; only the policy differs."""
    for k in LEPTON_DEPTHS:
        assert (
            s_k_membership(k, "B2_maslov_standard")
            == s_k_membership(k, "B2_single_radial_excitation")
        )


def test_sk_membership_candidate_c2_maslov_matches_c2():
    """C2_maslov_standard shares the C2 mode set; only the policy differs."""
    for k in LEPTON_DEPTHS:
        assert (
            s_k_membership(k, "C2_maslov_standard")
            == s_k_membership(k, "C2_eigenvector_weighted_B2")
        )


def test_locked_lepton_eigenvectors_are_orthonormal():
    """Eigenvectors of the locked block are orthonormal (rows of V.T)."""
    eigenvalues, eigenvectors = _locked_lepton_eigenvectors()
    assert len(eigenvalues) == 3
    assert len(eigenvectors) == 3
    for v in eigenvectors:
        assert len(v) == 3
        norm = sum(c ** 2 for c in v)
        assert math.isclose(norm, 1.0, abs_tol=1e-9)
    for i in range(3):
        for j in range(i + 1, 3):
            dot = sum(
                eigenvectors[i][a] * eigenvectors[j][a] for a in range(3)
            )
            assert math.isclose(dot, 0.0, abs_tol=1e-9)


def test_lepton_eigenvector_weights_sum_to_one_per_species():
    """|v_species,i|² over the depth basis sums to 1 by orthonormality."""
    weights = lepton_eigenvector_weights()
    assert set(weights.keys()) == set(LEPTON_DEPTHS)
    for k, ws in weights.items():
        assert math.isclose(sum(ws), 1.0, abs_tol=1e-9)
        for w in ws:
            assert 0.0 <= w <= 1.0


def test_lepton_eigenvector_weights_are_deterministic():
    """Two independent calls produce bit-identical weight sequences."""
    w1 = lepton_eigenvector_weights()
    w2 = lepton_eigenvector_weights()
    for k in LEPTON_DEPTHS:
        assert w1[k] == w2[k]


def test_lepton_eigenvector_weights_match_expected_dominance():
    """e is dominantly k=1, μ dominantly k=3, τ dominantly k=5."""
    weights = lepton_eigenvector_weights()
    # e: weight on depth-basis index 0 (depth=1) is the largest.
    e_w = weights[1]
    assert e_w[0] > e_w[1] > e_w[2]
    assert e_w[0] > 0.5
    # μ: weight on depth-basis index 1 (depth=3) is the largest.
    mu_w = weights[3]
    assert mu_w[1] > mu_w[0] > mu_w[2]
    assert mu_w[1] > 0.5
    # τ: weight on depth-basis index 2 (depth=5) is the largest.
    tau_w = weights[5]
    assert tau_w[2] > tau_w[0] and tau_w[2] > tau_w[1]
    assert tau_w[2] > 0.99


def test_s_k_weighted_modes_returns_normalized_weights_for_c_candidates():
    """s_k_weighted_modes returns weights summing to 1 for C1 and C2."""
    for cand in ("C1_eigenvector_weighted_B1", "C2_eigenvector_weighted_B2"):
        for k in LEPTON_DEPTHS:
            modes_with_weights = s_k_weighted_modes(k, cand)
            total_w = sum(w for _, w in modes_with_weights)
            assert math.isclose(total_w, 1.0, abs_tol=1e-9)


def test_s_k_weighted_modes_returns_unit_weights_for_unweighted_candidates():
    """A/B1/B2 use unit weights."""
    for cand in (
        "A_lowest_radial_per_l",
        "B1_single_angular_mode",
        "B2_single_radial_excitation",
    ):
        for k in LEPTON_DEPTHS:
            for mode, weight in s_k_weighted_modes(k, cand):
                assert weight == 1.0


def test_sk_membership_none_is_empty():
    """sk_candidate='none' yields an empty membership."""
    assert s_k_membership(1, "none") == []
    assert s_k_membership(5, "none") == []


def test_phi_radial_for_mode_l1_n0_is_finite_and_positive():
    """The Bohr-Sommerfeld radial action for the (l=1, n=0) mode resolves."""
    result = phi_radial_for_mode(l=1, n=0)
    assert result.status == "computed"
    assert result.phi is not None and math.isfinite(result.phi)
    assert result.phi > 0
    # Eigenfrequency for the lowest l=1 mode should be near 1 in geometric units.
    assert result.omega is not None and 0.5 < result.omega < 2.0


def test_phi_radial_from_sk_candidate_a_total_is_sum_of_modes():
    """Φ_radial(k) is the sum of Φ(l, n) over S(k)."""
    result = phi_radial_from_sk(5, "A_lowest_radial_per_l")
    assert result.status == "computed"
    assert result.total_phi is not None
    expected = sum(m.phi for m in result.modes)
    assert math.isclose(result.total_phi, expected, abs_tol=1e-12)
    assert len(result.modes) == 3   # l ∈ {1, 3, 5}


def test_radial_channel_wired_in_default_run():
    """Default run has no row blocked on radial_bulk_phase."""
    result = run_experiment()
    for row in result.rows:
        names = [t["name"] for t in row["terms"]]
        assert "radial_bulk_phase" in names
        radial = next(t for t in row["terms"] if t["name"] == "radial_bulk_phase")
        assert radial["status"] == "available"
        assert radial["value"] is not None
        assert "radial_bulk_phase" not in row["blocking_terms"]


def test_layer2_blocker_marks_implemented_candidate():
    """When candidate A is wired, blocker reports it as the implemented one."""
    result = run_experiment(sk_candidate="A_lowest_radial_per_l")
    blocker = result.sk_bridge_blocker
    assert blocker["implemented_candidate"] == "A_lowest_radial_per_l"
    impls = {
        c["name"]: c["implementation_status"] for c in blocker["candidates"]
    }
    assert impls["A_lowest_radial_per_l"] == "implemented"
    assert impls["B1_single_angular_mode"] == "open"
    assert impls["B2_single_radial_excitation"] == "open"
    assert impls["C1_eigenvector_weighted_B1"] == "open"
    assert impls["C2_eigenvector_weighted_B2"] == "open"
    assert impls["C1_maslov_standard"] == "open"
    assert impls["B2_maslov_standard"] == "open"
    assert impls["C2_maslov_standard"] == "open"
    assert impls["D0_overlap_phase"] == "open"
    assert impls["D1_potential_difference_phase"] == "open"
    assert impls["D2_symmetrized_momentum_phase"] == "open"


def test_phi_radial_uses_wkb_convention_label():
    """ModePhase carries the WKB convention label by default."""
    mp = phi_radial_for_mode(l=1, n=0)
    assert mp.convention == WKB_CONVENTION
    assert mp.maslov_correction == 0.0


def test_phi_grid_convergence_for_l1_n0():
    """Φ(l=1, n=0) is grid-stable across N: variation under ~1e-3."""
    table = phi_convergence_table(1, 0, Ns=(60, 80, 100, 120))
    phis = [row["phi"] for row in table]
    assert all(p is not None for p in phis)
    # All values within 1e-3 of the largest-grid value.
    ref = phis[-1]
    assert max(abs(p - ref) for p in phis) < 1e-3, (
        f"Φ(l=1, n=0) grid spread = {phis} (ref={ref})"
    )


def test_maslov_shift_is_additive():
    """A maslov_correction adds a constant to phi at fixed (l, n, N)."""
    base = phi_radial_for_mode(l=1, n=0)
    shifted = phi_radial_for_mode(l=1, n=0, maslov_correction=math.pi / 4.0)
    assert math.isclose(
        shifted.phi - base.phi, math.pi / 4.0, abs_tol=1e-12,
    )
    assert shifted.maslov_correction == math.pi / 4.0


def test_b1_b2_falsify_universality_under_wkb():
    """Both B1 and B2 break universal closure mod 2π under WKB convention."""
    for cand in ("B1_single_angular_mode", "B2_single_radial_excitation"):
        result = run_experiment(sk_candidate=cand)
        assert result.universality_check["universal"] is False, (
            f"{cand} unexpectedly universal: {result.universality_check}"
        )
        # Spread should be visibly larger than numeric tolerance.
        assert result.universality_check["spread"] > 1e-3


def test_c1_c2_falsify_universality_under_wkb():
    """Both C1 and C2 still break universality under WKB convention."""
    for cand in (
        "C1_eigenvector_weighted_B1",
        "C2_eigenvector_weighted_B2",
    ):
        result = run_experiment(sk_candidate=cand)
        assert result.universality_check["universal"] is False, (
            f"{cand} unexpectedly universal: {result.universality_check}"
        )
        assert result.universality_check["spread"] > 1e-3


def test_s_k_weighted_modes_normalized_for_c1_maslov_standard():
    """C1_maslov_standard inherits C1's normalized weights."""
    for k in LEPTON_DEPTHS:
        modes_with_weights = s_k_weighted_modes(k, "C1_maslov_standard")
        total_w = sum(w for _, w in modes_with_weights)
        assert math.isclose(total_w, 1.0, abs_tol=1e-9)


def test_c1_tightens_spread_relative_to_b1():
    """The eigenvector mixing of C1 reduces the spread compared to B1."""
    b1 = run_experiment(sk_candidate="B1_single_angular_mode")
    c1 = run_experiment(sk_candidate="C1_eigenvector_weighted_B1")
    assert c1.universality_check["spread"] < b1.universality_check["spread"], (
        f"Expected C1 spread ({c1.universality_check['spread']:.3e}) < "
        f"B1 spread ({b1.universality_check['spread']:.3e})"
    )


def test_c1_maslov_standard_emits_turning_point_metadata():
    """C1_maslov_standard reports per-mode turning-point counts and Δ."""
    result = run_experiment(sk_candidate="C1_maslov_standard")
    seen_modes = 0
    for row in result.rows:
        detail = row["radial_detail"]
        assert detail is not None
        for m in detail["modes"]:
            assert m["status"] == "computed"
            assert m["maslov_policy"] == "standard"
            # Per-mode shift = -π/2 * N_turning, exact equality.
            expected = -(math.pi / 2.0) * m["n_turning_points"]
            assert math.isclose(
                m["maslov_correction"], expected, abs_tol=1e-12,
            )
            seen_modes += 1
    assert seen_modes >= 3


def test_c1_maslov_standard_falsifies_universality_under_bohr_sommerfeld():
    """Maslov-shifted C1 still breaks universal closure mod 2π."""
    result = run_experiment(sk_candidate="C1_maslov_standard")
    assert result.universality_check["universal"] is False
    assert result.universality_check["spread"] > 1e-9


def _per_mode_turning_points(rows: list) -> list[int]:
    counts: list[int] = []
    for row in rows:
        detail = row["radial_detail"]
        assert detail is not None
        for m in detail["modes"]:
            counts.append(m["n_turning_points"])
    return counts


def _residues_in_pi(result) -> list[float]:
    return [v / math.pi for v in result.universality_check["per_lepton_mod_2pi"]]


def test_b2_maslov_standard_residues_are_deterministic():
    """B2_maslov_standard mod-2π residues are reproducible to 1e-12."""
    r1 = _residues_in_pi(run_experiment(sk_candidate="B2_maslov_standard"))
    r2 = _residues_in_pi(run_experiment(sk_candidate="B2_maslov_standard"))
    assert len(r1) == 3 and len(r2) == 3
    for a, b in zip(r1, r2):
        assert math.isclose(a, b, abs_tol=1e-12), f"{a} vs {b}"


def test_c2_maslov_standard_residues_are_deterministic():
    """C2_maslov_standard mod-2π residues are reproducible to 1e-12."""
    r1 = _residues_in_pi(run_experiment(sk_candidate="C2_maslov_standard"))
    r2 = _residues_in_pi(run_experiment(sk_candidate="C2_maslov_standard"))
    assert len(r1) == 3 and len(r2) == 3
    for a, b in zip(r1, r2):
        assert math.isclose(a, b, abs_tol=1e-12), f"{a} vs {b}"


def test_b2_maslov_standard_falsifies_universality():
    """B2_maslov_standard breaks universal closure mod 2π."""
    result = run_experiment(sk_candidate="B2_maslov_standard")
    assert result.universality_check["universal"] is False
    assert result.universality_check["spread"] > 1e-9


def test_c2_maslov_standard_falsifies_universality():
    """C2_maslov_standard breaks universal closure mod 2π."""
    result = run_experiment(sk_candidate="C2_maslov_standard")
    assert result.universality_check["universal"] is False
    assert result.universality_check["spread"] > 1e-9


def test_b2_modes_have_non_uniform_turning_point_count_if_solver_says_so():
    """
    The B2 ladder is the canonical place to look for differential N_turning:
    (l=1, n=0) crosses the centrifugal barrier; (l=1, n≥1) sits above it.
    If the solver reports non-uniform counts, the B2 Maslov candidate is in
    the differential regime; if uniform, it is in the degenerate regime.
    Either is valid — we only assert the test surfaces what the solver
    actually reports, so the regime is recorded in the run artifact.
    """
    result = run_experiment(sk_candidate="B2_maslov_standard")
    counts = _per_mode_turning_points(result.rows)
    assert all(c >= 0 for c in counts)
    # Each row contributes exactly one B2 mode (single-mode candidate).
    assert len(counts) == len(result.rows)
    # The ground mode of B2 (k=1 row, l=1 n=0) must have at least one
    # turning point — otherwise either the solver or the integration grid
    # has changed in a way the Maslov diagnostic should flag immediately.
    e_count = next(
        m["n_turning_points"]
        for row in result.rows
        if row["k"] == 1
        for m in row["radial_detail"]["modes"]
    )
    assert e_count >= 1, (
        "B2 ground mode (l=1, n=0) reported zero turning points — "
        "differential Maslov diagnostic is moot if this regresses."
    )


def test_no_maslov_baselines_are_preserved_by_new_candidates():
    """Adding B2/C2 Maslov variants does not perturb the existing baselines."""
    # Snapshot of pre-change residues for the no-Maslov candidates.
    expected = {
        "none": [0.0, 0.0, 0.0],
        "A_lowest_radial_per_l": [0.881876, 1.652620, 0.413097],
        "B1_single_angular_mode": [0.881876, 0.770744, 0.760477],
        "B2_single_radial_excitation": [0.881876, 1.994352, 0.999437],
        "C1_eigenvector_weighted_B1": [0.864195, 0.788396, 0.760506],
        "C2_eigenvector_weighted_B2": [1.059227, 1.817711, 0.998727],
    }
    for cand, ref in expected.items():
        residues = _residues_in_pi(run_experiment(sk_candidate=cand))
        for a, b in zip(residues, ref):
            assert math.isclose(a, b, abs_tol=1e-5), (
                f"{cand} residue drift: got {residues}, expected {ref}"
            )


def test_b2_maslov_standard_residues_relate_to_b2_via_per_mode_shift():
    """
    Each B2 row has a single mode, so the row's mod-2π residue must be
    (B2_residue + (−π/2)·N_turning_for_that_row) mod 2π exactly.
    """
    b2 = run_experiment(sk_candidate="B2_single_radial_excitation")
    b2m = run_experiment(sk_candidate="B2_maslov_standard")
    for r_base, r_shift in zip(b2.rows, b2m.rows):
        modes = r_shift["radial_detail"]["modes"]
        assert len(modes) == 1
        n_tp = modes[0]["n_turning_points"]
        delta = -(math.pi / 2.0) * n_tp
        expected = (r_base["available_total_mod_2pi"] + delta) % TAU
        assert math.isclose(
            r_shift["available_total_mod_2pi"], expected, abs_tol=1e-9,
        ), (
            f"row k={r_base['k']}: expected {expected}, "
            f"got {r_shift['available_total_mod_2pi']} (N_tp={n_tp})"
        )


def test_c2_maslov_standard_residues_relate_to_c2_via_weighted_shift():
    """
    Each C2 row's residue under the Maslov policy equals
    (C2_residue + Σ_i w_i · (−π/2) · N_turning_i) mod 2π.
    """
    c2 = run_experiment(sk_candidate="C2_eigenvector_weighted_B2")
    c2m = run_experiment(sk_candidate="C2_maslov_standard")
    for r_base, r_shift in zip(c2.rows, c2m.rows):
        modes = r_shift["radial_detail"]["modes"]
        delta = sum(
            m["weight"] * (-(math.pi / 2.0) * m["n_turning_points"])
            for m in modes
        )
        expected = (r_base["available_total_mod_2pi"] + delta) % TAU
        assert math.isclose(
            r_shift["available_total_mod_2pi"], expected, abs_tol=1e-9,
        ), (
            f"row k={r_base['k']}: expected {expected}, "
            f"got {r_shift['available_total_mod_2pi']}"
        )


def test_c1_maslov_standard_preserves_c1_spread_when_uniform():
    """
    If every C1 mode has the same turning-point count, the Maslov shift
    is a uniform offset and the spread mod 2π must equal the C1 spread
    mod 2π exactly (modular arithmetic preserves differences for uniform
    shifts smaller than a single 2π wrap).
    """
    c1 = run_experiment(sk_candidate="C1_eigenvector_weighted_B1")
    c1m = run_experiment(sk_candidate="C1_maslov_standard")
    # Collect per-mode N_tp across all rows; uniform iff all equal.
    counts: list[int] = []
    for row in c1m.rows:
        detail = row["radial_detail"]
        assert detail is not None
        for m in detail["modes"]:
            counts.append(m["n_turning_points"])
    if len(set(counts)) == 1:
        assert math.isclose(
            c1m.universality_check["spread"],
            c1.universality_check["spread"],
            abs_tol=1e-9,
        )


def test_run_comparison_covers_layer1_baseline_and_wired_candidates():
    """run_comparison runs Layer-1 + every wired candidate and labels each."""
    from experiments.closure_ledger.runner import run_comparison

    comparison = run_comparison()
    candidates_run = comparison["candidates_run"]
    assert "none" in candidates_run
    for cand in WIRED_CANDIDATES:
        assert cand in candidates_run

    by_cand = {row["candidate"]: row for row in comparison["status_table"]}
    assert by_cand["none"]["result"] == "PASS"
    for cand in WIRED_CANDIDATES:
        assert by_cand[cand]["result"] == "FAIL"


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


# --- Candidate D: operator-valued radial phase ----------------------------

def test_d_family_membership_matches_c1():
    """All D variants share the C1 depth-basis labels (l=1,3,5; n=0)."""
    for variant in D_OPERATOR_VARIANTS:
        for k in LEPTON_DEPTHS:
            assert (
                s_k_membership(k, variant)
                == s_k_membership(k, "C1_eigenvector_weighted_B1")
            )


def test_d_family_phase_matrix_is_hermitian():
    """Each D variant's Φ matrix must equal its transpose (real Hermitian)."""
    for variant in D_OPERATOR_VARIANTS:
        Phi, _omegas, _labels = _operator_phase_matrix(variant)
        for i in range(3):
            for j in range(3):
                assert math.isclose(
                    Phi[i][j], Phi[j][i], abs_tol=1e-12,
                ), f"{variant}: Phi[{i}][{j}]={Phi[i][j]} ≠ Phi[{j}][{i}]={Phi[j][i]}"


def test_d_family_phase_matrix_is_deterministic():
    """Repeated calls produce bit-identical matrices (lru_cache + pure compute)."""
    for variant in D_OPERATOR_VARIANTS:
        m1, _, _ = _operator_phase_matrix(variant)
        m2, _, _ = _operator_phase_matrix(variant)
        assert m1 == m2


def test_d0_diagonal_is_phase_scale_pi():
    """D0 diagonal Φ_ii = π · ⟨u_i|u_i⟩ = π by L²-normalization."""
    Phi, _, _ = _operator_phase_matrix("D0_overlap_phase")
    for i in range(3):
        assert math.isclose(Phi[i][i], math.pi, rel_tol=1e-9, abs_tol=1e-9)


def test_d1_diagonal_is_zero():
    """D1 diagonal Φ_ii = ⟨u_i|V_i − V_i|u_i⟩ = 0 by construction."""
    Phi, _, _ = _operator_phase_matrix("D1_potential_difference_phase")
    for i in range(3):
        assert math.isclose(Phi[i][i], 0.0, abs_tol=1e-12)


def test_d_family_off_diagonals_are_nonzero():
    """
    Each D variant must produce non-trivial off-diagonal phase coupling;
    otherwise the operator candidate collapses to scalar B1 modes and
    contributes no new physics relative to C1.
    """
    for variant in D_OPERATOR_VARIANTS:
        Phi, _, _ = _operator_phase_matrix(variant)
        off_diag_max = max(
            abs(Phi[i][j])
            for i in range(3) for j in range(3) if i != j
        )
        assert off_diag_max > 1e-3, (
            f"{variant}: off-diagonal max |Φ_ij|={off_diag_max} "
            "is below 1e-3; operator candidate is degenerate."
        )


def test_d_family_modes_carry_pair_metadata():
    """Each D-variant ModePhase entry exposes (l_i, n_i) and (l_j, n_j)."""
    for variant in D_OPERATOR_VARIANTS:
        for k in LEPTON_DEPTHS:
            r = phi_radial_from_sk(k, variant)
            assert r.status == "computed"
            assert len(r.modes) == 9   # 3×3 matrix
            for mp in r.modes:
                assert mp.pair_l is not None
                assert mp.pair_n is not None
                assert mp.operator_variant == variant


def test_d_family_total_phi_equals_quadratic_form():
    """RadialPhaseResult.total_phi equals v_species^T Φ v_species exactly."""
    _, eigenvectors = _locked_lepton_eigenvectors()
    for variant in D_OPERATOR_VARIANTS:
        Phi, _, _ = _operator_phase_matrix(variant)
        for sp_idx, k in enumerate(LEPTON_DEPTHS):
            v = eigenvectors[sp_idx]
            quad = sum(
                v[i] * Phi[i][j] * v[j]
                for i in range(3) for j in range(3)
            )
            r = phi_radial_from_sk(k, variant)
            assert math.isclose(r.total_phi, quad, abs_tol=1e-12), (
                f"{variant} k={k}: total_phi={r.total_phi}, v^T Φ v={quad}"
            )


def test_d_family_residues_are_deterministic():
    """D-family per-lepton residues are reproducible to 1e-12."""
    for variant in D_OPERATOR_VARIANTS:
        r1 = _residues_in_pi(run_experiment(sk_candidate=variant))
        r2 = _residues_in_pi(run_experiment(sk_candidate=variant))
        for a, b in zip(r1, r2):
            assert math.isclose(a, b, abs_tol=1e-12), (
                f"{variant}: {a} vs {b}"
            )


def test_d_family_runs_to_completion_through_runner():
    """Each D variant runs end-to-end and produces a Layer-2 result."""
    for variant in D_OPERATOR_VARIANTS:
        result = run_experiment(sk_candidate=variant)
        assert result.experiment_name == "closure_ledger.layer2"
        assert result.universality_check["per_lepton_mod_2pi"]
        # Spread fields populated.
        assert result.universality_check["spread"] >= 0
        assert result.universality_check["circular_spread"] >= 0


def test_circular_spread_caps_linear_spread():
    """Circular spread is always ≤ linear spread for any residue distribution."""
    for cand in WIRED_CANDIDATES + ("none",):
        result = run_experiment(sk_candidate=cand)
        uc = result.universality_check
        if not uc["per_lepton_mod_2pi"]:
            continue
        assert uc["circular_spread"] <= uc["spread"] + 1e-12, (
            f"{cand}: circular={uc['circular_spread']}, linear={uc['spread']}"
        )


def test_no_d_baselines_change_existing_candidate_residues():
    """Adding D-family does not perturb the previously-recorded residues."""
    expected = {
        "none": [0.0, 0.0, 0.0],
        "A_lowest_radial_per_l": [0.881876, 1.652620, 0.413097],
        "B1_single_angular_mode": [0.881876, 0.770744, 0.760477],
        "B2_single_radial_excitation": [0.881876, 1.994352, 0.999437],
        "C1_eigenvector_weighted_B1": [0.864195, 0.788396, 0.760506],
        "C2_eigenvector_weighted_B2": [1.059227, 1.817711, 0.998727],
        "C1_maslov_standard": [0.364195, 0.288396, 0.260506],
        "B2_maslov_standard": [0.381876, 1.994352, 0.999437],
        "C2_maslov_standard": [0.638760, 1.738289, 0.998617],
    }
    for cand, ref in expected.items():
        residues = _residues_in_pi(run_experiment(sk_candidate=cand))
        for a, b in zip(residues, ref):
            assert math.isclose(a, b, abs_tol=1e-5), (
                f"{cand} residue drift: got {residues}, expected {ref}"
            )
