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


# --- Dynamic-phase probe -----------------------------------------------

def test_dynamic_phase_probe_runs_to_completion():
    """The dynamic-phase probe builds a structured summary without errors."""
    from experiments.closure_ledger.dynamic_phase_probe import (
        BASELINE_RESIDUES_PI,
        run_probe,
    )

    summary = run_probe()
    assert "best_per_candidate" in summary
    for cand in BASELINE_RESIDUES_PI:
        assert cand in summary["best_per_candidate"]
        info = summary["best_per_candidate"][cand]
        assert "baseline_spread_rad" in info
        assert "best_per_mechanism" in info
        assert "overall_best" in info
        # Overall best spread must be ≤ baseline (the m=0 / α=0 grid point
        # always exists and reproduces the baseline).
        assert (
            info["overall_best"]["circular_spread_rad"]
            <= info["baseline_spread_rad"] + 1e-12
        )


def test_dynamic_phase_probe_natural_mechanisms_do_not_close_c1_or_d1():
    """
    Empirical claim: none of the four natural BAM loop phases close
    C1 or D1 mod 2π (within 1e-9). If a future repo change ever flips
    this — e.g. a new mechanism is added that DOES close — this test
    will fail and force a re-read of the verdict text in the probe.
    """
    from experiments.closure_ledger.dynamic_phase_probe import run_probe

    summary = run_probe()
    for cand, info in summary["best_per_candidate"].items():
        assert info["any_mechanism_closes"] is False, (
            f"{cand}: a natural loop phase unexpectedly closes the "
            f"residual. Best: {info['overall_best']}"
        )


def test_dynamic_phase_probe_baseline_grid_point_exists():
    """
    The probe's null hypothesis (zero dynamic phase) must be in every
    mechanism's parameter grid, so the baseline is reachable. If this
    invariant is dropped, the 'helps' verdict can become inconsistent.
    """
    from experiments.closure_ledger.dynamic_phase_probe import (
        _build_mechanisms,
    )
    for m in _build_mechanisms():
        deltas_at_each_grid_point = [
            [m.delta_k(k, p) for k in (1, 3, 5)]
            for p in m.parameter_grid
        ]
        any_zero = any(
            all(abs(d) < 1e-12 for d in deltas)
            for deltas in deltas_at_each_grid_point
        )
        assert any_zero, (
            f"{m.name} has no zero-Δ grid point; baseline is unreachable."
        )


# --- Geometric-Hamiltonian probe ---------------------------------------

def test_geometric_hamiltonian_probe_runs_to_completion():
    """The geometric-Hamiltonian probe builds a structured summary."""
    from experiments.closure_ledger.geometric_hamiltonian_probe import (
        run_probe,
    )

    summary = run_probe()
    assert "variants" in summary
    assert len(summary["variants"]) >= 5
    for v in summary["variants"]:
        assert "eigenvalues" in v and len(v["eigenvalues"]) == 3
        assert "eigenvectors" in v and len(v["eigenvectors"]) == 3
        assert "closure_spread_b1_rad" in v
        assert "closure_spread_d1_rad" in v


def test_geometric_hamiltonian_d1_identity_holds_to_machine_precision():
    """
    Verifies the structural identity ⟨u_i|V_j−V_i|u_j⟩ = (ω_j²−ω_i²)·⟨u_i|u_j⟩
    on the canonical N=80 Chebyshev grid. This is the operator-theoretic
    statement that GH_C and D1 derive their off-diagonals from.
    """
    import numpy as np
    from experiments.closure_ledger.geometric_hamiltonian_probe import (
        _build_radial_basis,
    )

    b = _build_radial_basis()
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            ovl = float(np.trapezoid(b.us[i] * b.us[j], b.rstar))
            lhs = float(np.trapezoid(
                b.us[i] * (b.Vs[j] - b.Vs[i]) * b.us[j], b.rstar,
            ))
            rhs = (b.omegas[j] ** 2 - b.omegas[i] ** 2) * ovl
            assert math.isclose(lhs, rhs, abs_tol=1e-9), (
                f"({i},{j}): lhs={lhs}, rhs={rhs}"
            )


def test_geometric_hamiltonian_no_variant_jointly_passes():
    """
    Empirical claim, recorded as a regression: no Tangherlini-matrix-element
    Hamiltonian both (a) reproduces observed lepton mass ratios within a
    factor of 2 AND (b) closes the closure ledger non-trivially. If a
    future variant ever flips this — flag it as the headline finding.
    """
    from experiments.closure_ledger.geometric_hamiltonian_probe import (
        run_probe,
    )
    summary = run_probe()
    for v in summary["variants"]:
        joint = (
            v["matches_observed_within_factor_of_2"]
            and (
                v["closes_b1_within_1e_9"]
                or (v["closes_d1_within_1e_9"]
                    and not v["is_trivial_d1_closure"])
            )
        )
        assert not joint, (
            f"{v['name']} unexpectedly satisfies both mass ratio AND "
            f"non-trivial closure: {v}"
        )


def test_geometric_hamiltonian_flags_trivial_d1_closure():
    """GH_A's identity eigenvectors × D1's zero diagonal → trivial closure."""
    from experiments.closure_ledger.geometric_hamiltonian_probe import (
        run_probe,
    )
    summary = run_probe()
    by_name = {v["name"]: v for v in summary["variants"]}
    assert by_name["GH_A_diagonal_omega_squared"]["is_trivial_d1_closure"] is True
    # Other variants have non-identity eigenvectors, so even if their D1
    # spread happened to vanish it would not be flagged trivial.
    for name in (
        "GH_B_potential_average",
        "GH_C_induced_hamiltonian",
        "GH_D_omega_squared_overlap",
        "GH_E_symmetric_momentum",
    ):
        assert by_name[name]["is_trivial_d1_closure"] is False


# --- Composed-Hamiltonian probe ----------------------------------------

def test_composed_hamiltonian_probe_runs_to_completion():
    """The composed-Hamiltonian probe builds a structured summary."""
    from experiments.closure_ledger.composed_hamiltonian_probe import (
        run_probe,
    )

    summary = run_probe()
    assert "variants" in summary
    assert len(summary["variants"]) >= 5
    for v in summary["variants"]:
        assert len(v["eigenvalues"]) == 3
        assert len(v["eigenvectors"]) == 3
        # Each variant carries the near-identity diagnostic.
        assert "eigenvector_max_off_diagonal" in v
        assert "is_near_identity_d1_regime" in v


def test_composed_hamiltonian_uses_only_geometric_constants():
    """
    Probe's composition formulas must reference only β = 50π and
    action_base = 2π as numerical constants — no fitted values.
    """
    from experiments.closure_ledger.composed_hamiltonian_probe import (
        BETA, ACTION_BASE,
    )
    assert math.isclose(BETA, 50.0 * math.pi)
    assert math.isclose(ACTION_BASE, 2.0 * math.pi)
    # 4·β / (2π) is the integer closure-quantum count (= 100).
    assert round(4 * BETA / (2 * math.pi)) == 100


def test_composed_hamiltonian_tau_eigenvalue_dominated_by_closure_quantum():
    """
    Across every variant, λ_τ should be within ~15% of 4β = 200π. This
    is the structural claim that the τ row mass comes overwhelmingly
    from the closure quantum, not from radial matrix elements.
    """
    from experiments.closure_ledger.composed_hamiltonian_probe import (
        run_probe, BETA,
    )
    target = 4 * BETA   # 200π
    summary = run_probe()
    for v in summary["variants"]:
        lam_tau = v["eigenvalues"][2]
        rel = abs(lam_tau - target) / target
        assert rel < 0.15, (
            f"{v['name']}: λ_τ = {lam_tau:.3f}, target = {target:.3f} "
            f"(4β = 200π). Relative diff {rel:.3f}."
        )


def test_composed_hamiltonian_no_variant_matches_observed_within_factor_5():
    """
    Empirical claim: even with the closure quantum included, no variant
    in the catalog reaches factor-5 agreement on the muon row. If a
    future variant ever does, it should flip this test and force a
    re-read of the verdict text.
    """
    from experiments.closure_ledger.composed_hamiltonian_probe import (
        run_probe,
    )
    summary = run_probe()
    for v in summary["variants"]:
        assert not v["matches_observed_within_factor_of_5"], (
            f"{v['name']} unexpectedly within factor 5 of observed: "
            f"{v['mass_ratios_predicted_to_e']}"
        )


def test_composed_hamiltonian_flags_near_identity_d1_regime():
    """HC_8 (large diagonal prefactor) should be flagged near-identity."""
    from experiments.closure_ledger.composed_hamiltonian_probe import (
        run_probe,
    )
    summary = run_probe()
    by_name = {v["name"]: v for v in summary["variants"]}
    # HC_8 has (2π)² ≈ 39.5 prefactor on ω², dominating off-diag D1
    # entries (which are O(0.5)). Eigenvectors near-identity → flag.
    assert by_name[
        "HC_8_action_squared_plus_closure_with_d1"
    ]["is_near_identity_d1_regime"] is True


# --- Pinhole-origin probe ----------------------------------------------

def test_pinhole_origin_probe_runs_to_completion():
    """The pinhole-origin probe builds a structured summary across categories."""
    from experiments.closure_ledger.pinhole_origin_probe import run_probe
    summary = run_probe()
    for cat in ("barrier_sum", "non_ground_mode", "matrix_element", "closure_quantum"):
        assert cat in summary["candidates_per_category"], (
            f"missing category {cat}"
        )
        assert summary["candidates_per_category"][cat], (
            f"category {cat} is empty"
        )
        assert cat in summary["best_per_category"]


def test_pinhole_origin_barrier_sum_within_3pct_of_gamma_lepton():
    """
    The QCD-style Σ_{l=1..5} V_max(l) on the Chebyshev grid must remain
    within 3% of γ_lepton = 22.5. This is the headline structural
    finding of the probe; if the radial solver or grid changes, the
    test will surface the drift.
    """
    from experiments.closure_ledger.pinhole_origin_probe import run_probe
    summary = run_probe()
    barrier = next(
        c for c in summary["candidates_per_category"]["barrier_sum"]
        if c["name"] == "Sum_l=1..5_V_max[chebyshev_N80]"
    )
    assert abs(barrier["pct_diff_lepton"]) < 3.0, (
        f"Σ V_max diverged: {barrier['value']} ({barrier['pct_diff_lepton']:+.3f}%)"
    )


def test_pinhole_origin_locked_baseline_reproduces_observed_masses():
    """
    The mass-sensitivity reference row must agree with the locked
    surrogate: γ = 22.5 should give m_μ and m_τ within 0.2% of PDG.
    """
    from experiments.closure_ledger.pinhole_origin_probe import run_probe
    summary = run_probe()
    locked = next(
        s for s in summary["mass_sensitivity"]
        if s["label"] == "locked_baseline_22.5"
    )
    assert locked["rel_err_pct"][3] < 0.2
    assert locked["rel_err_pct"][5] < 0.2


def test_pinhole_origin_geometric_gamma_breaks_mass_match_significantly():
    """
    Plugging the bare geometric γ = Σ V_max ≈ 22.008 into the locked
    block must significantly degrade the muon-mass match (>10% error).
    This is the empirical statement that the locked γ is fine-tuned at
    the percent level on top of its geometric origin.
    """
    from experiments.closure_ledger.pinhole_origin_probe import (
        _mass_sensitivity_for_gamma,
    )
    sens = _mass_sensitivity_for_gamma(22.0082)
    assert sens["rel_err_pct"][3] > 10.0, (
        f"Geometric γ unexpectedly preserved muon match: {sens}"
    )


# --- γ-offset probe ----------------------------------------------------

def test_gamma_offset_probe_runs_to_completion():
    """The γ-offset probe builds a structured candidate list."""
    from experiments.closure_ledger.gamma_offset_probe import run_probe
    summary = run_probe()
    assert "candidates" in summary
    assert len(summary["candidates"]) >= 5
    assert "best_by_gamma_closeness" in summary
    assert "best_by_muon_match" in summary


def test_gamma_offset_l_eq_0_extension_closes_offset_within_1pct():
    """
    Σ_{l=0..5} V_max(l) — the QCD pinhole with the 5D-specific l=0
    barrier added — must land within 1% of γ_lepton = 22.5. This is
    the headline finding of the probe.
    """
    from experiments.closure_ledger.gamma_offset_probe import run_probe
    summary = run_probe()
    cands = {c["name"]: c for c in summary["candidates"]}
    c = cands["extend_Sum_l_0to5_V_max"]
    assert abs(c["pct_diff_lepton"]) < 1.0, (
        f"l=0 extension drifted: γ = {c['value']} ({c['pct_diff_lepton']:+.3f}%)"
    )


def test_gamma_offset_l_eq_0_extension_recovers_muon_to_within_5pct():
    """
    Re-evaluating the locked block with γ = Σ_{l=0..5} V_max(l) must
    bring the muon-mass error down to < 5% (vs ~64% at the bare
    Σ_{l=1..5} value).
    """
    from experiments.closure_ledger.gamma_offset_probe import run_probe
    summary = run_probe()
    cands = {c["name"]: c for c in summary["candidates"]}
    c = cands["extend_Sum_l_0to5_V_max"]
    assert c["muon_mass_err_pct"] < 5.0, (
        f"l=0-extended γ does not recover the muon: "
        f"err = {c['muon_mass_err_pct']:.3f}%"
    )
    assert c["tau_mass_err_pct"] < 5.0


def test_gamma_offset_has_a_joint_winner():
    """
    The probe must return at least one (γ-within-1%, muon-within-5%)
    joint winner. Otherwise the offset is not closed by any candidate.
    """
    from experiments.closure_ledger.gamma_offset_probe import run_probe
    summary = run_probe()
    assert summary["joint_winners_within_1pct_gamma_and_5pct_muon"], (
        "no candidate jointly closes both γ and muon-mass criteria"
    )


# --- Quark β-origin probe ----------------------------------------------

def test_quark_beta_origin_probe_runs_to_completion():
    """The quark β-origin probe builds a structured candidate list."""
    from experiments.closure_ledger.quark_beta_origin_probe import run_probe
    summary = run_probe()
    assert summary["targets"]["N_quark"] == 466
    assert summary["targets"]["delta_N_quark_lepton_gap"] == 366
    assert summary["n_candidates_total"] > 100
    assert summary["principled_categories"]
    # Every principled category should have at least one near-miss.
    assert summary["best_per_category_for_466"]
    assert summary["best_per_category_for_366"]


def test_quark_beta_origin_no_principled_exact_match():
    """
    Empirical claim: no principled enumeration in this catalog produces
    N = 466 or ΔN = 366 exactly. A future positive result should flip
    this test and force a re-read of the verdict.
    """
    from experiments.closure_ledger.quark_beta_origin_probe import run_probe
    summary = run_probe()
    assert not summary["principled_exact_matches_466"], (
        f"unexpected principled exact match on N=466: "
        f"{summary['principled_exact_matches_466']}"
    )
    assert not summary["principled_exact_matches_366"], (
        f"unexpected principled exact match on ΔN=366: "
        f"{summary['principled_exact_matches_366']}"
    )


def test_quark_beta_origin_k5_squared_plus_one_pattern_brackets_targets():
    """
    `(k_5² + 1) = 26` near-misses bracket both targets within ±2:
    18·26 = 468 (vs N=466) and 14·26 = 364 (vs ΔN=366). This pattern
    is the principled candidate the probe surfaces; if it disappears
    after a structure change, the verdict must be re-read.
    """
    from experiments.closure_ledger.quark_beta_origin_probe import run_probe
    summary = run_probe()
    p_near_466 = summary["principled_near_matches_466_within_5pct"]
    p_near_366 = summary["principled_near_matches_366_within_5pct"]
    has_18x26 = any(c["value"] == 468 for c in p_near_466)
    has_14x26 = any(c["value"] == 364 for c in p_near_366)
    assert has_18x26, "18·(k_5²+1) = 468 missing from 466 near-miss list"
    assert has_14x26, "14·(k_5²+1) = 364 missing from 366 near-miss list"


# --- Quark β: focused boundary-correction probe -----------------------

def test_quark_beta_boundary_probe_runs_to_completion():
    """The boundary-correction probe builds a structured fit catalog."""
    from experiments.closure_ledger.quark_beta_boundary_probe import run_probe
    summary = run_probe()
    assert summary["targets"]["N_lepton"] == 100
    assert summary["targets"]["N_quark"] == 466
    assert summary["targets"]["delta_N"] == 366
    assert summary["all_fits"]
    assert summary["best_fit"]
    assert summary["robustness_snapshots_under_best_fit"]


def test_quark_beta_boundary_best_fit_is_k_5_with_plus_one_correction():
    """
    Headline finding: the cleanest joint decomposition is
        N_l =  20 · k_5     (δ = 0)
        N_q =  93 · k_5 + 1 (δ = +1)
        ΔN  =  73 · k_5 + 1 (δ = +1)
    The probe must continue to surface this as the best fit.
    """
    from experiments.closure_ledger.quark_beta_boundary_probe import run_probe
    summary = run_probe()
    best = summary["best_fit"]
    assert best["unit_name"] == "k_5"
    assert best["unit_value"] == 5
    assert best["m_lepton"] == 20 and best["delta_lepton"] == 0
    assert best["m_quark"] == 93 and best["delta_quark"] == 1
    assert best["m_delta"] == 73 and best["delta_gap"] == 1
    assert best["sum_abs_delta"] == 2
    assert best["deltas_in_natural_set"] is True
    assert best["deltas_consistent_across_quark_and_gap"] is True
    assert best["deltas_zero_for_lepton"] is True


def test_quark_beta_boundary_k5_squared_plus_one_is_strictly_worse():
    """
    The earlier near-miss `(k_5² + 1) = 26` must rank strictly below
    `k_5` on the joint-cleanness measure — it had Σ|δ| = 8 with
    inconsistent δ across targets, while `k_5` has Σ|δ| = 2 with
    consistent +1 boundary correction.
    """
    from experiments.closure_ledger.quark_beta_boundary_probe import run_probe
    summary = run_probe()
    fits = {f["unit_name"]: f for f in summary["all_fits"]}
    f_k5 = fits["k_5"]
    f_k5sq1 = fits["k_5² + 1"]
    assert f_k5["sum_abs_delta"] < f_k5sq1["sum_abs_delta"]
    assert f_k5["deltas_consistent_across_quark_and_gap"]
    # k_5²+1 has δ = (-2, +2) which are not consistent across q and ΔN.
    assert not f_k5sq1["deltas_consistent_across_quark_and_gap"]


def test_quark_beta_boundary_delta_invariance_under_documented_drifts():
    """
    Audit-corrected statement (was overstated in the original probe).
    Using the actual §8 ablation N values (not extrapolated), only ~33%
    of perturbations leave δ at +1; the structural reading is
    descriptively useful for the BASELINE locked value but does NOT
    survive per-species or anchor perturbations. The audit probe
    `quark_beta_robustness_audit.py` carries the rigorous treatment.
    Pin the actual rate so any future drift in the §8 numbers surfaces.
    """
    from experiments.closure_ledger.quark_beta_boundary_probe import run_probe
    summary = run_probe()
    rate = summary["delta_invariance_rate_under_documented_drifts"]
    # 4 of 12 logged ablations sit at δ = +1.
    assert math.isclose(rate, 4.0 / 12.0, abs_tol=0.01), (
        f"δ-invariance rate {rate:.3f} drifted from the documented "
        f"4/12 = 0.333; the §8 ablation table may have changed."
    )


# --- Quark β: structural decomposition sub-probe -----------------------

def test_quark_beta_decomposition_probe_runs_to_completion():
    """The decomposition probe builds a structured catalog of fits."""
    from experiments.closure_ledger.quark_beta_decomposition_probe import run_probe
    summary = run_probe()
    assert summary["targets"]["m_lepton"] == 20
    assert summary["targets"]["m_quark"] == 93
    assert summary["targets"]["m_delta"] == 73
    assert summary["n_decompositions"]["m_lepton"] >= 1
    assert summary["n_decompositions"]["m_quark"] >= 1


def test_quark_beta_decomposition_lepton_is_single_block():
    """
    The simplest lepton decomposition is `(k_5-1)·k_5 = 20` — a single
    structural block at complexity 1.
    """
    from experiments.closure_ledger.quark_beta_decomposition_probe import run_probe
    summary = run_probe()
    best = summary["best_decompositions"]["m_lepton"]
    assert best is not None
    assert best["nonzero_blocks_count"] == 1
    assert best["coeff_complexity"] == 1
    assert best["block_coeffs"]["(k_5-1)·k_5"] == 1


def test_quark_beta_decomposition_quark_contains_lepton_subpiece():
    """
    The cleanest joint reading places the lepton block as a strict
    sub-decomposition of the quark formula, mirroring the
    'minimal closure ⊂ shell-coupled closure' framing in §1.
    """
    from experiments.closure_ledger.quark_beta_decomposition_probe import run_probe
    summary = run_probe()
    pairs = summary["joint_sub_decomposition_pairs"]
    assert pairs, "no lepton ⊆ quark decomposition pair found"
    top = pairs[0]
    # The leading pair must include the (k_5-1)·k_5 lepton block.
    assert "(k_5-1)·k_5" in top["lepton_quark_overlap_blocks"]


def test_quark_beta_decomposition_boundary_origin_is_color_residue():
    """
    The C_color_residue_N_c_minus_2 candidate scores 3/3 on the three
    boundary-origin tests (predicts +1 for quarks, 0 for leptons,
    drift-invariant). The other two candidates score lower.
    """
    from experiments.closure_ledger.quark_beta_decomposition_probe import run_probe
    summary = run_probe()
    best = summary["best_boundary_origin"]
    assert best["name"] == "C_color_residue_N_c_minus_2"
    assert best["consistency_score"] == 3
    by_name = {o["name"]: o for o in summary["boundary_origin_analysis"]}
    assert by_name["A_Z2_partition_residue"]["consistency_score"] < 3
    assert by_name["B_l_zero_s_wave_closure"]["consistency_score"] < 3


def test_quark_beta_decomposition_full_identity_holds():
    """
    The full structural identity must reproduce N_quark = 466 exactly:
        N_q = ((k_5-1)·k_5 + 2·k_5·(k_5+2) + N_c) · k_5 + (N_c - 2)
            = 93 · 5 + 1
            = 466
    """
    k_5 = 5
    n_c = 3
    m_q = (k_5 - 1) * k_5 + 2 * k_5 * (k_5 + 2) + n_c
    delta_q = n_c - 2
    n_q = m_q * k_5 + delta_q
    assert m_q == 93
    assert delta_q == 1
    assert n_q == 466
    # Lepton parallel: m_l = (k_5-1)·k_5, δ_l = 0.
    n_c_lepton = 1   # colorless
    m_l = (k_5 - 1) * k_5
    delta_l = max(0, n_c_lepton - 2)   # spinor factorization yields 0
    assert m_l == 20
    assert delta_l == 0
    assert m_l * k_5 + delta_l == 100


# --- Quark β: §8 robustness audit --------------------------------------

def test_quark_beta_robustness_audit_runs_to_completion():
    """The robustness audit produces a structured verdict over §8 data."""
    from experiments.closure_ledger.quark_beta_robustness_audit import run_probe
    summary = run_probe()
    assert summary["k_5"] == 5
    assert summary["baseline_N"] == 466
    assert len(summary["ablation_decompositions"]) == 12
    assert "claim_verdicts" in summary
    assert summary["n_claims_true"] + summary["n_claims_partial"] + \
        summary["n_claims_false"] == len(summary["claim_verdicts"])


def test_quark_beta_robustness_audit_baseline_and_uniform_scale_pass():
    """
    Claim 1 must hold: the +1 boundary correction is preserved under
    the no-physics-change perturbations (baseline + uniform scale).
    This is the only structural claim that survives the audit.
    """
    from experiments.closure_ledger.quark_beta_robustness_audit import run_probe
    summary = run_probe()
    by_name = {c["name"]: c for c in summary["claim_verdicts"]}
    claim1 = by_name["Claim 1: δ = +1 in baseline / uniform-scale runs"]
    assert claim1["verdict"] == "TRUE"


def test_quark_beta_robustness_audit_perturbation_claims_fail():
    """
    Claims 2, 3, 4 must FAIL: under per-species and anchor perturbations,
    δ wanders across {-1, 0, +1, +2}, contradicting the prior probes'
    structural-invariance reading. This is the audit's headline
    correction.
    """
    from experiments.closure_ledger.quark_beta_robustness_audit import run_probe
    summary = run_probe()
    by_name = {c["name"]: c for c in summary["claim_verdicts"]}
    assert by_name[
        "Claim 2: δ = +1 across single-species perturbations"
    ]["verdict"] == "FALSE"
    assert by_name[
        "Claim 3: δ = +1 across anchor-species changes"
    ]["verdict"] == "FALSE"
    assert by_name[
        "Claim 4: overall δ-invariance rate ≥ 50%"
    ]["verdict"] == "FALSE"


def test_quark_beta_robustness_audit_plus_one_rate_is_one_third():
    """The audit's headline rate is ~1/3 (4 of 12 logged ablations)."""
    from experiments.closure_ledger.quark_beta_robustness_audit import run_probe
    summary = run_probe()
    assert math.isclose(summary["plus_one_rate"], 4.0 / 12.0, abs_tol=0.01)


# --- Quark β: sub-block stability -------------------------------------

def test_quark_beta_subblock_stability_runs_to_completion():
    """The sub-block stability probe runs end-to-end."""
    from experiments.closure_ledger.quark_beta_subblock_stability import run_probe
    summary = run_probe()
    assert "subblock_tests" in summary
    assert "modular_invariant_tests" in summary
    assert "preserved_modular_invariants" in summary
    assert summary["baseline_decomposition"]["total"] == 466


def test_quark_beta_subblock_no_constant_shift_reduces_variance():
    """
    Sanity check: every constant-shift subtraction of a baseline sub-block
    leaves the residual width invariant. (Proof: Var(N − c) = Var(N).)
    The probe must surface this — if it ever shows a strict reduction,
    the underlying arithmetic is wrong.
    """
    from experiments.closure_ledger.quark_beta_subblock_stability import run_probe
    summary = run_probe()
    widths = {t["residual_width"] for t in summary["subblock_tests"]}
    assert len(widths) == 1, (
        f"sub-block subtractions gave varying widths {widths}; expected "
        "all equal under constant-shift invariance"
    )


def test_quark_beta_subblock_only_invariant_is_parity():
    """
    The only modular invariant preserved across all 12 §8 ablations is
    `N_q ≡ 0 (mod 2)`. Other natural moduli (3, 4, 5, 10, 15, 25, 100)
    fail.
    """
    from experiments.closure_ledger.quark_beta_subblock_stability import run_probe
    summary = run_probe()
    invariants = summary["preserved_modular_invariants"]
    assert len(invariants) == 1
    assert invariants[0]["modulus"] == 2
    assert invariants[0]["baseline_residue"] == 0


def test_quark_beta_subblock_all_logged_N_values_are_even():
    """
    The Z₂ partition-class invariance: every logged §8 N value is even.
    This is the structural foundation of the N_q = 2·n_part reading.
    """
    from experiments.closure_ledger.quark_beta_subblock_stability import run_probe
    summary = run_probe()
    half = summary["half_partition_diagnostic"]
    assert half["all_N_even"] is True
    assert half["half_baseline"] == 233   # 466 / 2


# --- Closure-cycle action quantum probe (ℏ-origin sub-target #1) -------

def test_closure_cycle_action_probe_runs_to_completion():
    """The action-quantum probe builds a structured summary."""
    from experiments.closure_ledger.closure_cycle_action_probe import run_probe
    summary = run_probe()
    assert summary["p1_status"] in {"PASS", "FAIL"}
    assert "n_quanta_per_species" in summary
    assert summary["all_species_integer_quantized"] is True
    assert "hbar_conversion_check" in summary
    assert "layer_2_effect_analysis" in summary


def test_closure_cycle_p1_pass_with_specific_integer_counts():
    """
    Headline finding: each species' Layer-1 ledger sum is an integer
    multiple of 2π under the locked baseline, with the specific
    counts (electron=2, muon=4, tau=106) reading off as
    (k+1) closure passes plus the 100·(2π) τ-uplift quantum.
    """
    from experiments.closure_ledger.closure_cycle_action_probe import run_probe
    summary = run_probe()
    assert summary["p1_status"] == "PASS"
    assert summary["n_quanta_per_species"] == {
        "electron": 2,
        "muon": 4,
        "tau": 106,
    }


def test_closure_cycle_layer_2_candidates_break_integer_quantization():
    """
    Sanity check on the closure-ledger sweep's Layer-2 candidates: C1
    and D1 each shift the per-species φ/2π by non-integer amounts,
    confirming that a Layer-2 closure form has to contribute exactly
    an integer multiple of 2π per species (not just be universal mod
    2π) to preserve the action-quantum reading.
    """
    from experiments.closure_ledger.closure_cycle_action_probe import run_probe
    summary = run_probe()
    by_name = {e["layer_2_candidate"]: e for e in summary["layer_2_effect_analysis"]}
    assert by_name[
        "Layer 2 absent (ledger universal at 0 mod 2π)"
    ]["preserves_integer_quantization"] is True
    assert by_name[
        "C1 (eigenvector-weighted B1 modes)"
    ]["preserves_integer_quantization"] is False
    assert by_name[
        "D1 (operator-valued V_j-V_i Hermitian matrix element)"
    ]["preserves_integer_quantization"] is False


def test_closure_cycle_R_MID_self_consistency_under_compton_identification():
    """
    Under R_MID = ℏ/(m_e c), the closure cycle in Compton wavelengths
    equals N · 2π by construction. The check verifies this self-
    consistency (it does NOT predict ℏ; it just records the
    convention's tautology).
    """
    from experiments.closure_ledger.closure_cycle_action_probe import run_probe
    summary = run_probe()
    for c in summary["hbar_conversion_check"]:
        assert c["self_consistent"] is True
        # cycle length in Compton wavelengths == N · 2π
        expected = c["n_quanta"] * 2.0 * math.pi
        assert math.isclose(
            c["closure_cycle_physical_compton_wavelengths"],
            expected,
            rel_tol=1e-12,
        )


# --- Closed-orbit radial action probe (ℏ-origin sub-target #1) --------

def test_closed_orbit_radial_action_probe_runs_to_completion():
    """The closed-orbit radial action probe builds a structured summary."""
    from experiments.closure_ledger.closed_orbit_radial_action_probe import run_probe
    summary = run_probe()
    assert "radial_action_table" in summary
    assert "species_couplings" in summary
    assert "exact_quantum_reading" in summary
    # Coverage: 3 l-values × 4 n-values = 12 entries.
    assert len(summary["radial_action_table"]) == 12


def test_closed_orbit_integer_quantization_holds_at_high_n():
    """
    Headline P1-at-exact-level finding: S_full(l, n) / 2π → (n + 1)
    integer at high n. Numerically, max deviation at n ≥ 2 must be
    < 0.05 (essentially exact for excited states).
    """
    from experiments.closure_ledger.closed_orbit_radial_action_probe import run_probe
    summary = run_probe()
    assert summary["wkb_to_exact_max_dev_at_n_ge_2"] < 0.05


def test_closed_orbit_b2_radial_ladder_integer_counts():
    """
    Under the B2_radial_ladder coupling (l=1, n=(k-1)/2), each species
    lands at exactly one bound mode that gives an integer total cycle:
        electron + (1, 0)  →  N_total_exact = 3
        muon     + (1, 1)  →  N_total_exact = 6
        tau      + (1, 2)  →  N_total_exact = 109
    The tau and muon WKB deviations should be < 0.01 (essentially
    exact for excited states); the electron's is larger (0.118 at
    the ground state) but the EXACT integer is still 3.
    """
    from experiments.closure_ledger.closed_orbit_radial_action_probe import run_probe
    summary = run_probe()
    rows = summary["species_couplings"]["B2_radial_ladder (l=1, n=(k-1)/2)"]
    by_label = {r["label"]: r for r in rows}
    assert by_label["electron"]["n_total_exact"] == 3
    assert by_label["muon"]["n_total_exact"] == 6
    assert by_label["tau"]["n_total_exact"] == 109
    # WKB convergence is sharp for the excited-state couplings.
    assert by_label["muon"]["wkb_deviation_from_exact"] < 0.01
    assert by_label["tau"]["wkb_deviation_from_exact"] < 0.01


def test_closed_orbit_wkb_error_explains_c1_residue():
    """
    The earlier closure-ledger probes' residues at the ground state
    (l, n=0) match the WKB-to-exact deviation: C1 reported residues
    0.760π–0.882π, and our probe shows the WKB deviations at n=0 are
    0.118, 0.229, 0.240 (1 − 0.882, 1 − 0.771, 1 − 0.760). The
    earlier failures are now isolated to a WKB approximation issue.
    """
    from experiments.closure_ledger.closed_orbit_radial_action_probe import run_probe
    summary = run_probe()
    table = summary["radial_action_table"]
    by_ln = {(r["l"], r["n"]): r for r in table}
    # Match the previously-recorded WKB single-pass values to ~3 sig figs.
    assert math.isclose(
        by_ln[(1, 0)]["s_single_wkb"] / math.pi, 0.882, abs_tol=0.001
    )
    assert math.isclose(
        by_ln[(3, 0)]["s_single_wkb"] / math.pi, 0.771, abs_tol=0.001
    )
    assert math.isclose(
        by_ln[(5, 0)]["s_single_wkb"] / math.pi, 0.760, abs_tol=0.001
    )


# --- Hard-wall boundary verification (ℏ-origin sub-target #1c) --------

def test_hard_wall_boundary_verification_runs_to_completion():
    """The verification probe builds a structured summary."""
    from experiments.closure_ledger.hard_wall_boundary_verification import run_probe
    summary = run_probe()
    assert "endpoint_checks" in summary
    assert "boundary_hypothesis_comparison" in summary
    assert "t_fixed_point_argument" in summary
    assert summary["dirichlet_verified_numerically"] is True


def test_hard_wall_endpoints_vanish_to_machine_precision():
    """
    All eigenfunctions u(r) returned by `solve_radial_modes` vanish at
    BOTH grid endpoints (inner = R_MID + 5e-4, outer = R_OUTER − 5e-4)
    to machine precision. This is the numerical Dirichlet condition,
    a consequence of the eigensolver's [1:N, 1:N] interior-point slice.
    """
    from experiments.closure_ledger.hard_wall_boundary_verification import run_probe
    summary = run_probe()
    for e in summary["endpoint_checks"]:
        assert e["inner_is_zero"], f"inner not zero: {e}"
        assert e["outer_is_zero"], f"outer not zero: {e}"


def test_hard_wall_DD_hypothesis_fits_observed_at_high_n():
    """
    The DD (Dirichlet+Dirichlet) hypothesis has the smallest deviation
    from observed closed-orbit action at n ≥ 2. Alternative BCs (DN,
    ND, NN, soft+soft) are decisively rejected.
    """
    from experiments.closure_ledger.hard_wall_boundary_verification import run_probe
    summary = run_probe()
    by_name = {h["name"]: h for h in summary["boundary_hypothesis_comparison"]}
    dd = by_name["DD_dirichlet_both"]
    dn = by_name["DN_dirichlet_inner_neumann_outer"]
    nn = by_name["NN_neumann_both"]
    soft = by_name["soft_both_standard_bs"]
    # DD should fit at high n.
    assert dd["deviation_max_at_high_n"] < 0.01
    # All other BCs should have substantial deviation at high n.
    assert dn["deviation_max_at_high_n"] > 0.1
    assert nn["deviation_max_at_high_n"] > 0.1
    assert soft["deviation_max_at_high_n"] > 0.1


def test_hard_wall_t_fixed_point_forces_dirichlet():
    """
    T = iσ_y has T² = −I, and ψ = T·ψ admits only ψ = 0 as solution.
    This is the topological argument that forces Dirichlet at the
    throat — Dirichlet inner is physical, not just numerical.
    """
    from experiments.closure_ledger.hard_wall_boundary_verification import run_probe
    summary = run_probe()
    t_arg = summary["t_fixed_point_argument"]
    assert t_arg["T_squared_eq_minus_I"] is True
    assert t_arg["only_zero_solution_to_psi_eq_T_psi"] is True


def test_hard_wall_best_hypothesis_is_DD():
    """The data-best boundary hypothesis must be DD."""
    from experiments.closure_ledger.hard_wall_boundary_verification import run_probe
    summary = run_probe()
    assert summary["best_boundary_hypothesis"]["name"] == "DD_dirichlet_both"


def test_geometric_hamiltonian_locked_surrogate_matches_locked_lepton_eigenvectors():
    """
    The probe's reference eigensystem must agree with the existing C1
    eigenvectors used by the closure ledger. This catches accidental
    drift between the two computation paths.
    """
    import numpy as np
    from experiments.closure_ledger.geometric_hamiltonian_probe import (
        _locked_surrogate_eigensystem,
    )
    from experiments.closure_ledger.sk_bridge import (
        _locked_lepton_eigenvectors,
    )
    evals_probe, evecs_probe = _locked_surrogate_eigensystem()
    evals_sk, evecs_sk = _locked_lepton_eigenvectors()
    for a, b in zip(evals_probe, evals_sk):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9)
    # Eigenvectors are sign-ambiguous; compare |v|² element-wise.
    for vp, vs in zip(evecs_probe, evecs_sk):
        for p, s in zip(vp, vs):
            assert math.isclose(p ** 2, s ** 2, abs_tol=1e-9)
