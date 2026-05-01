"""
Layer 2 of the closure-phase ledger: structured blocker report.

This module emits a structured description of what the closure ledger
cannot compute, given the current state of the repo. Per the repo
audit:

- `lepton_spectrum.py` is an instanton-transition surrogate operating
  on depth labels k ∈ {1, 3, 5}. The radial-eigenmode parameters
  (`l`, `n_points`, `rs`, `r_outer`) are accepted for API compatibility
  but unused in `compute_knotted_lepton_spectrum`.

- The bridge from generation depth k to Tangherlini eigenmodes (l, n)
  — call it S(k) — does not yet exist in either sector.

This means:

- Φ_radial(k) is structurally blocked for the lepton sector.
- The 366-quanta gap prediction is downgraded to a hypothesis.
- The full closure-phase invariant cannot be computed end-to-end.

Defining S(k) is now the next theoretical target, not the next
implementation task. Three principled candidate maps are recorded
below for reference.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional


@dataclass
class CandidateSkMap:
    """A principled candidate definition of S(k) → {(l, n) modes}."""
    name: str
    formula: str
    physical_picture: str
    advantages: str
    open_questions: str
    implementation_status: str = "open"   # "implemented" | "open" | "deferred"

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class SkBridgeBlocker:
    """
    Structured Layer-2 status report for the S(k) → Φ_radial bridge.

    Originally a pure blocker (Layer 2 unimplemented). Now parameterized
    by `implemented_candidate`: when one of the candidates is wired, the
    verdict reports the wiring and the candidate's `implementation_status`
    flips to "implemented"; the other candidates continue to live here as
    open thesis-level alternatives. The downgraded P3 prediction is
    retained because the quark sector still requires its own S(k) bridge.
    """
    verdict: str
    evidence: list[str] = field(default_factory=list)
    candidates: list[CandidateSkMap] = field(default_factory=list)
    next_steps: list[str] = field(default_factory=list)
    downgraded_predictions: list[str] = field(default_factory=list)
    implemented_candidate: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "verdict": self.verdict,
            "evidence": list(self.evidence),
            "candidates": [c.to_dict() for c in self.candidates],
            "next_steps": list(self.next_steps),
            "downgraded_predictions": list(self.downgraded_predictions),
            "implemented_candidate": self.implemented_candidate,
        }


def layer2_blocker_report(
    implemented_candidate: Optional[str] = None,
) -> SkBridgeBlocker:
    """
    Construct the Layer-2 status report.

    If `implemented_candidate` matches one of the candidate names, that
    candidate is marked as implemented and the verdict reports the wiring;
    other candidates remain open thesis-level alternatives. Pass None or
    "none" to get the original "Layer 2 blocked" verdict.
    """
    candidates = [
        CandidateSkMap(
            name="A_lowest_radial_per_l",
            formula="S(k) = { (l, n=0) : l = 1, 3, ..., k }",
            physical_picture=(
                "Generation k is the sum of odd-l radial ground states up "
                "to angular harmonic l = k. Ties odd-k closure to odd-l "
                "angular content."
            ),
            advantages=(
                "Consistent with the non-orientable transport rule: "
                "closure must flip the Z₂ partition class, and odd-l "
                "modes are the natural candidates for partition-flipping "
                "angular harmonics."
            ),
            open_questions=(
                "Result with WKB radial-action convention: FALSIFIES "
                "universal closure mod 2π. The cumulative-odd-ground-mode "
                "interpretation of lepton depth is rejected; whether a "
                "different phase convention (Maslov, Bohr-Sommerfeld with "
                "two soft turning points) revives it is open."
            ),
        ),
        CandidateSkMap(
            name="B1_single_angular_mode",
            formula="S(k) = { (l = k, n = 0) }",
            physical_picture=(
                "Generation k is a single angular harmonic (l = k) in its "
                "radial ground state. One mode per generation, indexed by "
                "the angular quantum number alone."
            ),
            advantages=(
                "Cleanest single-mode-per-generation interpretation: one "
                "angular eigenstate per lepton family. Easy to falsify "
                "because the per-row Φ is a single integrated quantity."
            ),
            open_questions=(
                "Does Φ(l = k, n = 0) close to a universal value mod 2π "
                "under the chosen WKB convention? The B1 phases are the "
                "same numbers that already appear in the candidate-A "
                "decomposition, so candidate-A's failure mode constrains "
                "B1 directly."
            ),
        ),
        CandidateSkMap(
            name="B2_single_radial_excitation",
            formula="S(k) = { (l = 1, n = (k − 1) / 2) }",
            physical_picture=(
                "All generations share the lowest angular harmonic l = 1; "
                "depth labels successive radial excitations n = 0, 1, 2. "
                "Lepton depth is reinterpreted as radial-mode quantum "
                "number rather than angular content."
            ),
            advantages=(
                "Uses one fixed angular sector and a single radial ladder, "
                "matching the surrogate's `β · k²` uplift if and only if "
                "ω²(1, n) ≈ ω²(1, 0) + β · n² (testable directly). "
                "Maximally falsifiable: a single mode integral per row."
            ),
            open_questions=(
                "Does Φ(l = 1, n) close to a universal value mod 2π? "
                "Asymptotically Φ(1, n) → (n + 1) π by WKB, so for high n "
                "the residues converge to a parity pattern {π, 0, π, 0, …} "
                "rather than a single universal value — universality "
                "requires a Maslov-shifted convention or a convention "
                "that absorbs the (n + 1) π structure."
            ),
        ),
        CandidateSkMap(
            name="C1_eigenvector_weighted_B1",
            formula=(
                "Φ_radial(k) = Σ_i |v_species(k),i|² · Φ(l = k_i, n = 0), "
                "weights from the locked lepton generation block"
            ),
            physical_picture=(
                "The depth basis {1, 3, 5} is shared with the instanton "
                "surrogate; each species' radial phase is the squared-"
                "amplitude weighted sum of the B1 ground modes (l = k_i, "
                "n = 0) over the depth basis. The B1 hand-imposed "
                "single-mode bridge is the |v|² = δ_ij limit of this map."
            ),
            advantages=(
                "Derives weights from the existing lepton Hamiltonian "
                "without introducing fitted parameters. The eigenvector "
                "mixing lifts B1's degeneracy and produces tighter "
                "residues than any prior hand-imposed candidate."
            ),
            open_questions=(
                "Does the lifting close the residues to a single value "
                "mod 2π? If not, the eigenvector mixing alone does not "
                "supply the missing bridge under the WKB convention; "
                "test whether a different phase convention (Maslov, "
                "Bohr-Sommerfeld) revives universality given the "
                "tightened spread."
            ),
        ),
        CandidateSkMap(
            name="C2_eigenvector_weighted_B2",
            formula=(
                "Φ_radial(k) = Σ_i |v_species(k),i|² · Φ(l = 1, n = (k_i−1)/2), "
                "weights from the locked lepton generation block"
            ),
            physical_picture=(
                "Same eigenvector weights as C1 but selecting B2's "
                "single-l radial-excitation ladder. Each species mixes "
                "across radial excitations n ∈ {0, 1, 2} of l = 1."
            ),
            advantages=(
                "Tests whether the eigenvector mixing of B2's ladder "
                "(which under WKB asymptotes to (n+1)π) collapses the "
                "parity pattern across generations to a single value."
            ),
            open_questions=(
                "Same convention-dependence as C1; weight rows are the "
                "same in both candidates, only the per-mode Φ values "
                "differ. The two are independent tests of which mode "
                "ladder (B1 angular or B2 radial) the lepton "
                "Hamiltonian's eigenvectors actually align with."
            ),
        ),
        CandidateSkMap(
            name="C1_maslov_standard",
            formula=(
                "Φ_radial(k) = Σ_i |v_species(k),i|² · "
                "[ ∫√max(ω² − V_eff, 0) dr* − (π/2)·N_turning(l_i) ], "
                "S(k) modes as in C1, Maslov shift = −π/2 per detected "
                "classical turning point"
            ),
            physical_picture=(
                "Same eigenvector-weighted B1 mode set as C1, but each "
                "mode integral is corrected by the standard Bohr-Sommerfeld "
                "Maslov shift: −π/2 for every sign change of (ω² − V_eff) "
                "inside the tortoise grid. Hard-wall grid endpoints (where "
                "ω² > V_eff at the boundary) contribute no Maslov phase."
            ),
            advantages=(
                "Adds the textbook turning-point correction without "
                "introducing any fitted parameters. Tests whether the "
                "WKB-bare convention is the obstruction to closure: if "
                "the correction is uniform across species (e.g. every B1 "
                "mode has the same turning-point count) the spread is "
                "preserved and only the universal residue shifts; if the "
                "count varies across modes, the correction can in "
                "principle redistribute residues."
            ),
            open_questions=(
                "Empirical result: every B1 ground mode reports "
                "N_turning = 1 on the canonical tortoise grid, so the "
                "correction is a uniform −π/2 across the depth basis. "
                "Spread is preserved relative to C1; residues translate "
                "but do not collapse. Universality is not restored by "
                "the standard Maslov correction on this mode set."
            ),
        ),
        CandidateSkMap(
            name="B2_maslov_standard",
            formula=(
                "Φ_radial(k) = ∫√max(ω² − V_eff(l=1), 0) dr* "
                "− (π/2) · N_turning(l=1, n=(k−1)/2), "
                "S(k) modes as in B2"
            ),
            physical_picture=(
                "Same single-mode l=1 ladder as B2, with the standard "
                "Bohr-Sommerfeld Maslov shift (−π/2 per detected sign "
                "change of ω² − V_eff) applied to each mode integral."
            ),
            advantages=(
                "Tests whether the per-mode turning-point count varies "
                "across the radial ladder. The (l=1, n=0) ground mode "
                "sits below the centrifugal+grid-bounded barrier and has "
                "a soft inner turning point; higher excitations (l=1, "
                "n≥1) sit above the barrier and may have N_turning = 0. "
                "If so, the Maslov shift is differential across "
                "generations — exactly the structure that could "
                "redistribute residues mod 2π."
            ),
            open_questions=(
                "Does the differential Maslov correction collapse the "
                "(n+1)π parity pattern of B2 to a universal value, or "
                "does it merely permute the residues? The radial-detail "
                "table reports n_turning_points per mode so the regime "
                "can be read off directly."
            ),
        ),
        CandidateSkMap(
            name="C2_maslov_standard",
            formula=(
                "Φ_radial(k) = Σ_i |v_species(k),i|² · "
                "[ ∫√max(ω² − V_eff(l=1), 0) dr* "
                "− (π/2)·N_turning(l=1, n=i) ], "
                "S(k) modes as in C2"
            ),
            physical_picture=(
                "C2's eigenvector-weighted B2 ladder with the standard "
                "Maslov shift applied per mode. Each species' radial "
                "phase is the squared-amplitude weighted sum over the "
                "B2 ladder of (Φ_WKB − (π/2)·N_turning)."
            ),
            advantages=(
                "Combines the two structural levers that have any chance "
                "of breaking the ledger non-universality: eigenvector "
                "mixing across the depth basis (C-family) and a per-mode "
                "turning-point-dependent phase (Maslov). If the B2 modes "
                "have non-uniform N_turning, the eigenvector weights "
                "redistribute the differential Maslov shift across "
                "species."
            ),
            open_questions=(
                "Whether either the C-family weighting alone, the Maslov "
                "shift alone, or their composition supplies the missing "
                "bridge. C1_maslov_standard already showed the composition "
                "with the B1 mode set is degenerate (uniform N_turning); "
                "this candidate is the natural follow-on test on the B2 "
                "ladder where N_turning is expected to vary."
            ),
        ),
    ]

    impl = implemented_candidate if implemented_candidate not in (None, "none") else None
    if impl is not None:
        for c in candidates:
            if c.name == impl:
                c.implementation_status = "implemented"
        verdict = (
            f"Φ_radial(k) is wired in the lepton sector via candidate "
            f"`{impl}` using `geometrodynamics.tangherlini.radial.solve_radial_modes` "
            f"plus a Bohr-Sommerfeld integration of √max(ω² − V_eff, 0) "
            f"over the tortoise grid. The remaining candidates are kept "
            f"as open thesis-level alternatives. The quark-sector S(k) "
            f"bridge is still undefined, so the 366-quanta gap remains "
            f"S(k)-conditional."
        )
    else:
        verdict = (
            "Φ_radial(k) is structurally blocked for the lepton sector: "
            "lepton_spectrum.py operates on depth labels k ∈ {1, 3, 5} "
            "via an instanton-transition surrogate, with no map from "
            "generation depth to Tangherlini eigenmodes (l, n)."
        )

    return SkBridgeBlocker(
        verdict=verdict,
        implemented_candidate=impl,
        evidence=[
            "compute_knotted_lepton_spectrum accepts `l`, `n_points`, "
            "`rs`, `r_outer` for API compatibility but explicitly does "
            "not use them in the surrogate Hamiltonian.",
            "The locked lepton diagonal H_kk = action_base + "
            "resistance_scale·k² + res_diag(k) + pinhole(k∈{3,5}) + "
            "β·max(0, k−3)² has no eigenmode index and no Tangherlini "
            "potential evaluation.",
            "The quark residual sector (transport, pinhole, resistance) "
            "DOES read scalars off the tortoise grid — but those are "
            "single integrated quantities, not a per-generation "
            "S(k) → (l, n) mapping.",
        ],
        next_steps=(
            [
                "Read the universality_check on the lepton ledger: candidate "
                f"`{impl}` either preserves closure mod 2π (supports the "
                "candidate) or breaks it (falsifies the candidate; rerun with "
                "another wired candidate).",
                "Define a quark-sector S(k) bridge. The lepton implementation "
                "operates on (l, n) Tangherlini modes; the quark sector needs "
                "its own per-generation mode assignment before P3 (the "
                "366-quanta gap) becomes testable end-to-end.",
                "Investigate whether the candidate's per-mode Φ(l, n) values "
                "match the surrogate's β·k² uplift coefficients — this is the "
                "open-question line in candidate A's record.",
                "Update THESIS.md to reflect whichever way the universality "
                "check went: either close the ℏ open-problem bullet to a "
                "derived dimensionless invariant, or report which channel "
                "was overcredited.",
            ]
            if impl is not None
            else [
                "Define S(k) for the lepton sector. Three candidates listed "
                "in `candidates`; selection is thesis-level, not "
                "implementation-level.",
                "Once S(k) is fixed, the radial phase extractor becomes "
                "well-posed: Φ_radial(k) = Σ_(l,n)∈S(k) ∫ k_local(r*) dr*.",
                "Test the three falsification predictions (P1 lepton "
                "universality, P2 quark universality, P3 the 366-quanta gap) "
                "against the implemented S(k).",
                "Update THESIS.md to reflect the result: either close the "
                "ℏ open-problem bullet to a derived dimensionless invariant, "
                "or report which channel was overcredited.",
            ]
        ),
        downgraded_predictions=[
            "P3 (366-quanta gap as integrated bulk-mode contribution): "
            "downgraded from near-term falsification test to S(k)-"
            "conditional hypothesis. The number 366 = 466 − 100 is a "
            "structural prediction about what the bulk-coupling channel "
            "must contribute across the quark ladder relative to the "
            "lepton ladder, contingent on S(k) being defined for both.",
        ],
        candidates=candidates,
    )
