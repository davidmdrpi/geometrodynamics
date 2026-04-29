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
                "Does this map reproduce the locked lepton masses when "
                "Φ_radial(k) is added back into the diagonal Hamiltonian? "
                "If yes, the lepton instanton surrogate is a depth-block "
                "effective theory of this microscopic mode set."
            ),
        ),
        CandidateSkMap(
            name="B_fixed_total_quantum_number",
            formula="S(k) = { (l, n) : 2n + l = k }",
            physical_picture=(
                "Generation k is a single shell in the joint (l, n) "
                "spectrum, with each generation mapping to one mode "
                "rather than a sum."
            ),
            advantages=(
                "Clean single-mode-per-generation interpretation. Easier "
                "to falsify: each generation's Φ_radial value is a single "
                "eigenmode integral, with no sum-over-modes ambiguity."
            ),
            open_questions=(
                "Does the surrogate's β·k² term then read as the "
                "eigenfrequency ω(l, n) for the selected single mode? "
                "If so, this candidate gives a sharper microscopic "
                "interpretation than candidate A."
            ),
        ),
        CandidateSkMap(
            name="C_closure_coherent_superposition",
            formula="S(k) = { (l, n) : ω(l, n) satisfies S³ closure at depth k }",
            physical_picture=(
                "The mode set is determined by which (l, n) combinations "
                "actually close on the antipodal cavity at depth k — the "
                "S³ closure condition selects the membership of S(k)."
            ),
            advantages=(
                "Most physically principled: makes S(k) emergent from the "
                "antipodal closure condition rather than imposed. Connects "
                "channel 1 (closure) and channel 3 (eigenmodes) directly."
            ),
            open_questions=(
                "Requires the S³ closure phase condition to be spelled out "
                "as an explicit equation on ω(l, n) before the membership "
                "of S(k) becomes computable. This is the highest-leverage "
                "candidate but also the highest-effort to define."
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
