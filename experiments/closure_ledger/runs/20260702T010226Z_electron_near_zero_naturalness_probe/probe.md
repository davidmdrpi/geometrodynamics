# Attacking the fine-tuning: the electron near-zero — stabilized or dialed? (PR #194)

**Run:** 2026-07-02T01:02:26+00:00

The follow-up to #192/#193: is the surrogate's electron near-zero protected by a symmetry, an index, a seesaw, or an attractor — or genuinely dialed? Every candidate is tested. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | goal: stabilized or dialed? | **PASS** |
| T2 | `T2_anatomy` | anatomy: a 2.9% cancellation (diagonal vs repulsion) | **PASS** |
| T3 | `T3_symmetry_candidates_excluded` | chiral / reflection / index / seesaw — all EXCLUDED | **PASS** |
| T4 | `T4_bg_sensitivities` | BG sensitivities: only the near-zero is tuned | **PASS** |
| T5 | `T5_dialed_combination` | the ONE dialed combination + the #192 cross-check | **PASS** |
| T6 | `T6_monte_carlo_null` | Monte Carlo null: as rare as generic; no attractor | **PASS** |
| T7 | `T7_origin_and_guardrail` | origin: calibration imports the hierarchy; guardrail held | **PASS** |
| T8 | `T8_assessment` | DIALED, not stabilized; mechanism question sharply posed | **PASS** |

## Barbieri–Giudice sensitivities (only the near-zero is tuned)

| parameter | Δ(E_e) | Δ(E_μ) | Δ(E_τ) |
|---|---:|---:|---:|
| `phase_per_pass` | 0.0 | -0.0 | -0.0 |
| `transport_strength` | -57.07 | 0.261 | 0.001 |
| `hard_pinhole_gamma` | 17.93 | 0.459 | 0.032 |
| `resistance_scale` | 7.43 | 0.127 | 0.054 |
| `action_base` | 31.47 | 0.152 | 0.009 |
| `action_slope` | 30.92 | -0.136 | -0.001 |
| `k_uplift_beta` | 1.24 | 0.001 | 0.904 |

## The dialed combination

- n̂ ∝ ∇ln E_e: transport **-0.764**, base action **0.421**, slope **0.414**, pinhole **0.24** (global Δ = **74.7**)
- fiber-map contraction: **+71.1** — reproduces the #192 Berger λ-sensitivity (+71): the two probes see the same dial

## The Monte Carlo null

- P(|E_e| ≤ 0.1996) = **0.0767** vs linear-measure estimate **0.06** — same order; distribution FLAT through zero (no attractor, no repulsion)
- P(μ/e ≥ 206.7) = **0.0399** — the hierarchy lives only in the tuned sliver

## Verdict

**ELECTRON_NEAR_ZERO_IS_UNPROTECTED_CANCELLATION_TUNING_IMPORTED_BY_CALIBRATION_NO_STABILIZING_SYMMETRY_FOUND.** DIALED, NOT STABILIZED — and the dial is identified.

EXCLUSIONS. The near-zero survives no protection candidate: no chiral/sublattice sign conjugation exists (tr H = 736 ≠ 0), the spectrum is not reflection-symmetric, the near-zero eigenvector aligns with no structured/topological direction, and the smallness FLIPS SIGN under a ±2% transport change — cancellation between comparable terms, not a sign-stable seesaw suppression.

THE DIAL. The zero locus is codimension 1 and the one tuned combination is measured: n̂ = 76% transport vs 42% base action + 41% slope + 24% pinhole, with global tuning Δ = 74.7. The Monte Carlo null (log-uniform ±25%, 20000 samples) gives P(|E_e| ≤ observed) = 0.077 — matching the naive linear-measure estimate 0.060 with a FLAT distribution through zero: no attractor, no repulsion, exactly as rare as generic. Cross-check: the dialed direction contracted with the #192 fiber map gives +71.1, reproducing the measured Berger λ-sensitivity — #192 and this probe see one and the same tuned combination.

ORIGIN AND UPSHOT. The tuning is imported by the calibration itself: fitting μ/e = 206.77 with an O(10)-scale matrix forces |E_e| = E_μ/206.77. The surrogate CARRIES the hierarchy problem rather than solving it (and the #193 operator has no near-zero at all — the tuning lives entirely in the instanton/transport dynamics, exactly where the dialed direction points). The mechanism question is now sharp: new structure pinning a k=1 zero mode — an index, a chiral grading of winding sectors, or a geometric identity tying the 2π base action to the transport repulsion — is what a real solution requires; finding or refuting it in the throat geometry is the follow-up. Numerology guardrail held: det-zero roots measured against round constants and NOT matched (transport root vs 8π: 1.6% — rejected).
