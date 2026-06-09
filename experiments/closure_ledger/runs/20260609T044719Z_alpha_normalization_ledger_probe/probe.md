# Alpha normalization ledger for the gauge–matter coupling (PR #143)

**Run:** 2026-06-09T04:47:19+00:00

Consolidating ledger for the EM coupling normalisation α (the strength left input by #141/#142) — parallel to the bulk-scale ledger (#133). Separates what the geometry derives (the charge quantum, the 1/2π loop measure, the coupling structure, the running) from the one irreducible input: the value α ≈ 1/137 (the 137 problem). *(QFT on the classical throat, not quantum gravity.)*

- **How α enters**: A_EM = α·ℏc/2; vertex ∝ c₁²α; a = α/2π — one number
- **Charge quantum**: DERIVED: |c₁| = 1 (integer Hopf number, #58/#74)
- **Measure**: DERIVED: the 1/2π in a = α/2π (closure quantum, #74)
- **Running**: DERIVED (vacuum polarisation, #142); value α input (#105)
- **Value**: the one EM residual α ≈ 1/137 (the 137 problem), no clean closure match (#108)
- **Input budget**: {n_part, √σ/m_e, ε, α} (#104/#108) — α the EM residual
- **Open**: the value α (137 problem); EM normalisation; bulk scale (#133); flavor (#134)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | α normalization ledger for the gauge–matter coupling | **PASS** |
| T2 | `T2_how_alpha_enters` | how α enters: A_EM = α·ℏc/2, vertex ∝ c₁²α, a = α/2π — one number | **PASS** |
| T3 | `T3_charge_quantum_derived` | charge quantum derived: |c₁| = 1 (integer Hopf number, #58/#74) | **PASS** |
| T4 | `T4_two_pi_loop_measure_derived` | 1/2π measure derived: a = α/2π, the 2π = closure quantum (#74) | **PASS** |
| T5 | `T5_running_derived_value_input` | running derived (#142); value α(μ_0) input (#105) | **PASS** |
| T6 | `T6_value_is_one_em_residual` | value α ≈ 1/137 = one residual; no clean closure match (137 problem) | **PASS** |
| T7 | `T7_ledger_and_input_budget` | ledger: derived (quantum/measure/structure/running) vs input (value α) | **PASS** |
| T8 | `T8_assessment` | ALPHA_NORMALIZATION_LEDGER_..._VALUE_ONE_RESIDUAL | **PASS** |

## No clean closure match for α⁻¹ = 137.036 (the 137 problem)

| candidate | value | % off | needs ad-hoc term? |
|---|---:|---:|:---:|
| 2π | 6.283 | -95.41% | — |
| β_lepton = 50π | 157.08 | 14.63% | — |
| 2π·k₅² = 50π | 157.08 | 14.63% | — |
| k₅³ + 2π | 131.283 | -4.2% | — |
| 50π − 20 | 137.08 | 0.03% | ✗ (fit) |
| 4·k₅² + 37 | 137.0 | -0.03% | ✗ (fit) |
| 8π·k₅ | 125.664 | -8.3% | — |

The sub-% near-misses (`50π − 20`, `4·k₅² + 37`) each require an ad-hoc additive `O(20–37)` integer — fits, not derivations (the #107/#108 failure mode). No clean closure number lands near 137: α is plausibly irreducible, like √σ/m_e (#108).

## Verdict

**ALPHA_NORMALIZATION_LEDGER_CHARGE_QUANTUM_AND_2PI_MEASURE_DERIVED_VALUE_ONE_RESIDUAL.** THE EM COUPLING'S CHARGE QUANTUM, 1/2π MEASURE, STRUCTURE, AND RUNNING ARE DERIVED; ONLY THE VALUE α ≈ 1/137 IS INPUT — ONE EM RESIDUAL, THE 137 PROBLEM. PRs #141/#142 derived the gauge–matter coupling structure but left the strength α input; this ledger consolidates what α is, parallel to the bulk-scale ledger (#133).

HOW α ENTERS. Every EM observable is a function of α and the geometry: the amplitude A_EM = α·ℏc/2 (#105), the gauge–matter vertex strength ∝ c₁²α (#141), the Schwinger anomaly a = α/2π (#74). α is a single dimensionless number feeding the EM sector.

THE CHARGE QUANTUM IS DERIVED. The Hopf charge is the integer Hopf number |c₁| = 1 (#58/#74) — charge quantisation is topological, not input. The charge unit is geometric; only the strength is α.

THE 1/2π MEASURE IS DERIVED. In a = α/2π, the 2π is the BAM closure-quantum loop measure (#74) — derived. The geometry fixes the 1/2π and leaves only α as the input prefactor: BAM derives the measure, not the coupling.

THE RUNNING IS DERIVED; THE VALUE IS NOT. The RG flow of α (the vacuum polarisation, transverse by the #142 Ward identity) is derived structurally — BAM derives HOW α runs. The boundary value α(μ_0) ≈ 1/137 is the input — BAM does not derive WHERE it starts (the #105 classification, sharpened).

THE VALUE IS THE ONE EM RESIDUAL. A fit-independent scan of α⁻¹ = 137.036 against the BAM closure numbers (2π, k₅, β_lepton = 50π) finds NO clean match: the apparent near-misses (50π − 20 = 137.08, 4·k₅² + 37 = 137.0) each need an ad-hoc additive O(20–37) term — fits, not derivations, the same reverse-engineering failure mode #107/#108 documented for √σ/m_e. So α is plausibly irreducible, like √σ/m_e (#108): the EM sector contributes exactly one dimensionless residual, the value α.

TIE TO THE INPUT BUDGET. α joins the program's handful of dimensionless residuals — {n_part, √σ/m_e, ε, α} (#104/#108) — as the EM one. The charge quantum, the 1/2π measure, the coupling structure, and the running are derived; the value α is the single EM input, the 137 problem.

SCOPE. A consolidating/accounting ledger. It separates the derived EM structure from the one input (the value α) and shows α has no clean closure match. It does NOT derive α (the 137 problem stays open) or the EM normalisation absolutely; the α (#105/#108), bulk-scale (#133), and flavor (#134) residuals stand.
