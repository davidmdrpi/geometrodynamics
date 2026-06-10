# Residual-bracket synthesis and input-budget ledger (PR #150)

**Run:** 2026-06-10T04:27:08+00:00

The synthesis capstone for the program's input accounting: one categorized budget assembled from #104–#149, a keystone re-verified from every category, and the no-loose-knobs claim made checkable — zero inputs added across #144–#149, the budget constant since #104/#125 while the derived ledger grew. *(QFT on the classical throat, not quantum gravity.)*

- **Budget**: 1 anchor + {α, √σ/m_e} + {n_part, ε} + 2 brackets + flavor puzzle
- **Anchor**: G → ΔR = 0.52·R_MID (B4-mandatory, #52/#53); √6 fixed (#57)
- **Scans**: best principled: α⁻¹ −4.2%, √σ/m_e −5.4%; sub-% all ad-hoc (#107/#108/#143)
- **Brackets**: k·r_s ∈ (0, 0.0064–0.070] (#148); ε₃/ε₂ ∈ [1.32, 1.44] (#149)
- **No new knobs**: #144–#149: six probes, zero inputs added
- **Contrast**: lepton N = 4k₅² = 100 fully derived — the no-residual sector (#124)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | one categorized input budget; no-loose-knobs checkable | **PASS** |
| T2 | `T2_anchor_sector` | anchor: ΔR/R_MID = 0.52; √6 fixed; exactly one (B4) | **PASS** |
| T3 | `T3_universal_residuals_scans` | α and √σ/m_e scans re-run: best principled −4.2% / −5.4% | **PASS** |
| T4 | `T4_program_residuals` | n_part recycling circular ([764, 920] vs fixed 830); lepton 100 derived | **PASS** |
| T5 | `T5_bracketed_subresiduals` | brackets re-checked: k·r_s c ≈ 9.9; ε₃/ε₂ = 1.435 | **PASS** |
| T6 | `T6_the_input_budget_ledger` | THE LEDGER: #144–#149 added zero inputs; budget = #104/#125 | **PASS** |
| T7 | `T7_scope` | organizes, does not remove; what would change the ledger | **PASS** |
| T8 | `T8_assessment` | INPUT_BUDGET_ONE_ANCHOR_FOUR_RESIDUALS_TWO_BRACKETED_ZERO_NEW_KNOBS | **PASS** |

## The categorized input budget

| category | item | status | source PRs |
|---|---|---|---|
| ANCHOR (dimensionful) | G (→ ΔR unit) | mandatory (B4), relocatable | #52/#53/#57/#106/#133 |
| FIXED TUNING | √6 (RS flatness) | derived constant, not a knob | #57 |
| UNIVERSAL RESIDUAL | α ≈ 1/137 | structure/running derived; value scan-excluded | #74/#141–#147; #143 |
| UNIVERSAL RESIDUAL | √σ/m_e ≈ 830 | one-G repackaging derived; value scan-excluded | #106; #107/#108 |
| PROGRAM RESIDUAL | n_part = 233 | doubling topological (APS); value compensator | #97/#123/#125 |
| PROGRAM RESIDUAL | ε (ν compliance) | order-of-mag derived; window [2π, k₅√(2π)] | #89/#112 |
| BRACKETED SUB-RESIDUAL | k·r_s | (0, 0.0064–0.070] two-sided | #133/#148 |
| BRACKETED SUB-RESIDUAL | ε_n spread | [1.32, 1.44]/step, ~0.3%; power laws excluded | #113/#149 |
| UNIVERSAL OPEN PROBLEM | flavor puzzle | RG-invariant ⟹ not running; no theory derives it | #97/#107/#108/#134 |
| NO RESIDUAL (contrast) | lepton N = 4k₅² = 100 | structure AND value derived | #124 |

## Zero new inputs across the recent arc

| PR | inputs added |
|---|---:|
| #144 vacuum polarisation / running | 0 |
| #145 Z₁ = Z₂ charge non-renormalization | 0 |
| #146 charge form factor / geometric radius | 0 |
| #147 F₁/F₂ EM-arc capstone | 0 |
| #148 k·r_s bracket | 0 |
| #149 ε_n bracket | 0 |

## The universal-residual scans (re-run)

| α⁻¹ candidate | value | % off | ad-hoc? |
|---|---:|---:|:---:|
| 2π | 6.283 | -95.41% | — |
| β_lepton = 50π | 157.08 | 14.63% | — |
| k₅³ + 2π | 131.283 | -4.2% | — |
| 8π·k₅ | 125.664 | -8.3% | — |
| 50π − 20 (ad-hoc) | 137.08 | 0.03% | ✗ |
| 4k₅² + 37 (ad-hoc) | 137 | -0.03% | ✗ |

| √σ/m_e candidate | value | % off | ad-hoc? |
|---|---:|---:|:---:|
| 2π·k₅³ = β_lepton·k₅ | 785.4 | -5.4% | — |
| k₅⁴ | 625.0 | -24.7% | — |
| e^{2π} | 535.5 | -35.5% | — |
| (4/3)·k₅⁴ (ad-hoc) | 833.3 | 0.4% | ✗ |

## The bracketed sub-residuals (keystones re-checked)

| residual | bracket | upper bound | lower bound | keystone re-check |
|---|---|---|---|---|
| k·r_s (AdS scale, #133) | (0, 0.0064–0.07] | locked spectrum (#148) | static throat (#56/#57) | sensitivity c = 9.86 (≈ 9.86) |
| ε_n spread (flavor, #113) | ε₃/ε₂ ∈ [1.32, 1.44] (~0.3%) | oscillation data ÷ p (#149) | attribution endpoint | pure-bounce ε₃/ε₂ = 1.435 |

## Verdict

**INPUT_BUDGET_ONE_ANCHOR_FOUR_RESIDUALS_TWO_BRACKETED_ZERO_NEW_KNOBS.** THE INPUT BUDGET IS CATEGORIZED, KEYSTONED, AND CONSTANT: ONE GRAVITATIONAL ANCHOR, TWO UNIVERSAL RESIDUALS (α, √σ/m_e — SCAN-EXCLUDED FROM CLEAN CLOSURE MATCHES), TWO PROGRAM RESIDUALS (n_part, ε — STRUCTURE DERIVED), TWO BRACKETED SUB-RESIDUALS (k·r_s, THE ε_n SPREAD — BOXED TWO-SIDED), AND THE UNIVERSAL FLAVOR PUZZLE — WITH ZERO NEW INPUTS ADDED ACROSS #144–#149. BAM is not accumulating loose knobs; it carries a fixed, categorized budget while the derived ledger grows.

THE ANCHOR SECTOR. Exactly one dimensionful input is mathematically mandatory (B4 scale-freeness, #52), relocatable along G → λ_crit → σ → R_MID with ΔR/R_MID = 0.52 the invariant unit (#53/#133); the √6 RS tuning is a derived constant (#57), not a knob.

THE UNIVERSAL RESIDUALS. α and √σ/m_e: their structure is derived (the charge quantum, the 1/2π measure, the running, the full one-loop EM sector #141–#147; the one-G repackaging #106) but their values match no principled closure number — the re-run scans put the best principled candidates at −4.2% (k₅³ + 2π for α⁻¹) and −5.4% (2π·k₅³ for 830), with every sub-% match requiring an ad-hoc term (#107/#108). Plausibly irreducible — and shared with every current theory.

THE PROGRAM RESIDUALS. n_part: the even doubling N_q = 2·n_part is topological (APS index, #123) and §8-stable; the value drifts 216–255 across the ablations, so 4n_part − 100 ∈ [764, 920] against the fixed 830 — recycling the compensator is circular (#107). ε: order-of-magnitude derived (~R_c³), window [2π, k₅√(2π)] (#89/#112). The lepton contrast: N = 4k₅² = 100 and 3 generations fully derived, no residual (#124).

THE BRACKETED SUB-RESIDUALS. k·r_s ∈ (0, 0.0064–0.070]: bounded above by the locked spectrum (sensitivity re-checked here, c ≈ 9.9) and below by static-throat existence (#148). The ε_n spread ∈ [1.32, 1.44] per step, data-pinned to ~0.3%, with the χ-driven law and ALL single power laws excluded (#149). Neither is free: both are boxed by derived structure.

THE NO-LOOSE-KNOBS CLAIM. #144 (Π/running), #145 (Z₁ = Z₂), #146 (G_E), #147 (F₁/F₂ capstone), #148 (k·r_s), #149 (ε_n): six probes, ZERO inputs added — every one derived structure, re-verified keystones, or tightened a bracket. The budget today is the #104/#125 budget.

SCOPE. The synthesis organizes, it does not remove: deriving any residual would change the ledger, and none has been achieved. One anchor, a short categorized list of bounded residuals, the universal flavor puzzle — the honest ledger.
