# Bounce path-integral derivation of the channel-dominant saddle (PR #152)

**Run:** 2026-06-10T04:56:55+00:00

Derives the pair-tunneling rule #151 modelled: the channel-conversion vertex lives at the cavity mouths (the neck is single-channel, #88/#132), so the off-diagonal bounce amplitude is a two-path sum dominated by the cheaper tunneling segment — channel dominance. Confirmed by exact diagonalization of a controlled multi-channel double well, with the counterfactual (vertex inside the barrier) flipping the rule to factorized. The #151 β knob is retired: modelled → derived, no new input. *(QFT on the classical throat, not quantum gravity.)*

- **Decomposition**: A_nm = c_near·e^{−S_m} + c_far·e^{−S_n} ≍ O(1)·e^{−min(S)}
- **Exact test**: t/max constant (×<1.3); t/geo varies ×9 — channel-dominant
- **Counterfactual**: vertex in barrier ⟹ rule flips to factorized (t/geo constant)
- **Linearity**: t ∝ W₀ exactly — single conversion vertex
- **Closure**: #151 consequences derived-footed; β knob retired; budget unchanged
- **Open**: full 5D bounce path integral; anarchic prefactor; PMNS angles/CP

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | derive the #151 channel-dominant saddle (retire the β knob) | **PASS** |
| T2 | `T2_two_path_decomposition` | two-path decomposition: A_nm ≍ O(1)·e^{−min(S_n,S_m)} | **PASS** |
| T3 | `T3_controlled_model` | controlled model: 3 channels ×2000 splitting span; extraction faithful | **PASS** |
| T4 | `T4_exact_element_is_channel_dominant` | t/max constant, t/geo ×9 — exact element channel-dominant | **PASS** |
| T5 | `T5_counterfactual_vertex_in_barrier` | vertex in barrier ⟹ rule flips to factorized (location decides) | **PASS** |
| T6 | `T6_linearity_and_151_closure` | t ∝ W₀ (single vertex); #151 consequences derived-footed | **PASS** |
| T7 | `T7_ledger_and_scope` | ledger: rule derived; β retired; budget unchanged | **PASS** |
| T8 | `T8_assessment` | CHANNEL_DOMINANT_SADDLE_DERIVED_MOUTH_CONVERSION_COUNTERFACTUAL_FLIPS | **PASS** |

## The controlled model

| channel | barrier height | exact splitting | WKB action |
|---:|---:|---:|---:|
| 1 | 1.0 | 1.458e-08 | 15.44 |
| 2 | 0.55 | 1.063e-06 | 11.43 |
| 3 | 0.3 | 2.851e-05 | 8.4 |

## The rule test (mouth-localized conversion)

| pair | exact \|t\| | t/(Δ_max/2) | t/(Δ_geo/2) |
|---|---:|---:|---:|
| (1,2) | 8.905e-08 | 0.1675 | 1.43 |
| (1,3) | 2.922e-06 | 0.205 | 9.06 |
| (2,3) | 2.883e-06 | 0.2023 | 1.05 |

max-rule spread ×1.22 (constant — the O(1) conversion factor); geo-rule spread ×8.65 (= e^{|ΔS|/2}). The exact element is channel-dominant.

## The counterfactual (vertex inside the barrier)

| pair | t/(Δ_max/2) | t/(Δ_geo/2) |
|---|---:|---:|
| (1,2) | 0.00519 | 0.044 |
| (1,3) | 0.0016 | 0.071 |
| (2,3) | 0.0109 | 0.056 |

The rule flips: geo-rule spread ×1.59 (constant), max-rule spread ×6.79. The conversion-vertex location decides — and the #151 data exclusion of factorized corroborates mouth conversion.

## Verdict

**CHANNEL_DOMINANT_SADDLE_DERIVED_MOUTH_CONVERSION_COUNTERFACTUAL_FLIPS.** THE #151 CHANNEL-DOMINANT SADDLE RULE IS DERIVED: THE CONVERSION VERTEX LIVES AT THE MOUTHS, THE TWO-PATH DECOMPOSITION GIVES A_nm ≍ O(1)·e^{−min(S_n,S_m)}, THE EXACT MULTI-CHANNEL COMPUTATION CONFIRMS IT, AND THE COUNTERFACTUAL FLIPS THE RULE WHEN THE VERTEX MOVES INTO THE BARRIER — THE β KNOB IS RETIRED, MODELLED → DERIVED, NO NEW INPUT. #151 flagged "derive G_ij from the bounce path integral" as its lead open item; this probe supplies it.

THE TWO-PATH DECOMPOSITION. The channel-conversion vertex (the anarchic O(1) overlap, #91/#151) has support only where the overtone wavefunctions coexist — the cavity mouths; the neck interior is the single-channel tunneling region (#88/#132). So the off-diagonal amplitude is a two-path sum, A_nm = c_near·e^{−S_m} + c_far·e^{−S_n} (convert-then-tunnel ⊕ tunnel-then-convert), dominated by the cheaper segment: A_nm ≍ O(1)·e^{−min(S_n,S_m)} — channel dominance.

THE EXACT COMPUTATION. A controlled three-channel double well (WKB actions ≈ 15.4/11.4/8.4, splittings spanning ×2000) with mouth-localized coupling, solved exactly and Löwdin-projected: |t|/(Δ_max/2) is CONSTANT across pairs (spread ×1.22) while |t|/(Δ_geo/2) varies ×8.65 (= e^{|ΔS|/2}, the two-path prediction). The exact element follows the channel-dominant rule with the constant equal to the O(1) conversion factor.

THE COUNTERFACTUAL DECIDES. Move the vertex INSIDE the barrier and the rule flips: t/geo becomes the constant (spread ×1.59) while t/max varies ×6.79 — conversion mid-tunnel gives the factorized geometric mean. The vertex LOCATION decides the saddle rule — and the measured r₃₂ already excluded factorized (#151, 0.1th percentile), so the data independently corroborate mouth conversion, which is the BAM cavity/neck structure itself.

ONE VERTEX. t ∝ W₀ exactly: a single conversion — the two-segment saddle, not a multi-conversion path.

THE CLOSURE. The #151 chain now stands on derived footing: mouth conversion (geometry) ⟹ channel dominance (this probe) ⟹ measured r₃₂ natural at the ~77th percentile, large mixing 0.44, and the falsifiable m₁ ≈ 0.04 meV / Σm_ν ≈ 58.8 meV prediction. The β interpolation knob is retired; the #150 budget is unchanged (one modelling assumption REMOVED).

SCOPE. The 1D model is the controlled analogue of the neck; the full 5D bounce path integral, the anarchic prefactor distribution (the localized flavor residual), explicit PMNS angles, and CP phases remain open.
