# The v3+CP joint re-lock: the v4 candidate lock (PR #163)

**Run:** 2026-06-12T22:42:02+00:00

Realizes the #161 target state in a structured parameter set. The minimal v3 law fails by an exact two-equality no-go — and the breaking pattern is the partition asymmetry, concentrated on the minus block's d-row, exactly where #155 located the physical mixing. Three new targeted couplings (continuing the lock's own χ/η pattern) plus one retune complete the realization: the v3 masses inherited exactly, all nine flavor-CP observables at ≤ 1% at the derived φ_h = π/k₅, net predictive surplus +2, CP at zero parameters. The library migration is staged. *(QFT on the classical throat, not quantum gravity.)*

- **The no-go**: partition split ×1.424 + minus d-row ×2.0 break the v3 law (up block survives 5e-6)
- **The v4 lock**: v3 law + η_12^+, η_12^−, η_13^− (new) + η_35^− retune (+11.7%) + diag retunes
- **Verified**: masses exact (1e-15); nine flavor-CP observables ≤ 1% at φ_h = π/k₅
- **Counting**: +3 params, +5 independent observables (net +2); CP at 0 params
- **Migration**: four-step library follow-up staged; library untouched here

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | realize the #161 targets — the v4 candidate lock | **PASS** |
| T2 | `T2_targets_reproduced` | targets reproduced (residual 0.005; eigenvalues 1e-15) | **PASS** |
| T3 | `T3_minimal_law_no_go` | minimal-law no-go: two equalities broken; up block survives | **PASS** |
| T4 | `T4_v4_element_lock_verified` | v4 verified: masses exact; nine observables ≤ 1% | **PASS** |
| T5 | `T5_structured_realization` | three new targeted couplings + one retune (partition pattern) | **PASS** |
| T6 | `T6_predictive_counting` | +3 params, +5 independent observables; CP at 0 params | **PASS** |
| T7 | `T7_staged_migration` | library migration staged (four steps) | **PASS** |
| T8 | `T8_assessment` | V4_LOCK_REALIZED_THREE_TARGETED_COUPLINGS_NO_GO_PARTITION | **PASS** |

## The v4 structured couplings

| coupling | value | element | factor vs v3 |
|---|---:|---|---:|
| η_12^+ (new) | -0.1018 | H₊[12] | ×1.287 |
| η_12^− (new) | -0.2953 | H₋[12] | ×1.832 |
| η_13^− (new) | -0.2671 | H₋[13] | ×1.996 |
| η_35^− (retune) | 5.5861 | H₋[23] | ×1.111 |

Diagonal retunes — up: [0.0019, -0.0019, -0.0], down: [0.1128, -0.0647, -0.048] (inside the existing diagonal law's reach).

## The v4 lock, verified

| observable | ratio to observed |
|---|---:|
| V_us | ×1.0025 |
| V_cb | ×1.0024 |
| V_ub | ×0.9976 |
| V_td | ×1.008 |
| V_ts | ×1.0021 |
| J | ×1.0044 |
| (β, γ, α) | (22.29, 65.9, 91.81)° vs (22.2, 65.9, 91.9)° |
| sin δ | 0.889 vs 0.887 |

## Verdict

**V4_LOCK_REALIZED_THREE_TARGETED_COUPLINGS_MINIMAL_LAW_NO_GO_PARTITION.** THE #161 TARGETS ARE REALIZED AS THE v4 CANDIDATE LOCK: THE MINIMAL LAW FAILS BY AN EXACT TWO-EQUALITY NO-GO WHOSE BREAKING PATTERN IS THE PARTITION ASYMMETRY; THREE NEW TARGETED COUPLINGS PLUS ONE RETUNE COMPLETE THE REALIZATION; THE v3 MASSES ARE INHERITED EXACTLY AND ALL NINE FLAVOR-CP OBSERVABLES LAND AT ≤ 1% — NET PREDICTIVE SURPLUS +2, CP AT ZERO PARAMETERS.

THE NO-GO. The v3 off-diagonal law enforces (i) partition-symmetric transport magnitudes (H₊[12] = H₋[12]) and (ii) the dk = max degeneracy (H[13] = H[23]); the targets break both — a partition split of ratio 1.424 and a minus d–b enhancement ×1.996 — while in the up block, where data permits, the law survives EXACTLY (5e-6). Every required deviation sits in the partition-asymmetric sector, on the minus block's d-row: exactly where #155 located the physical mixing, and the sector where the v3 lock already carries targeted couplings.

THE v4 LOCK. Element level: twelve tabulated numbers + the derived phases reproduce the v3 masses exactly and all nine observables at ≤ 1% — the first complete flavor state of the program in one parameter set. Structured level: the v3 law + η_12^+ = −0.102, η_12^− = −0.295, η_13^− = −0.267 (new) + η_35^−: 5.0 → 5.586 (retune) + diagonal retunes within the existing diagonal law — the extension continues the lock's own targeted-coupling pattern.

THE COUNTING. +3 parameters for +5 independent observables (the other four of the nine follow from unitarity and the derived phase): net surplus +2; the CP sector costs zero parameters (φ_h = π/k₅ derived); the #150 budget unchanged.

THE MIGRATION. Staged in four steps (QuarkParams fields, the complexified transport with the derived default, the lock update, the regression re-baseline) — one dedicated follow-up PR; the library is untouched here.
