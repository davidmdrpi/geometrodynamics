# Three throat modes from `k_5`: `#gen = (k_5+1)/2 = 3`

**Run:** 2026-05-26T07:07:29+00:00

Closes the PR #70 follow-on "why exactly 3 throat modes" by deriving the count from the same `k_5 = 5` topological-charge primitive that gave `β_lepton = k_5²·(2π)` in PR #71. Lepton depths `{1, 3, 5}` are the odd integers in `[0, k_5]`; `#generations = (k_5+1)/2 = 3`. The cavity geometry independently shows ~3 throat-localized modes (PR #68's saturating crossover) — consistent.

- **Structural form**: `#generations = (k_5 + 1) / 2 = 3`
- **Lepton depths**: `[1, 3, 5]` (= e, μ, τ)
- **Shared with PR #71**: k_5 = 5 (β_lepton = k_5²·2π = 50π)
- **Cavity consistency**: PR #68 saturating crossover at ~3 throat modes
- **B4 caveat**: k_5 dimensionless integer; count structural/topological

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_odd_k_fermionic_selection` | lepton depths (1, 3, 5) all odd (recap #67) | **PASS** |
| T2 | `T2_topological_charge_cap` | γ-lock cap k ≤ k_5 = 5; odd subset = [1, 3, 5] | **PASS** |
| T3 | `T3_three_generations_count` | (k_5+1)/2 = 3 = 3 generations | **PASS** |
| T4 | `T4_k5_dependence` | k_5 = 3/5/7/9 → 2/3/4/5 generations | **PASS** |
| T5 | `T5_cavity_geometry_consistency` | cavity: 3 throat modes (#68 crossover) | **PASS** |
| T6 | `T6_unification_with_PR71` | same k_5 as #71: β=k_5²·2π=50π; #gen=(k_5+1)/2=3 | **PASS** |
| T7 | `T7_falsification_b4` | 4th lepton would require k_5≥7; observation 3 | **PASS** |
| T8 | `T8_assessment` | 3 generations from k_5 = 5 (Tangherlini) | **PASS** |

## T4: `k_5`-dependence of the generation count

| `k_5` | odd integers in `[0, k_5]` | `(k_5+1)/2` |
|---:|---|---:|
| 3 | `[1, 3]` | 2 |
| 5 | `[1, 3, 5]` | 3 ← **BAM** |
| 7 | `[1, 3, 5, 7]` | 4 |
| 9 | `[1, 3, 5, 7, 9]` | 5 |

## T5: Cavity-geometry consistency (#68)

| n | ⟨r⟩−R_MID | sector |
|---:|---:|---|
| 0 | 0.0212 | lepton (throat) |
| 1 | 0.0389 | lepton (throat) |
| 2 | 0.0440 | lepton (throat) |
| 3 | 0.0454 | shell (QCD) |
| 4 | 0.0460 | shell (QCD) |
| 5 | 0.0462 | shell (QCD) |

throat-localized lepton modes: 3; matches algebraic (k_5+1)/2 = 3: True

## Verdict

**THREE_GENERATIONS_FROM_K5.** #GENERATIONS = (k_5+1)/2 = 3. The three observed charged-lepton generations follow from the same Tangherlini topological-charge primitive (k_5 = 5) that derived β_lepton = k_5²·(2π) = 50π in PR #71.

STRUCTURAL FORM. Combining the odd-k fermionic selection (#67: charged leptons are spin-½ fermions ⟹ orientation-reversing closure ⟹ odd k) with the γ-lock l-cap (hbar_origin_status: the cross-species R_OUTER fixed point uses Σ V_max[0..k_5] = 22.5 with k_5 = 5), the lepton sector is the odd subset of [0, k_5]:

    lepton depths = {1, 3, 5} = (e, μ, τ),
    #generations = (k_5 + 1) / 2 = 3 .

k_5-DEPENDENCE. The count is locked to k_5: 2 (k_5=3), 3 (k_5=5, BAM), 4 (k_5=7), 5 (k_5=9). Observation of 3 charged leptons fixes k_5 = 5.

INDEPENDENT CAVITY CHECK (#68). The radial overtone ladder (l=1) shows the throat-localized modes n=0,1,2 (e,μ,τ) form the monotonic-rise lepton sequence, with n≥3 saturating into shell standing waves (the QCD/quark channel, #69). The cavity geometry independently supports ~3 throat-localized modes — consistent with the algebraic (k_5+1)/2 = 3. Honest note (from #68): the cavity count is a saturating crossover, not razor-sharp; the algebraic form is exact.

UNIFICATION. Both lepton-sector structural results come from the same k_5 primitive:
    β_lepton    = k_5²·(2π) = 50π   (PR #71, quadratic)
    #generations = (k_5+1)/2  = 3    (this PR,  linear)
The mass-scaling coupling and the generation count are the two faces of the Tangherlini topological charge. HONEST SCOPE. Closes the PR #70 follow-on by deriving the count from k_5; does NOT first-principles derive k_5 = 5 itself (the γ-lock R_OUTER fixed point in hbar_origin_status is the established input). B4: k_5 is a dimensionless integer; the count is structural/topological — scale-independent.

## What this leaves open

- **First-principles `k_5 = 5`.** Tied to the γ-lock cross-species fixed point at `R_OUTER ≈ 1.262`, `Σ V_max[0..5] = 22.5` (per `hbar_origin_status`); a deeper geometric derivation of why the `l`-range caps at 5 is future work.
