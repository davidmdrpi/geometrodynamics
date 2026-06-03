# The factorized BAM sector sum Z (PR #122)

**Run:** 2026-06-03T05:40:33+00:00

Assembles the validated ingredients (PRs #74, #115–#121) into the full BAM loop-measure sector sum, and shows it **factorizes** into a discrete Z₂-signed (topological) sum × a continuous η-phased (analytic) moduli integral — with the Z₂ grading **cancelling the leading UV**.

```
Z = Σ_{k odd, c₁∈ℤ, n_part}  (−1)^k  ∫₀^∞ (dL/L)  det^{−1/2}_matter(L) · e^{i(π/2)(1−2a)} · e^{−S_BAM}
    └──── discrete Z₂-signed (topological) ────┘ └──── continuous η-phased (analytic) ────┘
```

- **Factorization**: discrete Z₂-signed (topological) sum ⊗ continuous η-phased (analytic) integral; (−1)^k pulls out
- **No double-count**: U(1) η-phase ⟂ Z₂ sign (PR #121)
- **Graded UV**: θ_per − θ_anti ~ e^{−π²/t} → 0: Z₂ grading cancels the BC-independent Weyl UV term
- **Open**: absolute normalization (κ₅²/Λ₅); non-perturbative convergence; multi-loop

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | assemble factorized Z from PRs #74, #115–#121 | **PASS** |
| T2 | `T2_assembled_measure_formula` | the assembled measure formula (each factor + source PR) | **PASS** |
| T3 | `T3_discrete_sector_sum_signs_pull_out` | discrete sector sum; (−1)^k sector-constant pulls out | **PASS** |
| T4 | `T4_continuous_moduli_integral` | continuous moduli integral: dL/L + finite dets + η-phase | **PASS** |
| T5 | `T5_factorization` | Z = (discrete Z₂ sum) ⊗ (continuous integral); no double-count | **PASS** |
| T6 | `T6_z2_graded_uv_cancellation` | Z₂-graded UV cancellation: θ_per − θ_anti ~ e^{−π²/t} → 0 | **PASS** |
| T7 | `T7_scope` | scope: assembled; absolute normalization / non-pert open | **PASS** |
| T8 | `T8_assessment` | BAM_FACTORIZED_SECTOR_SUM_Z_..._GRADED_UV_CANCELS | **PASS** |

## The assembled measure (each factor validated)

| factor | meaning | source |
|---|---|---|
| `Σ_{k,c₁,n_part}` | closure-ledger sector sum | #115 |
| `(−1)^k` | discrete Z₂ orientation sign | #115/#118/#121 |
| `∫ (dL/L)` | gauge-fixed moduli measure (CKV = closure quantum) | #74/#117/#118 |
| `det^{−1/2}_matter(L)` | matter fluctuation det (finite, GY) | #116 |
| `det'(P) = L` | first-order ghost det | #117/#118 |
| `e^{i(π/2)(1−2a)}` | continuous η-phase (holonomy a) | #119/#121 |
| `e^{−S_BAM}` | leading bounce saddle | #87–#90 |

## The Z₂-graded UV cancellation

| t | θ_per | θ_anti | Weyl L/√(4πt) | θ_per − θ_anti |
|---:|---:|---:|---:|---:|
| 0.2 | 3.96333 | 3.96333 | 3.96333 | 0.0 |
| 0.1 | 5.60499 | 5.60499 | 5.60499 | 8.882e-16 |
| 0.05 | 7.92665 | 7.92665 | 7.92665 | 0.0 |
| 0.02 | 12.53314 | 12.53314 | 12.53314 | 1.776e-15 |

Each `θ ~ L/√(4πt) → ∞` as `t → 0` (UV divergent), but the Z₂-graded difference `θ_per − θ_anti ~ e^{−π²/t} → 0` is **UV-finite** — the BC-independent Weyl term cancels between the orientable (+) and Möbius (−) sectors.

## Verdict

**BAM_FACTORIZED_SECTOR_SUM_Z_DISCRETE_Z2_TIMES_CONTINUOUS_ETA_GRADED_UV_CANCELS.** THE BAM LOOP MEASURE ASSEMBLES INTO A FACTORIZED SECTOR SUM: A DISCRETE Z₂-SIGNED (TOPOLOGICAL) SUM TIMES A CONTINUOUS η-PHASED (ANALYTIC) MODULI INTEGRAL, WITH THE Z₂ GRADING CANCELLING THE LEADING UV. PRs #74 and #115–#121 validated every ingredient; this probe puts them together.

THE ASSEMBLED MEASURE. Z = Σ_{k odd, c₁∈ℤ, n_part} (−1)^k ∫_0^∞ (dL/L) det^{−1/2}_matter(L) · e^{i(π/2)(1−2a)} · e^{−S_BAM}, with the closure-ledger sector sum (#115), the discrete Z₂ orientation sign (−1)^k (#115/#118/#121), the gauge-fixed dL/L moduli measure whose 1/L is the closure-quantum CKV factor (#74/#117/#118), the finite matter determinant (#116), the first-order ghost determinant det'(P) = L (#117/#118), the continuous η-phase e^{i(π/2)(1−2a)} (#119/#121), and the leading bounce saddle e^{−S_BAM} (#87–#90).

WHY IT FACTORIZES. The discrete Z₂ orientation sign (−1)^k is a SECTOR-CONSTANT — it depends only on the winding parity, not on the continuous moduli L or holonomy a — so it pulls OUT of the continuous integral: Z = Σ_{discrete sectors} (−1)^k × [∫(dL/L) det η-phase e^{−S}]. Thus Z is a discrete Z₂-signed sum of continuous moduli integrals; the continuous part carries the η-phase (confined to the right half-circle, #121) and the finite determinants, the discrete part the orientation signs. They do not double-count (#121: the U(1)-valued η-phase and the Z₂ sign are independent).

THE Z₂-GRADED UV CANCELLATION. Grouped by orientation, Z is the Z₂-graded combination of the orientable (periodic, +) and Möbius (antiperiodic, −) contributions. The leading heat-kernel (Weyl) coefficient a_{−1/2} = L/√(4π) is a BULK quantity, independent of the boundary condition, so it is identical in both sectors and CANCELS in the graded difference: the heat traces θ_per(t) and θ_anti(t) each diverge as L/√(4πt) as t → 0, but their difference θ_per − θ_anti ~ e^{−π²/t} → 0 is UV-finite (the polynomial UV divergence cancels to all orders, leaving only exponentially small instanton terms). So the orientation Z₂ grading renders the bulk UV of the sector sum finite.

SCOPE. ASSEMBLED: the full factorized one-loop sector sum Z, every factor finite/validated (PRs #74, #116–#121), with the Z₂ grading cancelling the leading UV. OPEN: the absolute normalization (the bulk κ₅²/Λ₅ anchor, PR #112), the full non-perturbative convergence of the sector sum, and the multi-loop measure. The assembly organizes the structure; it does not fix the overall scale.
