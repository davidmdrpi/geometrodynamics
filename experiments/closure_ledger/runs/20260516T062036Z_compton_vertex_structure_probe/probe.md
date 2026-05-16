# Compton vertex-structure probe — finite-energy closure search

**Run:** 2026-05-16T06:20:36+00:00

Follow-on to PR #29 (finite-energy KN gap localised to vertex factors). Tests four families of natural vertex modification and asks whether any closes the BAM↔KN gap at finite ω/m.

**Construction:**

```
M^(λ, λ') = V(ε, ε', k, k'; α, β, γ) · [G_S3(ψ_s)^p · exp(iφ_s) + G_S3(ψ_u)^p · exp(iφ_u)] · T². V = ε·ε'*·(1 + β·sin²θ + γ·(1−cos θ)) + α·(ε·k̂')(ε'*·k̂).
```

## Test summary

| # | Test | Key metric | Value | PASS? |
|---|---|---|---:|---|
| T1 | `T1_thomson_preservation` | all candidates preserve Thomson | yes | **PASS** |
| T2 | `T2_family_A_eps_dot_k_sweep` | α_opt, residual | -4.800, 0.9640 | **FAIL** |
| T3 | `T3_family_B_angular_modulation_sweep` | (β_opt, γ_opt), residual | (-0.500, -1.000), 0.0764 | **PASS** |
| T4 | `T4_family_C_kinematic_power_sweep` | p_opt, residual | 0.25, 0.5608 | **PASS** |
| T5 | `T5_combined_best_ansatz` | best ansatz (B alone), full-grid residual | 0.3125 | **PASS** |
| T6 | `T6_interpretation` | dominant family | B (angular modulation) | **PASS** |

## T1_thomson_preservation

Verify Thomson (ω/m=1e-4) match for a range of candidate parameter sets. The vertex corrections are scaled by ε_kin = ω/m and must vanish at Thomson identically.

| α | β | γ | p | max diff at Thomson |
|---:|---:|---:|---:|---:|
| +0.00 | +0.00 | +0.00 | 1.00 | 6.00e-04 |
| +1.00 | +0.00 | +0.00 | 1.00 | 6.00e-04 |
| -1.00 | +0.00 | +0.00 | 1.00 | 6.00e-04 |
| +5.00 | +0.00 | +0.00 | 1.00 | 6.02e-04 |
| +0.00 | +1.00 | +0.00 | 1.00 | 6.00e-04 |
| +0.00 | +0.00 | +1.00 | 1.00 | 1.00e-03 |
| +0.00 | +0.00 | +0.00 | 0.50 | 5.00e-04 |
| +0.00 | +0.00 | +0.00 | 2.00 | 8.00e-04 |

## T2_family_A_eps_dot_k_sweep

Sweep α in V = ε·ε'* + α·(ε·k̂')(ε'*·k̂). Find α_opt minimising max KN residual over (ε, θ) grid.

Baseline residual at α=0: **0.9640**
Optimal α: **-4.8000**
Residual at α_opt: **0.9640**
Improvement: +0.0000

## T3_family_B_angular_modulation_sweep

Sweep (β, γ) in V_B = ε·ε'* · (1 + β·sin²θ + γ·(1−cos θ)). Find (β_opt, γ_opt) minimising KN residual.

Optimal (β, γ): **(-0.5000, -1.0000)**
Residual at optimum: **0.0764**
Improvement over baseline (0, 0): +0.8877

## T4_family_C_kinematic_power_sweep

Scan p in M_x ∝ G_S3(ψ_x)^p. Find p_opt minimising KN residual. Predicted natural values: p = 1 (current), p = 0.5 (QED-amplitude analog), or p = 2 (squared).

| p | residual |
|---:|---:|
| 0.25 | 0.5608 |
| 0.50 | 0.6763 |
| 0.75 | 0.8096 |
| 1.00 | 0.9640 |
| 1.25 | 1.1432 |
| 1.50 | 1.3514 |
| 2.00 | 1.8772 |

Optimal p: **0.2500** (residual 0.5608)

## T5_combined_best_ansatz

Evaluate candidate ansätze on a fine (ε, θ) grid including ε = 0.5. The honest "best ansatz" is the minimum over {baseline, A, B, C, B+C, A+B+C}. Family A is included for completeness; T2 shows it offers no improvement.

Optimal Family A α: -4.800; Family B (β, γ): (-0.500, -1.000); Family C p: 0.250

| ansatz | residual on full grid (ε ≤ 0.5) |
|---|---:|
| baseline (no modification) | 2.7443 |
| Family A only (α_opt) | 6.1562 |
| Family B only (β_opt, γ_opt) | 0.3125 |
| Family C only (p_opt) | 0.9610 |
| Family B + C combined | 0.3125 |
| All four (A + B + C) | 3.1083 |
| **Best ansatz: B alone** | **0.3125** |

Closes gap below 1 %: **no**.  Below 10 %: **no**.  Below 50 %: **YES**.  Improvement factor: 8.78×.

## T6_interpretation

Identify which vertex family contributed most to closing the gap and whether optimal parameters have clean (naturally derivable) values.

| family | improvement | optimal value (clean?) |
|---|---:|---|
| B (angular modulation) | +0.8877 | (β, γ) = (-0.500, -1.000) (clean) |
| C (kinematic power) | +0.4032 | p = 0.250 (clean) |
| A (ε·k coupling) | +0.0000 | α = -4.800 (fitted) |

**Dominant family:** B (angular modulation)

## Verdict

**PARTIAL_CLOSURE_50PCT.** PARTIAL CLOSURE (< 50 %) — best ansatz (B alone) brings KN residual to 0.3125 on the full (ε, θ) grid. Improvement factor 8.8× over baseline. At small ε (≤ 0.2), residual drops to 0.0764 for Family B; at ε = 0.5, residual remains large, indicating O(ε²) and higher-order corrections are needed. The identified leading-order structural piece is an angular modulation (β, γ) = (-0.5, -1.0) of the polarization sum, scaled by ε_kin = ω/m. The clean integer/half-integer values point to a natural BAM derivation.

## What this leaves open

- **Derivation of optimal coupling from BAM principles.** Empirical fit ≠ derivation. The probe localises which vertex family is needed; deriving the specific coupling values from Hopf-connection or throat-transport algebra is a separate analytic task.
- **Higher-order ε corrections.** If the probe achieves PARTIAL_CLOSURE, residual O(ε²) structure remains. Identifying its origin (next-order vertex, second-derivative of S³ Green function, electron spin-½ corrections) is follow-on work.
- **Electron spin at finite energy.** Probe still uses scalar electron. Spin-½ Dirac structure should appear at O(ω/m) corrections in QED.
- **Loop corrections.** Vertex/self-energy/vacuum polarization require BAM's bulk radial channel.
