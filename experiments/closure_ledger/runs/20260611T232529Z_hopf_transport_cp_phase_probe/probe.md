# Hopf-connection derivation of the quark CP phase (PR #158)

**Run:** 2026-06-11T23:25:29+00:00

Runs the #156 acceptance test from the Hopf connection — and relocates the mechanism, correcting #156: the partition-mixing CP was a non-unitarity artifact (unitarized J ≈ 0 for every phase form); the true home is the Hopf-fiber transport phase of the same-partition shell couplings, with orientation sign set by the Z₂ partition class. One parameter calibrated to J alone predicts the full unitarity triangle to ~1°, and the pure closure value π/k₅ reproduces all five CP observables uncalibrated (candidate, flagged per #107/#108). *(QFT on the classical throat, not quantum gravity.)*

- **Correction**: #156 partition-mixing J was non-unitarity artifact (unitarized J ≈ 0)
- **Relocation**: cos(phase·dk) = Re(e^{iφ_h·dk}); (H±) ∝ e^{±iφ_h·dk} (#63 orientation)
- **Calibrated**: φ_h* = 0.611 → (β,γ,α) ≈ (22.9, 65.8, 91.3)° vs (22.2, 65.9, 91.9)°
- **Candidate**: π/k₅ = 0.6283 (2.7% away): 5 observables, 0 parameters — flagged
- **Open**: derive the orientation-signed transport explicitly; confirm π/k₅

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | run the #156 acceptance test from the Hopf connection | **PASS** |
| T2 | `T2_partition_mixing_excluded` | EXCLUSION (corrects #156): unitarized J ≈ 0; ε* ×40 over unitarity | **PASS** |
| T3 | `T3_hopf_transport_relocation` | relocation: e^{±iφ_h·dk} per partition class; exact unitarity | **PASS** |
| T4 | `T4_one_parameter_full_triangle` | one parameter ⟹ full triangle to ~1°; β = 22° test PASSED | **PASS** |
| T5 | `T5_pure_pi_over_k5_candidate` | pure π/k₅: five observables, zero parameters (candidate, flagged) | **PASS** |
| T6 | `T6_budget_impact` | budget: #156 input conditionally returned; new ε bound | **PASS** |
| T7 | `T7_scope` | scope: orientation rule motivated, not yet derived | **PASS** |
| T8 | `T8_assessment` | QUARK_CP_FROM_HOPF_TRANSPORT_TRIANGLE_PREDICTED_PI_OVER_K5_CANDIDATE | **PASS** |

## The exclusion (corrects #156)

| φ_q(k) form | J raw | J_db raw | non-unitarity | J unitarized |
|---|---:|---:|---:|---:|
| linear k | 7.66e-06 | -3.17e-09 | 0.158 | -3e-09 |
| winding k² | 7.66e-06 | -3.19e-09 | 0.158 | -3.02e-09 |
| Casimir k(k+2) | -2.14e-06 | 8.8e-10 | 0.159 | 8.3e-10 |

First-row unitarity deficit at the #156 point: 0.1584 vs the experimental scale ~0.005 — excluded ×40.

## The triangle: calibrated and pure-candidate points

| quantity | calibrated (φ_h = J-fit) | pure π/k₅ (no fit) | observed |
|---|---:|---:|---:|
| φ_h | 0.61102 | 0.62832 | — |
| J/target | 1.0 | 0.969 | 1.0 |
| β (°) | 22.89 | 22.78 | 22.2 |
| γ (°) | 65.79 | 63.48 | 65.9 |
| α (°) | 91.32 | 93.75 | 91.9 |
| sin δ | 0.905 | 0.888 | 0.887 |
| max mass shift | 0.000868 | 0.000902 | < 0.016 |
| V_cb | 0.0377 | 0.0377 | 0.0418 (stiff pred 0.0377) |

## Verdict

**QUARK_CP_FROM_HOPF_TRANSPORT_TRIANGLE_PREDICTED_PI_OVER_K5_CANDIDATE.** THE QUARK CP PHASE LIVES IN THE HOPF-FIBER TRANSPORT OF THE SAME-PARTITION SHELL COUPLINGS — ONE PARAMETER CALIBRATED TO J ALONE PREDICTS THE FULL UNITARITY TRIANGLE TO ~1°, AND THE PURE CLOSURE VALUE π/k₅ REPRODUCES ALL FIVE CP OBSERVABLES UNCALIBRATED (CANDIDATE, FLAGGED) — WHILE THE #156 PARTITION-MIXING ROUTE IS EXCLUDED AND CORRECTED.

THE CORRECTION TO #156. The partition-mixing CP was an artifact: 16% charged-current non-unitarity (the u–d near-degeneracy), quartet-inconsistent invariants (J₁₂ ~ 7.7e-6 vs J_db ~ −3e-9), a unitarized core with J ≈ 0 for EVERY φ_q(k) form, and an ε* that violates first-row CKM unitarity ~×40. Partition mixing is doubly excluded as the CP origin. (The #156 ceiling identity and J = 0 baseline stand — they are |V| arithmetic.)

THE RELOCATION. The locked coupling −t·e^{−α·dk}·cos(phase·dk) is the real part of the Hopf transport factor e^{iφ·dk}; the two Z₂ partition classes traverse the fiber with OPPOSITE orientation (the #63 C-swap flips c₁) ⟹ (H±) ∝ e^{±iφ_h·dk}. Hermitian blocks ⟹ exactly unitary V with quartet-consistent J — genuine CP; the locked baseline is φ_h = 0.

ONE PARAMETER, THE FULL TRIANGLE. φ_h* = 0.61102 from J alone ⟹ (β, γ, α) = (22.89, 65.79, 91.32)° vs (22.2, 65.9, 91.9)° — max deviation 0.69°; sin δ = 0.905; masses shifted 9e-04; V_cb untouched; V_us moves toward the data. The #156 acceptance test is PASSED.

THE π/k₅ CANDIDATE. The calibration sits 2.7% from π/k₅ = 0.6283 — the χ = 0 fiber holonomy over the k₅ winding quanta. Pure and uncalibrated: J at 0.969 of target, (β, γ, α) = (22.78, 63.48, 93.75)°, sin δ = 0.888 vs 0.887 — five observables, zero parameters. Flagged per #107/#108 as a candidate closure identification pending an independent transport derivation (the #152 path).

BUDGET. No input consumed; #156's input conditionally returned (downgraded today to a one-parameter family with a principled candidate); a new independent bound emerges (partition mixing ε ≲ 0.004 from first-row unitarity). SCOPE: the orientation-sign rule is motivated (C-swap), not yet derived from explicit transport — the flagged next step.
