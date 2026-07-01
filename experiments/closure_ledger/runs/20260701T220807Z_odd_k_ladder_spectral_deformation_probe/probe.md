# The spectral deformation test: the {1,3,5} ladder on the Berger-squashed S³ (PR #192)

**Run:** 2026-07-01T22:08:07+00:00

Upgrades #183 from algebra to spectrum: the locked lepton Hamiltonian's geometric ingredients are rebuilt on the Berger sphere S³_λ (the #165 machinery) and the ladder is tracked as λ moves off 1. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | goal: algebra → spectrum; the pass/fail criteria | **PASS** |
| T2 | `T2_machinery` | machinery: #165 SU(2) Berger spectrum + locked λ=1 baseline | **PASS** |
| T3 | `T3_ingredient_map` | the ingredient map (declared before the sweep) | **PASS** |
| T4 | `T4_spectral_deformation` | the ladder deforms smoothly over a finite window | **PASS** |
| T5 | `T5_no_infinitesimal_breakdown` | NO infinitesimal breakdown (linear response at λ=1) | **PASS** |
| T6 | `T6_finite_boundary_and_fine_tuning` | λ_break = 0.986; μ/e sensitivity = 1/(1−λ_break) | **PASS** |
| T7 | `T7_robustness` | robust across all four ingredient maps | **PASS** |
| T8 | `T8_assessment` | structure protected; e–μ hierarchy metric-fine-tuned | **PASS** |

## The deformed ladder (action units; ratios electron-calibrated)

| λ | e level | μ level | τ level | μ/e | τ/μ |
|---:|---:|---:|---:|---:|---:|
| 0.99 | 0.0574 | 41.2008 | 688.2703 | 718.159 | 16.705 |
| 0.995 | 0.1286 | 41.2303 | 691.6267 | 320.601 | 16.775 |
| 1.0 | 0.1996 | 41.26 | 694.9832 | 206.679 | 16.844 |
| 1.01 | 0.3411 | 41.3198 | 701.6964 | 121.14 | 16.982 |
| 1.05 | 0.8991 | 41.5643 | 728.5516 | 46.229 | 17.528 |
| 1.1 | 1.5798 | 41.8818 | 762.1255 | 26.511 | 18.197 |
| 1.25 | 3.5204 | 42.9108 | 862.8725 | 12.189 | 20.109 |
| 1.5 | 6.4688 | 44.8578 | 1030.8377 | 6.935 | 22.98 |
| 2.0 | 11.5701 | 49.4568 | 1366.8588 | 4.275 | 27.637 |
| 3.0 | 19.9452 | 60.3843 | 2038.9991 | 3.028 | 33.767 |

## The boundary and the fine-tuning

- λ_break (electron level → 0): **0.985983** (a 1.4% squash — finite, not infinitesimal)
- electron level at round point: **0.19963** (vs μ 41.26, τ 694.98) — a near-zero
- d ln(μ/e)/dλ = **-70.92** vs −1/(1−λ_break) = −71.34 — the sensitivity IS the distance to the boundary
- d ln(τ/μ)/dλ = **0.822** — gentle; τ/μ is metric-robust

## Robustness (all ingredient maps)

| map | λ_break | d ln(μ/e)/dλ | 1/(1−λ_break) | identity | d ln(τ/μ)/dλ |
|---|---:|---:|---:|---|---:|
| default_all_fiber | 0.98598 | -70.92 | 71.34 | True | 0.822 |
| flip_resistance_blind | 0.98435 | -63.61 | 63.89 | True | 0.895 |
| flip_uplift_blind | 0.98573 | -69.68 | 70.06 | True | -0.081 |
| minimal_only_2pi_base | 0.96823 | -31.32 | 31.47 | True | -0.143 |

## Verdict

**ODD_K_LADDER_SPECTRUM_DEFORMS_SMOOTHLY_OVER_A_FINITE_BERGER_WINDOW_BUT_THE_MU_E_HIERARCHY_IS_METRIC_FINE_TUNED.** SPECTRUM, NOT JUST ALGEBRA — with a sharp caveat the algebra could not see.

THE PASS. The locked lepton Hamiltonian, rebuilt on the Berger-squashed S³_λ (fiber-riding ingredients × λ, connection/base ingredients fixed), keeps its {1,3,5} structure — three positive, ordered levels — over a FINITE window λ ∈ (0.9860, ≥3], and the mass ratios deform SMOOTHLY with a linear response at the round point. The ladder does not break at infinitesimal squash: the round metric is not smuggling in the sector structure, and the #183 protection claim is upgraded from algebraic invariants to the actual spectrum.

THE DISCOVERY. The protection does not extend to the HIERARCHY. The electron level at the round point is a near-zero (0.1996 in action units vs μ 41.3, τ 695) that crosses zero at a 1.4% fiber squash (λ_break = 0.98598); the μ/e log-sensitivity -70.9 equals −1/(1−λ_break) = −71.3 — the steepness IS the proximity to the spectral boundary — while τ/μ moves gently (+0.82). Robust across all four ingredient maps (flipped and minimal assignments). So the topology guarantees THREE GENERATIONS; the round metric tunes the e–μ hierarchy — the two claims are now separated by measurement, which is what the thesis needed.
