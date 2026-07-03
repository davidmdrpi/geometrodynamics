# The coupled 5D+soliton solve - companion probe (PR #203)

**Run:** 2026-07-03T03:48:30+00:00

The deliverable is `docs/coupled_5d_soliton_solve.md` - the confrontation of the exact #202 law with the locked #180 energetics. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the final register item; the edge is live | **PASS** |
| T2 | `T2_no_knobs_left` | no knobs left: exact law + locked solution | **PASS** |
| T3 | `T3_scales_extracted` | the scales from one solution (RMS, R*, q-core, Phi0) | **PASS** |
| T4 | `T4_confrontation` | the edge FIRED: band 4.6-6.7 vs 88.6 (m_e x15 over) | **PASS** |
| T5 | `T5_bound_and_nonclosability` | upper bound holds; gap not closable in weak field | **PASS** |
| T6 | `T6_nr_target` | the NR target: core contraction 13-45x (falsifiable) | **PASS** |
| T7 | `T7_honest_scope` | scope: locked, banded, no new fits; weak-field strained | **PASS** |
| T8 | `T8_assessment` | the thread closes onto one number | **PASS** |

## The confrontation

| sigma_mode | r_core | ratio |
|---|---|---:|
| RMS | r_q_half | 1.528 |
| RMS | rho_c | 1.175 |
| R*(A) | r_q_half | 6.018 |
| R*(A) | rho_c | 4.629 |
| R*(B) | r_q_half | 6.711 |
| R*(B) | rho_c | 5.163 |

(needed: 88.6 conv A / 206.8 conv B; m_e over-prediction x13.2)

## The binding sweep (the trend is the wrong way)

| M | RMS | r_q(half) | RMS/r_q | Phi(0) |
|---:|---:|---:|---:|---:|
| 2.75 | 1.5837 | 0.5 | 3.1674 | -2.5753 |
| 3.5 | 1.2733 | 0.8333 | 1.5279 | -4.2158 |
| 4.5 | 1.0379 | 0.8333 | 1.2455 | -6.8646 |

## Verdict

**WEAK_FIELD_COUPLED_SOLVE_OVERPREDICTS_M_E_BY_15X_GAP_NOT_CLOSABLE_IN_WEAK_FIELD_THE_NR_CORE_CONTRACTION_13_TO_45X_IS_THE_TARGET.** THE EDGE FIRED, AND THE RESULT IS CLEAN (the argument is in docs/coupled_5d_soliton_solve.md; this probe computes everything live from the locked solution).

THE CONFRONTATION. With no knobs left - the #202 law is exact, the #180 solution locked - the coupled weak-field solve gives sigma_mode/r_core = [4.63, 6.71] (pairing definitions) versus the required 88.6: m_e over-predicted by ~x13. The weak-field solve does not land the electron mass ratio.

THE STRUCTURE OF THE FAILURE. The direction is right - the true 5D core is the strong-field endpoint of the #179 runaway, smaller than the weak-field q-core, so the weak-field value is an UPPER BOUND on m_e, which holds. And the gap is provably not closable inside the weak-field model: the binding sweep gives RMS/r_q = [3.1674, 1.5279, 1.2455] - the ratio moves the WRONG way with binding strength. The missing factor is physics absent from the model.

THE TARGET. One number remains: the NR core contraction r_q/r_s = [13.2, 44.7] - with the failure mode stated (an O(1) contraction refutes the pairing mechanism as m_e's quantitative origin). The mass-ladder thread - from the #192 fine-tuning to here - closes onto a single dimensionless output of a well-posed GR computation, with every step upstream exact, measured, or bounded, and no knobs added anywhere.
