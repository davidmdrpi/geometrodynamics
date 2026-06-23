# The sensitivity audit: Jacobian rank, the forced core, the isolation dimension (PR #173)

**Run:** 2026-06-23T07:23:40+00:00

The dynamical inverse problem: vary the continuous geometry and MEASURE the observable response, to quantify the predictive content rather than assert it. *(QFT on the classical throat, not quantum gravity.)*

- **Isolation dimension**: rank(J) = 10 of 14 observables
- **Forced core**: 4 (entirely CKM unitarity relations; V = U₊†U₋ structural)
- **Masses**: fitted (no forced mass relation, quark or lepton)
- **Redundancy**: 10 compensator directions, dominated by the diagonal shifts
- **CP test**: deriving φ_h saves 0 effective inputs (CP direction already spanned)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal_inverse_problem` | the inverse problem (rank, forced core, isolation) | **PASS** |
| T2 | `T2_observables_and_inputs` | live observables vs free inputs (k₅-locks excluded) | **PASS** |
| T3 | `T3_jacobian_rank` | Jacobian rank: isolation dimension 10/14 (clean gap) | **PASS** |
| T4 | `T4_forced_core_ckm_unitarity` | forced core = 4 = CKM unitarity (largest zero-cost set) | **PASS** |
| T5 | `T5_masses_are_fitted` | the masses are fitted (no forced mass relation) | **PASS** |
| T6 | `T6_compensator_redundancy` | compensator redundancy 10 (the diagonal-shift over-completeness) | **PASS** |
| T7 | `T7_cp_at_zero_cost_test` | CP-at-zero-cost: φ_h adds no rank (honest) | **PASS** |
| T8 | `T8_assessment` | JACOBIAN_AUDIT_FORCED_CORE_4_CKM_UNITARITY | **PASS** |

## The quark singular-value spectrum (the rank gap)

| i | singular value |
|---|---:|
| 0 | 2.263e+01 |
| 1 | 1.185e+01 |
| 2 | 2.819e+00 |
| 3 | 1.474e+00 |
| 4 | 4.401e-01 |
| 5 | 3.989e-01 |
| 6 | 3.723e-02 |
| 7 | 1.244e-02 |
| 8 | 6.492e-08 ← rank gap |
| 9 | 1.969e-09 |
| 10 | 4.723e-11 |
| 11 | 1.126e-11 |

## Verdict

**JACOBIAN_AUDIT_FORCED_CORE_4_CKM_UNITARITY_MASSES_FITTED_REDUNDANCY_10.** MEASURED, NOT ASSERTED. The Jacobian of the live observables with respect to the free continuous geometry gives a definite, and honestly mixed, picture of the predictive content.

ISOLATION DIMENSION. rank(J) = 10 of 14 observables (quark 8, lepton 2), with a clean singular-value gap — the free knobs dial that many independent observable directions.

THE FORCED CORE. n_obs − rank = 4, entirely CKM combinations: the CKM UNITARITY relations. BAM's V = U₊†U₋ is exactly unitary (‖V†V−I‖ = 7e-16), so the 8 CKM observables lie on the 4-parameter unitary manifold and 4 relations are forced at zero input cost — the largest such set. A genuine structural prediction, but the standard unitarity, not a BAM-specific numerical relation.

THE MASSES ARE FITTED. The quark and lepton mass ratios carry no weight in the forced core — the ladder sets their values, but the knobs span them; there is no forced mass relation in the live observable set.

THE REDUNDANCY. n_inputs − rank = 10 compensator directions, dominated by the mass-preserving diagonal shifts (kernel share 0.68) — the over-completeness the program flagged qualitatively, now measured. The v4 quark parametrization is substantially over-complete.

CP AT ZERO COST. Adding φ_h as an input leaves the rank unchanged (8 → 8): the CP-phase direction is already spanned, so deriving φ_h saves no effective input — 'CP at zero parameters' is a counting economy, not a Jacobian reduction. An audit: a measurement, honest where it is not flattering.
