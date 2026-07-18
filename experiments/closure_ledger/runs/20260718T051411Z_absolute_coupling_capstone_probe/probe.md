# The absolute-coupling capstone (PR #225, FINAL)

**Run:** 2026-07-18T05:14:11+00:00

The deliverable is `docs/absolute_coupling_capstone.md` - canonical Hopf-KK normalization, geometric alpha, and the global no-retuning holdout. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the capstone question | **PASS** |
| T2 | `T2_canonical_chain` | the canonical chain; the #193 weld | **PASS** |
| T3 | `T3_charge_normalization` | charge = winding: exact flow + flux quantization | **PASS** |
| T4 | `T4_force_check` | the adiabatic-ramp force check | **PASS** |
| T5 | `T5_answer` | alpha = 4k^2/rho^2 exactly; the modulus REMAINS | **PASS** |
| T6 | `T6_stabilizer` | the EM cap selects rho* = 23.4; guardrail: no match | **PASS** |
| T7 | `T7_holdout` | the global no-retuning holdout passes | **PASS** |
| T8 | `T8_honest_scope` | honest scope | **PASS** |
| T9 | `T9_assessment` | assessment: the arc closes | **PASS** |

## Verdict

**Class:** `CANONICAL_HOPF_KK_NORMALIZATION_GIVES_ALPHA_EQUALS_4K2_OVER_RHO2_EXACTLY_A_CONTINUOUS_RADION_MODULUS_REMAINS_THE_EM_CAP_IS_ITS_DYNAMICAL_STABILIZER_AND_THE_GLOBAL_NO_RETUNING_HOLDOUT_PASSES`

ESTABLISHED (the argument is in docs/absolute_coupling_capstone.md).

THE CHAIN. The fiber KK tower (FD = closed form to 2e-11; continuum to 4e-04) welds to the #193 Berger spectra (identity exact) with the Wu-Yang half-charge; charge = fiber winding with spectral-flow slope 2k/R_f^2 to 3e-04 and EXACT flux quantization (4e-11); the adiabatic-ramp force check chirps at the predicted frequency to 3e-04.

THE ANSWER. alpha_k = 4 k^2 (l_P/R_f)^2 EXACTLY - and the kinematics is invariant under rescaling the fiber (machine zero: 0e+00): A CONTINUOUS RADION MODULUS REMAINS.  Uniqueness fails at the canonical level; alpha is not a kinematic pure number.

THE STABILIZER. The EM cap selects rho* = 2/sqrt(alpha) = 23.4125 (the Hopf fiber at ~23 Planck lengths); the #165 guardrail scan finds no closure-constant match (nearest e^pi at 1.2% - rejected): the selection is dynamical, correctly located outside canonical kinematics.

THE HOLDOUT. 8 keystone constants of the arc re-read from the committed ledgers and independently re-derived - the Bessel universal, the pi/2 step, the quarter wave, sqrt(mu_crit), the Rabi identity, the Weyl commutator, phi_h, beta_lepton - ALL ratios, roots, or topological phases, independent of rho: fixing rho* RETUNES NOTHING.  The arc closes with its books balanced.
