# The absolute-coupling capstone (PR #225, FINAL; revised)

**Run:** 2026-07-21T02:20:10+00:00

The deliverable is `docs/absolute_coupling_capstone.md` - canonical Hopf-KK normalization, geometric alpha, the Einstein-frame radion, and the alpha-dependent holdout. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the capstone question + the revision | **PASS** |
| T2 | `T2_canonical_chain` | the canonical chain; the #193 weld; the geometry map | **PASS** |
| T3 | `T3_charge_normalization` | charge = winding: flow + quantization + force check | **PASS** |
| T4 | `T4_einstein_frame` | the Einstein-frame radion: e^{-sqrt3 phi} F^2, V = 0 | **PASS** |
| T5 | `T5_answer` | alpha = 4k^2/rho^2; the modulus; the rank audit | **PASS** |
| T6 | `T6_stabilizer` | the EM cap selects rho*; leptons compatible (decoupling) | **PASS** |
| T7 | `T7_holdout` | the alpha-dependent holdout: one common alpha | **PASS** |
| T8 | `T8_honest_scope` | honest scope | **PASS** |
| T9 | `T9_assessment` | assessment: the arc closes | **PASS** |

## Verdict

**Class:** `CANONICAL_ALPHA_EQUALS_4K2_OVER_RHO2_THE_EINSTEIN_FRAME_RADION_IS_ITS_E_SQRT3_PHI_DILATON_THE_EM_CAP_FIXES_THE_ONLY_ALPHA_COUPLED_DIRECTION_AND_THE_ALPHA_DEPENDENT_HOLDOUT_PASSES_ON_ONE_COMMON_ALPHA`

ESTABLISHED (the argument is in docs/absolute_coupling_capstone.md).

THE CHAIN AND THE MAP. The fiber KK tower (FD = closed form to 2e-11; continuum to 4e-04) welds to the #193 Berger spectra with the Wu-Yang half-charge, and the fiber is EXPLICITLY the committed geometry: R_f = lambda R_u (quadrature 1e-07), half-radius base, and with R_u = r_h = 1 the cap converts 1 model unit = 23.41 l_P.

CHARGE = WINDING. Spectral-flow slope 2k/R_f^2 to 3e-04, EXACT flux quantization (4e-11), and the adiabatic-ramp force check to 3e-04.

THE EINSTEIN FRAME. Symbolically: b = -2a is the frame, -6a^2 the kinetic coefficient (canonical -1/2 at a = 1/(2 sqrt 3)), e^{3b phi} = e^{-sqrt 3 phi} the dilaton coupling, V_tree = 0 EXACTLY, and the dilaton exponent EQUALS the geometric-law exponent: alpha(phi) = alpha(0) e^{sqrt 3 phi} - the alpha law is the dilaton coupling.

THE ANSWER AND THE RANK AUDIT. alpha_k = 4 k^2 (l_P/R_f)^2 EXACTLY (rescale invariance 0e+00): A CONTINUOUS RADION MODULUS REMAINS.  Rank audit: before the cap rank 1 of 5 knobs and grad(alpha) SURVIVES in the flat space; after the cap rank 2 and grad(alpha) is annihilated (2e-16) - the cap fixes exactly the one alpha-coupled direction; the three remaining flats are alpha-decoupled.

THE STABILIZER AND THE LEPTONS. rho* = 2/sqrt(alpha) = 23.4125; guardrail: no closure-constant match (nearest e^pi at 1.2% - rejected). The implied KK scale 5.2e+17 GeV sits 5e+18 above m_mu: leptons are NOT fiber KK modes and need not be - COMPATIBLE BY DECOUPLING, with the #197 ladder slope 1/lambda = 1/R_f (4e-05) and beta_lepton/L_fiber = k5^2 exact.

THE ALPHA-DEPENDENT HOLDOUT. 12 committed constants that GENUINELY carry alpha - linear, quartic, inverse-quartic - all re-derive at machine zero (worst 1e-16); inverted, ONE COMMON ALPHA to 4e-16 relative. The arc closes with its books balanced.
