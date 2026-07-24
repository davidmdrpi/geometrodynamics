# Radion stabilization from the primordial EM-capped throat (PR #226)

**Run:** 2026-07-24T07:26:19+00:00

The deliverable is `docs/radion_stabilization.md` - V_eff(phi), alpha at the minimum, and the radion mass. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: execute #225's open dynamical problem | **PASS** |
| T2 | `T2_cap_imported` | the #55 cap imported; alpha-dependent re-read | **PASS** |
| T3 | `T3_radion_charges` | the radion charges: e^{+sqrt3 phi} vs e^{-sqrt3 phi} | **PASS** |
| T4 | `T4_potential` | the dyonic minimum; the no-go; m^2 = (35/12) V_min | **PASS** |
| T5 | `T5_alpha` | alpha* = sqrt(5 kappa/28); Dirac point 0.4226 | **PASS** |
| T6 | `T6_radion_mass` | the radion mass: E'' = 7 U_el*; anchored scales | **PASS** |
| T7 | `T7_arc_consistency` | the stabilizer lands in the reserved rank-audit slot | **PASS** |
| T8 | `T8_honest_scope` | honest scope | **PASS** |
| T9 | `T9_assessment` | assessment | **PASS** |

## Verdict

**Class:** `THE_HOPF_CHARGE_STABILIZES_THE_RADION_V_EFF_HAS_A_MINIMUM_WITH_ALPHA_STAR_EQUALS_SQRT_5KAPPA_OVER_28_ORDER_ONE_AT_THE_DIRAC_POINT_AND_M_PHI2_EQUALS_35_OVER_12_V_MIN_THE_OBSERVED_ALPHA_NEEDS_KAPPA_3E_MINUS_4_THE_NAMED_OPEN_PROBLEM`

ESTABLISHED (the argument is in docs/radion_stabilization.md).

THE CAP. The #55 ledger re-read at machine zero (A = alpha hbar c/2, R* = (A/2B)^{1/3}, U/mc^2 = alpha/2, E'' = 6B > 0) - the alpha-dependent holdout continues.

THE CHARGES. Fixed charge -> e^{+sqrt3 phi} (4e-06), fixed Hopf flux -> e^{-sqrt3 phi} (4e-06), Dirac ratio 1/(4 alpha^2) exact; primordial radius (#222) adds e^{a phi}: exponents 7/(2 sqrt3) and -5/(2 sqrt3) - ONE OF EACH SIGN.

THE POTENTIAL. The dyonic minimum exists (numeric = closed form to 1e-08); m_phi^2 = (35/12) V_min holds numerically (7e-09) and symbolically; WITHOUT the Hopf charge V is monotonic - runaway to decompactification: the #58 topology saves the vacuum.

ALPHA. alpha*^2 = 5 kappa/28: at the Dirac point alpha* = 0.4226 (order one, 58x observed); observed alpha needs kappa = 2.98e-04 (verified by minimization; guardrail: no match, nearest alpha/k5^2 at 2.1%); cohesion tilt < 0.2%.

THE MASS. U_mag/U_el = 7/5 at the minimum, E'' = 7 U_el* exactly (9e-09); anchored: 1.3e+16 GeV per throat at the observed-alpha minimum (fiber anchor; GUT-scale - the radion is heavy), 13 keV on the #55 Compton anchor; m_phi^2 = 16 pi n E''/m_P^2.

THE ARC. The stabilizer row = the reserved cap row (rank 2, grad alpha annihilated 2e-16); rho at the Dirac point 3.08, at observed alpha 23.41; #222's x210 primordial exclusion re-read. Nothing else in the arc moves.
