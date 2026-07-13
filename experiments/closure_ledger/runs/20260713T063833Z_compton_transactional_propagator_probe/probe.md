# The transactional Compton propagator (PR #213)

**Run:** 2026-07-13T06:38:33+00:00

The deliverable is `docs/compton_transactional_propagator.md` - the Feynman propagator, its two completion orderings, and their coherent relative phase from the frozen bulk. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: derive the imported time structure | **PASS** |
| T2 | `T2_transactional_construction` | G_F = Gbar + unique absorber response | **PASS** |
| T3 | `T3_two_orderings` | the two orderings -> the exact covariant pole | **PASS** |
| T4 | `T4_frozen_bulk_conditions` | static + closed: both conditions are geometric | **PASS** |
| T5 | `T5_coherent_phase` | the relative phase forced; deform test refutes others | **PASS** |
| T6 | `T6_compton_denominators` | the Compton denominators derived, not assumed | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | assessment | **PASS** |

## Verdict

**Class:** `FEYNMAN_PROPAGATOR_FROM_COMPLETE_HISTORIES_THE_TWO_ORDERINGS_ARE_OFFER_AND_CONFIRMATION_THE_COHERENT_PHASE_FORCED_BY_THE_FROZEN_BULK`

DERIVED (the argument is in docs/compton_transactional_propagator.md).

THE CONSTRUCTION. The static bulk forces the time-symmetric Wheeler-Feynman field; the closed S^3 makes complete absorption a geometric fact (antipodal refocus concentration 0.50 at t = pi R vs 1e-28 at mid-flight, full return at 2 pi R); the complete-history boundary condition is a well-posed 2x2 system (uniqueness error 0e+00) whose unique solution is the Feynman propagator, G_F = Gbar + (i/2w) cos(wt).

THE TWO ORDERINGS. The theta(+t)/theta(-t) split of G_F is offer and confirmation; each half-line transform is an OFPT energy denominator - individually regulator-dependent (truncation spread 0.50) and non-covariant (evenness violation 0.21) - and the coherent sum is the exact covariant pole (identity to 7e-15).

THE PHASE. The t -> -t isometry (tower identity to 0e+00) forces relative phase 0; deforming it breaks evenness and the pole form at O(1) and fails the Compton denominators (min error 0.07 at phi = pi/2). The s- and u-channel denominators the Compton arc assumed are recovered exactly (9e-15 / 1e-14): the i*epsilon prescription is, on this geometry, a theorem about complete histories.
