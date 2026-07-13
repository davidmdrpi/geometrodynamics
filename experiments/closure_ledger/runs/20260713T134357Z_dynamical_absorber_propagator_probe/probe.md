# The dynamical absorber (PR #214)

**Run:** 2026-07-13T13:43:57+00:00

The deliverable is `docs/dynamical_absorber_propagator.md` - the absorber promoted from a boundary condition to a degree of freedom. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: absorption as an outcome | **PASS** |
| T2 | `T2_exact_elimination` | exact elimination of the Gaussian absorber | **PASS** |
| T3 | `T3_derived_damping` | the damping derived: pole = live decay = rule | **PASS** |
| T4 | `T4_placement_adjudication` | the antipode is the impedance-matched absorber | **PASS** |
| T5 | `T5_recurrence` | finite bank = recurrence; the order of limits | **PASS** |
| T6 | `T6_213_closure` | the #213 epsilon filled in as a physical rate | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | assessment | **PASS** |

## Verdict

**Class:** `EPSILON_IS_THE_ABSORBER_RESPONSE_RATE_AND_THE_ANTIPODE_IS_THE_IMPEDANCE_MATCHED_ABSORBER_OF_THE_CLOSED_BULK`

DERIVED (the argument is in docs/dynamical_absorber_propagator.md).

THE PROMOTION. S = S_field + S_absorber + S_coupling, all Gaussian: the elimination is exact (resolvent identity to machine precision), and #213's imposed epsilon is now the computed damping rate of the absorber response - the complex pole, the live time-domain decay, and the golden rule agree (0.0499 / 0.0499 / 0.0500), scaling exactly as g^2 (slope 1.9996), with g -> 0 recovering the #213 kernel.

THE ADDRESS. The geometry adjudicates placement: the antipodal point absorber couples as Y_n ~ n - exactly cancelling the kinematic suppression - so EVERY tower mode damps at the same rate g^2/(4 pi); live residuals 0.014 (antipode, uniform across modes) vs 0.95 (generic point, near-zero couplings) vs 1.00 (uniform distribution - an exact selection rule confines it to the ground mode). The #166 focusing caustic is the impedance matching.

THE MEANING OF EPSILON. A finite bank revives at T_rec = 2 pi/dnu (linear in N, measured ratio 1.90): on the closed bulk, finite epsilon is Poincare recurrence, and the eps -> 0 limit is the continuum limit of the throat's internal spectrum taken first. The orderings' i*epsilon is filled in as a physical rate (quadrature error 3e-08), and the covariant pole form agrees to O(gamma^2).
