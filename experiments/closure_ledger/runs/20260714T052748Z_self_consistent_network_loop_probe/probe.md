# The self-consistent network loop (PR #217)

**Run:** 2026-07-14T05:27:48+00:00

The deliverable is `docs/self_consistent_network_loop.md` - the full two-mouth transfer system solved and the effective Green function derived before comparison with I+-. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: solve the loop, then compare | **PASS** |
| T2 | `T2_clock_rates` | clock-rate-correct traversal (redshift fixed) | **PASS** |
| T3 | `T3_transfer_system` | the transfer system; ring-validated eigenvalue | **PASS** |
| T4 | `T4_transaction_points` | Wigner-correct closure; the transaction points | **PASS** |
| T5 | `T5_state_evolution` | source and mouth state evolution | **PASS** |
| T6 | `T6_geff_vs_orderings` | G_eff first, then I+-: lines, passivity, epsilon | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | assessment | **PASS** |

## Verdict

**Class:** `THE_LOOP_SOLVED_SELF_CONSISTENTLY_G_EFF_DERIVED_COMPLETED_TRANSACTIONS_SIT_ON_THE_THROAT_RESONANCES_AND_THE_NOVIKOV_FIXED_POINT_IS_PASSIVE`

DERIVED (the argument is in docs/self_consistent_network_loop.md).

THE SYSTEM. Both barriers in one signal-flow system: the resolvent equals g/(1 - Lambda) exactly (9e-15); the interior x winding double series resums nested (1e-15); and the value-transport eigenvalue Lambda = t_net e^{i w D_loop} is validated against a first-principles ring spectrum (modes to 0.000; the opposite convention misses by 0.12) - correcting the #216 closure bookkeeping.

THE CLOSURES. Clock rates exact (redshift fixed; elastic confirmation forces rate-matched mouths). Wigner-correct closure lands the packet at 0.006 (uncorrected: 0.175 = tau_W), phase-aligned (-0.005). Above the barrier the mouth phase closes the carrier (Lambda = 0.999998); below it BOTH closures are solved by the throat alone and the completed transaction sits ON the interior resonance (|t_net|^2 = 0.969 at tau* = 5.1105, residual 5e-13).

THE STATES. Input-output rates match the derived linewidth (kappa = 0.0652 vs FWHM 0.0520) and the resonant storage time (31.1 vs 30.7); the source iteration converges at rate |Lambda| (0.2324 vs 0.2324) - marginal at completion; the cavity builds to steady state.

G_EFF, THEN I+-. Quartic line at the completed transaction (FWHM 0.0118 vs 0.0118 - group closure makes the resonance anomalously flat), Lorentzian at carrier-only closure (0.0485 vs 0.0456); passivity max|Lambda| = 0.999999 <= 1 (Novikov, no runaway); the O(Lambda) truncation IS the #216 assembly, the resummation renormalizes the confirmation weight to Lambda/(1-Lambda) (3e-14; at completion the weight reaches 5.0e+05 - the divergence IS the transaction pole), and the confirmation deficit 1-|Lambda| -> 0 at completion (2.0e-06 vs partial 0.11): the completed transaction is the #213 coherent limit, its deficit the #214 absorber damping in transactional form.
