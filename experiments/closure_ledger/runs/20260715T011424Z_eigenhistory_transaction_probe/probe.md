# The eigenhistory transaction (PR #218)

**Run:** 2026-07-15T01:14:24+00:00

The deliverable is `docs/eigenhistory_transaction.md` - the transaction as a homogeneous eigenhistory with the amplitude fixed by energy and state closure. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: the existence theorem | **PASS** |
| T2 | `T2_energy_closure_leg` | energy closure: lossless on resonance, reactive source | **PASS** |
| T3 | `T3_existence` | existence: IVT + the null space opens at the point | **PASS** |
| T4 | `T4_energy_dynamics` | energy conserved exactly; 1e4-pass persistence | **PASS** |
| T5 | `T5_state_closure` | state closure fixes a discrete amplitude spectrum | **PASS** |
| T6 | `T6_driven_correspondence` | the driven pole = the eigenhistory (Novikov populated) | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | assessment | **PASS** |

## Verdict

**Class:** `A_FINITE_SOURCE_INCLUSIVE_ENERGY_CONSERVING_SELF_CONSISTENT_WORMHOLE_TRANSACTION_EXISTS_THE_EIGENHISTORY_AMPLITUDE_FIXED_BY_ENERGY_AND_STATE_CLOSURE`

ESTABLISHED (the argument is in docs/eigenhistory_transaction.md).

THE FORMULATION. F = Lambda_tot(w, I) F, no source term: the transaction is a homogeneous, globally constrained eigenhistory. ENERGY CLOSURE: the throat is lossless exactly ON its interior resonance (unitarity identity, exact tier 2e-16; Tangherlini deficit 2.0e-04 = 10x the port flux error - solver precision), the source reactive (0e+00).

THE THEOREM. The IVT range argument guarantees the phase target 3.4201 is reachable (source range 12.0 > 2 pi); the amplitude solves to I* = 1.5944 with Lambda_tot - 1 = 4e-15; the null space of the transfer system opens EXACTLY there (singular value 6e-16 at the eigenhistory vs 0.14/0.27 detuned). A finite, source-inclusive, energy-conserving, self-consistent wormhole transaction EXISTS.

THE DYNAMICS. Energy conserved pass by pass (4e-16); the eigenhistory persists 10^4 passes (amplitude drift 2e-15); the interior mouth state finite (3.0x resonant enhancement); the physical tier decays at exactly its solver deficit. STATE CLOSURE fixes the amplitude: a discrete branch spectrum [1.5944455007897045, 16.899064794055423], perturbed amplitudes dephasing at d phi_s/dI (1.5250 vs 1.5337); the source's internal-state shift 3.4201 is part of the closed history.

THE CORRESPONDENCE. The #217 driven pole sits exactly at the eigenhistory (response 2.3e+14 there vs 1.9 detuned); weak driving shadows the homogeneous solution. The marginal Novikov point is populated: the completed transaction is an explicit self-consistent history.
