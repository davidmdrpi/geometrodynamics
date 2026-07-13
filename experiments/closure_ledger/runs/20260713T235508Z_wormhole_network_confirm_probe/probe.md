# The wormhole-network confirmation (PR #216)

**Run:** 2026-07-13T23:55:08+00:00

The deliverable is `docs/wormhole_network_confirmation.md` plus `geometrodynamics/transaction/network.py` - the advanced confirmation as an explicit, globally causal network traversal. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: a mechanism for the advanced half | **PASS** |
| T2 | `T2_traversal_ledger` | the ten-step ledger: locally forward, globally past | **PASS** |
| T3 | `T3_greybody_amplitude` | the complex greybody amplitude + Wigner delay | **PASS** |
| T4 | `T4_advanced_projection` | the projection is advanced; phase closure = comb | **PASS** |
| T5 | `T5_live_packet` | the live packet intersects the crossing | **PASS** |
| T6 | `T6_pole_structure` | the pole at closure; detuning = the #213 deform | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | assessment | **PASS** |

## The traversal ledger

| leg | global start | global end | local duration |
|---|---:|---:|---:|
| offer: source -> mouth A | 0.000 | 3.142 | 3.142 |
| throat | 3.142 | -3.142 | 0.800 |
| return: mouth B -> source | -3.142 | 0.000 | 3.142 |

## Verdict

**Class:** `THE_ADVANCED_HALF_IS_A_FUTURE_DIRECTED_NETWORK_TRAVERSAL_THE_CLOCK_OFFSET_MOUTH_IS_THE_ANALYTIC_CONTINUATION_AND_PHASE_CLOSURE_IS_THE_COHERENCE_CONDITION`

DERIVED (the argument is in docs/wormhole_network_confirmation.md).

THE MECHANISM. One retarded C-wave: source -> future antipode (t = 3.142) -> greybody transmission -> forward throat traversal -> emergence at t = -3.142 (the global PAST, via the frozen mouth offset Delta_BA = -7.083) -> forward return -> the original crossing (t_return = 0e+00). Every leg locally future-directed; the derived transfer factor U_BA = e^{i w Delta_BA} (identity to 3e-15).

THE PROJECTION IS ADVANCED. The absorption -> return kernel is Lambda(w) e^{-i w D'} with D' < 0 and |Lambda| = sqrt(T): the retarded rule continued to a negative exterior interval - the advanced kernel with the #215 greybody weight. Phase closure arg t + w Delta = 0 (mod 2 pi) selects a comb (spacing 2 pi/|Delta + tau_W| to 0%); the engine's retro_phase_match accepts the loop at the comb (1.000) and rejects it detuned (0.03). Live packet: emerges at -3.05 (< 0), intersects the crossing at 0.088 (the Wigner delay 0.085), energy fraction 0.9998 = mean T.

THE POLE. At closure I+ + Lambda I- is the covariant pole (|Lambda| = 1.0000, phase 0.000); detuning the network IS the #213 deform test (phi = w delta exactly; deviations match the #213 diagnostic to 2%); below the barrier the confirmation is partial with the coherent share (1+sqrt T)/2. The two Compton completion families are two globally causal wormhole-network path classes.
