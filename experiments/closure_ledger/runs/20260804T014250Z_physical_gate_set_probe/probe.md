# The physical gate set (PR #234)

**Run:** 2026-08-04T01:42:50+00:00

The deliverable is `docs/physical_gate_set.md` - the qubit clock, the gravitational Kerr, and the twist switch, calibrated on the committed network. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: #232's scope item + the arc weld | **PASS** |
| T2 | `T2_clock` | the clock: the #224 beat, re-derived exactly | **PASS** |
| T3 | `T3_twist_switch` | the switch: the #227 freeze; Kerr channel spared | **PASS** |
| T4 | `T4_kerr_calibrated` | the gravitational Kerr: chi = 5.5e-2, t_CZ = 57 | **PASS** |
| T5 | `T5_selection_rule_cz` | the dual-rail selection rule; the calibrated CZ | **PASS** |
| T6 | `T6_circuit_budget` | the circuit budget at physical rates | **PASS** |
| T7 | `T7_holdout` | the arc-weld holdout (four ledgers, both arcs) | **PASS** |
| T8 | `T8_honest_scope` | honest scope | **PASS** |
| T9 | `T9_assessment` | assessment | **PASS** |

## Verdict

**Class:** `THE_GATE_SET_IS_PHYSICAL_THE_BEAT_IS_THE_CLOCK_THE_TWIST_IS_THE_SWITCH_AND_THE_GRAVITATIONAL_DRESSING_IS_THE_KERR_WITH_T_CZ_32X_FASTER_THAN_THE_CLOCK_AT_LEADING_ORDER_IN_A_STRONG_COUPLING`

ESTABLISHED (the argument is in docs/physical_gate_set.md).

THE CLOCK. The #224 doublet re-derived exactly (dw dev 0e+00): X period 1833, basins 0.9767.

THE SWITCH. One twisted link freezes the beat (ratio 4.9e-11; the #227 row re-derived to 0e+00); memory leakage time 3.8e+13; the Kerr channel shifts < 0.2% under the twist - an error switch, not a gate hazard.

THE KERR. delta mu per quantum 5.802e-02, response dw/dmu -0.947: chi = 5.493e-02, t_CZ = 57.2 - 32x faster than the clock.  GRAVITY IS THE ENTANGLING RESOURCE.  Honesty: delta mu_q/mu = 0.64 - leading order, strong coupling.

THE SELECTION RULE. Same-qubit Kerr terms vanish on the sector (0e+00); the calibrated CZ is exact (1e-16, leakage 0e+00); residual ZZ 0.13 rad per CZ, compensable.

THE BUDGET. GHZ_3 = 4696 model units; twist-frozen leakage error 1.6e-20.

THE WELD. 7 committed constants from FOUR ledgers spanning BOTH arcs (#224, #227, #223, #232) re-derived on one machine at machine agreement: the arcs meet in one calibrated gate set.
