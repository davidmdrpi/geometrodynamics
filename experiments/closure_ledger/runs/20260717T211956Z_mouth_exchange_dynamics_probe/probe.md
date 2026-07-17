# Mouth-exchange dynamics (PR #224)

**Run:** 2026-07-17T21:19:56+00:00

The deliverable is `docs/mouth_exchange_dynamics.md` - the #223 successor: the exchange observables on the two-throat network. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: the five exchange observables | **PASS** |
| T2 | `T2_network` | the two-throat network; genuine mouth basins | **PASS** |
| T3 | `T3_p_other_mouth` | P_other(t) = P_max sin^2(dw t/2), exact + evolved | **PASS** |
| T4 | `T4_max_and_period` | P_max = localization weight; period = pi/dw | **PASS** |
| T5 | `T5_coupling_law` | J(r_s) -> (w r_s)^4: exchange frozen at the anchor | **PASS** |
| T6 | `T6_asymmetric_survival` | Rabi/Lorentzian survival: localization by asymmetry | **PASS** |
| T7 | `T7_complex_width` | Gamma = 0 closed; lead-independent QNM width open | **PASS** |
| T8 | `T8_honest_scope` | honest scope | **PASS** |
| T9 | `T9_assessment` | assessment | **PASS** |

## Verdict

**Class:** `THE_MOUTH_EXCHANGE_IS_A_TEXTBOOK_TWO_LEVEL_BEAT_WITH_PERIOD_PI_OVER_DW_ASYMMETRY_LOCALIZES_BY_THE_RABI_LAW_AND_THE_WIDTH_IS_ZERO_UNTIL_A_CHANNEL_OPENS`

ESTABLISHED (the argument is in docs/mouth_exchange_dynamics.md).

THE SYSTEM. The genuine mouth-to-mouth two-level system is the two-throat network (a single bridge ring has ONE exterior cavity - stated as a finding): the working interior doublet has interior fraction 0.965, half-half eigenmodes (0.500), and basins with 0.977 localization.

P_OTHER-MOUTH(t). The exact two-mode regional beat is P_max sin^2(dw t/2) IDENTICALLY (deviation 1e-16); direct evolution: transfer max 0.9768 at t/T = 1.0006, A depleted to 0.0000, half-period point at 0.501 of max.

THE TWO NUMBERS. P_max = 0.9768 (eigen) = 0.9768 (evolution) - limited only by the dressing tail; exchange period pi/dw = 1833 vs evolution 1834.

THE COUPLING LAW. J(r_s) exponents 3.26/3.54/3.68 rising toward the #223 limit 4: at the primordial anchor the period extrapolates to 1e+11 - mouth-to-mouth transfer is dynamically frozen.

THE SURVIVAL LAW. Clock-rate asymmetry: the Rabi identity to 1.0%, the Lorentzian to 0.1%; at a 1.5% bias survival = 0.986: LOCALIZATION BY ASYMMETRY.

THE WIDTH. Compact network: Gamma = 0 exactly (persistence). Open channel: QNM doublet Gamma = 4.0e-04/4.2e-04 (lead-independent to 0.0%); strong coupling J/(Gamma/4) = 8.3; Gamma_direct/Gamma_pair = 2.09 (the width shared by hybridization).
