# The Dirac-tower mass ladder: un-dialing the electron (PR #201)

**Run:** 2026-07-03T02:38:27+00:00

The first item of the #200 register: the electron level rebuilt on the #195 index-protected zero mode, the mouth coupling checked against the #185 throat overlap machinery, the #192/#194 diagnostics re-run. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the register item: rebuild the ladder, re-run the diagnostics | **PASS** |
| T2 | `T2_protected_foundation` | foundation: exact zero, forbidden lift, superselection | **PASS** |
| T3 | `T3_rebuilt_fit` | the rebuilt fit (uplift refit; eps1 convention band) | **PASS** |
| T4 | `T4_geometry_of_eps1` | eps1 IS O(1) geometry: neck aspect + soliton-kernel R* | **PASS** |
| T5 | `T5_undialing_bg` | BG 74.7 -> ~4.5; zero sign flips in 2000 +-25% draws | **PASS** |
| T6 | `T6_undialing_berger` | Berger -70.9 -> +4.1; NO lambda_break on (0, inf) | **PASS** |
| T7 | `T7_undialing_mc_and_impossibility` | MC 7.7% -> ~50% typical; hierarchy provably not pairing | **PASS** |
| T8 | `T8_assessment` | the ladder is fully natural; the dial is gone | **PASS** |

## The geometry of the fitted coupling

| convention | eps1 | l/a | R* | R*/RMS |
|---|---:|---:|---:|---:|
| A | 0.011285 | 4.484 | 5.015 | 3.937 |
| B | 0.004836 | 5.332 | 5.593 | 4.391 |

## The rebuilt Barbieri-Giudice table

| parameter | Delta(m_e) | Delta(m_mu) | Delta(m_tau) |
|---|---:|---:|---:|
| `c` | -4.484 | 0.0 | 0.0 |
| `o1` | 1.0 | 0.0 | 0.0 |
| `transport` | -0.011 | -0.011 | 0.001 |
| `pinhole` | 0.648 | 0.648 | 0.039 |
| `resistance` | 0.177 | 0.177 | 0.064 |
| `base` | 0.181 | 0.181 | 0.011 |
| `slope` | 0.01 | 0.01 | -0.001 |
| `beta` | 0.005 | 0.005 | 0.886 |

(worst Delta(m_e) = 4.484 vs the surrogate's 74.7; sign flips in 2000 draws: 0)

## Verdict

**ELECTRON_LEVEL_REBUILT_ON_THE_INDEX_PROTECTED_ZERO_MODE_FINE_TUNING_REMOVED_HIERARCHY_REMAINS_DYNAMICAL_AND_NATURAL.** THE #194 DIAL IS REMOVED (the argument is in docs/dirac_tower_mass_ladder.md; this probe re-runs the diagnostics).

THE REBUILD. The electron level sits on the #195 index-protected zero mode: bare level exactly zero, one-mouth lift forbidden, inter-sector mixing superselected - the surrogate's cancellation structure has no counterpart. The mass is the multiplicative chain eps1*o1*S1; fitting mu/e fixes eps1 = 0.0113 (conv A), and the fitted number IS O(1) geometry: a neck aspect l/a = 4.484-5.332, a mouth separation 3.937-4.391 soliton radii on the actual #180 profile (constrained, not derived).

THE DIAGNOSTICS COLLAPSE. Barbieri-Giudice: worst Delta(m_e) = 4.484 vs the surrogate's 74.7, zero sign flips in 2000 +-25% draws (vs flipping at +-2%). Berger: d ln m_e/d lam = 4.151 (= c - 1/3) vs -70.9, and NO lambda_break on (0, infinity) (vs 0.986). Monte Carlo: P(m_e <= observed) = 0.4963 - typical, smooth, no cliff (vs the 7.7% sliver).

THE DIVISION OF LABOR, FINAL FORM. The impossibility bound (pairing needs eps3/eps1 = 88.6; tunneling gives 1e-4) proves the inter-generation hierarchy cannot be mouth pairing: smallness = geometry (index + pairing), ratios = dynamics (the uplift - fitted, and natural per #194). After the rebuild every sensitivity in the ladder is O(few) or below: the ladder is fully natural. The remaining register item: derive the neck aspect from the 5D core solve.
