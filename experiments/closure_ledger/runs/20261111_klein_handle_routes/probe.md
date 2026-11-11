# Two routes to the sheet pairing, on equal terms

Round 14's lattice electrodynamics on the doubly reversed seam. Trapped flux and a matter current compared at matched total energy, over 12 draws at each of two sizes, with gauge-invariant observables only. Classical throughout; lattice units; no norm, detector law or statement about which route BAM uses.

Public freeze: `fd78693`. Seed `2026111101`. Long window `[0, 200]`.

| verdict field | value |
|---|---|
| numerical_validity | **PASS** |
| flux_route_requirement | **FLUX_ROUTE_NEEDS_FIELD_MOMENTUM** |
| energy_scaling | **MATTER_LINEAR_FLUX_THREE_HALVES** |
| short_time_ranking | **MATTER_CURRENT_DOMINATES_BY_INVERSE_SQRT_ENERGY** |
| flux_fraction | **OPTIMUM_AT_ONE_THIRD_IN_THE_FLUX** |
| coupling_dependence | **MATTER_ROUTE_IS_TRANSPORT_FLUX_ROUTE_IS_ELECTRIC** |
| electric_sectors | **ROUTES_OCCUPY_DIFFERENT_ELECTRIC_SECTORS** |
| interference | **NO_FIRST_ORDER_INTERFERENCE_RELATIVE_SIGN_PHYSICAL** |
| screening_rule | **COUNTEREXAMPLE** |
| screening_leading_order | **BOND_ORDER_FIXES_THE_SIGN_AND_LEADING_ORDER_OF_EARLY_FLUX_CHANGE** |
| reversal_requirements | **MATTER_ROUTE_NEEDS_ONLY_FIBRE_REVERSAL** |
| long_time_ranking | **FLUX_ROUTE_GAINS_AT_LONG_TIMES** |
| long_window_numerics | **STATISTICS_STABLE_UNDER_STEP_AND_INTEGRATOR** |
| support_dynamics | **EFFECTIVE_SUPPORT_ASSUMED** |
| orienting_datum | **NOT_DERIVED** |
| positive_response | **NOT_DERIVED** |
| event_law | **NOT_DERIVED** |

**Frozen prediction contradicted (screening_rule).** Frozen P8 required the early flux change to match -1/2 g^2 Phi (B_0/|C|) t^2 within 5 % at t = 0.1 for every resolved draw. One replicate draw misses by 6.4 %. The sign rule holds in every resolved draw. The formula is the leading order: its error falls as t^2 and is below 1 % at t = 0.04. The tolerance was mis-set, and so was my first explanation of the miss.

## Short window: sheet charge at matched energy

| size | energy | matter current | at rest | flux route | ratio F/M (median, range) |
|---|---:|---:|---:|---:|---|
| primary | 0.00025 | 3.707e-06 | 5.668e-06 | 3.579e-08 | 0.0054 (0.0023 to 0.0329) |
| primary | 0.001 | 1.483e-05 | 2.267e-05 | 2.863e-07 | 0.0108 (0.0045 to 0.0658) |
| primary | 0.004 | 5.932e-05 | 9.069e-05 | 2.290e-06 | 0.0215 (0.0091 to 0.1314) |
| primary | 0.016 | 2.373e-04 | 3.628e-04 | 1.831e-05 | 0.0426 (0.0181 to 0.2611) |
| replicate | 0.00025 | 4.408e-06 | 5.240e-06 | 7.809e-09 | 0.0023 (0.0006 to 0.0099) |
| replicate | 0.001 | 1.763e-05 | 2.096e-05 | 6.247e-08 | 0.0047 (0.0012 to 0.0198) |
| replicate | 0.004 | 7.054e-05 | 8.384e-05 | 4.998e-07 | 0.0094 (0.0024 to 0.0395) |
| replicate | 0.016 | 2.822e-04 | 3.353e-04 | 3.999e-06 | 0.0188 (0.0048 to 0.0788) |

| size | slope M | slope M0 | slope F | slope of ratio |
|---|---:|---:|---:|---:|
| primary | 1.0000 | 1.0000 | 1.4998 | 0.4972 |
| replicate | 1.0001 | 1.0000 | 1.5000 | 0.4998 |

## Long window (measurements, except the ranking)

| size | energy | rms M | rms F | rms B | ratio long | ratio short | F larger | flux retained (min) | current half-life |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| primary | 0.00025 | 7.213e-06 | 2.429e-06 | 3.861e-06 | 0.318 | 0.0054 | 0 of 12 | 0.9997 | 0.46 |
| primary | 0.001 | 2.969e-05 | 1.043e-05 | 1.510e-05 | 0.282 | 0.0108 | 0 of 12 | 0.9974 | 0.46 |
| primary | 0.004 | 1.246e-04 | 5.357e-05 | 6.024e-05 | 0.422 | 0.0215 | 0 of 12 | 0.9986 | 0.46 |
| primary | 0.016 | 4.656e-04 | 2.265e-04 | 2.774e-04 | 0.485 | 0.0426 | 0 of 12 | 0.9901 | 0.46 |
| replicate | 0.00025 | 6.229e-06 | 1.627e-06 | 2.800e-06 | 0.320 | 0.0023 | 0 of 12 | 0.9998 | 0.48 |
| replicate | 0.001 | 2.509e-05 | 8.001e-06 | 1.188e-05 | 0.341 | 0.0047 | 0 of 12 | 0.9996 | 0.48 |
| replicate | 0.004 | 1.011e-04 | 3.463e-05 | 4.492e-05 | 0.347 | 0.0094 | 0 of 12 | 0.9992 | 0.48 |
| replicate | 0.016 | 3.575e-04 | 1.388e-04 | 2.046e-04 | 0.352 | 0.0188 | 1 of 12 | 0.9928 | 0.48 |

| Check | Pass |
|---|---|
| P1 flux with the field at rest gives no sheet charge | True |
| P1 the flux still evolves, and a moving field does give one | True |
| P2 matter-current sheet charge is linear in energy | True |
| P2 flux-route sheet charge goes as energy to the three halves | True |
| P3 the ratio grows as the square root of energy | True |
| P3 the median ratio is below one at all four energies | True |
| P4 ratios to the one-third value are within 10 % | True |
| P5 matter route does not depend on the coupling | True |
| P5 flux route is first order in the coupling | True |
| P5 seam flux is g times the sheet charge | True |
| P6 matter current at rest has no loop flux, holonomy, loop current or swap-even field | True |
| P6 it does have sheet charge, and the flux route has loop flux | True |
| P7 additivity defect is second order in the flux | True |
| P7 reversing the flux alone changes the history | True |
| P7 reversing both orienting data is a symmetry | True |
| P8 sign of the early flux change is minus the sign of the bond order | True |
| P8 early flux change within 5 % of the formula (frozen tolerance) | False |
| post-freeze: the formula's error falls as t squared | True |
| post-freeze: the error is below 1 % at t = 0.04 | True |
| post-freeze: the rate of the bond order vanishes by parity | True |
| post-freeze: the sign rule holds at all three times | True |
| P9 matter current at rest is fixed by the time-preserving fibre-reversing deck | True |
| P9 its history stays in that fixed set and separates charge | True |
| P9 ordinary handle keeps the sheet charge at zero | True |
| P10 the long-window ratio exceeds the short-window ratio at every energy | True |
| long-window Gauss law holds | True |
| rms sheet charge is stable under a halved step | True |
| rms sheet charge agrees with DOP853 | True |

27 of 28 checks hold; 1 frozen prediction check contradicted, reported above.

Matched total energy is one definition of equal terms. Neither orienting datum is derived, and nothing here is a probability.
