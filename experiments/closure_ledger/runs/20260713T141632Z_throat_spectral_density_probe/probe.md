# The throat's spectral density (PR #215)

**Run:** 2026-07-13T14:16:32+00:00

The deliverable is `docs/throat_spectral_density.md` - the flat bank replaced by the geometry's own transmission spectrum. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: the bank retired | **PASS** |
| T2 | `T2_solver` | the greybody solver, flux-conserving | **PASS** |
| T3 | `T3_area_theorem` | IR pinned by the universal area theorem | **PASS** |
| T4 | `T4_photon_sphere` | UV pinned by the photon sphere | **PASS** |
| T5 | `T5_horizon_continuum` | the horizon is the continuum | **PASS** |
| T6 | `T6_quasimode_weld` | quasimodes obey the parameter-free transit law | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | assessment | **PASS** |

## The quasimode ladder vs the transit law

| w_q | gamma (measured) | T(w)/(2L) (predicted) | ratio |
|---:|---:|---:|---:|
| 0.23188 | 4.117e-04 | 4.143e-04 | 0.994 |
| 0.33472 | 1.476e-03 | 1.457e-03 | 1.013 |
| 0.43635 | 3.913e-03 | 3.623e-03 | 1.080 |

## Verdict

**Class:** `THE_BANK_RETIRED_THE_SPECTRAL_DENSITY_IS_THE_THROATS_GREYBODY_IR_TRANSPARENT_BY_THE_AREA_THEOREM_UV_BLACK_ABOVE_THE_PHOTON_SPHERE_AND_THE_HORIZON_IS_THE_CONTINUUM`

DERIVED (the argument is in docs/throat_spectral_density.md).

THE SPECTRUM FROM GEOMETRY. The 5D Tangherlini greybody T_l(w), flux-conserving across the band: the IR wing is the universal area theorem - sigma_s/A_h = 1.012 at w r_h = 0.04, extrapolating to 0.983, with T ~ w^3 (slope 3.02) - and the UV wing is the photon sphere: r_c = sqrt(2) r_h, b_c = 2 r_h, T = 1/2 at the eikonal (l+1)/b_c (ratios 1.04/1.02/1.01), T -> 1 above.

THE CONTINUUM. The tortoise depth grows as -(r_h/2) ln delta (coefficient -0.5005): level spacing -> 0, recurrence -> infinity - the #214 order of limits is a property of the horizon, not a choice.

THE WELD. Cavity quasimodes obey the parameter-free transit law gamma = T(w)/(2 L_cav): ratios 0.994/1.013/1.080 down the ladder, spacing pi/L_cav to 4%. The tower density (r_h/R = 0.1): per-pass absorption 1.7e-03 at n = 1 rising as n^3 (ratio 9.1/8 expected) to 1.000 at n = 24: IR-transparent, UV-black. Nothing about the absorber is chosen anymore - spectrum, address, and epsilon are all read off the frozen bulk.
