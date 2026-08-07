# Probe M — the charge-conserving apparatus and the complete pointer statistics (PR #240)

_Run 2026-08-06T06:53:22.208966+00:00 · 8/8 PASS_

**Q. Must the winding carrier be replaced by the spin doublet?**

**A. No — the 2.33 ceiling was an artifact too.**

## The apparatus conserves the *right* charge

| `n_max` | `‖[H, e^{2πiK/8}]‖` | `‖[H, K_integer]‖` |
|---:|---:|---:|
| 8 | **1.1e-15** | 9.2 |
| 16 | **1.1e-15** | 13.1 |
| 32 | **6.9e-15** | 18.5 |

Winding on a discrete `N_χ = 8` fiber is a **Z₈** charge. #239's integer-`K` test passed only because it excluded the wrap-around transitions.

## Every channel is detected — there is no no-click

The orbit under `Δk = ±2` is the 4-cycle `[-3, -1, 1, 3]`, and a winding Stern–Gerlach resolves all four (min separation 4.0σ, overlap 4.6e-02).

## The ceiling disappears

| `\|t\|` | `sign(k)` binning | #239 no-click binning | best of 16 |
|---:|---:|---:|---:|
| 0.4 | **2.134677** | 2.116652 | 2.134677 |
| 0.8 | **2.457397** | 2.305325 | 2.457397 |
| 1.2 | **2.758947** | 2.330532 | 2.758947 |
| 1.5 | **2.824508** | 2.327297 | 2.824508 |
| 1.9 | **2.828028** | 2.331306 | 2.828028 |

Tsirelson is 2.828427. *The 2.33 ceiling was the cost of throwing away channels the apparatus resolves.*

## Finite-carrier back-action (exact)

| `n̄` | max CHSH | coherent truncation |
|---:|---:|---:|
| 1 | 2.290317 | 0.0e+00 |
| 4 | 2.601231 | 1.0e-15 |
| 16 | 2.760053 | 2.5e-13 |

## Verdict

**THE_APPARATUS_IS_CHARGE_CONSERVING_UNDER_THE_Z8_OPERATOR_AND_THE_COMPLETE_POINTER_STATISTICS_REACH_TSIRELSON_BECAUSE_THE_LEAKED_CHANNELS_ARE_DETECTED_SO_THE_2_33_CEILING_WAS_AN_ARTIFACT_AND_THE_WINDING_CARRIER_NEED_NOT_BE_REPLACED**
