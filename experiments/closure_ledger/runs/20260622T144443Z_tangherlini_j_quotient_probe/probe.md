# Tangherlini J-quotient consistency and brane non-orientability (PR #169)

**Run:** 2026-06-22T14:44:43+00:00

Identifies the topological root of the non-orientable throat (#167): the antipodal (J) quotient is consistent but dimension-parity-split — orientable on the bulk, non-orientable on the brane mouth. *(QFT on the classical throat, not quantum gravity.)*

- **Bulk**: S³/antipodal = RP³ orientable (det +1)
- **Brane**: S²/antipodal = RP² non-orientable (det −1)
- **Consistency**: one free isometric involution; det = (−1)^{n+1}; metric descends
- **Realization**: J fixes the χ=π/2 brane and restricts to the S² antipodal map (#168)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | goal: the J-quotient parity split (bulk vs brane) | **PASS** |
| T2 | `T2_free_isometric_involution` | J is a free isometric involution (metric descends) | **PASS** |
| T3 | `T3_bulk_orientation_determinant` | bulk orientation determinant +1 → RP³ orientable | **PASS** |
| T4 | `T4_brane_angular_determinant` | brane angular determinant −1 → RP² non-orientable | **PASS** |
| T5 | `T5_dimension_parity_law` | the parity law det = (−1)^{n+1} (odd vs even) | **PASS** |
| T6 | `T6_tangherlini_realization` | Tangherlini realization: J fixes the χ=π/2 brane (#168) | **PASS** |
| T7 | `T7_spin_pin_consistency` | spin/Pin + even-k consistency (#63/#67) | **PASS** |
| T8 | `T8_assessment` | J_QUOTIENT_CONSISTENT_RP3_ORIENTABLE_RP2_NOT | **PASS** |

## The dimension-parity law  (orientation det = (−1)ⁿ⁺¹)

| n | orientation det | RPⁿ orientable? |
|---|---:|---|
| 1 | +1 | yes |
| 2 | -1 | no |
| 3 | +1 | yes |
| 4 | -1 | no |
| 5 | +1 | yes |

(bulk = S³, n=3 → orientable RP³; brane mouth = S², n=2 → non-orientable RP²)

## Verdict

**J_QUOTIENT_CONSISTENT_BULK_RP3_ORIENTABLE_BRANE_RP2_NON_ORIENTABLE.** CONSISTENT, AND IT EXPLAINS THE NON-ORIENTABLE THROAT. The antipodal (J) quotient is a single free isometric involution whose orientation determinant is (−1)^{n+1}, so it acts oppositely on the bulk and the brane mouth.

THE BULK. S³ / antipodal = RP³: the S³ antipodal map preserves orientation (determinant +1), so the BAM spatial bulk remains a consistent ORIENTABLE manifold — RP³ ≅ SO(3).

THE BRANE MOUTH. S² / antipodal = RP²: the S² antipodal map reverses orientation (determinant -1), so the throat mouth is the NON-ORIENTABLE RP² (a cross-cap). This is exactly PR #167's 'non-orientable throat gluing'.

THE CONSISTENCY. The bulk angular sphere (S³, odd) and the throat mouth (S², even) sit one dimension apart, on opposite sides of the parity (−1)^{n+1}: ONE antipodal quotient is therefore consistently orientable on the bulk and non-orientable on the mouth. In the #168 coordinates J = (χ,θ,φ)↦(π−χ,π−θ,φ+π) fixes the equatorial χ=π/2 brane and restricts to the S² antipodal map; the round-angular Tangherlini metric is J-invariant and descends. So the #167 non-orientable throat is the RP² cross-cap inside the orientable RP³ bulk of #168 — forced by the single-dimension drop, and the topological seat of the throat's Pin/spin-½ structure (#63/#67).
