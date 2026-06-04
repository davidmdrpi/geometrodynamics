# The 5D Tangherlini bulk lift (PR #127)

**Run:** 2026-06-04T04:27:43+00:00

Lifts the PR #116 Tangherlini cavity operator to its explicit 5D parent metric, verifies the throat is the boundary trace of a genuine D=5 vacuum geometry, and reconciles that bulk with the AdS₅ / Randall–Sundrum brane bulk (PR #57). Curvature is computed by a self-contained numerical GR routine (no symbolic algebra).

- **Metric**: ds² = −f dt² + f⁻¹dr² + r²dΩ₃², f = 1 − (rs/r)² (D=5)
- **Ricci**: R_μν = 0 (vacuum, Λ = 0), asymptotically flat
- **Kretschmann**: K = 72 rs⁴/r⁸, regular on cavity, singularity only at r=0
- **Horizon**: throat = S³ horizon at r = rs = R_MID (= Hopf base)
- **Hawking**: T_H = 1/(2π rs), carries the closure quantum 2π
- **AdS₅ lift**: f = 1 − rs²/r² + k²r², Einstein (Λ₅ = −6k²), interpolates neck → AdS₅/RS
- **Open**: exact AdS scale k (= κ₅²/Λ₅, PR #112); global brane-localised solution

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | lift the PR #116 cavity to its 5D bulk; reconcile with AdS₅/RS | **PASS** |
| T2 | `T2_metric_and_horizon` | 5D metric f = 1 − (rs/r)²; throat = S³ horizon (Hopf base) | **PASS** |
| T3 | `T3_ricci_flat_vacuum` | Ricci-flat vacuum R_μν = 0, Λ = 0 (asymptotically flat) | **PASS** |
| T4 | `T4_kretschmann_regular_on_cavity` | Kretschmann K = 72 rs⁴/r⁸ regular on cavity; sing. only at r=0 | **PASS** |
| T5 | `T5_cavity_potential_descends_from_D5` | cavity potential coefficients = D=5 reductions ⟹ k₅ = 5 | **PASS** |
| T6 | `T6_hawking_temperature_carries_2pi` | Hawking T_H = 1/(2π rs) carries the closure quantum 2π | **PASS** |
| T7 | `T7_ads5_rs_reconciliation` | AdS₅ lift f = 1 − rs²/r² + k²r² Einstein; interpolates; k open | **PASS** |
| T8 | `T8_assessment` | FIVE_D_TANGHERLINI_BULK_LIFT_RICCI_FLAT_VACUUM_CAVITY_REGULAR | **PASS** |

## Kretschmann scalar on the cavity (K = 72 rs⁴/r⁸, regular)

| r | K (numeric) | 72 rs⁴/r⁸ | match |
|---:|---:|---:|:---:|
| 1.05 | 48.7332 | 48.7324 | ✓ |
| 1.13 | 27.0835 | 27.0835 | ✓ |
| 1.26 | 11.3336 | 11.3336 | ✓ |

Finite across the whole cavity (72 at the throat r = rs down to ≈ 11.3 at R_OUTER); the only true curvature singularity is at r = 0, behind the throat. The throat r = rs is a coordinate (horizon) singularity.

## AdS₅ / RS reconciliation (Schwarzschild–Tangherlini–AdS₅)

| k (AdS) | Λ₅ = −6k² | max\|R_μν − (−4k²g)\| | Einstein? |
|---:|---:|---:|:---:|
| 0.1 | -0.06 | 1.14e-07 | ✓ |
| 0.3 | -0.54 | 1.63e-07 | ✓ |

`f = 1 − rs²/r² + k²r²` is Einstein with `Λ₅ = −6k²`, interpolating the Tangherlini neck (`k²r² → 0`, PR #116) to the AdS₅/RS asymptote (PR #57). On the cavity the AdS correction is `O(10⁻²)` for `k·rs ≲ 0.1`; the exact `k = κ₅²/Λ₅` stays open (PR #112).

## Verdict

**FIVE_D_TANGHERLINI_BULK_LIFT_RICCI_FLAT_VACUUM_CAVITY_REGULAR.** THE BAM THROAT LIFTS TO A GENUINE D=5 TANGHERLINI VACUUM BULK. PR #116 ran the Tangherlini cavity operator V = f[l(l+2)/r² + 3rs²/r⁴] (f = 1 − (rs/r)²) as a reduced radial object; this probe lifts it to the explicit Schwarzschild–Tangherlini (D=5) metric and verifies the throat is the boundary trace of a real 5D geometry.

THE 5D METRIC & HORIZON. ds² = −f dt² + f⁻¹dr² + r²dΩ₃² with f = 1 − (rs/r)^{D−3} = 1 − (rs/r)² (D=5, power D−3 = 2). The throat is the 5D horizon at r = rs = R_MID (f = 0); its spatial section is the round S³ — exactly BAM's Hopf base S¹ → S³ → S² (the spin/CPT angular base).

RICCI-FLAT VACUUM. R_μν = 0, Λ = 0 (verified numerically across the cavity; the residual is finite-difference noise). The throat parent is a genuine asymptotically-flat vacuum Einstein solution — distinct from the AdS₅ RS bulk (Λ₅ = −6k² < 0).

CAVITY CURVATURE-REGULAR. The Kretschmann scalar K = 72 rs⁴/r⁸ (numeric ≈ analytic to 1e-3) is finite on the whole cavity [R_MID, R_OUTER] — 72 at the throat down to ≈ 11.3 at R_OUTER. The only true curvature singularity is at r = 0, behind the throat; r = rs is a coordinate (horizon) singularity, not a curvature one.

THE CAVITY POTENTIAL DESCENDS FROM D=5. Both coefficients of the PR #116 potential are D=5 reductions of this metric: the centrifugal l(l+2) is the S³ Casimir l(l+D−3) (D−3 = 2), and the curvature term 3rs²/r⁴ = (D−2)/(2r)·f'(r) has coefficient D−2 = 3. So k₅ = D_bulk = 5 (PR #73) is realised as the genuine bulk dimension of the metric — the throat is not a 4D ansatz, it is the boundary of a real D=5 vacuum.

THE HAWKING PERIOD CARRIES 2π. Surface gravity κ = f'(rs)/2 = 1/rs ⟹ T_H = κ/2π = 1/(2π rs): the closure quantum 2π is the thermal/closure period (rs = R_MID = 1 ⟹ T_H = 1/2π).

ADS₅ / RS RECONCILIATION. The Schwarzschild–Tangherlini–AdS₅ metric f = 1 − rs²/r² + k²r² is Einstein with R_μν = −4k²g_μν, Λ₅ = −6k² (verified). It interpolates: near the throat k²r² → 0 and f → the pure Tangherlini neck (PR #116); far away f → k²r², the AdS₅ / RS asymptote (PR #57, the √6 tuning). On the BAM cavity the AdS correction k²r² is O(10⁻²) for k·rs ≲ 0.1, so the pure-Tangherlini cavity of PR #116 is the near-throat limit, good to ~1%.

SCOPE. This establishes the 5D bulk GEOMETRY of the throat and reconciles it with the AdS₅/RS bulk. It does NOT pin the AdS scale k (= κ₅²/Λ₅, the known open bulk ratio, PR #112), nor construct the full global brane-localised black-hole-on-brane solution, nor address bulk backreaction / quantum gravity beyond the classical metric.
