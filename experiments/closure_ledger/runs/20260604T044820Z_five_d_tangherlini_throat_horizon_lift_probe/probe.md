# Horizon-regular coordinate lift for the 5D Tangherlini throat (PR #128)

**Run:** 2026-06-04T04:48:20+00:00

PR #127 showed the throat r = rs is a *coordinate* (horizon) singularity, not a curvature one. This probe builds the horizon-regular coordinates (Eddington–Finkelstein, Kruskal–Szekeres) that remove it, make the throat crossing smooth, and exhibit the antipodal bifurcation structure that is the geometric home of BAM's throat ↔ antithroat C-swap. Curvature is computed by a self-contained numerical GR routine.

- **Coordinate singularity**: removable (g_rr = 1/f → ∞ but K = 72 rs⁴/r⁸ finite)
- **Eddington–Finkelstein**: ds² = −f dv² + 2 dv dr + r²dΩ₃², det g = −r⁶ sin⁴χ sin²θ regular
- **Proper distance**: finite √(2 rs Δr) = ε healing length (#112); tortoise r* → −∞
- **Surface gravity**: κ = 1/rs (κ·rs = 1); Kruskal F(rs) = 4 e⁻²; T_H = 1/(2π rs)
- **Maximal extension**: bifurcate Killing horizon U = V = 0; four regions
- **Antipodal identification**: (U,V) → (−U,−V) = throat ↔ antithroat C-swap (#63, #58)
- **Open**: nucleation rate (#58/#88); exact AdS scale k (#112); global brane solution (#127)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | construct horizon-regular coordinates (PR #127 coord. singularity) | **PASS** |
| T2 | `T2_coordinate_singularity_removable` | g_rr = 1/f → ∞ but K = 72 rs⁴/r⁸ finite ⟹ coordinate artifact | **PASS** |
| T3 | `T3_eddington_finkelstein_regular` | Eddington–Finkelstein regular: det g = −r⁶ sin⁴χ sin²θ, K finite | **PASS** |
| T4 | `T4_tortoise_vs_proper_distance` | tortoise r* → −∞ but proper √(2 rs Δr) finite (= ε healing length) | **PASS** |
| T5 | `T5_surface_gravity_and_kruskal_factor` | κ = 1/rs (κ·rs = 1); Kruskal F(rs) = 4 e⁻²; T_H = 1/(2π rs) | **PASS** |
| T6 | `T6_maximal_kruskal_extension` | maximal extension: bifurcate Killing horizon U = V = 0, 4 regions | **PASS** |
| T7 | `T7_antipodal_identification_is_C_swap` | antipodal (U,V)→(−U,−V) = throat ↔ antithroat C-swap (#63, #58) | **PASS** |
| T8 | `T8_assessment` | FIVE_D_TANGHERLINI_THROAT_HORIZON_REGULAR_ANTIPODAL_BIFURCATION | **PASS** |

## Eddington–Finkelstein metric is regular at the throat

| r | g_vv = −f | g_vr | det g | K (EF) | K analytic |
|---:|---:|---:|---:|---:|---:|
| 1.0 | -0.0 | 1.0 | -0.29904 | 72.0 | 72.0 |
| 1.13 | -0.21685 | 1.0 | -0.62259 | 27.0835 | 27.0835 |
| 1.26 | -0.37012 | 1.0 | -1.19661 | 11.3336 | 11.3336 |

At the throat g_vv = 0 but g_vr = 1, so the metric is nondegenerate (det g = −r⁶ sin⁴χ sin²θ ≠ 0) and the Kretschmann scalar in EF coordinates is 72 rs⁴/r⁸ — the same regular geometry, now with a finite metric.

## Tortoise (infinite) vs proper (finite) distance to the throat

| Δr = r − rs | tortoise r* | proper ∫dr/√f | √(2 rs Δr) |
|---:|---:|---:|---:|
| 0.1 | -0.4223 | 0.4615 | 0.44721 |
| 0.01 | -1.6417 | 0.14148 | 0.14142 |
| 0.001 | -2.7997 | 0.04429 | 0.04472 |

The throat is infinitely far in the tortoise/optical coordinate (r* → −∞) but a finite proper distance √(2 rs Δr) away — exactly the ε healing length √(2 rs ε) (PR #112).

## Verdict

**FIVE_D_TANGHERLINI_THROAT_HORIZON_REGULAR_ANTIPODAL_BIFURCATION.** THE 5D TANGHERLINI THROAT LIFTS TO HORIZON-REGULAR COORDINATES, AND ITS MAXIMAL EXTENSION IS THE ANTIPODAL STAGE OF BAM. PR #127 showed the throat r = rs is a coordinate (horizon) singularity, not a curvature one (K = 72 rs⁴/r⁸ finite); this probe builds the regular charts that remove it and exhibits the antipodal bifurcation structure.

THE COORDINATE SINGULARITY IS REMOVABLE. In Schwarzschild coordinates g_rr = 1/f → ∞ as r → rs, but K = 72 rs⁴/r⁸ is finite there: the divergence is a coordinate artifact.

EDDINGTON–FINKELSTEIN IS REGULAR. With the tortoise r* = r + (rs/2)ln|(r−rs)/(r+rs)| and v = t + r*, ds² = −f dv² + 2 dv dr + r²dΩ₃². At the throat g_vv = 0 but g_vr = 1, so det g = −r⁶ sin⁴χ sin²θ is finite and nonzero, and the Kretschmann scalar computed in EF coordinates is still 72 rs⁴/r⁸ — the same regular geometry, now with a nondegenerate metric.

THE THROAT IS FINITE PROPER DISTANCE AWAY. The tortoise distance to the throat is infinite (r* → −∞), but the proper radial distance ∫dr/√f ≈ √(2 rs(r−rs)) is finite — exactly the ε healing length √(2 rs ε) (PR #112).

SURFACE GRAVITY & KRUSKAL FACTOR. κ = f'(rs)/2 = 1/rs, so κ·rs = 1. The Kruskal conformal factor F = −f·e^{−2κr*} = (r+rs)²/r²·e^{−2r/rs} is finite/nonzero at the throat (F(rs) = 4 e⁻²) precisely because κ·rs = 1 makes the (r−rs)^{−κrs} factor cancel the simple zero of f. The Hawking temperature T_H = κ/2π = 1/(2π rs) carries the closure quantum (PR #127).

THE MAXIMAL EXTENSION. In Kruskal coordinates UV = −(1/κ²)e^{2κr*} → 0 at the throat: the bifurcate Killing horizon U = V = 0. The extension has four regions (exterior I, interior II, antipodal exterior III, white hole IV).

THE ANTIPODAL IDENTIFICATION = BAM's C-SWAP. The antipodal map (U, V, Ω) → (−U, −V, Ω_antipodal) is an isometry preserving UV (hence r) and mapping region I ↔ region III. This is the geometric realisation of BAM's throat ↔ antithroat identification — the C = inner/outer swap (PR #63) with c₁ → −c₁ (PR #58). The maximally-extended 5D Tangherlini horizon with its antipodal gluing is the geometric stage of "Bulk Antipodal Mechanics" itself.

SCOPE. This is the classical, maximally-extended vacuum geometry — the regular charts, the finite proper distance, the surface gravity, the antipodal bifurcation. It does NOT compute the dynamical throat ↔ antithroat nucleation amplitude (the bounce rate, PR #58/#88) — the lift provides that process's kinematic stage, not its rate. The exact AdS scale k (= κ₅²/Λ₅, PR #112) and the global brane-localised solution (PR #127) remain open.
