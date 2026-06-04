# Throat–antithroat dynamical nucleation rate on the horizon-regular 5D background (PR #132)

**Run:** 2026-06-04T23:19:46+00:00

Closes the #131 lead open item: the throat ↔ antithroat nucleation rate, placed on the horizon-regular 5D background and connected to the Majorana bounce arc (#87–#90). The nucleation is the odd antipodal instanton on a smooth Euclidean cigar; the geometry supplies the smoothness condition, the ln(1/ε) origin of the action, and the prefactor.

- **Smooth cigar**: β = 2π/κ = 2π rs (deficit 0); T_nuc = 1/(2π rs) = T_H (closure quantum)
- **Bounce**: odd antipodal instanton, region I↔III (#128), c₁→−c₁ (#63), ΔL=2 (#58)
- **Action**: S ∝ L*(ε) = (rs/2) ln(1/ε) + const (horizon tortoise divergence)
- **Rate**: Γ ~ det^{−1/2} e^{−S}; t ∈ [2π, k₅√(2π)] (#89), ε ~ R_c³ ⟹ S ≈ 15–18, m_ν ~ meV
- **Prefactor**: the #116 Tangherlini determinant 1.574370
- **Open**: exact ε; absolute scale κ₅²/Λ₅; precise S / m_ν (#88–#90, #112)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | throat↔antithroat nucleation rate on the regular 5D background | **PASS** |
| T2 | `T2_smooth_euclidean_cigar` | smooth Euclidean cigar: deficit 0 at β = 2π/κ = 2π rs (T_nuc = T_H) | **PASS** |
| T3 | `T3_antipodal_odd_instanton` | bounce = antipodal odd instanton: I↔III (#128), c₁→−c₁ (#63) | **PASS** |
| T4 | `T4_log_epsilon_horizon_tortoise_divergence` | action ln(1/ε) = horizon tortoise divergence (slope rs/2) | **PASS** |
| T5 | `T5_reduced_action_and_rate` | S ≈ 15–18, m_ν = m_D e^{−S} ~ meV (t window #89, ε ~ R_c³) | **PASS** |
| T6 | `T6_prefactor_is_tangherlini_determinant` | prefactor = the #116 Tangherlini determinant 1.574370 | **PASS** |
| T7 | `T7_scope` | scope: new (smoothness, ln(1/ε), prefactor) vs inherited residuals | **PASS** |
| T8 | `T8_assessment` | THROAT_ANTITHROAT_NUCLEATION_REGULAR_CIGAR_LOG_EPSILON_BOUNCE | **PASS** |

## The smooth Euclidean cigar (deficit 2π − κβ)

| β | deficit 2π − κβ | smooth? |
|---|---:|:---:|
| 6.28319 (β = 2π/κ) | 0.0 | ✓ |
| 5.02655 (0.8×) | 1.256637 | ✗ |
| 7.53982 (1.2×) | -1.256637 | ✗ |

Smooth (no conical defect) only at `β = 2π/κ = 2π rs`; `T_nuc = 1/β = 0.159155 = T_H` — the closure quantum 2π.

## The bounce action is the horizon tortoise divergence

| ε | L*(ε) |
|---:|---:|
| 0.01 | 1.8204 |
| 0.001 | 2.9785 |
| 0.0001 | 4.1304 |
| 1e-05 | 5.2818 |
| 1e-06 | 6.4331 |

`L*(ε) = (rs/2) ln(1/ε) + const`, asymptotic slope `0.5 = rs/2`. The exact-horizon limit `ε → 0` costs infinite tortoise length ⟹ `S → ∞`, `m_ν → 0` (rigid throat ⟹ massless ν, #88).

## Verdict

**THROAT_ANTITHROAT_NUCLEATION_REGULAR_CIGAR_LOG_EPSILON_BOUNCE.** THE THROAT ↔ ANTITHROAT NUCLEATION ON THE HORIZON-REGULAR BACKGROUND IS THE ODD ANTIPODAL INSTANTON ON A SMOOTH EUCLIDEAN CIGAR. The geometric throat arc (#127–#131) built the regular 5D background and identified the antipodal throat ↔ antithroat map but left the nucleation rate open; the Majorana bounce arc (#87–#90) had the bounce action on the EM/tortoise picture. This probe puts the nucleation on the regular background and supplies what the geometry contributes.

THE SMOOTH EUCLIDEAN CIGAR. Wick-rotating, the near-horizon Euclidean metric in the proper radius ρ = √(2 rs(r−rs)) is ds²_E ≈ dρ² + κ²ρ² dτ² with κ = f'(rs)/2 = 1/rs. This is the flat plane in polar coordinates (ρ, κτ) and is smooth — deficit 2π − κβ = 0 — iff the imaginary-time period is β = 2π/κ = 2π rs (Gibbons–Hawking). So the Euclidean throat closes off smoothly, the nucleation temperature is T_nuc = 1/β = 1/(2π rs) = T_H, and the period is the closure quantum 2π (#127).

THE BOUNCE IS THE ANTIPODAL ODD INSTANTON. The throat ↔ antithroat transition (the ΔL = 2 Majorana / pair-production channel, #58) is the region I ↔ III crossing of the maximal Kruskal extension (#128), mediated by the odd (c₁ → −c₁, the C-swap #63) instanton, with rate Γ ~ det^{−1/2} e^{−S}.

THE ACTION IS THE HORIZON TORTOISE LENGTH. The bounce tortoise length to the ε healing length is L*(ε) = (rs/2) ln(1/ε) + const (slope rs/2 = 0.5, verified to 4 digits). The reduced action S = (tension)·√(2μ E_c)·L*(ε) therefore grows as ln(1/ε): the exact-horizon limit ε → 0 costs infinite tortoise length ⟹ S → ∞, Γ → 0, m_ν → 0 — the "rigid throat ⟹ massless neutrino" of #88, now read off directly as the horizon tortoise divergence, regulated by the finite ε healing length (#112).

THE RATE. With the ΔL = 2 tension window t ∈ [2π, k₅√(2π)] ≈ [6.28, 12.53] (#89) and ε ~ R_c³ (#112), the #88–#90 chain gives S ≈ 15–18 and m_ν = m_D e^{−S} ~ few meV — the observed scale, retrodicted to order of magnitude.

THE PREFACTOR CLOSES THE ARC. The one-loop nucleation prefactor is the Tangherlini fluctuation determinant of PR #116, det(H)/det(H_free) = 1.574370. The geometric arc closes on itself: #116 is the bounce prefactor, #127/#128 the regular stage, #58/#87–#90 the bounce.

SCOPE. NEW here: the Euclidean smoothness β = 2π rs (the closure quantum), the ln(1/ε) as the horizon tortoise divergence (slope rs/2), the prefactor = the #116 determinant, the antipodal-instanton structure. INHERITED / open: the exact ε value, the absolute gravitational scale κ₅²/Λ₅, and hence the precise S and m_ν (#88–#90, #112). The rate is order-of-magnitude (meV), not pinned.
