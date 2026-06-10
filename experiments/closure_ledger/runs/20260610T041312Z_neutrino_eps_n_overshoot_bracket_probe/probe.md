# Neutrino log-bounce sensitivity audit and ε_n overshoot bracket (PR #149)

**Run:** 2026-06-10T04:13:12+00:00

The #148-pattern audit applied to the flavor sector's #113 overshoot residual: the log-bounce keystones re-verified, the hypersensitivity inverted into a sub-percent data bracket on the required ε_n profile, the Dirac-attribution endpoints bracketed, and the χ-driven power law excluded under both attributions. The residual is sharply localized — and plausibly belongs to the mixing/anarchy sector (#92). *(QFT on the classical throat, not quantum gravity.)*

- **Mechanism**: m_ν ∝ ε^p, p = 4.8 (L*(ε) slope r_s/2 re-fit)
- **Required profile**: ε₂/ε₁ = 1.352, ε₃/ε₂ = 1.435 (pure-bounce; ~0.1–0.3%)
- **Bracket**: ε₃/ε₂ ∈ [1.323, 1.435] (attribution endpoints)
- **Exclusion**: χ-driven χ₂/χ₃ = 2.487 ⟹ ×14–21 in mass (#113: ×28); no single q fits
- **Consistency**: normal ordering; Σm_ν ≈ 61 meV inside the program window
- **Open**: deriving the gentle profile (mixing/anarchy, #92); #112 anchor; m_D attribution

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the #148-pattern audit of the #113 overshoot residual | **PASS** |
| T2 | `T2_mechanism_keystones` | keystones: L* slope = r_s/2; e^{−cL*} ∝ ε^p identity (p = 4.8) | **PASS** |
| T3 | `T3_hypersensitivity_inverted` | inverted sensitivity: required ε ratios pinned to ~0.1–0.3% | **PASS** |
| T4 | `T4_dirac_attribution_bracket` | attribution bracket [1.32, 1.44]; χ-driven 2.49 far outside | **PASS** |
| T5 | `T5_power_law_exclusion` | no single power law in χ_n; q = 1 overshoots ×21 (#113: ×28) | **PASS** |
| T6 | `T6_required_profile_consistency` | normal ordering; Σm_ν ≈ 61 meV inside the program window | **PASS** |
| T7 | `T7_ledger_and_scope` | ledger: mechanism/bracket/exclusion derived; profile origin residual | **PASS** |
| T8 | `T8_assessment` | NEUTRINO_EPS_N_REQUIRED_PROFILE_PERCENT_BRACKETED_POWER_LAW_EXCLUDED | **PASS** |

## The hypersensitivity, inverted (errors ÷ p)

| required ratio | value | relative uncertainty |
|---|---:|---:|
| ε₂/ε₁ (pure-bounce) | 1.3525 | 0.00275 |
| ε₃/ε₂ (pure-bounce) | 1.4351 | 0.0039 |

## The attribution bracket vs the χ-driven prediction

| quantity | pure-bounce | #113 m_D | χ-driven |
|---|---:|---:|---:|
| required ε₂/ε₁ | 1.352 | 1.186 | 3.134 |
| required ε₃/ε₂ | 1.435 | 1.323 | 2.487 |

The χ-driven ε₃/ε₂ = 2.487 overshoots the pure-bounce requirement by ×14.0 in mass on the 3↔2 pair alone — against a sub-percent data bracket.

## No single power law ε_n ∝ χ_n^{−q}

| attribution | q₁₂ | q₂₃ | q ratio |
|---|---:|---:|---:|
| pure_bounce | 0.264 | 0.396 | 1.5 |
| md_113 | 0.149 | 0.307 | 2.06 |

Best single q = 0.323 still misses a mass ratio by ×1.38; the principled q = 1 overshoots m₃/m₂ by ×20.7 (the #113 headline, re-derived).

## Verdict

**NEUTRINO_EPS_N_REQUIRED_PROFILE_PERCENT_BRACKETED_POWER_LAW_EXCLUDED.** THE ε_n OVERSHOOT RESIDUAL IS NOW A SHARPLY LOCALIZED, DATA-BRACKETED NUMBER: THE LOG-BOUNCE MECHANISM STANDS, THE OSCILLATION DATA PIN THE REQUIRED COMPLIANCE PROFILE AT THE SUB-PERCENT LEVEL, THE REQUIRED SPREAD SITS IN [1.323, 1.435] PER GENERATION STEP AGAINST THE χ-DRIVEN 2.487 (×14–21 IN MASS; ×28 AT THE #113 VINTAGE), AND NO SINGLE POWER LAW IN χ_n FITS. #113 diagnosed the overshoot; this audit brackets it — the #148 pattern applied to the flavor sector.

THE KEYSTONES STAND. The bounce length is the horizon tortoise divergence L*(ε) = (r_s/2)ln(1/ε) (slope re-fit 0.500), and the logarithm converts the exponential mass-action relation into the power law m_ν ∝ ε^p (identity constant over three ε decades; p = 4.8 at the #112 constant).

HYPERSENSITIVITY INVERTED. Forward, the steepness amplified the χ_n spread into the ×28 overshoot (#113). Inverted, δln ε = δln m/p: the oscillation-data errors compress by 1/4.8 into ~0.1–0.3% on the required ε ratios — the residual is not fuzzy, it is a sharply localized number the geometry has not yet produced.

THE ATTRIBUTION BRACKET. The data fix only m_D,n·ε_n^p; bracketing the Dirac attribution between the pure-bounce and #113-implied endpoints puts the required ε₃/ε₂ in [1.323, 1.435]. The principled χ-driven χ₂/χ₃ = 2.487 sits far outside both endpoints — excluded by ×1.7–1.9 in ε, ×14–28 in mass, against a sub-percent bracket.

THE POWER-LAW EXCLUSION. The per-pair implied exponents in ε_n ∝ χ_n^{−q} disagree under BOTH attributions (ratio 1.5–2.1); the best single q = 0.323 still misses a mass ratio by ×1.38 (×1.07 in ε — ~25× the data bracket); q = 1 overshoots ×20.7 at this data vintage (#113 reported ×28). The spread is not a power law in the boundary stress.

CONSISTENCY. With the required profile: normal ordering (derived, #113) and Σm_ν ≈ 61 meV — inside the program's own 59–65 meV window. The residual is consistent, just not derived; it plausibly belongs to the mixing/anarchy sector (#92).

SCOPE. The audit brackets, it does not derive: the gentle profile's origin, the #112 anchor, and the m_D attribution remain open; the input budget is unchanged.
