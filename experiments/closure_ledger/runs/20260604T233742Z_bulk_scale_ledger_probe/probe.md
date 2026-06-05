# Bulk scale ledger for κ₅²/Λ₅ and ΔR normalization (PR #133)

**Run:** 2026-06-04T23:37:42+00:00

Consolidating ledger for the absolute bulk scale — the recurring κ₅²/Λ₅ residual (#57/#112/#127/#132). It counts the bulk dimensionful content, separates the scale modulus ΔR (the unit) from the residual, fixes √6 (the RS tuning), and bounds the one open AdS ratio.

- **Unit**: ΔR = R_OUTER − R_INNER = 0.52 R_MID (scale modulus, #52/#53)
- **Fixed tuning**: √6 = λ_crit κ₅²/√|Λ₅| ≈ 2.449 (RS flatness, #57)
- **Anchor**: G = κ₅²/ΔR³ (gravity strength, the one dimensionful anchor, #105/#106)
- **Open (bounded)**: k·R_MID = R_MID/L_AdS ≲ 0.1 (the κ₅²/Λ₅ residual, #112/#127)
- **Residual count**: one bounded dimensionless number

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | bulk scale ledger for κ₅²/Λ₅ and ΔR (#57/#112/#127/#132) | **PASS** |
| T2 | `T2_dimensionful_content` | dimensions: κ₅²[L³], Λ₅[L⁻²], k[L⁻¹], λ_crit[L⁻⁴], R_MID/ΔR[L] | **PASS** |
| T3 | `T3_delta_r_is_scale_modulus` | ΔR = scale modulus (the unit, #52/#53), not a residual | **PASS** |
| T4 | `T4_sqrt6_fixed_tuning` | √6 = λ_crit κ₅²/√|Λ₅| — the one fixed RS tuning (#57) | **PASS** |
| T5 | `T5_ads_scale_ratio_bounded` | open k·R_MID bounded ≲ 0.1 by the cavity correction (k r)² (#127) | **PASS** |
| T6 | `T6_consolidated_ledger` | ledger: {κ₅²,Λ₅} → {G} + {√6} + {k·rs bounded} + {ΔR unit} | **PASS** |
| T7 | `T7_scope` | scope: bounds the residual, does not pin it (= the #112 residual) | **PASS** |
| T8 | `T8_assessment` | BULK_SCALE_LEDGER_ONE_BOUNDED_ADS_RATIO_DELTAR_IS_SCALE_MODULUS | **PASS** |

## The open AdS ratio is bounded (cavity correction (k r)², #127)

| k·R_MID | cavity correction | L_AdS/R_MID |
|---:|---:|---:|
| 0.05 | 0.4% | 20.0 |
| 0.1 | 1.59% | 10.0 |
| 0.2 | 6.35% | 5.0 |

`k·R_MID ≲ 0.1` keeps the cavity correction `≲ 1.59%`, so `R_MID ≲ L_AdS/10`: the throat sits deep in the near-flat AdS region — why the pure-Tangherlini cavity (#116/#127) is a good approximation.

## The consolidated ledger

| category | content |
|---|---|
| **unit** | ΔR (scale modulus, #52/#53) — the length unit |
| **fixed tuning** | √6 = λ_crit κ₅²/√|Λ₅| (#57) |
| **anchor** | G = κ₅²/ΔR³ (gravity strength, the one dimensionful anchor, #105/#106) |
| **open (bounded)** | k·R_MID = R_MID/L_AdS ≲ 0.1 (the κ₅²/Λ₅ residual, #112/#127) |

The recurring `κ₅²/Λ₅` residual is **one bounded dimensionless number** (`k·R_MID ≲ 0.1`), not a multi-parameter freedom — isolated and bounded, though not pinned.

## Verdict

**BULK_SCALE_LEDGER_ONE_BOUNDED_ADS_RATIO_DELTAR_IS_SCALE_MODULUS.** THE BULK SCALE LEDGER REDUCES κ₅²/Λ₅ TO ONE BOUNDED DIMENSIONLESS NUMBER, WITH ΔR THE UNIT. The absolute bulk scale surfaced as an open residual at every step (#57/#112/#127/#132); this ledger counts the content and separates the scale modulus from the residual.

THE DIMENSIONFUL CONTENT. In D=5 (ℏ=c=1): κ₅²[L³] (the 5D Newton constant G₅), Λ₅[L⁻²] ⟺ k = √(|Λ₅|/6)[L⁻¹] (the AdS₅ inverse radius, L_AdS = 1/k), λ_crit = 6k/κ₅²[L⁻⁴] (the 4D brane tension), and the geometric lengths R_MID, ΔR[L].

ΔR IS THE SCALE MODULUS — THE UNIT, NOT A RESIDUAL. The B4 scale-modulus theorem (#52) proved BAM cannot derive an absolute unit from scale-free topology: exactly one external dimensionful anchor is required. ΔR = R_OUTER − R_INNER = 0.52 R_MID is that anchor — a proper, cosmologically-invariant length (#53) — and it sets the unit (R_MID = 1). The geometry ratios ΔR/R_MID = 0.52, R_OUTER/R_MID = 1.26 are fixed. So ΔR is units, not a residual.

√6 IS THE ONE FIXED DIMENSIONLESS TUNING. The Randall–Sundrum flatness condition is λ_crit κ₅²/√|Λ₅| = √6 ≈ 2.449 (the 6 is the AdS₅ curvature coefficient Λ₅ = −6k²) — one condition among (λ, Λ₅, κ₅), derived (#57).

THE OPEN BULK NUMBER IS BOUNDED. Once the unit (ΔR) and the tuning (√6) are fixed, the only remaining dimensionless bulk freedom is the AdS scale in throat units, k·R_MID = R_MID/L_AdS (= κ₅²/Λ₅ in throat units). It is not pinned, but it is BOUNDED: the cavity correction to the pure-Tangherlini background is (k r)² (#127), so k·R_MID ≲ 0.1 keeps it ≲ 1.6% on the cavity. Hence R_MID ≲ L_AdS/10 — the throat sits deep in the near-flat region of the AdS bulk, which is exactly why the pure-Tangherlini cavity (#116/#127) is a good approximation.

THE CONSOLIDATED LEDGER. {κ₅², Λ₅} ⟶ { G (gravity strength κ₅²/ΔR³ = the dimensionful anchor, #105/#106) } + { √6 (RS tuning, FIXED) } + { k·R_MID (AdS scale, OPEN but bounded ≲ 0.1) }, with ΔR the unit. So the recurring "κ₅²/Λ₅ mystery" is ONE bounded dimensionless number, not a multi-parameter freedom.

SCOPE. A consolidating/accounting ledger: it bounds and isolates the residual, it does not PIN k·R_MID (still open, = the #112 residual), nor add any new free parameter. Ties to the input budget (#104): G is the one dimensionful anchor; this audits its bulk origin.
