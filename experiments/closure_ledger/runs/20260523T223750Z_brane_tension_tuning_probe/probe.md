# Brane-tension / bulk-gravity tuning probe

**Run:** 2026-05-23T22:37:50+00:00

Derives the exact Randall–Sundrum-like fine-tuning of the throat brane tension from the Israel junction conditions, sharpening PR #56's parametric σ ∝ √|Λ₅|/κ₅ to an exact relation with the dimensionless factor √6.

- **Derived**: `λ_crit = √(6|Λ₅|)/κ₅² (RS fine-tuning)`
- **Tuning factor**: sqrt(6) ≈ 2.449
- **Meaning**: flat / static-throat condition (Λ₄ = 0)
- **B4 caveat**: one condition among (λ,Λ₅,κ₅); √6 derived, bulk scale is the anchor

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_israel_junction_pure_tension` | K_μν = −(κ₅²λ/6) h_μν (factor 1/6) | **PASS** |
| T2 | `T2_bulk_ads5_relation` | Λ₅ = −6k² ⟹ k = √(|Λ₅|/6) | **PASS** |
| T3 | `T3_fine_tuning_factor_sqrt6` | λ_crit κ₅²/√|Λ₅| = √6 ≈ 2.449 | **PASS** |
| T4 | `T4_flat_static_throat_condition` | Λ₄ ∝ λ²−λ_crit²; zero at tuning (flat throat) | **PASS** |
| T5 | `T5_cohesive_connection` | B=4πσ inherits √6 + bulk scale | **PASS** |
| T6 | `T6_b4_accounting` | one condition; √6 invariant; scale = anchor | **PASS** |
| T7 | `T7_warp_hierarchy` | warp w = e^{−kΔR} (RS hierarchy) | **PASS** |
| T8 | `T8_assessment` | RS tuning derived; value still the anchor | **PASS** |

## T1: Israel junction (pure tension) → factor 1/6

- K_μν coefficient derived: -0.500000
- expected −κ₅²λ/6: -0.500000
- 1/6 factor verified: True

## T2: Bulk AdS₅ relation

| k | Λ₅ = −6k² | k recovered | deviation |
|---:|---:|---:|---:|
| 1.00 | -6.0000 | 1.000000 | 0.0e+00 |
| 2.00 | -24.0000 | 2.000000 | 0.0e+00 |
| 3.50 | -73.5000 | 3.500000 | 0.0e+00 |

## T3: Fine-tuning factor √6

| k | κ₅² | Λ₅ | λ_crit=6k/κ₅² | √(6|Λ₅|)/κ₅² | factor |
|---:|---:|---:|---:|---:|---:|
| 1.00 | 1.00 | -6.000 | 6.0000 | 6.0000 | 2.449490 |
| 2.00 | 0.50 | -24.000 | 24.0000 | 24.0000 | 2.449490 |
| 3.50 | 2.00 | -73.500 | 10.5000 | 10.5000 | 2.449490 |

Dimensionless tuning factor = √6 = 2.449490.

## T4: Flat / static-throat condition

| λ/λ_crit | λ | Λ₄ ∝ λ²−λ_crit² | sign |
|---:|---:|---:|---|
| 0.50 | 6.0000 | -108.0000 | AdS(-) |
| 0.90 | 10.8000 | -27.3600 | AdS(-) |
| 1.00 | 12.0000 | +0.0000 | flat |
| 1.10 | 13.2000 | +30.2400 | dS(+) |
| 2.00 | 24.0000 | +432.0000 | dS(+) |

Tuned flat: True; over-tension dS: True; under-tension AdS: True.

## T5: Connection to the cohesive term (PR #56)

- λ_crit = 12.0000; σ (4D-effective image) = 12.0000
- B = 4πσ = 150.7964; inherits √6: True
- the tuning (flat brane) is the static-throat condition (PR #55): True

## T6: B4 accounting

| k scale | k | λ_crit | dimensionless factor |
|---:|---:|---:|---:|
| 1.0 | 2.00 | 12.0000 | 2.449490 |
| 2.0 | 4.00 | 24.0000 | 2.449490 |
| 0.5 | 1.00 | 6.0000 | 2.449490 |
| 10.0 | 20.00 | 120.0000 | 2.449490 |

One tuning condition among (λ, Λ₅, κ₅) → one dimensionful combination remains (the anchor); √6 is invariant under rescaling the bulk scale.

## T7: Warp hierarchy

Bulk depth ΔR (geometric) = 0.52

| k·L | warp w = e^{−kL} | hierarchy (orders) |
|---:|---:|---:|
| 1.0 | 3.679e-01 | 0.4 |
| 5.0 | 6.738e-03 | 2.2 |
| 10.0 | 4.540e-05 | 4.3 |
| 35.0 | 6.305e-16 | 15.2 |

## T8: Assessment

- **Tuning factor**: sqrt(6) ≈ 2.449
- **Tuned tension**: λ_crit = √(6|Λ₅|)/κ₅²
- **Meaning**: flat / static-throat condition (Λ₄ = 0)
- **Anchor**: bulk gravitational scale k = √|Λ₅/6| (one input)

## Verdict

**BRANE_TUNING_DERIVED.** BRANE TUNING DERIVED. The Randall–Sundrum-like fine-tuning of the throat brane tension is derived from the junction conditions, sharpening PR #56's parametric σ ∝ √|Λ₅|/κ₅ to an exact relation with a derived dimensionless factor.

DERIVATION. (1) The Z₂ Israel junction for a pure-tension 4D worldvolume brane (S_μν=−λh_μν) trace-reverses to K_μν = −(κ₅²λ/6) h_μν — the 1/6 factor. (2) The bulk AdS₅ vacuum equation gives Λ₅ = −6k², so k = √(|Λ₅|/6). (3) Staticity (a flat brane, K_μν = k h_μν) matches the junction at k = κ₅²λ/6, giving the fine-tuned tension λ_crit = 6k/κ₅² = √(6|Λ₅|)/κ₅². The dimensionless tuning factor is λ_crit κ₅²/√|Λ₅| = √6 ≈ 2.449.

MEANING. Detuning induces a 4D cosmological constant Λ₄ ∝ λ²−λ_crit², so Λ₄=0 ⟺ λ=λ_crit: the fine-tuning is the flat, static-throat condition (a stable particle with no induced 4D vacuum energy). Over-tension → de Sitter throat, under-tension → anti-de Sitter; the critically-tuned, static configuration is the cohesive equilibrium of PR #55.

CONNECTION + B4. The cohesive B = 4πσ (PR #56) is the 4D-effective image of λ_crit; it inherits the √6 factor and the bulk-gravity scale √|Λ₅|/κ₅². The fine-tuning is one condition among (λ, Λ₅, κ₅), so a net one dimensionful combination remains — the single anchor the B4 scale-modulus theorem requires. The dimensionless content (√6, the flatness condition) is derived; the absolute scale (k = √|Λ₅/6|, the bulk gravitational scale) is the one external input. The AdS₅ warp over the bulk depth ΔR gives an RS exponential hierarchy w = e^{−kΔR}.

## What this leaves open

- **The bulk-gravity scale (the anchor).** k = √|Λ₅/6| (equivalently Λ₅, κ₅) is the one dimensionful input; a genuine derivation needs a second fixed scale.
- **The BAM-throat junction from S_BAM.** The derivation uses the canonical RS Z₂ brane; matching to the exact BAM throat (Tangherlini interior + closure-quantum surface) is the follow-on.
- **Pair-production threshold.** 2 m_e c² at the lowest stable R* as a dynamical nucleation calculation.
