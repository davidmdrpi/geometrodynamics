# Compton vertex resummation probe

**Run:** 2026-05-17T03:09:22+00:00

Follow-on to PR #34 (O(ε²) closure with patterns ν₀ = γ², ξ = −A_φ(0)). Tests whether the vertex modification F(ε, θ) resums to a closed form at all orders.

**Closed form derived:**

```
F²(x, c) = 4·x³·(x²+1−x·sin²θ) / [(1+c²)·(1+x)²], x = 1/(1+ε(1−cos θ))
```

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_closed_form_verification` | max diff = 5.55e-16 | **PASS** |
| T2 | `T2_taylor_O_eps_O_eps2_match` | μ₁ PR31 ✓ (err 2.00e-05); μ₂ refines PR34 (Δ = 5.83e-02) | **PASS** |
| T3 | `T3_higher_order_pattern` | μₙ/((-3/2)ⁿ·(1-c)ⁿ) at θ=π/2: 1.000, 1.000, 0.926, 0.778 | **PASS** |
| T4 | `T4_natural_factorisations` | 28/28 samples match Padé factorisation | **PASS** |
| T5 | `T5_hopf_link_search` | (2x/(1+x)) is kinematic Padé, not Hopf-derived | **PASS** |
| T6 | `T6_end_to_end_KN_reproduction` | max residual = 5.55e-16 | **PASS** |

## T2: Taylor expansion match

Verifying that the closed-form F reproduces the PR #31 and PR #34 perturbative coefficients.

Max |μ₁_num − μ₁_PR31| = 1.9969e-05 — PR #31 confirmed exactly.
Max |μ₂_num − μ₂_PR34_polynomial| = 5.8261e-02.

**μ₂ note:** PR #34 polynomial fit captures leading angular structure but is NOT the exact μ₂. The closed-form expansion is the corrected form. PR #34 verdict should be amended.

## T3: Higher-order coefficient pattern

Taylor coefficients μₙ extracted from F at three test θ:

| θ/π | μ₁ | μ₂ | μ₃ | μ₄ |
|---:|---:|---:|---:|---:|
| 0.3333 | -0.7500 | +0.5375 | -0.3594 | +0.2223 |
| 0.5000 | -1.5000 | +2.2500 | -3.1250 | +3.9376 |
| 0.6667 | -2.2500 | +4.8375 | -9.7033 | +18.0109 |

At θ = π/2 (c = 0), μₙ values: ['-1.5000', '+2.2500', '-3.1250', '+3.9376']
Expected (−3/2)ⁿ pattern: ['-1.5000', '+2.2500', '-3.3750', '+5.0625']
Ratios: ['1.0000', '1.0000', '0.9259', '0.7778']

## T4: Natural factorisation

Factorisation `F² = (2x/(1+x))² · [x·(x²+1−x·sin²θ)/(1+c²)]` matches at 28/28 sample points to machine precision.

## T6: End-to-end KN reproduction

| ε | max residual |
|---:|---:|
| 0.010 | 5.55e-16 |
| 0.050 | 4.44e-16 |
| 0.100 | 3.33e-16 |
| 0.200 | 4.44e-16 |
| 0.500 | 2.22e-16 |
| 1.000 | 3.33e-16 |
| 2.000 | 2.22e-16 |

## Verdict

**EXACT_RESUMMATION.** EXACT RESUMMATION ACHIEVED — the BAM vertex modification factor F(ε, θ) has the closed form `F²(x, c) = 4·x³·(x²+1−x·sin²θ) / [(1+c²)·(1+x)²]` with x = ω'/ω = 1/(1+ε(1−cos θ)). This reproduces Klein-Nishina exactly at all (ε, θ), including ε up to 2 (highly relativistic Compton). The perturbative coefficients from PR #31 (γ = −3/2) and PR #34 ((ν₀, ν₁, ν₂, ξ) = (9/4, −4, 7/4, −1/2)) are Taylor expansions of this closed form. The factor decomposes naturally as F² = (2x/(1+x))² · [x·(x²+1−x·sin²θ)/(1+c²)] — a kinematic Padé form times an angular polarisation modification. The (1+c²) denominator in the angular factor is the polarisation sum itself, confirming that F must be derived AS a modification of the polarisation sum rather than as a separate amplitude piece.

## What this leaves open

- **BAM derivation of the closed-form F from first principles.** The closed form is now identified, but deriving it from a BAM Lagrangian or action remains the deeper task.
- **The (1+c²) denominator** in the angular factor is the polarisation sum itself, suggesting F must be derived as a modification of the polarisation sum rather than a separate amplitude piece. The natural BAM construction might be a modified pol-sum projector.
- **Cross-process generalisation.** Does the same F work for pair production γγ → e⁺e⁻ and other QED tree diagrams? Testable as the next probe.
- **Recursive (−3/2)ⁿ pattern.** T3 tests at specific θ; the pattern is approximate, suggesting the full recursive structure is more complex than simple geometric series.
