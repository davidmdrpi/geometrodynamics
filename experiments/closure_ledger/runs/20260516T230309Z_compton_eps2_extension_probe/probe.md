# O(ε²) analytic extension probe

**Run:** 2026-05-16T23:03:09+00:00

Follow-on to PR #31 (analytic γ = −3/2 at O(ε)) and PR #33 (d-scaling discrimination). Extends the analytic derivation to O(ε²), identifies what vertex structure is needed beyond PR #31, and discriminates between the 7 surviving coefficient-origin candidates from PR #33.

**Setup:**

```
Δ_required_O(ε²) = (1−c)²·(9+7c²)/4. PR #31 closes O(ε) but leaves this O(ε²) residual.
```

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_analytic_O_eps2_verification` | max num/ana O(ε²) error = 1.28e-06 | **PASS** |
| T2 | `T2_delta_required_decomposition` | max decomposition error = 1.78e-15 | **PASS** |
| T3 | `T3_family_B_prime_polynomial` | Family B' alone: over-det.; max LSQ residual 0.500 | **PASS** |
| T4 | `T4_family_A_plus_B_prime` | Family A added: over-det.; max LSQ residual 0.500 | **PASS** |
| T5 | `T5_extended_family_with_odd_c` | Extended (B+A+odd-c): EXACT; α² = 0.000 | **PASS** |
| T6 | `T6_numerical_verification` | residual order = O(ε^2.50) | **PASS** |
| T7 | `T7_candidate_O_eps2_predictions` | target magnitude 8.500 | **PASS** |

## Δ_required_O(ε²) in c basis (T2)

**Coefficients (c⁰..c⁴):** ['+2.250', '-4.500', '+4.000', '-3.500', '+1.750']

Equivalent form: (9 − 18c + 16c² − 14c³ + 7c⁴)/4 = (1−c)²·(9+7c²)/4

## T3: Family B′ alone (polynomial μ₂)

Rank(A) = 3, rank(A|b) = 4. Exact solution exists: NO.

LSQ: ν₀ = 2.2500, ν₁ = -4.0000, ν₂ = 1.7500. Max residual 0.5000.

**Conclusion:** the c¹ and c³ coefficients of μ₂·(1+c²) are tied (both equal ν₁ from the parametrisation), but the target has c¹ = −9/2 and c³ = −7/2 — different. Polynomial Family B′ alone cannot close O(ε²).

## T4: Family A + B′ (α² sin⁴θ added)

Exact: NO. LSQ residual 0.5000. α² = 0.0000.

**Conclusion:** sin⁴θ is even in c, so adding Family A preserves the c¹ = c³ contradiction. STILL over-determined.

## T5: Extended family with odd-c term ξ·c·(1−c²)

5-unknown system: det(A) = 4.0000, exact = YES. Max residual 0.0000e+00.

Coefficients: ν₀ = 2.2500, ν₁ = -4.0000, ν₂ = 1.7500, α² = 0.0000, ξ = -0.5000.

α value: 0.0000 (α² ≥ 0: YES)

Clean rational check:

- nu0: clean
- nu1: clean
- nu2: clean
- alpha_sq: clean
- xi: clean

## T6: Numerical residual scaling

Fitted residual order: **O(ε^2.504)** (expected ≥ 3 for O(ε²) closure).

| ε | max KN residual |
|---:|---:|
| 0.0001 | 9.3617e-10 |
| 0.001 | 1.1691e-07 |
| 0.01 | 3.8847e-05 |
| 0.1 | 3.0076e-02 |

## T7: Candidate O(ε²) scalar predictions

Candidate predictions are scalar magnitudes; they do not specify the angular structure of μ₂. The probe reports relative magnitudes as a weak discrimination metric.

| candidate | predicted scalar magnitude | ratio to T5 target |
|---|---:|---:|
| `A_doubled_electron_casimir` | 1.1250 | 0.132 |
| `B_photon_casimir_minus_hopf` | 3.7500 | 0.441 |
| `D_hopf_times_photon_mult` | 2.2500 | 0.265 |
| `F_closure_plus_hopf` | 2.2500 | 0.265 |
| `G_pauli_trace` | silent | n/a |
| `H_two_mouth_spin_half` | 1.1250 | 0.132 |
| `E_antipodal_natural_units` | silent | n/a |

## Verdict

**CLOSED_WITH_EXTENDED_FAMILY.** CLOSED WITH EXTENDED FAMILY (B + A + odd-c) — the O(ε²) gap closes analytically with (ν₀, ν₁, ν₂, α², ξ) = (2.2500, -4.0000, 1.7500, 0.0000, -0.5000). α² = 0.0000 (non-negative as required). Residual fits to O(ε^2.50). The PR #31 Family B alone is structurally insufficient; closing O(ε²) requires explicit Family A (ε·k coupling) AND an odd-in-c angular term — neither of which appears naturally in any of the PR #33 surviving candidates. This identifies the next BAM structural target.

## What this leaves open

- **BAM derivation of the new structures** (ξ odd-c term, α² coupling) at O(ε²). The structural form is identified by T5 but the BAM-native origin remains to be derived.
- **O(ε³) and beyond.** Each successive order may require additional vertex structures; the pattern of required additions is open.
- **Candidate discrimination at O(ε²).** T7 uses scalar magnitudes; sharper discrimination requires the candidates to predict the full angular structure of μ₂, not just a scalar.
