# Bulk-scale residual audit for k·r_s (PR #148)

**Run:** 2026-06-10T03:59:24+00:00

Makes the #133 bulk-scale bound quantitative: the #127 Tangherlini–AdS₅ background is put under the actual cavity operator, the locked closure precisions are converted into computed upper bounds on k·r_s, and the static-throat existence condition closes the bracket from below. The residual emerges two-sided — the #89 ε pattern: structure and bracket derived, value residual. *(QFT on the classical throat, not quantum gravity.)*

- **Background**: f_k = 1 − r_s²/r² + k²r² (#127); V_l = f[l(l+2)/r² + (3/2r)f′]
- **Scaling**: shifts ∝ (k·r_s)² (exponents 1.98–2.00)
- **Upper bound**: k·r_s ≲ 0.0064 (bridge 0.04%) / 0.0696 (pinhole 2.2%)
- **Lower bound**: k > 0: flat bulk ⟹ B = 4πσ = 0 ⟹ no throat equilibrium (#56/#57)
- **Bracket**: 0 < k·r_s ≲ 0.006–0.07 (the #89 two-sided pattern)
- **Open**: the value of k·r_s (κ₅²/Λ₅, #112/#133); global brane solution (#127)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | make the #133 bound quantitative; bracket k·r_s two-sided | **PASS** |
| T2 | `T2_interpolating_background` | V_l reduces to #116 at k = 0 (exact); pinhole lock residual reproduced | **PASS** |
| T3 | `T3_quadratic_scaling` | shifts ∝ (k·r_s)² — fitted exponents 1.98–2.00 | **PASS** |
| T4 | `T4_spectrum_bound` | spectrum bound: ≲ 0.0064 tight / 0.070 conservative (~16× tighter) | **PASS** |
| T5 | `T5_lower_bound_static_throat_requires_ads` | k = 0 ⟹ no throat equilibrium ⟹ k > 0: bracket two-sided | **PASS** |
| T6 | `T6_ledger_classification` | ledger: structure derived, value residual; budget unchanged | **PASS** |
| T7 | `T7_scope` | scope: bounds, does not derive; domain fixed (ΔR modulus) | **PASS** |
| T8 | `T8_assessment` | BULK_SCALE_K_RS_TWO_SIDED_BRACKET_SPECTRUM_BOUND_VALUE_RESIDUAL | **PASS** |

## The AdS correction enters at (k·r_s)²

| k·r_s | Δω(1,0)/ω | Δω(0,0)/ω | Δpinhole/pinhole |
|---:|---:|---:|---:|
| 0.02 | 0.003975 | 0.004075 | 0.001817 |
| 0.05 | 0.02466 | 0.02527 | 0.01136 |
| 0.1 | 0.09634 | 0.09861 | 0.04552 |
| 0.2 | 0.362 | 0.369 | 0.1834 |

Fitted log-log exponents: ω(1,0) `1.981`, ω(0,0) `1.98`, pinhole `2.001` — quadratic, as the #127 metric implies.

## The locks bound the residual

| lock (tolerance) | sensitivity c | bound on k·r_s |
|---|---:|---:|
| Compton bridge / ε-lock (0.04%) on ω(1,0) | 9.863 | 0.0064 |
| same tolerance on ω(0,0) | 10.108 | 0.0063 |
| pinhole γ-lock residual (2.2%) on Σ V_max | 4.544 | 0.0696 |

The #133 estimate was `≲ 0.1`; the tight bound `0.0064` is ~15.7× sharper. Combined with the lower bound (T5): **0 < k·r_s ≲ 0.006–0.07 (the #89 two-sided pattern)**.

## The flat-bulk limit has no throat (lower bound)

| B/B₀ | R* = (A/2B)^{1/3} |
|---:|---:|
| 1.0 | 0.63 |
| 0.1 | 1.357 |
| 0.01 | 2.924 |
| 0.001 | 6.3 |

As k → 0 the cohesive tension vanishes (#56/#57) and the equilibrium runs away — a static throat requires k > 0.

## Verdict

**BULK_SCALE_K_RS_TWO_SIDED_BRACKET_SPECTRUM_BOUND_VALUE_RESIDUAL.** THE BULK-SCALE RESIDUAL k·r_s IS BRACKETED TWO-SIDED BY THE PROGRAM'S OWN STRUCTURE: THE LOCKED CAVITY SPECTRUM BOUNDS IT ABOVE (k·r_s ≲ 0.0064 TIGHT, 0.0696 CONSERVATIVE — UP TO ~15.7× TIGHTER THAN THE #133 ESTIMATE) AND THE EXISTENCE OF THE STATIC THROAT BOUNDS IT BELOW (k > 0); ONLY THE VALUE STAYS RESIDUAL. #133 isolated the recurring κ₅²/Λ₅ residual as one bounded number with an order-of-magnitude estimate; this audit computes the bound.

THE BACKGROUND UNDER THE OPERATOR. The #127 interpolating metric f_k = 1 − r_s²/r² + k²r² carries the potential V_l = f[l(l+2)/r² + (3/2r)f'], which reduces to the #116 Tangherlini potential at machine precision when k = 0, and the k = 0 pinhole operator reproduces the documented −2.2% γ-lock residual — the audit stands on the locked machinery itself.

QUADRATIC SCALING, DERIVED. Every observable shifts as (k·r_s)² (fitted exponents 1.98–2.00 across a decade of k·r_s): the structure of the bound is geometry; only its tolerance comes from the locks.

THE SPECTRUM BOUND. With sensitivities c ≈ 9.9 (ω(1,0)) and ≈ 4.5 (pinhole), the locked precisions give k·r_s ≤ √(tol/c): the 0.04% Compton-bridge closure ⟹ ≲ 0.0064; the 2.2% pinhole lock ⟹ ≲ 0.070. The throat sits deep in the near-flat AdS region — why pure Tangherlini approximated everything so well, now quantified.

THE LOWER BOUND. k → 0 kills the throat: λ_crit → 0 (#57), B = 4πσ → 0 (#56), R* = (A/2B)^{1/3} → ∞, E(R) = A/R monotone (no minimum, computed). A static throat exists only for k > 0 — the bracket is two-sided, the #89 ε pattern.

LEDGER. Derived: the AdS sign, the quadratic scaling, the bracket. Residual: the value (κ₅²/Λ₅, #112). The input budget is unchanged — the one bounded bulk number is now bracketed by the program's own locks.

SCOPE. The audit bounds, it does not derive: the absolute normalisation stays open (#112/#133), the cavity domain is held fixed (ΔR modulus), and the ε-region horizon shift is neglected.
