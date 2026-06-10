# Bulk-scale residual audit for k·r_s (PR #148)

> **Framing (to avoid a category error).** QFT on the *fixed classical* throat
> geometry — geometry → fields, **not** quantum gravity. The AdS scale k is a
> parameter of the classical 5D bulk; the audit asks how much of it the locked
> matter spectrum tolerates.

The bulk-scale ledger (#133) reduced the recurring `κ₅²/Λ₅` residual to one
bounded dimensionless number — the AdS scale in throat units,
`k·r_s = R_MID/L_AdS` — with an order-of-magnitude estimate `≲ 0.1` read off
the #127 metric correction `(k·r)² ~ O(10⁻²)`. This PR makes the bound
**quantitative** and closes the bracket from below.

## Method

The #127 interpolating Schwarzschild–Tangherlini–AdS₅ metric
`f_k(r) = 1 − r_s²/r² + k²r²` (Einstein with `Λ₅ = −6k²`) is placed under the
actual cavity operator with the potential `V_l = f[l(l+2)/r² + (3/2r)f′]`,
which reduces to the #116 Tangherlini potential at **machine precision** when
k = 0. The cavity domain is held fixed (the geometry ratios are fixed by the
ΔR modulus, #133); the tortoise coordinate is rebuilt from `f_k` at each k. As
a cross-check, the k = 0 pinhole operator reproduces the documented γ-lock
residual (Σ V_max = 22.02, −2.1% off the locked 22.5).

## Results

- **Quadratic scaling, derived:** every locked spectral observable — ω(1,0),
  ω(0,0), the pinhole sum — shifts as `(k·r_s)²` (fitted log-log exponents
  1.98–2.00 across a decade). The structure of the bound is geometry; only
  its tolerance comes from the locks.
- **The spectrum bound:** with sensitivities `c = |Δobs/obs|/(k·r_s)²`
  computed (≈ 9.9 for ω(1,0), ≈ 4.5 for the pinhole), each locked precision
  implies `k·r_s ≤ √(tol/c)`:

  | lock | bound on k·r_s |
  |---|---|
  | Compton bridge / ε-lock closure (0.04%) on ω(1,0) | **0.0064** |
  | same tolerance on ω(0,0) | 0.0063 |
  | pinhole γ-lock residual (2.2%) on Σ V_max | 0.0696 |

  Both inside the #133 estimate; the tight bound is **~16× sharper**. The
  throat sits deep in the near-flat AdS region — *why* pure Tangherlini
  (#116/#127) approximated everything so well, now quantified by the locks
  themselves.
- **The lower bound:** `k → 0` kills the throat — `Λ₅ → 0` ⟹ the RS tension
  `λ_crit = √(6|Λ₅|)/κ₅² → 0` (#57) ⟹ the cohesive coefficient `B = 4πσ → 0`
  (#56) ⟹ `R* = (A/2B)^{1/3} → ∞` and `E(R) = A/R` is monotone (no minimum,
  computed). A static throat exists **only for k > 0**.

## The bracket

    0 < k·r_s ≲ 0.0064 (tight) … 0.070 (conservative)

— the same two-sided epistemic shape as the neutrino-compliance bracket
`ε ∈ [2π, k₅√(2π)]` (#89): **structure and bracket derived, value residual.**

## Ledger

- **Derived:** the AdS sign (static-throat necessity), the `(k·r_s)²`
  scaling, the two-sided bracket.
- **Residual:** the value of `k·r_s` (= `κ₅²/Λ₅` in throat units, #112).
- **Budget unchanged:** `{G anchor} + {n_part, √σ/m_e, ε, α}` + the flavor
  puzzle — the one bounded bulk number is now bracketed by the program's own
  locks instead of an estimate.

## Scope

The audit bounds the residual; it does **not** derive its value. The cavity
domain is held fixed (ΔR modulus); the `O(k²r_s²)` horizon shift inside the
ε-healing region is neglected; the lock-to-observable pairing maps each
closure to the operator it rides on (bridge ↔ ω(1,0), pinhole ↔ Σ V_max).

## Reproduce

```bash
python -m experiments.closure_ledger.bulk_scale_k_rs_audit_probe
# Verdict: BULK_SCALE_K_RS_TWO_SIDED_BRACKET_SPECTRUM_BOUND_VALUE_RESIDUAL
```
