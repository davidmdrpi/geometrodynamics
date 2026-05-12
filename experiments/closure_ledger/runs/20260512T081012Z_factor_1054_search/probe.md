# Closed-form search for the 1.054 factor

**Run:** 2026-05-12T08:10:12+00:00

**Targets** (ω(l=1, n=0)):

- Canonical R_OUTER = 1.26: **ω = 1.054727** (+5.473% from 1)
- γ-lock R_OUTER ≈ 1.2623 (cross-species fixed point): **ω = 1.053694** (+5.369% from 1)

**Candidates scanned:** 432 across four families (small rationals, roots of rationals, k₅-and-π forms, series-truncation `1 + ε/k₅^m`). A candidate is 'high-plausibility' if it uses only small integers (≤ 30) and / or k₅ = 5 directly.

## Best per family

| family | best formula | value | %Δ vs γ-lock | plausibility |
|---|---|---:|---:|---|
| `rational` | `157/149` | 1.053691 | -0.0003% | medium |
| `root` | `(75/52)^(1/7)` | 1.053714 | +0.0019% | medium |
| `k5_pi` | `10/9` | 1.111111 | +0.0378% | high |
| `series` | `1 + 7/k_5^3` | 1.056000 | +0.2188% | high |

## High-plausibility candidates within 0.1% of γ-lock

| formula | family | target | value | %Δ vs γ-lock |
|---|---|---|---:|---:|
| `10/9` | rational | omega_sq | 1.111111 | +0.0378% |
| `20/18` | rational | omega_sq | 1.111111 | +0.0378% |
| `10/9` | k5_pi | omega_sq | 1.111111 | +0.0378% |
| `(10/9)^(1/2)` | k5_pi | omega | 1.054093 | +0.0378% |

## Verdict

**No high-plausibility candidate matches within 0.01%.**

Overall best (no plausibility filter): `157/149` (rational) at -0.0003% — but typically with large integers (> 30) that lack obvious structural meaning.

Best high-plausibility near-miss: `10/9` at +0.0378% — suggestive but not exact.

**This is the THESIS.md-flagged clean negative result.** If the 1.054 factor has no closed form in natural BAM ingredients, the dimensional bridge to ℏ retains an irreducible structural constant that must be anchored externally (via m_e). BAM remains **dimensional-ratio-complete and dimensional-scale-incomplete**.

Caveat: 'closed form' depends on what's considered natural. The search restricts to (small integers, k₅, π); a wider net (e.g. involving the antipodal-closure quantum 100, ratios of Tangherlini ground-mode eigenvalues, etc.) might surface a more structural expression. The probe is exhaustive only within the catalog scanned.