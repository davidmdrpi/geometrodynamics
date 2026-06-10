# Neutrino log-bounce sensitivity audit and ε_n overshoot bracket (PR #149)

> **Framing (to avoid a category error).** QFT on the *fixed classical* throat
> geometry — geometry → fields, **not** quantum gravity. The bounce action
> lives on the classical throat; the audit asks how tightly the oscillation
> data pin the per-generation compliance.

PR #113 found that the generation-dependent healing length `ε_n` gets the
neutrino hierarchy **direction** right (`ε_n ∝ 1/χ_n` ⟹ normal ordering,
untuned) but overshoots its **magnitude** ×28: the steep bounce
(`m_ν ∝ ε^{4.8}`, #112) amplifies the natural ×8 χ_n spread into orders of
magnitude in mass. The spread stayed a residual. This PR is the #148-pattern
audit of that residual — bracket it with the program's own structure and the
data, and sharpen the failure mode into an exclusion.

## The mechanism keystones, re-verified

- The bounce length is the horizon tortoise divergence:
  `L*(ε) = (r_s/2)ln(1/ε) + const` — slope re-fit on the metric, 0.500
  (#88/#132).
- The logarithm turns the exponential mass-action relation into a power law:
  `e^{−c·L*(ε)} ∝ ε^{c·r_s/2}` — identity constant over three ε decades; at
  the #112 bounce constant, `p = ∂ln m_ν/∂ln ε = 4.8`.

## Hypersensitivity inverted: the data pin the required profile

Forward, the steepness produced the overshoot (×2 in ε → 2^4.8 ≈ 28× in
mass). Inverted, `δln ε = δln m / p`: the oscillation-data errors
(Δm²₂₁ ± 2.8%, Δm²₃₁ ± 1.1%) compress by 1/4.8 into **~0.3–0.4%** on the
required ε ratios. The residual is not fuzzy — it is a sharply localized
number the geometry has not yet produced (the same pleasant inversion as the
#148 spectrum bound on `k·r_s`).

## The Dirac-attribution endpoints

The data fix only the product `m_D,n · ε_n^p`. Bracketing the attribution
between the **pure-bounce** endpoint (uniform m_D) and the **#113-implied
Dirac growth** (m_D ratios 1.88, 1.48):

| quantity | pure-bounce | #113 m_D | χ-driven |
|---|---|---|---|
| required ε₂/ε₁ | 1.352 | 1.186 | 3.134 |
| required ε₃/ε₂ | 1.435 | 1.323 | 2.487 |

The required spread lies in **[1.32, 1.44]** per generation step; the
principled χ-driven prediction sits far outside both endpoints — excluded by
×1.7–1.9 in ε, i.e. **×14–21 in mass** at this data vintage (#113 reported
×28), against a sub-percent bracket.

## No single power law ε_n ∝ χ_n^{−q}

The per-pair implied exponents disagree under **both** attributions
(q₂₃/q₁₂ = 1.5 pure-bounce, 2.1 with the #113 m_D); the best single q (0.32,
minimax) still misses a mass ratio by ×1.38 — ×1.07 in ε, ~25× the data
bracket; q = 1 overshoots ×21. The spread is **not** a power law in the
boundary stress — #113's "accommodates, not predicts," now a scanned
exclusion.

## Consistency

With the required gentle profile: normal ordering (derived, #113) and
Σm_ν ≈ 61 meV — inside the program's own 59–65 meV window
(`cosmological_sigma_mnu_probe`) and below the cosmology bound. The residual
is consistent, just not derived; it plausibly belongs to the mixing/anarchy
sector (#92).

## Ledger and scope

- **Derived:** the mechanism keystones; normal ordering; the sub-percent
  data bracket on the required profile; the attribution endpoints; the
  power-law exclusion.
- **Residual:** the origin of the gentle ε_n profile (mixing/anarchy, #92).
- **Open:** deriving the profile; the #112 anchor (m₁ = 2.08 meV rides on
  ε ~ R_c³); the m_D attribution (bracketed, not derived). Input budget
  unchanged.

## Reproduce

```bash
python -m experiments.closure_ledger.neutrino_eps_n_overshoot_bracket_probe
# Verdict: NEUTRINO_EPS_N_REQUIRED_PROFILE_PERCENT_BRACKETED_POWER_LAW_EXCLUDED
```
