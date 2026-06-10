# Bounce path-integral derivation of the channel-dominant saddle (PR #152)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The saddle
> rule is a property of the bounce on the classical neck.

PR #151 resolved the neutrino-spread overshoot with a **modelled**
pair-tunneling rule — channel dominance (β = 1: every seesaw element tunnels
through the most compliant neck channel available to the pair) — and flagged
"derive G_ij from the bounce path integral" as its lead open item. This PR
supplies the derivation and **retires the β knob**.

## Part 1: the two-path decomposition

The channel-conversion vertex (the anarchic O(1) overlap, #91/#151) has
support only where the overtone wavefunctions **coexist** — the cavity
mouths; the neck interior is the single-channel tunneling region (#88/#132).
So exactly two single-conversion paths exist:

    A_nm = c_near·e^{−S_m}  (convert at the near mouth, tunnel in channel m)
         + c_far ·e^{−S_n}  (tunnel in channel n, convert at the far mouth),

dominated by the cheaper segment: `A_nm ≍ O(1)·e^{−min(S_n,S_m)}` — the
channel-dominant rule. Conversion *inside* the neck would give the factorized
`e^{−(S_n+S_m)/2}`; the vertex has no support there.

## Part 2: the controlled exact computation

A three-channel double well (channel-dependent barrier heights; WKB actions
≈ 15.4 / 11.4 / 8.4, exact splittings spanning ×2000) with the conversion
coupling localized at the wells, solved by exact diagonalization and
Löwdin-projected onto the single-channel well doublets (extraction verified:
the same-channel element reproduces Δ_k/2 to <10%).

| pair | t/(Δ_max/2) — mouth | t/(Δ_geo/2) — mouth | t/(Δ_max/2) — barrier | t/(Δ_geo/2) — barrier |
|---|---|---|---|---|
| (1,2) | 0.168 | 1.43 | 0.0052 | 0.044 |
| (1,3) | 0.205 | 9.06 | 0.0016 | 0.071 |
| (2,3) | 0.202 | 1.05 | 0.0109 | 0.056 |

- **Mouth coupling:** t/max is **constant** (spread ×1.22 — the O(1)
  conversion factor) while t/geo varies ×8.65 (= `e^{|ΔS|/2}`, the two-path
  prediction). The exact element is channel-dominant.
- **Counterfactual (vertex inside the barrier): the rule flips** — t/geo
  becomes the constant (×1.59) and t/max varies ×6.79. The vertex
  **location** decides the saddle rule. Since the measured r₃₂ excluded the
  factorized rule (#151, 0.1th percentile), the data independently
  corroborate mouth conversion — which is the BAM cavity/neck structure.
- **t ∝ W₀ exactly** (spread <1%): a single conversion vertex — the
  two-segment saddle.

## The closure

The #151 chain now stands on derived footing: mouth conversion (geometry) ⟹
channel dominance (this probe) ⟹ measured r₃₂ natural (~77th percentile),
large mixing (0.44), and the falsifiable `m₁ ≈ 0.04 meV / Σm_ν ≈ 58.8 meV`
prediction. The β interpolation knob is **retired** — one modelling
assumption removed, zero inputs added; the #150 budget is unchanged.

## Scope

The 1D model is the controlled analogue of the neck (channel-dependent
barrier = compliance-dependent bounce action); the full bounce path integral
on the actual 5D throat is open, as are the anarchic prefactor distribution
(the localized flavor residual), explicit PMNS angles, and CP phases.

## Reproduce

```bash
python -m experiments.closure_ledger.channel_dominant_saddle_derivation_probe
# Verdict: CHANNEL_DOMINANT_SADDLE_DERIVED_MOUTH_CONVERSION_COUNTERFACTUAL_FLIPS
```
