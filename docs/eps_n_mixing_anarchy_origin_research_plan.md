# Mixing/anarchy origin of the ε_n profile (PR #151)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The seesaw
> matrix lives in the overtone basis of the classical cavity.

PR #149 bracketed the ε_n overshoot and concluded the gentle required
profile "plausibly belongs to the mixing/anarchy sector (#92)". This PR
builds the test — and it comes back **positive on the measured ratio**, with
a falsifiable prediction on the unmeasured one.

## The model (all hierarchical inputs derived)

The Majorana seesaw matrix in the overtone basis:

    M_ij = m_D,i · m_D,j · c_ij · G_ij(β)

- **Dirac growth** = the #91 cavity floors (1.055, 1.974, 2.894) — which
  equal the #149 m_D-attribution endpoint to <1% (a verified identity
  between independently-introduced numbers);
- **channel suppressions** F_n = (ε_n/ε₁)^p with the χ-driven
  `ε_n ∝ 1/χ_n` (#113) through the p = 4.8 bounce (#112);
- **c_ij** = anarchic O(1) cross-channel overlaps (the #91 large-PMNS
  structure; seeded ensemble of 2000);
- **β** interpolates the pair-tunneling saddle rule:
  `G_ij = exp[(1−β)·(lnF_i+lnF_j)/2 + β·max(lnF_i, lnF_j)]` — factorized
  (β = 0) to **channel-dominant** (β = 1: every pair element tunnels through
  the widest neck available to the pair, the least-action saddle).

## Results

| β | median r₃₂ | observed (5.66) percentile | mixing indicator |
|---|---|---|---|
| 0.0 (factorized) | 113 | **0.1% — excluded** | 0.085 |
| 0.5 | 16.3 | 20.8% | 0.233 |
| 1.0 (channel-dominant) | 2.75 | **77% — natural** | 0.44 |

1. **The measured overshoot resolves.** β = 0 re-derives the #113/#149
   overshoot in matrix form (excluded at the 0.1th percentile); β = 1
   collapses the steep hierarchy out of the heavy pair — the observed
   m₃/m₂ is anarchy-typical.
2. **Large mixing emerges from the same rule** (0.085 → 0.44) — one
   structural choice moves both observables toward the data, consistent
   with the #91 cross-channel large-PMNS identification.
3. **The unmeasured ratio becomes a prediction.** r₂₁ stays ≈ 200 at β = 1;
   since m₁ is unmeasured this is not a misfit but a discriminator:
   **m₁ ≈ 0.04 meV, Σm_ν ≈ 58.8 meV** (mixing solution) vs m₁ = 2.08 meV,
   Σm_ν = 61.1 meV (the #112 uniform anchor) — distinguishable at ~1–2 meV
   cosmology precision, with m_ββ shifted accordingly; both keep normal
   ordering.

## The residual relocation

Before (#149): a three-number fine-tuned ε_n profile, bracketed but
underived. After: derived compliances + derived floors + **one discrete
saddle rule** (selected by the measured ratio) + an anarchic O(1) draw. The
ratios become percentile-natural statistics — the flavor puzzle's BAM face,
now localized to the draw. No new continuous knob; the #150 input budget is
unchanged.

## Scope

The saddle rule G_ij is modelled (motivated as least-action through the
widest neck; deriving it from the bounce path integral is open — the #136
posture). The mixing indicator is ν-side only; explicit PMNS angles and
CP/Majorana phases are open. The m₁ prediction carries the ensemble spread.

## Reproduce

```bash
python -m experiments.closure_ledger.eps_n_mixing_anarchy_origin_probe
# Verdict: EPS_N_FROM_CHANNEL_DOMINANT_ANARCHY_R32_NATURAL_M1_PREDICTED_LIGHT
```
