# Mixing/anarchy origin of the ε_n profile (PR #151)

**Run:** 2026-06-10T04:43:23+00:00

Tests #149's hypothesis that the gentle ε_n profile belongs to the mixing/anarchy sector. The seesaw matrix in the overtone basis — derived χ-driven compliances, derived cavity-floor Dirac growth, anarchic cross-channel overlaps — with the pair-tunneling saddle rule audited from factorized to channel-dominant. Result: channel dominance resolves the measured overshoot (natural at the ~75th percentile), grows large mixing from the same rule, and converts the unmeasured m₂/m₁ into a falsifiable m₁ prediction. *(QFT on the classical throat, not quantum gravity.)*

- **Model**: M_ij = m_D,i m_D,j·c_ij·G_ij(β); floors #91, χ-compliances #113, anarchic c #91
- **Selection**: β = 0 excluded (0.1th pct); β = 1 natural (~75th pct) on the measured r₃₂
- **Mixing**: indicator 0.085 → 0.43 with the same β (cross-channel #91)
- **Prediction**: m₁ ≈ 0.04 meV; Σm_ν ≈ 58.8 meV (vs uniform-anchor 61.1) — falsifiable
- **Relocation**: three-number profile → one saddle rule + anarchic O(1) draw; budget unchanged
- **Open**: derive the saddle rule; PMNS angles/CP; the anarchic draw (flavor puzzle)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | test the #149 mixing/anarchy hypothesis with derived inputs | **PASS** |
| T2 | `T2_model_and_derived_inputs` | cavity floors = the #149 m_D endpoint (<1%, verified identity) | **PASS** |
| T3 | `T3_beta_scan_measured_ratio` | measured r₃₂: β = 0 excluded (0.1th pct); β = 1 natural (~75th) | **PASS** |
| T4 | `T4_large_mixing_from_same_beta` | large mixing from the same β: indicator 0.085 → 0.43 | **PASS** |
| T5 | `T5_m1_prediction` | prediction: m₁ ≈ 0.04 meV, Σm_ν ≈ 58.8 meV (vs 61.1 uniform) | **PASS** |
| T6 | `T6_residual_relocation` | residual relocates: profile → saddle rule + anarchic draw | **PASS** |
| T7 | `T7_scope` | scope: saddle rule modelled; PMNS angles/CP open | **PASS** |
| T8 | `T8_assessment` | EPS_N_FROM_CHANNEL_DOMINANT_ANARCHY_R32_NATURAL_M1_PREDICTED_LIGHT | **PASS** |

## The β scan on the measured ratio (m₃/m₂)

| β (saddle rule) | median r₃₂ | observed r₃₂ | observed percentile |
|---:|---:|---:|---:|
| 0.0 | 113.12 | 5.66 | 0.1% |
| 0.5 | 16.3 | 5.66 | 20.8% |
| 1.0 | 2.75 | 5.66 | 77.2% |

β = 0 (factorized) re-derives the #113/#149 overshoot in matrix form — excluded; β = 1 (channel-dominant: every pair element tunnels through the widest neck) makes the observed ratio anarchy-typical.

## Large mixing emerges from the same rule

| β | median mixing indicator |
|---:|---:|
| 0.0 | 0.085 |
| 0.5 | 0.233 |
| 1.0 | 0.44 |

## The m₁ prediction (the unmeasured ratio becomes falsifiable)

| quantity | mixing solution (β = 1) | uniform anchor (#112) |
|---|---:|---:|
| m₁ (meV) | 0.042 | 2.08 |
| Σm_ν (meV) | 58.8 | 61.1 |

Median r₂₁ = 203.3 at β = 1; the discriminator is Σm_ν at ~1–2 meV cosmology precision (and m_ββ); both scenarios keep normal ordering.

## Verdict

**EPS_N_FROM_CHANNEL_DOMINANT_ANARCHY_R32_NATURAL_M1_PREDICTED_LIGHT.** THE #149 HYPOTHESIS TESTS POSITIVE: CHANNEL-DOMINANT ANARCHIC MIXING IN THE OVERTONE BASIS RESOLVES THE MEASURED OVERSHOOT, PRODUCES LARGE MIXING FROM THE SAME RULE, AND CONVERTS THE UNMEASURED m₂/m₁ INTO A FALSIFIABLE PREDICTION — THE ε_n PROFILE RESIDUAL RELOCATES TO ONE SADDLE RULE PLUS AN ANARCHIC DRAW. #149 bracketed the overshoot and pointed at the mixing/anarchy sector; this probe builds the test.

THE MODEL, DERIVED INPUTS. M_ij = m_D,i m_D,j·c_ij·G_ij(β): the Dirac growth is the #91 cavity floors — which equal the #149 m_D-attribution endpoint to <1% (verified identity) — the channel suppressions are the χ-driven compliances through the #112 bounce, the c_ij are the #91 anarchic cross-channel overlaps, and β interpolates the pair-tunneling saddle from factorized to channel-dominant (the widest neck available to the pair).

THE MEASURED RATIO SELECTS CHANNEL DOMINANCE. The only measured ratio is the heavy pair r₃₂ ≈ 5.7. Factorized (β = 0): ensemble median ≈ 113 — the #113/#149 overshoot re-derived in matrix form; observed at the 0.1th percentile, excluded. Channel-dominant (β = 1): the steep hierarchy collapses out of the heavy pair (every element shares the widest neck); median ≈ 3.0; observed at the ~75th percentile — natural.

ONE RULE, TWO OBSERVABLES. The same β that compresses the spread grows the mixing indicator 0.085 → 0.43 (small-mixing aligned → large-mixing anarchic) — consistent with the #91 cross-channel large-PMNS identification. The spread and the mixing are two faces of one saddle rule.

THE m₁ PREDICTION. The lightest channel keeps its suppression: median r₂₁ ≈ 203.3 at β = 1. m₁ is unmeasured, so this is a PREDICTION, not a misfit: m₁ ≈ 0.042 meV and Σm_ν ≈ 58.8 meV, against the #112 uniform anchor (m₁ = 2.08, Σ = 61.1) — a ~2 meV discriminator for next-generation cosmology, with m_ββ shifted accordingly; both keep normal ordering.

THE RESIDUAL RELOCATES. Before: a three-number fine-tuned profile (bracketed, #149). After: derived compliances + derived floors + one discrete saddle rule + an anarchic O(1) draw — the ratios become percentile-natural statistics rather than deterministic outputs (the flavor puzzle's BAM face, localized). No new continuous knob; the #150 budget is unchanged.

SCOPE. The saddle rule is modelled, not derived from the bounce path integral; the charged-lepton rotation, explicit PMNS angles, and CP phases are open; the anarchic draw itself is the localized residual.
