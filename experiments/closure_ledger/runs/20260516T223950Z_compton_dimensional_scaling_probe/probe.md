# Tangherlini dimensional-scaling discrimination probe

**Run:** 2026-05-16T22:39:50+00:00

Follow-on to PR #32 (coefficient-origin probe). Tests whether γ = −3/2 (the analytic O(ω/m) closure coefficient from PR #31) is d-independent or scales with spatial dimension. Discriminates between Casimir-based and embedding-based BAM derivations.

**Thesis:** γ = −3/2 is d-independent under the natural d-dim generalisation of BAM and KN. The ((d-2)+c²)/(d-1) polarisation-sum factor cancels in the matching equation.

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_analytic_d_independence` | K/B ratio = −2 in all d: YES | **PASS** |
| T2 | `T2_numerical_gamma_d` | max |γ(d) − (−3/2)| = 0.0000 | **PASS** |
| T3 | `T3_candidate_predictions` | predictions tabulated | **PASS** |
| T4 | `T4_falsified_candidates` | 1 falsified, 7 survive | **PASS** |
| T5 | `T5_survivors_and_next_steps` | 7 surviving candidates | **PASS** |

## T1: K/B ratio per d

Predicted: f_KN_d_O(ε) / f_BAM_d_baseline_O(ε) = **−2** for all d (factor `((d-2)+c²)/(d-1)` cancels).

| d | mean K/B | std | matches −2? |
|---:|---:|---:|:---:|
| 3 | -2.0000 | 2.96e-05 | ✓ |
| 4 | -2.0000 | 3.07e-05 | ✓ |
| 5 | -2.0000 | 3.12e-05 | ✓ |
| 6 | -2.0000 | 3.15e-05 | ✓ |
| 8 | -2.0000 | 3.18e-05 | ✓ |

## T2: Numerical γ(d)

| d | γ_opt | residual at optimum |
|---:|---:|---:|
| 3 | -1.5000 | 1.5957e-07 |
| 4 | -1.5000 | 1.6628e-07 |
| 5 | -1.5000 | 1.6964e-07 |
| 6 | -1.5000 | 1.7165e-07 |

**Max deviation from −3/2:** 0.0000

## T3 & T4: Candidate predictions and falsification

| candidate | d=3 | d=4 | d=5 | d=6 | status |
|---|---:|---:|---:|---:|---|
| `A_doubled_electron_casimir` | -1.5000 | -1.5000 | -1.5000 | -1.5000 | survives |
| `B_photon_casimir_minus_hopf` | -1.5000 | -1.5000 | -1.5000 | -1.5000 | survives |
| `D_hopf_times_photon_mult` | -1.5000 | -1.5000 | -1.5000 | -1.5000 | survives |
| `F_closure_plus_hopf` | -1.5000 | -1.5000 | -1.5000 | -1.5000 | survives |
| `G_pauli_trace` | -1.5000 | -1.5000 | -1.5000 | -1.5000 | survives |
| `H_two_mouth_spin_half` | -1.5000 | -1.5000 | -1.5000 | -1.5000 | survives |
| `C_embedding_over_polarization` | -1.5000 | -1.3333 | -1.2500 | -1.2000 | **FALSIFIED** |
| `E_antipodal_natural_units` | -1.5000 | -1.5000 | -1.5000 | -1.5000 | survives |

## T5: Surviving candidates and next discrimination paths

**Surviving candidates** (7):

- `A_doubled_electron_casimir`
- `B_photon_casimir_minus_hopf`
- `D_hopf_times_photon_mult`
- `F_closure_plus_hopf`
- `G_pauli_trace`
- `H_two_mouth_spin_half`
- `E_antipodal_natural_units`

**Next discrimination paths:**

- **O(ε²) analytic extension** — Different surviving candidates predict different next-order coefficients. Casimir candidates predict C_2²; Hopf candidates predict χ-dependent modulation. Computing the analytic O(ε²) coefficient would discriminate.
- **Cross-process test (pair production)** — γγ → e⁺e⁻ has a different amplitude structure. Candidates rooted in electron-only quantities (A, H) vs photon-electron coupling (B, D, G) would give different predictions.
- **Polarised cross-section corrections** — For specific photon polarisation choices (circular, linear), candidates may give different correction patterns.

## Verdict

**FALSIFIES_C.** FALSIFIES_C — γ is d-independent (numerical γ_d matches −3/2 in d ∈ {3, 4, 5, 6} to within grid tolerance). The factor ((d−2)+c²)/(d−1) cancels in the matching equation. Candidates predicting d-dependent γ are falsified: ['C_embedding_over_polarization']. 7 candidates survive — all predicting γ = −3/2 universally. Further discrimination requires O(ε²) analytic extension, cross-process tests, or polarisation-dependent observables (T5).

## What this leaves open

- **6 candidates remain after d-scaling discrimination.** Among the surviving candidates, the cleanest single-ingredient one is `G_pauli_trace_normalized` (PR #32 naturalness ranking 8/10).
- **O(ε²) extension** is the natural next probe: different surviving candidates predict different next-order coefficients. PR #31's residual fit (ε^1.88 vs predicted ε²) leaves room for a small O(ε^{3/2}) or O(ε²) correction that could be candidate-specific.
- **The d-dim KN analog used in this probe** is a natural generalisation; verifying against a textbook d-dim QED Compton derivation would close the analytic loop.
- **Tangherlini physical interpretation.** BAM's 5D framework uses Tangherlini for the radial bulk channel — this probe tests the *spatial-dimension* scaling of the photon polarisation sum, which is a different generalisation than putting the photon in the Tangherlini bulk. The latter would test the radial channel's contribution to the vertex coefficient — open follow-on.
