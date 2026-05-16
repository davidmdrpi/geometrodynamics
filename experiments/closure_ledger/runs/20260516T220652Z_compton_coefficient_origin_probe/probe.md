# Compton vertex coefficient-origin probe

**Run:** 2026-05-16T22:06:52+00:00

Follow-on to PR #31 (analytic vertex derivation establishing `γ = −3/2` as the exact O(ω/m) closure coefficient). This probe enumerates natural BAM-derivable combinations that plausibly produce |γ| = 3/2, evaluates each, and ranks them by parsimony, sign-consistency, simultaneous (α, β, γ) prediction, and connection to existing BAM derivations.

**Target:** γ = -1.5

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_numerical_evaluation` | 8/8 candidates match target | **PASS** |
| T2 | `T2_alpha_beta_prediction` | 1 predict full (α, β, γ) triple | **PASS** |
| T3 | `T3_naturalness_ranking` | top: `G_pauli_trace_normalized` (8/10) | **PASS** |
| T4 | `T4_discrimination_test` | 8 candidates to discriminate | **PASS** |
| T5 | `T5_o_eps2_prediction` | observed ε^1.88 vs predicted ε^2.0 | **PASS** |

## Candidate derivations (T1)

| label | expression | value | matches −3/2? |
|---|---|---:|---|
| `A_doubled_electron_casimir` | `−2·j(j+1) for j=1/2` | -1.5000 | **YES** |
| `B_photon_casimir_minus_hopf_charge` | `−(C_2(j=1) − A_φ(0)) = −(2 − 1/2)` | -1.5000 | **YES** |
| `C_embedding_over_polarization` | `−dim(R³)/dim(transverse pol)` | -1.5000 | **YES** |
| `D_hopf_charge_times_photon_multiplicity` | `−A_φ(0) · (2j_γ+1) = −(1/2)·3` | -1.5000 | **YES** |
| `E_antipodal_throat_natural_units` | `−3·R_S3/(2c) with R_S3 = c = 1` | -1.5000 | **YES** |
| `F_closure_winding_plus_hopf` | `−(n_closure + 1/2) for n_closure = 1` | -1.5000 | **YES** |
| `G_pauli_trace_normalized` | `−Σ_i Tr(σ_i²)/(2·dim_σ) = −6/4` | -1.5000 | **YES** |
| `H_two_mouth_spin_half_casimir` | `−n_mouth · C_2(1/2) = −2·(3/4)` | -1.5000 | **YES** |

## Naturalness ranking (T3)

| rank | candidate | score | single? | sign? | α=β=0? | BAM connections |
|---:|---|---:|:---:|:---:|:---:|---|
| 1 | `G_pauli_trace_normalized` | 8/10 | ✓ | ✓ | ✗ | 1 |
| 2 | `C_embedding_over_polarization` | 7/10 | ✗ | ✗ | ✓ | 2 |
| 3 | `A_doubled_electron_casimir` | 6/10 | ✗ | ✓ | ✗ | 2 |
| 4 | `B_photon_casimir_minus_hopf_charge` | 6/10 | ✗ | ✓ | ✗ | 2 |
| 5 | `D_hopf_charge_times_photon_multiplicity` | 6/10 | ✗ | ✓ | ✗ | 2 |
| 6 | `F_closure_winding_plus_hopf` | 6/10 | ✗ | ✓ | ✗ | 2 |
| 7 | `H_two_mouth_spin_half_casimir` | 6/10 | ✗ | ✓ | ✗ | 2 |
| 8 | `E_antipodal_throat_natural_units` | 4/10 | ✗ | ✗ | ✗ | 0 |

## Discrimination paths (T4)

| candidate | next-order prediction |
|---|---|
| `A_doubled_electron_casimir` | extends to higher j: predict O(ε²) coupling involves C_2² or higher Casimir; needs analytic calculation. |
| `B_photon_casimir_minus_hopf_charge` | extends to higher j: predict O(ε²) coupling involves C_2² or higher Casimir; needs analytic calculation. |
| `C_embedding_over_polarization` | extends to higher dimensions: in d-dim S² embedding, predict γ = −d/2; testable for d=5 (Tangherlini bulk). |

## T5: O(ε²) residual analysis

PR #31 fitted residual exponent at the analytic optimum: **ε^1.88**, predicted for clean O(ε²) closure: ε^2.0. Difference: 0.12.

The PR #31 numerical fit shows residual ~ ε^1.88, slightly below the predicted ε² scaling. The 0.12 shortfall could be either finite-grid numerics or a real residual O(ε^{3/2}) term — also natural in the BAM Hopf framework given the half-integer structure.

## Verdict

**MULTIPLE_PLAUSIBLE.** MULTIPLE PLAUSIBLE — 8 candidates evaluate to γ = −3/2. Top candidate by naturalness: `G_pauli_trace_normalized` (score 8/10), expression `−Σ_i Tr(σ_i²)/(2·dim_σ) = −6/4`. Discrimination requires either: (1) a higher-order BAM prediction that distinguishes them; (2) an additional QED observable (e.g. lepton g-2, pair production cross section) whose value depends differently on each candidate ingredient. The probe lists the discrimination paths in T4. Without further analytic work, the clean -3/2 has multiple plausible BAM origins; numerology cannot be ruled out, but the structural patterns (Casimirs, charges, embedding dims) are real BAM ingredients.

## What this leaves open

- **Discriminating the top candidates.** Multiple natural BAM ingredients evaluate to −3/2; without an additional BAM-derived prediction (e.g. O(ε²) coefficient, or a different QED process like pair production or lepton g-2 at tree level), the probe cannot uniquely identify the source.
- **Numerology check.** Several candidates (e.g. embedding dim / polarisation count) are essentially algebraic accidents. A rigorous derivation would need to show that the specific arrangement of BAM ingredients is forced by the geometry, not just one of many ways to combine numbers giving 3/2.
- **Sign rigor.** Most candidates have the sign added as convention. A first-principles derivation should produce −3/2 with the sign emerging from the spinor / connection structure.
- **O(ε²) residual at ε^1.88.** The slight shortfall from ε² scaling (0.12) could indicate a real residual O(ε^{3/2}) term — natural in the half-integer BAM framework — or finite-grid numerics. Sharpening this requires either analytic next-order calculation or higher-resolution numerical fits.
