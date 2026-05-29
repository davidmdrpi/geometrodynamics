# Generation spread + PMNS mixing sector (PR #91)

**Run:** 2026-05-29T03:26:17+00:00

PR #90 closed the neutrino-mass *scale* but left two residuals: the generation spread (a uniform bounce gives `m_ν ∝ m_D`, ×2.7) and the large PMNS mixing. This probe addresses both. **Headline:** large PMNS vs small CKM is the BAM **cross-channel** (leptons: throat-winding × cavity-resolving) vs **intra-channel** (quarks: shell × shell) distinction; the spread is widened by the overtone-dependent throat-neck coupling (PR #79 `χ_n`).

- **Identification**: large PMNS = cross-channel (charged throat-winding × neutrino cavity-resolving); small CKM = intra-channel (shell × shell); generation spread widened by overtone-dependent neck coupling (PR #79 χ_n decreasing with n)
- **Bare prediction**: normal ordering, m_ν ∝ m_D (cavity floors 1:1.87:2.74)
- **Spread mechanism**: higher-n less throat-coupled (χ_n↓) ⟹ less suppressed ⟹ heavier
- **Mixing dichotomy**: leptons cross-channel (large), quarks intra-channel (small)
- **Open**: precise spectrum (ε_n(χ_n) O(1)), explicit angles, CP/Majorana phases

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_recap_residuals` | uniform S ⟹ m_ν ∝ m_D (×2.7); PMNS large | **PASS** |
| T2 | `T2_overtone_mass_prediction` | generations = overtones; bare m_ν ∝ m_D (1:1.87:2.74) | **PASS** |
| T3 | `T3_spread_from_overtone_neck_coupling` | χ_n decreases with n ⟹ higher-n less suppressed ⟹ heavier | **PASS** |
| T4 | `T4_pmns_large_cross_channel` | PMNS large = cross-channel (k≠0 winding × k=0 cavity) | **PASS** |
| T5 | `T5_ckm_small_intra_channel` | CKM small = intra-channel (up & down both shell) | **PASS** |
| T6 | `T6_honest_data_context` | only Δm² measured; spread is a prediction not a fit | **PASS** |
| T7 | `T7_honest_scope` | spread direction + dichotomy structural; angles open | **PASS** |
| T8 | `T8_assessment` | PMNS_CROSS_CHANNEL_CKM_INTRA_CHANNEL_SPREAD_FROM_OVERTONE_COUPLING | **PASS** |

## T2–T3: Bare prediction and the spread mechanism

- Bare (uniform-bounce) prediction: `m_ν ∝ m_D` (cavity floors), normal ordering `1.00 : 1.87 : 2.74`
- Bare `m₂:m₃ = 1:1.47`; Δm²-implied (m₁≈0) `1:5.77` ⟹ spread needs widening

| n (gen) | cavity floor m_D | boundary stress χ_n |
|---:|---:|---:|
| 0 (1) | 1.055 | 0.304 |
| 1 (2) | 1.974 | 0.097 |
| 2 (3) | 2.894 | 0.039 |

`χ_n` decreases with `n`, so higher-overtone neutrinos are less throat-coupled ⟹ more compliant ⟹ less suppressed ⟹ **relatively heavier** — lifting `m₃` toward the observed spread. Right direction; exact magnitude needs the `ε_n(χ_n)` coefficient.

## T4–T5: The mixing dichotomy

| sector | basis 1 | basis 2 | channels | mixing |
|---|---|---|---|---|
| **leptons (PMNS)** | charged: throat-winding (k≠0) | neutrino: cavity-resolving (k=0) | **different** | **large** |
| **quarks (CKM)** | up: cavity-shell (k=0) | down: cavity-shell (k=0) | **same** | **small** |

Leptons mix **across** the throat-winding / cavity-resolving divide (PR #83's two channels) ⟹ large PMNS. Quarks mix **within** the cavity-shell channel ⟹ small CKM. This is the BAM-native reason `PMNS ≫ CKM`.

## Verdict

**PMNS_CROSS_CHANNEL_CKM_INTRA_CHANNEL_SPREAD_FROM_OVERTONE_COUPLING.** LARGE PMNS IS CROSS-CHANNEL, SMALL CKM IS INTRA-CHANNEL; THE GENERATION SPREAD IS WIDENED BY THE OVERTONE-DEPENDENT NECK COUPLING. PR #90 closed the neutrino-mass SCALE from bulk geometry but left a generation-uniform bounce (m_ν ∝ m_D, a ×2.7 spread) and the large PMNS mixing. This probe addresses both with BAM-native structure.

GENERATIONS = OVERTONES; THE BARE PREDICTION. The neutrino generations are the cavity radial overtones n = 0, 1, 2 (k = 0, PR #83). A generation-uniform bounce gives m_ν,g ∝ m_D,g = √(ω²(0,n)), the cavity floors, in the ratio 1 : 1.87 : 2.74 — normal ordering, a falsifiable prediction. (Only Δm² are measured; the absolute scale is unknown, so this is a prediction, not a fit.)

THE SPREAD: OVERTONE-DEPENDENT NECK COUPLING. The bounce suppression S_n grows with how strongly overtone n couples to the throat neck (PR #88). That coupling is exactly PR #79's Z₂-antisymmetric boundary stress χ_n = (T_inner−T_outer)/2, which DECREASES monotonically with n (≈ 0.304, 0.097, 0.039): higher overtones fill the cavity more uniformly and feel the mouth asymmetry less. So higher-n neutrinos are LESS throat-coupled ⟹ a more compliant neck ⟹ a smaller bounce ⟹ LESS suppression ⟹ relatively HEAVIER. This widens the spread in the right direction, lifting m₃ relative to m₂ (from the bare m₂:m₃ = 1:1.47 toward the Δm²-implied 1:5.8). The mechanism and direction are BAM-native; the exact magnitude needs the ε_n(χ_n) coefficient (O(1), not derived).

PMNS LARGE = CROSS-CHANNEL. The PMNS matrix is the overlap of the charged-lepton mass basis (throat-winding, k≠0, throat-localised) with the neutrino mass basis (cavity-resolving, k=0, shell-distributed) — DIFFERENT channels of the PR #83 unified operator, related only by the throat↔shell Z₂ / the +3 shift (PR #82). Two bases built from different quantum numbers (k vs n) are generically strongly misaligned ⟹ LARGE PMNS mixing.

CKM SMALL = INTRA-CHANNEL. Up-type and down-type quarks are BOTH cavity-shell modes (k=0, n≥3; PR #85 quadrant, PR #82 quark pair). Same channel ⟹ nearly aligned mass bases ⟹ SMALL CKM mixing. So the long-standing puzzle — why is lepton mixing large and quark mixing small — is the BAM cross-channel (leptons: throat-winding × cavity-resolving) vs intra-channel (quarks: shell × shell) distinction.

HONEST SCOPE. ESTABLISHED (BAM-native): generations are overtones, so the bare prediction is normal ordering with m_ν ∝ m_D (cavity-floor ratios 1:1.87:2.74); the spread is widened in the right direction by the overtone-dependent neck coupling (PR #79 χ_n); and large PMNS vs small CKM is the cross-channel vs intra-channel distinction. NOT established: the precise mass spectrum (the ε_n(χ_n) coefficient is O(1), not derived; the absolute scale is unmeasured — only Δm²), the explicit PMNS angles (need the cross-channel overlap integrals; large but not computed to specific angles here), and the CP / Majorana phases.

## What this leaves open

- **The precise mass spectrum** — the `ε_n(χ_n)` coefficient is `O(1)`, not derived; the absolute scale is unmeasured (only `Δm²`).
- **The explicit PMNS angles** — need the cross-channel overlap integrals (large, but not computed to specific angles here).
- **The CP / Majorana phases.**
