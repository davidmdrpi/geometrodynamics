# BAM throat-vertex loop probe — the Schwinger anomalous moment

**Run:** 2026-05-24T03:47:56+00:00

Constructs the one-loop throat-vertex correction in BAM-native terms (the throat dressing its moment by one S³-photon self-exchange) and shows it reproduces the Schwinger anomaly a = α/2π. A reconstruction using the tree-normalized BAM primitives — honest scope stated in T6.

- **Anomaly**: a = α/2π (Schwinger, one loop)
- **Loop**: throat self-exchange: S³ virtual photon at the throat-pinch vertex
- **Integral**: `I = ∫₀¹ 2z dz = 1 ⟹ F₂(0) = α/2π`
- **Scope**: reconstruction (tree-normalized primitives); 1/2π from S_BAM open
- **B4 caveat**: a, g dimensionless; α the coupling; scale = single anchor

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_bam_one_loop_vertex_structure` | throat self-exchange triangle; loop param α/π | **PASS** |
| T2 | `T2_virtual_photon_is_s3_exchange` | S³ photon → 1/(4πd) (dev 5e-03) | **PASS** |
| T3 | `T3_schwinger_feynman_parameter_integral` | I = ∫₀¹2z dz = 1.0000 → F₂(0)=α/2π | **PASS** |
| T4 | `T4_anomaly_vs_experiment` | a=α/2π=0.0011614 vs a_e (Δ 0.15%) | **PASS** |
| T5 | `T5_one_loop_g_factor` | g = 2(1+a) = 2.0023228 | **PASS** |
| T6 | `T6_honest_scope` | reconstruction; 1/2π from S_BAM open | **PASS** |
| T7 | `T7_b4_accounting` | a, g dimensionless; α coupling; scale = anchor | **PASS** |
| T8 | `T8_assessment` | throat loop reproduces Schwinger | **PASS** |

## T1: BAM one-loop vertex structure

- **external_lines**: in/out throat (electron) at the pinch vertex
- **virtual_photon**: S³ Green-function exchange (flat limit 1/q²)
- **vertex**: throat-pinch F² vertex (PRs #38–#41)
- **coupling**: α (Hopf charge structure)
- **loop_parameter**: α/π (coupling α × loop phase space 1/π)

## T2: Virtual photon = S³ exchange

| ψ | G(ψ) | G·4π d (→ 1) |
|---:|---:|---:|
| 1e-02 | 7.9195 | 0.995192 |
| 1e-03 | 79.5394 | 0.999522 |
| 1e-04 | 795.7367 | 0.999952 |

Flat limit → Coulomb 1/(4π d) = 1/q² propagator (max dev 4.8e-03).

## T3: Schwinger Feynman-parameter integral

- I (reduced ∫₀¹ 2z dz) = 1.000000
- I (simplex form) = 1.000000 (= 1)
- F₂(0) = (α/2π)·I = 0.00116141

## T4: a = α/2π vs experiment

- a = α/2π = 0.00116141
- a_e (measured) = 0.00115965
- relative difference = 0.15% (leading term; higher orders make up the remainder)

## T5: One-loop g

- g = 2(1 + a) = 2.00232282 (tree g=2 + Schwinger a)

## T6: Honest scope

BAM-native:
  - virtual photon = S³ Green-function exchange
  - vertex = throat-pinch F²
  - structure = throat dressing its moment by one self-exchange
Inherited from tree normalization (#35–#46):
  - the 1/q² propagator normalization
  - the vertex normalization (tree QED match, #35–#46)
  - the coupling α
Open (first-principles): 1/2π from a covariant S_BAM loop measure

## T7: B4 accounting

- a = 0.0011614 (dimensionless); α the coupling; scale = single anchor

## T8: Assessment

- anomaly: a = α/2π (one loop)
- integral: I = ∫₀¹ 2z dz = 1
- one-loop g: 2.0023228
- scope: reconstruction (tree-normalized primitives), not first-principles 1/2π
- remaining: 1/2π from covariant S_BAM loop; higher-order a_e; throat spinor

## Verdict

**SCHWINGER_RECONSTRUCTED.** SCHWINGER RECONSTRUCTED. The BAM throat-vertex loop reproduces the Schwinger anomalous magnetic moment a = α/2π, the one-loop capstone of the spin sector (tree g=2 in PR #61).

THE LOOP. The one-loop vertex correction, in BAM terms, is the throat dressing its magnetic moment by one virtual-photon self-exchange: external in/out throat lines at the throat-pinch vertex (PRs #38–#41), the virtual photon an S³ Green-function exchange whose flat limit G → 1/(4π d) is the 1/q² photon propagator (PRs #45–#46), and the Hopf coupling α. The loop expansion parameter is α/π.

THE COEFFICIENT. After the loop integral the anomalous form factor is F₂(0) = (α/2π)·I with the Feynman-parameter integral I = ∫dxdydz δ(x+y+z−1) 2z(1−z)/(1−z)² = ∫₀¹ 2z dz = 1 (the x,y simplex slice measure (1−z) cancels the 1/(1−z)). So a = F₂(0) = α/2π ≈ 0.0011614, matching the measured a_e = 0.00115965 at leading order (~0.15%; higher orders α², α³ make up the remainder), and g = 2(1 + a) = 2.00232….

HONEST SCOPE. This is a RECONSTRUCTION, not an independent first-principles derivation of 1/2π. The one-loop integrand is the QED vertex integrand expressed in BAM primitives that were normalized to QED at tree level (the S³ propagator, the throat-pinch vertex, the coupling α); the probe verifies the loop integral gives the Schwinger coefficient. What is BAM-native: the loop pieces (virtual photon = S³ exchange; vertex = throat pinch) and the structure (one self-exchange dressing the moment). What remains open: the 1/2π from a pure-geometry loop measure derived from S_BAM (the full covariant throat loop), plus the higher-order a_e series. B4: a and g are dimensionless; α is the coupling; the scale is the single anchor.

## What this leaves open

- **1/2π from a first-principles S_BAM loop measure.** The covariant throat loop derived from the action, rather than the QED-normalized integrand — the genuine derivation.
- **Higher-order a_e.** The α², α³, … terms (two- and three-loop).
- **The explicit throat spinor / vertex from S_BAM** (shared with #59–#61).
