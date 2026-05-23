# Self-consistent throat radius / finite self-energy probe

**Run:** 2026-05-23T19:13:26+00:00

Targets the remaining BAM scale anchor — the throat radius R_MID (≈ the invariant bulk separation ΔR) — recasting it from an imposed constant into a finite-self-energy stable equilibrium, honestly against the B4 scale-modulus theorem (PR #52).

- **Target**: remaining BAM scale anchor: throat radius R_MID (≈ ΔR)
- **Finite self-energy**: U_EM = α ℏc/(2 R_MID); throat caps the field; U/mc² = α/2
- **Equilibrium**: `E(R)=A/R+B R² → stable R*=(A/2B)^(1/3)`
- **B4 caveat**: absolute R* rides on one dimensionful coupling (or a relation to α)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_finite_self_energy_throat_cutoff` | throat caps field → U_EM finite; point charge diverges | **PASS** |
| T2 | `T2_s3_farfield_regular` | S³ Green regular at antipode (zero-mean background) | **PASS** |
| T3 | `T3_scale_free_obstruction` | EM-only 1/R monotone → no equilibrium (B4) | **PASS** |
| T4 | `T4_self_consistent_equilibrium` | E=A/R+BR² → stable R*=(A/2B)^⅓ (err 8e-06) | **PASS** |
| T5 | `T5_scale_modulus` | B→B/λ³ ⟹ R*→λR* (one modulus) | **PASS** |
| T6 | `T6_renormalization_free` | U_EM/mc² = α/2 = 0.0036 (finite, no UV div) | **PASS** |
| T7 | `T7_pure_em_relocation` | mc²=U_EM R-independent ⟹ g=2/α=274.1 | **PASS** |
| T8 | `T8_assessment` | equilibrium recasts anchor; value still one input | **PASS** |

## T1: Finite self-energy — throat caps the field

| point-charge cutoff (m) | self-energy (J) |
|---:|---:|
| 1e-13 | 1.154e-15 |
| 1e-15 | 1.154e-13 |
| 1e-17 | 1.154e-11 |
| 1e-19 | 1.154e-09 |

Throat (capped at R_MID = 3.862e-13 m): U_EM = 2.987e-16 J (finite).

## T2: S³ far-field regular

- G near antipode: ['-0.03799', '-0.03800', '-0.03800']
- analytic G(π) = (−1−½)/(4π²R) = -0.03800 (finite; zero-mean −½ background)

## T3: Scale-free obstruction (EM-only 1/R)

| R/λ_C | E=A/R (J) | dE/dR |
|---:|---:|---:|
| 0.5 | 5.974e-16 | -3.094e-03 |
| 1.0 | 2.987e-16 | -7.736e-04 |
| 2.0 | 1.494e-16 | -1.934e-04 |
| 4.0 | 7.468e-17 | -4.835e-05 |
| 8.0 | 3.734e-17 | -1.209e-05 |

dE/dR < 0 everywhere → no stationary point → no equilibrium from scale-free 1/R alone (the B4 obstruction).

## T4: Self-consistent equilibrium

- E(R) = A/R + B·R², A = 1.154e-28, B = 1.002e+09
- R* analytic = (A/2B)^⅓ = 3.8616e-13 m
- R* numeric (scan minimum) = 3.8616e-13 m (rel err 8.0e-06)
- d²E/dR² > 0 (stable minimum): True

## T5: Scale modulus — absolute R* rides on one coupling

| λ | B coupling | R* (m) | predicted λ·R₀ | deviation |
|---:|---:|---:|---:|---:|
| 1.0 | 1.002e+09 | 3.8616e-13 | 3.8616e-13 | 0.0e+00 |
| 2.0 | 1.252e+08 | 7.7232e-13 | 7.7232e-13 | 0.0e+00 |
| 0.5 | 8.013e+09 | 1.9308e-13 | 1.9308e-13 | 0.0e+00 |
| 10.0 | 1.002e+06 | 3.8616e-12 | 3.8616e-12 | 2.1e-16 |

Rescaling the cohesive coupling rescales R* linearly — the one scale modulus (B4-consistent).

## T6: Finite, renormalization-free

- U_EM (throat) = 2.987e-16 J
- rest energy m c² = 8.187e-14 J
- U_EM/(m c²) = 0.003649 = α/2 = 0.003649 (finite, small; no UV divergence)

## T7: Pure-EM relocation → g = 2/α

| R (m) | m c² (J) | g solved |
|---:|---:|---:|
| 1.931e-13 | 1.637e-13 | 274.0720 |
| 3.862e-13 | 8.187e-14 | 274.0720 |
| 7.723e-13 | 4.094e-14 | 274.0720 |

g is R-independent: True; g = 2/α = 274.0720. The balance fixes a geometry–α relation, not the length — the scale relocates to α.

## T8: Assessment

Relocation chain:
  - imposed R_MID (constants.py)
  - ΔR invariant geometric length (#53)
  - finite-self-energy stable equilibrium R*=(A/2B)^(1/3) (this probe)

- self-energy finite: True; equilibrium stable: True
- absolute value derived: False
- remaining prize: pin B (or α-relation) to a second fixed scale (e.g. Planck via closure quantum)

## Verdict

**SELF_CONSISTENT_THROAT_EQUILIBRIUM.** SELF-CONSISTENT THROAT EQUILIBRIUM. The remaining BAM scale anchor — the throat radius R_MID (≈ the invariant bulk separation ΔR) — is recast from an imposed constant into a finite-self-energy stable equilibrium.

FINITE SELF-ENERGY. A point charge has divergent Coulomb self-energy U = α ℏc/(2R) → ∞ as R → 0. The BAM throat removes the divergence geometrically: (short distance) the wormhole throat is the inner boundary (no r < R_MID), capping the field at the throat → U_EM = α ℏc/(2 R_MID) finite; (long distance) on compact S³ the Green function G(ψ) is regular at the antipode with a zero-mean background. The result is U_EM/(m c²) = α/2 ≈ 0.0036 — a finite, small correction with NO UV divergence (no renormalization), unlike QED.

SELF-CONSISTENT EQUILIBRIUM. With only the scale-free EM energy (∝1/R) there is no stationary point — a scale-free 1/R energy cannot fix R (the B4 obstruction in self-energy form). Adding a cohesive term (∝R²) gives E(R)=A/R+B R² with a unique stable minimum R*=(A/2B)^(1/3) (d²E/dR²>0): the throat radius is a stationary point, not imposed.

HONEST CAVEAT (B4). The absolute R* rides on the single dimensionful coupling B: rescaling B → B/λ³ sends R* → λ R* (one scale modulus). Equivalently, the BAM-native balance "rest energy = EM self-energy" is R-independent (both ∝1/R) and fixes a geometry–α relation g = 2/α ≈ 274, not the length — relocating the scale question to α. So the self-consistency recasts the anchor as a finite-self-energy equilibrium condition (genuine progress) and relates it to α, but it does NOT derive the absolute value — by the scale-modulus theorem, that needs one external dimensionful coupling.

Relocation chain: imposed R_MID → ΔR invariant geometric length (#53) → finite-self-energy stable equilibrium (this probe). The remaining prize: pin B (or the α-relation) to a second fixed scale, e.g. a closure-quantum relation to the Planck length.

## What this leaves open

- **The absolute value.** Still one dimensionful coupling (B, or a relation to α). A genuine derivation would pin it to a second fixed scale (e.g. a closure-quantum relation to the Planck length).
- **The cohesive term from first principles.** B·R² is a Poincaré-stress-like cohesion; deriving it from the BAM throat action is the follow-on.
- **Pair-production threshold.** Identified structurally (2 m_e c² = nucleating a throat–antithroat pair at the lowest stable R*); a dynamical nucleation calculation is future work.
