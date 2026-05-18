# Two-mouth flux/action derivation of the Padé factor

**Run:** 2026-05-18T01:38:44+00:00

Derives the caustic Padé factor `K(x) = 2x/(1+x)` from a concrete BAM throat model: two antipodal mouths on S³, closure quantum 2π, equal-action flux continuity.

## Model

- **topology**: two antipodal mouths on S³ connected by throat
- **closure_quantum**: action_base = 2π (BAM closure-ledger)
- **splitting_postulate**: equal-action: ω₁·τ₁ = ω₂·τ₂ = π (flux continuity)
- **derivation_chain**: τ_i = π/ω_i  →  T = π·(ω₁+ω₂)/(ω₁·ω₂)  →  ω_eff = 2·ω₁·ω₂/(ω₁+ω₂) = harmonic mean  →  K = ω_eff/ω₁ = 2x/(1+x)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_BAM_closure_quantum` | action_base = 6.283185 (target 6.283185); source = `repo` | **PASS** |
| T2 | `T2_equal_action_algebraic_derivation` | max |K_derived − 2x/(1+x)| = 2.22e-16 | **PASS** |
| T3 | `T3_numerical_two_segment_orbit` | max |K_numerical − Padé| = 8.82e-12 | **PASS** |
| T4 | `T4_alternative_splittings_rejected` | 1/5 match Padé; equal-action unique = True | **PASS** |
| T5 | `T5_series_impedance_equivalence` | max |2·admittance − H| = 0.00e+00 | **PASS** |
| T6 | `T6_K_squared_reproduces_F2_Padé` | max |K²·Q − F²| = 5.68e-14 | **PASS** |
| T7 | `T7_cross_process_analytic_continuation` | continuation max diff = 0.00e+00 | **PASS** |

## T1: BAM closure quantum

`action_base = 6.2831853072` (target 2π = 6.2831853072), source `repo`.

## T2: Equal-action algebraic derivation

| ω₁ | ω₂ | x | K_derived | K_Padé | diff |
|---:|---:|---:|---:|---:|---:|
| 1.0000 | 0.0100 | 0.0100 | 0.019802 | 0.019802 | 0.00e+00 |
| 1.0000 | 0.1098 | 0.1098 | 0.197873 | 0.197873 | 0.00e+00 |
| 1.0000 | 0.2096 | 0.2096 | 0.346561 | 0.346561 | 0.00e+00 |
| 1.0000 | 0.3094 | 0.3094 | 0.472583 | 0.472583 | 5.55e-17 |
| 1.0000 | 0.4092 | 0.4092 | 0.580755 | 0.580755 | 0.00e+00 |
| 1.0000 | 0.5090 | 0.5090 | 0.674619 | 0.674619 | 0.00e+00 |
| 1.0000 | 0.6088 | 0.6088 | 0.756837 | 0.756837 | 0.00e+00 |
| 1.0000 | 0.7086 | 0.7086 | 0.829451 | 0.829451 | 1.11e-16 |

## T3: Numerical two-segment orbit simulation

| x = ω₂/ω₁ | T (num) | T (analytic) | ω_eff (num) | H(ω₁, ω₂) | K (num) | K (Padé) | diff |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.100 | 34.557519 | 34.557519 | 0.181818 | 0.181818 | 0.181818 | 0.181818 | 5.71e-13 |
| 0.250 | 15.707963 | 15.707963 | 0.400000 | 0.400000 | 0.400000 | 0.400000 | 1.32e-12 |
| 0.500 | 9.424778 | 9.424778 | 0.666667 | 0.666667 | 0.666667 | 0.666667 | 2.30e-12 |
| 0.750 | 7.330383 | 7.330383 | 0.857143 | 0.857143 | 0.857143 | 0.857143 | 4.11e-12 |
| 1.000 | 6.283185 | 6.283185 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 3.16e-12 |
| 1.500 | 5.235988 | 5.235988 | 1.200000 | 1.200000 | 1.200000 | 1.200000 | 8.82e-12 |
| 2.000 | 4.712389 | 4.712389 | 1.333333 | 1.333333 | 1.333333 | 1.333333 | 3.84e-12 |
| 4.000 | 3.926991 | 3.926991 | 1.600000 | 1.600000 | 1.600000 | 1.600000 | 4.24e-12 |
| 10.000 | 3.455752 | 3.455752 | 1.818182 | 1.818182 | 1.818182 | 1.818182 | 4.48e-12 |

## T4: Alternative splittings rejected

| postulate | K values (x = 0.1, 0.5, 1.0, 2.0, 5.0) | max diff from Padé | matches? |
|---|---|---:|:---:|
| `equal_action` | 0.1818, 0.6667, 1.0000, 1.3333, 1.6667 | 0.00e+00 | True |
| `equal_time` | 0.5500, 0.7500, 1.0000, 1.5000, 3.0000 | 1.33e+00 | False |
| `equal_energy_squared` | 0.1089, 0.6000, 1.0000, 1.2000, 1.1538 | 5.13e-01 | False |
| `linear_ratio_swap` | 0.9182, 0.8333, 1.0000, 1.6667, 4.3333 | 2.67e+00 | False |
| `equal_root_action` | 0.3162, 0.7071, 1.0000, 1.4142, 2.2361 | 5.69e-01 | False |

**Equal-action is the unique matching postulate**: True

## T5: Series-impedance equivalence

| ω₁ | ω₂ | Z_total = 1/ω₁ + 1/ω₂ | 1/Z | 2·(1/Z) | H(ω₁, ω₂) | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 1.000 | 0.100 | 11.0000 | 0.0909 | 0.1818 | 0.1818 | 0.00e+00 |
| 1.000 | 0.500 | 3.0000 | 0.3333 | 0.6667 | 0.6667 | 0.00e+00 |
| 1.000 | 1.000 | 2.0000 | 0.5000 | 1.0000 | 1.0000 | 0.00e+00 |
| 1.000 | 2.000 | 1.5000 | 0.6667 | 1.3333 | 1.3333 | 0.00e+00 |
| 1.000 | 10.000 | 1.1000 | 0.9091 | 1.8182 | 1.8182 | 0.00e+00 |

## T6: K² in F² (cross-check with PR #38)

| x | cosθ | K | Q | K²·Q | F² closed | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0500 | -0.900 | 0.09524 | 0.02743 | +2.4881e-04 | +2.4881e-04 | 1.08e-19 |
| 0.0500 | -0.720 | 0.09524 | 0.03222 | +2.9223e-04 | +2.9223e-04 | 1.63e-19 |
| 0.0500 | -0.540 | 0.09524 | 0.03744 | +3.3957e-04 | +3.3957e-04 | 5.42e-20 |
| 0.0500 | -0.360 | 0.09524 | 0.04245 | +3.8501e-04 | +3.8501e-04 | 1.08e-19 |
| 0.0500 | -0.180 | 0.09524 | 0.04621 | +4.1913e-04 | +4.1913e-04 | 1.63e-19 |
| 0.0500 | -0.000 | 0.09524 | 0.04763 | +4.3197e-04 | +4.3197e-04 | 1.08e-19 |

## T7: Cross-process analytic continuation

| β | cosθ_CM | x_⊗ | K(x_⊗) | dist to pole x_⊗=−1 |
|---:|---:|---:|---:|---:|
| 0.100 | -0.700 | -1.1505 | +1.5286e+01 | 0.1505 |
| 0.100 | -0.300 | -1.0619 | +3.4333e+01 | 0.0619 |
| 0.100 | +0.300 | -0.9417 | -3.2333e+01 | 0.0583 |
| 0.100 | +0.700 | -0.8692 | -1.3286e+01 | 0.1308 |
| 0.300 | -0.700 | -1.5316 | +5.7619e+00 | 0.5316 |
| 0.300 | -0.300 | -1.1978 | +1.2111e+01 | 0.1978 |
| 0.300 | +0.300 | -0.8349 | -1.0111e+01 | 0.1651 |
| 0.300 | +0.700 | -0.6529 | -3.7619e+00 | 0.3471 |

**Pole at x_⊗ = -1.0**: Throat-closure breakdown at perpendicular scattering (β·cosθ = 0 with β > 0) in the BW frame — the two mouths face exact antipodal directions and the harmonic-mean rate formally diverges.

## Verdict

**PADÉ_DERIVED.** PADÉ DERIVED. The caustic Padé factor K(x) = 2x/(1+x) is derived from a concrete BAM throat model:
  (1) Two-mouth throat-pair on S³ with mouth frequencies ω₁, ω₂.
  (2) BAM closure quantum action_base = 2π (verified from the closure-ledger).
  (3) Equal-action splitting at the two mouths (flux continuity at the throat bottleneck): ω₁·τ₁ = ω₂·τ₂ = π.
  (4) Effective angular frequency = 2π/T_total = 2·ω₁·ω₂/(ω₁+ω₂) = harmonic mean.
  (5) Normalised to ω₁: K = 2x/(1+x).
Alternative splitting postulates (equal-time, equal-energy², linear-ratio-swap, equal-root-action) all fail to reproduce 2x/(1+x), so the equal-action postulate is the unique flux-continuous choice. A classical series-impedance derivation (Z = 1/ω) converges on the same harmonic mean, confirming the flux-continuity interpretation. The derived K(x) squared reproduces the Padé factor in the F² closed form (PR #38), and the analytic continuation to BW/annihilation kinematics has a pole at x_⊗ = -1 — the throat-closure breakdown locus at perpendicular scattering.

## What this leaves open

- **Why exactly equal action?** The equal-action postulate is *physically* the most natural (flux continuity at a bottleneck), and uniquely picks out the Padé factor — but a deeper BAM derivation tying it to a specific S³ throat action or Lagrangian is the remaining task.
- **Q(x, c) derivation (Hopf-helicity transport channel).** PR #38 identified the Hopf-fibre helicity-transport meaning of the Q factor; this probe only derives the K factor. A complementary "Hopf-fibre helicity transport" probe would close the F² derivation thread.
- **Loop corrections.** Tree-level only.
