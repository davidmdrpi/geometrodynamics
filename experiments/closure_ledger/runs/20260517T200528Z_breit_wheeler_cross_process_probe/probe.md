# Breit-Wheeler cross-process validation probe

**Run:** 2026-05-17T20:05:28+00:00

Cross-process test of the closed-form BAM Compton vertex factor F²(x, c) under standard QED Mandelstam crossing to pair production γγ → e⁺e⁻.

**Crossing map:**

```
s_C → u_BW, t_C → s_BW, u_C → t_BW
```

**Crossed variables at BW kinematics (β, c = cosθ_CM):**

```
x_⊗ = −(1−βc)/(1+βc),  c_⊗ = (2β²−β²c²−1)/(1−β²c²),  sin²θ_⊗ = 4β²(1−β²)sin²θ/(1−β²c²)²
```

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_mandelstam_crossing_identity` | max |M̄|²_BW − (−|M̄|²_KN_crossed) = 5.68e-14 | **PASS** |
| T2 | `T2_BAM_crossed_reproduces_BW` | max |M̄|²_BW_BAM − |M̄|²_BW_textbook = 2.13e-14 | **PASS** |
| T3 | `T3_individual_factor_inspection` | baseline range = [0.000e+00, 3.950e+01]; F² range = [-inf, -8.260e-04]; 0 baseline-negative samples; 78 F²-negative samples | **PASS** |
| T4 | `T4_total_cross_section_recovery` | ⟨ratio⟩ = 1.000000, rel-σ = 2.14e-16 | **PASS** |
| T5 | `T5_threshold_and_ultrarelativistic_limits` | ⟨BAM/textbook⟩ = 1.000000, rel-spread = 1.77e-14 | **PASS** |

## T1: Mandelstam crossing identity

Standard QED gives |M̄|²_BW(β, c) = −|M̄|²_KN(s,t,u) evaluated at the crossed Mandelstam triple. The minus sign comes from one crossed fermion line.

| β | cosθ | M²_KN (crossed) | M²_BW (direct) | diff |
|---:|---:|---:|---:|---:|
| 0.050 | -0.950 | -2.0100e+00 | +2.0100e+00 | 4.44e-16 |
| 0.050 | -0.831 | -2.0100e+00 | +2.0100e+00 | 0.00e+00 |
| 0.050 | -0.712 | -2.0100e+00 | +2.0100e+00 | 1.78e-15 |
| 0.050 | -0.594 | -2.0100e+00 | +2.0100e+00 | 8.88e-16 |
| 0.050 | -0.475 | -2.0100e+00 | +2.0100e+00 | 0.00e+00 |
| 0.050 | -0.356 | -2.0100e+00 | +2.0100e+00 | 0.00e+00 |
| 0.050 | -0.238 | -2.0100e+00 | +2.0100e+00 | 1.78e-15 |
| 0.050 | -0.119 | -2.0100e+00 | +2.0100e+00 | 1.33e-15 |

## T2: BAM crossed prediction vs textbook BW

BAM-predicted |M̄|²_BW = −2·(f_baseline · F²)/x_⊗² at the crossed variables vs textbook |M̄|²_BW. The conversion factor −2/x_⊗² is the standard Compton x²/2 relation between f_KN_invariant and |M̄|²/(8e⁴), with the fermion-crossing sign.

| β | cosθ | x_⊗ | c_⊗ | BAM M² | textbook M² | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.100 | -0.950 | -1.2099 | -0.9980 | +2.0404e+00 | +2.0404e+00 | 4.44e-16 |
| 0.100 | -0.831 | -1.1813 | -0.9938 | +2.0402e+00 | +2.0402e+00 | 0.00e+00 |
| 0.100 | -0.712 | -1.1534 | -0.9901 | +2.0401e+00 | +2.0401e+00 | 4.44e-16 |
| 0.100 | -0.594 | -1.1262 | -0.9870 | +2.0400e+00 | +2.0400e+00 | 1.78e-15 |
| 0.100 | -0.475 | -1.0997 | -0.9845 | +2.0398e+00 | +2.0398e+00 | 4.44e-16 |
| 0.100 | -0.356 | -1.0739 | -0.9825 | +2.0397e+00 | +2.0397e+00 | 1.78e-15 |
| 0.100 | -0.238 | -1.0487 | -0.9811 | +2.0397e+00 | +2.0397e+00 | 4.44e-16 |
| 0.100 | -0.119 | -1.0240 | -0.9803 | +2.0396e+00 | +2.0396e+00 | 8.88e-16 |
| 0.100 | +0.000 | -1.0000 | -0.9800 | +2.0396e+00 | +2.0396e+00 | 4.44e-16 |
| 0.100 | +0.119 | -0.9765 | -0.9803 | +2.0396e+00 | +2.0396e+00 | 4.44e-16 |
| 0.100 | +0.238 | -0.9536 | -0.9811 | +2.0397e+00 | +2.0397e+00 | 4.44e-16 |
| 0.100 | +0.356 | -0.9312 | -0.9825 | +2.0397e+00 | +2.0397e+00 | 1.33e-15 |

Max |BAM − textbook| = 2.1316e-14

## T3: Individual factor inspection at BW kinematics

f_baseline_C(x_⊗, c_⊗) range: [0.000e+00, 3.950e+01]
F²_C(x_⊗, c_⊗) range: [-inf, -8.260e-04]
Baseline-negative samples: 0; F²-negative samples: 78.

| β | cosθ | x_⊗ | c_⊗ | baseline | F² |
|---:|---:|---:|---:|---:|---:|
| 0.100 | -0.900 | -1.1978 | -0.9962 | +6.7915e-03 | -2.1551e+02 |
| 0.100 | -0.750 | -1.1622 | -0.9912 | +4.8249e-03 | -2.8555e+02 |
| 0.100 | -0.600 | -1.1277 | -0.9872 | +3.1631e-03 | -4.1005e+02 |
| 0.100 | -0.450 | -1.0942 | -0.9840 | +1.8250e-03 | -6.6917e+02 |
| 0.100 | -0.300 | -1.0619 | -0.9818 | +8.3302e-04 | -1.3804e+03 |
| 0.100 | -0.150 | -1.0305 | -0.9804 | +2.1417e-04 | -5.0562e+03 |
| 0.100 | -0.000 | -1.0000 | -0.9800 | +0.0000e+00 | -inf |
| 0.100 | +0.150 | -0.9704 | -0.9804 | +2.2741e-04 | -4.2232e+03 |
| 0.100 | +0.300 | -0.9417 | -0.9818 | +9.3927e-04 | -9.6298e+02 |
| 0.100 | +0.450 | -0.9139 | -0.9840 | +2.1851e-03 | -3.8982e+02 |

## T4: Total BW cross section

| β | textbook ∫ | BAM ∫ | ratio | textbook σ_BW |
|---:|---:|---:|---:|---:|
| 0.100 | +1.0098e-01 | +1.0098e-01 | +1.000000 | +3.1723e-01 |
| 0.300 | +3.2128e-01 | +3.2128e-01 | +1.000000 | +1.0093e+00 |
| 0.500 | +5.5394e-01 | +5.5394e-01 | +1.000000 | +1.7403e+00 |
| 0.700 | +6.8170e-01 | +6.8170e-01 | +1.000000 | +2.1416e+00 |
| 0.900 | +4.5217e-01 | +4.5217e-01 | +1.000000 | +1.4205e+00 |
| 0.950 | +2.8873e-01 | +2.8873e-01 | +1.000000 | +9.0688e-01 |
| 0.980 | +1.4887e-01 | +1.4887e-01 | +1.000000 | +4.6711e-01 |

## T5: Threshold and ultra-relativistic limits

### Threshold (β → 0): linear suppression

| β | textbook ∫ | BAM ∫ | textbook/β | BAM/β | ratio |
|---:|---:|---:|---:|---:|---:|
| 0.010 | +1.0001e-02 | +1.0001e-02 | +1.0001e+00 | +1.0001e+00 | +1.00000000 |
| 0.020 | +2.0008e-02 | +2.0008e-02 | +1.0004e+00 | +1.0004e+00 | +1.00000000 |
| 0.050 | +5.0124e-02 | +5.0124e-02 | +1.0025e+00 | +1.0025e+00 | +1.00000000 |
| 0.100 | +1.0098e-01 | +1.0098e-01 | +1.0098e+00 | +1.0098e+00 | +1.00000000 |

### Ultra-relativistic (β → 1): logarithmic falloff

| β | textbook ∫ | BAM ∫ | log env | textbook/env | BAM/env | ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 0.9500 | +2.8873e-01 | +2.8873e-01 | +3.5720e-01 | +8.0833e-01 | +8.0833e-01 | +1.00000000 |
| 0.9900 | +8.7707e-02 | +8.7707e-02 | +1.0534e-01 | +8.3263e-01 | +8.3263e-01 | +1.00000000 |
| 0.9990 | +1.5569e-02 | +1.5569e-02 | +1.5193e-02 | +1.0248e+00 | +1.0248e+00 | +1.00000000 |

Overall ratio statistics: ⟨BAM/textbook⟩ = 1.00000000, relative spread = 1.77e-14

## Verdict

**PROCESS_GENERAL_UNDER_CROSSING.** PROCESS-GENERAL UNDER CROSSING. The BAM Compton vertex factor F²(x, c) = 4·x³·(x²+1−x·sin²θ)/[(1+c²)(1+x)²], when expressed in Lorentz-invariant form and analytically continued via standard Mandelstam crossing (s_C → u_BW, t_C → s_BW, u_C → t_BW), exactly reproduces the Breit-Wheeler differential and total cross sections. The Compton-side construction is not a process-specific algebraic fit — it is the closed form of the invariant QED amplitude in Compton variables, and crossing carries it to BW automatically. The same closed-form F (with x_⊗ < 0 as the analytic continuation of ω′/ω) governs pair production.

## What this leaves open

- **Pair annihilation e⁺e⁻ → γγ.** A second crossed partner with the photons outgoing rather than incoming; should follow trivially from BW by reversing the crossing.
- **Bhabha / Møller scattering.** Both s- and t-channel diagrams together; the BAM construction would need to combine two crossed copies of the elementary vertex factor.
- **Loop corrections.** The closed-form F is tree-level. Whether the same algebraic structure persists at one loop is a separate (large) thread target.
- **BAM derivation of F from first principles.** Still open from PR #35 — the closed form is identified but not derived from a BAM Lagrangian / closure action.
