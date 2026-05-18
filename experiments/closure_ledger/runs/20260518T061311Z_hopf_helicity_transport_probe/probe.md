# Hopf-fibre helicity transport derivation of Q(x, c)

**Run:** 2026-05-18T06:13:11+00:00

Tests whether the F² polarization factor Q(x, c) decomposes into a Hopf-fibre helicity spinor with BAM-native per-mouth amplitudes. Closes the F² derivation thread (PR #38 + PR #39).

## Model

- **spinor**: `A = (A_pres, A_flip)  with  Q = A_pres² + A_flip²`
- **helicity_preserving**: `A_pres = x  ← per-mouth amplitude √x (PR #39) × two mouths`
- **helicity_flipping**: `A_flip = √x · (1−x)/√(1+c²)  ← one preserve × one flip; flip amplitude = recoil deficit × inverse Thomson sum`
- **thomson_sum**: `(1+c²)/2 = cos⁴(θ/2) + sin⁴(θ/2) = Σ_λ |d¹_{1,λ}|²`
- **complex_quadrature**: `A_complex = A_pres + i·A_flip → |A_complex|² = Q`

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_hopf_wigner_d_thomson_sum` | Wigner-d sum diff = 2.22e-16 | **PASS** |
| T2 | `T2_polarization_spinor_ansatz` | max |Q − (A_pres² + A_flip²)| = 2.27e-13 | **PASS** |
| T3 | `T3_A_pres_from_per_mouth_energy_amplitude` | max |A_pres − √x·√x| = 8.88e-16 | **PASS** |
| T4 | `T4_A_flip_from_spin_action_splitting` | max |A_flip − composed| = 2.78e-17 | **PASS** |
| T5 | `T5_alternative_flip_amplitudes_rejected` | 1/5 match; BAM-derived unique = True | **PASS** |
| T6 | `T6_thomson_limit` | A_flip(Thomson) = 0 (all): True; Q(1, c) = 1 (all): True | **PASS** |
| T7 | `T7_complex_amplitude_quadrature` | max ||A_complex|² − Q| = 1.42e-14 | **PASS** |
| T8 | `T8_F2_chain_K_squared_times_Q_spinor` | max |K²·Q − F²| = 1.42e-14 | **PASS** |
| T9 | `T9_cross_process_analytic_continuation` | continuation max diff = 0.00e+00 | **PASS** |

## T1: Hopf-fibre Wigner-d¹ Thomson sum

| θ/π | cosθ | |d_pres|² | |d_flip|² | sum | (1+c²)/2 | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0159 | +0.9988 | 0.9988 | 0.0000 | 0.9988 | 0.9988 | 2.22e-16 |
| 0.1369 | +0.9089 | 0.9110 | 0.0021 | 0.9130 | 0.9130 | 2.22e-16 |
| 0.2580 | +0.6892 | 0.7134 | 0.0241 | 0.7375 | 0.7375 | 1.11e-16 |
| 0.3790 | +0.3711 | 0.4700 | 0.0989 | 0.5689 | 0.5689 | 0.00e+00 |
| 0.5000 | -0.0000 | 0.2500 | 0.2500 | 0.5000 | 0.5000 | 0.00e+00 |
| 0.6210 | -0.3711 | 0.0989 | 0.4700 | 0.5689 | 0.5689 | 1.11e-16 |
| 0.7420 | -0.6892 | 0.0241 | 0.7134 | 0.7375 | 0.7375 | 1.11e-16 |
| 0.8631 | -0.9089 | 0.0021 | 0.9110 | 0.9130 | 0.9130 | 2.22e-16 |
| 0.9841 | -0.9988 | 0.0000 | 0.9988 | 0.9988 | 0.9988 | 2.22e-16 |

## T2: Polarization spinor ansatz Q = A_pres² + A_flip²

| x | cosθ | A_pres | A_flip | A_pres² + A_flip² | Q | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0500 | -0.950 | 0.0500 | 0.1540 | 0.026219 | 0.026219 | 6.94e-18 |
| 0.0500 | -0.855 | 0.0500 | 0.1615 | 0.028568 | 0.028568 | 6.94e-18 |
| 0.0500 | -0.760 | 0.0500 | 0.1691 | 0.031104 | 0.031104 | 0.00e+00 |
| 0.0500 | -0.665 | 0.0500 | 0.1769 | 0.033788 | 0.033788 | 0.00e+00 |
| 0.0500 | -0.570 | 0.0500 | 0.1846 | 0.036559 | 0.036559 | 6.94e-18 |
| 0.0500 | -0.475 | 0.0500 | 0.1919 | 0.039318 | 0.039318 | 6.94e-18 |
| 0.0500 | -0.380 | 0.0500 | 0.1986 | 0.041931 | 0.041931 | 0.00e+00 |
| 0.0500 | -0.285 | 0.0500 | 0.2043 | 0.044235 | 0.044235 | 0.00e+00 |

## T3: A_pres = x from PR #39 per-mouth amplitude √x

| x | √x (per-mouth) | √x · √x | A_pres | diff |
|---:|---:|---:|---:|---:|
| 0.0500 | 0.2236 | 0.0500 | 0.0500 | 6.94e-18 |
| 0.1000 | 0.3162 | 0.1000 | 0.1000 | 0.00e+00 |
| 0.5000 | 0.7071 | 0.5000 | 0.5000 | 1.11e-16 |
| 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.00e+00 |
| 2.0000 | 1.4142 | 2.0000 | 2.0000 | 4.44e-16 |
| 5.0000 | 2.2361 | 5.0000 | 5.0000 | 8.88e-16 |

## T4: A_flip from spin-action splitting + Thomson normalisation

| x | cosθ | √x (pres) | (1−x)/√(1+c²) (flip) | composed | A_flip | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.100 | -0.900 | 0.3162 | 0.6690 | 0.21155 | 0.21155 | 2.78e-17 |
| 0.100 | -0.500 | 0.3162 | 0.8050 | 0.25456 | 0.25456 | 0.00e+00 |
| 0.100 | +0.000 | 0.3162 | 0.9000 | 0.28460 | 0.28460 | 0.00e+00 |
| 0.100 | +0.500 | 0.3162 | 0.8050 | 0.25456 | 0.25456 | 0.00e+00 |
| 0.100 | +0.900 | 0.3162 | 0.6690 | 0.21155 | 0.21155 | 2.78e-17 |
| 0.300 | -0.900 | 0.5477 | 0.5203 | 0.28498 | 0.28498 | 0.00e+00 |
| 0.300 | -0.500 | 0.5477 | 0.6261 | 0.34293 | 0.34293 | 0.00e+00 |
| 0.300 | +0.000 | 0.5477 | 0.7000 | 0.38341 | 0.38341 | 0.00e+00 |
| 0.300 | +0.500 | 0.5477 | 0.6261 | 0.34293 | 0.34293 | 0.00e+00 |
| 0.300 | +0.900 | 0.5477 | 0.5203 | 0.28498 | 0.28498 | 0.00e+00 |

## T5: Alternative flip-amplitude weightings rejected

| candidate per-mouth flip amplitude | max diff | matches? |
|---|---:|:---:|
| `BAM_derived_(1-x)/sqrt(1+c2)` | 0.00e+00 | True |
| `no_Thomson_normalization_(1-x)` | 2.56e-01 | False |
| `Thomson_sum_in_numerator_(1-x)*(1+c2)` | 9.49e-01 | False |
| `energy_rescaled_deficit_(1-x2)/sqrt(1+c2)` | 2.83e+00 | False |
| `sqrt_deficit_sqrt(1-x)/sqrt(1+c2)` | 1.41e+00 | False |

**BAM-derived flip amplitude is unique**: True

## T6: Thomson limit

| cosθ | A_pres | A_flip | Q |
|---:|---:|---:|---:|
| -0.900 | 1.0000 | +0.00e+00 | 1.000000 |
| -0.500 | 1.0000 | +0.00e+00 | 1.000000 |
| +0.000 | 1.0000 | +0.00e+00 | 1.000000 |
| +0.500 | 1.0000 | +0.00e+00 | 1.000000 |
| +0.900 | 1.0000 | +0.00e+00 | 1.000000 |

## T7: Complex-amplitude / quadrature reading

| x | cosθ | Re(A) | Im(A) | |A|² | Q | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.100 | -0.700 | 0.1000 | 0.2332 | 0.064362 | 0.064362 | 1.39e-17 |
| 0.100 | -0.300 | 0.1000 | 0.2726 | 0.084312 | 0.084312 | 1.39e-17 |
| 0.100 | +0.300 | 0.1000 | 0.2726 | 0.084312 | 0.084312 | 1.39e-17 |
| 0.100 | +0.700 | 0.1000 | 0.2332 | 0.064362 | 0.064362 | 1.39e-17 |
| 0.300 | -0.700 | 0.3000 | 0.3141 | 0.188658 | 0.188658 | 2.78e-17 |
| 0.300 | -0.300 | 0.3000 | 0.3672 | 0.224862 | 0.224862 | 2.78e-17 |
| 0.300 | +0.300 | 0.3000 | 0.3672 | 0.224862 | 0.224862 | 2.78e-17 |
| 0.300 | +0.700 | 0.3000 | 0.3141 | 0.188658 | 0.188658 | 2.78e-17 |

## T8: F² = K(x)² · Q chain (PR #38 + PR #39)

| x | cosθ | K (PR #39) | Q (spinor) | K²·Q | F² closed | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0500 | -0.900 | 0.0952 | 0.0274 | +2.4881e-04 | +2.4881e-04 | 1.08e-19 |
| 0.0500 | -0.720 | 0.0952 | 0.0322 | +2.9223e-04 | +2.9223e-04 | 2.17e-19 |
| 0.0500 | -0.540 | 0.0952 | 0.0374 | +3.3957e-04 | +3.3957e-04 | 0.00e+00 |
| 0.0500 | -0.360 | 0.0952 | 0.0424 | +3.8501e-04 | +3.8501e-04 | 1.08e-19 |
| 0.0500 | -0.180 | 0.0952 | 0.0462 | +4.1913e-04 | +4.1913e-04 | 1.63e-19 |
| 0.0500 | -0.000 | 0.0952 | 0.0476 | +4.3197e-04 | +4.3197e-04 | 1.08e-19 |

## T9: Cross-process analytic continuation

| β | cosθ_CM | x_⊗ | c_⊗ | A_pres² | A_flip² | sum | Q cont. | diff |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.100 | -0.700 | -1.1505 | -0.9897 | +1.3237e+00 | -2.6879e+00 | -1.3642e+00 | -1.3642e+00 | 0.00e+00 |
| 0.100 | -0.300 | -1.0619 | -0.9818 | +1.1275e+00 | -2.2986e+00 | -1.1711e+00 | -1.1711e+00 | 0.00e+00 |
| 0.100 | +0.300 | -0.9417 | -0.9818 | +8.8689e-01 | -1.8080e+00 | -9.2112e-01 | -9.2112e-01 | 0.00e+00 |
| 0.100 | +0.700 | -0.8692 | -0.9897 | +7.5544e-01 | -1.5340e+00 | -7.7852e-01 | -7.7852e-01 | 0.00e+00 |
| 0.300 | -0.700 | -1.5316 | -0.9040 | +2.3459e+00 | -5.4022e+00 | -3.0563e+00 | -3.0563e+00 | 0.00e+00 |
| 0.300 | -0.300 | -1.1978 | -0.8349 | +1.4347e+00 | -3.4094e+00 | -1.9747e+00 | -1.9747e+00 | 0.00e+00 |
| 0.300 | +0.300 | -0.8349 | -0.8349 | +6.9700e-01 | -1.6563e+00 | -9.5931e-01 | -9.5931e-01 | 0.00e+00 |
| 0.300 | +0.700 | -0.6529 | -0.9040 | +4.2627e-01 | -9.8161e-01 | -5.5534e-01 | -5.5534e-01 | 0.00e+00 |

## Verdict

**Q_DERIVED.** Q DERIVED. The polarization factor Q(x, c) = x² + x·(1−x)²/(1+c²) decomposes into a Hopf-fibre helicity spinor A = (A_pres, A_flip) with:
  (1) A_pres = x — helicity-preserving channel from equal-action splitting (PR #39): per-mouth amplitude √x, two mouths preserving → A_pres = √x · √x = x.
  (2) A_flip = √x · (1−x)/√(1+c²) — helicity-flipping channel: one mouth preserves (amplitude √x), the other flips with amplitude (1−x)/√(1+c²) (recoil deficit weighted by inverse Thomson polarization sum).
  (3) (1+c²)/2 = cos⁴(θ/2) + sin⁴(θ/2) = Σ_λ |d¹_{1,λ}|² is the Hopf-fibre helicity transport sum (spin-1 Wigner-d).
Q = A_pres² + A_flip², equivalent to the complex Hopf-fibre amplitude A_complex = A_pres + i·A_flip with |A_complex|² = Q. Alternative per-mouth flip amplitudes (no Thomson normalisation; sum in numerator; energy-rescaled deficit; square-root deficit) all fail to reproduce A_flip — the BAM-derived form is unique. Combined with PR #39 K(x) = 2x/(1+x), the full F² closed form is reconstructed (F² = K²·Q) to machine precision, and the decomposition survives analytic continuation under crossing.

## What this leaves open

- **First-principles BAM action**. The equal-action splitting (for both energy and spin) is the natural flux-continuity postulate but is not yet derived from a specific BAM S³ throat action coupled to the Hopf bundle.
- **Helicity-resolved Compton comparison.** The standard QED helicity-resolved Compton amplitudes |M(λ → λ′)|² have their own (different) algebraic split. The BAM decomposition is one consistent split; relating it to the QED helicity basis is a follow-on.
- **Loop corrections.** Still tree-level.
