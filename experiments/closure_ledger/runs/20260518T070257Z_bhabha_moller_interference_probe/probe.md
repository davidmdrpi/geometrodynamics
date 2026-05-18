# Bhabha / Møller two-channel interference probe

**Run:** 2026-05-18T07:02:57+00:00

Tests whether the BAM scalar-intensity tree kernel can reproduce tree-level QED Bhabha/Møller, where two crossed kernels must be coherently summed with a Fermi-statistics-determined relative sign.

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_QED_bhabha_scalar_intensity` | QED Bhabha intensity tabulated at 5 CM angles | **PASS** |
| T2 | `T2_QED_moller_scalar_intensity` | QED Møller intensity tabulated at 5 CM angles | **PASS** |
| T3 | `T3_interference_signs` | Bhabha all-negative: True; Møller all-positive: True | **PASS** |
| T4 | `T4_construction_I_positive_sum` | Bhabha max residual = 6.419e+01; Møller max residual = 1.186e+01 | *falsified (documented)* |
| T5 | `T5_construction_II_fermi_sign` | Bhabha max residual = 1.221e+01; Møller max residual = 1.186e+01 | *falsified (documented)* |
| T6 | `T6_construction_III_compton_vertex` | Bhabha diag residual = 9.000e+00; Møller diag residual = 9.000e+00 | *falsified (documented)* |
| T7 | `T7_residual_analysis` | residuals identified: Dirac numerators + interference magnitude | **PASS** |
| T8 | `T8_required_BAM_extension` | 3 missing ingredients; 3 next-probe candidates | **PASS** |

## T1: QED Bhabha scalar intensity

| θ (deg) | s | t | u | M²_t | M²_s | interference | total |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 30 | 1.000 | -0.067 | -0.933 | 416.8461 | 0.8750 | -25.9904 | 391.7307 |
| 60 | 1.000 | -0.250 | -0.750 | 25.0000 | 0.6250 | -4.5000 | 21.1250 |
| 90 | 1.000 | -0.500 | -0.500 | 5.0000 | 0.5000 | -1.0000 | 4.5000 |
| 120 | 1.000 | -0.750 | -0.250 | 1.8889 | 0.6250 | -0.1667 | 2.3472 |
| 150 | 1.000 | -0.933 | -0.067 | 1.1539 | 0.8750 | -0.0096 | 2.0193 |

## T2: QED Møller scalar intensity

| θ (deg) | s | t | u | M²_t | M²_u | interference | total |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 30 | 1.000 | -0.067 | -0.933 | 416.8461 | 1.1539 | +32.0000 | 450.0000 |
| 60 | 1.000 | -0.250 | -0.750 | 25.0000 | 1.8889 | +10.6667 | 37.5556 |
| 90 | 1.000 | -0.500 | -0.500 | 5.0000 | 5.0000 | +8.0000 | 18.0000 |
| 120 | 1.000 | -0.750 | -0.250 | 1.8889 | 25.0000 | +10.6667 | 37.5556 |
| 150 | 1.000 | -0.933 | -0.067 | 1.1539 | 416.8461 | +32.0000 | 450.0000 |

## T3: Interference signs

Bhabha signs across angles: `[np.float64(-1.0), np.float64(-1.0), np.float64(-1.0), np.float64(-1.0), np.float64(-1.0)]`
Møller signs across angles: `[np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0)]`
Bhabha all negative: **True**; Møller all positive: **True**.

## T4: BAM Construction I (positive √-channel sum)

| θ | Bhabha QED interf | Bhabha BAM+ interf | Bhabha total resid | Møller QED interf | Møller BAM+ interf | Møller total resid |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | -25.9904 | +38.1964 | +64.1867 | +32.0000 | +43.8634 | +11.8634 |
| 60 | -4.5000 | +7.9057 | +12.4057 | +10.6667 | +13.7437 | +3.0770 |
| 90 | -1.0000 | +3.1623 | +4.1623 | +8.0000 | +10.0000 | +2.0000 |
| 120 | -0.1667 | +2.1731 | +2.3397 | +10.6667 | +13.7437 | +3.0770 |
| 150 | -0.0096 | +2.0096 | +2.0193 | +32.0000 | +43.8634 | +11.8634 |

## T5: BAM Construction II (Fermi sign by hand)

| θ | Bhabha QED interf | Bhabha BAM± interf | Bhabha total resid | Møller QED interf | Møller BAM± interf | Møller total resid |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | -25.9904 | -38.1964 | -12.2060 | +32.0000 | +43.8634 | +11.8634 |
| 60 | -4.5000 | -7.9057 | -3.4057 | +10.6667 | +13.7437 | +3.0770 |
| 90 | -1.0000 | -3.1623 | -2.1623 | +8.0000 | +10.0000 | +2.0000 |
| 120 | -0.1667 | -2.1731 | -2.0064 | +10.6667 | +13.7437 | +3.0770 |
| 150 | -0.0096 | -2.0096 | -2.0000 | +32.0000 | +43.8634 | +11.8634 |

## T6: BAM Construction III (Compton-F² vertex factor)

| θ | Bhabha QED M²_t | BAM M²_t | resid | Bhabha QED M²_s | BAM M²_s | resid |
|---:|---:|---:|---:|---:|---:|---:|
| 60 | 25.0000 | 16.0000 | -9.0000 | 0.6250 | 1.0000 | +0.3750 |
| 90 | 5.0000 | 4.0000 | -1.0000 | 0.5000 | 1.0000 | +0.5000 |
| 120 | 1.8889 | 1.7778 | -0.1111 | 0.6250 | 1.0000 | +0.3750 |

## T7: Residual analysis

### Bhabha interference vs geometric-mean (BAM ansatz)

| θ | |QED interference| | 2·√(M²_t · M²_s) BAM | ratio |
|---:|---:|---:|---:|
| 60 | 4.5000 | 7.9057 | 1.757 |
| 90 | 1.0000 | 3.1623 | 3.162 |
| 120 | 0.1667 | 2.1731 | 13.038 |

## T8: Required BAM extension

### Missing ingredients
1. Dirac spinors at the vertex. The Compton BAM kernel F² = K²·Q already integrates over Dirac structure implicitly via the closed form, but for multi-channel processes the per-channel amplitudes must carry explicit spinor weights.
2. Virtual-photon propagator structure beyond 1/q². Internal photon lines couple two BAM vertices; their transport on the throat may carry Hopf-fibre phase in addition to 1/q².
3. Channel sign rule from non-orientable throat topology. The README's channel-2 (Möbius / antipodal T = iσ_y) is the geometric origin of spin-½ and could in principle drive Pauli signs for identical-fermion interference (Møller t↔u). Whether it correctly predicts the Bhabha s↔t Wick sign is a deeper open question.

### Next-probe candidates
1. Dirac-trace BAM probe: derive (s²+u²) etc. from a geometric trace over per-vertex spinor configurations.
2. Möbius-sign probe: test whether T = iσ_y gives the Møller Pauli sign automatically.
3. Virtual-photon throat-fibre propagator probe: derive 1/q² + corrections from BAM throat-fibre dynamics.

## Verdict

**DIAGONAL_FAILS_TOO.** DIAGONAL_FAILS_TOO. Even the diagonal channel intensities for Bhabha/Møller do not match the BAM Compton-F² analytic continuation — these 4-fermion processes sit outside the Compton kernel's crossing-orbit class. A fundamentally new BAM tree kernel is needed for multi-fermion scattering, with explicit Dirac-trace structure at the vertex and a virtual-photon propagator beyond 1/q². The Compton/BW/ann triangle does NOT extend to Bhabha/Møller by Mandelstam crossing alone.

## What this leaves open

- **Dirac-trace structure in BAM**: capturing the QED γ-matrix algebra that produces (s²+u²), (u²+t²), (s²+t²) numerators requires either explicit Dirac spinors at vertices or a geometric proxy.
- **Möbius/antipodal Pauli signs**: BAM's `T = iσ_y` antipodal transport (README channel 2) is a candidate origin for the Møller t↔u Pauli sign. Whether it gives the right sign automatically is a separate probe target.
- **Virtual-photon throat-fibre propagator**: 1/q² alone does not reproduce QED; additional Hopf-fibre phase structure may be needed.
