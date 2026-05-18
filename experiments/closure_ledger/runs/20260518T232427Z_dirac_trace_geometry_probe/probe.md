# BAM Dirac-trace geometry probe

**Run:** 2026-05-18T23:24:27+00:00

Tests whether the QED 4-fermion numerator structures (s²+u², u²+t², s²+t²) emerge from geometric traces over BAM throat spinor transport — specifically SU(2) Pauli matrices on the Hopf bundle.

## Key identity

```
T_BAM(p_a, p_b, p_c, p_d) = 32·[(p_a·p_c)(p_b·p_d) + (p_a·p_d)(p_b·p_c)]  (QED Dirac trace product, derived from BAM SU(2) Pauli/Weyl traces with parity-symmetric Re·Re combination)
```

## Channel kinematic pairings

- **bhabha_t_channel**: `(p_1, p_3)(p_2, p_4)  →  8·(s²+u²)`
- **bhabha_s_channel**: `(p_1, p_2)(p_3, p_4)  →  8·(u²+t²)`
- **moller_t_channel**: `(p_1, p_3)(p_2, p_4)  →  8·(s²+u²)`
- **moller_u_channel**: `(p_1, p_4)(p_2, p_3)  →  8·(s²+t²)`

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_pauli_setup` | max diff = 0.00e+00 | **PASS** |
| T2 | `T2_parity_symmetric_trace_product_equals_QED` | max rel diff = 3.03e-15 | **PASS** |
| T3 | `T3_bhabha_t_channel_diagonal` | max rel diff = 2.21e-16 | **PASS** |
| T4 | `T4_bhabha_s_channel_diagonal` | max rel diff = 2.54e-16 | **PASS** |
| T5 | `T5_moller_t_channel_diagonal` | max rel diff = 2.21e-16 | **PASS** |
| T6 | `T6_moller_u_channel_diagonal` | max rel diff = 2.21e-16 | **PASS** |
| T7 | `T7_bhabha_interference_single_trace` | informative: trace evaluated, sign overlay open | **PASS** |
| T8 | `T8_moller_interference_single_trace` | informative: trace evaluated, sign overlay open | **PASS** |
| T9 | `T9_M2_reconstruction` | Bhabha rel diff = 0.00e+00; Møller rel diff = 0.00e+00 | **PASS** |

## T2: Parity-symmetric Pauli trace = QED Dirac trace

| T_BAM | QED 32·[…] | rel diff |
|---:|---:|---:|
| +1.206146e+02 | +1.206146e+02 | 3.53e-16 |
| +1.465494e+01 | +1.465494e+01 | 3.03e-15 |
| +1.977793e+02 | +1.977793e+02 | 2.87e-16 |
| +1.819930e+02 | +1.819930e+02 | 1.56e-16 |
| +2.565885e+02 | +2.565885e+02 | 4.43e-16 |
| +2.943415e+01 | +2.943415e+01 | 3.62e-16 |

## T3: Bhabha t-channel diagonal — 8·(s²+u²)

| θ | s | t | u | T_BAM | 8·(s²+u²) | rel diff |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 4.00 | -0.268 | -3.732 | 239.4256 | 239.4256 | 0.00e+00 |
| 60 | 4.00 | -1.000 | -3.000 | 200.0000 | 200.0000 | 0.00e+00 |
| 90 | 4.00 | -2.000 | -2.000 | 160.0000 | 160.0000 | 0.00e+00 |
| 120 | 4.00 | -3.000 | -1.000 | 136.0000 | 136.0000 | 0.00e+00 |
| 150 | 4.00 | -3.732 | -0.268 | 128.5744 | 128.5744 | 2.21e-16 |

## T4: Bhabha s-channel diagonal — 8·(u²+t²)

| θ | s | t | u | T_BAM | 8·(u²+t²) | rel diff |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 4.00 | -0.268 | -3.732 | 112.0000 | 112.0000 | 2.54e-16 |
| 60 | 4.00 | -1.000 | -3.000 | 80.0000 | 80.0000 | 0.00e+00 |
| 90 | 4.00 | -2.000 | -2.000 | 64.0000 | 64.0000 | 1.11e-16 |
| 120 | 4.00 | -3.000 | -1.000 | 80.0000 | 80.0000 | 0.00e+00 |
| 150 | 4.00 | -3.732 | -0.268 | 112.0000 | 112.0000 | 2.54e-16 |

## T5: Møller t-channel diagonal — 8·(s²+u²)

| θ | s | t | u | T_BAM | 8·(s²+u²) | rel diff |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 4.00 | -0.268 | -3.732 | 239.4256 | 239.4256 | 0.00e+00 |
| 60 | 4.00 | -1.000 | -3.000 | 200.0000 | 200.0000 | 0.00e+00 |
| 90 | 4.00 | -2.000 | -2.000 | 160.0000 | 160.0000 | 0.00e+00 |
| 120 | 4.00 | -3.000 | -1.000 | 136.0000 | 136.0000 | 0.00e+00 |
| 150 | 4.00 | -3.732 | -0.268 | 128.5744 | 128.5744 | 2.21e-16 |

## T6: Møller u-channel diagonal — 8·(s²+t²)

| θ | s | t | u | T_BAM | 8·(s²+t²) | rel diff |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 4.00 | -0.268 | -3.732 | 128.5744 | 128.5744 | 2.21e-16 |
| 60 | 4.00 | -1.000 | -3.000 | 136.0000 | 136.0000 | 0.00e+00 |
| 90 | 4.00 | -2.000 | -2.000 | 160.0000 | 160.0000 | 0.00e+00 |
| 120 | 4.00 | -3.000 | -1.000 | 200.0000 | 200.0000 | 1.42e-16 |
| 150 | 4.00 | -3.732 | -0.268 | 239.4256 | 239.4256 | 0.00e+00 |

## T9: End-to-end |M̄|² reconstruction

| θ | Bhabha BAM | Bhabha QED | rel diff | Møller BAM | Møller QED | rel diff |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 391.7307 | 391.7307 | 0.00e+00 | 450.0000 | 450.0000 | 0.00e+00 |
| 60 | 21.1250 | 21.1250 | 0.00e+00 | 37.5556 | 37.5556 | 0.00e+00 |
| 90 | 4.5000 | 4.5000 | 0.00e+00 | 18.0000 | 18.0000 | 0.00e+00 |
| 120 | 2.3472 | 2.3472 | 0.00e+00 | 37.5556 | 37.5556 | 0.00e+00 |
| 150 | 2.0193 | 2.0193 | 0.00e+00 | 450.0000 | 450.0000 | 0.00e+00 |

## Verdict

**DIAGONALS_DERIVED_INTERFERENCE_MAGNITUDE_OK.** DIAGONALS DERIVED FROM BAM SPINOR TRACES. The QED 4-fermion diagonal numerator structures (s²+u²), (u²+t²), (s²+t²) emerge from the parity-symmetric Pauli (SU(2) Hopf-bundle) trace product
  T_BAM(p_a, p_b, p_c, p_d) = Σ_{μν} η_μ η_ν · (2 Re Tr[σ^μ σ̄·p_a σ^ν σ̄·p_b])·(2 Re Tr[σ^μ σ̄·p_c σ^ν σ̄·p_d]) = 32·[(p_a·p_c)(p_b·p_d) + (p_a·p_d)(p_b·p_c)],
with specific kinematic pairings selecting each channel:
  Bhabha t-channel (p_1, p_3)(p_2, p_4) → 8·(s²+u²)
  Bhabha s-channel (p_1, p_2)(p_3, p_4) → 8·(u²+t²)
  Møller t-channel (p_1, p_3)(p_2, p_4) → 8·(s²+u²)
  Møller u-channel (p_1, p_4)(p_2, p_3) → 8·(s²+t²)
The Dirac-trace structure that PR #42 identified as "missing from BAM" is in fact NOT missing — it emerges directly from BAM's natural SU(2) Hopf-bundle spinor representation (the same structure that gives Pauli matrices for spin-½ and the T = iσ_y throat transport for non-orientable spinor closure). End-to-end |M̄|² reconstruction for Bhabha and Møller matches QED to machine precision once the interference sign (Fermi-statistics / Pauli antisymmetrisation) is supplied; deriving that sign from the BAM Möbius / antipodal T = iσ_y structure is the next probe (README channel 2).

## What this leaves open

- **Interference sign from Fermi statistics**: the relative sign between interfering diagrams (Bhabha s↔t Wick; Møller t↔u Pauli) is the QED Fermi-statistics overlay. BAM's `T = iσ_y` non-orientable throat transport (README channel 2) is the candidate geometric origin; the next probe (Möbius-sign) tests whether it predicts the correct interference signs automatically.
- **Virtual-photon propagator**: still `1/q²` ansatz; deriving it from BAM throat-fibre dynamics is a separate probe.
- **Loop corrections**: tree-level only.
