# Pair-annihilation crossing probe

**Run:** 2026-05-17T20:46:57+00:00

Closes the Compton/Breit-Wheeler/annihilation crossing triangle. PR #35 established the closed-form F²(x, c); PR #36 verified Compton → BW; this probe verifies Compton → annihilation and the full triangle loop closure.

**Crossing triangle:**

```
       Compton (γe → γe)
         /            \
        v              v
    BW (γγ→ee) ◀▶ ann (ee→γγ)
                 (T-reversal)
```

| edge | crossing |
|---|---|
| `C_to_BW` | s_C → u_BW, t_C → s_BW, u_C → t_BW (PR #36) |
| `C_to_ann` | s_C → u_ann, t_C → s_ann, u_C → t_ann (this PR) |
| `BW_to_ann` | T-reversal (same Mandelstam, same |M̄|²) |

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_mandelstam_crossing_ann` | max diff = 5.68e-14 | **PASS** |
| T2 | `T2_BAM_crossed_to_ann` | max |M̄|²_ann_BAM − textbook = 2.13e-14 | **PASS** |
| T3 | `T3_T_invariance_BW_equals_ann` | max |BAM-BW − BAM-ann| = 0.00e+00 | **PASS** |
| T4 | `T4_triangle_loop_closure` | label diff = 0.00e+00, amp diff = 0.00e+00 | **PASS** |
| T5 | `T5_total_cross_section` | ⟨BAM/textbook⟩ = 1.00000000, rel-spread = 4.44e-16 | **PASS** |
| T6 | `T6_threshold_and_ultrarelativistic_limits` | ⟨BAM/textbook⟩ = 1.00000000, rel-spread = 1.80e-14 | **PASS** |

## T1: Mandelstam crossing identity for annihilation

`|M̄|²_ann(β, c) = −|M̄|²_KN(s = u_ann, t = s_ann, u = t_ann)`.

| β | cosθ | M²_KN (crossed) | M²_ann (direct) | diff |
|---:|---:|---:|---:|---:|
| 0.050 | -0.950 | -2.0100e+00 | +2.0100e+00 | 4.44e-16 |
| 0.050 | -0.831 | -2.0100e+00 | +2.0100e+00 | 0.00e+00 |
| 0.050 | -0.712 | -2.0100e+00 | +2.0100e+00 | 1.78e-15 |
| 0.050 | -0.594 | -2.0100e+00 | +2.0100e+00 | 8.88e-16 |
| 0.050 | -0.475 | -2.0100e+00 | +2.0100e+00 | 0.00e+00 |
| 0.050 | -0.356 | -2.0100e+00 | +2.0100e+00 | 0.00e+00 |
| 0.050 | -0.238 | -2.0100e+00 | +2.0100e+00 | 1.78e-15 |
| 0.050 | -0.119 | -2.0100e+00 | +2.0100e+00 | 1.33e-15 |

## T2: BAM kernel crossed to annihilation

`|M̄|²_ann_BAM = −2·(f_baseline · F²)/x_⊗²` at the annihilation-crossed variables.

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

## T3: T-invariance — BAM-BW = BAM-ann

| β | cosθ | x_⊗ (BW) | x_⊗ (ann) | BAM M²_BW | BAM M²_ann | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.100 | -0.950 | -1.2099 | -1.2099 | +2.0404e+00 | +2.0404e+00 | 0.00e+00 |
| 0.100 | -0.831 | -1.1813 | -1.1813 | +2.0402e+00 | +2.0402e+00 | 0.00e+00 |
| 0.100 | -0.712 | -1.1534 | -1.1534 | +2.0401e+00 | +2.0401e+00 | 0.00e+00 |
| 0.100 | -0.594 | -1.1262 | -1.1262 | +2.0400e+00 | +2.0400e+00 | 0.00e+00 |
| 0.100 | -0.475 | -1.0997 | -1.0997 | +2.0398e+00 | +2.0398e+00 | 0.00e+00 |
| 0.100 | -0.356 | -1.0739 | -1.0739 | +2.0397e+00 | +2.0397e+00 | 0.00e+00 |
| 0.100 | -0.238 | -1.0487 | -1.0487 | +2.0397e+00 | +2.0397e+00 | 0.00e+00 |
| 0.100 | -0.119 | -1.0240 | -1.0240 | +2.0396e+00 | +2.0396e+00 | 0.00e+00 |

## T4: Triangle loop closure

- `C_to_BW`: (s,t,u) → (u,s,t)  perm "ust"
- `BW_to_ann`: (s,t,u) → (s,t,u)  identity (T-reversal)
- `ann_to_C`: (s,t,u) → (t,u,s)  perm "tus"  [inverse of C→ann]

Round-trip identity check at the Mandelstam label level and at the |M̄|²_KN amplitude level:

| original (s, t, u) | after C→BW | after BW→ann | after loop | label diff | amp diff |
|---|---|---|---|---:|---:|
| [2.0, -0.3, -0.7] | [-0.7, 2.0, -0.3] | [-0.7, 2.0, -0.3] | [2.0, -0.3, -0.7] | 0.00e+00 | 0.00e+00 |
| [3.0, -0.1, -1.9] | [-1.9, 3.0, -0.1] | [-1.9, 3.0, -0.1] | [3.0, -0.1, -1.9] | 0.00e+00 | 0.00e+00 |
| [5.0, -0.5, -3.5] | [-3.5, 5.0, -0.5] | [-3.5, 5.0, -0.5] | [5.0, -0.5, -3.5] | 0.00e+00 | 0.00e+00 |
| [10.0, -1.0, -8.0] | [-8.0, 10.0, -1.0] | [-8.0, 10.0, -1.0] | [10.0, -1.0, -8.0] | 0.00e+00 | 0.00e+00 |

## T5: Total annihilation cross section

| β | textbook ∫ | BAM ∫ | ratio | Dirac σ_ann |
|---:|---:|---:|---:|---:|
| 0.100 | +1.0098e+01 | +1.0098e+01 | +1.00000000 | +3.1723e+01 |
| 0.200 | +5.1810e+00 | +5.1810e+00 | +1.00000000 | +1.6276e+01 |
| 0.300 | +3.5697e+00 | +3.5697e+00 | +1.00000000 | +1.1215e+01 |
| 0.500 | +2.2158e+00 | +2.2158e+00 | +1.00000000 | +6.9610e+00 |
| 0.700 | +1.3912e+00 | +1.3912e+00 | +1.00000000 | +4.3706e+00 |
| 0.900 | +5.5824e-01 | +5.5824e-01 | +1.00000000 | +1.7537e+00 |
| 0.950 | +3.1993e-01 | +3.1993e-01 | +1.00000000 | +1.0049e+00 |
| 0.980 | +1.5501e-01 | +1.5501e-01 | +1.00000000 | +4.8637e-01 |

## T6: Threshold and ultra-relativistic limits

### Threshold (β → 0): σ ∼ π/β s-wave divergence

| β | textbook ∫ | BAM ∫ | textbook·β | BAM·β | ratio |
|---:|---:|---:|---:|---:|---:|
| 0.010 | +1.0001e+02 | +1.0001e+02 | +1.0001e+00 | +1.0001e+00 | +1.00000000 |
| 0.020 | +5.0020e+01 | +5.0020e+01 | +1.0004e+00 | +1.0004e+00 | +1.00000000 |
| 0.050 | +2.0050e+01 | +2.0050e+01 | +1.0025e+00 | +1.0025e+00 | +1.00000000 |
| 0.100 | +1.0098e+01 | +1.0098e+01 | +1.0098e+00 | +1.0098e+00 | +1.00000000 |

### Ultra-relativistic (β → 1): log envelope

| β | textbook ∫ | BAM ∫ | log env | ratio |
|---:|---:|---:|---:|---:|
| 0.9500 | +3.1993e-01 | +3.1993e-01 | +3.6636e+00 | +1.00000000 |
| 0.9900 | +8.9488e-02 | +8.9488e-02 | +5.2933e+00 | +1.00000000 |
| 0.9990 | +1.5600e-02 | +1.5600e-02 | +7.6004e+00 | +1.00000000 |

## Verdict

**TRIANGLE_CLOSES.** TRIANGLE CLOSES. The same closed-form Compton vertex factor F²(x, c) = 4·x³·(x²+1−x·sin²θ)/[(1+c²)(1+x)²], analytically continued via Mandelstam crossing, reproduces all three corners of the QED two-photon-two-fermion triangle: Compton (PR #35), Breit-Wheeler (PR #36), and annihilation (this PR). The Compton → BW and Compton → ann crossings share the same Mandelstam permutation (s,t,u) → (u,s,t); BW ↔ ann is T-reversal at the kernel level. The triangle loop π_C→BW ∘ π_BW→ann ∘ π_ann→C is identity, both at the Mandelstam label level and at the |M̄|²_KN amplitude level. The BAM tree kernel is process-general across the full crossing triangle.

## What this leaves open

- **Bhabha (e⁺e⁻ → e⁺e⁻) and Møller (e⁻e⁻ → e⁻e⁻).** Two-channel processes with both s- and t-channel diagrams interfering. The BAM kernel would need to combine two crossed copies coherently.
- **One-loop corrections.** Still tree-level only.
- **BAM first-principles derivation of F².** Still open from PR #35.
