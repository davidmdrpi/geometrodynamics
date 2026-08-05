# Probe J — which correlation tables does BAM's encoding admit? (PR #237)

_Run 2026-08-05T03:14:56.752961+00:00 · 8/8 PASS_

## Level A — the representation in the abstract

Exactly the unit-vector Gram matrices `E_xy = ⟨u_x, v_y⟩`, equivalently the TLM body. Verified both directions: fit vs predicate 60/60, and Gram-built tables never violate TLM (min slack +2.6e-04).

| body | fraction of the correlation cube |
|---|---:|
| local polytope | 66.83% |
| quantum (TLM) | 92.64% |
| no-signaling | 100% |

## Level B — what BAM implements

Read off the lattice: the T-matrix x–y block is exactly `c·I₂` with `c = −sin 2β`, so the achievable tables are exactly **`E_xy = c·cos(θ_x − φ_y)`**.

| `β` | `c` | `−sin 2β` | deviation from `c·I₂` |
|---:|---:|---:|---:|
| 0.2618 | -0.50000 | -0.50000 | 2.1e-15 |
| 0.3927 | -0.70711 | -0.70711 | 3.1e-15 |
| 0.5236 | -0.86603 | -0.86603 | 3.0e-15 |
| 0.7854 | -1.00000 | -1.00000 | 2.8e-15 |

| coverage of the quantum body | |
|---|---:|
| maximal preparation `c = ±1` (d=2 Gram) | **0.0%** (measure zero) |
| with the bridge preparation `β` free | **55.2%** |
| general `d = 3` | 100.0% |

Of the `c = 1` surface, 33.7% lies exactly on the Tsirelson boundary — which is why #236 found saturation there.

**An explicit unreachable quantum table:**

```
[[0.7424, -0.9971], [-0.5465, 0.4867]]
```

TLM slack +0.7400 (quantum) · `d=3` residual 4.3e-09 · BAM-family residual 4.6e-01

And the implemented readout has **identically zero marginals** (0.0e+00), so no biased-marginal behavior is representable at all.

## Verdict

**THE_REPRESENTATION_ADMITS_EXACTLY_THE_UNIT_VECTOR_GRAM_TLM_BODY_BUT_BAM_IMPLEMENTS_ONLY_C_COS_THETA_MINUS_PHI_WHICH_IS_MEASURE_ZERO_AT_MAXIMAL_PREPARATION_AND_ABOUT_HALF_THE_BODY_WITH_BETA_FREE_AND_HAS_ZERO_MARGINALS**
