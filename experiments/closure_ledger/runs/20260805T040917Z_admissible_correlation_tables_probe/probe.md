# Probe J — which correlation tables does BAM's encoding admit? (PR #237)

_Run 2026-08-05T04:09:17.612603+00:00 · 9/9 PASS_

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

## But convex mixing reaches the whole body

| generators | coverage | max residual |
|---:|---:|---:|
| 2000 | 99.2% | 1.2e-02 |
| 8000 | 100.0% | 0.0e+00 |
| 32000 | 100.0% | 0.0e+00 |

On 300 fresh random TLM tables: **98.7%**, support size at most **5** — the Carathéodory bound. PR box excluded (residual 1.18).

**The witness above decomposes into 5 `c=1` cosine tables** — residual 0.0e+00, reconstruction error 4.4e-16, weights [0.212356, 0.34472, 0.040881, 0.010237, 0.391806]. *So it was never operationally unreachable; #237's original ~45% headline was a single-shot statement.*

And the implemented readout has **identically zero marginals** (0.0e+00), so no biased-marginal behavior is representable at all — and **this restriction survives mixing**, so it is the one that stands.

## Verdict

**THE_REPRESENTATION_ADMITS_THE_GRAM_TLM_BODY_AND_SO_DOES_BAM_ONCE_MIXING_IS_ALLOWED_BECAUSE_THE_CONVEX_HULL_OF_THE_C_EQUALS_ONE_COSINE_TABLES_IS_THE_WHOLE_BODY_SO_THE_ONLY_SURVIVING_RESTRICTION_IS_IDENTICALLY_ZERO_MARGINALS**
