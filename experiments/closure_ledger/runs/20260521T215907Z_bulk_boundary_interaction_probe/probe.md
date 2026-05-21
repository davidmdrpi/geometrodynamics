# Bulk-boundary interaction probe

**Run:** 2026-05-21T21:59:07+00:00

Targets the B5′ residual (PR #50): unify the bulk radial modes (masses) and the boundary throat-pinch (F²) on the same footing. The same throat cavity produces both, via one bulk Green function and its boundary impedance.

## The bulk-boundary structure

```
one throat cavity → bulk Green function G(r,r′;ω) [poles = masses] + throat impedance Z(ω)=π/ω [series → K(x)=2x/(1+x)]
```

- **Unifies**: radial (masses) + throat (K)
- **Residual**: S³ (Q) channel; full F² = K²·Q

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_bulk_green_function_poles_are_masses` | G poles at masses: True | **PASS** |
| T2 | `T2_boundary_normal_derivatives` | throat couplings u'(throat) ≠ 0: True | **PASS** |
| T3 | `T3_throat_to_throat_response` | Π(ω) poles at masses: True | **PASS** |
| T4 | `T4_throat_impedance_to_K` | Z(ω)=π/ω series → K (max diff 2.2e-16) | **PASS** |
| T5 | `T5_one_cavity_two_outputs` | one cavity → masses (bulk) + K (boundary) | **PASS** |
| T6 | `T6_shared_substrate` | R_MID, hard-wall BC, closure quantum shared | **PASS** |
| T7 | `T7_b5prime_assessment` | radial+throat unified; S³ (Q) residual | **PASS** |

## T1: Bulk Green function poles = masses

| l | masses ω(l,n) | pole blow-up detected |
|---:|---|---|
| 1 | 1.0547, 1.9744, 2.8940, 3.8245 | [True, True, True] |
| 3 | 1.2191, 2.1412, 3.0219, 3.9223 | [True, True, True] |
| 5 | 1.3960, 2.3694, 3.2277, 4.0859 | [True, True, True] |

## T2: Boundary normal derivatives (throat couplings)

| l | u(throat) | u'(throat) couplings |
|---:|---:|---|
| 1 | 0.0 | -0.8185, -1.5049, +2.2217, +2.9486 |
| 3 | 0.0 | -0.9729, -1.6013, +2.2726, +2.9845 |
| 5 | 0.0 | -1.1629, -1.7793, +2.3690, +3.0423 |

## T3: Throat-to-throat response Π(ω)

| l | masses | throat couplings² | Π poles at masses |
|---:|---|---|---|
| 1 | 1.0547, 1.9744, 2.8940, 3.8245 | 0.6699, 2.2646, 4.9358, 8.6942 | [True, True, True] |
| 3 | 1.2191, 2.1412, 3.0219, 3.9223 | 0.9465, 2.5640, 5.1649, 8.9073 | [True, True, True] |

## T4: Throat impedance → K factor

| x | Z(1)=π | Z(x)=π/x | K series | 2x/(1+x) | diff |
|---:|---:|---:|---:|---:|---:|
| 0.10 | 3.1416 | 31.4159 | 0.1818 | 0.1818 | 2.8e-17 |
| 0.50 | 3.1416 | 6.2832 | 0.6667 | 0.6667 | 0.0e+00 |
| 1.00 | 3.1416 | 3.1416 | 1.0000 | 1.0000 | 0.0e+00 |
| 2.00 | 3.1416 | 1.5708 | 1.3333 | 1.3333 | 0.0e+00 |
| 5.00 | 3.1416 | 0.6283 | 1.6667 | 1.6667 | 2.2e-16 |
| 10.00 | 3.1416 | 0.3142 | 1.8182 | 1.8182 | 0.0e+00 |

## T5: One cavity, two outputs

- **Bulk output** (masses, l=1): 1.0547, 1.9744, 2.8940, 3.8245
- **Boundary output** (K at x=0.5): 0.6667
- Same cavity: R_MID = 1.0, R_OUTER = 1.26

## T6: Shared substrate

- R_MID = 1.0 (cavity inner wall = throat)
- Hard-wall BC: Dirichlet at throat (B3, from T²=−I)
- Closure half π = 3.141593; dwell time τ(1) = 3.141593

## T7: B5′ assessment

- **Unified**: radial (masses) + throat (K) via one bulk-boundary cavity
- **Residual**: combine with S³ (Q) channel for full F² = K²·Q

B5 progression:
  - **before_PR50**: reduction unconstructed
  - **PR50**: three-channel factorization; F² not a radial overlap
  - **this_probe**: radial+throat unified by bulk-boundary cavity; S³ residual

## Verdict

**BULK_BOUNDARY_FORMULATED.** BULK-BOUNDARY INTERACTION FORMULATED. The same throat cavity (r ∈ [R_MID, R_OUTER], Dirichlet hard walls) produces both outputs the B5′ residual asked to unify:

  BULK aspect — the mass spectrum as poles of the bulk Green function G(r,r′;ω) = Σ u_n(r)u_n(r′)/(ω²−ω_n²); the magnitude blows up at each radial mass ω(l,n).
  BOUNDARY aspect — the K factor as the series of throat dwell-time impedances Z(ω) = τ(ω) = π/ω: the Compton in/out photons see Z(ω), Z(ω′) in series → harmonic mean → K(x) = 2x/(1+x), to machine precision.

The link is the Dirichlet hard wall (B3): u_n(R_MID) = 0 but the normal derivative u_n′(R_MID) ≠ 0 is the throat coupling of each bulk mode, giving the throat-to-throat response Π(ω) = Σ[u_n′(R_MID)]²/(ω²−ω_n²) with poles at the masses. Bulk spectrum and boundary impedance share the substrate: R_MID, the hard-wall BC (B3, from T²=−I), and the closure quantum (B1, via the dwell time τ = π/ω).

This unifies the radial (mass) and throat (K) channels in ONE bulk-boundary structure — the master integral for those two channels. Residual: the S³ (Q, Hopf-fibre helicity) channel is not part of the throat cavity; the full vertex F² = K²·Q still combines the bulk-boundary (K + masses) with the S³ channel (Q). B5′ is narrowed from "three separate channels" to "radial+throat unified by the bulk-boundary cavity; S³ (Q) combined separately."

## What this leaves open

- **S³ (Q) combination.** The Hopf-helicity Q factor is the angular channel; combining it with the bulk-boundary (K + masses) to write the complete F² and the mass spectrum in one integral is the remaining piece of B5′.
- **B4 — dimensional bridge.** Unaffected; the single m_e anchor remains.
