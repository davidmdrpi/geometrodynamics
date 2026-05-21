# 5D → 4D radial reduction bridge probe

**Run:** 2026-05-21T04:13:37+00:00

Builds the Kaluza–Klein-like reduction of the 5D Tangherlini bulk and determines what it produces — connecting the mass (radial) and amplitude (F²) sub-threads, and honestly identifying the residual (B5).

## The reduction

```
Ψ(x^μ, r, Ω) = Σ ψ_{l,n}(x^μ)·u_{l,n}(r)·Y_l(Ω); integrate over (r, Ω) → three channels
```

Three channels:
  - **radial**: KK masses ω(l,n)
  - **s3**: c₁ = 1 + propagator 1/q²
  - **throat**: F²(x, c) form factor

**Central finding:** F²(x, c) is the throat form factor, NOT a radial overlap (overlaps are kinematics-independent constants)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_kk_reduction_basis` | Gram off-diag 1.7e-15; spectrum match 1.5e-04 | **PASS** |
| T2 | `T2_radial_channel_masses` | discrete KK mass spectrum: True | **PASS** |
| T3 | `T3_s3_channel_gauge_propagator` | c₁ = 1.0000; Coulomb residual 3.8e-06 | **PASS** |
| T4 | `T4_throat_channel_F2` | F² = K²·Q (max diff 2.8e-14) | **PASS** |
| T5 | `T5_F2_is_not_a_radial_overlap` | F² spread 10.50 (varies); overlaps constant — F² ≠ overlap | **PASS** |
| T6 | `T6_shared_substrate_consistency` | R_MID, 2π, T²=−I shared across 3 channels | **PASS** |
| T7 | `T7_b5_assessment` | 3-channel factorization; residual = master integral | **PASS** |

## T1: KK reduction basis

| l | ω (symmetric FD) | ω (Chebyshev) | Gram off-diag | spectrum diff |
|---:|---|---|---:|---:|
| 1 | 1.0547, 1.9744, 2.8940, 3.8245 | 1.0547, 1.9744, 2.8941, 3.8247 | 1.1e-15 | 1.5e-04 |
| 3 | 1.2191, 2.1412, 3.0219, 3.9223 | 1.2191, 2.1412, 3.0220, 3.9225 | 8.9e-16 | 1.5e-04 |
| 5 | 1.3960, 2.3694, 3.2277, 4.0859 | 1.3960, 2.3694, 3.2277, 4.0861 | 1.7e-15 | 1.5e-04 |

## T2: Radial channel → mass spectrum

| l | KK masses ω(l,n) | spacings | discrete |
|---:|---|---|:---:|
| 1 | 1.0547, 1.9744, 2.8941, 3.8247 | 0.9197, 0.9197, 0.9306 | True |
| 3 | 1.2191, 2.1412, 3.0220, 3.9225 | 0.9221, 0.8808, 0.9005 | True |
| 5 | 1.3960, 2.3694, 3.2277, 4.0861 | 0.9734, 0.8583, 0.8584 | True |

## T3: S³ channel → gauge coupling + propagator

Hopf charge c₁ = **1.000000** (error 1.29e-08); S³ Green function flat limit → Coulomb (residual 3.8e-06) → propagator 1/q² (PRs #45–#46).

## T4: Throat channel → F²

| x | cosθ | K²·Q | F² closed | diff |
|---:|---:|---:|---:|---:|
| 0.10 | -0.70 | 0.0021 | 0.0021 | 0.0e+00 |
| 0.10 | +0.00 | 0.0030 | 0.0030 | 4.3e-19 |
| 0.10 | +0.70 | 0.0021 | 0.0021 | 0.0e+00 |
| 0.50 | -0.70 | 0.1484 | 0.1484 | 2.8e-17 |
| 0.50 | +0.00 | 0.1667 | 0.1667 | 0.0e+00 |
| 0.50 | +0.70 | 0.1484 | 0.1484 | 2.8e-17 |

## T5: F² is NOT a radial overlap (central finding)

Radial overlaps (kinematics-independent constants):
  - `<u0|u0>` = +1.0000
  - `<u0|u1>` = -0.0000
  - `<u1|u1>` = +1.0000
  - `<u0|u2>` = -0.0000

F² across kinematics (varies strongly):
  - `F2(x=0.5,c=0)` = 0.1667
  - `F2(x=1.0,c=0)` = 1.0000
  - `F2(x=2.0,c=0)` = 10.6667
  - `F2(x=1.0,c=0.9)` = 1.0000

F² kinematic spread = **10.50** (overlaps have zero kinematic spread). F²(x,c) varies with the external kinematics, so it CANNOT be a radial overlap — it is the throat-channel form factor. The naive 'F² from radial integration' is falsified.

## T6: Shared-substrate consistency

| datum | radial channel | S³ channel | throat channel |
|---|---|---|---|
| `R_MID` | inner boundary of [R_MID, R_OUTER] | S³ Green function radius (unit) | throat-pinch location r = R_MID |
| `closure_quantum_2pi` | Bohr–Sommerfeld mode quantization | Hopf holonomy / winding | throat closure (K from 2π split) |
| `spin_structure_T2_minus_I` | hard-wall Dirichlet BC (B3) | RP³ = S³/Z₂ non-trivial spin structure | throat transport T = iσ_y (F² helicity) |

The three channels share R_MID, the closure quantum 2π, and the spin structure T² = −I — the bridge connecting the mass and amplitude sub-threads.

## T7: B5 assessment / scaffold status

Three channels of one reduction:
  - **radial**: r ∈ [R_MID, R_OUTER] → KK masses ω(l,n) [mass thread]
  - **s3_angular**: Ω ∈ S³ → c₁ = 1 + propagator 1/q² [gauge/prop]
  - **throat**: r → R_MID pinch → F²(x, c) form factor [vertex]

Scaffold status:
  - **B1**: closed (PR #49 winding θ-term)
  - **B2**: closed (PR #49 RP³ + spin structure)
  - **B3**: closed (hard-wall from spin structure)
  - **B4**: open (dimensional bridge, m_e anchor)
  - **B5**: substantially reduced (3-channel factorization; residual = master integral)

## Verdict

**BRIDGE_FACTORIZED.** BRIDGE FACTORIZED. The 5D → 4D radial reduction of the Tangherlini bulk action factorizes into three consistent channels, all reductions of one action on the same internal geometry:

  radial  (r ∈ [R_MID, R_OUTER]) → KK masses ω(l,n)    [mass sub-thread]
  S³      (Ω ∈ S³)               → c₁ = 1 + 1/q²        [gauge + propagator]
  throat  (r → R_MID pinch)      → F²(x, c) form factor [vertex]

The radial modes form an orthonormal KK basis (Gram ≈ I) with masses matching the Chebyshev solver; the S³ channel gives the Hopf charge c₁ = 1 and the Coulomb/1/q² propagator; the throat channel gives F² = K²·Q. All three share the substrate: R_MID (throat radius), the closure quantum 2π, and the spin structure T² = −I (the radial Dirichlet BC and the throat transport are the same datum). The mass and amplitude sub-threads are thereby structurally connected.

CENTRAL FINDING (T5): F²(x, c) is NOT a radial overlap integral. Radial overlaps ∫u_m u_n dr are kinematics-independent constants (orthonormality δ_mn); F²(x, c) varies strongly with the scattering kinematics (x, c). So F² is the throat-channel form factor, not a KK overlap — the naive "F² from radial integration" is falsified. The reduction connects the sub-threads structurally (shared substrate, three channels of one action) but does NOT unify masses and F² into a single master integral.

B5 is substantially reduced: from "the reduction map is unconstructed" to "the reduction factorizes into three consistent channels, with one residual — a single master integral producing masses AND F² together." Scaffold status: B1, B2, B3 closed; B5 reduced to its residual; B4 (dimensional bridge / m_e anchor) remains.

## What this leaves open

- **The master integral (B5 residual).** A single covariant reduction integral producing the mass spectrum AND the F² vertex shape together. The three-channel factorization shows they are all reductions of one action, but writing them as one integral requires treating the throat-pinch (boundary) dynamics and the bulk radial modes on the same footing — open.
- **B4 — dimensional bridge.** Unaffected; the single m_e anchor (ℏ = m_e·R_MID·c) remains the last external input.
