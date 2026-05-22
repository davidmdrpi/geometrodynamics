# Master integral probe — S³ Hopf Q-channel integration

**Run:** 2026-05-22T03:32:07+00:00

Closes the B5′ residual (PR #51): integrate the S³ Hopf Q-channel into the bulk-boundary master functional, so a single functional yields the mass spectrum AND the full F²=K²·Q vertex.

## The master functional

```
ℳ(ω;x,c) = G_C(r,r′;ω) ⊗ 𝒢_{S³}(Ω,Ω′) on warped product C × S³ → poles=masses, throat boundary=K, S³ Hopf=Q; vertex residue = F²=K²·Q
```

- **Unifies**: radial (masses) + throat (K) + S³ (Q)
- **Remaining barrier**: B4 — dimensional bridge (single m_e anchor)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_separable_master_kernel` | separable; poles present: True | **PASS** |
| T2 | `T2_radial_times_s3_casimir_mass_spectrum` | product spectrum (mono l: True, mono n: True) | **PASS** |
| T3 | `T3_throat_boundary_to_K` | throat → K (max diff 2.2e-16) | **PASS** |
| T4 | `T4_s3_hopf_reduction_to_Q` | S³ Hopf → Q (max diff 2.8e-14) | **PASS** |
| T5 | `T5_master_vertex_residue_is_F2` | vertex = F²=K²·Q (max diff 2.1e-14) | **PASS** |
| T6 | `T6_one_functional_masses_and_F2` | one ℳ → masses + F² together | **PASS** |
| T7 | `T7_product_geometry_is_the_mechanism` | product geometry (sep. of vars: True) | **PASS** |
| T8 | `T8_shared_substrate` | R_MID, 2π, T²=−I shared across channels | **PASS** |
| T9 | `T9_b5prime_assessment` | radial+throat+S³ unified; B5′ closed | **PASS** |

## T2: Mass spectrum = radial ladder × S³ Casimir

| l (S³ Casimir) | ω(l,0) | ω(l,1) | ω(l,2) | ω(l,3) |
|---:|---:|---:|---:|---:|
| 1 | 1.0547 | 1.9744 | 2.8940 | 3.8245 |
| 3 | 1.2191 | 2.1412 | 3.0219 | 3.9223 |
| 5 | 1.3960 | 2.3694 | 3.2277 | 4.0859 |
| 7 | 1.5587 | 2.6054 | 3.4853 | 4.3126 |

## T3: Throat boundary → K

| x | Z(1)=π | Z(x)=π/x | K series | 2x/(1+x) | diff |
|---:|---:|---:|---:|---:|---:|
| 0.10 | 3.1416 | 31.4159 | 0.1818 | 0.1818 | 2.8e-17 |
| 0.50 | 3.1416 | 6.2832 | 0.6667 | 0.6667 | 0.0e+00 |
| 1.00 | 3.1416 | 3.1416 | 1.0000 | 1.0000 | 0.0e+00 |
| 2.00 | 3.1416 | 1.5708 | 1.3333 | 1.3333 | 0.0e+00 |
| 5.00 | 3.1416 | 0.6283 | 1.6667 | 1.6667 | 2.2e-16 |
| 10.00 | 3.1416 | 0.3142 | 1.8182 | 1.8182 | 0.0e+00 |

## T4: S³ Hopf reduction → Q

| x | cosθ | A_pres | A_flip | Q Hopf | Q target | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.050 | -0.900 | 0.0500 | 0.1579 | 0.027431 | 0.027431 | 0.0e+00 |
| 0.050 | -0.600 | 0.0500 | 0.1822 | 0.035680 | 0.035680 | 0.0e+00 |
| 0.050 | -0.300 | 0.0500 | 0.2035 | 0.043899 | 0.043899 | 0.0e+00 |
| 0.050 | -0.000 | 0.0500 | 0.2124 | 0.047625 | 0.047625 | 0.0e+00 |
| 0.050 | +0.300 | 0.0500 | 0.2035 | 0.043899 | 0.043899 | 0.0e+00 |
| 0.050 | +0.600 | 0.0500 | 0.1822 | 0.035680 | 0.035680 | 0.0e+00 |
| 0.050 | +0.900 | 0.0500 | 0.1579 | 0.027431 | 0.027431 | 0.0e+00 |
| 0.200 | -0.900 | 0.2000 | 0.2659 | 0.110718 | 0.110718 | 2.8e-17 |

Hopf helicity-sum identity max diff: 3.33e-16

## T5: Master vertex residue = F² = K²·Q

| x | cosθ | K | Q | K²·Q | F² closed | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0500 | -0.900 | 0.0952 | 0.0274 | 0.000249 | 0.000249 | 1.1e-19 |
| 0.0500 | -0.750 | 0.0952 | 0.0314 | 0.000285 | 0.000285 | 1.1e-19 |
| 0.0500 | -0.600 | 0.0952 | 0.0357 | 0.000324 | 0.000324 | 1.1e-19 |
| 0.0500 | -0.450 | 0.0952 | 0.0400 | 0.000363 | 0.000363 | 1.1e-19 |
| 0.0500 | -0.300 | 0.0952 | 0.0439 | 0.000398 | 0.000398 | 5.4e-20 |
| 0.0500 | -0.150 | 0.0952 | 0.0466 | 0.000423 | 0.000423 | 5.4e-20 |

**Max |vertex − F²| over the (x,c) grid: 2.13e-14**

## T6: One functional, masses AND F²

- **Mass output** (l=1): 1.0547, 1.9744, 2.8940, 3.8245
- **Mass output** (l=3): 1.2191, 2.1412, 3.0219, 3.9223
- **Vertex output**:
  - F2(x=0.5,c=0) = 0.166667
  - F2(x=1,c=0) = 1.000000
  - F2(x=2,c=0.5) = 9.955556
- masses ok: True; vertex matches F²: True

## T7: Product geometry is the mechanism

- Radial Gram off-diagonal: 7.10e-16 (orthonormal modes per l)
- Radial Gram diagonal error: 1.11e-16
- Separation of variables: True
- S³ Casimir shift (l=1→5): 0.3412, 0.3950, 0.3337 (monotone: True)

## T8: Shared substrate across all three channels

| substrate | radial (mass) | throat (K) | S³ (Q) |
|---|---|---|---|
| `R_MID` | cavity inner wall | throat radius (dwell) | S³ radius scale |
| `closure_quantum_2pi` | mode normalization | dwell time tau=pi/omega (half=3.141593) | Hopf holonomy pi cos chi (pole=3.141593) |
| `T2_eq_minus_I` | Dirichlet hard wall (B3) | throat node u(R_MID)=0 | helicity-flip epsilon (A_flip) |

- dwell half (π) = Hopf holonomy at pole (π): True

## T9: B5′ assessment

- **Unified**: radial (masses) + throat (K) + S³ (Q) in one functional
- **Remaining**: B4 — dimensional bridge (single m_e anchor)

B5 progression:
  - **before_PR50**: reduction unconstructed
  - **PR50**: three-channel factorization; F² not a radial overlap
  - **PR51**: radial+throat unified by bulk-boundary cavity
  - **this_probe**: S³ (Q) integrated → masses AND F² from one ℳ; B5′ closed

## Verdict

**MASTER_INTEGRAL_COMPLETE.** MASTER INTEGRAL COMPLETE. The S³ Hopf Q-channel is integrated into the bulk-boundary master functional. A SINGLE separable functional ℳ = G_C ⊗ 𝒢_{S³} on the warped-product internal geometry C × S³ (C = radial cavity [R_MID, R_OUTER]) yields all three channels:

  RADIAL — the mass spectrum as the ω-poles of ℳ, the product spectrum ω(l,n) (radial ladder n × S³ Casimir l, the latter entering as the centrifugal barrier of the warped product); monotone in both l and n.
  THROAT — the K factor as the throat-boundary dwell-time impedance Z(ω)=π/ω of the in/out photons in series → harmonic mean → K(x)=2x/(1+x).
  S³ HOPF — the Q factor as the Hopf-fibre helicity spinor (A_pres=x, A_flip=√x(1−x)/√(1+c²)) → Q=x²+x(1−x)²/(1+c²), with (1+c²)/2 = cos⁴(θ/2)+sin⁴(θ/2) the Hopf-fibre helicity sum.

The vertex residue of the SAME ℳ — throat boundary (K) dressed by the S³ Hopf reduction (Q) — reproduces the closed-form F²(x,c)=K(x)²·Q(x,c) to machine precision over the (x,c) grid, while its poles give the mass spectrum: masses AND the full F² vertex from ONE object — the master integral the B5′ residual asked for.

The F²=K²·Q factorization is NOT a failure to unify — it is the direct consequence of the product internal geometry: separation of variables Ψ=Σ u_{l,n}(r)𝒴_l(Ω) makes the internal Green function a sum of factor products. All three channels share the substrate: R_MID (cavity wall = throat = S³ scale), the closure quantum 2π (dwell time τ=π/ω for K; Hopf holonomy π cos χ for Q), and T²=−I (Dirichlet wall for the masses; throat node; helicity-flip ε for Q). B5′ is closed; B4 (the single m_e dimensional anchor) is the only surviving barrier.

## What this leaves open

- **B4 — dimensional bridge.** Unchanged: the single m_e anchor (ℏ = m_e·R_MID·c) sets the absolute MeV scale. The master integral is dimensionless; B4 is orthogonal.
- **First-principles internal action.** The master functional uses the established channel reductions (radial Sturm–Liouville, throat dwell-time, Hopf helicity); writing all three as one covariant 5D Lagrangian density with the throat boundary term is the natural follow-on.
