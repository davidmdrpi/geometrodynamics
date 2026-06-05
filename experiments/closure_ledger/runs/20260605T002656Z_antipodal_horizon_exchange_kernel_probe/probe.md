# Antipodal-horizon exchange kernel from the throat boundary data (PR #135)

**Run:** 2026-06-05T00:26:56+00:00

Builds the matter-sector exchange kernel — the two-point Green's function / resolvent of the matter cavity operator (#116) with the antipodal horizon boundary data (#129). It is the BAM matter propagator: a reciprocal, unitary kernel assembled as a sum over the stable modes (#130), with l-channels graded by the antipodal sign (−1)^l. (The gauge-sector photon kernel 1/q² is the separate PR #42–#44 probe.)

- **Kernel**: K_l(r,r';ω) = (H_l − ω²)^{-1}, antipodal BC (#129)
- **Spectral**: K_l = Σ_n ψ_n ψ_n/(ω_n² − ω²), poles = real #130 spectrum
- **Reciprocity**: K_l(r,r') = K_l(r',r) (self-adjoint ⟹ symmetric)
- **Unitarity**: antipodal real poles ⟹ unitary; absorbing complex poles ⟹ lossy (#130)
- **Parity grading**: l-channels carry (−1)^l (antipodal exchange sign, #129/#134, C-swap #63)
- **Open**: interacting/multi-loop kernel; absolute normalisation (#133); flavor residuals (#134)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | build the matter exchange kernel from the antipodal data (#129) | **PASS** |
| T2 | `T2_kernel_is_cavity_resolvent` | kernel = cavity resolvent K_l = (H_l − ω²)^{-1}, self-adjoint | **PASS** |
| T3 | `T3_spectral_representation_mode_sum` | spectral rep K_l = Σ ψψ/(ω_n²−ω²); poles = #130 spectrum | **PASS** |
| T4 | `T4_reciprocity_symmetric_kernel` | reciprocity K_l(r,r') = K_l(r',r) (~1e-14) | **PASS** |
| T5 | `T5_unitary_vs_lossy_kernel` | antipodal real poles (unitary) vs absorbing complex (lossy, #130) | **PASS** |
| T6 | `T6_angular_antipodal_parity_grading` | angular channels antipodal-parity-graded (−1)^l (#129/#134) | **PASS** |
| T7 | `T7_scope` | scope: free/one-loop kernel; interacting kernel / normalisation open | **PASS** |
| T8 | `T8_assessment` | ANTIPODAL_HORIZON_EXCHANGE_KERNEL_UNITARY_RECIPROCAL_PARITY_GRADED | **PASS** |

## Spectral representation: poles = the stable #130 spectrum

| l | poles ω² (kernel) | mode-sum vs resolvent |
|---:|---|---:|
| 0 | [1.363, 10.581, 28.856] | 4.4e-14 |
| 1 | [5.269, 19.046, 41.99] | 4.6e-16 |
| 2 | [2.017, 11.583, 29.868] | 5.2e-15 |
| 3 | [6.728, 20.581, 43.522] | 3e-16 |

The exchange kernel `K_l = Σ_n ψ_n(r)ψ_n(r')/(ω_n² − ω²)` is a sum over the stable modes; its poles are the real #130 spectrum. (Mode sum equals the matrix resolvent to ~1e-14.)

## Angular antipodal-parity grading of the kernel

| l | antipodal sign (−1)^l | under C-swap |
|---:|---:|---|
| 0 | 1 | symmetric |
| 1 | -1 | antisymmetric |
| 2 | 1 | symmetric |
| 3 | -1 | antisymmetric |

Each l-channel of `K(x,x') = Σ_l K_l(r,r';ω) C_l(Ω·Ω')` carries the antipodal sign `(−1)^l` under the throat ↔ antithroat exchange — the same `(−1)^l` that fixed the boundary condition (#129).

## Verdict

**ANTIPODAL_HORIZON_EXCHANGE_KERNEL_UNITARY_RECIPROCAL_PARITY_GRADED.** THE ANTIPODAL-HORIZON MATTER EXCHANGE KERNEL IS THE ANTIPODAL-BC CAVITY RESOLVENT — RECIPROCAL, UNITARY, AND ANTIPODAL-PARITY-GRADED. PR #129 fixed the boundary data the antipodal horizon imposes and PR #130 the stable spectrum; this probe builds the object that boundary data defines, the matter two-point exchange kernel. (The gauge-sector photon kernel 1/q² is the separate PR #42–#44 probe.)

THE KERNEL = THE CAVITY RESOLVENT. For each angular channel l the exchange kernel is the resolvent K_l(r,r';ω) = ⟨r|(H_l − ω²)^{-1}|r'⟩ of the matter cavity operator H_l = −d²/dr*² + V_l (#116) with the antipodal boundary data of #129 — Neumann for even l, Dirichlet for odd l, Dirichlet shell wall. The operator is exactly self-adjoint.

SPECTRAL REPRESENTATION = EXCHANGE OF THE STABLE MODES. Because the antipodal operator is self-adjoint, the kernel is the mode sum K_l = Σ_n ψ_n(r)ψ_n(r')/(ω_n² − ω²), with poles at the real normal-mode spectrum of #130. The exchange kernel is a propagator assembled as a sum over the stable exchanged modes, with no decaying contribution (mode sum = matrix resolvent to ~1e-14).

RECIPROCITY. The antipodal boundary data makes H_l self-adjoint, so the kernel is symmetric, K_l(r,r') = K_l(r',r) (reciprocal exchange; |K − Kᵀ|/|K| ~ 1e-14).

UNITARY VS LOSSY — THE BOUNDARY DATA DECIDES. The antipodal (real) boundary data gives a Hermitian H_l with REAL poles ⟹ an undamped, unitary exchange kernel; the absorbing horizon (ingoing BC) would give a non-Hermitian operator with COMPLEX poles ⟹ a lossy kernel (#130). So the antipodal horizon boundary data is exactly what makes the matter exchange kernel unitary — the propagator-level face of the unitary mirror (#129) and the global CPT/unitarity (#64).

ANGULAR ANTIPODAL-PARITY GRADING. The full kernel factorises K(x,x') = Σ_l K_l(r,r';ω) C_l(Ω·Ω') with C_l the S³ zonal harmonic; under the throat ↔ antithroat exchange Ω'→AΩ', C_l(−Ω·Ω') = (−1)^l C_l(Ω·Ω') (#129/#134), so each l-channel of the exchange kernel carries the antipodal sign (−1)^l — even-l symmetric, odd-l antisymmetric under the C-swap (#63) — the same (−1)^l that fixed the boundary condition.

SCOPE. This constructs the FREE / one-loop matter exchange kernel on the fixed antipodal background — the propagator of the S_BAM fluctuation measure (#115–#122) with the #129 boundary data. It does NOT include the interacting / multi-loop kernel (vertices, self-energy) or fix the absolute normalisation; the bulk-scale (#133) and flavor (#134) residuals stand.
