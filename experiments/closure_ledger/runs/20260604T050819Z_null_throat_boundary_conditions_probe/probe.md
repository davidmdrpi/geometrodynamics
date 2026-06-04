# Null throat boundary conditions for wave transport on the 5D horizon (PR #129)

**Run:** 2026-06-04T05:08:19+00:00

Derives the boundary condition the null throat (the 5D horizon) imposes on transported waves. The answer is the BAM-native antipodal condition — l-parity-graded (even-l Neumann, odd-l Dirichlet) — which makes the throat a unitary mirror, not an absorbing horizon. Builds on the PR #128 antipodal Kruskal structure and the PR #116 cavity operator.

- **Vanishing potential**: V_l ∝ f → 0 at throat; null modes ψ ~ e^{±iωr*}
- **Antipodal BC**: even-l Neumann (ψ'=0), odd-l Dirichlet (ψ=0) from Y_l parity (−1)^l
- **Unitarity**: real BC ⟹ zero throat flux (unitary mirror); ingoing BC ⟹ flux −ω (absorbing)
- **Spectrum**: real, positive, discrete; even-l (N) vs odd-l (D) distinct
- **Open**: full QNM spectrum (complex ω); throat↔antithroat nucleation rate (#58/#88)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | derive the null throat BC for wave transport on the 5D horizon | **PASS** |
| T2 | `T2_vanishing_potential_null_modes` | V_l ∝ f → 0 at the throat ⟹ near-horizon null modes e^{±iωr*} | **PASS** |
| T3 | `T3_three_candidate_bcs` | three candidate BCs: ingoing/absorbing, reflective, antipodal | **PASS** |
| T4 | `T4_antipodal_map_fixes_l_parity_bc` | Y_l(−x) = (−1)^l Y_l ⟹ even-l Neumann, odd-l Dirichlet | **PASS** |
| T5 | `T5_unitary_mirror_flux_conservation` | unitary mirror: real BC ⟹ zero flux; ingoing ⟹ flux −ω | **PASS** |
| T6 | `T6_real_discrete_l_graded_spectrum` | real discrete spectrum; even-l (N) vs odd-l (D) distinct | **PASS** |
| T7 | `T7_scope` | scope: kinematic BC established; QNM / nucleation open | **PASS** |
| T8 | `T8_assessment` | NULL_THROAT_BC_ANTIPODAL_L_PARITY_UNITARY_MIRROR | **PASS** |

## The antipodal map fixes the BC by l-parity (Y_l(−x) = (−1)^l Y_l)

| l | antipodal parity Y_l(−x)/Y_l(x) | (−1)^l | throat BC |
|---:|---:|---:|---|
| 0 | 1.0 | 1 | Neumann (ψ'=0, antinode) |
| 1 | -1.0 | -1 | Dirichlet (ψ=0, node) |
| 2 | 1.0 | 1 | Neumann (ψ'=0, antinode) |
| 3 | -1.0 | -1 | Dirichlet (ψ=0, node) |
| 4 | 1.0 | 1 | Neumann (ψ'=0, antinode) |

## The exterior cavity spectrum is real, discrete, and l-graded

| l | antipodal BC | lowest ω² |
|---:|:---:|---|
| 0 | N | [1.368, 10.629, 28.99] |
| 1 | D | [5.269, 19.047, 41.995] |
| 2 | N | [2.025, 11.632, 30.004] |
| 3 | D | [6.728, 20.582, 43.527] |

Real, positive, discrete (a unitary bound cavity); even-l (Neumann-at-throat) and odd-l (Dirichlet-at-throat) families give distinct spectra — the wave-transport face of BAM's even-k/odd-k Z₂ structure (#67/#121).

## Verdict

**NULL_THROAT_BC_ANTIPODAL_L_PARITY_UNITARY_MIRROR.** THE NULL THROAT IMPOSES THE ANTIPODAL, l-PARITY-GRADED BOUNDARY CONDITION, MAKING IT A UNITARY MIRROR. PR #128 built the horizon-regular charts and identified the antipodal map as BAM's throat ↔ antithroat C-swap; this probe derives the boundary condition that antipodal null throat imposes on transported waves.

THE POTENTIAL VANISHES AT THE HORIZON. The wave equation □Φ = 0, separated as Φ = e^{−iωt} Y_l(Ω) ψ(r)/r^{3/2}, gives −d²ψ/dr*² + V_l ψ = ω²ψ with V_l = f[l(l+2)/r² + 3rs²/r⁴] (PR #116). Since V_l ∝ f → 0 at the throat, the near-horizon equation is −ψ = ω²ψ: the modes are the pure null phases ψ ~ e^{±iωr*} (the ingoing/outgoing null rays of PR #128).

THREE CANDIDATE BCs. (1) ingoing/absorbing ψ ~ e^{−iωr*} (flux lost into the horizon, non-unitary); (2) reflective wall (Dirichlet/Neumann, real discrete spectrum — the cavity of #116); (3) antipodal (BAM-native): the #128 identification Φ(U,V,Ω) = Φ(−U,−V,Ω_antipodal).

THE ANTIPODAL MAP FIXES THE BC BY l-PARITY. The S³ scalar harmonics carry antipodal parity Y_l(−x) = (−1)^l Y_l(x) (degree-l harmonic polynomials on ℝ⁴; verified). Single-valuedness of the untwisted scalar under the antipodal identification then forces the radial function to carry the compensating (−1)^l across the throat: EVEN l ⟹ Neumann ∂ψ(throat) = 0 (an antinode), ODD l ⟹ Dirichlet ψ(throat) = 0 (a node). A twisted/Möbius field flips even ↔ odd, connecting to the Z₂ orientation grading (#67/#121).

THE THROAT IS A UNITARY MIRROR. Both antipodal BCs are REAL, so the Klein–Gordon/Wronskian flux j ∝ Im(ψ* ∂ψ) through the throat VANISHES — a perfect unitary mirror, no net flux lost. This is the sharp contrast with the ingoing/absorbing horizon BC, whose flux j = −ω|A|² ≠ 0 carries probability into the hole. The antipodal throat conserves flux: what falls toward it on one sheet re-emerges on the antipodal sheet (global CPT/unitarity, #64).

THE SPECTRUM IS REAL, DISCRETE, AND l-GRADED. On the exterior [R_MID+ε, R_OUTER] with the shell wall (Dirichlet at R_OUTER) and the antipodal throat BC, −d²/dr*² + V_l has a real, positive, discrete spectrum (a unitary bound cavity), with even-l (Neumann) and odd-l (Dirichlet) families giving distinct spectra — the wave-transport face of BAM's even-k/odd-k structure (#67/#121).

SCOPE. This establishes the kinematic BC structure of classical wave transport across the null throat: the vanishing potential, the l-parity-graded antipodal BC, the unitary-mirror flux, the real discrete spectrum. It does NOT solve the full quasinormal spectrum or the dynamical throat ↔ antithroat nucleation amplitude (#58/#88); the absorbing-vs-antipodal choice is fixed by the BAM antipodal postulate (#128), shown here to be self-consistent and unitary.
