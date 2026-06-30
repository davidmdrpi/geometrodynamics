# Two-throat Hartree–Fock sandbox: direct plus exchange terms (PR #187)

**Run:** 2026-06-30T01:11:47+00:00

Convolves the #186 hardened overlap kernels with an interaction to build the overlap-normalized two-throat Hartree–Fock energy `E±(R) = (J(R) ± K_ex(R))/(1 ± S²)` — direct plus exchange — from the actual #180 throat-soliton orbitals. *(QFT on the classical throat, not quantum gravity.)*

- **Structure**: E±(R) = (J(R) ± K_ex(R))/(1 ± S²) : overlap-normalized; S=⟨φ_a|φ_b⟩
- **Numerators**: J, K_ex positive/decaying, J ≥ K_ex (unnormalized HF numerators)
- **Ordering**: fermion branch E₋ below boson E₊ at finite R — FOR A REPULSIVE V (attractive reverses)
- **Pauli**: at coincidence the antisymmetric state is the zero vector (forbidden); contact V → E₋=0 at finite R

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | assemble the overlap-normalized two-throat HF energy | **PASS** |
| T2 | `T2_hartree_fock_structure` | the HF structure: E± = (J ± K_ex)/(1 ± S²) | **PASS** |
| T3 | `T3_overlap_and_numerators` | the overlap S(R) and the numerators J(R), K_ex(R) | **PASS** |
| T4 | `T4_overlap_normalized_energies` | the normalized energies: E₋ below E₊ for a repulsive V | **PASS** |
| T5 | `T5_pauli_zero_vector_at_coincidence` | Pauli: the antisymmetric state is the zero vector (forbidden) | **PASS** |
| T6 | `T6_controls_and_convergence` | controls + convergence (far-vanishing; grid-convergent) | **PASS** |
| T7 | `T7_honest_scope` | honest scope (a sandbox) | **PASS** |
| T8 | `T8_assessment` | TWO_THROAT_HF_OVERLAP_NORMALIZED | **PASS** |

## Overlap, numerators, and the normalized energies

| R | S(R) | J (direct num.) | K_ex (exch. num.) | E₊ = (J+K_ex)/(1+S²) | E₋ = (J−K_ex)/(1−S²) |
|---:|---:|---:|---:|---:|---:|
| 1.0 | 0.7911 | 0.03101 | 0.02375 | 0.03368 | 0.0194 |
| 2.0 | 0.4093 | 0.01727 | 0.00579 | 0.01975 | 0.01379 |
| 3.0 | 0.1529 | 0.0086 | 0.00071 | 0.0091 | 0.00808 |
| 4.0 | 0.0454 | 0.0044 | 5e-05 | 0.00445 | 0.00436 |

For the **repulsive** screened `V`, the fermion branch `E₋` (the Pin⁻-selected antisymmetric sector) sits below the boson `E₊` at every finite separation; the gap closes as the overlap `S → 0`. At coincidence the antisymmetric state is the **zero vector** (Pauli-forbidden), not a zero-energy state.

## Verdict

**TWO_THROAT_HARTREE_FOCK_OVERLAP_NORMALIZED_DIRECT_PLUS_EXCHANGE_FERMION_LOWER_FOR_REPULSIVE_V_ANTISYM_FORBIDDEN_AT_COINCIDENCE.** ASSEMBLED — THE OVERLAP-NORMALIZED TWO-THROAT HARTREE–FOCK ENERGY. E±(R) = (J(R) ± K_ex(R))/(1 ± S²) from the GR soliton orbitals.

OVERLAP + NUMERATORS. The non-orthogonal throats have overlap S(R) = ⟨φ_a|φ_b⟩ (1 at contact → 0 far apart); the direct J(R) and exchange K_ex(R) numerators are positive and decay, the direct dominating (K_ex/J = 0.77 at R = 1).

NORMALIZED ENERGIES, FERMION LOWER (REPULSIVE V). For the repulsive screened interaction the antisymmetric branch E₋ = (J − K_ex)/(1 − S²) lies below the symmetric E₊ = (J + K_ex)/(1 + S²) at every finite separation (the gap closing as the overlap dies) — the GR-selected Pin⁻ state is the lower one; an attractive V would reverse this.

PAULI — THE ZERO VECTOR. At coincidence the antisymmetric state is the ZERO VECTOR (both J − K_ex → 0 and 1 − S² → 0): two identical throats in the same orbital are Pauli-FORBIDDEN, not a zero-energy state. (For a contact V the antisymmetric energy is zero at all finite R.)

CONTROLS. Both numerators and the overlap vanish at far separation (orthogonal, distinguishable), and the energies are grid-convergent (J(2) to 0.0%, K_ex(2) to 0.0%). The exchange interaction and the Pauli physics of throat matter, from the GR-derived kernel.
