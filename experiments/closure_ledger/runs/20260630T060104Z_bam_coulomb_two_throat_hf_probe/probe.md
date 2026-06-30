# The BAM Coulomb-photon kernel for the two-throat HF: replacing the Yukawa stand-in (PR #190)

**Run:** 2026-06-30T06:01:04+00:00

Replaces the screened-photon (Yukawa) stand-in of #187/#189 with the genuine BAM Coulomb-photon kernel — the unscreened `1/(4πd) ⟷ 1/q²` (the flat limit of the S³ Green function, #42–#44), regulated by the Hockney open-boundary solver — and recomputes the two-throat direct + exchange energies. *(QFT on the classical throat, not quantum gravity.)*

- **Kernel**: unscreened BAM Coulomb photon V(d)=1/(4πd) ⟷ 1/q² (flat limit of the S³ Green fn)
- **Regulator**: Hockney open-boundary Coulomb (validated ~0.2% on the Gaussian self-energy)
- **Long-ranged**: the direct J(R) → 1/(4πR) (the Coulomb tail); the exchange stays short-ranged
- **Physics robust**: fermion-lower for the repulsive Coulomb; zero-vector Pauli state at coincidence

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | replace the Yukawa stand-in with the BAM Coulomb-photon kernel | **PASS** |
| T2 | `T2_bam_coulomb_photon_kernel` | the kernel: 1/(4πd) ⟷ 1/q²; the S³ Green-function flat limit | **PASS** |
| T3 | `T3_isolated_system_regulator` | the regulator: the isolated-system Coulomb, validated (~0.2%) | **PASS** |
| T4 | `T4_long_ranged_direct_short_exchange` | long-ranged direct: J(R) → 1/(4πR); exchange stays short-ranged | **PASS** |
| T5 | `T5_187_physics_robust_under_kernel` | the #187 physics robust: fermion-lower for the repulsive Coulomb | **PASS** |
| T6 | `T6_pauli_zero_vector_at_coincidence` | Pauli at coincidence: the antisymmetric state is the zero vector | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | BAM_COULOMB_PHOTON_KERNEL_REPLACES_YUKAWA | **PASS** |

## The recomputed energies (unscreened BAM Coulomb)

| R | direct J | exchange K_ex |
|---:|---:|---:|
| 0.0 | 0.0627 | 0.0627 |
| 1.0 | 0.05356 | 0.03833 |
| 2.0 | 0.03765 | 0.00963 |
| 3.0 | 0.02646 | 0.00123 |
| 4.0 | 0.01995 | 0.0001 |
| 6.0 | 0.0133 | 0.0 |

The direct `J → 1/(4πR)` (long-ranged Coulomb tail; `J(6)` ratio 1.0031), the exchange stays short-ranged; and the fermion branch `E₋` sits below the boson `E₊` at every finite separation for the repulsive photon (the #187 result survives).

## Verdict

**BAM_COULOMB_PHOTON_KERNEL_REPLACES_YUKAWA_STANDIN_LONG_RANGED_DIRECT_TWO_THROAT_HF_PHYSICS_ROBUST.** REPLACED — THE BAM COULOMB-PHOTON KERNEL, AND THE #187 PHYSICS SURVIVES. The screened Yukawa stand-in is gone.

THE KERNEL. The unscreened Coulomb 1/(4πd) ⟷ 1/q² (the photon propagator, #42–#44) — the flat limit of the S³ Green function (G·4πs = 0.9571 → 1) — regulated by the Hockney open-boundary solver (validated against the Gaussian self-energy to ratio 0.9992, ~0.1%).

LONG-RANGED DIRECT. The direct energy is now correctly long-ranged — J(6) = 0.0133 ≈ 1/(4π·6) (ratio 1.0031) — while the exchange stays short-ranged; the screened stand-in lacked this Coulomb tail.

PHYSICS ROBUST. With the overlap-normalized E± = (J ± K_ex)/(1 ± S²), for the repulsive photon the antisymmetric (Pin⁻) branch E₋ sits below the symmetric E₊ at every finite separation (fermion-lower survives), and at coincidence the antisymmetric state is the zero vector (Pauli-forbidden: J = K_ex, S → 1). The statistics are geometric, not an artifact of the stand-in. SCOPE: the S³ curvature corrections and the self-consistent SCF with the Coulomb kernel are follow-ups.
