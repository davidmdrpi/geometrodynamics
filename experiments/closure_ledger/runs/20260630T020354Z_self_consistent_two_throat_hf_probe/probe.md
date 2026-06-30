# Self-consistent two-throat Hartree–Fock relaxation (PR #189)

**Run:** 2026-06-30T02:03:54+00:00

Relaxes PR #187's rigid two-throat orbitals self-consistently — a genuine HF SCF (self-interaction-free Fock operator) lets the two same-spin throats deform in each other's direct + exchange field, lowering the energy to a self-consistent variational fixed point. *(QFT on the classical throat, not quantum gravity.)*

- **Convergence**: imaginary-time HF SCF, monotone descent to a fixed point (robust across seeds)
- **Relaxation**: the self-consistent energy lies below the rigid 1D #187-style reference (~2.5%)
- **Deformation**: the two-throat density polarizes in the mean field (RMS shift, fidelity<1)
- **Exchange**: turning off the Fock −K (consistent control) raises the energy — the exchange lowers it

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | relax the orbitals self-consistently in the direct + exchange field | **PASS** |
| T2 | `T2_hartree_fock_mean_field` | the HF mean field: F_i = h + J_{≠i} − K_{≠i} (self-interaction-free) | **PASS** |
| T3 | `T3_scf_convergence_and_robustness` | convergence: monotone descent to a fixed point, robust across seeds | **PASS** |
| T4 | `T4_relaxation_lowers_energy` | relaxation: the energy is below the rigid 1D #187-style reference | **PASS** |
| T5 | `T5_orbitals_deform_in_mean_field` | the orbitals deform: the density relaxes (RMS shift, fidelity<1) | **PASS** |
| T6 | `T6_exchange_field_lowers_energy` | the exchange field: turning off −K (consistent control) raises the energy | **PASS** |
| T7 | `T7_honest_scope` | honest scope (a 1D sandbox SCF; a variational fixed point) | **PASS** |
| T8 | `T8_assessment` | SELF_CONSISTENT_TWO_THROAT_HF_RELAXED | **PASS** |

## The relaxation

| quantity | value |
|---|---|
| energy: rigid 1D ref → relaxed | -3.80798 → -3.9046 (2.537% lower) |
| SCF convergence (final ΔE) | -4.4e-08 (monotone) |
| seeded-restart spread | 7.6e-03 (robust fixed point) |
| density RMS: rigid → relaxed | 2.6458 → 3.0554 |
| density fidelity (rigid vs relaxed) | 0.9783 |
| exchange lowering (E_Hartree − E_HF, consistent control) | 0.56703 |

## Verdict

**SELF_CONSISTENT_TWO_THROAT_HF_RELAXATION_CONVERGES_LOWERS_ENERGY_ORBITALS_DEFORM_IN_DIRECT_PLUS_EXCHANGE_FIELD.** RELAXED — A SELF-CONSISTENT TWO-THROAT HARTREE–FOCK VARIATIONAL FIXED POINT. The orbitals deform in each other's direct + exchange field.

CONVERGENT + ROBUST. The imaginary-time HF relaxation (with a self-interaction-free Fock operator) lowers the energy monotonically from -3.80798 to -3.9046, settling to a fixed point (ΔE = -4.4e-08 → 0) that is robust across seeded restarts (spread 7.6e-03).

RELAXATION. The self-consistent energy (-3.9046) lies below the rigid 1D #187-style reference (-3.80798) by 2.54% — the variational gain from optimizing the orbital shapes.

DEFORMED. The two-throat density polarizes in the mean field — RMS width 2.6458 → 3.0554, fidelity to the rigid density 0.9783 < 1.

EXCHANGE WORKS. With the consistent self-interaction-free control, turning off the non-local Fock exchange raises the energy by 0.56703 (E_HF = -3.9046 vs -3.33756): the same-spin exchange hole of #185–#188 substantially lowers the energy in the self-consistent mean field. SCOPE: a 1D sandbox SCF, a variational fixed point (not certified the global ground state); the full 3D self-gravitating two-throat solve is the follow-up.
