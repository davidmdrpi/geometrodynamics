# Self-consistent two-throat Hartree–Fock relaxation (PR #189)

**Run:** 2026-06-30T01:50:30+00:00

Relaxes PR #187's rigid two-throat orbitals self-consistently — a genuine HF SCF lets the two same-spin throats deform in each other's direct + exchange field, lowering the energy to its variational minimum. *(QFT on the classical throat, not quantum gravity.)*

- **Convergence**: imaginary-time HF SCF, monotone energy descent to a fixed point
- **Relaxation**: the self-consistent energy lies below the rigid #187 value (~2.5%)
- **Deformation**: the two-throat density polarizes in the mean field (RMS shift, fidelity<1)
- **Exchange**: turning off the Fock −K raises the energy — the exchange lowers it

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | relax the orbitals self-consistently in the direct + exchange field | **PASS** |
| T2 | `T2_hartree_fock_mean_field` | the HF mean field: F = h + V_H − K (the Fock operator) | **PASS** |
| T3 | `T3_scf_convergence` | convergence: the SCF lowers the energy monotonically to a fixed point | **PASS** |
| T4 | `T4_relaxation_lowers_energy` | relaxation: the self-consistent energy is below the rigid value | **PASS** |
| T5 | `T5_orbitals_deform_in_mean_field` | the orbitals deform: the density relaxes (RMS shift, fidelity<1) | **PASS** |
| T6 | `T6_exchange_field_lowers_energy` | the exchange field: turning off −K raises the energy | **PASS** |
| T7 | `T7_honest_scope` | honest scope (a 1D sandbox SCF) | **PASS** |
| T8 | `T8_assessment` | SELF_CONSISTENT_TWO_THROAT_HF_RELAXED | **PASS** |

## The relaxation

| quantity | value |
|---|---|
| energy: rigid → relaxed | -3.80798 → -3.9046 (2.537% lower) |
| SCF convergence (final ΔE) | -4.4e-08 (monotone) |
| density RMS: rigid → relaxed | 2.6458 → 3.0554 |
| density fidelity (rigid vs relaxed) | 0.9783 |
| exchange lowering (E_Hartree−E_HF) | 0.61613 |

## Verdict

**SELF_CONSISTENT_TWO_THROAT_HF_RELAXATION_CONVERGES_LOWERS_ENERGY_ORBITALS_DEFORM_IN_DIRECT_PLUS_EXCHANGE_FIELD.** RELAXED — A SELF-CONSISTENT TWO-THROAT HARTREE–FOCK GROUND STATE. The orbitals deform in each other's direct + exchange field.

CONVERGENT. The imaginary-time HF relaxation lowers the energy monotonically from -3.80798 to -3.9046, settling to a self-consistent fixed point (ΔE = -4.4e-08 → 0).

RELAXATION. The self-consistent energy (-3.9046) lies below the rigid #187-style value (-3.80798) by 2.54% — the variational gain from optimizing the orbital shapes.

DEFORMED. The two-throat density polarizes in the mean field — RMS width 2.6458 → 3.0554, fidelity to the rigid density 0.9783 < 1 — the throats are no longer rigid.

EXCHANGE WORKS. Turning off the non-local Fock exchange (Hartree only) raises the energy by 0.61613 (E_HF = -3.9046 vs -3.28847): the same-spin exchange hole of #185–#188 substantially lowers the energy in the self-consistent mean field. SCOPE: a 1D sandbox SCF; the full 3D self-gravitating two-throat solve is the follow-up.
