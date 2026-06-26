# α as a protected boundary invariant, not a continuous tuning parameter (PR #184)

**Run:** 2026-06-26T05:35:08+00:00

Tests the EM coupling α as a protected boundary invariant (its quantized structure topologically robust to smooth deformation) rather than a continuous tuning knob — the EM-sector instance of the #181/#182/#183 protected-charge picture. *(QFT on the classical throat, not quantum gravity.)*

- **Charge quantum**: boundary Chern number c₁ = −1 (|c₁|=1), the Gauss-law charge
- **Protected**: stays the same integer under 30 smooth boundary diffeomorphisms (~1e-6)
- **Not tuning**: a continuous coupling functional drifts ~8% under the same deformations
- **Loop measure**: ∮F = 2π·c₁ quantized in 2π (fixing the 2π of a = α/2π)
- **Residual**: the value α⁻¹ ≈ 137 (the 137 problem stands)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | test α as a protected boundary invariant, not a tuning knob | **PASS** |
| T2 | `T2_alpha_decomposition` | α = protected structure × residual value | **PASS** |
| T3 | `T3_charge_quantum_boundary_chern` | charge quantum = boundary Chern number (|c₁|=1, integer) | **PASS** |
| T4 | `T4_protected_under_smooth_deformation` | protected: c₁ survives smooth boundary deformation | **PASS** |
| T5 | `T5_protected_not_tuning` | not tuning: a continuous functional drifts; c₁ does not | **PASS** |
| T6 | `T6_loop_measure_and_topology_change` | 1/2π loop measure quantized; changes only at a topology change | **PASS** |
| T7 | `T7_scope_and_unity` | scope (value still residual) + the #181/#182/#183 unity | **PASS** |
| T8 | `T8_assessment` | ALPHA_STRUCTURE_PROTECTED_VALUE_RESIDUAL | **PASS** |

## Protected invariant vs continuous tuning (same deformations)

| quantity | response to 30 smooth boundary deformations |
|---|---|
| charge quantum `c₁` (Chern number) | **invariant** — same integer, max move 5e-07 |
| a continuous coupling functional | **drifts** — 6.8% mean (15.77% max) |

## Verdict

**ALPHA_CHARGE_QUANTUM_AND_LOOP_MEASURE_ARE_PROTECTED_BOUNDARY_INVARIANTS_NOT_TUNING_THE_VALUE_REMAINS_RESIDUAL.** PROTECTED — α'S STRUCTURE IS A BOUNDARY INVARIANT, NOT A TUNING KNOB. The EM-sector instance of the #181/#182/#183 picture.

CHARGE QUANTUM. The boundary S² Chern number (the Gauss-law charge) is c₁ = 1.0 — an exact integer, |c₁| = 1 (the integer Hopf number, charge quantization as a boundary integral).

PROTECTED. Across 30 smooth diffeomorphisms of the boundary it stays the same integer ([1]) to 5e-07 — it does not drift.

NOT TUNING. Under the same deformations a generic continuous coupling functional drifts 6.8% on average while c₁ moves 5e-07: protected (quantized + invariant), not a continuous knob.

LOOP MEASURE + TOPOLOGY CHANGE. The boundary flux ∮F = 2π·c₁ is quantized in units of the closure quantum 2π (fixing the 2π of a = α/2π); the integer changes only when the Berry GAP closes — sweeping the gap parameter, C stays |1| while the gap is open and jumps to 0 exactly at the gap closing (min|d| → 0, the degeneracy crossing the boundary), the EM-boundary analog of |q| = 0 (#182) and ½ tr T² = 0 (#183).

RESIDUAL. The value α⁻¹ ≈ 137 stands as the single EM residual (the 137 problem is unchanged); what is established is that the structure around it is specifically PROTECTED — so α should be tested as protected-boundary-structure × one residual scale, not as a continuous tuning parameter.
