# Probe C — bulk extension, APS inflow, Pin⁻ mouth (PR #230)

_Run 2026-07-28T23:23:21.190117+00:00 · 5/5 PASS_

## APS per boundary spin structure

| `ε` | `η(0)` | `h` | boundary defect `−(h+η)/2` |
|---:|---:|---:|---:|
| +1 | +0.2500 | 0 | -0.1250 |
| -1 | -0.2500 | 0 | +0.1250 |

**Ω₃^Spin = 0 ⟹ both extend.** Inflow does **not** select — both sectors remain admissible.

## The Pin⁻ mouth: ABK in Z₈

| mouth A | mouth B | total mod 8 | bounds? |
|---:|---:|---:|---|
| +1 | +1 | 2 | no |
| +1 | -1 | 0 | **yes** |
| -1 | +1 | 0 | **yes** |
| -1 | -1 | 6 | no |

Only opposite Pin⁻ structures bound — deriving the repo's `ThroatDefect` assertion. But it is a *relative* rule: it does not fix the absolute sector.

## Verdict

**BOTH_SPIN_STRUCTURES_EXTEND_SO_INFLOW_DOES_NOT_SELECT_BUT_ABK_FORCES_THE_TWO_MOUTHS_OF_A_THROAT_TO_CARRY_OPPOSITE_PIN_MINUS_STRUCTURES**
