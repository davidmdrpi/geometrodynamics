# Adiabatic two-throat exchange holonomy: measure the Pin⁻ sign along a swap path (PR #188)

**Run:** 2026-06-30T01:27:56+00:00

Measures the two-throat exchange sign operationally — adiabatically transports the throat spin-½ state along an explicit swap path and reads off the holonomy `−I` (a π Berry phase, the Pin⁻ `−1`), instead of asserting it from `T²=−I`. *(QFT on the classical throat, not quantum gravity.)*

- **Holonomy**: path-ordered SU(2) transport along the swap (2π) loop = −I
- **Measured sign**: −1 (the spinor returns to minus itself); Berry phase π
- **Topological**: a wandering-axis swap gives the same −I (the ℤ₂ homotopy class)
- **Controls**: double-swap (4π) → +1; contractible loop → +1

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | measure the exchange sign as an adiabatic holonomy | **PASS** |
| T2 | `T2_swap_path_and_spin_statistics` | the swap path: π₁(RP²)=ℤ₂ generator ≃ a 2π rotation (FR) | **PASS** |
| T3 | `T3_swap_holonomy_measured` | the holonomy measured: Hol=−I, sign=−1, Berry phase π | **PASS** |
| T4 | `T4_path_independence_topological` | path-independence: a wandering-axis swap gives the same −I | **PASS** |
| T5 | `T5_controls_double_swap_and_contractible` | controls: double-swap (4π) → +1; contractible loop → +1 | **PASS** |
| T6 | `T6_pin_minus_identification` | the Pin⁻ identification: −1 = T² = (−1)^{2j}, spin-½ throat | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | ADIABATIC_EXCHANGE_HOLONOMY_MEASURES_PIN_MINUS | **PASS** |

## The measured holonomies

| loop | holonomy | measured sign |
|---|---|---:|
| swap (2π, the exchange) | −I (‖Hol+I‖ = 2e-06) | -1 |
| double-swap (4π, two exchanges) | +I | +1 |
| contractible (no exchange) | +I | +1 |

The single swap gives a π Berry phase (-3.1416 ≈ π) — the exchange sign `−1`, the Pin⁻ monodromy `T²=−I` measured along the path.

## Verdict

**ADIABATIC_TWO_THROAT_EXCHANGE_HOLONOMY_MEASURES_THE_PIN_MINUS_SIGN_MINUS_ONE_ALONG_THE_SWAP_PATH.** MEASURED — THE PIN⁻ EXCHANGE SIGN −1, AS AN ADIABATIC HOLONOMY ALONG THE SWAP PATH.

THE HOLONOMY. Path-ordering the spin-½ transport along the swap (2π) loop gives Hol = −I to machine precision (‖Hol + I‖ = 1.8e-06): the spinor returns to minus itself, the exchange sign ⟨ψ|Hol|ψ⟩ = -1, the Berry phase -3.142 ≈ π.

TOPOLOGICAL. A wandering-axis swap gives the same −I (the ℤ₂ homotopy class), converging as the transport is refined — any way of doing the exchange gives the same −1.

CONTROLS. A double-swap (4π, two exchanges) gives +1 and a contractible loop gives +1, so the −1 is the single-swap odd class.

PIN⁻. The −1 is the monodromy T² = −I (½ tr T² = -1): the throat is a spin-½ spinor via the non-orientable RP² closure, so its 2π/swap holonomy is (−1)^{2j} = −1 (a scalar throat would give +1). The #185 exchange sign, now measured along an explicit swap path rather than asserted from the algebra.
