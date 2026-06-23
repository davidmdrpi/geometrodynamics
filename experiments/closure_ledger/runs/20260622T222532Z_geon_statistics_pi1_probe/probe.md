# π₁ of the two-mouth configuration space and FR-homotopy survival for the Pin⁻ throat (PR #171)

**Run:** 2026-06-22T22:25:32+00:00

Replaces the **orientable** Finkelstein–Rubinstein citation of #170 with the correct **geon-statistics** framework, and checks whether the −1 exchange sign survives the non-orientable Pin⁻ mouth. *(QFT on the classical throat, not quantum gravity.)*

- **Configuration space**: π₁: σ²=e (no 3D braiding) + R_i (spin) + τ_i (non-orientable reversal)
- **Spinorial**: 2π rotation = −1 (Pin⁻; Friedman–Sorkin)
- **New ingredient**: Pin⁻ reflection² = −1 (RP² admits Pin⁻ only; orientable FR lacks this)
- **Achiral**: non-orientable ⇒ own mirror image ⇒ handedness hypothesis met
- **Result**: exchange = −1 (Fermi) survives; conditional on the cited Dowker–Sorkin slide hypothesis

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | replace orientable FR with the geon-statistics check | **PASS** |
| T2 | `T2_pi1_two_mouth_configuration_space` | π₁ of the two-mouth config space (σ²=e, R_i, τ_i) | **PASS** |
| T3 | `T3_single_geon_spinorial` | single geon is spinorial (2π = −1, Pin⁻) | **PASS** |
| T4 | `T4_pin_minus_reflection_squared` | the new ingredient: Pin⁻ reflection² = −1 (not Pin⁺) | **PASS** |
| T5 | `T5_achirality_from_non_orientability` | non-orientable ⇒ achiral ⇒ a hypothesis met | **PASS** |
| T6 | `T6_exchange_sign_survives` | exchange = −1 (Fermi) survives, conditional on slide | **PASS** |
| T7 | `T7_honesty_and_scope` | honest scope (computed vs the cited geon theorem) | **PASS** |
| T8 | `T8_assessment` | FR_HOMOTOPY_SURVIVES_PIN_MINUS_GEON_FERMI_CONDITIONAL | **PASS** |

## The ingredient orientable FR lacks

| structure | reflection² | RP² admits it? |
|---|---:|---|
| Pin⁻ | -1 | yes (w₂+w₁²=0) |
| Pin⁺ | +1 | no (w₂≠0) |

(the non-orientable exchange carries a reflection; Pin⁻ makes it square to −1)

## Verdict

**FR_HOMOTOPY_SURVIVES_PIN_MINUS_MOUTH_GEON_FERMI_CONDITIONAL_ON_EXCHANGEABILITY.** THE −1 SURVIVES — ON THE CORRECT FOOTING. PR #170's orientable Finkelstein–Rubinstein citation is replaced by the geon-statistics framework, and the Pin⁻ mouth still gives Fermi statistics, conditional on one cited hypothesis.

THE CONFIGURATION SPACE. π₁ of the unordered two-mouth configuration space of an asymptotically-flat 3-slice has the exchange σ with σ² = e (≥3 dimensions ⇒ symmetric group, no braiding; only ±1 statistics), the per-geon 2π rotation R_i, and — because the mouth is non-orientable — an orientation-reversing loop τ_i.

SPINORIAL. The single geon's 2π rotation acts as −I (4π = +I): spinorial, the Pin⁻ holonomy of #170 and Friedman–Sorkin's spin-½ from gravity.

THE NEW INGREDIENT. The orientable FR argument never sees the orientation reversal a non-orientable exchange carries. RP² admits Pin⁻ only, and in Pin⁻ a reflection squares to −1 (computed -1; Pin⁺ would give +1). The reversal therefore contributes a consistent −1.

ACHIRALITY. Non-orientability makes the geon its own mirror image, meeting the geon spin–statistics theorem's handedness hypothesis automatically.

THE RESULT. Spinorial + Pin⁻ reflection² = −1 + achiral ⇒ the exchange is the −1 irrep of σ: a FERMION. CONDITIONAL on the Dowker–Sorkin exchangeability ('slide') hypothesis — holding for identical asymptotically-flat throats, cited not derived. The geon literature has spin–statistics violation examples, so this is a genuine check BAM passes (Pin⁻ + achiral + exchangeable), not an automatic result. The remaining honest gap is that hypothesis and the field-theory mapping class group, not the spinor sign or the reflection algebra.
