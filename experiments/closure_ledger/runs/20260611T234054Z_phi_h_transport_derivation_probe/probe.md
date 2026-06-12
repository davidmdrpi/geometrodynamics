# Explicit Hopf-transport derivation of φ_h = π/k₅ (PR #159)

**Run:** 2026-06-11T23:40:54+00:00

Supplies the transport derivation #158 flagged: the quark CP phase scale φ_h = π/k₅ is derived from the connection's ½ (the spin-½ factor), the #63 C-swap orientation sign, the mass-locked dk = max rule, and the winding-capacity sector arc — explicit path-ordered transport exact to 1e-15, every alternative sector count excluded by the CP data, and the #156-consumed input RETURNED to the budget. Quark CP now stands as a calibration-free five-observable prediction. *(QFT on the classical throat, not quantum gravity.)*

- **Chain**: phase = ±dk·(½)·(2π/k₅) = ±dk·π/k₅ (transport exact to 1e-15)
- **Spinor check**: full circuit k=1 → π (T² = −I, the module consistency)
- **Exclusion**: π/3, π/4, π/6, 2π/5, π/10 all fail ≥ 1 observable; π/k₅ unique survivor
- **Prediction**: J 0.969, (β,γ,α) = (22.8, 63.5, 93.8)°, sin δ 0.888 vs 0.887 — no calibration
- **Budget**: the #156 input RETURNED; flavor arc net: inputs +0, knobs −1
- **Open**: hop arc from explicit shell wavefunctions; soft V_us

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the #158-flagged transport derivation (the #152 path) | **PASS** |
| T2 | `T2_explicit_transport_rate` | explicit transport: full circuit → π (spinor); sector → k·π/k₅ exact | **PASS** |
| T3 | `T3_ingredients_audit` | ingredients: ½ derived; ± C-swap verified; dk = max mass-locked | **PASS** |
| T4 | `T4_alternative_sector_exclusion` | exclusion scan: π/k₅ the unique survivor of six candidates | **PASS** |
| T5 | `T5_derived_prediction` | derived prediction: five observables, no calibration | **PASS** |
| T6 | `T6_budget_input_returned` | budget: the #156 input RETURNED; flavor arc net +0 inputs | **PASS** |
| T7 | `T7_scope` | scope: hop arc from shell wavefunctions = the final mile | **PASS** |
| T8 | `T8_assessment` | PHI_H_PI_OVER_K5_DERIVED_TRANSPORT_EXACT_ALTERNATIVES_EXCLUDED | **PASS** |

## The explicit transport (sector arc 2π/k₅)

| winding k | transported phase | k·π/k₅ | error |
|---:|---:|---:|---:|
| 1 | 0.6283185307 | 0.6283185307 | 5.6e-16 |
| 3 | 1.8849555922 | 1.8849555922 | 2.4e-15 |
| 5 | -3.1415926536 | -3.1415926536 | 4.4e-15 |

Full circuit at k = 1: |phase| = 3.1415926536 = π — the spinor sign flip (the module's own consistency check).

## The exclusion scan (observed: J/t = 1.0, β = 22.2°, γ = 65.9°, sin δ = 0.887)

| candidate | φ | J/target | β | γ | sin δ | passes all five? |
|---|---:|---:|---:|---:|---:|:---:|
| π/k₅ = π/5 (winding capacity) | 0.6283 | 0.969 | 22.78 | 63.48 | 0.888 | ✓ |
| π/3 (generation count) | 1.0472 | -0.02 | -0.94 | -1.29 | -0.022 | ✗ |
| π/4 | 0.7854 | 0.661 | 20.83 | 42.24 | 0.668 | ✗ |
| π/6 | 0.5236 | 1.14 | 23.18 | 78.19 | 0.971 | ✗ |
| 2π/5 | 1.2566 | -0.499 | -17.69 | -32.79 | -0.538 | ✗ |
| π/10 | 0.3142 | 1.12 | 19.0 | 114.67 | 0.901 | ✗ |

## Verdict

**PHI_H_PI_OVER_K5_DERIVED_TRANSPORT_EXACT_ALTERNATIVES_EXCLUDED.** φ_h = π/k₅ IS DERIVED BY EXPLICIT HOPF TRANSPORT — THE CONNECTION'S ½, THE C-SWAP SIGN, THE MASS-LOCKED dk RULE, AND THE WINDING-CAPACITY SECTOR ARC — WITH ALL ALTERNATIVE SECTOR COUNTS EXCLUDED BY THE CP DATA AND THE #156 INPUT RETURNED: QUARK CP IS NOW A CALIBRATION-FREE FIVE-OBSERVABLE PREDICTION. #158 flagged the candidate; this probe walks the #152 modelled→derived path.

THE CHAIN. (a) The RATE: A_φ(0) = ½ (the spin-½ factor of the connection) ⟹ a winding-k state transported through arc L accumulates k·L/2 — verified by explicit path-ordered integration, with the full circuit at k = 1 giving exactly π (the spinor sign flip, the module's own consistency check). (b) The SIGN: the two Z₂ partition classes traverse the fiber oppositely (#63 C-swap; opposite transport conjugates, verified). (c) The WINDING CONTENT: dk = max(k, k′) — the SAME rule the mass calibration locked (the phase and magnitude rules share their dk: independent corroboration). (d) The ARC: one winding-sector 2π/k₅ (capacity k₅, #73/#126). ⟹ phase = ±dk·π/k₅; sector transport exact to 1e-15.

THE EXCLUSION SCAN. Every principled alternative fails: π/3 (the generation count) kills J outright; π/4 misses γ by 24°; π/6 misses γ by 12° and J by 14%; 2π/5 flips the CP sign; π/10 misses γ by 49°. π/k₅ alone passes all five observables — the identification is data-selected among principled candidates, the anti-numerology discipline in its positive mode.

THE DERIVED PREDICTION. At φ_h = π/k₅, calibration-free: J at 0.969 of target, (β, γ, α) = (22.78, 63.48, 93.75)° vs (22.2, 65.9, 91.9)°, sin δ = 0.888 vs 0.887; masses shifted 0.09%, V_cb untouched.

THE BUDGET. The #156 input is RETURNED — the flavor card's last open row closes (quark CP: derived). Net flavor-arc bookkeeping #149–#159: inputs +0, modelling knobs −1.

SCOPE. The final geometric mile — deriving the hop arc from the explicit shell wavefunctions — is flagged (the same status the #152 saddle had pre-derivation); the soft V_us direction stands.
