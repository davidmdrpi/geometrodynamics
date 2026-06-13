# The v4 library migration: the flavor-CP lock lands in geometrodynamics/qcd (PR #164)

**Run:** 2026-06-13T08:14:26+00:00

Performs the migration PR #163 staged: the v4 flavor-CP lock moves from probe-local code into the calibrated library, additively over the FROZEN v3 lock. At the φ_h = 0 default the v3 surface is bit-for-bit reproducible (real Hamiltonian, real CKM, no CP), so every PR #155–#162 probe is untouched. The v4 lock inherits the v3 mass spectrum exactly (the Hopf holonomy is a pure mixing phase — the #158 relocation) and realizes the complete nine-observable flavor-CP dataset at ≤ 1% from `extract_ckm_matrix` at the derived φ_h = π/k₅. *(QFT on the classical throat, not quantum gravity.)*

- **The surface**: 6 new QuarkParams fields + extract_ckm_matrix + LOCKED_QUARK_PARAMS_V4 (φ_h = π/k₅), all default-off
- **Default-off**: v3 lock real H + real CKM (no CP); #155–#162 untouched
- **Mass inheritance**: holonomy stripped; v4 masses ≡ v3 (1e-9)
- **The nine observables**: ≤ 1% from the library at φ_h = π/k₅, unitary
- **Counting**: +3 params, +5 independent observables (net +2); CP at 0 params

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | perform the staged migration into the library | **PASS** |
| T2 | `T2_public_surface` | new public surface: fields + extract_ckm_matrix + v4 lock | **PASS** |
| T3 | `T3_default_off_v3_bit_reproducible` | default-off: v3 bit-reproducible (real H, real CKM) | **PASS** |
| T4 | `T4_mass_inheritance` | mass inheritance: v4 masses ≡ v3 (holonomy stripped) | **PASS** |
| T5 | `T5_nine_observables` | nine observables from the library, ≤ 1%, unitary, at π/k₅ | **PASS** |
| T6 | `T6_counting_and_regression` | counting (net +2) + additive regression decision | **PASS** |
| T7 | `T7_program_position` | flavor-CP sector closed in library; residuals remain | **PASS** |
| T8 | `T8_assessment` | V4_LOCK_MIGRATED_LIBRARY_DEFAULT_OFF_NINE_OBSERVABLES | **PASS** |

## The nine observables, from the library

| observable | ratio to observed |
|---|---:|
| V_us | ×1.0025 |
| V_cb | ×1.0024 |
| V_ub | ×0.9976 |
| V_td | ×1.008 |
| V_ts | ×1.0021 |
| J | ×1.0044 |
| (β, γ, α) | (22.29, 65.9, 91.81)° vs (22.2, 65.9, 91.9)° |
| sin δ | 0.889 vs 0.887 |
| unitarity | dev 6.7e-16 |

## Verdict

**V4_LOCK_MIGRATED_LIBRARY_DEFAULT_OFF_MASSES_INHERITED_NINE_OBSERVABLES.** THE v4 FLAVOR-CP LOCK IS MIGRATED INTO THE LIBRARY, ADDITIVELY OVER THE FROZEN v3 LOCK: v3 IS BIT-FOR-BIT REPRODUCIBLE, v4 INHERITS THE MASSES EXACTLY AND REALIZES THE NINE OBSERVABLES AT ≤ 1% AT THE DERIVED φ_h = π/k₅.

THE SURFACE. Six new QuarkParams fields (the Hopf phase φ_h, three targeted couplings, the per-shell diagonal retunes), the extract_ckm_matrix() reader, and LOCKED_QUARK_PARAMS_V4 at φ_h = π/k₅ — all default-off.

DEFAULT-OFF. At φ_h = 0 the v3 Hamiltonian is exactly real and the CKM is a real rotation with no CP; every PR #155–#162 probe pins to the frozen v3 lock and is untouched. The migration is additive, not a re-baseline.

MASS INHERITANCE. The holonomy is a pure mixing phase (the #158 relocation): extract_physical_spectrum strips φ_h, so the v4 lock inherits the v3 spectrum to 3e-09.

THE NINE OBSERVABLES. extract_ckm_matrix at the v4 lock returns the magnitudes and J at ≤ 1%, (β, γ, α) = (22.29, 65.9, 91.81)° and sin δ = 0.889, unitarily — the complete flavor state as a library call.

THE COUNTING. Unchanged from #163: +3 parameters for +5 independent observables (net +2); the CP sector costs zero parameters; the #150 budget is unchanged. The four-step migration is complete.
