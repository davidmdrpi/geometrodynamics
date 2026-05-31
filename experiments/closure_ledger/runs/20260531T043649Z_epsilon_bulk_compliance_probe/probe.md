# Can ε be computed from bulk compliance, or only inferred from the meV scale? (PR #112)

**Run:** 2026-05-31T04:36:49+00:00

PRs #87–#90 INFERRED the chargeless-throat healing length `ε` by demanding the observed `S = ln(m_D/m_ν) ≈ 15–18` — `ε` read back from the meV scale it is meant to predict. This probe asks whether `ε` can instead be COMPUTED from the bulk compliance (the elastic/nucleation response of the throat neck), with no neutrino input. **Answer: a genuine partial — the smallness yes, the precise value no.**

- **Computed (meV-free)**: ε ~ R_c³ ~ 10⁻² (sub-throat); m_ν ≈ 2.1 meV output at t = k_5√(2π)
- **Derives**: the exponential smallness of m_ν (sub-throat ⟹ large S)
- **Residual**: the precise ε — m_ν ∝ ε^{4.8}, O(1) ambiguity spans 2–108 meV
- **Obstruction**: absolute compliance normalization = unpinned κ₅²/Λ₅ (the one anchor)
- **One line**: the smallness is derived from bulk compliance; the exact value is not

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_question_and_current_status` | ε currently inferred via S = ln(m_D/m_ν) | **PASS** |
| T2 | `T2_define_bulk_compliance` | compliance: ε = ℓ²/2rs, ℓ ~ R_c = 2σ/ρ (electron-calibrated) | **PASS** |
| T3 | `T3_compute_eps_from_compliance` | ε ~ R_c³ 0.011, Δ³ 0.018, R_c²/2 0.025 — sub-throat, meV-free | **PASS** |
| T4 | `T4_chain_closes_meV_free` | ε=R_c³ + t=k_5√(2π) ⟹ S≈16.85 ⟹ m_ν≈2.1 meV (retrodiction) | **PASS** |
| T5 | `T5_steep_sensitivity_precision_not_geometric` | m_ν ∝ ε^{4.8}: O(1) ambiguity spans 2–108 meV (not geometric) | **PASS** |
| T6 | `T6_obstruction_unpinned_bulk_normalization` | absolute compliance = unpinned κ₅²/Λ₅; only √6 fixed | **PASS** |
| T7 | `T7_reclassification_smallness_derived_value_residual` | smallness DERIVED (meV-free); precise value RESIDUAL | **PASS** |
| T8 | `T8_assessment` | EPSILON_BULK_COMPLIANCE_DERIVES_SMALLNESS_PRECISE_VALUE_STAYS_RESIDUAL | **PASS** |

## The chain closes meV-free (ε = R_c³, t = k_5√(2π))

`S ≈ 16.85` ⟹ predicted m_ν spectrum:

| gen | m_D (keV) | m_ν predicted (meV) |
|---|---:|---:|
| 1 | 43 | 2.08 |
| 2 | 80 | 3.86 |
| 3 | 118 | 5.7 |

The meV **scale** (2.1–5.7 meV) is an **output** of bulk geometry — no neutrino input. A sub-throat `ε ≪ 1` makes `L*` large ⟹ `S` large ⟹ `m_ν = m_D·e^{−S}` exponentially small: **the smallness is derived.** _Scale only_ — uniform `S` gives `m_ν ∝ m_D` (×2.7), so the observed generation spread (up to `√Δm²_31 = 50 meV`) is the separate PR #90/#91 residual, not addressed by `ε`.

## The catch: m_ν ∝ ε^{4.8} (precision not geometric)

| ε candidate | value | m_ν (meV) |
|---|---:|---:|
| R_c³ | 0.0110 | 2.1 |
| Δ³ | 0.0176 | 20.4 |
| R_c²/2 (ℓ=R_c) | 0.0247 | 107.8 |

`d ln m_ν/d ln ε ≈ 4.75` — a factor ~2 in `ε` spans m_ν by ~×51.3. Landing on 'few meV' still selects `ε ≈ R_c³` via the observed scale.

## Verdict

**EPSILON_BULK_COMPLIANCE_DERIVES_SMALLNESS_PRECISE_VALUE_STAYS_RESIDUAL.** BULK COMPLIANCE DERIVES THE SMALLNESS OF ε (AND HENCE m_ν), BUT NOT ITS PRECISE VALUE. PRs #87–#90 inferred the chargeless-throat healing length ε by demanding the observed S = ln(m_D/m_ν) ≈ 15–18 — ε read back from the meV scale it is meant to predict. This probe asks whether ε can instead be COMPUTED from the bulk compliance, with no neutrino input. The answer is a genuine partial: the order of magnitude yes, the precise value no.

THE COMPLIANCE COMPUTATION (meV-free). The compliance cutoff is a sub-throat healing length ε = ℓ²/(2rs), with ℓ ~ R_c = 2σ/ρ the critical-bubble scale. σ = 1/(12π) and ρ = 3/(4π) are fixed by the ELECTRON rest-energy calibration (PR #58), so R_c = 2/9 ≈ 0.222 is a pure number from the charged-throat geometry — no neutrino mass anywhere. The candidate compliances are all sub-throat, O(10⁻²): R_c³ ≈ 0.011, Δ³ ≈ 0.018, R_c²/2 ≈ 0.025. So ε's order of magnitude and its decisive sub-throat character ARE computable from bulk compliance, independent of the meV scale.

THE CHAIN CLOSES meV-FREE. With ε = R_c³ and the winding-edge tension t = k_5√(2π) ≈ 12.53 (the PR #89 closure quantity, also meV-free), S = t²·P0·L*(R_c³) ≈ 16.85, so m_ν(gen 1) = m_D·e^{−S} ≈ 2.1 meV with m_D ≈ 43 keV the electron-anchored cavity-floor Dirac mass. The meV scale comes out as an OUTPUT — a genuine retrodiction. This structurally DERIVES the neutrino's lightness: a sub-throat healing length (ε ≪ 1) makes L* large, hence S large, hence m_ν exponentially small.

THE CATCH: PRECISION IS NOT GEOMETRIC. The bounce action is steep in ε — d ln m_ν/d ln ε = t²·P0·rs/2 ≈ 4.8 — so the O(1) ambiguity among the healing-length candidates blows up: R_c³ → 2 meV, Δ³ → 20 meV, R_c²/2 → 108 meV, a factor ~50 spread from a factor ~2 in ε. Landing precisely on "few meV" still SELECTS ε ≈ R_c³ using the observed scale; there is no first-principles reason to prefer R_c³ over R_c²/2.

THE OBSTRUCTION. A true compliance = 1/stiffness needs the absolute bulk gravitational stiffness (κ₅², Λ₅). BAM fixes only the dimensionless RS-tuning combination λ_crit κ₅²/√|Λ₅| = √6 (PR #57); κ₅² and Λ₅ separately are absorbed into the single dimensionful anchor and are not pinned. So ε is computable only up to an O(1) factor — the same residual status everything tied to the one anchor shares.

VERDICT. ε is upgraded from "inferred from the meV scale" to "bulk-geometric to order of magnitude": the SMALLNESS is derived from the electron-calibrated nucleation compliance (meV-independent), and m_ν ~ meV is a retrodiction at the winding-edge tension. The PRECISE value stays a residual — m_ν ∝ ε^{4.8}, the O(1) ambiguity spans 2–108 meV, and the absolute normalization is the unpinned κ₅²/Λ₅. The smallness is derived from bulk compliance; the exact value is not.

## What this changes, and what stays open

- **Changed:** `ε` moves from "inferred from the meV scale" to "bulk-geometric to order of magnitude." The neutrino's exponential lightness is now DERIVED from the electron-calibrated nucleation compliance (sub-throat `ε`), and `m_ν ~ meV` is a retrodiction — not an input.
- **Open:** the PRECISE `ε`. Because `m_ν ∝ ε^{4.8}`, the O(1) ambiguity (`R_c³`/`Δ³`/`R_c²/2`) spans `m_ν ~ 2–108 meV`; the absolute compliance normalization is the unpinned `κ₅²/Λ₅` (the one anchor). The exact value is still residual.
