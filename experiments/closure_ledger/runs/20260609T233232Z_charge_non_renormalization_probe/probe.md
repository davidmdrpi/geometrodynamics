# Charge non-renormalization: Z₁ = Z₂ on the antipodal cavity (PR #145)

**Run:** 2026-06-09T23:32:32+00:00

Closes the renormalization-constant triangle opened by #136 (Z₂) and #144 (Z₃): the Ward identity Z₁ = Z₂ is computed on the antipodal cavity — the q = 0 photon insertion equals minus the self-energy slope term by term, so the dressed charge equals the bare Hopf charge c₁ exactly and universally across matter species. Only the universal photon √Z₃ (the #144 vacuum polarisation) renormalizes charge, and the protection is exact charge conservation at the unitary antipodal throat. *(QFT on the classical throat, not quantum gravity.)*

- **Triangle**: Z₂ (#136 Σ) · Z₁ (#141/#142 vertex) · Z₃ (#144 Π); e = (Z₂/Z₁)·√Z₃·e₀
- **Ward**: Λ(0) = −c₁Σ' (term by term); Σ' = -0.014670, Z₂ = 0.985543
- **Dressed charge**: F(0) = Z₂(c₁ + Λ) = c₁ exactly; universal across species
- **Consequence**: e = √Z₃·e₀ — the running of α is purely vacuum polarisation (#144)
- **Protection**: exact charge conservation at the unitary throat (#58/#141/#142)
- **Open**: q ≠ 0 form factors; higher loops; normalisation (#133); flavor (#134)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | close the renormalization triangle: compute Z₁ = Z₂ | **PASS** |
| T2 | `T2_charged_model` | charged χ (Dirichlet) × neutral φ (Neumann); pole below threshold | **PASS** |
| T3 | `T3_z2_wavefunction_renormalization` | Z₂ < 1 computed two ways (spectral sum vs finite difference) | **PASS** |
| T4 | `T4_ward_identity_z1_eq_z2` | Λ(0) = −c₁Σ′ term by term ⟹ Z₁ = Z₂ (machine precision) | **PASS** |
| T5 | `T5_dressed_charge_exact_and_universal` | F(0) = c₁ exact + universal (Z₂ varies, charge does not) | **PASS** |
| T6 | `T6_charge_violation_counterfactual` | charge-violating vertex ⟹ F(0) ≠ c₁ — conservation is the protection | **PASS** |
| T7 | `T7_ledger_and_scope` | ledger: derived identity/universality vs modelled coupling vs input α | **PASS** |
| T8 | `T8_assessment` | CHARGE_NON_RENORMALIZATION_Z1_EQ_Z2_DRESSED_CHARGE_EXACT_UNIVERSAL | **PASS** |

## The Ward identity, computed

| quantity | value |
|---|---:|
| Σ′(s₀) — spectral sum | -0.01466957 |
| Σ′(s₀) — finite difference | -0.01466957 |
| Z₂ = 1/(1 − Σ′) | 0.985543 |
| Λ(0) — q=0 charged-line insertion | 0.01466957 |
| −c₁Σ′(s₀) | 0.01466957 |
| residual Λ + c₁Σ′ (vs spectral) | 0.0 |
| residual Λ + c₁Σ′ (vs finite diff) | -1.49e-12 |

## The dressed charge is exact and universal

| species | Z₂ | F(0) − c₁ |
|---|---:|---:|
| χ(l=1)·φ(l=0), g=1.0 | 0.985543 | -1.11e-16 |
| χ(l=3)·φ(l=0), g=1.0 | 0.987631 | 0.0 |
| χ(l=1)·φ(l=2), g=0.7 | 0.995636 | 0.0 |
| χ(l=1)·φ(l=0), g=0.5 | 0.996346 | 0.0 |

## Counterfactual: a charge-violating vertex shifts the charge

| internal charge c′ | F(0) − c₁ |
|---:|---:|
| 0.8 | -0.002891 |
| 0.5 | -0.007229 |
| 0.0 | -0.01446 |

Z₂ varies by sector but the conserving dressed charge never moves; break charge conservation at the vertex and it does — the protection is the unitary throat's exact charge conservation (#58/#141/#142).

## Verdict

**CHARGE_NON_RENORMALIZATION_Z1_EQ_Z2_DRESSED_CHARGE_EXACT_UNIVERSAL.** THE WARD IDENTITY Z₁ = Z₂ HOLDS ON THE ANTIPODAL CAVITY TO MACHINE PRECISION: THE DRESSED CHARGE EQUALS THE BARE HOPF CHARGE c₁ EXACTLY AND UNIVERSALLY ACROSS SPECIES — ONLY THE PHOTON √Z₃ (#144) RENORMALIZES CHARGE, AND THE PROTECTION IS EXACT CHARGE CONSERVATION AT THE UNITARY THROAT. PR #144 completed the one-loop two-point sector (Z₂ from #136, Z₃ from #144); this probe closes the triangle by computing Z₁ = Z₂.

THE CHARGED MODEL. A charged field χ (the odd-l Dirichlet tower, Hopf charge c₁ = 1, #129/#58) couples to a neutral field φ (the even-l Neumann tower) through the charge-conserving cubic triple-overlap vertex (#136/#137/#141); the external pole sits below the lowest threshold, so the constants are real.

Z₂ ≠ 1, COMPUTED TWO WAYS. Σ'(s₀) from the analytic spectral sum and from a finite difference of the #136 self-energy agree to ~1e-9; Z₂ = 1/(1 − Σ') ≈ 0.986 at g = 1 — the dressed particle is genuinely renormalized, so there is something real to cancel.

THE WARD IDENTITY, COMPUTED. The q = 0 photon insertion on the internal charged line doubles the propagator: Λ(0) = Σ c₁|g|²/(s₀ − s_nm)² = −c₁Σ'(s₀) term by term (machine precision; the neutral line contributes zero). This is the #142 Ward–Takahashi identity Γ(p,p) = ∂S⁻¹/∂p at one loop — equivalently Z₁ = Z₂.

THE DRESSED CHARGE IS EXACT AND UNIVERSAL. F(0) = Z₂(c₁ + Λ(0)) = c₁ exactly (machine precision), across species with different towers and couplings: Z₂ varies by sector, F(0) never moves. Each sector's self-interaction cancels out of its own charge — why different generations (k ∈ {1,3,5}, #71) carry exactly the same |c₁| = 1 — and the charge renormalization collapses to e = √Z₃·e₀: the running of α is PURELY the #144 vacuum polarisation.

THE COUNTERFACTUAL. A charge-violating vertex (c' ≠ c₁) breaks the cancellation: F(0) ≠ c₁ and species-dependent. Charge non-renormalization rests on exact charge conservation at the throat — Σc₁ = 0 from the unitary antipodal mirror (#58/#141/#142); the absorbing throat leaks charge and loses the protection. The same single postulate as stable matter, gauge invariance, and the massless photon (#129/#130/#142/#144).

SCOPE. One loop on the fixed background; modelled cubic coupling (the #136 posture); the q = 0 charge-operator limit only. q ≠ 0 form factors, higher loops, the absolute normalisation (#133), and the flavor residuals (#134) stand; α(μ₀) stays the one EM input (#143).
