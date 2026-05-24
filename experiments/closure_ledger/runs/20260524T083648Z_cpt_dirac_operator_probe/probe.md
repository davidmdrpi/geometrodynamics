# Explicit CPT / Dirac-spinor operator on the BAM throat

**Run:** 2026-05-24T08:36:48+00:00

Builds the explicit gamma-matrix operators C, P, T on the throat Dirac 4-spinor (the inner/outer mouth doubling) and the composite Θ = CPT, closing the remaining note of PR #64. Result: Θ ∝ γ⁵ and Θ² = −I.

- **Spinor**: throat Dirac 4-spinor = inner/outer mouth doubling
- **C**: `iγ²γ⁰ (defining C⁻¹γ^μC=−(γ^μ)ᵀ; the #63 inner/outer swap)`
- **P**: `γ⁰ (P²=+I)`
- **T**: `γ¹γ³ K (T²=−I; the B2 iσ_y)`
- **Θ = CPT**: `γ⁰γ¹γ²γ³=−iγ⁵ (∝ γ⁵); Θ_m²=−I, antiunitary Θ²=+I`
- **Consistency**: Θ⁻¹γ^μΘ=−γ^μ → j^μ→−j^μ (the #64 sign table)
- **B4 caveat**: dimensionless constant matrices; group facts; scale-independent

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_throat_dirac_4spinor` | Clifford {γ^μ,γ^ν}=2η; γ⁵ mixes the mouths | **PASS** |
| T2 | `T2_charge_conjugation` | C=iγ²γ⁰: C⁻¹γ^μC=−(γ^μ)ᵀ, C²=−I (#63) | **PASS** |
| T3 | `T3_parity` | P=γ⁰: γ⁰γⁱγ⁰=−γⁱ, P²=+I | **PASS** |
| T4 | `T4_time_reversal` | T=γ¹γ³K: T²=−I (B2 iσ_y) | **PASS** |
| T5 | `T5_theta_is_gamma5` | Θ=γ⁰γ¹γ²γ³=−iγ⁵ (∝ γ⁵) | **PASS** |
| T6 | `T6_theta_squared` | Θ_m²=−I; antiunitary Θ²=+I (−1 = T², B2) | **PASS** |
| T7 | `T7_consistency_falsification_b4` | Θ⁻¹γ^μΘ=−γ^μ → j^μ→−j^μ (#64) | **PASS** |
| T8 | `T8_assessment` | throat = Dirac spinor with standard CPT | **PASS** |

## T2: C = iγ²γ⁰ (charge conjugation)

- defining relation C⁻¹γ^μC = −(γ^μ)ᵀ: True
- matrix square C_m² = −I: True

## T3: P = γ⁰ (parity)

- spatial γⁱ flip (γ⁰γⁱγ⁰=−γⁱ): True; temporal γ⁰ fixed: True; P²=+I: True

## T4: T = iγ¹γ³ K (time reversal)

- T² = −I: True (the B2 iσ_y signature: True)

## T5: Θ = CPT ∝ γ⁵

- Θ_m = γ⁰γ¹γ²γ³ (total spacetime inversion) = −iγ⁵: True
- anticommutes with every γ^μ (∝ γ⁵): True

## T6: Θ² (matrix −I; operator +I)

- matrix Θ_m² = −I: True
- antiunitary operator Θ² = Θ_m Θ_m* = +I: True ((CPT)²=+1, a symmetry)
- fermionic −1 double cover is T²=−I (B2): True

## T7: Consistency / falsification / B4

- Θ⁻¹γ^μΘ = −γ^μ (current j^μ → −j^μ, recovering #64): True
- γ⁵ anticommutes with every γ^μ: True
- C ↔ inner/outer swap (#63): True; T ↔ B2 iσ_y: True

## T8: Assessment

- C: iγ²γ⁰ (C⁻¹γ^μC=−(γ^μ)ᵀ, C_m²=−I); P: γ⁰ (P²=+I); T: γ¹γ³ K (T²=−I)
- Θ = CPT: γ⁰γ¹γ²γ³=−iγ⁵ (∝ γ⁵); Θ² = Θ_m²=−I (matrix); antiunitary Θ²=+I; fermionic −1 = T²=−I (B2)
- remaining: the throat spinor as an S_BAM bulk solution; P vs antipodal Z₂ at the spinor level

## Verdict

**CPT_OPERATOR_CONSTRUCTED.** CPT OPERATOR CONSTRUCTED. The explicit gamma-matrix CPT operator on the throat Dirac spinor is built, closing the remaining note of PR #64.

THE SPINOR. The throat Dirac 4-spinor is the doubling of the inner/outer mouth 2-spinors, Ψ=(Ψ_inner,Ψ_outer), with the Dirac-representation Clifford algebra and γ⁵=[[0,I],[I,0]] mixing the mouths.

THE OPERATORS. C = iγ²γ⁰ (charge conjugation) satisfies the defining relation C_m⁻¹γ^μC_m=−(γ^μ)ᵀ with C²=−I — the inner/outer swap (#63, c₁→−c₁) on the spinor. P = γ⁰ (parity, γ⁰γⁱγ⁰=−γⁱ, P²=+I). T = iγ¹γ³ K (time reversal, antiunitary) with T²=−I — the fermionic signature, the same T²=−1 as the B2 iσ_y spin structure.

THE COMPOSITE. CPT inverts all four spacetime axes (P inverts space, T inverts time), so its spinor matrix is the total-inversion product Θ_m = γ⁰γ¹γ²γ³ = −iγ⁵ — proportional to the chiral matrix γ⁵, characterized by anticommuting with every γ^μ. The matrix squares to Θ_m²=(−iγ⁵)²=−I, but the CPT OPERATOR is antiunitary (C, P unitary; T antiunitary), so its square is Θ²=Θ_m Θ_m*=+I — (CPT)²=+1, consistent with CPT being a symmetry; the fermionic −1 spinor double cover is carried by T²=−I (the 2π rotation / RP³ spin structure, B2), not by CPT². On the 4-current j^μ=Ψ̄γ^μΨ, Θ_m⁻¹γ^μΘ_m=−γ^μ (γ⁵ anticommutes with every γ^μ), realizing j^μ(x)→−j^μ(−x) — the #64 CPT sign table at the operator level.

So the throat is a Dirac spinor carrying the standard CPT operator, realized geometrically: the inner/outer mouth doubling, C = the #63 swap, T = the B2 iσ_y. B4: the operators are dimensionless constant matrices (Θ∝γ⁵, T²=−I are group facts) — scale-independent. Remaining: the throat spinor as an explicit S_BAM bulk solution, and disentangling P from the antipodal Z₂ at the spinor level.

## What this leaves open

- **The throat spinor from S_BAM.** The explicit Dirac spinor as a bulk solution of the action, with the mouth doubling derived rather than posited.
- **P vs the antipodal Z₂.** Disentangling spatial parity from the RP³ deck transformation (B2) at the spinor level.
