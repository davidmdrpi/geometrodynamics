# Explicit CPT / Dirac-spinor operator on the BAM throat

Closes the "remaining" note of the CPT-assembly probe (PR #64), which
assembled CPT at the level of the observable sign-table plus the single
spinor fact `T² = −I`. This probe builds the **explicit gamma-matrix
operators** `C, P, T` on the throat Dirac 4-spinor and the composite
`Θ = CPT`, computing `Θ` and `Θ²` directly. The clean result: `Θ ∝ γ⁵`
(the chiral matrix) and `Θ² = −I` (the fermionic double cover).

## The throat Dirac 4-spinor

The wormhole has two mouths — inner (`r < R_MID`) and outer (`r > R_MID`)
— each carrying a 2-component Hopf/SU(2) spinor (`T = iσ_y`, B2). The
throat Dirac spinor is their **doubling**,

```
Ψ = (Ψ_inner, Ψ_outer)ᵀ ∈ ℂ⁴ ,
```

with the Dirac gamma matrices in the Dirac representation,
`γ⁰ = diag(I, −I)`, `γⁱ = [[0, σⁱ],[−σⁱ, 0]]`, `γ⁵ = iγ⁰γ¹γ²γ³ =
[[0, I],[I, 0]]`. The inner/outer mouth structure is the upper/lower
Dirac block; the inner/outer swap (PR #63) acts on it.

## The discrete operators (explicit matrices)

  - **C — charge conjugation.** `Ψ → Ψ^c = C_m Ψ*`, `C_m = iγ²γ⁰`. It
    satisfies the defining relation `C_m⁻¹ γ^μ C_m = −(γ^μ)ᵀ` and
    `C² = −I`. This is the inner/outer swap (PR #63, `c₁ → −c₁`) realized
    on the spinor.

  - **P — parity.** `Ψ(t, x) → γ⁰ Ψ(t, −x)`, with
    `γ⁰ γⁱ γ⁰ = −γⁱ`, `γ⁰ γ⁰ γ⁰ = γ⁰`; `P² = +I`.

  - **T — time reversal.** `Ψ(t, x) → T_m Ψ*(−t, x)`, `T_m = γ¹γ³`
    (antiunitary, `K` = complex conjugation); `T² = −I` — the fermionic
    signature, the same `T² = −1` as the B2 `iσ_y` spin structure.

## The composite Θ = CPT

CPT inverts **all four** spacetime axes (P inverts space, T inverts
time), so its spinor matrix is the total-inversion product

```
Θ_m = γ⁰γ¹γ²γ³ = −i γ⁵ ,
```

which is **`∝ γ⁵`** (the chiral matrix), characterized by anticommuting
with every `γ^μ`. The matrix squares to `Θ_m² = (−iγ⁵)² = −I`, **but** the
CPT *operator* is antiunitary (C, P unitary; T antiunitary — one
antiunitary factor), so its square is

```
Θ² = Θ_m Θ_m* = +I ,
```

i.e. `(CPT)² = +1` on the spinor, consistent with CPT being a symmetry.
The fermionic `−1` spinor double cover is carried by `T² = −I` (the `2π`
rotation / the RP³ spin structure, B2), **not** by `(CPT)²`. On the
4-current `j^μ = Ψ̄γ^μΨ`, `Θ_m⁻¹ γ^μ Θ_m = −γ^μ` (γ⁵ anticommutes with
every `γ^μ`), realizing `j^μ(x) → −j^μ(−x)` — the #64 CPT sign table at
the operator level.

(Note: naively multiplying the matrix parts `C_m P_m T_m` gives `γ⁰γ⁵`,
not `γ⁵`, because that product mishandles the field conjugations of the
antiunitary factors. The correct CPT matrix is fixed by its defining
property — inverting all four axes / anticommuting with every `γ^μ` —
which is the total-inversion product `γ⁰γ¹γ²γ³ ∝ γ⁵`.)

## B4 accounting

The operators are **dimensionless constant matrices**; `Θ ∝ γ⁵`,
`Θ² = −I`, `C² = −I`, `T² = −I` are group facts. The spinor structure is
geometric — independent of the single anchor `m_e` (B4-consistent).

## Tests

  T1. **Throat Dirac 4-spinor.** Inner/outer mouth doubling; the Dirac
      gammas (Dirac rep), `γ⁵ = [[0,I],[I,0]]`.
  T2. **C = iγ²γ⁰.** `C_m⁻¹ γ^μ C_m = −(γ^μ)ᵀ` (the robust defining
      relation), matrix `C_m² = −I` (the #63 inner/outer swap on the
      spinor).
  T3. **P = γ⁰.** `γ⁰ γⁱ γ⁰ = −γⁱ`, `P² = +I`.
  T4. **T = γ¹γ³ K.** `T² = −I` (fermionic; the B2 `iσ_y` structure).
  T5. **Θ = CPT ∝ γ⁵.** `Θ_m = γ⁰γ¹γ²γ³ = −iγ⁵` (total spacetime
      inversion; anticommutes with every `γ^μ`).
  T6. **Θ² (matrix −I; operator +I).** `Θ_m² = −I`; antiunitary
      `Θ² = +I` ((CPT)²=+1, a symmetry); the fermionic `−1` is `T² = −I`
      (B2).
  T7. **Consistency / falsification / B4.** `Θ_m⁻¹ γ^μ Θ_m = −γ^μ`
      (the current `j^μ → −j^μ`, recovering #64); C ↔ #63, T ↔ B2; the
      operators are dimensionless.
  T8. **Assessment.**

## Verdict structure

  - **CPT_OPERATOR_CONSTRUCTED** (expected): the explicit gamma-matrix
    `C = iγ²γ⁰` (`C⁻¹γ^μC = −(γ^μ)ᵀ`, `C_m² = −I`), `P = γ⁰` (`P² = +I`),
    `T = γ¹γ³ K` (`T² = −I`) give the CPT operator `Θ = γ⁰γ¹γ²γ³ = −iγ⁵`
    (`Θ ∝ γ⁵`, total spacetime inversion), with matrix `Θ_m² = −I` but
    antiunitary `Θ² = +I` ((CPT)²=+1, a symmetry; the fermionic `−1` is
    `T² = −I`, B2) and `Θ⁻¹γ^μΘ = −γ^μ` (the #64 current sign at the
    operator level). The throat is a Dirac spinor with the standard CPT
    operator, realized geometrically (the inner/outer mouth doubling,
    C = the #63 swap, T = the B2 `iσ_y`).

  - **OPERATOR_INCONSISTENT**: the defining relations fail, `Θ` is not
    `∝ γ⁵`, or `Θ² ≠ ±I`.

## What this leaves open

  - **The throat spinor from `S_BAM`.** The explicit Dirac spinor as a
    bulk solution of the action (not just the algebra), with the mouth
    doubling derived rather than posited.
  - **P vs the antipodal `Z₂`.** Disentangling spatial parity from the
    RP³ deck transformation (B2) at the spinor level.

## Cross-references

  - `docs/cpt_assembly_research_plan.md` — the CPT sign table (#64).
  - `docs/charge_conjugation_swap_research_plan.md` — C = inner/outer
    swap (#63).
  - `docs/topological_discrete_sector_research_plan.md` — `T = iσ_y`,
    `T² = −I` (B2).
  - `geometrodynamics/embedding/transport.py` — `T = iσ_y`.
  - `experiments/closure_ledger/cpt_dirac_operator_probe.py` — this probe.
