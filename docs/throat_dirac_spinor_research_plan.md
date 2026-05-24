# Throat Dirac spinor from S_BAM

Closes the open note of the explicit-CPT-operator probe (PR #65), which
*posited* the throat Dirac 4-spinor as the inner/outer mouth doubling.
This probe **derives** that doubling from the BAM radial structure via the
**Dirac/SUSY factorization** — the throat spinor is the square root of the
repo's second-order radial operator — together with the wormhole's
two-sidedness. It also disentangles parity `P` from the antipodal `Z₂`,
the other #65 open note.

## Honest scope

This derives the doubling **structure** (why the throat field is a
4-component Dirac spinor, what the components are), computably, from the
established radial operator. It does **not** write the full closed-form
bulk Dirac spinor of the complete `S_BAM` (all sectors, the S³ angular
coupling included) — that remains the radial/throat sector only, stated
plainly.

## The Dirac factorization

The closure-ledger radial operator is the second-order Schrödinger-like

```
H = −d²/dr*² + V_tangherlini(r, l) ,   spectrum  ω²(l, n) .
```

A Dirac field is first-order; squaring the Dirac operator gives a
second-order operator. Conversely, `H` factorizes as a perfect square
(SUSY QM / the radial Dirac factorization):

```
H − E₀ = A†A ,   A = d/dr* + W(r*) ,   W = −ψ₀'/ψ₀ ,
```

with `E₀ = ω₀²` the lowest eigenvalue and `ψ₀` its mode. Equivalently the
Riccati identity `V − E₀ = W² − W'`. **`A` is the first-order radial Dirac
operator**, and `H = A†A + E₀` is its square. So the Dirac structure is
literally the square root of the BAM radial operator — a 2-component
object (`A` maps between two sectors) where `H` was 1-component. That is
the doubling.

## The two partner sectors (the two mouths)

`A†A` (`= H − E₀`) and its partner `AA†` share their **nonzero**
spectrum (SUSY isospectrality) — the two Dirac components carry the same
mass ladder `ω(l, n)`. These are the two wormhole **mouths**: the field
has support on both sides of the throat (`r < R_MID` inner,
`r > R_MID` outer), joined at the throat by the B3 hard wall / the odd
extension `u(2R_MID − r) = −u(r)` (the PR #63 reflection).

## The 4-spinor (derived, not posited)

Combining the two structures:

```
throat Dirac 4-spinor = [2 : Dirac square-root (A/A†, the two mouths)]
                      ⊗ [2 : SU(2) spin (T = iσ_y, B2)]
                      = Ψ_inner ⊕ Ψ_outer .
```

The 4-component throat spinor of PR #65 is **derived**: the factor-2 from
the Dirac square root of the radial operator, the factor-2 from the B2
spin structure. The inner/outer mouths are the two SUSY-partner sectors.

## P vs the antipodal Z₂ (the other #65 note)

These act on **different factors** and are now disentangled:

  - **Parity `P = γ⁰`** acts on the **radial/Dirac** components — the
    inner/outer reflection `r → 2R_MID − r` across the throat (the #63
    radial swap direction), exchanging the mouth/chirality blocks.
  - **The antipodal `Z₂` (B2)** acts on the **S³ angular base** — the
    deck transformation of `RP³ = S³/Z₂` (`σ: p → −p`), independent of
    the radial sector.

`P` is radial (Dirac), the antipodal `Z₂` is angular (the base) — distinct
operations on distinct tensor factors.

## B4 accounting

The spinor structure (component counting, the factorization, the
isospectrality) is **dimensionless and geometric**; the mass scale
`ω(l, n)` rides on the single anchor `m_e` (`ω·R_MID` is dimensionless).
Independent of the anchor's value.

## Tests

  T1. **Dirac factorization.** `H − E₀ = A†A`, `V − E₀ = W² − W'`
      (`W = −ψ₀'/ψ₀`) verified in the interior — `H` is the square of the
      first-order `A`.
  T2. **A reproduces the BAM spectrum.** the discrete `A†A + E₀` low
      spectrum = `ω²(l, n)` (the closure-ledger ladder).
  T3. **SUSY isospectrality (the two mouths).** nonzero spec(`A†A`) =
      nonzero spec(`AA†`) — the two partner sectors share `ω(l, n)`.
  T4. **Inner/outer mouths + B3 odd extension.** the wormhole
      two-sidedness; `u(2R_MID − r) = −u(r)` joins them at the throat
      (#63).
  T5. **The 4-spinor derived.** `4 = 2 (Dirac square-root / mouths) × 2
      (SU(2) spin, B2)` = `Ψ_inner ⊕ Ψ_outer` — PR #65's spinor derived.
  T6. **P vs antipodal Z₂.** `P = γ⁰` radial (Dirac); antipodal `Z₂`
      angular (the S³ base) — distinct factors.
  T7. **Honest scope / B4.** structure derived; full closed-form bulk
      spinor (S³ coupling) open; dimensionless/geometric.
  T8. **Assessment.**

## Verdict structure

  - **THROAT_DIRAC_DERIVED** (expected): the throat Dirac 4-spinor
    structure is derived — the radial operator `H = −d²/dr*² + V` is a
    perfect square `A†A + E₀` (`A` the first-order radial Dirac operator,
    `V − E₀ = W² − W'`), its two SUSY-partner sectors (`A†A`, `AA†`)
    isospectral on the nonzero spectrum are the two wormhole mouths
    (joined by the B3 odd extension, #63), and `4 = 2 (mouths) × 2 (SU(2)
    spin, B2)` gives `Ψ_inner ⊕ Ψ_outer` — PR #65's posited spinor
    derived. Parity (`γ⁰`, radial) and the antipodal `Z₂` (angular) are
    disentangled. The full closed-form bulk spinor (with the S³ coupling)
    remains open.

  - **DIRAC_STRUCTURE_INCOMPLETE**: the factorization fails, the
    partners are not isospectral, or the doubling does not assemble.

## What this leaves open

  - **The full bulk Dirac spinor from `S_BAM`.** The radial/throat sector
    is derived here; the S³ angular spinor harmonics and the complete
    coupled solution are the next step.
  - **The superpotential `W` in closed form.** `W = −ψ₀'/ψ₀` is computed
    numerically; a closed form (from the Tangherlini `f(r)`) would sharpen
    it.

## Cross-references

  - `docs/cpt_dirac_operator_research_plan.md` — the posited 4-spinor
    (#65).
  - `docs/charge_conjugation_swap_research_plan.md` — the inner/outer
    swap / B3 odd extension (#63).
  - `docs/topological_discrete_sector_research_plan.md` — `T = iσ_y`,
    the antipodal `Z₂` (B2).
  - `geometrodynamics/tangherlini/radial.py` — `V_tangherlini`, the
    radial operator.
  - `experiments/closure_ledger/throat_dirac_spinor_probe.py` — this
    probe.
