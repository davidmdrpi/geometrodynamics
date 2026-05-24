# Charge conjugation from the inner/outer swap

The standing THESIS open problem: *"Charge conjugation from inner/outer
swap. Promote the C-symmetry from a postulate to a geometric statement —
that swapping `r < R_MID ↔ r > R_MID` in the throat eigenmodes flips the
sign of the integrated Hopf curvature."* This probe makes that statement
precise and checks it: the inner/outer reflection across the throat is an
involution under which the throat eigenmodes are odd and the integrated
Hopf curvature (the charge `c₁`) flips sign — so **C = the inner/outer
swap**, geometric rather than postulated.

## The swap

The two regions of the wormhole are `r < R_MID` (inner) and `r > R_MID`
(outer), with `R_INNER = R_MID − ΔR` and `R_OUTER = R_MID + ΔR`
symmetric about the throat. The inner/outer swap is the **reflection
across the throat**

```
S : r ↦ 2 R_MID − r ,
```

which fixes `R_MID`, exchanges `R_INNER ↔ R_OUTER`, and is an involution
(`S² = id`). It is exactly the reflection of the B3 hard-wall odd
extension: the radial throat modes are continued antisymmetrically across
the throat (`u(2R_MID − r) = −u(r)`), so the eigenmodes are **odd** under
`S`.

## The charge is the integrated Hopf curvature

The Hopf connection `A_φ = ½ cos χ` has curvature `F = dA = −½ sin χ
dχ∧dφ`, and the charge is the integrated curvature (the first Chern
number),

```
c₁ = (1/2π) ∮ F = ±1
```

(`geometrodynamics/hopf/chern.py`), with the sign set by the orientation
of the mouth 2-surface.

## Why the swap flips the charge

The wormhole mouth's induced orientation is set by its **outward
normal** `n̂ = ±r̂`: the outer mouth has `n̂ = +r̂`, the inner mouth
`n̂ = −r̂`. The two normals point oppositely, so the inner and outer mouths
carry **opposite** induced orientation, and

```
c₁(inner) = −c₁(outer) .
```

The swap `S` exchanges the mouths (and reverses `dr`), reversing the
mouth orientation and flipping the integrated curvature: `c₁ → −c₁`. The
two orientations are exactly the `c1_chiphi = −1` and `c1_phichi = +1` of
`compute_c1`. Combined with the modes being odd under `S`, the swap takes
a throat (`c₁ = +1`) to its **antithroat** (`c₁ = −1`) — the
charge-conjugate partner used in pair production (PR #58) and the
antipodal `Z₂` (B2).

## So C is geometric

Charge conjugation — particle → antiparticle — is realized as the
inner/outer swap: `C = S`, with `C: c₁ → −c₁`, `C² = id`. It is no longer
a postulate but the throat-reflection involution, consistent with the
antipodal `Z₂` / `T = iσ_y` (B2) and the pair-production antithroat
(#58).

## B4 accounting

`c₁ = ±1` is a **dimensionless topological integer** (the first Chern
number); `C` is a discrete geometric involution. Neither carries a scale
— independent of the single anchor `m_e` (B4-consistent).

## Tests

  T1. **Charge = integrated Hopf curvature.** `c₁ = (1/2π)∮F = ±1`
      (`compute_c1`).
  T2. **The inner/outer swap is an involution.** `S: r ↦ 2R_MID − r`
      fixes `R_MID`, exchanges `R_INNER ↔ R_OUTER`, `S² = id`; the throat
      modes are odd under `S` (the B3 antisymmetric extension).
  T3. **The swap reverses the mouth orientation.** The outward normal
      `n̂ = ±r̂` is opposite for inner/outer; the swap reverses it.
  T4. **The integrated curvature flips.** `c₁ → −c₁` (the two
      orientations `c1_chiphi = −1`, `c1_phichi = +1`).
  T5. **C = swap (the geometric statement).** `C: throat → antithroat`,
      `c₁ → −c₁`; tie to the pair-production antithroat (#58) and the
      antipodal `Z₂` (B2).
  T6. **C² = id; discrete-symmetry consistency.** The swap is an
      involution; consistent with `T = iσ_y` / the double cover (CPT).
  T7. **Falsification / B4.** A swap-invariant charge would break C;
      BAM flips `c₁`. `c₁` is a topological integer, scale-independent.
  T8. **Assessment.**

## Verdict structure

  - **C_IS_INNER_OUTER_SWAP** (expected): charge conjugation is the
    inner/outer reflection `S: r ↦ 2R_MID − r` — an involution fixing the
    throat, under which the eigenmodes are odd (B3) and the integrated
    Hopf curvature flips sign (`c₁ → −c₁`, the two mouth orientations
    `∓1`), taking a throat to its antithroat. C is promoted from a
    postulate to a geometric statement, consistent with the antipodal
    `Z₂` (B2) and the pair-production antithroat (#58).

  - **C_NOT_GEOMETRIC**: the swap does not flip `c₁`, or is not an
    involution — C would remain a postulate.

## What this leaves open

  - **The full CPT theorem from `S_BAM`.** C (this probe), P (parity),
    and T (`T = iσ_y`, B2) as one geometric CPT statement on the throat.
  - **C on the Dirac spinor.** The explicit action of `S` on the throat
    spinor (`ψ → C ψ̄ᵀ`), beyond the charge sign.

## Cross-references

  - `docs/pair_production_threshold_research_plan.md` — the antithroat
    (`c₁ = −1`) in pair production (#58).
  - `docs/topological_discrete_sector_research_plan.md` — the antipodal
    `Z₂` / `T = iσ_y` (B2).
  - `geometrodynamics/hopf/chern.py` — `compute_c1` (the Chern number).
  - `geometrodynamics/constants.py` — `R_MID, R_INNER, R_OUTER`.
  - `experiments/closure_ledger/charge_conjugation_swap_probe.py` — this
    probe.
