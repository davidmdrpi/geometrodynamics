# Hard-wall boundary derivation probe — research plan

Step (2) of the scaffold closure programme. PR #49 (topological/discrete
sector) promoted barriers B1 (closure quantum) and B2 (antipodal Z₂) to
action data, reducing the BAM effective-action scaffold (PR #49 — the
five-barrier map) from 5 barriers to 3. This probe targets the next one:

  - **B3 — hard-wall boundary conditions.** In the scaffold, the
    Dirichlet condition at the throat (`ψ = 0`) was "topologically
    imposed, not derived from `δS/δg`."

The PR #49 result makes B3 tractable: it identified `T = iσ_y`,
`T² = −I` as the **non-trivial spin structure** on `RP³ = S³/Z₂`. The
hard-wall BC should now follow as a *consequence* of that spin
structure (single-valuedness of the spinor at the throat fixed point),
**absorbing B3 into the topological sector** rather than imposing it
separately.

## The derivation chain

```
non-trivial RP³ spin structure (PR #49)
   →  throat transport T = iσ_y,  T² = −I  (eigenvalues ±i)
   →  single-valuedness at the throat fixed point:  ψ = T·ψ
   →  T²·ψ = T·ψ = ψ  but  T²·ψ = −ψ  ⟹  ψ = −ψ  ⟹  ψ(throat) = 0
   →  Dirichlet (hard wall) at the throat
```

The argument: at the throat, the non-orientable identification glues
the two sides via the spinor transport `T`. A field continuous across
the throat must be invariant, `ψ = T·ψ`. Because `T² = −I` (the
non-trivial spin structure, not a free choice), `T` has eigenvalues
`±i` and **no `+1` eigenvector**, so the only invariant spinor is
`ψ = 0`. The Dirichlet condition is forced by the spin structure, not
imposed.

This argument already appears in the prior
`hard_wall_boundary_verification` probe; what is new here is the
**anchoring to the PR #49 topological sector** (the spin structure is
now action data, so the BC is a derived consequence, not an
independent input) and the **concrete realization in the radial
solver**: the Tangherlini mode solver extends each radial mode
antisymmetrically across the throat,

```
u_full(r) = [ −u(reflected) , u ]    (odd across R_MID)
```

which is exactly the `T`-odd transport (spinor flips sign across the
throat) producing the Dirichlet node `u(R_MID) = 0`.

## Predictions

### P1. T eigenvalue structure forbids a +1 invariant

`T = iσ_y` has eigenvalues `±i`; no nonzero spinor satisfies `T·ψ = ψ`.

### P2. Single-valuedness forces ψ = 0 at the throat

`ψ = T·ψ` with `T² = −I` ⟹ `ψ = 0` (the T-fixed-point argument).

### P3. The radial solver realizes the T-odd transport

The Tangherlini mode solver's antisymmetric extension across the throat
(`u_full = [−u_reflected, u]`) is the `T`-odd transport; the throat node
`u(R_MID) = 0` is the Dirichlet condition, with both grid endpoints at
Dirichlet.

### P4. Alternatives inconsistent with the spin structure

Neumann (`u ≠ 0`, `u' = 0` at throat) requires a nonzero `T`-invariant
spinor — impossible under `T² = −I`. Only Dirichlet is consistent.

### P5. The derived BC gives the discrete bulk spectrum

Solving the radial modes with the spin-structure-derived Dirichlet BC
yields the discrete `ω` spectrum that feeds the lepton/quark ladder
(consistency with the closure-ledger thread).

## Tests

  T1. **Spin structure → T eigenvalues** (P1): `T = iσ_y` (PR #49
      non-trivial spin structure); eigenvalues `±i`; no `+1`
      eigenvector.
  T2. **T-fixed-point → Dirichlet** (P2): numerically verify
      `ψ = T·ψ` ⟹ `ψ = 0` under `T² = −I`.
  T3. **Radial solver odd extension** (P3): the mode solver's
      `u_full` is odd across the throat; `u(R_MID) = 0`; both
      endpoints Dirichlet.
  T4. **Alternatives ruled out** (P4): Neumann/Robin require a
      `T`-invariant nonzero spinor — none exists.
  T5. **Discrete spectrum** (P5): solve the radial modes; verify the
      discrete `ω` ladder (hard-wall Bohr–Sommerfeld spacing).
  T6. **Consistency with `hard_wall_boundary_verification`**: the prior
      probe showed DD wins over DN/ND/NN; this probe supplies the
      topological-sector derivation. Cross-check the T²=−I argument.
  T7. **B3 promotion / barrier reduction**: the hard-wall BC is a
      consequence of the PR #49 spin structure, not an independent
      imposition; scaffold barriers 3 → 2.

## Verdict structure

  - **HARD_WALL_DERIVED**: P1–P5 hold; the throat Dirichlet BC follows
    from the PR #49 non-trivial spin structure (`T² = −I`) via
    single-valuedness, realized concretely in the radial solver's
    odd extension. B3 is absorbed into the topological sector; scaffold
    barriers drop 3 → 2.

  - **HARD_WALL_PARTIAL**: the throat BC derives but the outer
    boundary requires a separate (cavity) condition not tied to the
    spin structure.

  - **HARD_WALL_FAILS**: the spin-structure argument does not force
    Dirichlet, or the solver does not realize it.

## What this leaves open (residual 2 barriers)

  - **B4 — dimensional bridge.** The discrete spectrum is computed in
    geometric units (R_MID = 1); the absolute MeV scale needs the
    `m_e` anchor (`ℏ = m_e·R_MID·c`).
  - **B5 — 5D → 4D reduction.** The radial spectrum (this BC) and the
    F² vertex still live in separate sub-threads; the reduction map
    connecting them is unconstructed (the largest gap).

Note: the **outer** boundary at `R_OUTER` is the cavity wall; the
closure-ledger thread fixes `R_OUTER` by the cross-species γ-lock
fixed point. Whether the outer Dirichlet is also a spin-structure
consequence or a distinct cavity condition is noted but not the focus
here (the throat BC is the one B3 names).

## Cross-references

  - PR #49: `topological_discrete_sector_probe` — RP³ non-trivial spin
    structure (`T² = −I`), the anchor for this derivation.
  - PR #49 scaffold: `bam_effective_action_scaffold_probe` — the
    five-barrier map naming B3.
  - `experiments/closure_ledger/hard_wall_boundary_verification.py` —
    prior probe (DD wins; T²=−I argument).
  - `geometrodynamics/tangherlini/radial.py` — the mode solver whose
    odd extension realizes the throat Dirichlet node.
  - `geometrodynamics/embedding/transport.py` — `T = iσ_y`.
  - `experiments/closure_ledger/hard_wall_boundary_derivation_probe.py`
    — this probe.
