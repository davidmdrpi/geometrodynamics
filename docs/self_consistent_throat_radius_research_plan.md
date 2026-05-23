# Self-consistent throat radius / finite self-energy probe

Targets the **remaining BAM scale anchor**: the throat radius `R_MID`
(≈ the invariant bulk separation `ΔR`), which is currently *imposed*
(`constants.py`). The THESIS open problem "Self-consistent throat radius"
asks for `R_MID` to be determined **dynamically as the equilibrium
throat radius**, with finite self-energy and the pair-production
threshold as the lowest stable configuration.

This probe builds that equilibrium and reports — honestly, against the
B4 scale-modulus theorem (PR #52) — exactly what it does and does not
fix.

## The B4 constraint (what is and isn't possible)

The B4 audit proved the closure-ledger machinery is **scale-free**: it
can fix dimensionless ratios but not an absolute length. So a
self-energy equilibrium built from scale-free ingredients (every term
`∝ 1/R`) has **no** stationary point in `R` — it can only fix a
dimensionless balance. To pin an absolute `R*`, the energy functional
must contain **two terms with different `R`-scaling**, at least one
carrying a dimensionful coupling. The probe makes this explicit rather
than hiding it.

## Finite self-energy — the BAM regularization

A point charge has divergent Coulomb self-energy
`U = e²/(8πε₀ r) → ∞` as `r → 0`. The BAM throat removes the divergence
**geometrically** from both ends:

  - **Short distance:** the wormhole throat is the inner boundary
    (`f(r) = 1 − (rs/r)² = 0` at `r = rs = R_MID`); there is no `r <
    R_MID`. The field is capped at the throat, so the self-energy is
    finite, `U_EM = α ℏc/(2 R_MID)` — the cutoff is the throat radius
    itself, not an arbitrary `ε`.

  - **Long distance:** on compact `S³` the Green function
    `G(ψ) = ((π−ψ)cot ψ − ½)/(4π²R)` is regular at the antipode
    (`G(π)` finite) with a zero-mean neutralizing background (the `−½`).
    No far-field divergence.

So the throat self-energy is finite — `U_EM/(m c²) = α/2 ≈ 0.0036`, a
small finite correction with **no UV divergence** (no renormalization
needed), in contrast to QED's divergent electron self-energy.

## The self-consistent equilibrium

Write the throat energy as a competition between the EM self-energy
(repulsive, `∝ 1/R`, wants to expand) and a cohesive term (`∝ R²`, wants
to shrink):

```
E(R) = A/R + B·R² ,    A = g·α ℏc/2 ,   B = cohesive stiffness
```

Stationarity gives a **unique stable equilibrium**:

```
dE/dR = −A/R² + 2B R = 0   ⟹   R* = (A / 2B)^{1/3} ,   d²E/dR² > 0.
```

The throat radius is then a self-consistent stationary point, not an
imposed constant. **But** `R*` depends on `A/B`: rescaling `B → B/λ³`
sends `R* → λ R*` — the one scale modulus. The absolute value rides on
the single dimensionful coupling `B` (or `A/B`), exactly as B4 requires.

## The pure-EM relocation

If instead of an ad-hoc cohesive term one demands the BAM-native balance
"rest energy = EM self-energy" (`m c² = U_EM`), both sides scale as
`1/R`, so the balance is **R-independent** — a dimensionless identity
`g = 2/α`, relating the geometric self-energy factor to the
fine-structure constant rather than fixing `R`. This is the B4
obstruction in self-energy language: the self-consistency relocates the
scale question to `α` (still one input), it does not fix the length.

## What this advances

The throat radius is realized as a **finite-self-energy stable
equilibrium**, replacing "imposed `R_MID`" with "`R_MID` = stationary
point of the throat self-energy." The self-energy is finite (the throat
caps the field — no UV divergence). The relocation chain sharpens:

```
imposed R_MID  →  ΔR invariant geometric length (#53)
               →  finite-self-energy equilibrium R* = (A/2B)^{1/3} (this probe)
```

Each step makes the anchor more physical; none derives its absolute
value — by the theorem, that needs one external dimensionful coupling
(here the cohesive stiffness `B`, equivalently a relation to `α`).

## Tests

  T1. **Finite self-energy (short distance).** Point charge diverges as
      cutoff→0; the throat caps the field at `R_MID` → `U_EM` finite.
  T2. **S³ far-field regular.** `G(ψ)` finite at the antipode; zero-mean
      `−½` background — no long-distance divergence.
  T3. **Scale-free obstruction.** EM-only `E(R)=A/R` is monotone → no
      equilibrium from `1/R` alone (B4 in self-energy form).
  T4. **Self-consistent equilibrium.** `E(R)=A/R+B R²` has a unique
      stable minimum `R*=(A/2B)^{1/3}`; verify numerically.
  T5. **Scale modulus.** `B → B/λ³` ⟹ `R* → λ R*`: absolute `R*` rides
      on one dimensionful coupling (B4-consistent).
  T6. **Finite, renormalization-free.** `U_EM/(m c²) = α/2 ≈ 0.0036` —
      a finite correction, no UV divergence (contrast QED).
  T7. **Pure-EM relocation.** `m c² = U_EM` is R-independent (both
      `∝1/R`) ⟹ fixes `g = 2/α`, not `R` — relocation to `α`.
  T8. **Assessment.** Throat radius = finite-self-energy stable
      equilibrium (not imposed); absolute value still one anchor (B4);
      progress = explicit equilibrium + finite self-energy + α-relocation.

## Verdict structure

  - **SELF_CONSISTENT_THROAT_EQUILIBRIUM** (expected): the throat radius
    is realized as a finite-self-energy stable stationary point, and the
    self-energy is finite (throat caps the field; no UV divergence). The
    absolute value rides on one dimensionful coupling, consistent with
    the B4 theorem; the self-consistency recasts the anchor as an
    equilibrium condition (and relates it to `α`), without deriving the
    absolute value.

  - **EQUILIBRIUM_FAILS**: no stable stationary point exists, or the
    self-energy is not finite.

## What this leaves open

  - **The absolute value.** Still one dimensionful coupling (`B`, or a
    relation to `α`). A genuine derivation would pin it to a second
    fixed scale (e.g. a closure-quantum relation to the Planck length).
  - **The cohesive term from first principles.** `B·R²` is a Poincaré-
    stress-like cohesion; deriving it from the BAM throat action (rather
    than positing the scaling) is the follow-on.
  - **Pair-production threshold.** Identified structurally (`2 m_e c²` =
    nucleating a throat–antithroat pair at the lowest stable `R*`); a
    dynamical nucleation calculation is future work.

## Cross-references

  - `docs/maslov_dimensional_bridge_research_plan.md` — the B4 audit /
    scale-modulus theorem (#52).
  - `docs/delta_r_scale_modulus_research_plan.md` — ΔR invariant (#53).
  - `docs/scaffold_closure_release_note.md` — scaffold closure summary.
  - `geometrodynamics/transaction/s3_geometry.py` — `s3_green_potential`.
  - `geometrodynamics/hopf/connection.py` — Hopf curvature `|F|=½sin χ`.
  - `geometrodynamics/tangherlini/radial.py` — `f(r)=1−(rs/r)²` throat.
  - `experiments/closure_ledger/self_consistent_throat_radius_probe.py` —
    this probe.
