# Deriving the cohesive `B·R²` term from the BAM throat action

Follows the self-consistent throat-radius probe (PR #55), which balanced
the EM self-energy (`A/R`, repulsive) against a cohesive term (`B·R²`,
posited as "Poincaré-stress-like") to give the equilibrium throat radius
`R* = (A/2B)^{1/3}`. The `B·R²` was assumed. This probe **derives** it:
identifies it as the throat **brane tension**, derives the `R²` power,
shows it is uniquely selected by power-counting, ties `B` to the bulk
gravity sector, and is honest (per the B4 theorem) that the *value* of
the coupling is still the one dimensionful anchor.

## The derivation

The throat, observed in 4D spacetime, is a 2-surface (the wormhole
mouth). The most general local surface action is a derivative expansion
in the induced metric `h`:

```
S_throat = −σ ∫ √h d²ξ  +  (1/16πG₂) ∫ √h R₂ d²ξ  +  (bending) + …
```

  - **Tension term** `−σ∫√h`: for a 2-sphere mouth of radius `R`,
    `∫√h = Area = 4πR²`, so the static energy is `E_tension = σ·4πR²`
    → **the cohesive `B·R²` with `B = 4πσ`.**
  - **Intrinsic-curvature term** `∫√h R₂`: `R₂ = 2/R²`, `√h R₂·Area =
    (2/R²)(4πR²) = 8π` — Gauss–Bonnet, `R`-independent (topological), no
    `R`-dependence.
  - Higher (bending `∝R⁰`, …) are subleading.

So the **leading cohesive `R`-dependence is exactly `σ·4πR²`** — the
brane tension. `B = 4πσ`.

## Why `R²` is unique (power-counting)

Each term of `S_BAM` evaluated on a throat of radius `R` has a definite
`R`-power:

| term in `S_BAM` | physical origin | `R`-scaling |
|---|---|---|
| `¼ F²` (Coulomb) | EM self-energy outside the capped throat | `1/R` (repulsion) |
| `ψ̄(iΓD−m)ψ` | radial zero-point `ℏω ∝ 1/R` | `1/R` |
| Israel junction `[K]` | induced surface stress of the glued throat | `R¹` |
| `R₅/2κ₅` (Einstein–Hilbert) | bulk curvature on the throat | `R¹` |
| **`L_throat` brane tension** | **constant surface tension `σ·Area`** | **`R²`** |
| `Λ₅` (cosmological) | vacuum "bag" energy in the throat volume | `R³` |

The cohesive partner to the leading `1/R` repulsion that yields a
*stable* minimum at the smallest non-trivial power is the `R²` brane
tension. `R²` is the unique signature of a **constant** surface tension —
distinct from the induced junction tension (`R¹`), curvature/EH (`R¹`),
and the cosmological bag (`R³`).

## The induced junction tension scales as `R¹`, not `R²`

A check that the `R²` term must be a *fundamental* brane tension (a
constant `σ` in `L_throat`), not the *induced* Israel/Lanczos junction
tension of the Tangherlini geometry. For a symmetric thin shell at radius
`a` in `f(r) = 1 − (rs/r)²`, the Israel surface energy density is
`σ_Israel(a) = −√(f(a))/(2πa)`, so the surface energy is

```
E_Israel(a) = σ_Israel · 4πa² = −2a √(1 − (rs/a)²)  ∝  a¹   (a ≫ rs).
```

Computing the local log-log slope confirms `d ln|E_Israel|/d ln a → 1`,
**not 2**. So the induced junction tension is `R¹`; the cohesive `R²`
term is a fundamental brane tension in `L_throat`.

## Where the coupling lives (and the B4 caveat)

`B = 4πσ` carries dimension `[σ] = energy/area`. In a brane-world /
thin-shell embedding the throat tension is set by the bulk gravity
sector — Randall–Sundrum-like, `σ ∝ √(|Λ₅|)/κ₅` (a dimensionless factor
times the bulk cosmological constant and gravitational coupling). So the
`R²` form and the identity (`B = 4πσ`, brane tension) are **derived**;
the *value* of `σ` (equivalently `Λ₅, κ₅`) is the single dimensionful
anchor — exactly the one external input the B4 scale-modulus theorem
(PR #52) says is mandatory. Rescaling `σ → σ/λ³` sends
`R* = (A/8πσ)^{1/3} → λ R*`.

## Tests

  T1. **Brane-tension area term.** Constant surface tension on the
      2-sphere mouth gives `E = σ·4πR²` → `B = 4πσ` (the `R²` form).
  T2. **Power-counting uniqueness.** The `S_BAM` term `R`-powers
      (`1/R, 1/R, R, R, R², R³`); `R²` ⟺ constant surface tension only.
  T3. **Induced junction tension is `R¹`.** `E_Israel(a) = −2a√(1−(rs/a)²)`
      from `f(r) = 1−(rs/r)²` has log-log slope → 1, not 2 → the `R²`
      term is fundamental, not induced.
  T4. **Dimensional consistency.** `[B] = [4πσ] = energy/area`;
      `R* = (A/8πσ)^{1/3}` has dimension length.
  T5. **B4 scale-modulus.** `σ → σ/λ³` ⟹ `R* → λ R*`: `σ` is the single
      dimensionful coupling.
  T6. **Reproduces PR #55.** `E(R) = A/R + 4πσ R²` has the stable
      minimum `R* = (A/8πσ)^{1/3}`; matches the PR #55 equilibrium.
  T7. **Coupling from the bulk gravity sector.** `σ ∝ √(|Λ₅|)/κ₅`
      (RS-like): the anchor is the bulk cosmological/gravitational scale.
  T8. **Assessment.** The `R²` term is derived as the throat brane
      tension (form + power-counting + junction discriminator); `B=4πσ`;
      value still the one anchor (B4).

## Verdict structure

  - **COHESIVE_TENSION_DERIVED** (expected): the `B·R²` term is the
    throat brane tension `σ·Area = 4πσR²` — its `R²` power derived from a
    constant surface tension and uniquely selected by power-counting
    (the induced Tangherlini junction tension is `R¹`, computed). `B =
    4πσ` with `σ` set by the bulk gravity sector. The form and identity
    are derived; the value of `σ` is the single dimensionful anchor
    (B4-consistent).

  - **DERIVATION_FAILS**: the `R²` power is not reproduced, or the
    junction tension is not distinguishable from a fundamental tension.

## What this leaves open

  - **The value of `σ` (the anchor).** Still one dimensionful coupling
    (equivalently `Λ₅, κ₅`); deriving it needs a second fixed scale.
  - **The RS-like tuning from the BAM action.** `σ ∝ √(|Λ₅|)/κ₅` is the
    brane-world form; deriving the exact dimensionless factor from
    `S_BAM`'s junction conditions is the follow-on.
  - **Pair-production threshold.** `2 m_e c²` at the lowest stable `R*`
    as a dynamical nucleation calculation.

## Cross-references

  - `docs/self_consistent_throat_radius_research_plan.md` — the `A/R +
    B·R²` equilibrium (#55).
  - `docs/maslov_dimensional_bridge_research_plan.md` — the B4 theorem.
  - `docs/bam_scaffold_status.md` — barrier ledger; `S_BAM` candidate.
  - `geometrodynamics/tangherlini/radial.py` — `f(r) = 1 − (rs/r)²`.
  - `experiments/closure_ledger/cohesive_tension_derivation_probe.py` —
    this probe.
