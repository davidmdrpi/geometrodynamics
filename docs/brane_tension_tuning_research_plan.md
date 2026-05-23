# Brane-tension / bulk-gravity tuning probe

Follows the cohesive-tension derivation (PR #56), which identified the
throat cohesive term as the brane tension `B·R² = σ·4πR²` and tied it
parametrically to the bulk gravity sector, `σ ∝ √|Λ₅|/κ₅`. This probe
**derives the exact Randall–Sundrum-like fine-tuning** from the Israel
junction conditions: the dimensionless tuning factor is **√6**, the
tuning is the **flat / static-throat condition** (zero induced 4D
cosmological constant), and the absolute scale remains the single
dimensionful anchor (B4).

## The derivation

A `Z₂`-symmetric pure-tension throat (3-brane worldvolume) in an `AdS₅`
bulk. Three ingredients:

  1. **Israel junction (pure tension).** For surface stress-energy
     `S_μν = −λ h_μν`, the `Z₂` junction `2(K_μν − K h_μν) = −κ₅² S_μν`
     gives, after trace-reversal in a 4D worldvolume,

     ```
     K_μν = −(κ₅² λ / 6) h_μν .
     ```

  2. **Bulk `AdS₅`.** The warped vacuum Einstein equation fixes the AdS
     inverse radius `k` (warp `a(y) = e^{−k|y|}`):

     ```
     Λ₅ = −6 k²        ⟹   k = √(|Λ₅| / 6) .
     ```

  3. **Staticity (flat brane).** A flat, static brane has extrinsic
     curvature `K_μν = k h_μν`. Matching to (1):

     ```
     k = κ₅² λ / 6     ⟹   λ_crit = 6k / κ₅² = √(6 |Λ₅|) / κ₅² .
     ```

The **dimensionless tuning factor** is therefore

```
λ_crit · κ₅² / √|Λ₅|  =  √6  ≈ 2.449 ,
```

sharpening PR #56's parametric `∝` to the exact `√6` (the `6` is the
`AdS₅` curvature coefficient `Λ₅ = −6k²`).

## The tuning is the flat / static-throat condition

Detuning the tension induces a 4D cosmological constant on the throat,

```
Λ₄  ∝  λ² − λ_crit²  ,
```

so `Λ₄ = 0 ⟺ λ = λ_crit`: the fine-tuning is exactly the condition for a
**flat, static throat** (a stable particle with no induced 4D vacuum
energy). Over-tension (`λ > λ_crit`) gives `Λ₄ > 0` (de Sitter throat),
under-tension `Λ₄ < 0` (anti-de Sitter). The cohesive equilibrium of
PR #55 is the critically-tuned, static configuration.

## Connection to the cohesive term and B4

The throat cohesive tension `σ` (PR #56, `B = 4πσ`) is the 4D-effective
image of this RS-tuned brane tension `λ_crit`; it inherits the `√6`
factor and the bulk-gravity scale `√|Λ₅|/κ₅²`. The fine-tuning is **one
condition** relating the three couplings `(λ, Λ₅, κ₅)`, so a net **one**
dimensionful combination remains — the single anchor the B4 scale-modulus
theorem (PR #52) requires. The dimensionless content (`√6`, the flatness
condition) is derived; the absolute scale (`k = √|Λ₅/6|`, the bulk
gravitational scale) is the one external input. Rescaling `k → k/λ_s`
rescales the throat radius.

## Hierarchy (warp factor)

The `AdS₅` warp `a(y) = e^{−k|y|}` over the bulk depth `L` (the BAM bulk
separation `ΔR`, #53) gives an exponential hierarchy `w = e^{−kL}`
between the bulk and throat scales — the RS hierarchy mechanism, here
relating the bulk-gravity scale to the throat geometry.

## Tests

  T1. **Israel junction (pure tension).** Trace-reversal in a 4D
      worldvolume gives `K_μν = −(κ₅² λ/6) h_μν` (the `1/6` factor).
  T2. **Bulk `AdS₅`.** `Λ₅ = −6k²` ⟹ `k = √(|Λ₅|/6)`.
  T3. **Fine-tuning factor √6.** `λ_crit = 6k/κ₅² = √(6|Λ₅|)/κ₅²`;
      `λ_crit κ₅²/√|Λ₅| = √6 ≈ 2.449`.
  T4. **Flat / static-throat condition.** `Λ₄ ∝ λ² − λ_crit²`; zero at
      `λ_crit`, `Λ₄>0` over-tension, `Λ₄<0` under-tension.
  T5. **Connection to the cohesive term.** `B = 4πσ` inherits `√6` and
      the bulk scale; the tuning = the static-throat condition (PR #55).
  T6. **B4 accounting.** One tuning condition among `(λ, Λ₅, κ₅)` → one
      dimensionful combination remains; `√6` derived, scale is the anchor.
  T7. **Warp hierarchy.** `w = e^{−kL}` over the bulk depth `L = ΔR`.
  T8. **Assessment.** RS-like tuning derived (factor `√6`, flatness
      condition); cohesive tension fixed up to the bulk scale; value the
      one anchor (B4).

## Verdict structure

  - **BRANE_TUNING_DERIVED** (expected): the RS-like fine-tuning of the
    throat brane tension is derived from the Israel junction + bulk
    `AdS₅` equations: `λ_crit = √(6|Λ₅|)/κ₅²`, the dimensionless factor
    `√6`, and the tuning is the flat/static-throat condition (zero
    induced `Λ₄`). The cohesive `B = 4πσ` inherits the factor and the
    bulk-gravity scale; the absolute value remains the single
    dimensionful anchor (B4-consistent).

  - **TUNING_FAILS**: the junction algebra does not give the `1/6`/`√6`
    factors, or the flatness condition does not vanish at `λ_crit`.

## What this leaves open

  - **The bulk-gravity scale (the anchor).** `k = √|Λ₅/6|` (equivalently
    `Λ₅, κ₅`) is the one dimensionful input; deriving it needs a second
    fixed scale.
  - **The BAM-throat junction from `S_BAM`.** The derivation uses the
    canonical RS `Z₂` brane; matching to the exact BAM throat
    (Tangherlini interior + closure-quantum surface) is the follow-on.
  - **Pair-production threshold.** `2 m_e c²` at the lowest stable `R*`
    as a dynamical nucleation calculation.

## Cross-references

  - `docs/cohesive_tension_derivation_research_plan.md` — `B = 4πσ`,
    `σ ∝ √|Λ₅|/κ₅` (#56).
  - `docs/self_consistent_throat_radius_research_plan.md` — the
    equilibrium (#55).
  - `docs/maslov_dimensional_bridge_research_plan.md` — the B4 theorem.
  - `docs/delta_r_scale_modulus_research_plan.md` — `ΔR` (the bulk depth).
  - `experiments/closure_ledger/brane_tension_tuning_probe.py` — this
    probe.
