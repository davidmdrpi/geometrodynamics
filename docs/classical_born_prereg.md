# Pre-registration: derive or refute the Born rule from classical BAM measure

**Status: frozen before any code.** Fourth round of the finite-mouth chain,
after `docs/mouth_spin_frame_prereg.md` (`6bc4306`). The spin-transport arc
is frozen at `7bd7ecc`; this round does not touch it, nor masses, QED
vertices, composition or CHSH. Success criteria may change afterwards only
by an explicit correction note.

## 0. The question, and what is not the question

A classical wave has intensity proportional to amplitude squared. That does
not explain why a single detection is stochastic with Born frequencies. The
question is:

> Can BAM derive why **individual outcomes** are sampled with probability
> `|A|²` from its classical global dynamics and a classically derived
> ensemble measure?

The microscopic classical state is what the previous rounds supply: a point
of the mouth spin-frame bundle `q ∈ S³ = P_Spin(S²)`, i.e. a Bloch direction
`x = q⁻¹iq ∈ S²` and a fibre phase `φ ∈ S¹` (the tangent frame rotates by
`2φ`), together with the Pin⁻ sector `ε` and the classical field. A
preparation `Ψ` fixes `x = x_Ψ` and leaves `φ` unresolved. A detector with
analyzer axis `a ∈ S²` is a deterministic map `D_a(X) ∈ {+1, −1}`. No
projector, no `|amplitude|²`, no random number drawn with quantum weights,
no tensor-product state.

## 1. Three theorems, fixed before coding

* **Theorem 1 (intensity is not probability).** Obtaining a `cos²`
  decomposition of classical field energy is *not* evidence for the outcome
  law. Any result whose only Born-shaped quantity is an intensity fraction
  is classified `CLASSICAL_INTENSITY_ONLY_NO_OUTCOME_PROBABILITY`.
* **Theorem 2 (invariant measure).** The ensemble measure must be
  flow-invariant or the discrete/global analogue: Haar on the fibre `S¹`,
  Haar on `SU(2) = S³`, Liouville, microcanonical. No distribution may be
  tuned to the Born answer. A basin shape tuned to the answer counts as a
  tuned distribution.
* **Theorem 3 (deterministic detector).** `D_a(X)` is a function of the
  complete classical state and boundary conditions. Any stochasticity must
  come from unresolved geometry or measure.

## 2. Classification, done before coding (H1)

By rotational covariance a deterministic detector depends only on the
relative configuration of `a` and the frame `(x, e₁, e₂)`: the polar angle
`θ = ∠(a, x)` and the azimuth `ψ` of `a` in the tangent plane at `x`.
Fibre Haar makes `ψ` uniform (the `2φ` double cover is invisible to a uniform
measure). Hence

```
P(+|θ) = (arc measure of {ψ : D(θ, ψ) = +1}) / 2π  ≡  f(θ)
```

* **H1a.** Analyzer reversal `D_{−a} = −D_a` (`θ ↦ π−θ`, `ψ ↦ ψ+π`) forces
  `f(π−θ) = 1 − f(θ)` and nothing else; complementarity of orthogonal
  outcomes (`x ↦ −x` with the lifted frame `(e₁, −e₂)`) adds
  `D(θ, ψ) = D(θ, π−ψ)`. Born `cos²(θ/2)`, the straight line `1 − θ/π`, and
  the step `Θ(π/2 − θ)` all satisfy these. **Symmetry plus Haar does not
  select Born.** *Falsifier:* a symmetry in the list that excludes the line
  or the step.
* **H1b.** For every `f` with `f(π−θ) = 1−f(θ)` there is a symmetric basin
  map realising it (`D = +1 iff |ψ| < π f(θ)`). Born is the basin
  `|ψ| < π cos²(θ/2)`: realisable, hence not excluded, but a *specific basin
  shape* that must come from a dynamical coupling to count as derived.

## 3. The dynamical candidates, with pre-computed numbers (H2)

Each candidate is a detector built from an actual classical coupling, with
its induced `f(θ)` under fibre Haar computed as an exact arc measure. Born
means `max_θ |f(θ) − cos²(θ/2)| < 1e-3`.

* **C1 — linear functional of the frame** (a torque or gradient force
  coupling `a` to a fixed combination of `x, e₁, e₂`):
  `D = sign(α cos θ + ρ sin θ cos(ψ − ψ₀))`. Induced
  `f = 1` for `cot θ > ρ/α`... explicitly
  `f(θ) = (1/π) arccos(−(α/ρ) cot θ)` clipped to `[0, 1]`. This family has
  plateaus at `0` and `1` (or is identically `1/2` when `α = 0`).
  *Pre-computed:* the best member (`α/ρ ≈ 0.80`) misses Born by
  `max|f − Born| = 0.109`. **Not Born.**
* **C2 — quadratic / intensity functional** (the classical Malus
  intensities `|⟨a,±|q⟩|² = cos²(θ/2), sin²(θ/2)`): these are
  `φ`-independent, so `D = sign(I₊ − I₋) = sign(cos θ)`: the step. The
  intensity fraction *is* `cos²(θ/2)` — Theorem 1 applies. **Intensity,
  not probability.**
* **C3 — two-harmonic coupling** (the detector sees the spinor at `e^{iφ}`
  and the frame at `e^{2iφ}`): `D = sign(A cos(φ+δ₁) + B cos(2φ+δ₂))` with
  the natural geometric weightings `(A, B) ∈ {(cos θ/2, sin θ/2), (sin θ,
  cos θ), …}`. *Pre-computed:* `max|f − Born| ≥ 0.50` for all four natural
  weightings. **Not Born.**
* **C4 — the repository's own detector** (Probe K's winding-number
  Stern–Gerlach with the fibre-`U(1)` setting): winding-diagonal, so
  `φ`-independent: the step, and the abelian `span{I, σ_z}` of PR #238.
  **Not Born for a superposed preparation.**
* **C5 — the Archimedes route.** Let the detector carry an unresolved
  direction `y ∈ S²` with Haar measure and register the sign of the total
  polarisation along `a`: `D = sign(a·(x + κ y))`. Since `a·y` is uniform on
  `[−1, 1]` for Haar `y` (Archimedes' hat-box theorem),
  `f(θ) = clip((1 + cos θ / κ)/2)`. *Pre-computed:* `κ = 1` gives Born
  **exactly** (`max|f − Born| = 0`; Monte Carlo `0.9772 / 0.7699 / 0.2920`
  against `0.9777 / 0.7702 / 0.2919` at `θ = 0.3, 1, 2`); `κ = 0.9, 1.1`
  miss by `0.050, 0.046`; `κ = 0.5, 2` by `0.25`. **Born, iff (i) the
  hidden variable is Haar on the base `S²`, not on the fibre, and (ii) the
  coupling weight is exactly `κ = 1`.**

*What BAM supplies.* For a prepared mouth the Bloch direction is fixed and
only the fibre is unresolved: Haar on `S¹`, which by H1 and C1–C4 does not
produce Born from any coupling in the repository. Haar on `S²` is Haar on
all of `SU(2)`, the measure of an **unprepared** mouth; the only candidate
carrier is the detector's own mouth (a second unit spin with unresolved
orientation). `κ = 1` is the statement that the detector registers the
total polarisation `x + y` with equal weights. Neither (i) nor (ii) is
derived anywhere in the repository. This is Bell's 1964 / Kochen–Specker
1967 single-spin hidden-variable construction, which is known to reproduce
the single-spin Born rule and known to be a local model.

## 4. Verdict rule

* `BORN_RULE_DERIVED_FROM_CLASSICAL_BAM_MEASURE` — a coupling already in
  the repository, under fibre Haar, gives `max|f − Born| < 1e-3` with no
  tuned weight.
* `CLASSICAL_INTENSITY_ONLY_NO_OUTCOME_PROBABILITY` — the only Born-shaped
  quantity produced is an intensity fraction (C2).
* `BORN_REQUIRES_AN_IMPORTED_MEASURE_OR_DETECTOR_LAW` — Born is reached
  only by importing a measure the preparation does not supply (Haar on
  `S²`) and/or a detector law with a tuned weight (`κ = 1`) or a tuned
  basin.
* `GLOBAL_BOUNDARY_DYNAMICS_GENERATES_CONTEXTUAL_BORN_MEASURE` — the
  ensemble for `X` depends on the setting `a` through a derived global
  boundary problem and yields Born.

Expected, stated in advance: **`BORN_REQUIRES_AN_IMPORTED_MEASURE_OR_DETECTOR_LAW`**,
with the precise content that the *only* classical route to Born found is
C5, which is outcome **C** of the round's typology (deterministic hidden
outcome, probabilities from ignorance of `y`), and which therefore places
the pair problem under Bell's theorem unless measurement independence fails
geometrically. Outcome D is not expected; a candidate for it must be a
computed setting-dependence of the ensemble, not prose.

## 5. Controls

* **Basin control.** The tautological basin `|ψ| < π cos²(θ/2)` reproduces
  Born under fibre Haar; it is reported as `tuned` and must not be counted.
* **Measure control.** Replace fibre Haar by a non-invariant weight (e.g.
  `∝ 1 + cos ψ`): the C1 family changes; a result that survives only for
  a non-invariant weight fails Theorem 2.
* **Weight control.** C5 with `κ ∈ {0.5, 0.9, 1.1, 2}` must miss Born by
  the pre-computed amounts.
* **Reversal control.** Every reported `f` must satisfy
  `f(π−θ) = 1 − f(θ)` to `1e-12`.

## 6. Dependency ledger to be printed

```
f(θ) under fibre Haar  =  f( rotational covariance [derived], Haar on S¹ [derived: unresolved phase],
                             basin shape [from the coupling: derived for C1-C4; tuned for the basin control] )
C5 Born                =  ( Haar on S² [chosen: detector mouth unprepared], κ = 1 [chosen], sign(a·(x+κy)) [chosen] )
```

## 7. Deferred, explicitly

Composition (`Γ₁ × Γ₂` versus `H₁ ⊗ H₂` for an opposite-sector pair) and
CHSH are the next stage and are not computed here.
