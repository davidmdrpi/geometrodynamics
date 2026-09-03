# Pre-registration: does the BAM two-boundary classical solution induce measurement dependence at the source, without an imported measure or signalling?

**Status: frozen before any code.** Fifth round of the finite-mouth chain,
after `docs/classical_born_prereg.md` (merged in PR #279). Success criteria
may change afterwards only by an explicit correction note.

## Correction note (post-implementation, before commit of the module)

P1 stated that `ρ(x|a,b)` and `ρ(x|a′,b′)` have distinct great-circle
supports, hence total variation `1`. That is true for **non-coplanar**
settings (computed: `TV = 1`). In the standard Bell configuration all
analyzers lie in one plane, so every setting pair shares the same great
circle and the setting dependence is in the **density** on it, not in the
support: total variation `0.018` between the pairs `(0, 1)` and `(0.5, 1)`,
`0.063` between `(0, π/4)` and `(π/2, −π/4)`. Measurement independence
fails in both cases; the correction is to the described mechanism in the
coplanar case. No criterion changes.

## 0. What the repository's two-boundary structure is

The repository has one global, two-boundary principle: the closed-history
axiom of `history/closure.py` — *only globally closed, conservation-respecting
histories become events*, with **phase closure `≡ 0` or `π` (mod `2π`)**
around every loop. Its Bell computation (`enumerate_bell_branches`) weights
each outcome branch by `|amplitude|² × closure_weight`, where the amplitude
is the singlet amplitude imported from `bell.analyzers` and the closure
weight is a **Gaussian in the phase mismatch with a chosen width
`σ = 0.6`**. Both factors are imported; that computation derives nothing
about measurement dependence and is audited here as such. The transactional
handshake (`transaction/handshake.py`) is phenomenological by its own
docstring. The MTY network (PR #216/#276) needs a traversable throat the
GR rounds closed.

This round therefore builds the closure principle **without** the imported
factors: no amplitude, no Gaussian, no width. The classical state is the
one the spin-frame rounds supply.

## 1. The model (chosen, stated first)

* **Source variable.** The created pair's frame direction `x ∈ S²` (the
  Bloch direction of the mouth spin-frame; the partner carries `−x` through
  the Pin transport). Prior: Haar on `S²` (isotropy of an unprepared pair;
  invariant). *Chosen: that the pair's direction is physical, not gauge.*
* **The closed history.** `source → A → B → source`, a loop of frame
  directions on the Bloch sphere: `x → u → v → x` with `u = s_A a`,
  `v = s_B b`, the outcome signs being **boundary data** of the history,
  and each leg the geodesic (parallel transport). *Chosen: that detection
  with setting `a` and outcome `s_A` realigns the frame to `s_A a` along
  the geodesic; this is the classical content of "absorption into the
  analyzer".*
* **The closure phase.** The spin-½ geometric phase of the loop,
  `Ω/2` with `Ω` the solid angle of the geodesic triangle:

  ```
  Ω(x; u, v) = 2 atan2( x·(u×v),  1 + x·u + u·v + v·x )  ≡  2 atan2(N, D)
  ```

  Derived: it is the holonomy of the spin-frame bundle around the loop,
  the object the previous rounds identified as the Hopf/Spin(2) phase.
* **Closure.** `Ω/2 ≡ 0 or π (mod 2π)`, i.e. `Ω ≡ 0 (mod 2π)`, i.e.
  `N = 0`: **`x` lies on the great circle `Γ` through `a` and `b`.** The
  condition depends on both settings jointly. No width parameter.
* **The conditioned measure.** Haar on `S²` conditioned on `N = 0` is the
  coarea measure on `Γ`, density `1/|∇Ω| = |D| / (2|u×v|)` (derived:
  `∇Ω = 2∇N/D` on the zero set, `|∇_{S²}N| = |u×v|` there). This is the
  `ε → 0` limit of "uniform on `{|Ω mod 2π| < ε}`", which is the
  independent construction used to check it.
* **Outcome probabilities.** `P(s_A, s_B | a, b) ∝ ∫_Γ |D| dσ / (2|u×v|)`:
  the relative closure measure of the four histories. *This is not a
  deterministic detector: for `x ∈ Γ` several histories close, and the
  outcome is which closed history is realised. That is the D-type
  structure by construction, and it is stated as such.*

## 2. Pre-computed predictions (closed form, derived by hand)

With `γ = ∠(a, b)`, `c = cos(γ/2)`, `s = sin(γ/2)`, and
`∫₀^{2π}|k + cos τ| dτ = 2πk + 4√(1−k²) − 4k arccos k`:

```
W(+,+) = W(−,−) = 2 + c(π−γ)/s          W(+,−) = W(−,+) = 2 + sγ/c

E(γ) = [c²(π−γ) − s²γ] / [2 sin γ + c²(π−γ) + s²γ]
```

* **P1 (measurement dependence).** `ρ(x|a,b)` is supported on `Γ(a,b)`;
  `ρ(x|a′,b) ≠ ρ(x|a,b)` whenever `Γ` changes. Total-variation distance
  between the conditioned measures for two setting pairs is `1` (disjoint
  supports up to two points) — measurement independence is violated
  maximally. *Falsifier:* none needed; this is the definition doing its
  work, and it is reported as such, not as a discovery.
* **P2 (no signalling).** `P(s_A|a,b) = 1/2` for every `b`, exactly, because
  `W(+,+) + W(+,−) = W(−,−) + W(−,+)` (the weight is invariant under
  `(x,u,v) → (−x,−u,−v)`). *Falsifier:* a marginal depending on the
  distant setting at the `1e-12` level.
* **P3 (the correlation).** `E(γ)` as above; numerically
  `E(0.3) = 0.82095`, `E(π/4) = 0.53557`, `E(1) = 0.39850`,
  `E(π/2) = 0`, `E(2) = −0.30350`, with `E(π−γ) = −E(γ)`. It is **not**
  `cos γ` (`0.95534, 0.70711, 0.54030, 0, −0.41615`). *Falsifier:* a
  quadrature or Monte Carlo value off the closed form by more than `1e-4`
  (quadrature) or `3e-3` (Monte Carlo, `2×10⁶` samples).
* **P4 (Bell violation).** At `(a, a′, b, b′) = (0, π/2, π/4, −π/4)`,
  `S = 3E(π/4) − E(3π/4) = 4E(π/4) = 2.1423 > 2`. The maximum over
  settings is predicted **strictly below `2√2 = 2.8284`** and is to be
  reported. *Falsifier:* `S_max ≥ 2√2`.
* **P5 (the signed variant is Born).** Replacing `|D|` by the signed `D`
  gives `E = (c²−s²)/(c²+s²) = cos γ` exactly and `S = 2√2`. Identity:
  `D/4 = Re Tr(P_x P_u P_v)`, the real part of the Bargmann invariant of
  the three spin-½ projectors (checked to `2e-16`). So **the quantum
  correlation is the signed closure measure**, which assigns negative
  weight to the arc of `Γ` where `D < 0` (angular size `γ` for like
  outcomes). A signed measure is not a probability: using it imports the
  quantum rule. *Class:* analytic identity; the quasi-probability reading is
  recorded, not adopted.
* **P6 (the strict-zero variant).** Requiring `Ω = 0` rather than
  `Ω ≡ 0 mod 2π` (dropping the `π` branch of the closure axiom) keeps only
  the `D > 0` arc: `E(1) = 0.46495`, `S = 2.4649` at the standard angles.
  Still not `cos γ`. Reported for completeness; the repository's axiom
  includes the `π` branch.
* **P7 (window convergence).** Uniform sampling of `x ~ Haar(S²)` with
  the window `|Ω mod 2π| < ε` converges to the coarea `E(1) = 0.39850`:
  pre-computed `0.3729, 0.3923, 0.3973, 0.3972` at
  `ε = 0.4, 0.2, 0.1, 0.05` (`2×10⁶` samples). This is the independent
  construction of the conditioned measure and confirms `|D|`, not `D`.

## 3. Controls

* **Loop-topology control.** A history that does not link both detectors
  (`x → u → x`, the two-leg loop) has zero solid angle for every `x`:
  closure is automatic, `ρ(x|a,b) = ρ(x)`, and no correlation is
  selected. The measurement dependence comes entirely from the loop
  passing through both settings.
* **Local-detector control.** The previous round's C5 detectors
  (`sign(a·(x + y_A))`, `sign(b·(−x + y_B))`, independent `y_A, y_B` Haar)
  on the same source prior give a local model: `S ≤ 2` (predicted
  `E = −cos γ /2`... the value is to be computed, the bound is the
  prediction). *Falsifier:* `S > 2` from a local model.
* **Width control.** The repository's Gaussian closure weight with
  `σ = 0.6` replaces the sharp condition; `E` changes and depends on `σ`.
  A result that needs a chosen `σ` is imported.
* **Sign control.** `v → −v` (the partner's antipode) flips `E → −E`,
  singlet versus triplet; `S` is unchanged.
* **Symmetry control.** Marginals `1/2` at random settings; `E(π−γ) = −E(γ)`.

## 4. Verdict rule

* `CLOSURE_INDUCES_SETTING_DEPENDENT_SOURCE_MEASURE_NO_SIGNALLING_NOT_BORN` —
  P1, P2, P4 hold and P3 gives `E ≠ cos γ` with `S_max < 2√2`. **Expected.**
* `CLOSURE_REPRODUCES_QUANTUM_CORRELATIONS` — the positive coarea measure
  gives `E = cos γ` to `1e-4`.
* `CLOSURE_INDUCES_NO_MEASUREMENT_DEPENDENCE` — `ρ(x|a,b) = ρ(x)`.
* `CLOSURE_SIGNALS` — a marginal depends on the distant setting.

A narrower verdict is permitted if the mathematics shows the trichotomy
mis-posed.

## 5. What "without an imported measure" will mean

The prior is Haar (invariant); the conditioning is the coarea limit of a
window, parameter-free; the closure axiom is the repository's own. What is
**chosen** is the model of detection as geodesic realignment to the outcome
direction, and that the outcome signs are history boundary data rather than
functions of the source state. Every headline is stated conditionally on
those two choices.

## 6. Dependency ledger to be printed

```
ρ(x|a,b)  =  ρ( Haar on S² [invariant prior; pair direction physical: chosen],
                closure axiom Ω ≡ 0 mod 2π [repository axiom],
                geodesic-realignment detection model [chosen],
                coarea conditioning [derived; window limit] )
E(γ)      =  E( ρ(x|a,b) [above], outcome signs as boundary data [chosen: D-type] )
Born      =  signed closure measure [identity: Re Bargmann]  — not a probability; imported if used
```
