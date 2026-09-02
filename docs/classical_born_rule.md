# Derive or refute the Born rule from classical BAM measure

*Fourth round of the finite-mouth chain. Pre-registered in
`docs/classical_born_prereg.md` (`7ff6e41`); module
`geometrodynamics/bulk/mouth_measurement.py`; tests
`tests/test_mouth_measurement.py`; probe `classical_born_probe.py`. The
spin-transport arc is frozen at `7bd7ecc`. Composition and CHSH are deferred
and not computed here.*

## Verdict

Pre-registered label: `BORN_REQUIRES_AN_IMPORTED_MEASURE_OR_DETECTOR_LAW`.
Scope-correct reading, adopted after review and recorded as a correction
in the pre-registration:

`CURRENT_BAM_PREPARATION_AND_DETECTOR_DYNAMICS_DO_NOT_DERIVE_BORN`

The pre-registered label reads as a no-go against every classical BAM
detector, and the round did not establish that. What it established:

```
Spin-frame symmetry does not select the Born basin.
Classical Malus intensity does not explain outcome probability.
The existing BAM winding detector does not derive Born.
A classical detector hidden-variable model reproduces Born exactly,
but its physical detector ensemble and coupling are not yet derived.
```

With the classical state the previous rounds actually supply — a point of
the mouth spin-frame bundle, Bloch direction fixed by the preparation, fibre
phase unresolved with Haar measure — no deterministic detector built from a
coupling in the repository produces Born outcome frequencies. The one
classical route that does (C5) needs a hidden variable Haar-distributed on
the base `S²`, which a prepared mouth does not have, and a coupling weight
tuned to exactly `κ = 1`. It is Bell's 1964 / Kochen–Specker 1967
single-spin hidden-variable model: outcome typology **C**, deterministic
hidden outcome with probabilities from ignorance.

## 1. The question

Classical intensity is quadratic in amplitude; that is not the Born rule.
The question is whether **individual outcomes** are sampled with `|A|²`
frequencies from a classically derived ensemble and a deterministic
detector. Three theorems were fixed in advance: intensity is not
probability; the measure must be invariant, and a basin tuned to the answer
counts as a tuned measure; the detector must be deterministic.

## 2. The classification (H1)

Rotational covariance reduces any deterministic detector to a function
`D(θ, ψ)` of the polar angle `θ = ∠(a, x)` and the azimuth `ψ` of `a` in the
tangent frame at `x`. Under the fibre `q ↦ e^{iφ}q` the frame turns by `2φ`
(`mouth_spin_frame`), so `ψ` shifts by `2φ` (checked: `θ` fixed to `1e-12`,
`|Δψ| = 2φ` to `1e-10`) and fibre Haar makes `ψ` uniform. Therefore

```
P(+|θ) = f(θ) = arc measure of {ψ : D(θ, ψ) = +1} / 2π.
```

Analyzer reversal forces `f(π−θ) = 1 − f(θ)`; composed with complementarity
of orthogonal outcomes it adds `D(θ, ψ) = D(θ, π−ψ)`, and nothing else.
Born `cos²(θ/2)`, the straight line `1 − θ/π` and the step all satisfy the
first (residuals `≤ 2e-16`), and each is realised by the basin centred at
`ψ = π/2`, `|ψ − π/2|_circ < π f(θ)`, which has both detector-level
symmetries (checked directly: zero violations on a `720`-point circle
across `θ`) and the required arc measure (residuals `1e-4`). *Correction:*
the first version used the basin `|ψ| < πf`, centred at `0`, which lacks
the `π−ψ` symmetry; the conclusion was unaffected, the constructive proof
was not what the pre-registration claimed. **Symmetry plus fibre Haar do
not select Born. The basin shape decides, and the basin shape must come
from a coupling.**

*On the fibre ensemble.* The fibre is the Spin(2) spin-frame fibre. If its
coordinate is gauge, averaging over it is gauge averaging and no detector
can obtain stochastic outcomes from it at all; if it is physical, Haar is
the natural invariant measure but its emergence from preparation dynamics
is not shown anywhere in the repository. The round grants the most natural
invariant ensemble, and Born still does not follow.

## 3. The couplings (H2)

| candidate | coupling | induced `f(θ)` | `max_θ |f − cos²(θ/2)|` | class |
|--|--|--|--|--|
| C1 | sign of a linear functional of the frame, `sign(α cos θ + ρ sin θ cos(ψ−ψ₀))` (torque / gradient force) | `arccos(−(α/ρ) cot θ)/π` clipped: plateaus at `0`, `1` | `0.109` at best (`α/ρ = 0.80`); closed form against the arc measure `3e-5` | analytic identity |
| C2 | classical Malus intensities `|⟨a,±|q⟩|²` | phase-independent: the step `sign(cos θ)` | `0.500` | analytic identity; **Theorem 1**: the intensity fraction *is* `cos²(θ/2)` |
| C3 | spinor (`e^{iφ}`) plus frame (`e^{2iφ}`) harmonics, four natural weightings | irregular | `≥ 0.500` | numerically converged |
| C4 | the repository's winding Stern–Gerlach with the fibre-`U(1)` setting (PR #238) | winding-diagonal, phase-independent: the step | `0.500` | analytic identity |
| C5 | detector carries a Haar-distributed direction `y ∈ S²`; `sign(a·(x + κy))` | `clip((1 + cos θ/κ)/2)` | `κ=1`: **`0`**; `κ = 0.9, 1.1`: `0.050, 0.046`; `κ = 0.5, 2`: `0.25` | analytic identity (Archimedes); conditional |

C5 works because `a·y` is uniform on `[−1, 1]` for Haar `y` (Archimedes'
hat-box theorem; Kolmogorov distance `1.5e-3` on `4×10⁵` samples), so
`P(a·x + κ a·y > 0) = (1 + cos θ/κ)/2`. Monte Carlo at `κ = 1`:
`0.9776, 0.7701, 0.2913` against Born `0.9777, 0.7702, 0.2919` at
`θ = 0.3, 1, 2`.

On C4, one more fact from the repository itself: its own measurement round
supplies outcome statistics for a superposed preparation from an explicitly
*Born* throat ensemble (`measurement_sector_probe.sg_run`: "the winding
Stern–Gerlach with a Born throat ensemble"). The repository's detector does
not derive Born; it is fed Born.

## 4. What BAM supplies, and what C5 needs

A prepared mouth has its Bloch direction fixed; its only unresolved variable
is the fibre phase, Haar on `S¹`. Under that measure, H1 and C1–C4 show that
no coupling in the repository gives Born, and the only way to get it is the
tuned basin `|ψ| < π cos²(θ/2)` (basin control: reproduces Born to `1e-4`;
classified **tuned**).

C5 needs Haar on `S²`, which is the base marginal (the pushforward under
the Hopf map) of Haar on `SU(2) ≅ S³`: the measure of an **unprepared**
mouth (checked: the pushforward of Haar `q ∈ S³` gives `a·y` uniform to a
Kolmogorov distance `< 1.5e-2` on `2×10⁴` samples). The only candidate
carrier is the detector's own mouth, a second unit spin with unresolved
orientation, and the coupling must register the total polarisation `x + y`
with equal weights. Neither the detector-mouth model nor `κ = 1` is derived
in the repository; both are exactly the ingredients of the classic
single-spin hidden-variable construction.

**C5 is an open derivation route, not a closed one.** If the detector
contains an identical unprepared mouth, isotropy gives `y` Haar on `S²`
without further input; if system and detector mouths enter a symmetric
classical polarisation coupling, equal normalisation could force `κ = 1`.
Then the single-system Born rule would be derived from apparatus microstate
ignorance. It would remain a deterministic classical completion of one spin
and would not touch Bell's theorem.

## 5. Controls

* **Basin control.** The tautological Born basin reproduces Born (`1e-4`);
  it is reported as tuned and not counted.
* **Measure control.** A non-invariant fibre weight with the detector's
  period, `1 + 0.8 cos(2φ + 1)`, changes the C1 curve; any result that needed
  such a weight would fail Theorem 2. *Correction note:* the first version
  of this control used `1 + cos φ`, which has twice the period of every
  frame-coupled detector and therefore integrates to zero against it
  (change `9e-5`); the control could not fire. The weight was changed to the
  detector's period. This changed the control, not a criterion.
* **Weight control.** C5 with `κ ≠ 1` misses by the pre-computed amounts.
* **Reversal control.** Every reported curve satisfies `f(π−θ) = 1 − f(θ)`
  to `1e-12`.

## 6. Typology and consequence

Outcome **C**: for every complete classical state `(x, y)` the outcome
`D_a` is fixed; probabilities come from ignorance of `y`. Outcome **D**
(a setting-dependent ensemble from a derived global boundary problem) was
not found, and nothing in the repository supplies a computation of it. This
matters for what comes next: a C-type model for a pair with independent
detector variables is a local hidden-variable model and falls under Bell's
theorem. The repository's own PR #206 record (`max CHSH = 2` for a single
classical field with local dynamics) is that theorem applied. Escaping it
would require measurement independence or outcome independence to fail
*geometrically*, which is the content the D branch would have to supply.

## 7. Dependency ledger

```
f(θ) under fibre Haar  =  f( rotational covariance [derived],
                             Haar on S¹ [natural invariant measure; gauge-or-physical fork and
                                         preparation derivation open],
                             basin shape [from the coupling: derived for C1-C4; tuned for the control] )
C5 Born                =  ( Haar on S² = base marginal of Haar on S³ [open route: identical unprepared
                                                                   detector mouth + isotropy],
                            κ = 1 [open route: symmetric polarisation coupling],
                            D = sign(a·(x + κy)) [chosen] )
```

## 8. Single-system Born and Bell have separated

Even if C5 is derived, it is a deterministic classical completion of one
spin, exactly the account Bell's 1964 paper grants a single spin-½ before
proving the two-particle no-go. The pair problem must then explain why the
joint global solution is not `λ = (λ_source, y_A, y_B)` with `A = A(a, λ)`,
`B = B(b, λ)` under a setting-independent distribution; if it is, CHSH is
capped at `2`. The next target is therefore sharper than "composition":

> Does the BAM global boundary-value problem, solved with `(a, b)` as
> boundary data, produce an admissible classical ensemble
> `ρ(λ | a, b) ≠ ρ(λ)`?

Not an assumed retrocausal weight, not a transaction product: the classical
boundary problem actually solved, and the resulting measure compared with
the setting-independent one. If `ρ(λ|a,b) = ρ(λ)`, BAM with deterministic
detectors is in Bell's local class regardless of the mouth topology. If
not, the questions become whether the dependence arises from two-time
boundary conditions without signalling and whether it reproduces the
quantum joint probabilities.

## 9. What is not claimed

Nothing about composition (`Γ₁ × Γ₂` versus `H₁ ⊗ H₂`) or CHSH; nothing
about the field sector beyond the spin-frame degree of freedom; nothing
about whether a global boundary-value formulation of BAM could produce a
setting-dependent ensemble. Those are the next questions, and the last one
is the only classical route left that is not already closed by this round
plus Bell's theorem.
