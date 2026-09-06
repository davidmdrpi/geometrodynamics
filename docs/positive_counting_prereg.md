# Pre-registration: can positive classical history counting reach the quantum correlation?

**Status: frozen before any code for the open questions.** Ninth round,
following PR #283 (`b0e372f`). Prompted by a pattern across rounds 5–8 and
#283: every concrete mechanism anyone has built lands *near but not on*
quantum mechanics — positive counting `2.1423`, quadratic readouts `3.3941`
and `3.7712`, oriented sum `2√2` available but unselected. This round asks
whether that near-miss is systematic.

## 0. The class, stated sharply

A **positive local counting model** on the closure geometry assigns to each
outcome sector `(s_A, s_B)` the weight

```
W_s = ∫_Γ Φ(D_s(x)) dσ,    Φ ≥ 0,
D_s(x) = 1 + x·u + u·w + x·w,   u = s_A a,  w = −s_B b,
```

with `Γ` the closure circle and a sector prior `π_s`. This is the class the
repository's own mechanisms live in: round 5's positive branch is `Φ = |·|`,
and #283's canonical-equilibrium low-temperature limit is the same `Φ = |·|`.

**What the class excludes, stated so the result is not overread.** Weights
that are nonlinear functionals of the integral rather than integrals of a
function of `D` — round 7's `(∫_Γ D dσ)²` is *not* in the class, though
`∫_Γ D² dσ` is. Weights depending on `x` other than through `D_s(x)`.
Sector-dependent `Φ`. Any of these escapes a no-go proved here.

## 1. Results established before freezing (analytic oracle)

Verified numerically before this file was committed; frozen as oracles.

* **O1 — the reduction.** On `Γ`, with `t = 1 + u·w`,

  ```
  D_s(x) = t + √(2t) cos ψ,     ψ uniform in arclength.
  ```

  Checked as an equality of laws (sorted residual `2.7e-5`, grid-limited)
  over 16 `(γ, sector)` cases.
* **O2 — the sector pair.** `t_like = 1 − cos γ`, `t_unlike = 1 + cos γ`, so
  `t_like + t_unlike = 2` identically.
* **O3 — one function of one variable.** Everything depends on `Φ` only
  through `W(t) = ∫₀^{2π} Φ(t + √(2t) cos ψ) dψ`, and with equal priors

  ```
  E(γ) = [W(τ) − W(2−τ)] / [W(τ) + W(2−τ)],    τ = 1 − cos γ.
  ```

  `Φ = |·|` reproduces round 5 (`E(1) = −0.398496650`) and `Φ = id`
  reproduces round 6 (`E = −cos γ`), both to `1e-9`.
* **O4 — the quantum target as a functional equation.** `E = −cos γ` for all
  `γ` is equivalent to `W(τ)/W(2−τ) = τ/(2−τ)`, i.e. to

  ```
  W(t) = t · G(t)   with   G(t) = G(2−t).
  ```

  In particular `W(t) ∝ t` suffices, and is what `Φ = id` achieves
  (`∫_Γ D dσ = 2πt`, round 6).
* **O5 — marginals are automatic.** `Γ` is invariant under `x → −x`, which
  exchanges `(s_A,s_B)` with `(−s_A,−s_B)` while preserving `D`. Every model
  in the class therefore has marginals exactly `1/2` with paired priors, so
  **no-signalling constrains nothing here** and cannot be cited as support.
* **O6 — a complex-analytic form.** `D + 1/2 = ½|1 + r e^{iψ}|²` with
  `r = √(2t)`, so `W(t)/2π` is the mean of the radial function
  `Ψ(|z|²) = Φ((|z|²−1)/2)` over the circle of radius `r` centred at `1`.

## 2. The two questions

**Q1 — does positivity bound CHSH at all?** Compute
`sup CHSH` over `Φ ≥ 0`.

*Prediction, frozen:* the supremum is **4**, the algebraic maximum, and is
approached by a narrow bump at `v = 1 + √2`. Reason: `W(t) > 0` iff
`v ≤ t + √(2t)`, whose right side increases in `t`, so a bump at the value
`t + √(2t)` takes at `t = 1` gives `W = 0` for `t < 1` and `W > 0` for
`t > 1`, hence `E(γ) = −sgn(cos γ)` and `S = 4` at the standard angles.
*Falsifier:* the construction fails, or the attained supremum is below `4`.

**Q2 — is the quantum law attainable?** Is there `Φ ≥ 0` with
`W(t) = t·G(t)`, `G(t) = G(2−t)`? Equivalently (O6): is there a nonnegative
radial `Ψ` whose mean over the circle of radius `r` centred at `1` equals
`r²/2` for every `r ∈ (0,2)`?

Permitted verdicts:
`QUANTUM_LAW_ATTAINABLE_BY_POSITIVE_COUNTING` (an explicit `Φ ≥ 0` exhibited),
`QUANTUM_LAW_UNATTAINABLE_IN_THIS_CLASS` (with a dual infeasibility
certificate), or `UNRESOLVED_NUMERICALLY`.

## 3. Rules fixed in advance

1. **A small residual is not attainability.** Nonnegative least squares
   returning a near-zero residual does **not** discharge Q2; only an explicit
   `Φ ≥ 0`, verified to reproduce `E = −cos γ` to `1e-6` at independent
   angles, does.
2. **A large residual is not a no-go.** A negative answer requires a **dual
   certificate**: a signed measure `λ(r)` with `∫λ(r) dμ_r(y) ≤ 0` for every
   `y` in the support while `∫λ(r)·(r²/2) dr > 0`. Report the certificate and
   verify it independently of the primal solve.
3. **The discretisation must be shown not to drive the answer.** Refine the
   grid in both `y` and `r` and report the trend.
4. Neither outcome may be preferred. A positive answer to Q2 would make the
   quantum branch reachable by counting and would *weaken* round 6's fork; a
   negative answer is a no-go that does not by itself select any mechanism.
   Both are reportable results.
5. Nothing in rounds 5–8 or #283 may change. A test asserts the round-5,
   round-6 and round-7 numbers are unmoved.

## 4. Expected verdict, stated in advance

For Q1, `4`. For Q2 I have no expectation and am not recording one; the
class is large enough to exceed Tsirelson (Q1), which makes it a genuinely
open question whether it can *hit* the quantum value.

If Q1 gives `4` and Q2 is negative, the statement is that this model class
can **exceed** the quantum bound but cannot **attain** the quantum law —
which would explain the observed near-misses as structural rather than
accidental, and would be the round's result.

## 5. Dependency ledger to be printed

Including at minimum: the class restriction of §0, the sector prior, the
geodesic itinerary and realignment inherited from rounds 5–6, and the
conditioning variable inherited from round 8.
