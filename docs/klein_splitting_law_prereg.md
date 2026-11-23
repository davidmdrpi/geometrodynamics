# Pre-registration: the tunnelling law for Klein splitting on an Alice handle

Publish this file before any held-out splitting is computed. Baseline: `main`
at `85c5f68`. The model, solver and in-sample archive are PR #355 at
`07d58de` (`geometrodynamics/bulk/alice_handle_ports.py`,
`experiments/closure_ledger/runs/20261118_alice_handle_ports/`). No random
numbers are used.

**Prior state.** #355's frozen P4 split statement failed. The lowest Klein
pair at `a = 0.3`, `m = 3` splits by `Δ/ω̄ = 1.01e-8`, below the frozen
`1e-6`. Its expectation `exp(−2m asinh(cot a))` was off by a factor that grows
with `m`. This file freezes a derived law for that splitting and tests it
where it has not been looked at.

## How this hypothesis was found

The law below was formulated after reading #355's archive. Freeze and
implementation are by one author. The paths taken, in order, all on the six
in-sample points `a ∈ {0.3, 0.5}`, `m ∈ {1, 2, 3}`:

1. **A pattern in the data.** With `W = ∫√(m²/ρ² − ω̄²) dℓ` over the whole
   forbidden region (neck plus the near-rim band), `Δ·e^W` was 0.70–0.714 for
   every lowest pair. That is a pattern, not a derivation.
2. **Two derivations tried and set aside.**
   * *Landau–Lifshitz in the Liouville variable* `ψ = √ρ w` gave
     measured/predicted from 0.73 to 0.99. The extra terms of that potential
     are not small in the neck. The conformal variable below avoids them
     exactly.
   * *Landau–Lifshitz in the conformal variable*, with a harmonic-well
     prefactor, gave 1.10–1.12.
3. **The frozen form.** Herring's formula with the exact band state and WKB
   only through the neck. Its in-sample ratios are tabulated below. Nothing
   in it was fitted.

No held-out geometry has been solved. The predictions below are the formula
evaluated by `experiments/closure_ledger/klein_splitting_law_predictions.py`,
which does not call the solver.

## Derivation

**Conformal coordinate.** On the loop of #355, put `du = dℓ/ρ`. Then
`(1/ρ)(ρ w′)′ = w_uu/ρ²`, and each azimuthal sector obeys exactly

    w_uu + (ω² ρ(u)² − m²) w = 0,

with no first-derivative term. On the band `ρ = sin θ = sech u`, so the
band is a Pöschl–Teller well whose `n = 0` state is `w = sin^m θ` with
`ω² = m(m+1)`. On the neck `s = b sinh v` and `ρ = b cosh v`. At `ω = 0`,
both regions are solved exactly by `e^{±mu}`. The WKB exponent in `u`,

    ∫ √(m² − ω²ρ²) du = ∫ √(m²/ρ² − ω²) dℓ,

is therefore exact at zero frequency on both pieces.

**One well, one barrier, a Bloch phase.** The band between its turning
points is the only well on the loop. The neck is the barrier, together with
the band's forbidden strips at the rims. A torus twist `c` multiplies the
field by `e^{imc}` across rim B, which is a Bloch phase `φ = mc` on a ring of
one cell. To leading order in tunnelling, then, `ω²(φ) = ω₀² − 2t cos φ`.
#355's P3, confirmed to `3e-14`, makes the Klein pair the union of `φ = 0`
and `φ = π`, so `Δ(ω²) = 4t`.

**Herring's formula.** The weight of the eigenproblem in `u` is `ρ²`. The
Wronskian identity between the single-well state `ψ` and its image one cell
over then gives

    t = 2 ψ(u*) |ψ_u(u*)| / ∫ ρ² ψ² du,

with `u*` the neck midpoint.

**The single-well state.** Take `ψ = sin^m θ` on the band, which is exact
there. Its normalisation is `∫ sin^{2m+1}θ dθ = √π Γ(m+1)/Γ(m+3/2) ≡ N_m`.
The band is truncated at the rims, which drops a relative `O(sin^{2m+2} a)`.
From the rim, where `ψ = sin^m a`, continue `ψ` into the neck as the decaying
WKB solution with `κ = √(m² − ω²ρ²)`. Then
`ψ ψ_u |_{u*} = κ_rim sin^{2m}a e^{−W_neck}`, with

    W_neck = ∫_{−S}^{S} √(m²/(s² + b²) − ω²) ds,   κ_rim = √(m² − ω² sin² a).

**The law.** With `ω = √(m(m+1))` and `Δ = Δ(ω²)/(2ω)`:

    H(a, m) = 4 κ_rim sin^{2m}(a) e^{−W_neck} / (ω N_m).

Nothing is fitted. The old expectation is the neck factor alone at `ω = 0`.
The factor `sin^{2m} a` is the band's forbidden strips. It adds
`ln(1/sin² a)` per unit `m` to the decay: 2.43 at `a = 0.3` and 1.47 at
`a = 0.5`. #355 measured an excess over its expectation of 2.57 and 1.51.

**Expected error.** The WKB step in the neck has a relative error of order
`κ′/κ² ~ 1/m`. So `r ≡ Δ_measured/H` should behave as `r − 1 ≈ c(a)/m`,
with `c` growing with `a`. Both the neck correction and the band truncation
grow with `b = sin² a`.

**Corollary.** P4 of #355 confirmed that `L` first changes sign at
`π/(2Δ)`. The law therefore predicts the handle's reversal time,
`π/(2H)`, from the geometry alone.

## In-sample values (#355's archive, not a test)

`Δ` is `|ω₊ − ω₋|` of the lowest cosine/sine pair, `n = 0`.

| a | m | W_neck | H | r at N=2000 | r at N=4000 | m(r − 1) |
|---|---|---|---|---|---|---|
| 0.3 | 1 | 3.6791 | 4.2490e-3 | 1.0508 | 1.0509 | 0.051 |
| 0.3 | 2 | 7.4091 | 1.3186e-5 | 1.0244 | 1.0245 | 0.049 |
| 0.3 | 3 | 11.1389 | 3.4482e-8 | 1.0155 | 1.0160 | 0.048 |
| 0.5 | 1 | 2.4362 | 3.1357e-2 | 1.1821 | 1.1821 | 0.182 |
| 0.5 | 2 | 5.0282 | 8.5775e-4 | 1.0819 | 1.0819 | 0.164 |
| 0.5 | 3 | 7.6179 | 1.8834e-5 | 1.0525 | 1.0526 | 0.158 |

`Δ` itself spans six decades here, while `m(r − 1)` stays near 0.05 at
`a = 0.3` and near 0.16–0.18 at `a = 0.5`.

## Held-out predictions

The geometries were chosen so that the predicted `Δ/ω̄ ≥ 1e-8`, the smallest
split #355 resolved (to `5e-4` relative across `N`). `a = 0.2` and `m = 4`
are outside the in-sample range. `a = 0.4` lies inside it.

| a | m | W_neck | predicted Δ = H | predicted Δ/ω̄ |
|---|---|---|---|---|
| 0.2 | 1 | 4.5558 | 8.4421e-4 | 5.97e-4 |
| 0.2 | 2 | 9.1331 | 4.9980e-7 | 2.04e-7 |
| 0.4 | 1 | 3.0080 | 1.3262e-2 | 9.38e-3 |
| 0.4 | 2 | 6.1111 | 1.3728e-4 | 5.60e-5 |
| 0.4 | 3 | 9.2134 | 1.1766e-6 | 3.40e-7 |
| 0.5 | 4 | 10.2072 | 3.8280e-7 | 8.56e-8 |

Diagnostics, reported but not gated: `a = 0.6`, `m = 1…5`, predicted
`Δ = 5.9819e-2, 3.8333e-3, 1.8530e-4, 8.2061e-6, 3.4798e-7`. This angle is
outside the in-sample range of `c(a)`, so its tolerance is not known.

## Measurement

Use #355's solver at `07d58de`, unchanged. For each geometry, build
`Loop(a, N)` and `Sector(loop, m, KLEIN)`. Take the lowest eigenvalue of the
cosine sector and of the sine sector, as paired in `p4_p5`. These are the
Rayleigh-quotient eigenvalues #355 uses. Compute `Δ = |ω₊ − ω₋|` at
`N = 2000` and `N = 4000`.

**Reproduction check, before any held-out number is read.** The
implementation must reproduce the six in-sample `Δ` of #355's archive to a
relative `1e-6` at both sizes. If it does not, the run stops and is reported
`INCOMPLETE`.

## Gates

* **G1, the law (structural).** At every gated geometry,
  `0 < m(r − 1) ≤ 0.25`. That is, `H` underestimates `Δ`, by at most
  `25%/m`.
* **G2, the error ordering (structural).**
  * `c(0.2) < c(0.4)`, where `c` is `m(r − 1)` averaged over the gated `m`
    at that angle.
  * Every gated `m(r − 1)` at `a = 0.4` lies in `[0.048, 0.182]`, the
    in-sample range spanned by `a = 0.3` and `a = 0.5`.
* **G3, resolution (tolerance).** For every gated geometry,
  `|Δ_2000/Δ_4000 − 1| ≤ 1e-3`.

## Verdicts

* `CONFIRMED` if G1, G2 and G3 hold.
* `COUNTEREXAMPLE` if G1 or G2 fails at both resolutions. The report names
  which gate failed.
* `UNRESOLVED_NUMERICALLY` if G3 fails, or if a G1 or G2 miss shrinks from
  `N = 2000` to `N = 4000`. One further doubling to `N = 8000` is allowed,
  and nothing else changes.
* `INCOMPLETE` if the reproduction check fails.
* `NO_VERDICT` for any case not named here.

Tolerances are not adjusted after the fact, and a miss is reported as a miss.
The report also gives, without gates:
* the `a = 0.6` diagnostics;
* `r` for the excited pairs `n = 1, 2`, using the exact band states
  `P_{m+n}^m(cos θ)`;
* the predicted reversal times `π/(2H)`.

## Scope

* A classical scalar on a fixed 2D surface, as in #355. Not a GR solution
  and not on `S³`.
* The splitting is classical wave tunnelling through a geometric barrier.
  The resemblance to a quantum tunnel splitting is mathematical. Nothing
  here is claimed to be quantum.
* A `CONFIRMED` verdict would explain #355's P4 counterexample and give the
  reversal time in closed form. It would not reopen #355's second gate, the
  dynamical reflection axis, which its freeze tied to a `CONFIRMED` verdict
  there.
