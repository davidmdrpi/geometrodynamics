# Pre-registration: the exact-neck law for Klein splitting on an Alice handle

Publish this file before any held-out splitting is computed.

* **Baseline.** Branch `claude/hopeful-babbage-ps703o` at `c40907f`, which is
  PR #362 on top of `main` at `1d29e3e`.
* **Model and solver.** PR #355's `geometrodynamics/bulk/alice_handle_ports.py`,
  unchanged from `07d58de`.
* **Prior archives.** `experiments/closure_ledger/runs/20261118_alice_handle_ports/`
  (#355) and `experiments/closure_ledger/runs/20261123_klein_splitting_law/`
  (#362).

No random numbers are used.

**Prior state.** #362 froze the Herring law `H` at `56055ab`, which takes a WKB
step through the neck. It was `CONFIRMED` on six held-out cells, with `H`
underestimating the splitting by 1–10%. The `a = 0.6` diagnostics reached
33%. #362's write-up attributed the error to the WKB step. It said that solving
the neck exactly should remove the error. This file freezes that law and tests
it on cells no one has computed.

## How this hypothesis was found

Freeze and implementation are by one author. All 17 cells whose splitting
has been measured are in-sample here:

* the six cells of #355;
* #362's six gated cells;
* #362's five `a = 0.6` diagnostics.

The steps, in order:

1. **The review of #362 made the substitution post hoc.** It kept `H`'s band
   state, frequency and normalisation, and integrated the neck equation in
   place of the WKB step. It reported these worst discrepancies:
   * 0.234% on #362's gated cells;
   * 1.02% at `a = 0.6`.

   The review did not specify the neck solution, so steps 2 and 3 had to fix
   one.
2. **Forward integration from the rim, set aside.** Starting at the rim with
   the band's value and slope and integrating to the neck midpoint does not
   give a usable law. For `m ≥ 2` the growing neck solution dominates, and
   `Δ_measured/prediction − 1` is near −1.
3. **Three far-rim conditions, compared in-sample.** The neck solution was
   normalised to the band value at the near rim. Three conditions at the far
   rim were tried:

   | condition at the far rim | worst \|r − 1\|, #362 gated | worst \|r − 1\|, `a = 0.6` | worst \|r − 1\|, #355 cells |
   |---|---|---|---|
   | `D = 0` (frozen) | 0.233% | 1.017% | 0.49% |
   | `D′/D = −κ_rim` | 0.50% | 2.9% | 1.25% |
   | `D′/D = −m cos a` | 0.48% | 2.2% | 1.11% |

   * For `m ≥ 2` the three agree to `3e-4` or better. They differ most at
     `m = 1`, where the neck action is smallest.
   * Only `D = 0` reproduces the review's figures. The review gives the
     first as 0.234%.
   * `D = 0` is the standard choice for a state localised in one well: it
     carries no component into the image well. It was also the best of the
     three in-sample. That is a selection on known cells, and it is disclosed
     here.
4. **The tolerance was set after the in-sample residuals were seen.** It is
   twice their envelope, plus a resolution floor (see Gates).

No held-out geometry has been solved. The predictions below come from
`experiments/closure_ledger/klein_exact_neck_predictions.py`, which does not
call the finite-element solver.

## The law

The setup is as in `docs/klein_splitting_law_prereg.md`:
* the conformal coordinate `du = dℓ/ρ`;
* one well, the band, and one barrier, the neck with the band's forbidden
  strips;
* the torus twist as a Bloch phase;
* `Δ(ω²) = 4t`, with Herring's `t = 2ψψ_u/∫ρ²ψ²` taken at the neck midpoint.

**The neck in its own conformal coordinate.** Put `s = b sinh v`, with
`b = sin² a`. Then `ρ = b cosh v`, `du = dv`, and the neck is
`v ∈ [−L, L]`, where `L = asinh(cot a)` and `b cosh L = sin a`. The sector
equation there is exactly

    D_vv = (m² − ω² b² cosh² v) D.

**The single-well state.**
* On the band, `ψ = sin^m θ`, exact at `ω² = m(m+1)`.
* In the neck, `ψ = sin^m(a) D(v)`. Here `D` solves the neck equation at the
  same `ω`, with `D(L) = 0` at the far rim and `D(−L) = 1` at the near rim.
* `ψ` is continuous at the rim. Its slope there, `D′(−L)`, matches the band's
  `−m cos a` to 0.008–1.2% on the in-sample cells. The residual kink is one
  source of error.
* On the in-sample cells the kink and `r − 1` are close in size. For
  example they are 1.18% and 1.02% at (0.6, 2), and 0.51% and 0.49% at
  (0.5, 2). This was noticed after the fact and is not used.

**The law.** With `ω = √(m(m+1))` and `N_m = √π Γ(m+1)/Γ(m+3/2)`:

    E(a, m) = 4 sin^{2m}(a) |D(0) D′(0)| / (ω N_m).

`H` is the same expression with the WKB value
`D(0)D′(0) ≈ −κ_rim e^{−W_neck}`. Nothing is fitted.

**Limit.** As `a → 0`, `E` tends to #360's leading law `C(m) tan^{4m}(a/2)`.
The relative difference is `O(a⁴)`: 7.8e-6 at `a = 0.04`, `m = 1`.

**Expected error, empirical.** No error order is derived. Three
approximations remain:
* the slope kink at the rim;
* the normalisation `N_m`, which is the full-sphere integral; it ignores both
  the truncated band caps and the neck's own weight;
* Herring's leading order in the tunnelling.

In-sample, `|r − 1| ≤ 0.11 a⁴` at N = 4000 in all 17 cells, where
`r = Δ_measured/E`. At small `a` this is comparable to the solver's
resolution floor of about `7e-4`. `r − 1` is
positive at N = 4000 in 16 of the 17 cells. The exception is (0.2, 2), at
−1.4e-4, where `|Δ₂₀₀₀/Δ₄₀₀₀ − 1|` is 7.2e-4.

## In-sample values (#355 and #362 archives, not a test)

| a | m | E | r, N=2000 | r, N=4000 | r_H − 1, N=4000 | \|r − 1\| / \|r_H − 1\| | \|r − 1\| / bound |
|---|---|---|---|---|---|---|---|
| 0.2 | 1 | 8.6186e-4 | 1.00009 | 1.00016 | 0.0211 | 0.008 | 0.12 |
| 0.2 | 2 | 5.0494e-7 | 0.99914 | 0.99986 | 0.0101 | 0.091 | 0.65 |
| 0.3 | 1 | 4.4612e-3 | 1.00087 | 1.00089 | 0.0509 | 0.017 | 0.34 |
| 0.3 | 2 | 1.3503e-5 | 1.00037 | 1.00052 | 0.0245 | 0.021 | 0.20 |
| 0.3 | 3 | 3.5029e-8 | 0.99961 | 1.00013 | 0.0160 | 0.025 | 0.15 |
| 0.4 | 1 | 1.4558e-2 | 1.00233 | 1.00233 | 0.1003 | 0.023 | 0.38 |
| 0.4 | 2 | 1.4348e-4 | 1.00186 | 1.00191 | 0.0471 | 0.040 | 0.31 |
| 0.4 | 3 | 1.2115e-6 | 1.00081 | 1.00098 | 0.0306 | 0.032 | 0.16 |
| 0.5 | 1 | 3.6938e-2 | 1.00351 | 1.00351 | 0.1821 | 0.019 | 0.26 |
| 0.5 | 2 | 9.2350e-4 | 1.00489 | 1.00491 | 0.0819 | 0.060 | 0.36 |
| 0.5 | 3 | 1.9769e-5 | 1.00271 | 1.00279 | 0.0526 | 0.053 | 0.21 |
| 0.5 | 4 | 3.9689e-7 | 1.00154 | 1.00171 | 0.0386 | 0.044 | 0.13 |
| 0.6 | 1 | 7.9640e-2 | 1.00111 | 1.00112 | 0.3329 | 0.003 | 0.04 |
| 0.6 | 2 | 4.3114e-3 | 1.01016 | 1.01017 | 0.1361 | 0.075 | 0.38 |
| 0.6 | 3 | 1.9992e-4 | 1.00636 | 1.00640 | 0.0858 | 0.075 | 0.24 |
| 0.6 | 4 | 8.6817e-6 | 1.00408 | 1.00417 | 0.0624 | 0.067 | 0.15 |
| 0.6 | 5 | 3.6395e-7 | 1.00269 | 1.00286 | 0.0489 | 0.058 | 0.11 |

The last two columns take the worse of the two resolutions, comparing `r`
and `r_H` at the same resolution.

## Held-out predictions

* All four angles are new and lie inside the in-sample range of `a`.
* `m` runs as high as the condition predicted `Δ/ω̄ ≥ 1e-8` allows. That
  threshold is the smallest split #355 resolved.
* The `E/H − 1` column is the fraction of the splitting the law attributes
  to the neck.

| a | m | predicted Δ = E | predicted Δ/ω̄ | H | E/H − 1 | bound on \|r − 1\| |
|---|---|---|---|---|---|---|
| 0.25 | 1 | 2.1249e-3 | 1.50e-3 | 2.0559e-3 | 0.0335 | 1.78e-3 |
| 0.25 | 2 | 3.0666e-6 | 1.25e-6 | 3.0173e-6 | 0.0163 | 1.78e-3 |
| 0.35 | 1 | 8.3890e-3 | 5.93e-3 | 7.8336e-3 | 0.0709 | 4.00e-3 |
| 0.35 | 2 | 4.7691e-5 | 1.95e-5 | 4.6145e-5 | 0.0335 | 4.00e-3 |
| 0.35 | 3 | 2.3234e-7 | 6.71e-8 | 2.2732e-7 | 0.0221 | 4.00e-3 |
| 0.45 | 1 | 2.3761e-2 | 1.68e-2 | 2.0984e-2 | 0.1323 | 9.20e-3 |
| 0.45 | 2 | 3.8202e-4 | 1.56e-4 | 3.6062e-4 | 0.0593 | 9.20e-3 |
| 0.45 | 3 | 5.2604e-6 | 1.52e-6 | 5.0643e-6 | 0.0387 | 9.20e-3 |
| 0.45 | 4 | 6.7961e-8 | 1.52e-8 | 6.6059e-8 | 0.0288 | 9.20e-3 |
| 0.55 | 1 | 5.5167e-2 | 3.90e-2 | 4.4471e-2 | 0.2405 | 1.93e-2 |
| 0.55 | 2 | 2.0630e-3 | 8.42e-4 | 1.8789e-3 | 0.0980 | 1.93e-2 |
| 0.55 | 3 | 6.6050e-5 | 1.91e-5 | 6.2144e-5 | 0.0629 | 1.93e-2 |
| 0.55 | 4 | 1.9820e-6 | 4.43e-7 | 1.8941e-6 | 0.0464 | 1.93e-2 |
| 0.55 | 5 | 5.7432e-8 | 1.05e-8 | 5.5391e-8 | 0.0368 | 1.93e-2 |

**Diagnostics, reported but not gated.**
* `a = 0.65` and `a = 0.7`, for `m = 1…5`. These are outside the in-sample
  range of `a`.

  | a | m = 1 | m = 2 | m = 3 | m = 4 | m = 5 |
  |---|---|---|---|---|---|
  | 0.65 | 1.1161e-1 | 8.5025e-3 | 5.5567e-4 | 3.3970e-5 | 2.0040e-6 |
  | 0.7 | 1.5231e-1 | 1.5918e-2 | 1.4317e-3 | 1.2029e-4 | 9.7475e-6 |

* The excited pairs `n = 1, 2` at every gated and diagnostic cell. They use
  the same construction with the band state `P_{m+n}^m(cos θ)`, its
  normalisation `2(ℓ+m)!/((2ℓ+1)(ℓ−m)!)`, and `ω² = ℓ(ℓ+1)`, where
  `ℓ = m + n`.
  * #362's WKB form failed near `ω = m/sin a` (`r` up to 2.8) and was not
    defined beyond it.
  * The exact neck needs no turning-point treatment at the rim, so `E` is
    evaluated at every cell.
  * At the larger apertures the predicted excited splittings are not small;
    `E` reaches 1.1 at (0.7, 1, n = 2). There the tunnelling expansion
    behind Herring's formula does not apply, which is one reason these
    pairs are not gated.
* The predicted reversal times `π/(2E)`.

## Measurement

The measurement is as in #362:
* Use #355's solver at `07d58de`, unchanged.
* For each geometry, take the lowest cosine-sector and sine-sector
  eigenvalues of `Sector(Loop(a, N), m, KLEIN)`, paired as in `p4_p5`.
* `Δ = |ω₊ − ω₋|`, at `N = 2000` and `N = 4000`.

**Reproduction check, before any held-out number is read.** The
implementation must reproduce all 17 in-sample `Δ` to a relative `1e-6` at
both sizes. These are #355's six and #362's eleven, against their own
archives. If it does not, the run stops and is reported `INCOMPLETE`.

## Gates

`r = Δ_measured/E` and `r_H = Δ_measured/H`, at each resolution.

* **G1, accuracy (tolerance).** At every gated cell,
  `|r − 1| ≤ 0.2 a⁴ + 1e-3`.
  * The `a⁴` term is about twice the in-sample envelope, `0.11 a⁴`.
  * The constant term is the solver's resolution floor.
  * The worst in-sample value of `|r − 1|/bound` is 0.65.
* **G2, the neck is the error (structural).** At every gated cell,
  `|r − 1| ≤ 0.2 |r_H − 1|`. That is, the exact neck removes at least 80% of
  `H`'s discrepancy. The worst in-sample ratio is 0.091.
* **G3, resolution (tolerance).** For every gated cell,
  `|Δ₂₀₀₀/Δ₄₀₀₀ − 1| ≤ 1e-3`.

## Verdicts

* `CONFIRMED` if G1, G2 and G3 hold at both resolutions.
* `COUNTEREXAMPLE` if G1 or G2 fails at both resolutions. The report names
  the gate and the cells.
* `UNRESOLVED_NUMERICALLY` if G3 fails, or if a G1 or G2 miss shrinks from
  `N = 2000` to `N = 4000`, that is, appears only at `N = 2000`. One further
  doubling to `N = 8000` is allowed for the affected cells, and nothing else
  changes.
* `INCOMPLETE` if the reproduction check fails.
* `NO_VERDICT` for any case not named here, including a miss that appears
  only at `N = 4000`.

Tolerances are not adjusted after the fact, and a miss is reported as a
miss. The report also gives, without gates:
* the signed `r − 1` at every cell;
* `max |r − 1|/a⁴`;
* the `a = 0.65` and `a = 0.7` diagnostics;
* the excited pairs;
* the predicted reversal times.

## Scope

* A classical scalar on a fixed 2D surface, as in #355 and #362. Not a GR
  solution and not on `S³`.
* The splitting is classical wave tunnelling through a geometric barrier.
  Nothing here is claimed to be quantum.
* A `CONFIRMED` verdict would show that, to the stated tolerance, the WKB
  neck step was the whole of `H`'s error. It would not derive the remaining
  residual, and it would not reopen #355's second gate.
