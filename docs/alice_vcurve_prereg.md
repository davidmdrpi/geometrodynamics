# Pre-registration: the window-length V-curve, aperture scaling and envelope control for the Alice-handle exchange

Publish this file before any detector run on the held-out cells (0.5, 2) and (0.35, 2).

* **Baseline.** Branch `claude/hopeful-babbage-ps703o` at `d767a19`, which is
  PR #370 on `main`.
* **Model.** The model, detector, readings and defect of
  `docs/alice_exchange_prereg.md` (#370), with the window length τ and the
  switching envelope as parameters. They are implemented in
  `experiments/closure_ledger/alice_vcurve_probe.py`.
* **Solver.** #355's solver, unchanged from `07d58de`.

## Prior state

#370 found that a counter-rotating dashpot, which absorbs only the +m channel
it reads, keeps the second reading complementary. Its defect, the largest
`|u + v − 1|` over all doublet states, is 5e-5 to 2e-2.

It also found a floor. At (0.3, 2) the defect was 2.8e-4 at η = 0.5, against the
frozen 1e-4. A post-hoc diagnostic showed the floor falling about 400-fold as the
window length went from 4 to 32.

## How this hypothesis was found (in-sample only)

All exploration used the cells already run in #370, (0.3, 1) and (0.3, 2), at
N = 2000 and K = 150. Its data are archived with the result. In order:

1. **sin², sin⁴ and a C^∞ bump at τ = 4 to 64, (0.3, 2), η = 0.5.**
   * sin⁴ beat sin² by 8× at τ = 8 and 15× at τ = 16.
   * The bump was **worse** than sin² at every τ up to 32.
   * So the smoothness of the edges is not the controlling variable.
2. **Which frequency matters.** I compared the measured defect with each
   window's area-normalised Fourier weight at two candidate offsets: the
   doublet's gaps to the strongly coupled modes, and 2ω̄.
   * The weight at **2ω̄** tracked the defect for all three envelopes, with
     defect ÷ |F̂| between 0.30 and 0.44 wherever switching dominates.
   * 2ω̄ is the counter-rotating coupling: a detector tuned to the +m channel at
     ω̄ couples it to the conjugate of the −m channel at 2ω̄. The rotating-frame
     two-level model of #370 drops that term.
3. **A fine τ scan on (0.3, 2).** The sin² weight vanishes at
   `2ω̄τ = 2πk`, i.e. `τ_k = kπ/ω̄`.
   * Between τ = 6 and 11 on a 0.2 grid, the measured minima fell at 6.4,
     7.6–7.8, 9.0 and 10.2, against the predicted 6.41, 7.70, 8.98 and 10.26.
   * Each minimum was 0.15–0.23 of the nearby peak. A residual of about 5e-5
     remains at the nulls.
4. **Long windows on (0.3, 1), and η = 0.9 on both cells.**
   * The defect rises again on the right, following the two-level exchange
     defect.
   * At η = 0.5 the V on (0.3, 1) is 2.7e-3, 8.2e-4, 8.1e-4, 1.2e-3, 2.9e-3 and
     9.7e-3 at τ = 8, 16, 32, 64, 128 and 256.
5. **A two-part fit.** `E = A(η)|F̂(2ω̄)| + c_x(η)·E_two-level` was fitted to the
   53 points at η = 0.5 and the 10 at η = 0.9.
   * Frozen as A = 0.33 and 0.86, and c_x = 1.3 and 1.1, for η = 0.5 and 0.9.
   * On the sin² doubling grids, measured ÷ predicted is 0.77–1.91 at η = 0.5
     and 0.89–1.48 at η = 0.9.
   * It misses by up to 3–5× at a spectral null, at sin⁴ with τ = 4, and in the
     crossover where both parts are comparable.

The predictor, the held-out design and the predictions were all fixed after
this exploration and before any held-out run.

## The predictor (frozen in `alice_vcurve_predictions.py`, which calls no solver)

    E_pred(τ, envelope, η) = A(η) |F̂_env(2ω; τ)| + c_x(η) E_2L(Δ, τ, envelope, η)

* `ω = √(m(m+1))`.
* `F̂` is the envelope's Fourier transform over `[0, τ]`, divided by its area.
* `E_2L` is #370's rotating-frame two-level defect, with the same envelope,
  τ, depletion and calibration.
* Δ is the exact-neck splitting.

It predicts the shape of the V:
* **Left branch.** For sin², the C¹ envelope, `|F̂|` decays as τ⁻³ with
  oscillations, and it vanishes at `τ_k = kπ/ω`.
* **Right branch.** The two-level exchange defect is roughly linear in Δτ,
  steepening as Δτ approaches 1.

## Held-out design

* **Cells.**
  * `a₁ = (0.5, 2)`, with Δ_E = 9.24e-4.
  * `a₂ = (0.35, 2)`, with Δ_E = 4.77e-5.
  * Both have m = 2, so ω is shared and Δ falls about 19× between them.
  * Neither has been run with a detector. (0.5, 2) is an in-sample cell of the
    exact-neck law only.
* **V-curve grid.** τ_j = (k + ½)π/ω for k = 2, 5, 11, 24, 49, 99 and 199, i.e.
  3.21, 7.05, 14.75, 31.42, 63.49, 127.61 and 255.87. These are sin² sidelobe
  peaks, so no grid point sits on a spectral null.
  * Envelope sin², at η = 0.5 and 0.9, on both cells.
* **Null scan.** `a₁`, sin², η = 0.5, τ from 6.0 to 11.0 in steps of 0.1.
* **Envelope control.** `a₂`, η = 0.5, τ = 8 and 16, with envelopes sin², sin⁴
  and bump.
* **Numerics.** N = 2000 and K = 150 throughout.
  * At each V minimum, the run is repeated with K = 80 and with N = 4000.
  * Each run also reports the undetected defect and the energy balance.

**Frozen predictions,** from `alice_vcurve_predictions.py`:

| cell | η | 3.21 | 7.05 | 14.75 | 31.42 | 63.49 | 127.61 | 255.87 | predicted minimum |
|---|---|---|---|---|---|---|---|---|---|
| a₁ | 0.5 | 8.0e-3 | 6.8e-4 | 1.17e-4 | **1.09e-4** | 2.1e-4 | 4.6e-4 | 1.0e-3 | 1.09e-4 at 31.42 |
| a₁ | 0.9 | 2.1e-2 | 1.8e-3 | **4.8e-4** | 6.6e-4 | 1.3e-3 | 2.7e-3 | 5.8e-3 | 4.8e-4 at 14.75 |
| a₂ | 0.5 | 8.0e-3 | 6.5e-4 | 7.2e-5 | 1.23e-5 | **1.12e-5** | 2.1e-5 | 4.2e-5 | 1.12e-5 at 63.49 |
| a₂ | 0.9 | 2.1e-2 | 1.7e-3 | 2.0e-4 | **5.2e-5** | 6.9e-5 | 1.3e-4 | 2.7e-4 | 5.2e-5 at 31.42 |

* **Envelope control at a₂, η = 0.5.**

  | τ | sin² | sin⁴ | bump |
  |---|---|---|---|
  | 8 | 3.0e-4 | 3.6e-5 | 3.7e-4 |
  | 16 | 5.7e-5 | 3.5e-6 | 2.1e-4 |

  The predicted ratios to sin² are 0.117 and 1.22 at τ = 8, and 0.061 and 3.64
  at τ = 16.
* **Predicted nulls at a₁.** τ = 6.413, 7.695, 8.978 and 10.260.

## Gates

* **G1, the predictor.** At every V-grid point, measured ÷ predicted lies in
  [1/3, 3].
* **G2, V shape.** For each cell and η:
  * the measured minimum over the grid is interior, neither the first nor the
    last point;
  * it lies at the predicted minimum's grid point or an adjacent one.
* **G3, aperture scaling.** For both η:
  * `E_min(a₂) ≤ E_min(a₁)/3`;
  * `E_min(a₂) ≤ 1e-4`.
* **G4, spectral nulls.** For each predicted null `τ_k` at a₁:
  * the scan's lowest value within ±0.64 of `τ_k` lies within 0.15 of `τ_k`;
  * that value is at most 0.4 times the scan's highest value in the same
    interval.
* **G5, envelope control.**
  * At each τ, the measured sin⁴/sin² and bump/sin² ratios are within a factor
    of 2.5 of the predicted ratios.
  * At τ = 16, sin⁴ ≤ 0.2 × sin², and bump > sin², so the smoother envelope is
    worse.
* **G6, numerics.**
  * The undetected defect is at most 1e-8, and the energy balance holds to 1e-8.
  * At the V minima, K = 80 against K = 150 agree within `max(1e-6, 10%)`, and
    N = 2000 against N = 4000 within `max(1e-6, 20%)`.

## Verdicts

* `CONFIRMED` if G1–G6 hold.
* `COUNTEREXAMPLE` if any of G1–G5 fails. The report names each failing gate and
  point.
* `UNRESOLVED_NUMERICALLY` if only G6 fails.
* `NO_VERDICT` for any case not named here.

Tolerances are not adjusted after the fact, and every point is reported.

## Scope

* As in #370, this is a classical scalar on #355's fixed surface, with a detector
  whose selectivity is put in by hand.
* The predictor's coefficients are empirical, fitted in-sample. The 2ω̄
  mechanism is identified, not derived from first principles.
* No probability law is derived, and nothing here is claimed to be quantum.
