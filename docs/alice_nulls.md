# The residual at the spectral nulls of the Alice-handle exchange: result

This tests `docs/alice_nulls_prereg.md`, frozen at `f03b034`.

* **Archive.** `experiments/closure_ledger/runs/20261126_alice_nulls/`, which
  holds `probe.json`, `probe.md`, `frozen_predictions.json`, and the in-sample
  exploration with its scripts.
* **Code.** `alice_nulls_run.py`, committed at `be69c01` before the run, and
  `alice_nulls_predictions.py`, which is part of the freeze.

## Verdict

**`CONFIRMED`.** All six gates pass on the held-out cells (0.45, 2) and
(0.55, 2). No primary point lies outside the domain. The largest peak drain
rate is 0.462ω̄, at k = 5 and η = 0.9.

| cell | η | τ = 6.41 | 10.26 | 16.67 | 25.65 |
|---|---|---|---|---|---|
| (0.45, 2) | 0.5 | 1.92e-4 (1.06) | 8.98e-5 (1.04) | 6.67e-5 (0.95) | 6.46e-5 (0.85) |
| (0.45, 2) | 0.9 | 1.18e-3 (1.04) | 4.07e-4 (1.02) | 2.46e-4 (0.93) | 2.56e-4 (0.86) |
| (0.55, 2) | 0.5 | 4.47e-4 (1.13) | 3.49e-4 (1.07) | 3.34e-4 (0.96) | 3.51e-4 (0.85) |
| (0.55, 2) | 0.9 | 1.84e-3 (1.11) | 1.12e-3 (1.05) | 1.09e-3 (0.93) | 1.35e-3 (0.86) |

Each entry is the full-field defect exactly at a sin² null, with the ratio of
measured to predicted in brackets.

* **N1, null values.** Measured ÷ predicted is 0.85–1.13 at all 16 points,
  against ×2.
* **N2, no radial-doublet contribution.** The doublet alone agrees with the full
  150-mode field to within 1.02%, against 3%.
* **N3, mismatch.**
  * The strip mismatch is `|ε|/(Δ/ω̄)` = 3.93 and 3.96. The in-sample value was
    3.9, and the gate was [3, 5].
  * The mismatch share, doublet minus identical-profile, is within 0.84–1.10 of
    `K·Δ/ω`, against ×2.
* **N4, the intrinsic law.** Identical-profile defect minus exchange, times
  `(ωτ)³`, is within 0.99–1.06 of C(η), against ×2.
* **N5, the depletion boundary.** The largest reachable depletion matches the
  in-sample values to within 0.0009, against ±0.005, at both cells and all three
  windows.
* **N6, numerics.**
  * K = 80 and K = 150 agree to 2e-7 relative, and N = 2000 and N = 4000 to
    6e-6.
  * The energy balance holds to 2e-11, and the undetected defect is 4e-11.

## What it means

The defect of the handle's exchange under a channel-selective, depleting
detector now has a complete description in this model. Every term was predicted
before the run on cells never measured:

    E = A(η)|F̂(2ω; τ)| + c_x(η) E_2L + C(η)/(ωτ)³ + K(η) Δ/ω

* **Switching, first order.** The counter-rotating coupling at 2ω̄, weighted by
  the window's spectrum. It sets the V-curve's left branch and vanishes at the
  sin² nulls.
* **Exchange during the window.** #370's rotating-frame two-level defect. It
  sets the right branch.
* **Intrinsic, second order in the drain rate.** It falls as `(ω̄τ)⁻³` with no
  spectral nulls. Measured in isolation, it matches the law to within 6%.
* **Mismatch.** The detector strip sees the two doublet modes differ by
  `ε ≈ 3.9 Δ/ω̄`, a tunnelling correction to the mode profiles. So it absorbs a
  little of the channel it should not. The term is roughly independent of τ and
  proportional to Δ.

The other radial doublets play no part: the doublet alone reproduces the full
field to 1%. Three of the four terms vanish as Δ → 0 or τ → ∞, and the fourth,
switching, vanishes at the nulls. So on this description the defect can be made
as small as wanted, inside the domain, by a small aperture with a long window
placed at a null.

**The domain is a property of the detector, not the handle.**
* The largest reachable depletion is the same at four apertures, to within
  0.001.
* It saturates where the medium's peak drain rate reaches about 1.3–1.9ω̄, and
  the medium then drags the field instead of absorbing it.
* The frozen validity rule, a peak drain rate of at most 0.5ω̄, kept every
  primary point inside.

## One systematic trend (not gated)

At both cells and both depletions, measured ÷ predicted falls from about 1.1 at
τ = 6.4 to 0.85 at τ = 25.7. The long windows are where the exchange term
dominates. So the four terms adding exactly, or the in-sample `c_x`, slightly
overestimate there. The trend is well inside the ×2 tolerance, and is recorded
for a future refinement.

## Scope

* A classical scalar on #355's fixed surface. The exchange is the handle's own.
  The detector's channel selectivity, a counter-rotating medium tuned to the
  doublet, is put in by hand.
* The coefficients A, c_x, C and K are empirical, fitted in-sample.
  * The intrinsic term's `(ωτ)⁻³` law is measured, not derived.
  * Its mechanism is presumably a second-order, non-rotating-wave process. It is
    not identified here.
* No probability law is derived: the measure on doublet states that would give
  Born weights remains an input. Nothing here is claimed to be quantum.

## Addendum (after the verdict): a closed form for the intrinsic term

Proposed in review after the `CONFIRMED` run, then checked with
`exploration/closed_form_check.py` (archived as `closed_form_check.json`).

**Origin.** The intrinsic term is not a separate second-order process. It is the
same first-order counter-rotating drive at 2ω̄ that gives the switching term,
`ȧ₋ = ½γ(t) a₊*(t) e^{2iω̄t}`, acting through a **depleting** +m amplitude.

**The calculation.**
* Write `g(t) = ½γ(t) e^{−½∫γ}`, with γ the sin² drain and `Γ = −ln(1 − η)`.
* Integrate the drive by parts at an exact sin² null, where `e^{2iω̄τ} = 1`.
* The sin² boundary terms no longer cancel, because `a₊(τ) = √(1 − η) ≠ a₊(0)`.
* Summing the even-derivative series gives the factor `(1 − k⁻²)⁻¹`.
* At balanced preparation the defect is `C_int/(ω̄τ)³`, with

      C_int(η, k) = (π²/4) (−ln(1 − η)) (1 − √(1 − η)) / (1 − k⁻²)

* At small η this is about `(π²/8)Γ²`, which is why the term looked second order
  in the drain rate.

**The check.** It has no free parameter. It was evaluated in the doublet-only
model with identical strip profiles, so there is no mismatch term, and with the
splitting cut 10⁴-fold, so the exchange part is negligible. It covers
(0.3, 2), (0.45, 2) and (0.55, 2), η = 0.1, 0.25, 0.5, 0.75 and 0.9, and
k = 5, 8, 13 and 20.
* Measured ÷ closed form is **0.9987–1.0215**, and within 0.2% wherever k ≤ 13.
* At η = 0.75 the closed form gives 1.7815, 1.7374, 1.7205 and 1.7146 for
  k = 5, 8, 13 and 20. The measured values agree to within 0.03% at k ≤ 13.
* The η = 0.25 and 0.75 points never entered any fit. The frozen in-sample
  coefficients, C = 0.51 and 3.95, are this formula's values at k ≈ 7–8.

**What it changes.** The empirical C(η) of the predictor is replaced by a
leading-order closed form in the counter-rotating drive, neglecting Δ. It has
no free parameter, and this check confirms it to 0.1–0.2%.

**The long-window drift.** Measured ÷ predicted falls from about 1.1 to 0.85 as
τ grows.
* It is not depletion during the window: the predictor's E_2L already depletes
  the +m amplitude.
* After the closed-form intrinsic term is subtracted at (0.5, 2), the exchange
  part needs a coefficient of 1.10 and 1.01 at k = 20, against the frozen
  c_x = 1.3 and 1.1.
* So c_x was inflated. When the two-term V-curve predictor was fitted, it
  absorbed the then-unmodelled terms, and it overshoots where exchange
  dominates.
* At short windows the inferred coefficient rises to as much as 2.1, because
  the terms add coherently rather than as magnitudes.
* A refit with the closed-form intrinsic term and coherent addition is the
  natural refinement. It is not done here.
