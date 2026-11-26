# V-curve, aperture scaling and envelope control for the Alice-handle exchange: result

This tests `docs/alice_vcurve_prereg.md`, frozen at `291ab5f`.

* **Archive.** `experiments/closure_ledger/runs/20261125_alice_vcurve/`, which
  holds:
  * `probe.json` and `probe.md`;
  * `frozen_predictions.json`;
  * the in-sample exploration in `exploration/`;
  * `posthoc_null_depth.json`.
* **Code.** `alice_vcurve_run.py`, committed at `a2a9fb4` before the run;
  `alice_vcurve_probe.py`; and `alice_vcurve_predictions.py`, which is part of
  the freeze.

## Verdict

**`COUNTEREXAMPLE`, through G4 alone: the spectral nulls are where predicted but
shallower than frozen.** The other five gates pass.

The table gives the largest `|u + v − 1|` for sin² at N = 2000 and K = 150, with
the frozen prediction in brackets.

| cell | η | τ = 3.21 | 7.05 | 14.75 | 31.42 | 63.49 | 127.61 | 255.87 |
|---|---|---|---|---|---|---|---|---|
| (0.5, 2) | 0.5 | 8.40e-3 (8.0e-3) | 7.25e-4 (6.8e-4) | 2.05e-4 (1.2e-4) | **1.66e-4** (1.1e-4) | 2.18e-4 (2.1e-4) | 3.82e-4 (4.6e-4) | 8.20e-4 (1.0e-3) |
| (0.5, 2) | 0.9 | unreachable | 1.85e-3 (1.8e-3) | **5.94e-4** (4.8e-4) | 6.95e-4 (6.6e-4) | 1.25e-3 (1.3e-3) | 2.52e-3 (2.7e-3) | 5.34e-3 (5.8e-3) |
| (0.35, 2) | 0.5 | 8.26e-3 (8.0e-3) | 5.93e-4 (6.5e-4) | 6.96e-5 (7.2e-5) | 1.45e-5 (1.2e-5) | **1.14e-5** (1.1e-5) | 1.76e-5 (2.1e-5) | 3.32e-5 (4.2e-5) |
| (0.35, 2) | 0.9 | unreachable | 1.51e-3 (1.7e-3) | 1.82e-4 (2.0e-4) | **5.09e-5** (5.2e-5) | 6.49e-5 (6.9e-5) | 1.23e-4 (1.3e-4) | 2.45e-4 (2.7e-4) |

* **G1, the predictor, passes.** At every measurable point, measured ÷ predicted
  is 0.79–1.76, against the frozen [1/3, 3].
* **G2, the V shape, passes.** All four minima are interior and fall exactly at
  the predicted grid point.
* **G3, aperture scaling, passes.** Going from (0.5, 2) to (0.35, 2) lowers the
  minimum 14.6× at η = 0.5, to 1.14e-5, and 11.7× at η = 0.9, to 5.09e-5. Both
  are below 1e-4.
* **G4, spectral nulls, fails.** All four nulls are where predicted: the minima
  fall at 6.4, 7.7, 9.0 and 10.3, against 6.413, 7.695, 8.978 and 10.260.
  Measured against the nearby peak, they are 0.29, 0.31, 0.40 and **0.50**.
  The last exceeds the frozen 0.4.
* **G5, envelope control, passes.**
  * sin⁴ ÷ sin² is 0.136 at τ = 8 and 0.151 at τ = 16, against predicted 0.117
    and 0.061. At τ = 16 that is 2.49× the prediction, just inside the
    allowed 2.5.
  * bump ÷ sin² is 1.76 and 3.40, against predicted 1.22 and 3.63.
  * At τ = 16 the C^∞ bump is worse than sin², as predicted.
* **G6, numerics, passes.**
  * K = 80 and K = 150 agree to 6e-6 relative at the minima, and N = 2000 and
    N = 4000 to 8e-5 relative.
  * The energy balance holds to 5.1e-12. With no detector, the defect is
    1.7e-11.

**Unreachable points.** At τ = 3.21 the selective detector cannot remove 90% of
the +m channel on the weak-coupling branch. It reaches at most 0.873 at a₁ and
0.874 at a₂. A stiffer medium drags the field instead of absorbing it, as #370
found. The prereg did not anticipate this.
* The runner first stopped at that point.
* It was fixed in `26927be` to record the largest reachable depletion rather
  than crash. No gate or tolerance changed.
* Had nothing else failed, these two points would have made the verdict
  `NO_VERDICT`.

## What it means

The floor found in #370 is explained, and its two-sided structure was predicted
before the held-out run.

* **Switching.** On the left of the V, the defect is set by the window's
  spectral weight at **2ω̄**. That is the counter-rotating coupling between the
  +m channel and the conjugate of the −m channel, which the rotating-frame
  model drops.
  * Its frozen coefficient predicts held-out points to within 0.79–1.76.
  * It also predicts where the sin² spectral nulls fall in τ, and the ranking of
    envelopes.
  * Smoother is not better as such. sin⁴, which is C³, beats sin² because its
    weight at 2ω̄ is smaller. The C^∞ bump loses, because at the same support
    it is narrower and carries more weight at 2ω̄.
* **Exchange.** On the right of the V, the defect is the exchange running
  during the window. It follows #370's two-level model, scaled by about 1.1–1.3.
* **Aperture.** Δ falls exponentially with the aperture while ω̄ stays fixed,
  so the best window lengthens and the best defect falls.
  * At (0.35, 2), a channel-selective detector removing 50% or 90% of the
    channel it reads leaves the complementary reading intact to 1.1e-5 and
    5.1e-5.
  * By the predictor, smaller apertures go lower still, until the residual
    below sets the limit.

## What failed (post hoc, `posthoc_null_depth.json`)

The predictor expected nulls of depth 0.05–0.19. The measured nulls are
shallower because a residual of 1.8–2.8e-4 remains at each null.

Neither predicted term accounts for that residual: the exchange part there is
only 2–3e-5. The in-sample cell (0.3, 2) showed a smaller residual of about
5e-5, which the prereg noted. The frozen depth threshold of 0.4 was set from
in-sample depths of 0.15–0.23. At the larger held-out aperture the residual is
about 4× larger, and the last null failed.

So the predictor is incomplete. A third contribution, which depends on the
aperture, sets how deep the nulls go. Its origin is not identified. Candidates
include:
* coupling through the other radial doublets, whose own counter-rotating
  frequencies `ω_n + ω̄` the window does not null;
* the cross term between switching and exchange.

Identifying it needs its own test.

## Procedure notes

* **The design came from in-sample exploration.**
  * It used #370's cells, (0.3, 1) and (0.3, 2), only.
  * The proposed C¹-against-C^∞ envelope control would have failed for the reason
    above, so the frozen version predicts the ranking instead.
  * All exploration data are archived.
* **The runner was fixed after it stopped.**
  * The held-out run stopped at the first unreachable point.
  * The runner was fixed and the run resumed from cached points, with no gate
    changed.
  * The archived `seconds` covers only the resumed run.
* **The PR's branch.** This work is on #370's branch, so it lands in that PR.

## Scope

* As in #370, a classical scalar on #355's fixed surface. The detector's
  selectivity is put in by hand.
* The predictor's coefficients are empirical, and the null residual is
  unexplained.
* No probability law is derived, and nothing here is claimed to be quantum.
