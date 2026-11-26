# Pre-registration: the residual at the spectral nulls of the Alice-handle exchange, and the detector's domain of validity

Publish this file before any detector run on the held-out cells (0.45, 2) and (0.55, 2).

* **Baseline.** Branch `claude/hopeful-babbage-ps703o` at `fd4e2f1`, which is
  PR #370.
* **Model.** That of `docs/alice_vcurve_prereg.md` (freeze `291ab5f`).
* **Solver.** #355's solver, unchanged from `07d58de`.

## Prior state

The V-curve test, `docs/alice_vcurve.md`, came out `COUNTEREXAMPLE` through G4.
The sin² nulls of the switching term fell where predicted, but they were
shallower. A residual of 1.8–2.8e-4 remained at (0.5, 2), which the two-part
predictor did not contain. The candidate offered there was the counter-rotating
coupling through the other radial doublets, at `ω_n + ω̄`.

## How this hypothesis was found (in-sample only)

All exploration used the cells already measured: (0.5, 2), (0.35, 2) and
(0.3, 2). It is archived in
`experiments/closure_ledger/runs/20261126_alice_nulls/exploration/`, with the
scripts that produced it. In order:

1. **Radial truncation, which refutes the candidate.**
   * Setup: at (0.5, 2), at the four null minima and three peaks of the archived
     scan, the detector was coupled to the doublet alone, then to the doublet plus
     the ℓ ≈ 4, 6, 8, 10 and 12 doublets, then to the full 150-mode basis with
     and without ℓ ≈ 4 or ℓ ≈ 6.
   * The equator-symmetric strip couples only to equator-even doublets; its
     couplings to the odd ones are below 1e-10.
   * Result: every variant gives the same defect to three digits, for example
     1.90e-4 at τ = 9.0.
   * **The other radial doublets own none of the residual. It lives inside the
     doublet.**
2. **Dissecting the doublet.** In the doublet-only model I changed two things
   separately: making the two mode profiles identical inside the strip, and
   cutting Δ tenfold.
   * At (0.5, 2) the profiles differ by `ε = w_AA/w_BB − 1 = 1.5e-3`. Making
     them identical removes about half of the residual.
   * At (0.3, 2), where ε = 2.1e-5, neither change matters. A residual that
     does not depend on the mismatch remains.
3. **The intrinsic part, measured exactly at the nulls `τ_k = kπ/ω̄` with
   identical profiles at (0.3, 2).**
   * It is close to quadratic in `ln(1/(1 − η))`, i.e. second order in the
     drain rate, and it falls as `(ω̄τ)⁻³`.
   * `(I − exchange)·(ω̄τ)³` is 0.50–0.52 at η = 0.5 and 3.88–4.04 at η = 0.9,
     over k = 5 to 13.
4. **The mismatch part,** actual minus identical-profile doublet at (0.5, 2).
   * It is about 1.0e-4 at η = 0.5 and 1.0–2.2e-4 at η = 0.9, roughly
     independent of τ.
   * The mismatch itself is set by the splitting: `ε/(Δ/ω̄)` = 3.96 at
     (0.5, 2) and 3.88 at (0.3, 2). So the mismatch term is written as
     `K(η)·Δ/ω̄`, with Δ from the exact-neck law.
5. **Combining the terms.** Adding the four terms reproduces every in-sample
   null value to 0.85–1.09×. That covers 8 points at (0.5, 2) and 8 at (0.3, 2).
   * The fitted coefficients are C = 0.51 and 3.95, and K = 0.26 and 0.42, for
     η = 0.5 and 0.9.
6. **The depletion boundary.** At τ = 2 to 9 with sin² and sin⁴, the largest
   reachable depletion is the same at (0.5, 2) and (0.35, 2) to three digits.
   For sin² it is 0.720, 0.797, 0.877 and 0.943 at τ = 2, 2.5, 3.21 and 4.
   * So the boundary is a property of the detector–wave coupling, a function of
     ω̄τ and the envelope, not of the handle.
   * The depletion stops rising when the peak drain rate is about 1.3–1.9ω̄:
     the medium then drags the field instead of absorbing it.
   * Reaching η = 0.9 needs a peak drain rate of 0.93ω̄ at τ = 4, 0.50ω̄ at
     τ = 6 and 0.42ω̄ at τ ≥ 7.

## The predictor (frozen in `alice_nulls_predictions.py`, which calls no solver)

    E = A(η)|F̂(2ω; τ)| + c_x(η) E_2L + C(η)/(ωτ)³ + K(η) Δ/ω

* `ω = √(m(m+1))`, and Δ is the exact-neck splitting.
* A and c_x are unchanged from `291ab5f`.
* C = 0.51 and 3.95, and K = 0.26 and 0.42, for η = 0.5 and 0.9.

At an exact sin² null the switching term vanishes. What remains is the
exchange, the intrinsic term and the mismatch term.

## Domain of validity

The detector is a dashpot absorbing a wave of frequency ω̄. When its peak drain
rate approaches ω̄, it stops absorbing cleanly and drags the field. The rule
frozen here is:

> **A point lies in the domain only if its calibrated peak drain rate is at most `0.5 ω̄`.**

* Points outside the domain are reported but not scored.
* All primary points below are predicted to lie inside. The tightest is k = 5
  at η = 0.9, which needs about 0.46ω̄.
* If any primary point falls outside, the verdict is `NO_VERDICT`, and the
  design gap is reported.

## Held-out design

* **Cells.** (0.45, 2) and (0.55, 2), with m = 2. Their exact-neck Δ is 3.82e-4
  and 2.06e-3. Neither has been run with a detector.
* **Primary points.**
  * Exact nulls `τ_k = kπ/ω̄_FEM` for k = 5, 8, 13 and 20.
  * Envelope sin², at η = 0.5 and 0.9.
  * Each is evaluated with the full field (N = 2000, K = 150), with the
    doublet alone, and with the doublet with identical strip profiles.
* **Boundary map.** At each cell, sin², τ = 2.5, 3.21 and 4.0, the largest
  depletion reachable on the weak-coupling branch.
* **Numerics.**
  * Every primary point is repeated at K = 80.
  * The k = 8 points are repeated at N = 4000.
  * The energy balance and the undetected defect are also reported.

**Frozen predictions** (`frozen_predictions.json`). The table gives the total
predicted defect at each null.

| cell | η | τ = 6.413 | 10.260 | 16.673 | 25.651 |
|---|---|---|---|---|---|
| (0.45, 2) | 0.5 | 1.81e-4 | 8.61e-5 | 6.99e-5 | 7.63e-5 |
| (0.45, 2) | 0.9 | 1.14e-3 | 4.00e-4 | 2.64e-4 | 2.97e-4 |
| (0.55, 2) | 0.5 | 3.96e-4 | 3.25e-4 | 3.48e-4 | 4.12e-4 |
| (0.55, 2) | 0.9 | 1.67e-3 | 1.07e-3 | 1.18e-3 | 1.56e-3 |

* **Mismatch term.** 4.05e-5 and 6.55e-5 at (0.45, 2), and 2.19e-4 and 3.54e-4
  at (0.55, 2), for η = 0.5 and 0.9.
* **Intrinsic term.** The same at both cells: 1.32e-4 at τ = 6.413 for η = 0.5,
  and 1.02e-3 for η = 0.9.

## Gates

* **N1, null values.** At all 16 primary points, the full-field defect is
  within `[1/2, 2]` of the predicted total.
* **N2, no radial-doublet contribution.** At all 16 points, the doublet-only
  defect agrees with the full field to within 3%.
* **N3, mismatch attribution.**
  * At both cells, ε/(Δ_FEM/ω̄) lies in [3, 5].
  * At every η = 0.5 point, doublet-only minus identical-profile is within
    `[1/2, 2]` of `K(0.5)·Δ/ω`.
* **N4, the intrinsic law.**
  * Take the identical-profile doublet at (0.45, 2), with k = 5 and 8, at both
    η, and subtract the predicted exchange term.
  * The result times `(ωτ)³` must be within `[1/2, 2]` of C(η).
* **N5, the depletion boundary.** At both cells and the three boundary windows,
  the largest reachable depletion equals the in-sample value to within ±0.005.
  The in-sample values are 0.797, 0.877 and 0.943.
* **N6, numerics.**
  * K = 80 and K = 150 agree within `max(1e-6, 5%)`.
  * N = 2000 and N = 4000 agree within `max(1e-6, 10%)`.
  * The energy balance holds to 1e-8, and the undetected defect is at most
    1e-8.

## Verdicts

* `NO_VERDICT` if any primary point falls outside the domain.
* Otherwise:
  * `CONFIRMED` if N1–N6 hold;
  * `COUNTEREXAMPLE` if any of N1–N5 fails, naming the gate and the point;
  * `UNRESOLVED_NUMERICALLY` if only N6 fails;
  * `NO_VERDICT` for any case not named here.

Tolerances are not adjusted after the fact.

## Scope

* As before, a classical scalar on #355's fixed surface, with a detector whose
  selectivity is put in by hand.
* The intrinsic term's `(ωτ)⁻³` law and its coefficient are empirical. A
  candidate mechanism is a second-order, non-rotating-wave process. It is not
  derived here.
* No probability law is derived, and nothing here is claimed to be quantum.
