# Pre-registration: does the Alice handle's mode exchange keep the second reading complementary when a detector depletes the first?

Publish this file before any detector run on the handle.

* **Baseline.** `main` at `a68964f`.
* **Geometry and solver.** PR #355's `geometrodynamics/bulk/alice_handle_ports.py`, unchanged from `07d58de`.
* **Splitting law.** `docs/klein_exact_neck.md` (#362).
* **Randomness.** None is used.

## Question

#365 showed that zero doubles and zero misses at every power split require
complementary readings, u + v = 1, on every history. #367 obtained u + v = 1
from a variational selection on a Berger sphere. #367 and #368 then showed that
every detector load tested so far breaks that selection at first order in the
load.

The handle offers a different route. In the circulating basis, the lowest Klein
doublet in sector `m` is a two-level system. Write z for the +m amplitude and w
for the −m amplitude. The doublet evolves as `ω̄·1 + (Δ/2)σ_x`. Over a time
`π/Δ` it exchanges z and w up to phases. This is the angular-momentum reversal
of #355, and it is #365's map J realised by the geometry.

Suppose a port reads the +m intensity, u at one visit and v at a visit `π/Δ`
later. Then u + v = 1 for every state, from the linear dynamics alone, with no
selection principle. A detector that depletes only the channel it reads removes
energy from z. After the exchange, the second reading is the old w, which the
detector never touched.

**Hypothesis H.** On the actual handle, a local detector that absorbs only the
+m circulating channel keeps the second reading complementary to the first. The
remaining defect vanishes as the share of the exchange that happens during the
detection window, `Δτ`, goes to zero. A non-selective detector with the same
depletion breaks complementarity at first order.

In the two-mode truncation, H is close to true by construction. What is tested
is whether the full field obeys it, with:
* a spatially localised detector;
* a finite window during which the exchange is already running;
* leakage into the other radial modes;
* the frequency spread of the doublet;
* the work a moving medium can do on the field.

## Model

**Field and cells.**
* The field is the real scalar of #355 in azimuthal sector `m` with Klein
  gluing, discretised as `Sector(Loop(a, N), m, KLEIN)`.
* The cells are `(a, m)` = (0.5, 1), (0.4, 1), (0.3, 1) and (0.3, 2). Their
  predicted `Δτ` runs from 0.30 to 1.1e-4.
* Every cell is computed at N = 2000 and 4000.

**The doublet and the two channels.**
* The doublet is the lowest cosine-sector (A) mode and the lowest sine-sector
  (B) mode, identified as in `p4_p5`.
* Its mean frequency is `ω̄`, and `Δ = |ω_A − ω_B|` from the Rayleigh-quotient
  eigenvalues.
* The doublet phase space is 4-dimensional: the two modal coordinates and their
  velocities.
* On it, define the energy form E and the angular-momentum form L, where
  `L = π m ẋᵀ(diag(mass) ⊗ J)x` as in `Sector.angular_momentum`.
* The **+m channel** is the 2-dimensional eigenspace of `E⁻¹L` with positive
  eigenvalue, and the −m channel is the negative one. They are E-orthogonal.
* `E₊(s)` is the energy of the E-orthogonal projection of a doublet state s onto
  the +m channel.

**The detector.**
* It occupies a strip of the band around the equator. Its weight is
  `w(θ) = cos²(π(θ − π/2)/(2h))` for `|θ − π/2| < h`, with `h = 0.35`, and zero
  elsewhere.
* It is active on `0 < t < τ`, with `τ = 8`, at strength `γ(t) = γ₀ sin²(πt/τ)`.
* The force on band node i, with `(A, B)` in the φ frame, is
  `F_i = −γ(t) mass_i w(θ_i) [(Ȧ − Ω_d B), (Ḃ + Ω_d A)]`.

**Selective detector: a counter-rotating dashpot.**
* Here `Ω_d = ω̄`, the doublet mean frequency measured on the same mesh. This is
  a viscous medium rotating at angular velocity `−ω̄/m`.
* In its frame a −m wave at `ω̄` is at rest, so it feels no force. A +m wave
  sees frequency `2ω̄` and is absorbed.
* Away from `ω̄`, the medium can also do work on the field. This is the
  Zel'dovich effect, and it is counted in the energy balance.

**Non-selective control: a static dashpot.**
* `Ω_d = 0`, the same strip and the same window.

**Calibration.**
* For each configuration, `γ₀` is set by bisection so that the pure +m doublet
  state at t = 0 loses a fraction η of the field energy over the window.
* The selective detector is run at η = 0.1, 0.5 and 0.9, and the control at
  η = 0.5.

## Evolution and readings

* **Window.** The field is expanded in the lowest K modes of the sector, with
  K = 150 primary and 300 as the check. The damped linear equation is
  integrated over `0 ≤ t ≤ τ` with DOP853 (rtol 1e-11, atol 1e-13), together
  with the power absorbed by the detector.
* **After the window.** Each mode evolves exactly.
* **First reading.** `u = E₊(s_ref(τ/2))`: the +m doublet energy at the window
  centre on the undamped reference evolution. This is what an ideal
  instantaneous reading at the window centre would see.
* **Second reading.** `v = E₊(s(t₂))`: the +m doublet energy of the detected
  evolution at `t₂ = τ/2 + π/Δ`.
* **Normalisation.** Both are relative to the initial energy, which is 1.
  Initial states range over the whole doublet phase space.
* **Defect.** Both readings are quadratic forms in the initial state. The
  reported defect is the largest `|u + v − 1|` over all normalised doublet
  states, an eigenvalue problem on four basis evolutions.
* **Also reported:**
  * the leakage, meaning the energy outside the doublet at `t₂`;
  * the absorbed energy `D₁` against `η·u`;
  * the defect when `t₂` uses the exact-neck prediction for Δ instead of the
    measured one.

## Two-level expectation (not a gate)

`experiments/closure_ledger/alice_exchange_predictions.py` evaluates the
rotating-frame two-level model:
* `ż = i(Δ/2)w − (g/2)z` and `ẇ = i(Δ/2)z`;
* the same window and calibration;
* Δ from the exact-neck law.

It predicts the following largest defects:

| a | m | Δτ | η = 0.1 | η = 0.5 | η = 0.9 |
|---|---|---|---|---|---|
| 0.5 | 1 | 2.96e-1 | 8.1e-5 | 1.07e-3 | 6.85e-3 |
| 0.4 | 1 | 1.17e-1 | 1.6e-5 | 3.46e-4 | 2.45e-3 |
| 0.3 | 1 | 3.57e-2 | 3.1e-6 | 9.7e-5 | 7.2e-4 |
| 0.3 | 2 | 1.08e-4 | 7.5e-9 | 2.8e-7 | 2.2e-6 |

The two-level model omits several effects:
* the absorption, and the gain from the rotating medium, on the −m channel at
  the offsets `±Δ/2`;
* leakage into the other radial modes;
* the imperfect overlap of the A and B profiles.

The gates below are therefore set as margins above this expectation, while
keeping them well apart from the control.

## Gates

At both N = 2000 and 4000:

* **G1, protection.** For the selective detector at every cell and every η,
  the largest defect is at most `0.05·η`.
* **G2, the limit.** At (0.3, 2), where `Δτ ≈ 1e-4`, the selective detector's
  largest defect is at most `1e-4` for every η.
* **G3, the control breaks it.** For the static dashpot at η = 0.5, the largest
  defect is at least `0.25` at every cell.
* **G4, numerics.**
  * With no detector, the defect is at most `1e-8`.
  * The K = 150 and K = 300 defects agree to within `max(1e-5, 10%)`.
  * The N = 2000 and N = 4000 selective defects agree to within
    `max(1e-4, 20%)`.
  * The absorbed work balances the change in field energy to within 1e-8 of
    the initial energy.

## Verdicts

* `CONFIRMED` if G1–G4 hold.
* `COUNTEREXAMPLE` if G1, G2 or G3 fails at both N. The report names the gate
  and the cell.
* `UNRESOLVED_NUMERICALLY` if G4 fails, or if a G1–G3 miss appears at one N
  only. One further refinement, K = 600, is allowed, and nothing else changes.
* `NO_VERDICT` for any case not named here.

Tolerances are not adjusted after the fact, and a miss is reported as a miss.

## Scope

* A classical scalar on #355's fixed 2D surface, with a detector put in by hand.
  The detector's selectivity, a counter-rotating medium tuned to the doublet,
  is an input and is not derived from the geometry.
* The exchange is the handle's own, through the Klein gluing and tunnelling.
* The test concerns complementarity of the readings under depletion. It does
  not derive a probability law: the uniform measure on doublet states that
  would give Born weights remains an input.
* It does not model detector durability, a threshold record or a GR solution.
  Nothing here is claimed to be quantum.
