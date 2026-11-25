# The Alice handle's mode exchange under detector depletion: result

This tests `docs/alice_exchange_prereg.md`, frozen at `0fb3b6e`.

* **Archive.** `experiments/closure_ledger/runs/20261125_alice_exchange/`
  (`probe.json`, `probe.md`, and `posthoc.json` for the diagnostics below).
* **Code.**
  * `experiments/closure_ledger/alice_exchange_probe.py`, committed at
    `811aaa8` before the run;
  * `alice_exchange_predictions.py`, which is part of the freeze;
  * `alice_exchange_posthoc.py`, written after the verdict.
* **Solver.** PR #355's `alice_handle_ports.py`, unchanged from `07d58de`.

## Verdict

**`COUNTEREXAMPLE`, through G2.** G1, G3 and G4 pass at both N.

The handle's exchange does protect complementarity against a detector that
depletes only the channel it reads: the defect is 1e-3 to 1e-2, against 0.5 for
a non-selective detector. But the frozen claim that the remaining defect
**vanishes as Δτ → 0** fails at the frozen window length. At (0.3, 2), where
`Δτ ≈ 1e-4`, the defect is 2.8e-4 at η = 0.5 and 8.2e-4 at η = 0.9, against the
frozen 1e-4.

The table gives the largest `|u + v − 1|` over all doublet states at N = 4000 and
K = 150. N = 2000 and K = 300 agree to every digit shown.

| a | m | Δ (FEM) | Δτ | selective η = 0.1 | η = 0.5 | η = 0.9 | static η = 0.5 |
|---|---|---|---|---|---|---|---|
| 0.5 | 1 | 3.707e-2 | 0.297 | 1.30e-3 | 7.54e-3 | 2.00e-2 | 0.500 |
| 0.4 | 1 | 1.459e-2 | 0.117 | 7.17e-4 | 4.14e-3 | 1.08e-2 | 0.500 |
| 0.3 | 1 | 4.465e-3 | 0.036 | 4.69e-4 | 2.71e-3 | 6.99e-3 | 0.500 |
| 0.3 | 2 | 1.351e-5 | 1.1e-4 | 4.73e-5 | 2.78e-4 | 8.17e-4 | 0.500 |

* **G1, protection, passes.** The selective defect is at most `0.05η` in every
  cell. The worst cell uses 0.45 of that bound, at (0.5, 1), η = 0.9.
* **G2, the limit, fails at both N.** It fails at η = 0.5 and 0.9.
* **G3, the control breaks it, passes.** The static dashpot gives 0.500 to
  three digits in every cell. That is first order in the depletion, as
  predicted: the detector reads u but removes energy from both channels.
* **G4, numerics, passes.**
  * With no detector, the defect is at most 6.2e-11. The exchange is exact.
  * Absorbed work and field energy balance to 3.0e-13.
  * K = 150 and K = 300 agree exactly. N = 2000 and 4000 agree to within 1e-4,
    as frozen.

## What it means

The test compares the handle with the earlier detector models.
* In #367 and #368, every detector coupling tested broke complementarity at
  first order in its strength.
* Here, a detector that removes up to 90% of the energy in the channel it reads
  leaves the complementary reading intact to between 5e-5 and 2e-2.
* The reason is structural. Over `π/Δ` the Klein doublet swaps its ±m channels,
  so the second reading sees the channel the detector never touched.
* The non-selective control shows what is at stake: absorbing from both
  channels breaks complementarity by half the depletion.

The detector's record is also faithful. Its absorbed energy is proportional to
the first reading: `|D₁ − η u|` is at most 9.9e-3 at (0.5, 1) and 2.8e-4 at
(0.3, 2) for η = 0.5. So the depleting detector records u while the
exchange carries `1 − u` to the second visit.

What failed is the frozen account of the remaining defect. It assumed the only
residual was the exchange running during the window, which is proportional to
Δτ. The full-field defect exceeds that two-level expectation by a factor of 3 to
6,000. The excess is largest where Δτ is smallest, so something else sets a
floor.

## Where the floor comes from (post hoc, not a result)

These diagnostics were run after the verdict, at N = 2000, K = 150 and η = 0.5.
They are `alice_exchange_posthoc.py`, archived in `posthoc.json`.

**Not the strip's size or place.** Changing the strip half-width from 0.15 to
1.0 changes the defect by under 1% at (0.3, 2) and by 4% at (0.3, 1).

**Not a profile mismatch in the strip.** At (0.3, 2), the energy-normalised A and
B amplitudes inside the strip differ by only 1.2e-5 to 1.7e-5, far below the
defect.

**It is the detector's switching.** Lengthening the window, while keeping its
sin² shape and the depletion η = 0.5, gives:

| window τ | (0.3, 2) defect | (0.3, 2) leakage | (0.3, 1) defect | (0.3, 1) Δτ |
|---|---|---|---|---|
| 4 | 1.42e-3 | 4.1e-3 | 1.78e-2 | 0.018 |
| 8 (frozen) | 2.78e-4 | 2.0e-5 | 2.71e-3 | 0.036 |
| 16 | 5.07e-5 | 5.6e-8 | 8.22e-4 | 0.071 |
| 32 | 3.67e-6 | 2.5e-9 | 8.14e-4 | 0.143 |

* At (0.3, 2) the defect falls by about 400 between τ = 4 and τ = 32. The
  energy left in other modes after the window falls by six orders of magnitude.
* A detector switched on and off over a few wave periods excites the other modes
  non-adiabatically. Those modes are not tuned to the medium's rotation, so the
  medium acts on their −m parts too.
* At (0.3, 1) the defect stops falling at τ = 32, where Δτ has grown to 0.14 and
  the exchange running during the window takes over.

The picture these diagnostics suggest is that the residual has two parts:
* a switching part, which falls rapidly as `ω̄τ` grows;
* an exchange part, which grows with Δτ.

If so, a selective detector protects complementarity to arbitrary accuracy only
inside the window `Δ ≪ 1/τ ≪ ω̄`: slow compared with the wave period, fast
compared with the exchange. For the m = 2 cell that window spans five decades.
This is a candidate for a new pre-registration, with its own held-out cells and
window lengths. It is not a result.

## Diagnostics

* **Timing from the exact-neck law.** Placing the second visit at `π/Δ_E`
  instead of the measured `π/Δ` adds the expected first-order timing defect.
  The law's 0.05–0.35% error in Δ gives total defects of 8.3e-3 at (0.5, 1) and
  6.9e-4 to 8.6e-4 at (0.3, 2), for η = 0.5.
* **Leakage.** After the frozen window, the energy outside the doublet is about
  2e-5 of the initial energy in every cell.
* **Selectivity of the medium.** In the smoke test, at γ₀ = 1, the
  counter-rotating dashpot removed 87% of the +m channel and changed the −m
  channel by about 1e-4.
  * Absorption is not monotonic in the medium's strength. A very stiff medium
    drags the field with it and can feed energy in, which is the Zel'dovich
    effect.
  * The probe therefore calibrates to the smallest strength that reaches each η.

## Procedure notes

* **The freeze commit carried a failing test.** `0fb3b6e` included a test whose
  threshold was too tight for the two-level scaling. It was fixed in `8ff9df4`,
  which changed no frozen file.
* **The probe was smoke-tested on a non-gated cell only.** The cell was
  (0.35, 1) at N = 400 and 600. That showed a defect about 15× the two-level
  expectation before the gated run. The gates were not changed.
* **Calibration.** It takes the smallest detector strength at which the +m
  channel, averaged over its two phases, loses η of the energy. That implements
  the prereg's "bisection" on the weak-coupling branch.

## Scope

As frozen.
* A classical scalar on #355's fixed surface, with a detector put in by hand.
* The detector's selectivity, a counter-rotating medium tuned to the doublet, is
  an input.
* The exchange is the handle's own.
* No probability law is derived: the measure on doublet states that would give
  Born weights remains an input. Nothing here is claimed to be quantum.
