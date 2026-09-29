# R3 refocusing-resonance sign gate: results

Date: 2026-09-29.

**Review qualification:** see the [dated audit of the freeze and run](r3_refocusing_resonance_review.md).
The registered UNRESOLVED result is retained. Restarted segments have not
been certified to shadow a single Einstein trajectory; a frequency-shift
sign would establish a local direction, not a resonance crossing or closure.
The new linear-only control reproduces finite-window estimator bias.

- Freeze: [`2e984ac`](r3_refocusing_resonance_prereg.md), pushed 00:46:58 UTC.
- Correction: [`6b55c5b`](r3_refocusing_resonance_prereg_correction.md), pushed 00:47:41 UTC, before any code.
- Implementation: `7b1f454`, committed while the run was in progress.
- Archive: `experiments/closure_ledger/runs/20260929_r3_resonance/r3_resonance.json`. It binds the SHA-256 of `r3_resonance.py`, the probe and `nonlinear_supported_tt.py`.

## 1. Registered verdict: UNRESOLVED

Gate **N3 (convergence) failed in all 10 amplitude/polarisation cells**. As a
consequence, **N5 (scaling) had no resolved pair** at eps <= .04. By the frozen
rules, the verdict is UNRESOLVED, and it is not relabelled.

The other recorded diagnostic gates passed everywhere (these are not a
full-state trajectory replay):
- **N1.** Hamiltonian residual <= 8e-11.
- **N2.** |A-1| <= 2.3e-3. The phase-plane radius stays within .84–1.28 eps.
- **N4.** At eps = .01, |rho - rho0| = 1.5e-4, below the 1e-3 limit.

No bracket widening from the correction note was ever used: every window
straddled at its initial bracket.

rho0 = 1.4846664 (frac = .48467). This is the branch nearer the measured rho,
and it agrees with the WKB value 1.4857. The resonance target 3/2 lies above
it, so reaching the resonance requires **c > 0**.

## 2. Measured rotation numbers

Columns:
- **rho(48):** rotation number from the weighted Birkhoff average over K = 48 periods.
- **rho(24):** the same average over the first 24 periods.
- **\|prim - sec\|:** difference between the two integrator tolerances.
- **err:** the registered N3 error estimate.
- **Delta:** rho(48) - rho0.
- **c:** Delta/eps^2.

| eps | pol | rho (48) | rho (24) | \|prim - sec\| | err | Delta | c | resolved |
|---|---|---|---|---|---|---|---|---|
| .01 | + | 1.484518 | 1.484160 | 1.1e-10 | 3.6e-4 | -1.48e-4 | -1.48 | no |
| .01 | - | 1.484516 | 1.484156 | 7.6e-11 | 3.6e-4 | -1.51e-4 | -1.51 | no |
| .02 | + | 1.484151 | 1.483732 | 4.3e-10 | 4.2e-4 | -5.15e-4 | -1.29 | no |
| .02 | - | 1.484128 | 1.483705 | 2.7e-10 | 4.2e-4 | -5.39e-4 | -1.35 | no |
| .04 | + | 1.482696 | 1.482079 | 6.2e-10 | 6.2e-4 | -1.97e-3 | -1.23 | no |
| .04 | - | 1.482499 | 1.481860 | 5.5e-10 | 6.4e-4 | -2.17e-3 | -1.35 | no |
| .08 | + | 1.476735 | 1.476118 | 2.0e-8 | 6.2e-4 | -7.93e-3 | -1.24 | no |
| .08 | - | 1.475024 | 1.474536 | 3.1e-9 | 4.9e-4 | -9.64e-3 | -1.51 | no |
| .16 | + | 1.455268 | 1.455415 | 5.5e-10 | 1.5e-4 | -2.94e-2 | -1.15 | yes (descriptive) |
| .16 | - | 1.441267 | 1.441306 | 2.8e-9 | 3.9e-5 | -4.34e-2 | -1.70 | yes (descriptive) |

rho does not cross 3/2 anywhere on the ladder. It moves away from 3/2 at
every amplitude.

## 3. Why N3 failed: a design error in the freeze

The two DOP853 tolerance settings agree to about 2e-8. The failure comes entirely
from |rho(48) - rho(24)|. This section is a post-hoc diagnostic
(`experiments/closure_ledger/r3_posthoc_beat.py`), not part of the registered
test.

- **A near-resonant beat limits this finite estimator.** rho0 lies near the 1:2
  resonance of the clock-section map: 2 rho0 mod 1 = .9693. The phase
  increments therefore carry a beat with period 1/(1 - .9693) = 32.6 clock
  periods. That is longer than the 24-period half-window and suggests a
  convergence problem, but is not by itself an error bound. The review's
  independent linear control gives a half-window discrepancy of 3.381e-4
  and a 48-period bias of -2.405e-5 for this angle and initial phase.
- **The problem was foreseeable.** It follows from the linear value rho0 ≈ 1.485
  and from the target sitting exactly at the 1:2 resonance. It should have
  been caught before the freeze. It was not, and the error is mine.

## 4. What the data indicate, without a verdict

The following is descriptive. It is not a registered outcome.

- **Every cell has the same sign.** Delta < 0 in all 10 cells, and in the
  unweighted raw means as well.
- **The descriptive coefficient varies.** c lies between -1.15 and -1.70.
  These unresolved estimates do not establish an asymptotic coefficient.
- **The raw-mean shifts grow quadratically.** Successive ratios are 3.7, 3.7,
  3.4 and 3.8 per doubling of eps.
- **Only eps = .16 clears the resolution criterion.** There |Delta| is 200x and
  1100x the error estimate, and it is at the wrong sign.
- **Implication if it holds.** The nonlinear shift points away from the
  refocusing resonance, which is the frozen FAIL direction. A certified
  leading coefficient of this sign would show initial motion away from the
  target, without excluding a later turn or closing the R3 route.

This is not adopted as a result. At eps <= .04 the shifts are the same size
as the estimator bias, so the registered test cannot decide.

## 5. Proposed prospective extension (not run)

A replacement must also address the trajectory, phase-coordinate and
decision-logic issues in the dated review before new nonlinear runs. Two
possible estimator improvements, neither guaranteed to resolve the test:
1. **More periods.** Use K = 256, with a half-window of 128, about four beat
   periods. Validate the error rather than assuming this duration suffices.
2. **Cancel the bias directly.** Measure Delta against a quasi-linear run at
   eps = .001 that uses the same estimator and the same K.

Any extension needs a separate prospective specification. A finite-amplitude
reference subtraction can retain nonlinear bias and does not cure restart
discontinuities. The historical decision rule also needs its interpretation
limited to the direction of a resolved shift.

## 6. Reproduction

    python -m pytest -q tests/test_r3_resonance.py
    python -m experiments.closure_ledger.r3_review_controls             # archive replay + linear controls
    python -m experiments.closure_ledger.r3_resonance_probe --workers 4 --output /tmp/r3-new-run.json  # historical protocol; not a corrected experiment
    python -m experiments.closure_ledger.r3_posthoc_beat

The six tests cover:
- constraint-satisfying section data;
- invariance of the LRS subsystem;
- agreement of the linearised section map with #310's tensor n=2 map (trace to 1e-5);
- the Birkhoff estimator on a synthetic quasi-periodic signal;
- re-scoring of the archive and its source binding;
- sensitivity of the score to tampered increments.

A smoke test at eps = .16 over 8 periods ran before the registered run, to
check the code executed. It is disclosed here and no decision used it.
