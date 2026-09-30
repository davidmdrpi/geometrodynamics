# R3 beyond the local LRS coefficient: results

Date: 2026-09-29.
- Specification and code: [`8584346`](r3_extension_prereg.md). It became public with the merge push at 04:58:16 UTC, before any cubic five-polarisation jet or new-amplitude circle.
- Archives: `experiments/closure_ledger/runs/20260929_r3_extension/` (`part_A.json`, `part_B.json`, `result.json`). They bind the SHA-256 of all six sources.
- Re-score test: `tests/test_r3_extension.py`.
- Authenticated replay (added after the #317 review): `python -m experiments.closure_ledger.r3_extension_replay --full` returns VERIFIED, with the largest re-evaluated invariance residual 4.9e-12. See [the correction note](r3_extension_prereg_correction.md).

## Summary

| part | question | frozen label |
|---|---|---|
| **B** | leading twist for all five homogeneous n=2 polarisations | **SOME_POLARISATION_SHIFTS_TOWARD** |
| **A** | does the LRS family turn at large amplitude? | **NO_TURN_IN_FAMILY** (27 circles through a = .3044; continuation beyond that failed numerically, unresolved) |
| D | which closure rationals does the LRS family pass? | descriptive; candidate rationals only, periodic histories unverified |
| C | inhomogeneous n >= 3 tensor, vector and scalar sectors | NOT_TESTED |

**The R3 refocusing closure is not ruled out.** #316's SHIFT_AWAY_FROM_TARGET
holds for every *linear* polarisation. For the LRS family it continues to
hold at large amplitude, through the last accepted circle at a = .3044.
Beyond that point, continuation failed numerically and the question is
unresolved. But the homogeneous n=2 tensor also has complex (elliptical)
polarisations, and some of them shift *toward* the resonance. By the
frozen continuation rule, that calls for a separate prospective crossing
and closure test. It does not establish a crossing.

## Part B: all five polarisations

Every check B1–B9 passes:

| check | value |
|---|---|
| B1 fixed point | 0 |
| B2 linear | theta = theta0; off-block 0; block spread 7.8e-15 |
| B3 dissipative quartic part | 1.3e-13 (twist scale about 4) |
| B4 SO(3) covariance | 1.0e-14 |
| B5 nu(LRS) against #316 | -0.9501180968942944, vs -0.9501180968942993 |
| B6 spread over linear polarisations | 1.2e-14 |
| B7 Richardson convergence | 3.8e-13 |
| B8 smallest relative divisor; \|mu^k - 1\| | .988; >= .192 |
| B9 constraint jet | 2.5e-8 on coefficients up to 5.0e5 |

**Twist range.** Over the unit sphere of C^5, nu runs from **-1.8441** to
**+0.8480** turns per unit action. The Sym^2 bounds (+2.87, -4.10) are
consistent with this but not tight. So no certificate of a uniform sign
was possible, as expected once the maximum is positive.

Descriptive, computed after the label was fixed:

| polarisation | nu | gravitational angular momentum |
|---|---|---|
| any linear (real) b, including LRS | -0.950118 | none |
| **diagonal circular**, W ∝ diag(1, ω, ω²), ω = e^{2πi/3} (the maximum) | **+0.847972** | none |
| helicity m = 1 about an axis, (xz) + i(yz) | +0.174947 | yes |
| helicity m = 2 about an axis, diag(1,-1,0) + i(xy) (the minimum) | -1.844127 | yes |

**The maximiser is simple to describe.** It is a shape oscillation that
cycles the anisotropy through the three principal axes: beta(t) ∝
Re(e^{it} diag(1, ω, ω²)). It lies inside the exact diagonal Bianchi IX
subsystem, with M and L diagonal. It carries no angular momentum, so the
momentum-constraint argument of the specification §2 is not even needed
for it.

**The frame-rotating wave shifts away most strongly.** The m=2 helicity
state, the natural "circularly polarised" homogeneous wave, has the
strongest shift away.

**Perturbative crossing estimate.** For the maximiser,
I* ≈ (1/2 - .484666)/.848 = 0.018. This comes from the leading coefficient
alone. It is not a crossing claim.

## Part A: LRS at large amplitude

27 invariant circles were accepted on the registered ladder,
a = .004 to .3044 (I = 2.3e-5 to 0.128). Each has a residual <= 5e-12,
a Fourier tail <= 3.3e-11, and constraint residuals <= 1e-10.
- **frac(rho) decreases monotonically**, from .484645 to .358709. The
  largest step in omega is -5.7e-5, so it never increases.
- **Continuation then failed.** The registered retries went as follows:
  - a = .2560 failed and its retry at .2348 passed;
  - .3320 failed and its retry at .3044 passed;
  - .3620 failed (integrator: step size underflow), and its retry at .3320
    failed with "no real clock velocity".

**Scope of the label (corrected 2026-09-29, after the #317 review).**
NO_TURN_IN_FAMILY covers the 27 accepted circles through a = .3044. It
says nothing beyond them. The endpoint is a failure of **numerical
continuation**, not a demonstrated physical boundary. The final retry's
initial predictor (the a = .3044 circle rescaled to .3320) has positive
q'^2 at every grid point, between 2.114 and 2.244. The predictor for the
failed a = .3620 step is also positive everywhere, between 1.907 and 2.140. So the constraint
failure occurred at a later Newton trial state. It does not show that
section data or invariant circles cease to exist. An earlier version of
this section claimed that they "do not exist beyond a ≈ .33". That claim
is withdrawn.

At the last accepted circle the rotation has moved about 0.126 turns away
from the target. Whether the family continues, and whether it turns,
beyond a = .3044 is unresolved.

## Part D: closure conditions reached

The LRS family's frac(rho) range [.358709, .484666] contains the
rationals with q <= 12:
**5/11, 4/9, 3/7, 5/12, 2/5, 3/8 and 4/11**.

These are **candidate closure rationals only: no periodic history has
been verified** (corrected 2026-09-29, after the #317 review). The code
lists the rationals that lie between sampled rotation numbers. It does
not solve P^q(z) = z.

An earlier version of this section, and §4 of the specification, invoked
Poincaré–Birkhoff to assert closed histories at each rational and closure
actions dense in the family. Those claims are withdrawn. That theorem
needs an invariant annulus of an area-preserving map, with twist
conditions on its boundaries. Those hypotheses have not been established
for this four-dimensional return map over the measured range: the centre
manifold, its induced area form, and the boundary twist were not checked.

For example, 2/5 lies between the rotation numbers sampled at a = .2348
and a = .2792. Whether a period-5 orbit exists there is untested. The
natural low-order closures (0, 1/4, 1/3, 1/2, 2/3, 3/4, 1) all lie
outside the sampled range of the linearly polarised family.

## Part C: not tested

- **Inhomogeneous n >= 3 tensor, vector and scalar sectors.** They need
  nonlinear inhomogeneous GR, which this code does not evolve.
- **Scalar n = 2.** It is linearly hyperbolic (#310).
- **Homogeneous quartet rotations.** Neutral, with nothing to twist.

No conclusion here extends to them.

## What follows

1. **For the elliptical polarisations, the R3 question is open again.**
   #316 and #317 apply to linear polarisations only. That now includes
   large amplitude in the LRS case.
2. **Next step under the frozen rule:** a separate prospective freeze for
   the diagonal Bianchi IX subsystem, where the toward-shifting mode lives
   exactly. It would:
   - continue the diagonal-circular family, which is a relative
     equilibrium at quartic order, perturbed at higher order by the
     three-fold symmetry of the axes;
   - bracket any crossing of frac = 1/2, with error bounds;
   - test the full-state closure P^2(z) = z by multiple shooting.
3. **Selection would still not follow.** Even a verified closed history
   would not establish selected or quantised action. That still needs a
   physical reason for the closure condition. Part D lists candidate
   closure rationals only. It does not show that closure is dense.

## Disclosures

- **Aborted first launch.** Both parts were first started a few seconds
  before I noticed that the push of `8584346` had been rejected: #315 had
  been merged into this branch in the meantime. I killed both processes
  within about a minute. They had written no archive and printed no
  output. I merged the remote branch, confirmed the freeze was public
  (04:58:16 UTC), and restarted both parts unchanged (04:58:38 UTC). The
  merge touched none of the six bound sources.
- **Scoring order.** Part B finished first. It was scored with the frozen
  `score_b` before Part A ended. No rule or code changed between the two.
- **Post-hoc material.** The named-polarisation table and the I*
  estimate were computed after the label was fixed. They are descriptive.
