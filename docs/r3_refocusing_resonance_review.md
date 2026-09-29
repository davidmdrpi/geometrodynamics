# Dated review of the R3 freeze and its completed run

Date: 2026-09-29. Reviewed freeze: `2e984accf42cd649ae54bbba8427e909a33924c3`.
The review also inspected correction `6b55c5b`, implementation `7b1f454`,
and results `acd38b2`, which already existed when this review began.
This is retrospective. Neither preregistration nor the original producer,
archive, or registered **UNRESOLVED** verdict is changed. No new nonlinear
amplitude scan was performed. Any replacement experiment needs a new freeze.

## Findings that affect the experiment

### 1. A sign test establishes a local direction, not a resonance crossing

The freeze's PASS sentence, "The resonance is reached at finite amplitude,"
does not follow from the sign of c. For example, with d=3/2-rho0>0,

    rho(epsilon) = rho0 + (d/2) epsilon^2/(1+epsilon^2)

has positive quadratic coefficient but never reaches 3/2. The estimated
sqrt(d/c) may also lie outside the range in which the expansion is accurate.
A negative c establishes the leading direction away from the target; it
does not exclude a later turn or close the whole R3 route. Appropriate
future labels are `SHIFT_TOWARD_TARGET`, `SHIFT_AWAY_FROM_TARGET`, and
`UNRESOLVED`, with crossing and closure assessed separately. No historical
label is re-scored using these proposed meanings.

### 2. Rational phase advance does not establish full-state return

An averaged rotation number of 3/2 need not imply that every tensor state
is inverted at each section or that A, A', q', x and x' all return at the
second section. Those implications require additional dynamical structure.
Even for a periodic orbit, inversion is not guaranteed when x -> -x is not
a symmetry. The appropriate closure residual is P^2(z)-z for the full,
constraint-reduced section state z, with the coordinate-time accumulator
excluded and any allowed gauge/field identification specified in advance.

Nor are Hamiltonian periodic histories necessarily isolated in amplitude.
The elementary harmonic oscillator has a continuous family at one fixed
period; two uncoupled oscillators with frequency ratio 3:2 retain continuous
amplitude families at that winding. A nondegenerate twist/continuation
argument and physically justified boundary conditions are needed here.

The freeze also conflates the free conformal control with #310's supported
tensor. The supported n=2 map has trace about -1.990725, not -2, and its
refocusing operator is not exactly -I. The ratio 3/2 is a candidate target
inspired by the free control, not an established nonlinear identification.

### 3. Restart kicks have no demonstrated trajectory error budget

`track` changes A and re-solves q0' after every four accepted periods. This
keeps the Hamiltonian constraint small but produces a sequence of segments
with state jumps. It does not by itself prove a nearby single trajectory
or an invariant centre manifold. The review archive check finds later A
kicks as large as 4.12490e-7. Their small size alone bounds neither their
effect on phase nor accumulated bias in Delta/epsilon^2.

The saved record contains phase increments, sampled extrema and scalar
constraint residuals, but no segment start/end states or q0' jumps. Those
data can reproduce the score, not independently verify the flow, continuity
at restarts, or phase reconstruction. Accordingly the new replay reports
`diagnostic_replay=VERIFIED`, `trajectory_shadowing=NOT_ESTABLISHED`.
The latter is an evidence limitation, not a proof that no shadowing orbit
exists. A centre-stable finite-time straddle is also not automatically a
bounded two-sided centre-manifold orbit.

The bisection implementation differs from the freeze: `classify` stops at
9*pi+.5 and, absent an exit, reads the final A. It does not stop at the ninth
clock crossing. It treats solver failure as collapse, conflating numerical
failure with a physical branch. The widening implementation multiplies by
4 and aborts before trying the stated cap when the next step exceeds it
(.16 -> .64 skips .32; .000256 -> .001024 skips .001). No widening occurred
in the saved run, so that last discrepancy did not change its recorded
outcome. These facts belong in a replacement protocol, not a silent edit
to the archived producer.

### 4. The finite estimator has a measurable linear bias

The independent **linear-only** control uses the same clock phase, angle,
48 periods and sampling as the completed nonlinear run. Linear phase
unwrapping selects the Floquet branch without consulting a nonlinear run.
DOP853/256 samples and RK45/512 samples give:

| Quantity | Linear control |
|---|---:|
| rho0 from the elliptic trace, branch fixed by linear winding | 1.484666408416 |
| weighted rho(24) | 1.484304250755 |
| weighted rho(48) | 1.484642354781 |
| rho(48)-rho0 | -2.40536e-5 |
| abs[rho(48)-rho(24)] | 3.38104e-4 |
| difference between integrators in rho(48) | < 5e-14 |

Thus the original N3 failure is already plausible in the linear estimator;
the beat explanation is supported by a direct control. But "24 periods
cannot converge" is not a general theorem: the error also depends on phase
coordinate and oscillatory coefficients. Replacing theta by the smooth
degree-one angle theta+.2 sin(2 theta) changes this finite rho(48) to
1.484868348658, although it leaves the limiting rotation number unchanged.
Finite coordinate independence in the freeze is therefore false.

Weighted Birkhoff superconvergence requires hypotheses including smooth
quasiperiodicity and a Diophantine rotation vector; it is not guaranteed
at the rational target or for kicked trajectories. See Das and Yorke,
[Theorem 1.1](https://arxiv.org/html/1506.06810v5). A longer window or an
epsilon=.001 subtraction may help, but neither is guaranteed to resolve
the nonlinear coefficient. Bias changes with amplitude and may not cancel.
The historical "secondary integrator" is the same DOP853 method at looser
tolerances; it is a tolerance check, not an independent integration method.

### 5. Even powers of the signed turning displacement were not derived

For b0=diag(1,1,-2)/sqrt(6), tr(b0^3)=-1/sqrt(6). The LRS potential lacks
an x -> -x symmetry. A normal-form frequency may be smooth in an invariant
action I, but the initial coordinate displacement is not I: generically
I=a epsilon^2+s_pol b epsilon^3+... . Thus rho(epsilon) can contain a
polarisation-dependent cubic term. The O(epsilon^4) remainder in the freeze
needs proof or replacement by a controlled expansion in I. Comparing both
signs is useful; identifying equal signed magnitudes with equal-action
orbits is not justified. The finite-difference N5 slope check alone does
not supply that derivation or a bound on the remainder.

## R2 claims and scratch-script corrections

Odd maps S3 -> S3 have odd degree when the normalized map is defined
([Hatcher, Algebraic Topology, section 2.B](https://pi.math.cornell.edu/~hatcher/AT/ATpage.html)).
The parity restriction is sound. A generic isolated antipodal pair of
simple zeros can change the degree by ±2; simultaneous pairs, degenerate
zeros and the homogeneous collapse need not have exactly that jump.
At a zero, the normalized map and its degree are undefined.

Continuity rules out a universal identity J=lambda*N on a connected family
crossing an N jump if J remains continuous there and lambda is a fixed
nonzero constant. It does not rule out asymptotic, event-conditioned or
restricted-family observables for which those assumptions fail. The scratch
ansatz does not evolve Einstein data and contains no absorbed-action ledger.
It cannot close R2 or make R3 the only remaining mechanism. Similarly #314
excludes its registered homogeneous plateau family, not all R1 mechanisms.

The original event-time equality was inserted by the code, not measured:
it set tau=epsilon*lambda/sqrt(3). For the actual static ansatz
R(tau)=-(sqrt(3)/2)sin(2 tau), the central simple root is

    tau = (1/2) arcsin(2 epsilon lambda/sqrt(3)),
    |2 epsilon lambda/sqrt(3)| < 1.

Hence tau/epsilon is constant only to leading order. For lambda=1 and
epsilon=.001,.01,.05 it is .5773503975, .5773631000, .5776715014.
The Hopf control has |R x+epsilon Jx|^2=R^2+epsilon^2 for unit x and an
orthogonal complex structure J; its absence of zeros at nonzero epsilon
follows directly, without a grid search.

`r2_degree_intervals.py` originally failed with `FileNotFoundError: deg.py`.
The successor scripts now use explicit imports, defer work until invoked,
use the exact central-branch time, and label their output `KINEMATIC_ONLY`.
They reject sampled zeros before normalizing and include the identity and
Hopf controls. Their finite root search and grid quadrature still do not
certify root completeness or integer degree. The original pre-freeze bytes
remain available at `2e984ac`; the repaired scripts are not represented as
the code originally run before that freeze.

## Concrete replacement requirements before another nonlinear run

1. Derive the constrained return map or normal form and its amplitude/action
   convention. Fix the linear rotation branch from linear dynamics alone.
2. Replace kicks with a specified collocation/multiple-shooting/invariant-map
   solve, or bound their effect with a variational shadowing analysis. Freeze
   continuity, residual, restart/precision and frequency-error thresholds;
   store the full section and pre/post-segment states. Constraint residuals
   must be monitored along the segments as well as at their endpoints.
3. Validate the estimator at the actual near-resonant linear frequency,
   vary an admissible phase coordinate, and use a genuinely different
   integration method. Keep numerical failure distinct from a collapsing
   solution, and implement the stopping section actually registered.
4. Report a resolved local shift direction first. Claim a crossing only
   after a bracketed continuation with error bounds; claim closure only
   after the full return residual passes independently. Discrete history
   families and a selected action require further physical justification.

No numerical thresholds for a new nonlinear experiment are retroactively
chosen here. The current scientific verdict remains **UNRESOLVED**.

## Reproduction and preservation

The original nonlinear archive fingerprint is pinned in the review checker.
Missing/duplicate cells, short or nonfinite arrays, invalid bounds, missing
restart records, and jointly altered evidence/results cannot authenticate.
This authenticates the historical record; it does not recreate missing
trajectory evidence or certify the old producer for future data.

```sh
python -m experiments.closure_ledger.r3_review_controls
python -m experiments.closure_ledger.r3_prefreeze.r2_degree_intervals --grid 6 --seeds 8
python -m pytest -q tests/test_r3_resonance.py tests/test_r3_review_controls.py
```

The small-grid command is a reproducibility smoke check, not a physical
degree measurement. The linear-control record is
`experiments/closure_ledger/runs/20260929_r3_review/review.json`. No nonlinear
result is substituted or re-scored under looser thresholds.
