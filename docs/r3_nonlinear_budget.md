# Nonlinear diagonal detuning with a fixed corrective-reset budget

The prospective hypothesis **failed**: all four controlled cases exceeded the
pre-stated reset budget after one two-return step, independently in DOP853
and Radau. No over-budget reset was applied. The aggregate label is
`CONTROLLED_NONLINEAR_BOUND_FAILED`.

This is a failure of this particular linear tuning/controller, coordinate
tube and budget for these four initial conditions. The post-hoc review below
shows that the design confounds nonlinear detuning with amplified initial
tuning error. It does not discriminate intrinsic robustness of the nearby
family. Autonomous nonlinear stability and dynamical action selection remain
`NOT_ESTABLISHED`.

## Prospective commitment

The [specification](r3_nonlinear_budget_prereg.md) was published in commit
`75bed15f1afc7c932462d2721c5e76d3c2546458` before these trajectories.
The cases, tolerances, controller, budgets, horizons and labels were not
changed after observing outcomes. The experiment follows #320 at
`80d894a9304b9c1cbe4db51b955865a0083d6f4b`; #321 is stacked on that branch.
The measured sources and specification are hashed in `provenance.json`.
Only the prose specification was published before trajectories. The producer
and dynamics implementation first appeared in measurement commit `6504e1c`;
their hashes bind what ran but do not constitute a pre-run code freeze.
For example, the departure interval `1e-6` implements the specification's
"short initial departure" without having been fixed numerically in that
specification. This limits the strength of the prospective provenance claim.

The finite-power result in #320 excluded the physical hyperbolic pair and
conditioned away marginal action drift. This study instead evolves the
unchanged full nonlinear 29-state equations with diagonal initial data,
records actual drift, and charges for external interventions. It makes no
off-diagonal or inhomogeneous nonlinear claim.

## Norm, fine tuning and budgets

On descending clock sections the coordinate vector is

\[
w=(A,p_A,x,p_x,y,p_y,q'_0),\qquad
p_A=-6A',\quad (p_x,p_y)=A^2(x',y').
\]

All distances and reset costs use its unweighted Euclidean norm in the
dimensionless normalization of the model. This is a declared coordinate
norm, not an invariant physical distance. The Hamiltonian-completed clock
momentum is included in every cost. Distance is minimized over the periodic
cubic interpolation of #319's family, while the fitted phase is separately
unwrapped; phase fitting never changes the state.

The initial scale coordinates are tuned to annihilate the inherited linear
hyperbolic components. This is deliberate fine tuning, not dynamical
selection or construction of a nonlinear centre manifold. A trial loop's
canonical action offset defines each input detuning; that loop is not proved
invariant and its action is not a measured conserved action of the trajectory.

Let `d0` be initial distance to the family. The frozen tube radius is `2*d0`,
each proposed reset must cost at most `0.02*d0`, and the sum of applied reset
magnitudes must remain at most `0.25*d0`. Costs cannot cancel. Initial tuning
is reported separately. Resets are external changes to the state, not GR
forces. Both the pre-reset and proposed post-reset distances must satisfy
the tube bound. A budget or tube failure in either integrator counts as a
resolved negative only when independently confirmed by the other integrator.

## Horizon tied to twist

Quadratic interpolation of the authenticated #319 values of `I(omega)` gives
`d omega/dI = 4.908582988640943` at `omega=pi`. The prescribed number of
two-return steps is

\[
N=\left\lceil\frac{1}{2|(d\omega/dI)\Delta I|}\right\rceil,
\]

giving 102 or 204 steps, approximately one radian of predicted drift.
A positive outcome additionally required at least 0.5 radian of measured
unwrapped drift by that horizon. Because the trial loop is not proved
invariant, the empirical twist alone cannot certify adequate exposure.

The hyperbolic multiplier is about 7242 per two returns. The linear estimate
`log10(2*d0) - N*log10(7242)` is about -395.54 or -789.55: unforced shadowing
over the proposed horizons would require an unstable initial component of
roughly `10^-396` or `10^-790` in the chosen normalization, within that linear
estimate. This is not precision delivered by the present double-precision
calculation and not a nonlinear stability theorem.

## Measured result

Numbers below use DOP853; Radau independently gives the same terminal
reason and step for every controlled case.

| Input action offset | `d0` | Target steps | Initial tuning norm | First proposed cost / `d0` | Cost / per-step cap |
|---:|---:|---:|---:|---:|---:|
| -0.0010 | 0.0073373182 | 102 | 0.0000268002 | 1.5031722 | 75.16 |
| -0.0005 | 0.0036448931 | 204 | 0.0000133088 | 0.7470815 | 37.35 |
| +0.0005 | 0.0035992626 | 204 | 0.0000131335 | 0.7370501 | 36.85 |
| +0.0010 | 0.0071546738 | 102 | 0.0000260985 | 1.4631299 | 73.16 |

Every controlled case terminates at step 1 with `CONTROL_BUDGET_EXCEEDED`.
Each first proposal exceeds both the per-step and total budget. Applied
cumulative reset cost is zero. Initial tuning is the only intervention
actually performed. Pre-reset distances are respectively 1.8290, 1.2699,
1.2433 and 1.7671 times `d0`, within the tube at this step. The recorded
failure is therefore the cost of this controller's proposed correction,
not a demonstrated need for every possible controller to spend that amount.

The largest controlled discrepancy between integrators over compared states,
distances and costs is `1.06e-10`, below the case-specific limits
`3.60e-7` to `7.34e-7`. Sampled constraint residuals on these completed
controlled histories are at most `3.60e-13`, below the `1e-8` gate.
Both solvers used `rtol=2e-12`, `atol=2e-14` and `max_step=0.025`.

All eight unforced runs reach `NUMERICALLY_UNRESOLVED` on step 2 under the
registered two-return scoring rule. However, each archive contains a valid
first return of that step, already outside the tube, independently confirmed
by both integrators. The original report overlooked these completed partial
histories when it said the failed attempts were not evidence of escape.
That statement is corrected here. The second return is not reached; the
archive does not establish its subsequent fate or a collapse singularity.
The registered labels remain unchanged, and the aggregate negative still
rests on the controlled budget failures.

No case completed its target horizon. Measured drift after the first step
has magnitude only 0.00336 to 0.00813 radian. Early failure was an explicit
stopping rule, so it falsifies the specified bounded-control hypothesis;
it cannot support a long-horizon or autonomous claim.

## Review follow-up: what the negative result measures

This section is post hoc, prompted by the
[independent review](https://github.com/davidmdrpi/geometrodynamics/pull/321#issuecomment-6051205277).
No new trajectories were run. The specification, measured sources, archives
and registered labels are unchanged.

**Quadratic cost and poor discrimination.** The first proposed reset costs
divided by `d0^2` are 204.8667, 204.9666, 204.7781 and 204.4999 in case order.
Both integrators reproduce this near-quadratic scaling. It is consistent
with an O(d0^2) residual from linear initial tuning that grows along the
hyperbolic direction before the first correction. The coefficient near 205
is an effective two-return correction-cost coefficient in this norm: it
already includes evolution and controller projection. It is neither a
direct measurement of centre-manifold curvature nor a coefficient to
multiply by 7242 again. Establishing the manifold mechanism separately
requires nonlinear invariant-manifold calculations.

The review's feasibility estimate is compelling **within that local scaling
model**. Combining `d0 ≈ 7.2*abs(delta_I)` and `cost ≈ 205*d0^2` with the
256-step cap gives `abs(delta_I) >= 3.979e-4` and a first cost near
`0.587*d0` at the smallest permitted detuning, about 29 times the per-step
cap. This exposes why reducing the tested detuning modestly would not rescue
the design. If the same cost recurred at every step, the cumulative cost
would be approximately `150*d0`, roughly independent of detuning.

Those are extrapolations from four first-step measurements, not a proof of
impossibility for every admissible detuning or a measured cumulative cost.
No reset was applied, and no controlled later step was observed. The
registered negative stands, but should not be used as evidence against
nearby invariant tori. A future design needs a feasibility/pilot stage
disclosed separately before freezing both its specification and producer.

**Correction timing.** The reviewer reports that pre-step A corrections of
about `3.4e-7` at +0.001 and `8.8e-8` at -0.0005 reduce the subsequent
proposed reset costs to about `6.9e-5*d0` and `3.3e-5*d0`, respectively.
These are the reviewer's post-hoc one-step diagnostics, not independently
rerun here or incorporated into our registered evidence. They support
investigating initialization and timing rather than interpreting the failed
delayed controller as intrinsic nonlinear fragility. They do not establish
budget compliance, stability or a predictable positive outcome over the
102/204-step horizons.

**Valid intermediate unforced escape.** The added archive-only diagnostic
validates the completed first return within each failed second step:
connection to the previous state, time ordering, chart, constraints and
descending section. It checks DOP853/Radau endpoint agreement and computes
the global range of the reference spline's A coordinate, including all
interior extrema: `[1.00043661784, 1.00050323387]`.

| Input action offset | A at the third descending return (DOP853) | Lower bound on distance / d0 from A alone | Maximum endpoint discrepancy between integrators |
|---:|---:|---:|---:|
| -0.0010 | 0.8962493834 | 14.1996 | 5.88e-10 |
| -0.0005 | 0.9736829327 | 7.3400 | 2.84e-10 |
| +0.0005 | 0.9743871616 | 7.2374 | 5.64e-10 |
| +0.0010 | 0.9013404345 | 13.8506 | 1.22e-9 |

Every lower bound exceeds the tube radius `2*d0`; this conclusion does not
depend on finding a nearest phase. Sampled constraint residuals in these
partial returns are at most `7.43e-13`. The agreement is numerical, not exact
identity. This is resolved finite-time departure of these linearly tuned
unforced initial conditions, consistent with the inherited hyperbolic
instability. It does not establish instability of motion restricted to a
nonlinear centre manifold. The post-hoc observation supplements, rather than
replaces, the registered `NUMERICALLY_UNRESOLVED` full-step status.

**Related evidence and next question.** The review points to
[`r3_breaking.md`](https://github.com/davidmdrpi/geometrodynamics/blob/claude/geometrodynamics-qft-audit-vpktax/docs/r3_breaking.md).
That report finds two resonances unbroken at its registered resolution,
while its post-hoc analysis suggests a small harmonic-10 breaking signal
at the LRS 2/5 resonance. It explicitly does not identify an extra integral
or establish exact integrability. We inspected the report; this follow-up
does not independently replay that separate experiment.

A useful next study would test a candidate extra integral or a normal-form
remainder on held-out amplitudes and phases, with pre-stated accuracy and
failure thresholds and both code and specification published before the
validation runs. A merely successful one-step predictive correction or
an in-sample normal-form fit would not suffice. Two nearly unbroken
resonances do not guarantee a positive full-horizon controller result.

## Evidence and verification

The evidence directory is
`experiments/closure_ledger/runs/20261005_r3_nonlinear_budget/`.
Four gzip/base64 JSON case files contain initial loops and states, completed
sampled full-state histories, proposed corrections, application flags and
terminal records. The encoding is lossless and decoded by the producer's
`read` function. `result.json` is a summary, not the source of decisions.

The manifest SHA-256 is
`4ea7db2a0f0dc002bb1a4d0f7a12375d9c9f5a38d719d395741efdb4d51e44e1`.
From the repository root, after installing the project dependencies:

```bash
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.r3_nonlinear_budget_replay
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.r3_nonlinear_budget_review
pytest -q tests/test_r3_nonlinear_budget.py tests/test_r3_normal_quotient.py
```

Replay authenticates parent evidence and measured sources and recomputes
initial loop actions, tuning, horizons, phase fits, distances, proposed
corrections, costs, constraints and terminal decisions from archived states.
It is a diagnostic replay, not a fresh integration of the differential
equations; independent integration was performed during production.
Tests include reachable tube, per-step budget, cumulative budget and drift
failures, a synthetic permitted outcome, clock-completion cost, refused
over-budget application, altered histories, fabricated pass flags, corrupted
bytes, changed labels and missing/reordered cases or returns.

Production used Python 3.12.14, NumPy 2.5.3 and SciPy 1.18.1. The 21 new
tests pass there and with NumPy 2.3.5 / SciPy 1.17.0; the new and parent
normal-quotient suites pass 38 tests together in the production environment.
The review follow-up adds five tests for the partial-return evidence and
tamper rejection; the combined suite now has 43 tests.

The portable replay was added after measurement. The original producer's
replay regenerated initial states before validating phase drift. Under the
alternate libraries, roundoff in that regeneration moved the initial fitted
phase by up to `3.10e-9` radians, exceeding its `1e-9` drift-reproduction
tolerance. The portable replay independently verifies the regenerated setup
against the archive, checks the fitted phase of the actual recorded state,
then replays from that recorded initial state and phase. It does not loosen
budget, state, cost, constraint or terminal-decision thresholds. Tests reject
altered initial states and phases. The original measured producer, dynamics,
specification, archives and their hashes remain unchanged.

Further controller or initialization changes must be a separately specified
experiment, with these failures retained.
