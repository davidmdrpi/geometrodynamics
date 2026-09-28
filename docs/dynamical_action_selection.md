# Step 4: nonlinear action control and receiver-test readiness

**The available nonlinear control retains a continuous action family. The full
receiver action-selection experiment remains NOT_READY.** This distinction is
the result, not an omitted test quietly counted as a pass.

The requested step 4 is to sweep preparation amplitude and duration and test
whether transferred action approaches robust selected values under changes in
preparation, resolution and receiver definition. Merged #312 supplies nonlinear
initial regions and instantaneous rates, but no evolved localized receiver or
absorbed-action observable. Its profile-overlap problem also remains relevant.
This PR runs a canonical-action control using an existing exact nonlinear
Einstein/quartet reduction. It does not evolve #312's two regions.

## What is measured

Use the nonspherical **spatially homogeneous** four-field Einstein evolution
from `nonlinear_supported_tt.py`, on its expanding FRW-connected branch with
initial A=1.15 and scalar phase zero. This is not the static initial geometry
of #312. There are no spatially separate emitting and absorbing bodies.
The quartet remains assumed matter. No new damping, field, detector threshold,
action unit or Planck constant is inserted.

For a closed loop C of constraint-satisfying initial preparations, evolve
every member independently and measure

    I = Vol(S3)/(2pi) integral_C Theta,
    Theta = -6 A' dA + q'.dq + tr(Pi dM),
    Pi = (H/2) L M^-1, H=A^2-q.q/6, Vol(S3)=2pi^2.

This is **canonical circulation over a preparation loop**. It is neither the
orbital action of one periodic trajectory nor action absorbed by a receiver.
The extra state coordinate eta is a clock, not a new canonical degree of
freedom. Primes in Theta are the stored conformal-time velocities even though
histories are sampled at common proper times. The chosen full S3-volume
normalization is explicit; using the quotient volume would halve these
numbers, without converting the continuous family into a discrete one.

The frozen loop changes relative shape displacement and velocity phase:
initial_data(cos(theta) U, -sin(theta) V, epsilon), with the noncommuting U,V
specified in the freeze. The inherited Hamiltonian and momentum constraints
are completed for every member, including the responsive scalar momentum.
At initial time A and q do not vary around the loop. Direct substitution gives

    Theta(d/dtheta) = H0 epsilon^2 sin(theta)^2,
    I0 = pi^2 H0 epsilon^2, H0=1.1975.

Thus a continuum of initial actions is allowed analytically. Conservation of
full canonical circulation by smooth constrained Hamiltonian flow is the
prospective null prediction. The experiment tests its nonlinear implementation
and shows how sector-only measurements can hide parts of the action ledger.
It is not a search that could establish a new universal action scale.

## Public freeze and results

The [freeze](dynamical_action_selection_prereg.md) was published as
`9f18137441ecb4d9158792e8d4eef11079f0ddea` with draft
[PR #313](https://github.com/davidmdrpi/geometrodynamics/pull/313) created at
2026-09-28 04:22:15 UTC, before the new implementation and numerical runs.
The analytic continuous-action prediction was disclosed before measurements.

All five amplitudes and proper-time samples t=0,.25,.5,1,2 were retained.
The table reports raw action in the stated units, without normalizing the
preparations to the same input action. Values at t=2 use 64 loop nodes.

| Amplitude | Initial/full action | Shape action at t=2 | Scale + scalar at t=2 |
|---|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| .01 | .001181885127 | .001181827311 | 5.78156e-8 |
| .02 | .004727540508 | .004726607509 | 9.32999e-7 |
| .04 | .018910162032 | .018894743162 | 1.54189e-5 |
| .08 | .075640648130 | .075366588397 | 2.74060e-4 |

At epsilon=.08 and t=2 the missing shape circulation is composed of scale
2.6767995e-4 and scalar 6.3797797e-6. These sector contributions are signed,
chart-dependent quantities, not positive local gravitational energies or
absorption events. Changing from shape-only to shape+scalar to the full
canonical ledger changes the apparent redistribution. It is not a change of
physical receiver definition.

Full-action amplitude log-slopes are 2 to about 9e-14 at every sampled time:
doubling amplitude multiplies action by four. The shape slopes at t=2 range
from 1.99594 to 1.99979. Nonlinear scalar and scale contributions scale roughly
as epsilon^4 at small amplitude and have smooth finite-amplitude corrections;
there is no sampled nonzero plateau in these readouts. The exact initial
formula and Hamiltonian circulation identity, rather than four amplitude
samples alone, establish the continuous full-action family in this sector.

| Numerical validity check | Result |
|---|---:|
| Largest relative full-circulation drift, DOP853 | 3.58e-14 |
| Largest sampled normalized constraint residual, DOP853 | 1.23e-15 |
| Largest sampled normalized constraint residual, RK45 | 4.32e-13 |
| Largest N32 to N64 sector difference / initial action | 6.13e-14 |
| Largest RK45 to DOP853 sector difference / initial action | 7.56e-12 |
| Zero-amplitude action | exactly zero numerically |

All seven registered numerical gates pass. Loop grids 16,32,64 use nested
trajectory subsets; they refine phase-loop integration, not the time solver.
The independent RK45 run repeats all 64 trajectories at epsilon=.08. Constraint,
determinant and symmetry maxima are over the saved time samples; positive
chart checks also occur during RHS evaluation. These are numerical checks,
not rigorous global error bounds.

## What this does and does not settle

Two verdicts are deliberately separate:

- `PASS_CANONICAL_CONTINUUM_CONTROL`: the nonlinear constrained dynamics
  preserve a continuous preparation-loop action, with a resolved sector ledger.
- `NOT_READY_FOR_RECEIVER_ACTION_SELECTION`: no operational receiver transfer
  is measured, so no action-selection milestone is awarded or rejected globally.

Canonical circulation conservation cannot exclude discrete coarse-grained
receiver responses, special globally selected histories, nontrivial mouth
boundary conditions or sectors outside this homogeneous family. Likewise,
sector redistribution cannot establish an absorbed quantum. Mode discreteness,
antipodal focusing, initial metric cross-response and raw curvature norms do
not fill that missing measurement.

Observation times in this run vary evolution duration. They do **not** vary
physical preparation-pulse duration. No physical source pulse exists in this
homogeneous experiment. That part of the requested step 4 remains outstanding,
as do receiver-definition controls; the machine-readable readiness verdict
lists them explicitly.

## Reproduction and next decisive test

Raw endpoint states for every finest-grid trajectory, including the RK45
control, are stored losslessly as gzip/base64 JSON in
[`states.json.gz.b64`](../experiments/closure_ledger/runs/20260928_action_selection/states.json.gz.b64).
The diagnostic report, gates and source/archive hashes are in
[`action.json`](../experiments/closure_ledger/runs/20260928_action_selection/action.json).
Decode with `decode_states` in the runner, or base64-decode, gzip-decompress
and parse JSON. No unsafe object/pickle deserialization is used.

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.action_selection_probe
python -m pytest -q tests/test_action_selection.py
# Recompute all action/constraint diagnostics without rerunning trajectories:
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.action_selection_probe \
  --replay experiments/closure_ledger/runs/20260928_action_selection/states.json.gz.b64 \
  --output-dir /tmp/action-replay
```

A zero process exit means numerical control gates passed, not that receiver
action selection was tested successfully. The readiness verdict remains in
the output regardless of that exit code.

The next decisive experiment needs constraint-monitored nonspherical evolution
of separated source perturbations and a responding localized receiver. Define
absorbed action from the retained canonical action and an operational readout,
with field/metric controls and a complete transfer ledger; then sweep actual
preparation amplitude and pulse duration and change receiver definitions.
Any nonzero selected value must persist under those controls without fitting
a unit or rounding outcomes. The #312 seed-overlap and phase-dependent active
density findings belong in that future freeze. No new trapped-object or
quantum milestone is claimed here.
