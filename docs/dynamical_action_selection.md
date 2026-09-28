# Poincaré–Cartan implementation check and receiver-test readiness

**This run checks implementation of the Poincaré–Cartan integral invariant.
It supplies no evidence for or against receiver action selection. The full
receiver experiment remains NOT_READY.** The quadratic amplitude scaling was
built into the preparation loop and is carried forward by the theorem.

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
and checks the sum of its chart-dependent sector integrals. These quantities
must not be summarized later as physical action transfer between sectors.
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
there is no sampled nonzero plateau in these readouts. Both the quadratic full-action formula and its conservation were analytic
inputs. Their numerical agreement is an integrator, canonical-normalization
and constraint-implementation check, not a discovered property of the matter
model or evidence about action selection.

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

- `INTEGRAL_INVARIANT_IMPLEMENTATION_CHECK`: the seven numerical gates confirm
  implementation of a known Hamiltonian invariant. `numerical_check_passed`
  records their pass/fail status. The original frozen label
  `PASS_CANONICAL_CONTINUUM_CONTROL` is retained only as
  `registered_numerical_verdict`, not as a scientific selection finding.
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

## Review response: theorem and mechanism preflight

The [review](https://github.com/davidmdrpi/geometrodynamics/pull/313#issuecomment-5863520331)
is correct that the primary check is predetermined for an exact smooth
Hamiltonian implementation. Numerical errors or a faulty implementation can
fail the gates, but a pass does not discriminate physical selection mechanisms.
For canonical Theta=p dq and omega=dq wedge dp=-dTheta, with
contraction_X omega=dH, Cartan's identity gives

    Lie_X Theta = d[Theta(X)-H],
    d/dt integral_C(t) Theta = 0.

The extended invariant uses Theta-H dt. For the common proper-time sampling
here, H is the Hamiltonian constraint C, not the coefficient A^2-Q/6 also named
H in the code. On C=0, the lapse-dependent generator NC gives the same
characteristic flow and the closed-loop invariant is unchanged. This is the
Poincaré–Cartan result, not specific to the quartet potential. See Sussman and
Wisdom, [Structure and Interpretation of Classical Mechanics, sections 5.3
and 5.5.1](https://mitp-content-server.mit.edu/books/content/sectbyfn/books_pres_0/9579/sicm_edition_2.zip/chapter005.html).

The review's Liouville concern is relevant but requires a narrower statement.
A regular finite-dimensional Hamiltonian phase space preserves its symplectic
volume. It cannot compress an open basin of nonzero volume into a bounded,
lower-dimensional attractor. A receiver-plus-bulk Hamiltonian truncation is
subject to this constraint after gauge reduction. The full infinite-dimensional
Einstein field theory does not come with a finite Lebesgue phase-volume measure
merely by analogy; the assumptions and any truncation must be stated.

This excludes a full-phase-space dissipative-attractor explanation under those
assumptions. It does not exclude focusing of projected coordinates, finite-time
subsystem plateaus, resonant or isolated periodic solutions, or selected global
boundary-value histories. A projection can lose preparation information while
other coordinates retain it. For an asymptotic receiver-relaxation hypothesis,
one must derive an effective reservoir mechanism and its duration of validity;
an irreversible law must not be silently added to the closed theory. For the
volume theorem, projected focusing and the finite-volume assumptions of
recurrence see David Tong, [Classical Dynamics, section 4.2](https://davidtong.org/teaching/classical-dynamics/dynhtml/S4).

The prior bulk results demand a return-flux and instability check, not an
assumption of an ideal sink: #311 measured approximate tensor refocusing for
specific linear preparations, including a failed low-frequency timing gate.
It did not show that every field component in an interacting nonlinear bulk
returns completely at time pi. #310's scalar growth constrains the experiment's
usable duration; growth is neither an established reservoir nor a universal
proof against a transient receiver response. Spatial compactness alone also
does not supply the bounded finite phase-volume assumptions of recurrence.

Topology and boundary conditions can select classical modes or winding sectors
without quantization. They do not alone establish a universal nonzero unit of
absorbed action. A relation such as I=n hbar cannot be inserted as its own
explanation. Conversely, the review's assertion that every possible boundary
selection of action requires prior quantization is stronger than established:
a concrete classical boundary-value problem would have to be derived and
solved. Nothing in this PR supplies one.

Before constructing the expensive receiver experiment, its freeze must name
which mechanism and which operational selection claim it tests:

| Proposed mechanism | Required derivation and falsifier |
|---|---|
| Effective receiver relaxation | Derive the receiver/bulk split from the retained action; measure outgoing and returning flux, information/energy retained in the bulk, relaxation time, and the instability budget. Reject a sink description when return flux spoils its frozen accuracy over the stated duration. |
| Finite-time resonance or plateau | Specify a physical readout, duration and basin of preparations; test robustness to conjugate initial data and receiver definition. Do not call a narrow resonance or thresholded display an action quantum. |
| Global topological/boundary selection | State the actual classical global conditions and derive the allowed histories and action normalization. Test continuous amplitude deformations within each allowed sector; parity or winding labels alone cannot pass. |

This is a required mechanism preflight, not three mechanisms supplied by this
model. No candidate is established here. If none is specified with a distinct
falsifiable observable, the readiness result stays NOT_READY and another
invariant-only calculation is not a substitute for step 4. The original freeze
is unchanged; these are post-review requirements for a future experiment.

## Reproduction and next decisive test

Raw endpoint states for every finest-grid trajectory, including the RK45
control, are stored losslessly as gzip/base64 JSON in
[`states.json.gz.b64`](../experiments/closure_ledger/runs/20260928_action_selection/states.json.gz.b64).
The original diagnostic report is preserved byte-for-byte as
[`action_initial.json`](../experiments/closure_ledger/runs/20260928_action_selection/action_initial.json).
The current report relabels interpretation and adds mechanism readiness; its
raw states, numerical rows, thresholds and gates are unchanged.
The current diagnostic report, gates and source/archive hashes are in
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
