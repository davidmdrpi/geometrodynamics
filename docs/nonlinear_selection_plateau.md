# Finite-time nonlinear action-readout selection

This experiment tests a concrete projected selection hypothesis: can the
existing nonlinear metric–quartet parametric coupling produce a common,
nonzero change in a scalar action readout across a broad set of preparations?
The criterion is fixed before the run and can fail with a correct integrator.
It is distinct from #313's Poincaré–Cartan implementation check.

Review clarification: the weak-amplitude response has a quadratic analytic
null model. Its observed phase dependence already defeats the all-phase
positive-plateau requirement in that regime. The useful additional numerical
finding is the absence of saturation in the finite-amplitude range tested.
The original freeze omitted this explicit null analysis; the section below
is a retrospective derivation, not an amended preregistration.

## Mechanism and observable

The coupled model is the inherited spatially homogeneous, nonspherical
Einstein–quartet system on its expanding branch, initially A=1.15. The scalar
quartet responds dynamically to the metric and scale. Its equation is

    q''+Omega^2 q=0,
    Omega^2=tr(M^-1)+(r+tr(L^2))/6.

Time-dependent Omega supplies parametric coupling. There is no imposed sink,
friction, apparatus threshold or new field. A finite-time plateau of a
projected observable would not contradict preservation of full Hamiltonian
phase volume. No such plateau is guaranteed by a conservation theorem.

Two snapshot oscillator-action readouts are fixed, with V3=2pi^2:

    J_ref=V3 (q'.q'+4q.q)/4,
    J_inst=V3 (q'.q'+Omega^2 q.q)/(2 Omega).

The first uses background frequency 2; the second uses the instantaneous
quadratic scalar equation with geometry and velocities held fixed. They are
specified canonical readouts in the inherited conformal chart, not a measured
unit of absorbed action, nonlinear closed-orbit actions, or gauge-independent
local detector outputs. The evolution is sampled in proper time, while q'
remains its stored conformal canonical momentum.

The measured response is J(t)-J(0) for each history, including its own initial
constraint correction. This subtraction is essential: the round quartet
already has J_bg=3V3/4=14.8044066016. A nearly unchanged background action is
not a dynamically generated selected value. No amplitude normalization,
integer rounding or fitted action scale enters plateau detection.

Independent work accumulators use the exact matter equations:

    dJ_ref/dt=V3 (4-Omega^2) q.q'/(2A),
    dJ_inst/dt=V3 Omega' [q.q-q'.q'/Omega^2]/(2A).

Omega' is evaluated from the metric equations, and both accumulated changes
are checked against direct readouts. This validates the matter-equation
accounting; it is not a spatial gravitational flux or recoil ledger.

## Frozen test

The [specification](nonlinear_selection_plateau_prereg.md) was published as
`aa7ff3d19d61c8555502eb34cb0021aa13411698` with draft
[PR #314](https://github.com/davidmdrpi/geometrodynamics/pull/314) opened at
2026-09-28 06:06:00 UTC, before implementation and measurements.

The scan includes zero and eight amplitudes .02*1.5^k, four shape
position/velocity phases and three quartet phases: 108 preparations in total.
Every constraint is completed before evolution. Window means are measured
around proper times 1,2,4,8, with the latter three used in the primary gate.
Both action readouts and all phase pairs remain separate in that gate.

A candidate needs four consecutive nonzero amplitude levels (a factor 3.375
range) whose responses all exceed 1e-6 J_bg, stay within a common 10% range,
and have adjacent log-amplitude slopes no larger than .1 in magnitude. All
five candidate windows are reported. This tests a common plateau over the
registered phase family; it does not exclude a smaller phase-specific basin.

The numerical gates concern constraints, chart positivity, the work identities,
round-background constancy, late-window quadrature and an independent RK45
control. They are separate from the physical selection criterion. Numerical
failure makes the mechanism verdict inconclusive, rather than negative.

## Results

The original run is **INCONCLUSIVE_NUMERICAL_FAILURE**: one preparation,
epsilon=.34171875, theta=0, phi=pi/2, had a late-window quadrature difference
of 1.04378327e-6 J_bg in both integrators, just above the 1e-6 limit. All
other numerical gates passed. No case or failed result was removed.

The [targeted refinement freeze](nonlinear_selection_plateau_refinement_prereg.md)
was published as `e1abc93f0accb1e49793fb4792f0cf5780fcb161` before repeating
that preparation in both integrators with .025-spaced window samples. Its
.025 versus .05 difference is 2.70086e-8 J_bg. The other histories retain their
original passing error checks. Original and refined records are separate.

**Refined verdict: NO_ROBUST_PLATEAU_IN_REGISTERED_FAMILY.** All refined
numerical gates pass. All five candidate amplitude windows still fail all
three selection criteria: resolved positive response, a common value and
flat amplitude dependence. Each candidate contains 288 signed readouts
(4 amplitudes x 12 phase pairs x 3 times x 2 conventions); 151–166 fail the
positive/resolved requirement. The largest resolved positive-response slope
in each candidate is about 2.04, 2.14, 2.49, 3.34 and 6.20 respectively,
far above the .1 flatness limit. These are not fitted quantum spacings.

Representative means at t=8, fixing shape phase theta=0:

| Amplitude | Quartet phase | Delta J_ref | Delta J_inst |
|---|---:|---:|---:|
| .02 | 0 | +.00100307 | +.00051492 |
| .02 | pi/2 | -.00101799 | -.00053005 |
| .10125 | 0 | +.02568801 | +.01357508 |
| .10125 | pi/2 | -.02603826 | -.01394107 |
| .34171875 | 0 | +.31141955 | +.21760866 |
| .34171875 | pi/2 | -.31497207 | -.22106264 |

The phase reversal and amplitude dependence falsify the registered robust
locking hypothesis. Late-time slowing of a response on the expanding branch
would not rescue it: a time plateau still retains preparation and readout
dependence. This is a negative result for this specified family and basin,
not a theorem against smaller phase-specific resonances or another GR model.

## Review follow-up: analytic null and what the scan adds

This section responds to the [review of commit 6d24d9a](https://github.com/davidmdrpi/geometrodynamics/pull/314#issuecomment-5864633430).
It uses the existing equations and archived measurements, with no new
trajectories, changed gates or changes to either frozen specification.

Write M=exp(2 beta), beta=epsilon b+O(epsilon^2), with tr(b)=0,
and use primes for conformal-time derivatives. The exact frequency in this
model expands as

    tr(M^-1) = 3 + 2 epsilon^2 tr(b^2) + O(epsilon^3),
    r = 6 - 8 epsilon^2 tr(b^2) + O(epsilon^3),
    tr(L^2) = epsilon^2 tr(b'^2) + O(epsilon^3),
    Omega^2 = 4 + epsilon^2 w2 + O(epsilon^3),
    w2 = (2/3)tr(b^2) + (1/6)tr(b'^2).

There is no first-order scalar-frequency forcing. On the round background,
q_b=(sqrt(3)/2)cos(2 eta+phi)e0 and
q_b'= -sqrt(3)sin(2 eta+phi)e0. Substitution in the exact work identities
therefore gives the two leading changes, including subtraction of the
initial readout:

    Delta J_ref = epsilon^2 C_ref + O(epsilon^3),
    C_ref = (3 V3/8) integral_0^eta w2(s;theta,phi) sin(4s+2phi) ds,
    Delta J_inst = epsilon^2 C_inst + O(epsilon^3),
    C_inst = (3 V3/32) integral_0^eta w2'(s;theta,phi) cos(4s+2phi) ds.

Equivalently, C_inst=C_ref+(3 V3/32)[w2 cos(4 eta+2phi)]_0^eta.
At fixed proper time these coefficients use the background eta(t), and a
window mean averages them over the same proper-time window as the data.
The constraint completion's transverse quartet components start at order
epsilon^2 and do not change these leading quadratic readouts.

The phase dependence needs one qualification to the review's shorthand
epsilon^2 cos(2phi) null. The linear shape equation is

    b'' + (H_b'/H_b)b' + (8+2 Q_b/H_b)b = 0,
    Q_b=(3/4)cos^2(2 eta+phi), H_b=A_b^2-Q_b/6,
    b(0)=cos(theta)U, b'(0)=-sin(theta)V.

Thus w2 itself depends on phi through the responding geometry. A kernel
held independent of phase would yield a linear combination of cos(2phi)
and sin(2phi), with an exact sign reversal under phi -> phi+pi/2; a pure
cosine requires an additional vanishing sine coefficient. Neither exact
antisymmetry nor an exact zero at pi/4 is a symmetry theorem for this
backreacting system. The archive shows an approximate sign reversal, which
is sufficient to fail the registered all-phase gate. Quadratic scaling
alone already excludes a nonzero amplitude plateau wherever the leading
coefficient is nonzero and perturbation theory is accurate. Its accuracy
at a particular finite amplitude is a numerical question.

For theta=0 and the t=8 window, the archived J_ref responses give:

| epsilon | Delta J, phi=0 | Delta J / epsilon^2, phi=0 | Adjacent log-slope, phi=0 | [Delta J(0)+Delta J(pi/2)] / Delta J(0) |
|---:|---:|---:|---:|---:|
| .02000000 | .00100307 | 2.507671 | — | -.014881 |
| .03000000 | .00225673 | 2.507480 | 1.999812 | -.014804 |
| .04500000 | .00507685 | 2.507089 | 1.999615 | -.014636 |
| .06750000 | .01141978 | 2.506398 | 1.999320 | -.014288 |
| .10125000 | .02568801 | 2.505765 | 1.999377 | -.013635 |
| .15187500 | .05786771 | 2.508786 | 2.002972 | -.012689 |
| .22781250 | .13161949 | 2.536093 | 2.026699 | -.012135 |
| .34171875 | .31141955 | 2.666909 | 2.124044 | -.011407 |

These are post hoc diagnostics, not new selection gates or a fit of an
action unit. Values come from `plateau.json`'s `cases[*].diagnostics.window_means[3][0]`,
with the DOP853 replacement from `refinement.json` for the final pi/2 row.
The slope is log[J(epsilon_i)/J(epsilon_(i-1))]/log(1.5).
The nearly constant coefficient confirms approximate quadratic behavior;
calling it exact epsilon^2 scaling would overstate the data. The last two
slopes steepen to 2.027 and 2.124 rather than approach the frozen .1 limit.
This finite-amplitude observation is the added information beyond the
weak-amplitude null. It does not establish behavior beyond epsilon=.34171875.

Any further plateau experiment must state and test its analytic response
null prospectively, and explain what mechanism could depart from it in the
chosen regime. This retrospective correction does not repair the omission
in the original freeze or convert these data into evidence for selection.

| Numerical check | Observed bound |
|---|---:|
| Original DOP853 sampled constraint residual | 7.20e-14 |
| Original DOP853 readout/work error / J_bg | 4.77e-15 |
| Original DOP853 determinant error | 7.78e-15 |
| Refined targeted quadrature difference / J_bg | 2.71e-8 |
| Largest remaining original quadrature difference / J_bg | 8.91e-7 |
| Largest independent-integrator mean difference / J_bg | 1.39e-12 |

The precise numerical checks, including RK45 and all chart minima, are in
the records. Bounds are sampled diagnostics, not rigorous continuum estimates.
Window refinement changes quadrature sampling, with a fresh solve for the two
targeted histories; it does not introduce a new spatial resolution.

## Limits and reproduction

This probes finite-time parametric locking within an existing homogeneous
sector. It has no localized emitter, mouth, independent spatial receiver,
physical preparation pulse, spatial return flux or detector model. Evolution
duration is not preparation-pulse duration. The experiment neither measures
quantized absorption nor tests the complete antipodal wormhole proposal.
#312's seed-overlap limitation is not solved by replacing it with this homogeneous
model; this is a separately labeled mechanism test. No claim about an ESU
reservoir or phase-averaged gravitational attraction follows from it.

The archived states and two work accumulators are lossless gzip/base64 JSON;
the report contains every signed time series, window mean, candidate gate,
error estimate and source hash. Both are in
`experiments/closure_ledger/runs/20260928_selection_plateau/`.

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.selection_plateau_probe
python -m pytest -q tests/test_selection_plateau.py
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.selection_plateau_probe \
  --replay experiments/closure_ledger/runs/20260928_selection_plateau/states.json.gz.b64 \
  --output-dir /tmp/plateau-replay
```

Exit 0 means the numerical evidence supports a candidate or a registered
negative result; it does not mean selection occurred. Exit 1 means the
numerical result is inconclusive. Read `selection_verdict`, not only the
process exit. The localized-receiver verdict is always `NOT_TESTED` here.

Regression controls accept a synthetic nonzero plateau and reject quadratic
scaling, near-zero flatness, phase-dependent values, missing cases and solver
failures. They also independently differentiate the readouts to check the
work-rate formulas. None substitutes for the measured selection criterion.

Replay the separately frozen refinement without overwriting the original run:

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.selection_plateau_refinement \
  --replay experiments/closure_ledger/runs/20260928_selection_plateau/refinement_states.json \
  --output-dir /tmp/plateau-refinement-replay
```

The next scientific choice should change or derive a mechanism, not re-label
this failed plateau test as selection. A localized receiver experiment would
still need its own action readout and spatial exchange accounting.

Two other questions raised by review remain open, with prerequisites:

- **Discrete events:** for n=phi/|phi|, an integer degree needs a closed
  oriented spatial domain mapped to S3, or a receiver region with boundary
  conditions that define such a degree (for example, a constant boundary
  map collapsed to a point). A general finite receiver-window integral is
  not integer-valued. Boundary transport must be distinguished from a
  change through |phi|=0, where n is undefined. Even a well-defined integer
  event count does not make the action per event universal; that is a
  separate amplitude, pulse-duration and receiver-definition test. The
  S3/RP3 identification must also be specified for the proposed field map.
- **Global history selection:** a boundary-value experiment needs explicit
  physical boundary conditions and a demonstrated discrete solution set.
  Periodicity alone does not establish discrete action values or remove
  continuous families and free scales.

This PR supplies neither calculation. Its homogeneous null also cannot
establish a general no-go theorem for localized self-gravitating receivers.
Claims about soliton existence, event stability or universal continuous
absorption require analysis of those configurations and their boundary
conditions. They are not consequences of this scan.
