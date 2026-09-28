# Finite-time nonlinear action-readout selection

This experiment tests a concrete projected selection hypothesis: can the
existing nonlinear metric–quartet parametric coupling produce a common,
nonzero change in a scalar action readout across a broad set of preparations?
The criterion is fixed before the run and can fail with a correct integrator.
It is distinct from #313's Poincaré–Cartan implementation check.

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
