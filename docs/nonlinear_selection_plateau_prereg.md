# Prospective test: finite-time nonlinear action-readout locking

Baseline: merged #313, 3c8f0d73b27a1076dc550629d8ef14fd09e17736.
Publish this freeze and open a draft PR before implementation or measurements.

## Mechanism and its falsifier

Test whether the existing nonlinear metric/quartet parametric coupling makes
a nonzero scalar action-readout change insensitive to preparation amplitude,
phase and observation duration. This is a finite-time projected plateau
hypothesis, not a full-phase-space attractor or a conserved-loop theorem.
The homogeneous scalar quartet is the responding subsystem; the metric and
scale evolve self-consistently. No bath, dissipation, new field, rounding or
fitted action unit is added. There is no localized receiver or spatial flux.
A positive result is only a candidate plateau in these readouts, not quantized
absorbed action. A negative result rejects this specified mechanism/family,
not every mechanism in classical GR or an unimplemented mouth geometry.

Use nonlinear_supported_tt.py unchanged, kappa=a=1, departure=.15. This is
its expanding homogeneous branch, not the ESU packet background or #312 data.
The exact scalar equation is q''+Omega^2 q=0, with
Omega^2=tr(M^-1)+(r+tr(L^2))/6. Primes mean conformal time. Proper time obeys
d/dt=(1/A)d/deta. Nonlinear variation of Omega can change the projected scalar
readouts, so their response is not fixed by Poincare-Cartan conservation.

## Readouts and balance checks

Fix full S3 volume V3=2pi^2 and reference frequency omega0=2. Define

    J_ref=V3 (q'.q'+4 q.q)/4,
    J_inst=V3 (q'.q'+Omega^2 q.q)/(2 Omega).

Each is the sum of oscillator ellipse areas divided by 2pi at a snapshot,
using respectively the background calibration or the instantaneous quadratic
scalar Hamiltonian with geometry/velocities held fixed. Neither is the action
of a closed nonlinear orbit or an absorbed-action observable. Their difference
is a readout convention to be tested, not a known physical equivalence.
There is no universal action unit in these definitions.

For each preparation use Delta J(t)=J(t)-J(0); do not normalize by amplitude,
input action or per-run maxima before plateau detection. The matched round
background has constant J in both conventions and is checked independently.
The initial constraint correction is subtracted with its own J(0), never
miscounted as a dynamically generated change.

Integrate independent work accumulators along with the inherited evolution:

    dJ_ref/dt=V3 (4-Omega^2) q.q'/(2 A),
    dJ_inst/dt=V3 Omega' [q.q-q'.q'/Omega^2]/(2 A).

Compute Omega' by differentiating tr(M^-1), r and tr(L^2) using the inherited
metric equations, not by differentiating the readout numerically. Check the
accumulators against endpoint changes. This is a matter-equation consistency
check, not a gravitational radiation-flux ledger or a selection verdict.
Retain Hamiltonian/momentum constraints and the six inherited energy terms.

## Frozen preparation and numerical schedule

Let U=diag(1,-1,0)/sqrt(2), W_01=W_10=1/sqrt(2), V=U+W/2.
Complete initial data using
initial_data(cos(theta)U,-sin(theta)V,epsilon,departure=.15,phase=phi).

Amplitudes: zero and epsilon=.02*1.5^k, k=0,...,7.
Shape phases theta=0,pi/4,pi/2,3pi/4.
Quartet phases phi=0,pi/4,pi/2.
Retain all 108 preparations, even if a domain or numerical check fails.
There are three distinct zero-amplitude backgrounds; duplicating them across
shape phases is allowed for complete bookkeeping. No phase averaging may
hide preparation sensitivity in the primary plateau criterion.

Evolve t=0..8.25 with DOP853, rtol=1e-12,atol=1e-14,max_step=.02.
Save t=0 and every .05 within [T-.25,T+.25] for T=1,2,4,8
(45 saved times per history). Measure those late-window means by Simpson integration; compare dt=.05 and .1 using the same
centered windows (the coarse window has endpoints on the .05 grid).
Repeat all largest-amplitude phases and the three round backgrounds using
RK45, rtol=1e-10,atol=1e-12,max_step=.01. No other integrator reruns unless
needed to diagnose a registered numerical failure.

Stop and record any failed solver, nonfinite state, A<=0,H<=0,M not positive,
or Omega^2<=0. Numerical validity requires normalized sampled constraints
<1e-8, det M and symmetry errors <1e-7, and readout/work residual <1e-8 J_bg.
Here J_bg=3 V3/4 is a fixed dimensional reference, not a fitted quantum.
Background readout change must be <1e-9 J_bg. Late-window quadrature changes
and independent-integrator differences must be <1e-6 J_bg. Refinement errors
remain visible even when a selection gate fails.

## Plateau criterion, fixed before seeing results

For each contiguous four-amplitude window (amplitude ratio 3.375), require:

1. Every window mean Delta J at T=2,4,8, for all 12 phase combinations and
   both readouts, exceeds the fixed resolution floor 1e-6 J_bg.
2. Their full range is <=10% of their common median. Do not fit multiple
   plateaus, an integer spacing or a shifted baseline.
3. For each readout/time/phase, every adjacent log-amplitude slope of positive
   Delta J within the window has absolute value <=.1. Points below the floor
   cannot pass by looking flat.

Report all five candidate windows and every criterion, plus per-phase/readout
amplitude ranges, signed raw changes, and logarithmic slopes where resolved.
A common candidate must survive all registered durations and readouts.
Record the T=1 window as transient context, not an additional plateau gate.
This deliberately tests a preparation-robust basin across the listed phases;
it does not exclude smaller phase-specific basins outside that definition.

If numerics pass and no window passes, return
NO_ROBUST_PLATEAU_IN_REGISTERED_FAMILY. If a window passes, return
CANDIDATE_FINITE_TIME_PLATEAU_NOT_QUANTIZATION. If numerics fail or any case is
missing, return INCONCLUSIVE_NUMERICAL_FAILURE. Never replace a failed gate
with a conservation-law pass. A plateau in time alone, especially freezing
under expansion, cannot pass without amplitude and phase robustness.

## Outputs and limits

Publish raw endpoint states and work integrals, all readouts/windows/errors,
source hashes, replay and tests. Preserve the original #312/#313 records.
The selected mechanism is finite-time parametric locking only. Physical pulse
duration, localized receiver definition, spatial return flux and gravitational
recoil remain absent; full receiver action selection is still not established.
#312's negative/phase-dependent active density does not by itself fix a force
or sink. This expanding homogeneous test does not claim to test the ESU's
multi-transit reservoir or its phase-averaged gravitational interaction.
