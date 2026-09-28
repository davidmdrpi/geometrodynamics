# Step 4: prospective dynamical action-selection audit

Baseline: merged #312, f193bd1a9b6eb12cd9eb75d2b23a973b5ef688f1.
Publish this specification and open a draft PR before new implementations or
runs. The analytic predictions below are disclosed prospective design work.

## Intended question and readiness

Step 4 asks whether transferred action, swept over preparation amplitude and
duration, approaches robust dynamically selected values surviving changes in
preparation, numerical resolution and receiver definition. A spatial eigenmode
index, phase winding, squared Weyl norm or normalized packet fraction is not
an absorbed action measurement. Do not insert hbar, round to integers, fit an
action unit, or normalize every preparation to unit action.

The merged repository has linear finite-packet propagation (#311) and
nonlinear two-region initial data/first time jets (#312), but no evolved
localized receiver, separated-source causal transfer or gravitational recoil
ledger. Its full receiver-selection verdict must remain NOT_READY until those
prerequisites and an operational absorbed-action observable are supplied.
This experiment is an executable nonlinear canonical-action control, not a
substitute receiver, a new matter model or a completion of the full step 4.

## Available exact nonlinear dynamics

Use the existing homogeneous, nonspherical four-field Einstein reduction in
geometrodynamics/waves/nonlinear_supported_tt.py, without changing its action,
constraint completion or evolution. This model has no distinct spatial source
and receiver. Variables are A,A',q,q',M,L=M^-1 M'/2 with det M=1 and
H=A^2-Q/6, Q=q.q, kappa=a=1. Use departure=.15 and scalar phase=0.
Its canonical one-form per unit S3 volume is

    Theta = -6 A' dA + q'.dq + tr(Pi dM), Pi=(H/2)L M^-1.

For a CLOSED LOOP of preparations C, define I=Vol(S3)/(2pi) integral_C Theta,
Vol(S3)=2pi^2. This is canonical circulation over a family of solutions, NOT
the action of a closed orbit of one solution or action absorbed by a detector.
The distinction is mandatory in every verdict. Separately record signed
scale, scalar and shape contributions in this fixed conformal canonical chart.
Their redistribution is chart-dependent and is not regional momentum recoil.

Let U=diag(1,-1,0)/sqrt(2), W_01=W_10=1/sqrt(2), other W entries zero,
V=U+W/2. For loop parameter theta use the inherited constraint completion

    initial_data(cos(theta)*U, -sin(theta)*V, epsilon,
                 departure=.15, phase=0).

U,V are STF, noncommuting, tr(U V)=1. The scalar momentum response is retained.
At initial time A,q are theta-independent and H0=1.15^2-1/8=1.1975. Thus

    integral_C Theta = pi H0 epsilon^2,
    I0 = pi^2 H0 epsilon^2.

This exact continuous initial family is an analytic null prediction, not a
fitted scale. Smooth Hamiltonian evolution preserves the full closed-loop
circulation. It need not preserve each sector contribution. Conservation is
not a theorem excluding discrete receiver responses or special globally
selected histories outside this homogeneous Cauchy-data family.

## Frozen schedules and checks

Amplitudes epsilon=0,.01,.02,.04,.08. Evolve in inherited proper time t/a,
record t=0,.25,.5,1,2. This varies observation/evolution duration, NOT a
physical source pulse duration; that missing preparation control remains in
the readiness checklist. Phase-loop grids N=16,32,64, endpoint excluded.
Use Fourier differentiation in theta and trapezoidal loop integration.
Keep all scalar/scale variables and constraints. No damping or projection.
DOP853 inherited tolerance (rtol=1e-12,atol=1e-14,max_step=.02).
Reuse nested theta trajectories across N, noting this is quadrature refinement,
not an independent time integrator. At epsilon=.08, N=64, repeat with RK45,
rtol=1e-10,atol=1e-12,max_step=.01. Reject domain violations and failed solves.

Record normalized constraint maximum and minima of A,H,eigenvalue(M); det M
and M-self-adjointness of L residuals. At finest N require constraints <1e-9,
det/symmetry errors <1e-8, positive domain, initial I agreeing with the exact
prediction within 1e-8 relative, and maximum full-circulation drift <1e-6
relative to I0. Zero-amplitude action must be <1e-12 absolute. Sector loop
quadrature differences N32->64 and RK45->DOP853 must be <1e-5 I0 (for nonzero
amplitudes). Publish N16->32 differences without an extra convergence gate.

Report I_shape, I_scalar, I_scale and their changes at each duration; also
I_shape+I_scalar and full I. Compare shape-only apparent losses/gains to the
full ledger. Test epsilon-scaling using raw actions, with no per-run amplitude
normalization in the primary data. Report adjacent log-slopes for each sector
only where both values are safely nonzero; do not turn near-zero projections
into evidence of selection. No positive selection verdict is available from
these proxy readouts. Amplitude doubling predicts a factor 4 in full I, not
a nonzero plateau. Numerical agreement verifies the continuum control.

## Required outputs and interpretation

Publish source hashes, reproducible runner, raw finest-loop endpoint states,
all amplitudes/durations, quadrature and integrator controls, tests and report.
Preserve failures. Distinguish numerical validity, the continuum-action finding,
and NOT_READY_FOR_RECEIVER_ACTION_SELECTION. The latter requires actual
preparation-duration sweeps and receiver-definition controls as well as evolved,
constraint-monitored localized transfer; changing homogeneous sector readouts
does not meet that prerequisite. #312's overlapping-profile review limitations
and its failed coarse-grid gate remain unchanged. No new quantum milestone.
