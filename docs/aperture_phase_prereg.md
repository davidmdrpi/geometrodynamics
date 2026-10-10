# Preregistration: controlled Berger phase-scaling family

Date: 2026-10-10. This specification, implementation and producer are to be
published together before any scheduled propagation. It follows the accepted
review of #323, including its portable source validation. No new propagation
in this parameter box has been used to choose the schedule or thresholds.

## Question and two falsifiable claims

At fixed aperture/wavelength, pulse cycles and dimensionless port strength,
does finite-lead capture retain a substantial fraction of the round reference
as deformation-induced phase spread increases? Does integrated capture obey
the often-suggested inverse-phase scaling over this declared finite family?

These are different claims and receive separate labels. Inverse-phase decay
can pass while retention fails; a plateau can pass retention and fail the law.
The inverse-phase law is a phenomenological hypothesis motivated by the
stationary-phase discussion in #323, **not a derived theorem for integrated
lead capture**. A peak-intensity estimate does not imply that hypothesis.
No inference to physical throat scales or infinite-frequency asymptotics is
registered. The largest carrier is 48 and the largest proxy is about 42.4 rad.

## Geometry, wave action and exact reduction

Use the same conformal scalar, ultrastatic Berger S3 action and two coherent
antipodal distributed ports as #323, with R=1/pi. There is no evolving GR
background, feedback reinjection, excised mouth, motion or corrective kick.
The eigenfrequency squares are

    k_lm = pi^2 [l(l+2) + (b^-2 - 1)m^2 + (4-b^2)/3].

The compact aperture profile is (1-(theta/a)^2)^4 for theta<a, physically
L2-normalized. Character quadrature gives coefficients A_l, with the omitted
squared norm retained as a diagnostic rather than renormalized away.

At exact antipodes, W_B=(-1)^l W_A. Dark Gram rows never get excited from zero
initial field; +/-m have equal frequency and parallel port vectors. Keep one
row for m>=0, with W_A=sqrt(gamma A_l^2 d_m/(l+1)), d_0=1 and d_m=2 otherwise.
This is an exact degeneracy reduction of #323, not a ray approximation.

    q''+Kq+WW^T q'=2W a_in;  b_out=W^T q'-a_in
    E=(|q'|^2+q^T Kq)/2;  dE/dt=|a_in|^2-|b_out|^2.

All initial states vanish. Source A is driven and receiving B has no incoming
wave. Capture means outgoing B lead energy divided by complete incident A
energy, including any prompt reflection. Record reflected and remaining bulk
energy as well. No flux is counted twice or discarded as an assumed loss.

## Controlled family and coordinate norm

Primary footprint f=a*w=4.8. Squashes b=(.8,.9,.98,1,1.02,1.1,1.2), carriers
w=(12,24,48). The secondary transducer family f=7.2 uses b=(.8,1) at all three
carriers. Round baselines are separate at every footprint/carrier.

- a=f/w; coordinate wavelength lambda/R=2*pi/w, so a/lambda=f/(2*pi).
- gamma=8*w/12; gamma/(pi*w)=2/(3*pi) is fixed.
- Source halfwidth h=3/w; source cos(pi*t/(2h))^4 cos(pi*w*t) on |t|<h.
  This holds three carrier cycles over the full source duration and fixes
  relative spectral bandwidth. Its analytic amplitude bound is one.
- Coordinate cap size is defined using the round angle theta. These are not
  equal physical geodesic-radius mouths as b changes.
- Every case uses t in [-h,1.75]. The fixed horizon is in round-transit units,
  not fitted to an outcome. The maximal metric stretch is 1.2, so this covers
  a first transit plus the longest halfwidth .25 and a .30 margin. This is a
  declared observation window, not a proof that all capture has finished.
  An extension to 2.25 at the most strongly squashed, highest-frequency primary
  case reports the later increment separately without changing any main gate.

For scale s=w/12, coarse L=40s, dt=1/(512s^2); fine L=56s,
dt=1/(1024s^2). Resolving midpoint phase drift motivates dt falling faster
than 1/w: its leading single-mode phase error over fixed time is O(omega^3 dt^2).
Keeping L/w fixed resolves the shrinking aperture. Paired convergence is
checked everywhere. At b=.8,w=48,f=4.8 add separate time refinement
(dt=1/32768), mode refinement (L=288), and the extended horizon at fine dt.
There are 54 paired trajectories plus three controls, exactly 57 cases in
producer order. No adaptive expansion, repeated selection or tuned horizon.

Waveform errors use ||out_x-interp(out_y)||_2 / ||incoming_x||_2 over both
ports and the full common sampled window. Interpolation is used for comparison
only. Final-state reconstruction uses unweighted Euclidean ||delta q||_2 +
||delta q'||_2. No fitted phase, amplitude rescaling or peak alignment is allowed.

## Exact phase diagnostics and specified fit

At one round transit,

    delta_phi_lm = sqrt(k_lm) - pi*(l+1).
    Phi = pi*abs(b^-2-1)*w/2.

Phi is a proxy, not the exact phase of the packet. Report the mean, standard
deviation and squared coherence of exp(i delta_phi) using fixed weights
W_A^2 |F(pi*(l+1))|^2, normalized to sum one, where F is the analytic Fourier
integral of the compact source. These are the free round-source modal energy
weights; they are not inferred from the receiving waveform and are not a
replacement for the coupled dynamics. The transform follows by expanding
cos^4 into its constant, first and second harmonics. It is tested against
independent quadrature.

Fit separately for f=4.8 and f=7.2, **only on the preselected b=.8 slice**:
log(capture_nonround/capture_round)=intercept-beta*log(Phi), ordinary unweighted
least squares on w=(12,24,48). Here Phi=(10.60,21.21,42.41), a factor-four
interval. These three deterministic points are not a statistical sample;
no population confidence interval or asymptotic exponent is claimed.
Other b values are measured cross-checks, not extra fit points or grounds to
select a better exponent after seeing results.

## Gates and failure rules

All 57 cases must have finite arrays, exact scheduled grid/configuration,
zero initial energy and nonnegative energy to 1e-12. Source formula comparison
uses absolute tolerance 1e-14, analytic amplitude one; inactive port and zero
outside compact support must be exact. Physics reconstruction uses the
archived source, not an overwritten locally evaluated source.

Numerical validity requires every applicable gate:

- field-plus-lead ledger error/input <1e-9; final energy consistency <1e-9;
- independently reconstructed port samples <1e-8 absolute; reconstructed
  energy/input <1e-8 and final-state norm <1e-8;
- aperture omitted squared norm <1e-4;
- premature B energy/input <1e-5 before -h+min(1,b)*(pi-2a)/pi;
- all paired and separate waveform comparisons <.03;
- coarse/fine fitted beta differs by <.1 for both footprints.

If any numerical gate fails, both claims are NUMERICALLY_UNRESOLVED. A
measurement with failed convergence is not physical evidence for either label.

The phase-coverage gate requires, for each footprint at b=.8, weighted exact
phase standard deviation >=5 rad at w=48 and its w=48/w=12 ratio >=3. If this
fails despite numerical validity, both labels are PHASE_RANGE_INSUFFICIENT.
This criterion prevents a large proxy alone from being called large sampled
phase spread. Passing it still does not establish an asymptotic limit.

With numerical validity and phase coverage established:

1. **Retention hypothesis:** every fine nonround case with Phi>=10 must retain
   >=.1 of its matching round capture. Any such case below .1 gives
   FAILED_IN_DECLARED_FAMILY; otherwise SUPPORTED_IN_DECLARED_FAMILY.
2. **Inverse-phase hypothesis:** both preselected footprint fits must have
   beta in [.5,1.5] and maximum absolute log residual <=.25. Any fit outside
   either bound gives FAILED_IN_DECLARED_FAMILY; otherwise
   SUPPORTED_IN_DECLARED_FAMILY. This tolerance tests a broad approximate
   inverse law, not exact beta=1. A failure is not relabelled as evidence for
   a different law selected after measurement.

There is no shifted-return success gate, disconnected-port mechanism claim,
or action-selection claim in this study. Those do not measure the proposed
scaling. R3 propagation and closed-feedback history remain untested.

## Fine tuning, budget and provenance

Homogeneous Berger geometry, conformal scalar coupling, compact polynomial
coherent ports, selected f values, gamma/w and pulse shape are explicit
assumptions. Scaling those quantities together is a controlled mathematical
family, not a derived physical evolution of mouths. It removes some
confounding from #323 but not all dependence on the transducer or geometry.
The corrective-kick budget is exactly zero; no kicks, timing fit or output
renormalization are permitted. There is no time-dependent metric work in
this static family. The input/lead/bulk energy ledger is mandatory.

Known information before freeze: #323 captures and retention, the review's
phase proxy, and its observation of coherent cancellation for the larger
aperture. Tests before freeze use algebra, cap quadrature, synthetic scorer
inputs, analytic source transforms and tiny evolution examples outside the
production schedule. They do not measure any scheduled capture or exponent.

Production refuses to overwrite an existing run directory, binds all source
hashes to the published freeze, records run-start/run-finish UTC and per-case
simulation timestamps, and archives every time grid, input/output history,
energy history and final modal state. NPZ timestamps are not used as run
evidence. A SHA-256 manifest authenticates all files; a later replay receives
its independently published pinned hash, reconstructs the field from archived
port forces, remeasures capture, and reproduces the two declared labels. No
new propagation is required for replay. Failed or unresolved outcomes are
published with the same evidence requirements as successes.
