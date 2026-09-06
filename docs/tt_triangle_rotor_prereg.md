# Preregistration: does the existing TT field contain the triangle rotor?

Publish this file before implementing or running the new experiment.
Base: PR #286 at cd7f300, including #284 and #285.
This follows the user's request to derive the triangle rotor from an existing
BAM field mode. It does not presuppose that such a reduction succeeds.

## Existing model and scope

Use the homogeneous, trace-free symmetric tensor beta in waves/backreaction.py.
Its five real components obey beta_ddot + (8/a^2) beta = S, where S is the
appropriately normalized projected matter stress. This is the repository's
linear TT response on a fixed 3+1 ESU with spatial S3, not its nonlinear 4+1
spherical Tangherlini model and not a completed resolved mouth solution.
Keep that distinction explicit. The existing shear_sources routine integrates
stress rather than dividing by the spatial volume; track physical source
normalization instead of silently identifying its output with kappa times
an average.

The candidate family is
    beta = A (n n^T - I/3), |n|=1, A != 0.
A is a tensor amplitude, not the analyzer setting. n and -n represent the
same field: the angular orbit is RP2. A=0 is a singular apex at which the
director cannot be recovered. An orientation is relational to a physical
source/apparatus frame, not an absolute coordinate-frame observable.

This family is a proposed embedding, not an already validated field reduction.
No new matrix field, shear modulus, rigid-shape constraint, fitted potential,
quantum probability, Born rule or canonical equilibrium is to be added to
make it succeed.

## Design analysis already known before freezing

The inherited frequency is omega^2=8/a^2. Substituting the candidate into the
Frobenius kinetic term gives
    tr(beta_dot^2) = (2/3) A_dot^2 + 2 A^2 |n_dot|^2.
Consequently a restricted action has a round angular metric at fixed A.
This identity alone is not a consistent truncation of the full field equations.

Analytic inspection suggests a normal obstruction: with P=I-n n^T, the
trace-free part within the two-dimensional transverse plane is
    N_n(M) = P M P - tr(P M P) P/2.
For v=n_dot tangent to n, the omitted part of the field equation is expected
to be
    N_n(beta_ddot + omega^2 beta - S)
      = 2 A [v v^T - |v|^2 P/2] - N_n(S).
Its first term has Frobenius norm sqrt(2) |A| |v|^2.
This is a pre-calculation analytic prediction, not an unanticipated numerical
discovery. It predicts failure of an autonomous rotating uniaxial family in
the source-free TT model. It does not exclude a specially constrained or
driven solution, other tensor shapes, other field modes or all BAM completions.

## Q1: action, normalization and field identity

Derive the quadratic tensor action from the homogeneous Einstein-Hilbert ADM
action, including the spatial volume and gravitational coupling. Use lapse
one, fixed ESU radius a, trace-free beta and spatial metric a^2 exp(2 beta).
Check the curvature expansion independently of the remembered oscillator
frequency. Obtain the source term from variation of the matter action.
A homogeneous isotropic supporting stress contributes no anisotropic source
in this approximation. Track the limits of this frozen-background expansion.

Report the restricted amplitude/angular kinetic metric and potential, and the
actual field-source coupling. Distinguish an overall action normalization
from the dynamical coefficients: the response equation alone cannot fix
Liouville normalization. Verify the frequency against the repository's
Cartan/Einstein derivation and its tensor constant.

Use two independent checks of the pullback metric: matrix velocities and
finite differences of the field embedding. Verify covariance under a common
SO(3) change of tensor/director/source frame. Check the antipodal field
identity and the amplitude-zero singularity.

## Q2: full equations, not just projected equations

Derive both Euler-Lagrange equations of the restricted action, allowing A
to vary. Evaluate their proposed motion in ALL five tensor equations.
Derive the transverse normal projector and test it against explicit
orthonormal basis components. Report tangential, amplitude and normal
residuals separately.

For S=0 and A!=0 predict that nonzero v generates a nonzero normal residual
and exits the axisymmetric cone. Constant director with oscillating A is
the invariant control, away from A=0. Treat frozen A separately: it requires
a radial reaction even before the normal obstruction is considered.
A projected trajectory or a fitted constraining stress is a control, never
a successful derivation from the unconstrained model.

For general symmetric trace-free S, report the necessary normal-stress
condition. Exhibit a manufactured stress that makes an assigned constant-A,
constant-speed great-circle trajectory solve the full equation, and label
that stress as added control data. Its existence prevents an overbroad
claim that a rotating axisymmetric field is impossible with every source.
Do not claim that the repository's computed scalar sources supply that
feedback unless that identification is actually derived.

## Frozen numerical checks

Use dimensionless radius a=1 and kappa=1 for reported numerical action
coefficients; report their analytic scaling. The tensor perturbations must
remain small.

Primary source-free data:
    A0=0.01, A_dot0=0, n0=e_z, v0=0.4 e_x.
Initialize beta and beta_dot from the candidate embedding.
Propagate the unconstrained tensor exactly as a harmonic oscillator and
independently with DOP853 at rtol=1e-11, atol=1e-13.

Measure distance to the FULL uniaxial cone, minimizing over director and
signed amplitude. Do not compare only with a frozen director or a single
projected trajectory. Use eigenvalue minimization and verify it against a
separate rotational invariant/eigenvalue formula. Evaluate t=0.04, 0.02,
0.01, 0.005. Predict
    distance = |A0| |v0|^2 t^2 / sqrt(2) + O(t^4)
for these A_dot0=0 data. Require step-halving ratios between 3.8 and 4.2 and
the final scaled coefficient within 1% of the prediction.
The control v0=0 must remain uniaxial to 1e-12 absolute error before A=0.

Repeat the primary motion with A0=0.005 and 0.0025, scaling beta_dot with A.
The relative distance must remain unchanged to 1e-9. This tests whether
leakage is leading order in field amplitude, rather than a discarded
nonlinear correction. Also report speed dependence for v0 magnitudes
0.2, 0.4, 0.8; no nonzero-speed result may be relabelled invariant because
the absolute perturbation is small.

For algebraic tensor/projector/metric and covariance checks use seed 20260907
with 32 random unit axes, tangent velocities, signed nonzero amplitudes and
symmetric trace-free sources. Target scaled residual 1e-10 or better.
Finite-difference metric checks target 1e-7; exact-flow/ODE agreement 1e-10.
Verify conserved free energy and rotational Noether charge [beta,beta_dot].
For the manufactured rotating solution require full equation residual below
1e-10, and explicitly report its nonzero radial and normal forcing.

## Q3: consequence for selection and readout

If Q2 obstructs the rotor, do not assign a canonical rotor ensemble or infer
Phi. The induced metric on a constrained submanifold is not the probability
law of the unconstrained five-component field. Discuss remaining shape,
amplitude and source-preparation variables. Gauge-invariant TT at linear
order does not turn an arbitrary global director into a local operational
readout without a relational apparatus.

Compute the field's existing anisotropic coupling from the matter variation.
On the candidate family its angular dependence should be quadrupolar,
A n^T S n. Distinguish that from a derived analyzer-itinerary frame-restoring
energy. A local tensor contraction along a physical apparatus direction m
is even under n -> -n; this does not by itself protect source readability
because informative even observables already exist. No complete BAM pointer
or future-setting intervention is claimed by this calculation.

Success means a justified field reduction OR a verified scoped obstruction.
A numerical failure without a certificate/analytic argument is unresolved,
not a no-go. Use separate fields for experimental check status and physical
reduction verdict; the probe must exit nonzero when a required check fails,
even if the expected physical verdict is an obstruction.

## Deliverables and provenance

New reusable module, independent tests, archived values and criterion statuses,
derivation with explicit assumptions, status-audit update and a draft PR.
Keep existing field solvers, round 9, pointer-spread archives and all historical
freezes unchanged. Record any implementation correction or numerical redesign.
There is no preselected CHSH value and no attempt to adjust a coupling to
produce the cubic counting function.
