# Pre-registration: nonlinear completion of the supported tensor sector

Baseline: merged main `934dabc1ede2d8e84c7585e6778c0dc201243faf` (#298).
Baseline tree: `5c6d0054b882e7faad0785818e91018a357cf9f6`.
Branch: `codex/nonlinear-supported-tt-prereg`. Seed: `2026091319`.

This is the proposed milestone experiment: **do the linear tensor histories
of #295–#298 belong to nearby solutions of the full Einstein–matter equations?**
Publish this document before implementing, symbolically checking, solving
constraints, or integrating the new nonlinear model. Formulas and likely
outcomes below are scoping predictions, not results of an experiment.
Preserve this freeze if any prediction fails; publish the correction and
its consequences separately. The milestone number must not affect verdicts.

## 1. Prior state and separate questions

The repository establishes four-component odd scalar support (#294), its
linear homogeneous scalar instability (#296), the fully supported linear
TT equation on an evolving FRW family (#297), and an all-phase turning
certificate with complete linear asymptotic transport (#298). These results
do not establish a finite-amplitude Einstein–matter solution with tensors.

There are four distinct targets:

- D: an exact homogeneous reduction of the full field equations, including
  all four scalar equations and Hamiltonian/momentum constraints;
- C: nonlinear constraint completion of specified linear tensor data;
- N: finite-amplitude evolution, nonlinear response, and recovery of #297
  as the tangent of the completed family;
- F: a controlled future-persistence statement for a stated neighborhood
  of admissible data, if a suitable bound can be proved.

The prior is positive for D and for a nontrivial family under C and N.
Completing arbitrary noncommuting tensor displacement and velocity may
require a second-order matter response. The future prior is persistence
for sufficiently small data on an expanding branch separated from the ESU,
with positive effective gravitational coupling and asymptotically frozen
shape. This is a conjecture to establish or leave unresolved, not a license
to call long integrations a global existence theorem.

A positive D/C/N result is already a nonlinear existence result in a
specified invariant sector. F is a stronger, separate claim. A failure of
a restricted matter ansatz is not a no-go for the four-field theory, still
less for every BAM matter sector. A sampled success is not an open-set proof.

## 2. Matter, metric and gauge budget

Use exactly the four real conformal scalars and action of #294–#298:

    I = integral sqrt(-g) [(R-2 Lambda)/(2 kappa)
        - (1/2) sum_I ((grad phi_I)^2 + R phi_I^2/6)],
    Lambda = 3/(2 a^2), kappa>0, a>0.

No extra field, perfect fluid, prescribed stress, damping, drive, or
constitutive relation is added. The original four fields may respond.
There is no assumption of the FRW radiation law away from isotropy.

Use the inherited invariant coframe e^i on the unit S3 and a general
positive symmetric shape matrix M of determinant one:

    ds^2 = A^2 [-n^2 d eta^2
               + M_ij (e^i + N^i d eta)(e^j + N^j d eta)],
    M = exp(2 beta), beta^T=beta, tr beta=0.

Retain lapse n and all three homogeneous shifts N^i until their equations
have been varied or independently recovered from the full Einstein tensor.
Only then impose n=1, N^i=0. Proper time satisfies dt=A d eta in that gauge.
Time-dependent diagonalization is not free: its shift and the corresponding
matter transformation must be retained. Do not eliminate momentum equations
by choosing a diagonal metric before checking them.

Let S_i be the real 4x4 quaternion derivative matrices already defined in
`geometrodynamics/waves/reciprocal_scalar_tt.py`. The primary candidate is

    B(q) = q_0 I_4 + sum_{i=1}^3 q_i S_i,
    phi_I(eta,x) = [B(q(eta)) x]_I / A(eta).

These are four real fields, not four new fields on top of the old four.
The inherited background is q=(q_b,0,0,0). Every component remains odd
under x -> -x. Allowing q_i to evolve changes the existing support's
configuration. It is not an assumed mechanism selecting that configuration.

**P1 (closure prediction):** the quaternion identities imply
B^T B=|q|^2 I, and summed stress is homogeneous in the invariant frame.
The scalar equations preserve this four-dimensional coefficient space for
a homogeneous metric. Verify the left/right action and shift signs using
the repository's actual S_i, not an interchangeable quaternion convention.

The more restrictive common-amplitude ansatz B=q_b I is a control, not the
only admissible support. A general real 4x4 coefficient matrix is a diagnostic
extension if quaternionic closure fails: its spatially varying stress cannot
simply be projected back onto a homogeneous metric. Report the generated
spatial multiplets if the proposed reduction is not exact.

## 3. P2: exact action and independent field-equation derivation

At zero shift and n=1, put

    Q=|q|^2, F=1/kappa-Q/(6 A^2), H=A^2/kappa-Q/6,
    L=(1/2) M^{-1} M', tr L=0,
    r(M)=2 [2 tr(M^{-1})-tr(M^2)], det M=1.

Here H is a kinetic coefficient, not the Hubble rate. Use a different
symbol/name for the latter in code. The predicted reduced Lagrangian per
unit-S3 volume, after the appropriate curvature boundary term, is

    L_red = -3 A'^2/kappa + |q'|^2/2
            + (H/2) tr(L^2) + (H/2) r(M)
            - (Q/2) tr(M^{-1}) - Lambda A^4/kappa.

With lapse n and zero shift, the first three kinetic terms carry 1/n and
the last three potential terms carry n. Derive the covariant velocities
with shift before extracting the three momentum constraints. This displayed
gauge-fixed expression is a prediction for the evolution action; by itself
it does not supply the missing shift equations.

Derive all of the following separately:

1. The exact spatial curvature, extrinsic-curvature terms, scalar gradient
   and kinetic terms, including nonminimal curvature coupling and its
   boundary contribution. Do not use a quadratic expansion in beta.
2. The lapse and shift constraints and their propagation identities.
3. The evolution equations for A, q and all independent shape variables.
4. Full coordinate curvature, every improved scalar stress, and all four
   Klein–Gordon equations evaluated independently of the reduced equations.

The two routes must agree off shell. The full metric continuation must be
positive definite at finite beta; the earlier I+2 epsilon beta device was
only a first-variation continuation and is insufficient here.

The isotropic limit must recover q_b''+4q_b=0 and

    A'^2 + A^2 - Lambda A^4/3 - kappa (q_b'^2+4q_b^2)/6 = 0.

The linear TT tangent must recover exactly #297:

    (H_b beta')' + K_b beta = 0,
    H_b=A_b^2/kappa-q_b^2/6,
    K_b=8 A_b^2/kappa+2 q_b^2/3.

Retain physical versus normalized canonical momenta and the fixed volume
2 pi^2. No reuse of #298's normalized pair on an enlarged phase space without
deriving its normalization and canonical embedding is allowed.

## 4. P3: nonlinear constraints and initial-data completion

Set eta=t=0 initially. The reference expanding data are

    A_b=a(1+d), A_b'=(A_b^2-a^2)/(sqrt(2) a),
    q_b=a sqrt(3/(4 kappa)) cos(delta),
    q_b'=-2 a sqrt(3/(4 kappa)) sin(delta).

For fixed STF matrices U,V prescribe

    M_0=exp(2 epsilon U),
    M_0'=2 epsilon M_0^{1/2} V M_0^{1/2}.

This is a symmetric tangent with tr(M_0^{-1} M_0')=0 and the desired linear
beta displacement/velocity (U,V). It is not claimed that beta'=epsilon V
exactly at finite amplitude when U and V do not commute.

Hold A_0, q_0=q_b and q_0'=q_b' fixed, and first try the following explicit
completion prescription. With D_b=q_b^2+q_b'^2>0, introduce three unknowns
xi_i and set

    q_i = -q_b' xi_i / D_b,
    q_i'=  q_b  xi_i / D_b, i=1,2,3.

Solve all four constraints for (A_0',xi_1,xi_2,xi_3), choosing the expanding
solution continuous from (A_b',0,0,0). The denominator is nonzero at both
field and field-velocity zeros. It uses only variables of the original four
fields; no division by q_b or q_b' alone is permitted.

**P3a:** near epsilon=0, the initial-data constraint Jacobian with respect
to these four unknowns is nonsingular for d>0. Establish its rank and its
phase dependence analytically. If true, the implicit function theorem
supplies local families for each specified U,V and phase. State precisely
which parameters range over an open neighborhood in the constraint surface
and which are fixed preparation choices; do not infer openness in the full
inhomogeneous field theory.

**P3b:** A_0'-A_b'=O(epsilon^2) and xi=O(epsilon^2). Noncommuting U,V may
carry a quadratic gravitational momentum source proportional, in the
appropriate frame, to the axial part of [U,V]. Derive its coefficient and
which momentum components it actually sources. The matter configuration
response is predicted to provide the compensating momentum. Do not adopt
that prediction because a numerical constraint solver finds a small residual.

The deliberately rigid control q_i=q_i'=0 tests whether a momentum
obstruction remains for noncommuting U,V. A failure of this control does not
falsify P3a. Conversely, success for commuting U,V does not prove completion
for noncommuting data. A spatial metric Killing field alone is insufficient
to infer an obstruction for the full matter background, whose individual
fields can transform under that spatial symmetry.

If this four-unknown prescription fails, report separately (i) a proven
constraint obstruction and its exact assumptions, (ii) a singular chosen
parameterization, or (iii) an unresolved numerical solve. A residual, a
failed Newton iteration, or insufficient rank in this prescription is not
by itself a nonexistence proof for the underlying theory.

## 5. P4: finite-amplitude response and recovery of linear transport

For completed data evolve the scale, support and shape together. The
reference FRW trajectory is a control, not an externally held background.
Predict A-A_b, q_0-q_b, and transverse q_i to start at second order in
amplitude; predict the first-order beta history to be #297's supported TT
solution. Compute the leading scalar and matter-configuration responses,
including whether any claimed vector response is physical or a coordinate
rotation. Use the full constraints to account for that distinction.

Compare trajectories at equal proper elapsed time and, separately, at equal
A/a where monotonicity is verified. A nonlinear conformal endpoint can shift
by O(epsilon^2); using #298's fixed eta_star as the nonlinear endpoint would
manufacture a divergence. Record the clock conversion explicitly.

For a quantity Z use matched +epsilon and -epsilon data. The odd first
variation [Z(epsilon)-Z(-epsilon)]/(2 epsilon) must tend to the predicted
linear tangent; the even second variation
[Z(epsilon)+Z(-epsilon)-2Z(0)]/(2 epsilon^2) estimates the quadratic response.
Derive its inhomogeneous second-variation equations independently and
compare them to the finite-amplitude difference. Do not infer an exponent
from a single amplitude or fit away an initial constraint correction.

A complete Hamiltonian flow must preserve its derived symplectic form on
the physical reduced phase space after constraints and gauge are accounted
for. A tensor-only submap need not be symplectic when the support exchanges
energy or momentum with it. This round must not claim a nonlinear two-column
asymptotic tensor map by discarding the additional canonical variables.
The finite-time milestone is exact reduction plus local existence for
constraint-completed data and independently verified evolution, not merely
an approximately balanced prescribed stress.

## 6. P5: future persistence requires a proof beyond integration

The future conjecture is deliberately separate from finite-time completion.
For a stated sufficiently small neighborhood of admissible data on a fixed
d>0 expanding branch, seek a bound guaranteeing:

    A -> infinity only at infinite proper time,
    F>0 throughout, M stays positive definite,
    M -> a finite positive M_infinity,
    proper Hubble rate -> sqrt(Lambda/3).

The expected decays are physical scalar amplitudes O(A^{-1}) and proper
trace-free expansion O(A^{-2}), with a possible faster decaying component.
Finite anisotropic limiting shape is compatible with physical curvature
and shear decaying; do not require M_infinity=I. These rates are predictions
to derive, not assumptions with which to close the evolution equations.

A positive F verdict needs a closed continuation/energy estimate, a proved
regular compactified system, or validated bounds with an explicit tail
argument, covering a specified neighborhood of constraint-satisfying data.
State the norm, domain, smallness condition and dependence on d and phase.
A theorem conditional on remaining in a regular chart does not establish
that the solution remains in that chart. General results about a perfect
fluid cannot be transferred to these nonminimally coupled fields without
checking all hypotheses.

A finite-time singularity claim likewise needs more than a failed solver.
Distinguish loss of F>0, collapse, gauge/coordinate failure, and numerical
stiffness. If neither continuation nor a genuine obstruction is established,
F is UNRESOLVED; retain finite-time D/C/N results. No scalar/vector or
inhomogeneous stability conclusion follows outside the verified sector.

## 7. Frozen numerical design and independent controls

Use the STF basis E_0,...,E_4 of reciprocal_scalar_tt.py. Primary normalized
(U,V) pairs are (E0,E0), (E0,E1), (E0,E2), (E2,E3), and
((E0+E3)/sqrt(2),(E1+E4)/sqrt(2)). Include U=0,V=E0 and U=E0,V=0 controls.
Record the commutator norms rather than trusting the labels. Add eight
seeded random unit-norm STF pairs; retain every generated pair.

Primary initial-data grid: d in {.05,.15,.30}, delta=j*pi/8 for j=0,...,7,
a=kappa=1, epsilon in {0,+/-.01,+/-.02,+/-.04,+/-.08}, all fifteen pairs.
Field and field-velocity zeros are included. Constraint completion must
cover the full grid even if integration of failed data is disallowed.
Repeat representative pairs (E0,E1) and (E0,E2) at d=.15, delta in {0,pi/4,pi/2},
a in {.7,1,2}, kappa in {.4,1}, epsilon in {0,+/-.02} to test units and clocks.

For finite-time evolution use the seven explicit pairs, d=.15,
delta in {0,pi/4,pi/2}, and epsilon in {0,+/-.01,+/-.02,+/-.04,+/-.08}.
Sample dimensionless proper time t/a at {0,.25,.5,1,2,4,8}. Include all eight
random pairs at delta=pi/4, epsilon in {0,+/-.02}, with the same times.
Use the representative scale/coupling controls through t/a=2.
For tail diagnostics use the two representative pairs, delta in
{0,pi/4,pi/2}, and epsilon in {0,+/-.02,+/-.08}, through t/a=16, with
additional samples at 10,12,16. These are diagnostics, never a proof of F.

Derive the constraints symbolically. Solve their local branch to normalized
residual <1e-11 and test Jacobian rank separately. Use DOP853 at
(rtol,atol)=(1e-10,1e-12) and (1e-12,1e-14), with max_step<=.02 in t/a.
If stiffness requires another solver, record the change and independently
cross-check overlapping trajectories. Do not silently drop cases.

For on-shell full-field verification use at least the seven explicit pairs,
delta in {0,pi/4,pi/2}, epsilon in {-.08,+.08}, at t/a={0,.5,2}, and
three coordinate points (.83,1.07,.61), (1.11,.72,1.27), (.58,1.32,.94).
For off-shell identities use twenty seeded positive-metric, arbitrary-jet
cases as well as isotropic and diagonal controls. Compare all independent
Einstein components and all four scalar equations pointwise, not their
spatial averages or only their TT projections.

Normalize field-equation residuals by max(1,norm of the individual terms)
in dimensionless a=kappa=1 variables, with explicit unit conversion in
scale controls. Require exact symbolic identities where tractable and
normalized full-field errors <1e-8 using analytic jets. If numerical jets
are needed, show their convergence independently instead of relaxing this
gate. Constraint propagation and final solver comparisons must be <1e-8.
Report raw absolute errors, denominators and relative errors.

For first/second variation comparisons, report all amplitude levels and
errors against independently integrated variational equations. On intervals
where errors are above the solver floor, the central differences should
show second-order convergence. Require relative finest-amplitude errors
<1e-3 through t/a=2 for both nonzero variations. Resolve cancellation with
higher precision or leave the affected coefficient unresolved; do not
replace an unresolved coefficient by zero. Tail-rate plots cannot pass F.

Required negative/limiting controls:

- zero tensor amplitude recovers the exact evolving FRW support;
- the linear limit recovers #297 and the appropriate finite-interval #298 map;
- freezing q and A during nonlinear evolution is tested against full fields;
- dropping the nonminimal curvature term fails against full improved stress;
- omitting the momentum equation is tested on noncommuting data;
- the common-amplitude-only control is separated from responsive matter;
- wrong conformal/proper clock conversion separates from the correct one;
- a tensor-only symplectic assertion cannot substitute for full reduction;
- epsilon=0 and pure-displacement/velocity controls avoid vacuous scaling.

## 8. Verdicts, failure dependencies and deliverables

| Gate | Targets |
|---|---|
| conventions_and_units | D, C, N, F |
| exact_field_closure | D, C, N, F |
| action_full_field_agreement | D, C, N, F |
| constraint_derivation | D, C, N, F |
| constraint_completion | C, N, F |
| constraint_propagation | N, F |
| linear_recovery | N, F |
| quadratic_response | N |
| finite_time_evolution | N, F |
| future_continuation_bound | F |
| negative_controls | D, C, N, F |
| scope | D, C, N, F |
| failure_paths | D, C, N, F |

Missing/false gates and malformed/nonfinite dependent evidence make the
corresponding target UNRESOLVED, preserving unrelated verified targets.
Unknown gate schemas cannot pass. Test this dependency table itself,
including affirmative evidence, missing evidence, injected numerical
failures and actual CLI failures overwriting stale success artifacts.

Potential affirmative labels are EXACT_HOMOGENEOUS_REDUCTION_VERIFIED,
LOCAL_CONSTRAINT_COMPLETED_FAMILIES_VERIFIED,
FINITE_AMPLITUDE_RESPONSE_VERIFIED, and
FUTURE_PERSISTENCE_PROVED_ON_STATED_DOMAIN. An obstruction must have its
assumptions and certificate attached; it must not be emitted on a failed
construction or an unverified Jacobian. A failure of a prediction does not
require withdrawing independent verified results.

Deliver a separate derivation/results document, an implementation, meaningful
tests, and an archived deterministic probe with --output-dir (--output may
be an alias). Archive the constraint choices, Jacobians, full residuals,
raw trajectories, convergence data, and the exact scope of any future bound.
Record the published freeze SHA in the implementation and run. Keep this
file untouched after publication; changes in assumptions belong in the
results document or a separately identified prospective extension.

The intended milestone is a defensible nonlinear existence or obstruction
result for this specified classical field model. No Phi selection,
quantization, operational causality result, nonlinear rotor, generic
inhomogeneous stability, or preferred preparation is inferred. The already
proved ESU scalar instability is not erased by persistence on a departing
expanding branch.
