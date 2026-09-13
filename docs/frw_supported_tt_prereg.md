# Pre-registration: supported TT response along the exact FRW family

Baseline: `8ac4fd3e100b2ab75f91c08b4fae9ef605885294` (#296).
Branch: `codex/frw-supported-tt`. Seed: `2026091207`.
Publish this document before implementation, symbolic verification, finite
variations or integration. The formulas below are precomputed analytic
predictions from #296's post-review scoping, not new numerical discoveries.
Preserve this freeze; record any correction separately.

## 1. Question and scope

Derive the fully supported **homogeneous five-component TT block** on #296's
exact evolving four-scalar FRW solution. Does that sector close at linear
order, and does its operator agree with action, full stress and both clocks?
This is not the general tensor-harmonic tower, a coupled rotor construction,
a nonlinear persistence result, or an adiabatic/Bogoliubov experiment.

Use exactly the independent real conformal components of #294–#296 with
xi=1/6 and no additional matter, constitutive law or damping. In conformal
time eta the background is

    ds^2=A(eta)^2(-deta^2+gamma), phi_I=q(eta)x_I/A(eta),
    q=a sqrt(3/(4 kappa)) cos(2eta+delta), Lambda=3/(2a^2),
    A'=sigma (A^2-a^2)/(sqrt(2)a), sigma=+1 or -1,
    w=(A-a)/(A+a)=w0 exp(sigma sqrt(2) eta),
    w0=epsilon_background/(2+epsilon_background).

The background departure epsilon_background is independent of the tensor
variation parameter epsilon. It is not expanded in this round. The local
chart requires A>0 and F=1/kappa-q^2/(6A^2)>0. Do not integrate across a
scale pole, A=0 or F=0. The ESU is an additional limiting control, not the
only background tested. The homogeneous scalar instability established by
#296 remains; zero scalar perturbation data defines the sector here.

Use the invariant coframe and real orthonormal STF basis of #295:

    g_ij=A^2 [exp(2 epsilon beta)]_ij, tr beta=0,
    beta=sum_{B=1}^5 b_B E_B, tr(E_B E_C)=delta_BC.

## 2. P1: scalar and constraint closure

For each of the four fields, predict delta(Box-R/6)phi_I=0 for this TT
variation, using Hess x_I=-gamma x_I, transversality and zero trace. Include
time-dependent A and q, not merely the static Hessian contraction. Show that
delta phi=delta phi'=0 initial data remain zero by linear KG uniqueness.
Do not set the stress perturbation to zero as a consequence.

All linear Hamiltonian, momentum, spatial trace and scalar curvature
variations are predicted to vanish. The five spatial TF equations must all
be satisfied on the derived tensor evolution. Check the full mixed Einstein
residual, not a projection that drops an unwanted component. If another
sector is sourced, report failure of closure rather than discarding it.

## 3. P2: action and full-stress predictions

Write Q=q^2/A^2, F=1/kappa-Q/6, h=A'/A. The expected quadratic action is

    L2=Vol(S3_unit)/2 [M tr(beta'^2)-K tr(beta^2)],
    M=A^2 F=A^2/kappa-q^2/6,
    K=8A^2/kappa+(2/3)q^2,
    (M beta')'+K beta=0.

Derive this with R3=A^-2[6-8 epsilon^2 tr(beta^2)+...],
sum |grad phi_I|^2=(Q/A^2) tr exp(-2 epsilon beta), and the ADM curvature.
Here tr K_extrinsic=3 A'/A^2, not zero. Demonstrate that it and its weighted
boundary term are tensor-independent in the exponential fixed-volume
parameterization, even when F depends on time. Assert that the Lambda term
has no tensor quadratic variation. A traceless linear metric continuation
has a different quadratic volume and is not an action substitute.

Independently predict the full mixed spatial TF variations

    delta G = A^-2 [beta''+2h beta'+8 beta],
    delta T = -2Q beta/A^2 + Q' beta'/(6A^2) + (Q/6) delta G,
    delta(G-kappa T) = kappa/A^4 [M beta''+M' beta'+K beta].

The mixed components remove the background pressure times delta g. The
derivative Q' includes the expansion contribution -2hQ. Terms proportional
to beta, beta' and beta'' must be varied separately, with additional
noncommuting STF data. No division by q or q' is permitted.

## 4. P3: propagation and controls

Use normalized coefficients m=kappa M/a^2, k=kappa K/a^2. The pair
(b,p=m b') has generator [[0,1/m],[-k,0]], a constant rescaling of the
physical canonical momentum. With y=sqrt(m)b the independent normal form is

    y''+[k/m - m''/(2m) + m'^2/(4m^2)] y=0.

Verify both transformations at the initial and final endpoints; a varying
momentum normalization must not be mistaken for damping. In proper time
dt=A deta the same equation is

    d/dt(A M db/dt)+(K/A)b=0.

Compare finite-interval fundamental matrices from conformal canonical,
normal-form and independent proper-time/background integrations. This is
transport on a nonperiodic background: do not classify the eigenvalues of
an arbitrarily chosen finite-interval map as Floquet stability multipliers.
Do not infer an adiabatic invariant, preferred frequency basis, particle
production, a rotor speed, quantum amplitudes or Phi from a determinant.

At A=a, m and k must reproduce #295's f and g, including its pi/2 period map.
At epsilon_background -> 0, transport on a fixed interval must approach the
coupled ESU map, not the bare map. At q=0 the instantaneous stress variation
vanishes; this is not an all-time decoupling. Setting the whole support
response to zero is allowed only as an operator negative control, since that
would not source the same background.

Required negative controls: bare FRW acceleration beta''=-2h beta'-8 beta;
wrong sign of the M' term; omit -2hQ in Q'; 5% wrong on-equation acceleration;
and the linear-volume substitution in the Lambda quadratic term. Use phases
where each control is active. The bare canonical potential 8-A''/A must not
be accepted as the supported normal form merely because it is familiar.

## 5. Frozen numerical domain and independent route

Use epsilon_background=0.1,0.3,-0.1, sigma=+1,-1, delta=0,0.31,pi/4,
a=0.7,1,2 and kappa=0.4,1. Primary transport interval: eta in [0,0.8].
All selected cases must remain in A>0 and m>0.2; check this explicitly.
Use the same backgrounds for off-shell variations at eta=0,0.37,pi/4,0.8,
including field/velocity zeros and at least three generic spatial points.
It is acceptable to use a deterministic representative subset of the
Cartesian product for expensive curvature checks, but report the selections.
Every tensor basis direction and all three time-jet coefficients are required.

Use independent coordinate metric jets, full curvature and the summed
off-shell improved stress. A linear continuation I+2 epsilon beta is allowed
for this first variation only. Compare conformal lapse N=A against proper
lapse N=1 after transforming every background, field and tensor time jet.
Check all equations, baseline stresses, parity and clock covariance.

Central variations epsilon=.002,.001,.0005 must match predicted linear
coefficients to normalized error <2e-4. Require halving ratios 3.5..4.5 where
the finer normalized error exceeds 1e-10. Baseline normalized Einstein/KG
residual and full finite-continuation clock disagreement must be <1e-9.
Integrate with DOP853 at (rtol,atol)=(1e-10,1e-12),(1e-12,1e-14),
max_step<=pi/100 in conformal time. Require relative map agreement <1e-8,
canonical determinant/symplectic error <1e-8, composition and inverse-map
errors <1e-8, and proper-clock endpoint error <1e-9. Compare ESU convergence
at epsilon_background=1e-2,1e-3,1e-4 without fitting an asymptotic spectrum.

## 6. Evidence gates and stopping rule

Required gate names:

    exact_background, action_derivation, full_stress_response,
    scalar_constraint_closure, independent_geometry, two_clocks,
    static_limit, canonical_transport, negative_controls, failure_paths.

Machine verdicts: linear_tensor_sector, supported_operator, clock_agreement,
finite_interval_transport. Require every gate and finite, well-shaped
independent map evidence. Missing/false gates, missing identities, malformed
or nonfinite reports and exceptions must produce UNRESOLVED and a failing
CLI exit rather than retain a success artifact. Exercise failure paths from
a genuinely passing baseline. Use --output-dir (optional --output alias).

Stop after the equation, closure and finite-interval transport are established
or refuted. Full nonlinear stability, tensor-driven secondary scalar/vector
sources, general tensor harmonics, admissible director histories, adiabatic
action, Bogoliubov bases, preparation selection, Phi and operational causality
remain outside this round. Correcting the equation is the prerequisite for
those later questions, not their answer.
