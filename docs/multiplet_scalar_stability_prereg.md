# Pre-registration: low scalar modes of the four-component ESU

Baseline: `9a0bfbdbc553e85a4bef5632ee5184d760db95dd` (#295).
Branch: `codex/multiplet-scalar-stability`. Seed: `2026091106`.
Publish this freeze before implementation, symbolic verification or numerical
measurements. The homogeneous formulas below are analytic predictions from
scoping; they are not experimental discoveries. Preserve this document if a
prediction fails and record corrections separately.

## 1. Question and inherited assumptions

Does the background underlying #294–#295 have a physical growing homogeneous
scalar perturbation? What becomes of the degree-1 scalar sector after gauge,
constraints and the inherited antipodal restriction are accounted for?
This takes priority over secondary scalar/vector sources quadratic in the
tensor perturbation and over redoing #292's driven rotor.

Use precisely the four independent real conformal fields of #294:

    I = integral sqrt(-g) [(R-2 Lambda)/(2 kappa)
        - sum_I ((grad phi_I)^2 + R phi_I^2/6)/2],
    phi_I = sqrt(3/(4 kappa)) cos(2t/a + delta) x_I,
    sum_I x_I^2=1, Lambda=3/(2a^2).

No fluid perturbation law, extra field, damping or junction is added.
The radiation-like background stress alone does not determine perturbations.
The prior is at least one growing homogeneous mode (Eddington-type); an
elliptic physical homogeneous block would contradict this prior.

The degree in this document labels the metric scalar harmonic under the
diagonal spatial/internal SO(4), not the degree of each perturbed field.
This round is linear in perturbation amplitude, except for an exact FRW
continuation used to check the homogeneous constraint and physical growth.
It does not classify all representations of four arbitrary field perturbations.

## 2. P1: homogeneous reduction, prediction to re-derive

In conformal time eta allow g=A(eta)^2(-deta^2+gamma) and
phi_I=q(eta)x_I/A(eta). Derive the reduced equations from the summed full
improved stress and Klein–Gordon equations, independently of any fluid model.
The predicted equations and conserved conformal energy are

    q'' + 4q = 0, E=(q'^2+4q^2)/2,
    rho=E/A^4, p=rho/3,
    A'^2+A^2-Lambda A^4/3=kappa E/3,
    A''+A-(2 Lambda/3)A^3=0.

At A=a, E0=3a^2/(2 kappa), q0=a sqrt(3/(4 kappa)) cos(2eta+delta).
Writing A=a(1+epsilon r), q=q0+epsilon z predicts

    r''=2r, z''+4z=0, delta E=q0'z'+4q0 z=0.

Do not count an unconstrained energy perturbation as a physical linear mode.
The surviving oscillator perturbation should be a phase/time-origin mode;
distinguish a fixed-clock phase change from the common time-translation gauge
orbit. The scale perturbation is invariant under linear time shifts because
A0'=0; test proper spatial volume as a physical observable.

The ESU Hamiltonian constraint is degenerate at first order in metric data.
A growing eigenvector must therefore be tangent to an actual constraint-
satisfying family, not just solve linear evolution. At fixed E=E0 the predicted
exact constraint is

    A'^2=(A^2-a^2)^2/(2a^2).

Its two first-order branches A'=+/- (A^2-a^2)/(sqrt(2)a) are candidate exact
continuations of the growing/decaying directions. Verify full Einstein–matter
residuals, positivity of A and F near the background, and proper time
dt=A deta. A linear exponential is not a claim of exponential growth for all
finite amplitudes or all times.

## 3. P2: degree 1, open reduction rather than a prior verdict

For Y=d.x, |d|=1, use the general diagonal-SO(4) scalar-type matter pair

    delta phi_I = u(eta) x_I Y + v(eta) D^a x_I D_a Y
               = (u-v) x_I Y + v d_I.

Include scalar lapse, scalar shift and spatial trace before fixing gauge.
Hess(Y)=-gamma Y, so its trace-free Hessian vanishes identically. Never divide
by l(l+2)-3 or infer a lapse-slip relation from that zero tensor. Demonstrate
the gauge fixing and any residual transformations. A round spatial metric
and zero shift is a candidate gauge; establish whether it is complete.
Derive Hamiltonian, momentum, spatial Einstein and all four KG equations.
Only evolve a reduced physical system after verifying the discarded equations
and constraint propagation. If reduction is singular, use regular patches or
report the physical classification unresolved. No division by a field that
crosses zero may masquerade as a dynamical obstruction.

The displayed matter perturbations are antipodally even, while the background
fields are odd; the metric scalar Y is odd. Report two distinct statements:
(i) the scalar-type degree-1 block on the unrestricted S3 cover; (ii) whether
it is admissible if componentwise odd fields and an antipodally invariant
metric are required. The latter metric restriction must be named explicitly,
not silently imported as a field equation. An excluded sector has no physical
period map in that restricted preparation class; it is not an elliptic mode.
An SO(4)-covariant completeness argument for the displayed scalar-type pair
is required; do not claim it covers vector-type or unrelated internal modes.

## 4. Period maps and independent verification

Use the full field period T_eta=pi as the primary Floquet interval; the fields
change sign after pi/2. For comparison with #295 report the homogeneous metric
block also over pi/2, explaining any internal-sign identification before using
a half-period matter map. In conformal background time eta=t/a the homogeneous
prediction is eigenvalues exp(+/-sqrt(2) T); the proper-time rate is sqrt(2)/a.
This is an analytic prediction, not a frozen measured decimal.

Use DOP853 at (rtol,atol)=(1e-10,1e-12) and (1e-12,1e-14), max_step<=pi/100.
Compare maps with an independent analytic matrix exponential when available,
and canonical determinant/symplectic constraints. Require map difference
<1e-8 in relative Frobenius norm, determinant error <1e-7 and analytic-map
relative error <1e-8. Test phases 0, 0.31, pi/4 and radii 0.7,1,2 with kappa
0.4,1. Use one-period maps for classification; long integrations are checks.
Report neutral/Jordan blocks separately from exponential growth. No exponent
from a lapse or time-origin mode licenses physical instability.

Use independent coordinate metric jets with the existing full-curvature and
off-shell improved-stress routines. Check all four fields, Hamiltonian,
momentum, spatial trace and TF equations at non-special spatial points and
both field and velocity zeros. Verify arbitrary off-shell variations as well
as solutions. Central variations at epsilon=.002,.001,.0005 must agree with
the derived linear coefficients to normalized error <2e-4; halving should be
consistent with second order when above a 1e-10 numerical floor. Exact FRW
solutions must have normalized full residual <1e-9. Include a wrong acceleration,
a nonzero linear energy perturbation and a pure-gauge degree-1 control.

## 5. Gates and interpretation

Required gate names are fixed:

    field_reduction, homogeneous_constraints, exact_continuation,
    clock_and_gauge, dipole_reduction, dipole_parity,
    independent_geometry, period_maps, negative_controls, failure_paths.

All ten gates are required for the combined affirmative report. Separately
retain the evidence for homogeneous and dipole claims so an unresolved dipole
does not erase a valid homogeneous theorem. Machine verdict fields:
homogeneous_physical_block, homogeneous_constraint_completion,
dipole_cover_block, dipole_restricted_admissibility. A failed, missing, unknown,
malformed or nonfinite required item must not yield an affirmative combined
verdict. Test failure paths starting from genuinely passing evidence. Do not
infer stability from positivity of F, a small residual or bounded gauge data.

A verified admissible growing homogeneous mode refutes linear stability of
this background even if the tensor block remains elliptic. It does not refute
existence of the background, #294's equal-stress family or #295's sector result.
It establishes neither a measure on preparations nor a measure-zero viability
claim, an inevitable trapped surface, nonlinear instability of every nearby
history, nor a failure of two-boundary selection. A stable result in these two
blocks would still not prove full stability. Vector response, tensor-driven
secondary scalar/vector sources, a coupled rotor, O(s^4) persistence, Phi
selection and the causality gate remain outside this round.
