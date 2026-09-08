# Freeze: ESU support response and the first cubic scalar reaction

Baseline: main `080c1cc58cd7a01fc619ff5d933c5fdfb5b8cf33`, including #289.
Publish this document before implementing or measuring the new coefficient.
Seed `2026090713`. Claude's separate task is the necessity of parity purity;
this branch assumes the same degree-3 source as #289 and does not answer it.

## Decision question and scope

Does the omitted scalar metric response cancel the retained homogeneous TT
reaction at order s^3 for the specified preparation? First target the initial
proper-time acceleration coefficient, which is narrower than a full evolution.
Then derive the linear scalar-metric/support response with propagated scalar
Einstein constraints, keeping its constitutive dependence explicit.

The existing ESU channel assumes an unspecified perfect-fluid support.
Its background equations fix

    kappa rho_0 = 3/a^2 - Lambda,
    kappa p_0 = Lambda - 1/a^2,
    H_f := rho_0+p_0 = 2/(kappa a^2).

They do not fix a pressure derivative, entropy response, or a BAM microscopic
realization of that fluid. Use the named class: perfect fluid, no anisotropic
stress, no exchange with the conformal scalar, adiabatic linear perturbations
delta p_f = c_s^2 delta rho_f. Treat c_s^2 as a parameter, not a fitted or
derived constant. The initial coefficient may be independent of it; a later
history generally is not. No support field, canonical ensemble or Phi is added.

This is a perturbative scalar-sector completion at order s^2 plus a projected
order-s^3 force. Inhomogeneous tensor responses are not evolved. No complete
Einstein–matter history, finite-scalar-multiplet closure of the corrected wave
equation, source-local apparatus, or causality/readout verdict is claimed.

## Gauge, clock, and preparation

Use Newtonian gauge for l>=2,

    ds^2 = -(1+2 alpha) dt^2 + (1-2 psi) gbar_ij dx^i dx^j,

with zero scalar shift/shear. Set the spatially homogeneous lapse alpha_0=0;
retain the homogeneous spatial scale response. This fixes the background
clock convention. The physical fluid proper derivative is D_U = U^mu nabla_mu.
Do not identify coordinate acceleration with D_U^2 phi.

Leading scalar: phi_1=A(t)Y, A=s cos(omega t), omega=4/a, with the inherited
normalized coherent Y=sqrt(8/V) Re[(x0+i m.x)^3], m=(1,2,3)/sqrt(14),
V=2 pi^2 a^3. Initial support delta rho_f=delta p_f=0 and fluid velocity=0.
All free tensor data, including homogeneous TT, are zero. The initial metric
is constraint-completed, not zero: psi(0)=-2u(0) from #289 and psi_dot(0)=0.
The initial slice is CMC; Newtonian gauge evolution is not a sequence of the
rigid CMC slices in #289.

All sources and added metric coefficients are order s^2. Corrections to the
leading scalar are order s^3 and first alter the matter stress at order s^4.
The computed leading-order source and linear gravitational response are thus
appropriate for this coefficient. Do not integrate a projected cubic force
and call it the complete corrected scalar field.

## Analytic predictions before measurement

Write Y^2=sum h_l Y_l, lambda_l=l(l+2)/a^2, L_l=3/a^2-lambda_l. The source
has degrees 0,2,4,6, with the exact powers established in #289:

    V ||h_l||^2 = 1, 27/25, 1/5, 201/175.

The improved stress gives

    rho_phi,l = [s^2 omega^2/2 - lambda_l A^2/12] h_l,
    J_phi,l = -A A_dot h_l/6,       j_phi = grad J_phi,
    p_phi = rho_phi/3.

Define the scalar anisotropic-stress potential by
pi_phi^S = (Hess - gbar Delta/3) Pi. Source conservation predicts

    Pi_l = -(3 J_phi_dot,l + rho_phi,l)/(2 L_l)
         = -A^2(omega^2-lambda_l/12) h_l/(2 L_l),  l>=2.

Check Pi independently by projecting the pointwise improved spatial stress
on scalar STF Hessian harmonics. Do not replace the full anisotropic stress
by this scalar projection; its other tensor components remain present.

The linearized Einstein equations predict

    2(Delta+3/a^2)psi = kappa(rho_f+rho_phi),
    -2 grad psi_dot = kappa(j_f+j_phi),
    psi_l-alpha_l = kappa Pi_l, l>=2,
    psi_ddot-psi/a^2 - Delta(psi-alpha)/3 = kappa(p_f+p_phi)/2.

Thus the scalar-sector evolution reduces, for each of l=0,2,4,6, to

    psi_ddot,l = [1/a^2+c_s^2 L_l] psi_l
        + kappa(1/3-c_s^2) rho_phi,l/2 - kappa lambda_l Pi_l/3,

where the last term is absent at l=0. Recover the support from

    rho_f,l = 2 L_l psi_l/kappa-rho_phi,l,
    J_f,l = -2 psi_dot,l/kappa-J_phi,l, l>=2.

The spatially constant J is immaterial and set to zero. Independently evolve
the fluid conservation equations

    rho_f_dot = -Delta J_f + 3 H_f psi_dot,
    J_f_dot = -c_s^2 rho_f - H_f alpha, l>=2,

and verify both Einstein constraints as unused equations. Compare against the
reduced oscillator system, rather than reporting algebraic reconstruction
residuals alone. The l=0 fluid energy equation remains active; no homogeneous
Euler-potential equation is imposed.

The free-support scalar frequency is
Omega_l^2=[c_s^2(l(l+2)-3)-1]/a^2. This predicts the usual l=2 threshold
c_s^2=1/5 and an unstable homogeneous scale mode. Do not silently remove the
mean to obtain stability. Finite-time coefficients need not be stable forever.
The standard fluid-only threshold is a cross-check, not a new discovery:
Barrow et al., https://arxiv.org/abs/gr-qc/0302094, equation (15).

## The cubic coefficient and why the clock matters

The predicted coordinate-time wave-operator correction, evaluated on phi_1, is

    F_t = 2 alpha phi_1,tt + (alpha_dot+3 psi_dot) phi_1,t
          + 2 psi Delta phi_1 + grad(alpha-psi).grad phi_1
          - delta R^(4) phi_1/6,
    delta R^(4) = 4 Delta psi+12 psi/a^2-6 psi_ddot-2 Delta alpha
                = kappa(1-3 c_s^2) rho_f.

The trace equality must be independently checked against the geometric
curvature expression. The fluid proper derivative differs from coordinate
acceleration by lapse and material-advection terms. On the initial slice the
fluid is at rest and has zero pressure gradient, hence geodesic acceleration,
and the scalar curvature perturbation vanishes. Therefore the initial proper
force is predicted to be

    F_proper(0) = 2 psi(0) Delta phi_1(0)
                 - grad psi(0).grad phi_1(0).

This is also obtained from D_U^2 phi = Delta_g phi-R^(4)phi/6 on that slice.
Project this local scalar correction onto the fixed initial Y with the leading
round volume. This is a stated modal diagnostic, not an operational pointer.
Changing the projection measure at order s^2 changes this order-s^3 correction
only at higher order; it does not turn the whole projected field amplitude
into a claimed measured record.

From the exact harmonic powers, the frozen predictions are

    integral Y F_proper(0) dV = -(7976/875) kappa s^3/(V a^2),
    integral Y F_t(0) dV      = +(55096/875) kappa s^3/(V a^2).

Their different signs are a clock/advection distinction, not two physical
predictions for one observable. Derive both and validate the conversion.
The initial result is independent of c_s^2 within the stated initial-data and
perfect-fluid class. Entropy/pressure perturbations at preparation, anisotropic
support, different free gravitational data, or another field preparation are
outside this certificate.

The induced homogeneous TT has beta_ind(0)=beta_ind_dot(0)=0, so its cubic
force initially vanishes. A nonzero proper coefficient therefore excludes
identical scalar/TT cancellation for this preparation; it does not prove that
their sum never vanishes at a later time, or for all BAM supports/preparations.
Use #289's exact induced TT solution for the later comparison. At later times
report the coordinate-gauge modal force explicitly as such. Higher odd scalar
degrees can be sourced; no corrected 16-mode closure is asserted.

## Frozen implementation and verification

Primary a=kappa=1, s=0.02, c_s^2=1/3, t in [0,2] at 201 times. Controls:
c_s^2=0,1/5,1; amplitudes 0.01 and 0.005; radii 0.7 and 2 for normalization.
The small amplitude and finite interval are chosen before measurement to
avoid promoting the homogeneous ESU instability into a finite-amplitude
perturbative evolution. Report maximum pointwise |alpha|,|psi|, requiring
both below 0.05 for the finite-amplitude primary/control demonstrations.

1. Derive the linearized Einstein equations from a metric/connection
   variation in an S3 coordinate chart. Require exact symbolic residuals
   for the Hamiltonian, momentum, spatial trace, scalar STF and curvature.
2. Check scalar rho,p,J and Pi against inherited improved stress on S3 rules
   (8,16) and (12,24), including nonzero momentum times. Scaled residual
   below 1e-9, with scaled=Euclidean error/max(1, reference norm).
3. Integrate the reduced metric equation and independently the fluid
   conservation equations plus spatial Einstein evolution with DOP853 at
   (rtol,atol)=(1e-10,1e-12) and (1e-12,1e-14). Require reconstructed support
   agreement, unused Hamiltonian/momentum residuals and tolerance refinement
   below scaled 1e-8. Include the homogeneous mode and all even source modes.
4. Verify delta R^(4)=kappa(1-3 c_s^2)rho_f to scaled 1e-9. Recover the
   fluid-only frequency and l=2 stability threshold algebraically.
5. Check F_t against the scalar wave operator on an exponential metric
   N=exp(epsilon alpha), g=exp(-2 epsilon psi)gbar, using its full ADM scalar
   curvature. Centered differences at epsilon=0.04,0.02,0.01 and Richardson
   extrapolation must agree to scaled 1e-7. Keep nonminimal curvature terms.
6. Verify the two initial coefficients to scaled 1e-9 by independent
   pointwise quadrature, exact harmonic powers and proper-time conversion.
   Amplitude-halving force ratios equal eight to 1e-8; metric ratios equal
   four. A coordinate-time coefficient cannot substitute for the proper one.
7. Derive or verify #289's induced TT force independently; record both
   scalar and TT modal forces. A failed witness proves no cancellation
   theorem. The claimed obstruction requires the nonzero exact initial
   proper coefficient with the frozen tensor initial data.
8. A missing or failed required gate must render stable UNRESOLVED verdicts,
   name the failures, overwrite stale reports with valid JSON and exit 1.

Preserve the freeze. Corrections, if needed, are separate notes. Report
separate verdicts for assumed support class, propagated scalar constraints,
initial proper coefficient, preparation-scoped cancellation, full evolution,
BAM support selection and physical readout. An unspecified support class is
not silently promoted to a derived BAM matter model.

Deliver one reusable module, probe, tests, compact archive and derivation.
Update the audit with #289's metric-size limitation and this round's scoped
coefficient. Leave the parity-necessity question for Claude. Open the next PR.
