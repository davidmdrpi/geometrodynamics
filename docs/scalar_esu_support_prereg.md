# Freeze: can the existing real conformal scalar supply the ESU support?

Baseline: `6631aa3057b16ee96c2a35f45fc43a848d65fca4` (main including #291).
Branch: `codex/scalar-esu-support`. Seed: `2026090915`.
This file is to be published before implementation or measurements. The
analytic predictions below are predictions, not results discovered by a scan.
Preserve this file; record any correction in the implementation report.

## 1. Question and exact scope

Can the *one existing real conformal scalar*, with no other matter, replace
the unspecified supporting fluid in the round 3+1 Einstein static universe?
Require the background Einstein equations pointwise for all time. No spatial,
temporal, ensemble, or multiplet average may replace a pointwise equation.

Use metric signature (-+++), radius a>0, kappa>0 and

    I = integral sqrt(-g) [(R-2 Lambda)/(2 kappa)
                          - (grad phi)^2/2 - R phi^2/12].

This is the curvature coupling and improved stress already used in
`geometrodynamics/waves/two_wave.py` and `waves/backreaction.py`. No new
potential, scalar species, phenomenological fluid or support equation of
state is added. Lambda is a background parameter whose required relation
to a must be reported, not fitted away in each Einstein component.

Primary admissibility condition: phi(t,-x)=-phi(t,x), imposed as the BAM
odd-field sector, not derived by this round (see #291). The mathematical
class is all smooth real fields on complete round S3, not merely finitely
many harmonics. Distinguish this conditional admissibility from the question
without an antipodal restriction. Include a homogeneous even-field control.

This is ESU *bulk* support, not 5D throat support: `bulk/source_audit.py`
and `shells/junction.py` answer different dimensional/geometric questions.
There is no claim about multiple scalars, complex fields with independent
components, singular fields, non-round metrics, averaged backgrounds,
additional stresses, nonminimal couplings other than 1/6, or modified gravity.

## 2. Frozen analytic predictions

### P1. Pointwise isotropy has a global consequence

On a fixed round spatial slice, the spatial trace-free part of the improved
stress is

    T_ij^TF = (2/3) (D_i phi D_j phi)^TF
              - (1/3) phi (D_i D_j phi)^TF.

No use of the scalar wave equation or averaging is needed for this identity.
On any connected open component where phi != 0, put h=1/phi. Isotropy is
equivalent there to (D_i D_j h)^TF=0. Derive by curvature commutation that

    h = A + B dot x

locally, with A and B constant on the connected component (unit embedding
coordinate x, physical metric a^2 dOmega_3^2). Prove the global step, rather
than assuming division by phi is valid at its zeros: the displayed h is
bounded on the complete sphere, so phi=1/h cannot approach zero on the
boundary of a nonzero component. Thus that component is closed as well as
open. Smoothness also excludes a pole of 1/h. Every nontrivial isotropic
slice should therefore be nowhere zero and have |A|>|B|.

A continuous real odd function on S3 has a zero. P1 should exclude every
nonzero odd isotropic slice, hence every nontrivial smooth odd scalar
history supplying the ESU. A slice on which phi is identically zero is not
excluded by this spatial argument; the all-time odd condition must be used
to finish the history argument. The zero history cannot supply the required
enthalpy 2/(kappa a^2). This is the predicted obstruction, conditional on
the explicitly stated class. Failure to establish its global step leaves
the general-class verdict UNRESOLVED, regardless of numerical scans.

Do not infer existence of a full Einstein-scalar history from the local
spatial form 1/(A+B dot x); momentum, density and time equations remain.

### P2. A positive even-field control

Without odd admissibility, predict an exact homogeneous solution

    phi_0(t) = sqrt(3/kappa) cos(t/a + delta),
    Lambda = 3/(2a^2),
    rho_0 = 3/(2 kappa a^2),   p_0 = rho_0/3.

Derive this from the scalar equation and all Einstein components, using
the actual improved stress, including at phi_0=0 and dot(phi_0)=0.
The coefficient multiplying R/2 in the action is

    F = 1/kappa - phi_0^2/6 >= 1/(2 kappa) > 0.

Check the scalar kinetic coefficient after the regular conformal change
g_E=(kappa F)g, namely

    K_E = 1/(kappa F) + 3/(2 kappa) (F_phi/F)^2.

This is a regularity/kinetic-sign check, not a proof of global stability.
The homogeneous solution is even and is not an admissible positive witness
for the primary odd-sector question. Its phase is chosen. Background
existence does not select a statistical measure.

### P3. Radiation background does not imply a perfect-fluid response

For the control, use scalar perturbations in longitudinal gauge

    ds^2 = -(1+2 alpha)dt^2 + a^2(1-2 psi)dOmega_3^2,
    phi = phi_0 + chi.

Derive delta rho, the momentum potential J with j=D J, delta p and the
scalar anisotropic-stress potential Pi, including their metric-dependent
terms. Predict the trace-free equation for degrees l>=2:

    Pi = -phi_0 chi/3 + phi_0^2 (psi-alpha)/6,
    (1-kappa phi_0^2/6)(psi-alpha) = -kappa phi_0 chi/3.

Thus the generic perturbation has nonzero anisotropic stress even though
the on-shell stress is traceless. Do not identify delta p=delta rho/3
with the zero-anisotropic-stress closure assumed in #290.

Derive and evolve the linearized scalar equation and spatial Einstein
equation; retain Hamiltonian and momentum constraints as independent
propagation checks. Treat at least degrees 2 and 3, over a complete scalar
period, including the two kinds of background turning point. Constraint
solves must not divide by phi_0 or dot(phi_0). This response belongs to the
even-field control, not to an admissible odd background. No claim about
homogeneous or degree-1 perturbation stability is licensed by these runs.

If support and signal are the same scalar, T[phi_0+epsilon chi] has cross
terms linear in epsilon. Explicitly report that #289--#292's zero-background
scalar power counting cannot be transferred unchanged to this control.

## 3. Mandatory checks, fixed before implementation

1. `stress_identity`: derive P1 from the action's improved stress and check
   against the inherited stress implementation on nontrivial odd harmonics.
2. `global_obstruction`: supply the curvature-commutation calculation, the
   local reciprocal classification and the component-boundary proof. Test
   the local identity numerically on 1/(A+B dot x), including nonconstant
   positive configurations. A collection of failed odd candidates alone
   cannot satisfy this gate.
3. `homogeneous_background`: exact symbolic scalar/Einstein residuals and
   independent pointwise full-stress checks of P2 over a full period; radius
   controls a=0.7,1,2 and kappa=0.4,1. Reversed Lambda and amplitude multiplied
   by 1.1 must fail at least one background equation.
4. `admissibility`: distinguish smooth odd-sector exclusion, even control,
   zero field and reciprocal functions with poles. Enforce F>0 and check
   K_E>0. A positive even control must never yield odd-sector existence.
5. `response_variation`: derive all four scalar stress responses. Compare
   them with direct metric/connection/Ricci evaluation and the improved
   stress under small finite perturbations, not just with algebraically
   rearranged copies. Use central differences at two steps and generic
   field/metric time jets, including off-shell scalar jets. Verify the
   off-shell trace identity T^mu_mu=phi(Box-R/6)phi.
6. `constraint_propagation`: solve initial Hamiltonian and momentum
   constraints and independently evolve the scalar/spatial Einstein
   equations for l=2,3, primary a=kappa=1 over 2 pi a. Check both constraints,
   the trace identity and the anisotropic equation throughout. Repeat at
   (rtol,atol)=(1e-10,1e-12) and (1e-12,1e-14). Include a=0.7,kappa=0.4 as
   a response scaling control. No projected-back-to-constraint integration.
7. `nonfluid_control`: show nonzero Pi in at least one constraint-compatible
   response. Forcing Pi=0 must leave a nonzero trace-free Einstein residual
   for that same response. No conclusion about every possible perturbation.
8. `scope_and_order`: record the even background's excluded parity, the
   chosen phase/Lambda relation, the first-order stress cross term, and the
   missing triangle, preparation, measure and physical-readout maps.
9. `fail_closed`: every required check is addressed by exact name. Deleting
   or failing any one produces UNRESOLVED physical verdicts, overwritten
   JSON/Markdown and a nonzero process exit. Unknown true keys cannot pass.

Exact symbolic identities must simplify to zero. Independent spatial
stress agreement tolerance: 1e-10 relative to a stated nonzero scale;
reciprocal-isotropy residual: 1e-10 normalized. Finite-variation response
error: 1e-6 normalized, using epsilon=1e-3 and 5e-4 (the observed order is
reported; no convergence claim from a residual alone). Initial and evolved
constraint/Einstein residuals: 1e-8 normalized; refined trajectories: 1e-7.
Report absolute as well as relative residuals and floors at zeros.

No search objective is used. If an analytic prediction needs correction,
preserve the freeze and state the correction before reporting verdicts.
Numerical failure is UNRESOLVED; a scoped exclusion requires the proof.

## 4. Deliverables and separate verdicts

Add an isolated module, probe, tests, derivation and reproducible JSON/Markdown
archive. Update the audit with separate fields for:

- odd-sector exact ESU support;
- unrestricted homogeneous control;
- kinetic/regularity admissibility of that control;
- its constitutive response and propagated constraints;
- BAM support selection;
- triangle map, Phi selection and causality gate.

Even if P1--P3 all hold, BAM support selection is NOT_DERIVED and the last
three maps/gates remain open. The result would exclude one concrete way of
filling the support gap and identify why a familiar even-field workaround
does not inherit either odd admissibility or the assumed fluid response.

Sources already in the baseline: `docs/esu_support_response.md`,
`docs/parity_solvability.md`, `geometrodynamics/waves/two_wave.py`,
`geometrodynamics/waves/backreaction.py`,
`geometrodynamics/waves/reciprocal_scalar_tt.py`,
`geometrodynamics/waves/scalar_tt_constraints.py`,
`geometrodynamics/waves/esu_support_response.py`,
`geometrodynamics/bulk/source_audit.py`, `geometrodynamics/shells/junction.py`.
