# Pre-registration: field-derived normal stress for a driven TT rotor

Baseline main `c08f46a56def92edc3ab1ad1146973ab90c40aa2` (#290 merged).
Publish this file before implementation or measurements. Seed `2026090814`.
No numerical search or evaluation of the new candidate has preceded this
freeze. The identities and candidate below are analytic predictions to be
checked, with corrections recorded rather than the freeze rewritten.

## Question and scope

Does an existing real, degree-3 conformal scalar supply the normal stress
needed for a nonzero, constant-amplitude, uniformly rotating uniaxial
**homogeneous TT projection**, with admissible leading constraint data?

Use the action and normalization of `waves/reciprocal_scalar_tt.py`, the
embedding and projector of `bulk/tt_triangle_rotor.py`, and the inherited
improved stress. The manufactured drive in the original rotor round is a
control only. No freely specified time-dependent stress is evidence for a
field realization.

Order counting is part of the question: phi=s phi_1+O(s^3), beta=O(s^2),
with fixed a,kappa and finite t/a. The free scalar fixes stress at O(s^2).
The scalar feedback of both TT and scalar/support metric sectors first
changes stress at O(s^4). Thus a leading normal-balance test is meaningful
without pretending #290 supplied a full corrected history. Other metric
harmonics can be sourced at O(s^2); they do not mix into this homogeneous
TT projection in the linear metric equation. Their evolution and all
O(s^4) gravitational and matter corrections remain outside this round.

Use the same conditional background support class as #290: perfect fluid
with no anisotropic stress, initially zero support perturbations. Its
constitutive law is not derived from BAM. Check CMC Hamiltonian and all ten
momentum compatibility conditions on the prepared slice. Leading free
source compatibility will also be checked along the history; no CMC gauge
is imposed at every later time and no constraint propagation is inferred
merely from these compatibility integrals.

## P1. The field source and the equations that must hold

Let E_mu be the inherited orthonormal STF basis, C=V/kappa,
V=2 pi^2 a^3, omega_T^2=8/a^2, omega_s=4/a. The inherited action gives

    S(q)=sum_mu (q^T F_mu q) E_mu / C,
    beta_ddot + omega_T^2 beta = S(q).

Retain the factor 1/C: `ReciprocalModel.source` is the integrated action
source, not the acceleration source. Its normalization and sign must agree
with an independent integral of the full improved spatial stress at a=1.

For beta=A(nn^T-I/3), P=I-nn^T, v=n_dot, the normal residual is

    R_N=N_n(S)-2A(vv^T-|v|^2 P/2),
    N_n(M)=PMP-tr(PMP)P/2.

For a constant A and great-circle speed Omega, the required source has
radial, angular and normal components

    A(omega_T^2-3 Omega^2)(nn^T-I/3),
    0,
    2A(vv^T-Omega^2 P/2).

All three components and the full five-component equation must agree over
an interval. An instantaneous equality, equality of norms only, or forced
projection back to the uniaxial cone is insufficient.

## P2. An analytic frequency restriction, with its exception

A general free degree-3 scalar is q=u cos(omega_s t)+w sin(omega_s t).
Its quadratic source has only frequencies 0 and 2 omega_s:

    S=S_0+S_c cos(2 omega_s t)+S_s sin(2 omega_s t).

A uniform great-circle rotor has beta=B_0+B_c cos(2 Omega t)
+B_s sin(2 Omega t), with B_c,B_s nonzero for A!=0. Its required source is

    omega_T^2 B_0
    +(omega_T^2-4 Omega^2)[B_c cos(2 Omega t)+B_s sin(2 Omega t)].

Linear independence of these frequencies on an open interval implies
Omega=omega_s, unless omega_T^2-4 Omega^2=0. At the exception
Omega=omega_T/2 a constant source can sustain the motion; it need not have
the scalar's oscillation frequency. This is a necessary frequency condition,
not an existence proof in either branch. In particular the original speed
0.4/a is excluded for this uniform, constant-amplitude, single-multiplet
leading-order class. Nonuniform rotation and varying A are outside this
frequency theorem. Record this as an analytic prediction, not a fitted scan.

## P3. Explicit constant-source candidate, before measurement

On unit embedding coordinates x=(x0,x1,x2,x3), define Haar-normalized
harmonic polynomials

    f=sqrt(8) Re[(x0+i x3)^3],
    g=sqrt(8) Re[(x1+i x2)^3],
    phi(t,x)=s [f(x) cos(omega_s t)+g(x) sin(omega_s t)]/sqrt(V).

The predicted identities in the inherited invariant frame are

    <f^2>=<g^2>=1, <fg>=0,
    <D_i f D_j f>=<D_i g D_j g>=diag(3,3,9),
    <D_(i f D_j) g>=0.

Here D is the dimensionless unit-sphere derivative, and <> is normalized
Haar. The cross identity is for the symmetric tensor. Prove or refute these
by exact monomial moments; an eigensolver's floating-point coefficients are
not an exact certificate.

If they hold, the field supplies the constant acceleration source

    S=k Q_z,  Q_z=e_z e_z^T-I/3,  k=6 s^2/(C a^2).

The proposed tensor solution is

    Omega=omega_T/2=sqrt(2)/a,
    A=-2k/omega_T^2=-3s^2/(2C),
    n=(cos(Omega t),sin(Omega t),0),  beta=A Q_n.

A is signed: this is the negative-amplitude, oblate branch of the existing
uniaxial cone. Deleting that branch would change the original field family.
Use beta(0)=A Q_x and beta_dot(0)=A Omega(e_x e_y^T+e_y e_x^T).
These are **chosen nonzero free tensor data**. They differ from #290's zero
tensor preparation; #290's initial no-cancellation certificate is not
transferred to this preparation. The source itself is obtained from the
scalar; the tensor orbit is not claimed to be an attractor or selected by
the scalar for arbitrary initial data.

The joint scalar/tensor history need not close: their predicted frequencies
have irrational ratio omega_s/Omega=2 sqrt(2). A constant source and a
periodic tensor do not establish a periodic full history or a triangle map.

## P4. Independent constraint compatibility

For this odd multiplet, the Hamiltonian dipole is zero at O(s^2). That says
nothing by itself about all momentum charges. Construct all six ambient
rotations x_A partial_B-x_B partial_A and all four gradient dipoles grad x_A.
Evaluate their charges from the full improved j_i=-T_0i and derive the
candidate's vanishing charges analytically. The polynomial separation of f
and g suggests all cross Killing overlaps <g K f> vanish; verify this rather
than substituting only the inherited three invariant generators.

Negative controls must include (i) parity-pure data with vanishing three
invariant-generator charges but a nonzero omitted rotation charge, and
(ii) phi in degree 1, phidot in degree 2, with zero Hamiltonian dipole and
rotation charges but a nonzero gradient conformal-Killing charge. Use the
complete ten-charge condition to reject them. These do not modify #291.

## Frozen implementation and gates

Primary a=kappa=1, s=0.02, over one tensor period pi/Omega at 401 times.
Controls: s=0.01,0.005; radii 0.7,2 with kappa=0.4 for source and tensor
normalization. Improved stress hardcodes a=1, so use it only there and
check general-radius scaling through the independently derived action.
Spatial rules (8,16) and (12,24); DOP853 at (rtol,atol)=(1e-10,1e-12)
and (1e-12,1e-14). Define scaled error as Euclidean error divided by
max(1, reference norm); divide forces by k and tensors by |A| first.

Required gates, all mandatory and named in a canonical list:

1. Exact polynomial norms, eigenvalues, gradient moments and cross-source
   identities; explicit zero residuals using rational monomial moments.
2. Action source versus full improved stress on both spatial rules at
   0, pi/(8 omega_s), pi/(4 omega_s), pi/(2 omega_s), pi/omega_s;
   normalized tensor error <1e-9. Keep nonzero scalar momentum times.
3. Hamiltonian dipole and all ten momentum charges on both rules, divided
   by s^2, <1e-9 at the same times; exact candidate charge certificate.
4. Both omitted-charge negative controls are rejected, with the missed
   charge >1e-3 after unit modal normalization and the purportedly zero
   charges <1e-9. Constraint gating must not accept three charges as ten.
5. Constant source, radial/angular/normal balance and full tensor equation
   over the 401 times, normalized residuals <1e-9.
6. Independent integration of all five tensor components with the freely
   evolved scalar source (no uniaxial projection); tensor trajectory and
   distance to the full uniaxial cone, divided by |A|, <1e-8; ODE refinement
   <1e-8. Check nonzero director motion from the tensor, modulo n~-n.
7. Necessity controls: zero source, source sign reversal, and a doubled
   initial director speed violate the full equation by >1e-3 in the fixed
   normalized units. Also perturb the scalar's second-quadrature amplitude
   by 10 percent while retaining the proposed orbit: failure expected.
8. Source and tensor amplitude-halving ratios are four, to 1e-9; radius
   controls satisfy the predicted normal/full balance to 1e-9. Maximum
   Frobenius norm of the tensor <0.05 for all frozen demonstrations.
9. Exact frequency restriction including the constant-source exception;
   do not infer that either admissible frequency alone proves feasibility.
10. Removing or failing any required gate renders every physical verdict
    UNRESOLVED, names failures, overwrites stale JSON/Markdown reports and
    exits nonzero. A failed candidate is not a no-go for other scalar data.

A successful candidate may yield `LEADING_ORDER_FIELD_SUPPORTED_ROTOR` and
`CHOSEN_CONSTRAINT_COMPATIBLE_PREPARATION`. Failure yields `UNRESOLVED` on
existence, not `NO_ROTOR`. Keep separate verdict fields for normal balance,
constraint compatibility, homogeneous tensor evolution, full Einstein-matter
evolution (`NOT_ESTABLISHED`), preparation selection (`NOT_DERIVED`),
triangle map (`NOT_DERIVED`) and Phi selection (`NOT_DERIVED`). Preserve this
freeze and all previous rounds; archive exact certificates and measured
checks, clearly distinguishing each from a premise.
