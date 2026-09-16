# Localized four-scalar bulk/handle initial data: prospective specification

Date: 2026-09-16. Scientific baseline: #304, `7f2acb3`, with its Python 3.10
compatibility correction recorded separately. #300 and the existing odd
four-scalar support supply the action. No new experiment in this document
has been implemented or run before publication. The algebra below consists
of prospective analytic predictions, to be verified.

## Milestone and exclusions

Construct a smooth compact handle joined to a region retaining the existing
round-S3 scalar profile, with NONZERO scalar momentum and solved metric and
extrinsic curvature. Quantify the bulk's geometric departure from round S3
as the gluing is moved into smaller antipodal caps. This tests the missing
localized connection, not the vacuum product geometry of #304 again.

The four existing real conformal scalars are retained. They have not been
derived from pure vacuum geometry. No additional matter, supporting shell,
constitutive law, detector, kick or probability rule is supplied. Their
collar profile and initial normal momenta below are chosen free initial
data, not dynamically selected preparations or a formation mechanism.
Initial constraints are the scope: no crossing evolution can be claimed
from their success. Failure to construct the desired branch is a recorded
outcome, not a general impossibility theorem.

## The exact action and field transport

Set kappa=a=1, Lambda=3/2 as in the radius-one four-scalar support:

    I=int sqrt(-g)[(R-2Lambda)/2 - (1/2)sum (dphi)^2
                                              -R sum phi^2/12].

Use g_E=f g_J, f=1-sum(phi^2)/6. The Einstein-frame field-space metric and
potential, retaining all four fields, are

    G_AB=delta_AB/f+phi_A phi_B/(6f^2), U=Lambda/f^2.

Work on -L<=s<=L, n in S2 with (L,n)~(-L,-n), conformal metric
gbar=ds^2+dOmega^2. Cutting two antipodal caps out of S3 gives the cylindrical
bulk chart sin(chi)=sech(s), cos(chi)=tanh(s). The collar replaces the caps'
connection. Geometry must be solved, not declared round by this chart.

Impose the field transition phi(L,n)=-phi(-L,-n), and the same transition
for normal momenta. This is an explicit Z2 line-bundle sector of the original
sign-symmetric action. It is an extra global sector choice: it is not an
ordinary single-valued scalar on the quotient, is not derived by the
constraints, and must not be hidden in a sign convention. All local stresses
and action densities are single-valued. The timelike normal is unchanged.

## Smooth bulk-to-collar fields

Put q=sqrt(3)/2 and

    phi=q (sin theta(s), cos theta(s) n),
    d_phi/d_theta=q (cos theta(s), -sin theta(s) n),
    tilde_Pi=p(s) d_phi/d_theta, Pi_E=psi^-6 tilde_Pi.

Then sum phi^2=3/4, f=7/8, and n_E f=0 on the initial slice. This is an
initial-data restriction, not a constraint on subsequent field evolution.
The full four-field dynamics may change the norm after this slice.

Let theta be odd and, for s>=0,

    theta(0)=0, theta'=sech(s) c(s),
    c=1 for s<=L-1, c=0 for s>=L-1/2.

In the transition set z=2(s-L+1), b(z)=exp(-1/z) for z>0 and zero otherwise,
c=b(1-z)/(b(z)+b(1-z)). Thus the fields agree EXACTLY with the old round-S3
profile in |s|<=L-1, and are smooth and constant in theta near the seam.
Compute theta by high-accuracy integration, retaining the analytic theta'
for the equations. Do not approximate the smooth cutoff by a discontinuity.

Choose antiperiodic normal-momentum coefficient and throat integration data

    p(s)=eta sech(L)^3 sin(pi s/(2L)),
    C=eta sech(L)^3/8.

These make the physical momentum scale comparable as the caps shrink; they
are a registered preparation family, not a prediction of quantization.

## Coupled constraints

Take gamma_E=psi^4 gbar and K_E=psi^-2 diag(2a,-a Omega_AB), tr K_E=0.
For alpha=q^2/f=6/7, the predicted momentum equation is

    2a'=-alpha p theta', a(L)=C.

Its source is odd on the double interval and has zero integral, as required
by the radial conformal Killing field. A deliberately even momentum p is
the incompatible control; projecting away its mean silently is prohibited.

The predicted Einstein-frame Hamiltonian equation is

    psi'' - [2-alpha(theta'^2+2cos(theta)^2)] psi/8
           +[6a^2+alpha p^2] psi^-7/8 + U psi^5/4 = 0.

Solve on [0,L] with psi'(0)=psi'(L)=0, psi>0. Reflection supplies the other
half. The nonlinear elliptic solve, not the choice of conformal chart,
determines the neck and bulk metric. The exact uncut reference at eta=0 is

    theta_ref=asin(tanh(s)), psi_ref=f^(1/4) sqrt(sech(s)), a=p=0.

Verify that reference by independent substitution. It fails the finite-L
Neumann seam condition; treating it as an already glued solution is a
negative control.

## Independent physical/Jordan-frame verification

On this slice f is spatially constant and n f=0, hence

    gamma_J=gamma_E/f, K_J=K_E/sqrt(f), Pi_J=sqrt(f) Pi_E.

Derive the ADM projections of the improved stress from the ORIGINAL action.
The predicted physical constraints, valid for these constant-norm/tangent
data only, reduce to

    f [R_J+(tr K_J)^2-|K_J|^2] = |Pi_J|^2+|Dphi|_J^2+2Lambda,
    f [D_j K_J^j_i-D_i tr K_J] = -sum Pi_J d_i phi.

Check them from coordinate metric/field/K finite differences at off-grid
points, independently of the conformal ODE contractions. Retain both local
frames, all four components, scalar momentum and gradients in raw evidence.
Discarding G_AB's nonminimal factors is the wrong-theory control.

Report r_J(s)=psi^2/sqrt(f), section areas 4pi r_J^2, and relative bulk metric
factor r_J/sech(s) on |s|<=1. For a minimal seam compute future null expansions
with K=-L_n gamma/2. Report their signs. A marginal or anti-trapped seam is
not a traversability result. In the regular f>0 Einstein frame,
R_kk=G_AB(k.dphi_A)(k.dphi_B)>=0; checking this null-convergence identity
does not prove every global traversability claim false, but it prevents
inventing a negative-energy support channel in this action.

## Frozen schedule, controls and acceptance criteria

L=3.5,4.5,5.5; eta=0,.1,.3; all nine cases, no amplitude cherry-picking.
Use solve_bvp or an equivalent collocation BVP method on initial grids
129,257,513, tolerances 1e-6,1e-8,1e-10 respectively, maximum 40000 nodes.
Use the positive reflected reference guess

    psi_guess=f^(1/4)[sqrt(sech(s))+sqrt(sech(2L-s))],

and its analytic derivative. Integrate theta and a with relative/absolute
tolerances 1e-12/1e-14. Each case starts from the stated guess, with no
unrecorded branch search. Record convergence failures and nonpositive
solutions as failures; do not clamp psi or relax tolerance. A new numerical
branch search requires a prospective addendum. Archive BVP polynomial
coefficients so off-grid equations can be rechecked without rerunning.

Off-grid scalar equation check: 1001 points on [0,L], omit no collar points;
use derivatives of the stored collocation polynomial. Mesh-refinement
comparison: common 1001 points. Physical finite-difference checks at the
finest (L=5.5,eta=.3) data, s/L=(.07,.23,.51,.79,.93) and
theta_angle=(.43,.91,1.47,2.13), phi_angle=.37; h=1e-3,5e-4,2.5e-4.
Also evaluate the seam directly from the one-sided smooth data.

Gates:

1. Exact reference, Einstein-frame kinetic and constraint-projection
   identities agree symbolically; f=7/8 and kinetic eigenvalues positive.
2. All momentum solves meet max off-grid residual 1e-9. Nonzero eta gives
   a nonzero correction a-C and nonzero scalar current. The even-p control's
   compact-source integral has absolute value >1e-7 at L=4.5,eta=.3 and is
   rejected. An omitted correction fails the momentum equation by >1e-7.
3. All 27 BVP solves converge with psi>0 and boundary derivative residual
   <1e-9. Finest off-grid Hamiltonian residual <1e-7 and medium-to-fine
   max relative psi difference <1e-6. Any miss is retained.
4. Metric/K and sign-bundle scalar/momentum seam value and first-derivative
   mismatch <1e-7. Omitting the scalar transition sign gives mismatch >.1.
5. Independent Jordan-frame Hamiltonian and momentum normalized residuals
   <1e-5 at finest h. Last refinement ratios in [2.5,5.5] when preceding
   error exceeds 1e-8. The wrong-theory control omitting f has a physical
   residual >10 times the correct maximum.
6. Localization: for each eta, maximum |r_J/sech(s)-1| on |s|<=1 decreases
   as L increases, and is <.1 at L=5.5. Require 0<r_J(L)/r_J(0)<.15 at
   L=5.5, and a genuine local section minimum: r_J(L-.05)>r_J(L),
   r_J(L-.1)>r_J(L). Failure leaves valid global data distinct from a
   localized-mouth claim.
7. At L=4.5,eta=-.3, the medium-grid metric agrees with eta=.3 to 1e-8,
   K and Pi reverse signs, and null expansions transform consistently.
   Numerical null contractions of G stay nonnegative. All are diagnostics
   of the chosen action, not positive traversability gates.
8. Evidence integrity: all cases and physical source data required; missing,
   duplicated, nonfinite, corrupted, mismatched-frame or wrong-sign data
   must withdraw the dependent affirmative verdict. Recompute verdicts
   from saved data; stale output must not survive a failed CLI run.

`FOUR_SCALAR_HANDLE_CONSTRAINT_DATA` requires 1-5,7,8.
`LOCALIZED_BULK_MOUTH_INITIAL_DATA` additionally requires 6.
An evolved crossing, reciprocal momentum-transfer event, traversability,
topology/sector selection, discrete action and quantum statistics remain
UNESTABLISHED regardless of those gates. There is no threshold or event
detector in this milestone.

## Sources and deliverables

Inherited action and kinetic matrix: `docs/odd_multiplet_support.md`,
`waves/odd_multiplet_support.py`, `waves/scalar_esu_support.py`.
Constraint benchmark: #304's `waves/mouth_momentum.py`.
For the standard conformal method and matter scaling see
https://arxiv.org/abs/gr-qc/0610045 and https://arxiv.org/abs/2106.15027;
for multifield conformal transformations see
https://doi.org/10.1103/PhysRevD.81.084044.
The calculation here must derive its own nonminimal factors, rather than
apply a minimally coupled scalar formula in the Jordan frame.

Deliver a reproducible module/probe, full raw evidence, tests and a results
document. Publish this freeze before any of its experimental checks. Keep
initial-data existence distinct from a physical crossing mechanism.
