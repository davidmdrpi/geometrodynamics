# Nonlinear completion of supported tensor histories

This round tests whether #297–#298's homogeneous linear tensor histories are
tangents to actual solutions of the four-field Einstein equations. The
scientific freeze is `a067782`; the pre-implementation provenance addendum is
`9f922b6`. See [the freeze](nonlinear_supported_tt_prereg.md) and
[the review-provenance record](nonlinear_supported_tt_provenance.md).
No frozen prediction, grid, or threshold has been changed.

**Result:** the quaternionic four-field system closes exactly, and all five
homogeneous TT displacement and velocity directions admit nearby nonlinear
constraint-satisfying histories. An explicit continuation argument gives
future persistence on a stated domain and, by continuous dependence, for
sufficiently small perturbations of each departing FRW reference.

The original numerical freeze passed **11/13** gates: its .01 central
amplitude differences missed the 1e-3 accuracy cutoff. Its N/F machine
verdicts remain UNRESOLVED. A separately published prospective extension,
[`75c7ab4`](nonlinear_supported_tt_refinement_prereg.md), specified .005 and
.0025 before they were run. That refinement passes **13/13**, with measured
second-order convergence and all the original non-amplitude gates retained.
The original failure is not relabeled as a pass.

## 1. Exact four-field closure

The original four real conformal scalars are retained, with
`Lambda=3/(2a^2)`. In the repository's invariant coframe write

    g=A^2[-n^2 d eta^2+M_ij(e^i+N^i d eta)(e^j+N^j d eta)],
    M=M^T>0, det M=1,
    phi=B(q)x/A, B(q)=q0 I+sum qi S_i.

The same S_i define the coframe `e^i=(S_i x).dx`. Their anticommutators are
`S_i S_j+S_j S_i=-2 delta_ij I`. Consequently `B^T B=Q I`, where Q=q.q,
and every component of phi remains antipodally odd. The opposite quaternion
algebra does not give homogeneous currents in this coframe; it is a
negative control, not an interchangeable convention.

For the homogeneous volume form every invariant frame vector is
divergence-free. Thus `Delta_M x=-tr(M^-1)x`, directly from the Clifford
identity. The summed field square, spatial gradient tensor, kinetic term
and momentum current are homogeneous. The curvature scalar is homogeneous
as well. Each scalar equation therefore stays within this four-dimensional
quaternionic coefficient space. This is an exact equivariant reduction,
not a projection of a spatially varying stress onto its mean.

At zero shift and conformal lapse one define

    H=A^2/kappa-Q/6, F=H/A^2,
    L=M^-1 M'/2, tr L=0,
    r=2[2 tr(M^-1)-tr(M^2)], t_M=tr(M^-1), ell=tr(L^2).

L is self-adjoint in M: ML=(ML)^T. It need not be symmetric in the fixed
coframe. A prime is a conformal-time derivative throughout. The exact
reduced action per unit-S3 volume is

    L_red=-3A'^2/kappa+q'.q'/2+H tr(L^2)/2
          +H r/2-Q tr(M^-1)/2-Lambda A^4/kappa.

Conformal invariance of the scalar action removes the q/A cross terms;
the Einstein curvature boundary term gives the displayed scale kinetic
term. The conformal spatial volume is constant. Its curvature scalar is
`r+tr(L^2)`, with no trace-extrinsic-curvature boundary contribution.
The spatial curvature follows from the diagonal invariant metric and
rotation covariance, so its formula applies to every positive symmetric M.

Lapse variation is performed before gauge fixing: the kinetic terms have
factor 1/n and the potential terms have factor n. For the shift, let
`D_N=2 cross(N)`. The covariant velocities are

    M'_cov=M'-[D_N,M], B'_cov=B'-B sum N^i S_i.

These follow from the repository's frame brackets. They produce all three
momentum equations. Choosing a diagonal metric at the outset would have
hidden those equations.

## 2. Evolution, constraints and full-field verification

The resulting equations are

    A''=-A(r+ell)/6+2 Lambda A^3/3,
    q''=-[t_M+(r+ell)/6]q,
    M'=2ML,
    (HL)'=H STF[(-4+Q/H)M^-1-4M^2].

The constraints are the energy

    E=-3A'^2/kappa+q'.q'/2+H ell/2-H r/2+Q t_M/2
      +Lambda A^4/kappa=0

and three momentum components

    C_i=H tr[(L-L^T) cross(e_i)]-j_i=0,
    j_i=tr[B(q')^T B(q) S_i]/4.

Both E and C_i are conserved by the conformal equations. For example,
`r'=-8 tr[(M^-1+M^2)L]`, `t_M'=-2 tr(M^-1 L)`, and the q and L equations
cancel every term in E'. For momentum, the force in `(HL)'` is symmetric,
so the antisymmetric part of HL is constant. The scalar equation has one
common real coefficient for all four q components, making j_i constant.
Thus completed initial constraints propagate; they are not repeatedly
projected back to zero during integration.

The independent route constructs full coordinate metric jets at finite M,
Christoffel symbols, Ricci curvature, all four improved stresses, and all
four scalar equations. It uses neither the reduced stress nor a TT-only
projection. In the frame (A d eta,A e^i), define residuals of the scale,
scalar and shape evolution equations as R_A,R_q,R_L respectively. The
independent full Einstein residual has

    Einstein_residual^0_0=kappa E/A^4,
    Einstein_residual^0_i=-kappa C_i/A^4,
    Einstein_residual^i_0=kappa (M^-1 C)^i/A^4,
    STF(Einstein_residual^i_j)=kappa R_L/A^4,
    tr(Einstein_residual)=-6 R_A/A^3+kappa q.R_q/A^4.

The scalar residual is `-B(R_q)x/A^3`. These identities are checked off
shell, including independent accelerations, and on shell at multiple
spatial points. Proper-time jets transform every derivative before the
curvature is computed. They independently reproduce the conformal route.

The kinetic Legendre map is regular for A>0,H>0,M>0. Before imposing the
constraints its canonical one-form is

    p_A dA+q'.dq+tr(Pi dM),
    p_A=-6A'/kappa, Pi=(H/2)L M^-1,

with the determinant constraint on M understood. The action has the fixed
overall factor `Vol(S3_unit)=2pi^2`. Symplectic conservation is a statement
on the complete constrained/gauge-reduced phase space. No isolated nonlinear
tensor submap is declared symplectic after eliminating the matter response.

## 3. Closed-form initial-data completion

For each prescribed STF pair U,V and amplitude epsilon set

    R=exp(epsilon U), M=R^2, L=epsilon R^-1 V R.

This realizes exactly the frozen M and M' initial data. Put

    g_i=tr[(L-L^T) cross(e_i)],
    D_b=q_b^2+q_b'^2, H_b=A_b^2/kappa-q_b^2/6,
    c=q_b'^2/(6D_b^2).

The freeze's completion ansatz gives the **exact** current j_i=xi_i and
`H=H_b-c|xi|^2`. The momentum equations consequently reduce to

    xi=H g, H+c|g|^2 H^2=H_b.

The positive root continuous from the reference solution is

    H=2H_b/[1+sqrt(1+4cH_b|g|^2)], xi=H g.

It is regular at g=0, at q_b=0, and at q_b'=0. The remaining Hamiltonian
constraint fixes the positive expanding A' root when its radicand is
positive. Near any frozen reference with d>0 it is positive by continuity.
The Jacobian in variables (A',xi) at epsilon=0 is exactly

    diag(-6A_b'/kappa,-1,-1,-1).

Its determinant is `6A_b'/kappa>0`, independent of matter phase. The implicit
function theorem therefore gives local nonlinear initial-data families for
**all five homogeneous TT displacement and velocity directions**, including
noncommuting pairs. This is an existence statement, not dimension counting
or an interpretation of a small Newton residual.

The leading momentum source is

    g_i=-2epsilon^2 tr([U,V] cross(e_i))+O(epsilon^3).

The responsive matter changes at second order, exactly as frozen. The
rigid-support control has xi=0 and fails for the noncommuting examples. Its
failure cannot be elevated to an obstruction for the original four fields.
A spatial Killing field alone is not a symmetry of each individual
background scalar. Omitting their allowed response would create a false
integrability obstruction.

## 4. Independent second-order response

Let `beta=epsilon b+epsilon^2 c2+...`,
`A=A_b+epsilon^2 a2+...`, and
`q=(q_b+epsilon^2 u2,epsilon^2 v2)+...`, initially at fixed conformal time.
Set h=H_b'/H_b and k=K_b/H_b, with `K_b=8H_b+2q_b^2`. The first variation is
precisely #297's equation `b''+h b'+k b=0`. The independently derived
second-variation equations are

    a2''=(-1+3A_b^2/a^2)a2
          -A_b[tr(b'^2)-8tr(b^2)]/6,
    u2''+4u2=-q_b[(2/3)tr(b^2)+tr(b'^2)/6],
    v2''+4v2=0,
    c2''+h c2'+k c2=(-40+2q_b^2/H_b) STF(b^2).

The initial scalar-amplitude displacement and velocity corrections vanish.
The transverse initial v2 and v2' are the frozen xi prescription with
`xi2=H_b g2`. The Hamiltonian correction is

    a2'(0)=kappa[H_b tr(V^2)+K_b tr(U^2)]/(12 A_b'),
    a2(0)=c2(0)=c2'(0)=0.

The quadratic tensor force is retained as well as the scale and support
responses. The noncommuting part of `L2=c2'-[b,b']` accounts for the momentum
source, rather than being dropped by replacing L with a symmetric matrix.

For equal proper elapsed time, let `T2'=a2`, `T2(0)=0`. Then
`eta2=-T2/A_b`. The scale and common-scalar quadratic responses become
`a2+A_b' eta2` and `u2+q_b' eta2`. Transverse q and the quadratic tensor
response receive no second-order time-shift term because their backgrounds
vanish. Central finite-amplitude differences test these converted equations.
The nonlinear endpoint is not fixed to #298's reference eta_star.

## 5. Future continuation with explicit bootstrap margins

Use dimensionless `A/a`, `sqrt(kappa)q/a`, and t/a in this section, writing
them as A,q,t. Thus Lambda=3/2 and kappa=1. Consider admissible expanding
initial data at any A_0>=32 with

    ||M_0-I||_F <= 1/10, ||L_0||_F <= 1/10,
    |q_0|+|q_0'| <= 4, det M_0=1, M_0 L_0 symmetric,
    E=C_i=0.

Bootstrap the larger region

    ||M-I||_F <= 1/5, ||L||_F <= 1, |q|+|q'| <= 8.

Here M has eigenvalues in [4/5,6/5]. Elementary matrix estimates give
`|r|<=15`, `tr(M^-1)<=4`, and `0<=ell<=1`; the lower bound on ell uses
M-self-adjointness. With h_0=49/50 and c_0=2/3,

    h_0 A^2 <= H <= A^2,
    c_0 A^2 <= A' <= A^2.

The first follows from Q<=64. For the second, divide the Hamiltonian
constraint by 3A^4:

    (A'/A^2)^2 = 1/2 + (|q'|^2 + H ell - H r + Q tr(M^-1))/(6A^4).

The curvature contribution has absolute value at most 15/(6A^2).
Retaining the nonnegative kinetic and Q terms gives the lower estimate.
For the upper estimate, use H ell<=A^2, |H r|<=15A^2, and
`|q'|^2+Q tr(M^-1)<=64+256=320`. Thus

    1/2 - 15/(6A^2) <= (A'/A^2)^2
      <= 1/2 + 16/(6A^2) + 320/(6A^4).

The conservative coefficient 16 is 1+15, including the absolute curvature
bound as well as shear; dropping curvature would require its positive sign.
The sign A'>0 cannot change under the lower bound.

The remaining conformal duration is at most `1/(c_0 A_0)<=3/64`. The scalar
frequency has absolute value at most 7 in its second-order equation, hence
`Y=|q|+|q'|` obeys the upper Dini derivative inequality `Y'<=7Y`. Since
`exp(z)<=1/(1-z)` for 0<=z<1,

    Y <= 4/[1-7/(c_0 A_0)] <= 256/43 < 8.

For delta=1/5, the shape force satisfies

    ||STF[(-4+Q/H)M^-1-4M^2]||_F <= 14 delta.

This uses `||M^-1-I||_F<=(5/4)||M-I||_F`,
`||M^2-I||_F<=(11/5)||M-I||_F`, and the contraction of STF in Frobenius norm.
Integrating the exact equation for HL gives

    ||L(A)||_F <= A_0^2 ||L_0||_F/(h_0 A^2)
                 +14delta(A-A_0)/(h_0 c_0 A^2),
    integral ||L||_F d eta
      <= ||L_0||_F/(3h_0 c_0 A_0)
         +14delta/(6h_0 c_0^2 A_0^2)
      <= 265/100352.

Consequently

    ||M-I||_F <= 1/10+2(1+delta)*265/100352 < 1/5,
    ||L||_F <= (1/10)/h_0+14delta/(4h_0 c_0 A_0) < 1.

Every bootstrap boundary is strictly improved. The archive's exact rational
certificate records each margin, including the force Lipschitz estimate.
On any finite-A interval the equations remain in a compact regular chart,
so ordinary continuation applies. A diverges only at the finite conformal
endpoint, and

    dt/dA=A/A' >= 1/A

shows that it takes **infinite proper time**. The above estimates also yield
`L=O(A^-1)` and integrable shape velocity. M therefore converges to a positive
M_infinity, physical shear `L/A=O(A^-2)`, and physical fields `B(q)x/A=O(A^-1)`.
The Hamiltonian constraint gives `A'/A^2 -> 1/sqrt(2)`; restoring units this
is the proper Hubble limit `sqrt(Lambda/3)`. Also kappa F tends to one.
This is a closed continuation argument, not an assumption that F stays
positive or an extrapolation of numerical tail slopes.

Every reference expanding FRW history with fixed d>0 reaches this tail
region strictly: M=I, L=0 and `|q|+|q'|<=sqrt(15)/2<4`. Local existence and
continuous dependence on the compact reference interval imply an open
neighborhood of its constraint-completed initial data reaches the interior
of the tail region as well. Together with the explicit constraint Jacobian,
this establishes nonlinear histories for sufficiently small tensor data in
all five directions, existing for all future proper time in this sector.

The smallness condition is explicit on the A>=32 entry surface. We do not
compute a numerical epsilon radius for its pullback to each initial d and
phase, and numerical entry by a sampled trajectory is not a rigorous bound
for a whole initial-amplitude grid. Compactness permits a common sufficiently
small neighborhood over a compact d>0/phase set; it does not establish
persistence arbitrarily close to d=0 with a uniform radius. This result does
not remove the ESU instability or establish inhomogeneous stability.

## 6. Numerical results, failure gates and scope

The pre-review targeted suite reported **214 passed**, including 61 new tests. It
covers the nonzero shift variation, independent off-shell accelerations,
field/velocity-zero completion, phase and unit changes, the prospective
refinement, source-hash provenance, every missing/false gate, malformed raw
evidence, and actual CLI failures overwriting stale success artifacts.
The review update's milestone suite reports **67 passed**, including six
new constraint-tampering cases at both refined amplitudes and three times.

| Check | Largest observed error / outcome |
|---|---:|
| all 3,240 initial constraint solves | 1.19e-16 normalized |
| full coordinate field agreement, 378 on-shell and 20 off-shell cases | 1.79e-15 |
| 108 scale/coupling/clock controls | 3.07e-15 |
| constraint propagation on 213 trajectories | 9.16e-14 normalized |
| constraint propagation on 812 smaller-amplitude refinement states | 5.38e-15 normalized |
| final solver comparisons, each state block | 7.38e-15 relative |
| independent linear-operator/FRW controls | 2.33e-14 |
| original first-variation error, epsilon=.01 | 1.176e-3: fails 1e-3 |
| original second-variation error, epsilon=.01 | 1.099e-3: fails 1e-3 |
| refined first-variation error, epsilon=.0025 | 7.37e-5 |
| refined second-variation error, epsilon=.0025 | 6.89e-5 |
| refined consecutive error ratios | 3.9927 to 4.0000 |

The rigid noncommuting control has momentum residual 0.0312077; allowing
the frozen matter response reduces it to numerical zero. Both coordinate
clocks agree. The frozen original run still exits with status 1 because its
two amplitude-difference gates fail; the prospective refinement exits 0.
A test run must not suppress the original nonzero exit code and describe
it as a fully passing original freeze.

The raw [original run](../experiments/closure_ledger/runs/20260913_nonlinear_supported_tt/probe.json.gz)
and [original gate report](../experiments/closure_ledger/runs/20260913_nonlinear_supported_tt/probe.md)
are preserved separately from the
[refinement run](../experiments/closure_ledger/runs/20260913_nonlinear_supported_tt/refinement.json.gz)
and [refinement verdict](../experiments/closure_ledger/runs/20260913_nonlinear_supported_tt/refinement.md).
The JSON runs are losslessly gzip-compressed for repository distribution;
their decompressed bytes reproduce exactly, including the original hash
recorded before refinement. Probe commands still write plain JSON, and the
refinement command accepts either plain or gzip-compressed original input.
The pre-refinement source snapshots accompany the archive so the hashes in
the extension freeze remain traceable after validation code was strengthened.
The post-freeze source changes added input validation, defensive copies of
cached certificates, and raw-evidence gate checks; the evolution equations,
integrator, tolerances and archived trajectories are unchanged. Re-scoring
the frozen original evidence with the current code reproduces all 13 saved
gate values and the complete original verdict exactly: still 11/13, with
N/F unresolved. This equality is also checked in the test suite.

The review update additionally checks Hamiltonian and all three momentum
constraints on every refinement state, including all 812 states at epsilon
.005/.0025, using the original normalized 1e-8 propagation cutoff. This
strengthens the existing `constraint_propagation` gate for the extension;
it does not alter the original report. A velocity-only corruption invisible
to the variation observables now fails that gate and withdraws N/F alone.
Both refinement accuracy and ratio checks use the registered `0<t/a<=2`
window; constraint checks include t=0 and the later samples through t/a=8.

The deterministic archive and gate report accompany this document. All
frozen cases, failed gates, raw constraints, coordinate-field residuals,
trajectories, solver comparisons and variation errors are retained. The
scientific status follows the frozen dependency table, including a separate
quadratic-response gate. No observed failure is converted to an obstruction
or erased by changing a threshold.

Reproduce with:

```sh
python -m experiments.closure_ledger.nonlinear_supported_tt_probe \
  --output-dir experiments/closure_ledger/runs/20260913_nonlinear_supported_tt
# The original command exits 1: retain its failed verdict.
python -m experiments.closure_ledger.nonlinear_supported_tt_refinement_probe \
  --original experiments/closure_ledger/runs/20260913_nonlinear_supported_tt/probe.json \
  --output-dir experiments/closure_ledger/runs/20260913_nonlinear_supported_tt
python -m pytest -q tests/test_nonlinear_supported_tt.py
```

This result concerns the homogeneous quaternionic four-field sector. It
adds no new matter fields but it does allow a new configuration of the four
existing fields. It does not select that preparation, construct a driven
rotor, derive Phi or quantization, or resolve operational causality.
