# Follow-up freeze: bound the scalar–TT model's omitted constraints

Baseline: PR #289 at `67940cb32d26ffed0c24634348ed9896d3964d64`, following
the COMMENT review against that commit. Commit and publish this document
before implementing or measuring the new response. Preserve the original
reciprocal-action freeze and archive. Seed: `2026090712`.

## Question, inherited assumptions, and exclusions

How large a metric response do the omitted Einstein constraints require for
the existing smooth degree-3 scalar, and does that bound the omitted scalar
backreaction? These are separate questions. First solve the **linearized
constraints on individual round-ESU slices**. A sequence of such solutions
is not an Einstein–matter evolution.

Use the inherited conformal scalar stress, radius a, kappa = 8 pi G, and
smooth whole S3. `waves/initial_data.py` already derives

    (Delta + 3/a^2) u = -kappa delta_rho / 4,
    g_ij = (1 + 4u) gbar_ij + 2 beta_ij + higher orders

(the tensor term is expressed in the background orthonormal frame).
Hold the supporting matter's normal-frame delta_rho and delta_j at zero
**on each independently considered slice** for the scalar-only particular
solution. This is an explicit assumption inherited from the initial-data
channel, not a consequence of separate scalar stress conservation and not
a claim that the support can stay rigid during evolution. Keep the spatial
mean; do not silently retune the background or cosmological constant.

For momentum choose constant mean curvature, delta K = 0, with the ADM
convention K_ij = -1/2 Lie_normal g_ij and j_i = -T_0i. Choose the minimal
longitudinal particular solution. Free TT data are separate. Spatial
conformal gauge and kernel coefficients are fixed explicitly below; lapse,
shift, matter equation of state, and an evolution of the added metric are
not selected. No pointer experiment or probability law follows here.

## Analytic predictions to verify, not numerical discoveries

### Hamiltonian constraint and a finite bound

At the retained order Kbar = 0 makes the extrinsic-curvature terms quadratic.
Delta Y_l = -lambda_l Y_l, lambda_l = l(l+2)/a^2. The constraint operator
has a four-dimensional l=1 kernel. A single degree-3 multiplet has odd phi
and phi_dot; its quadratic energy is even with degrees 0,2,4,6. Thus the
dipole solvability condition holds exactly. Set homogeneous kernel
coefficients to zero. Smooth even sources need no mouth excision to invert
this operator; that does not cure the singular-mouth problem of initial_data.

For rho_l the orthonormal harmonic coefficients, including the mean,

    u_l = -kappa rho_l / [4(3/a^2-lambda_l)].

With either physical L2 norms or volume-normalized RMS norms,

    ||u|| <= kappa a^2 ||rho|| / 12,
    kappa a^2 ||rho_perp|| / 180 <= ||u_perp||
        <= kappa a^2 ||rho_perp|| / 20.

The lower bound uses this finite degree-6 source. The upper inhomogeneous
bound holds for any compatible mean-zero source with l>=2, after removing
the l=1 kernel. Do not apply either bound to an incompatible dipole source.
If supporting matter is allowed, rho becomes rho_phi + delta_rho_support:
give a parameterized upper bound and the reverse-triangle lower bound.
Unknown support with no norm bound supplies no universal scalar-only bound
on the total response; cancellation and reinforcement are both possible
as constraint data. Such choices are not asserted to obey a matter evolution.

### Leading standing wave and the momentum constraint

To measure at one consistent order, use the free leading scalar

    phi(t,x) = A(t) Y(x), A(t) = s cos(omega t),
    omega = 4/a, s = 0.2,
    Y = normalized Re[(x0 + i m.x_vec)^3], m=(1,2,3)/sqrt(14).

Let h_l be the harmonic coefficients of Y^2. The improved stress predicts

    rho = (A_dot^2 + omega^2 A^2) Y^2/2 + A^2 Delta(Y^2)/12,
    rho_l = [s^2 omega^2/2 - lambda_l A^2/12] h_l,
    j = grad J, J = -A A_dot Y^2/6.

Derive these from the inherited improved stress, and test against its full
pointwise implementation. A generic q,p in the multiplet need not have a
gradient momentum density; do not extend this special reduction to it.

For K^L = L grad w, where (LX)_ij = nabla_i X_j + nabla_j X_i
- (2/3) gbar_ij div X, the momentum equation gives

    div L grad w = (4/3) grad(Delta + 3/a^2)w = kappa grad J,
    w_l = 3 kappa J_l / [4(3/a^2-lambda_l)] (l=2,4,6),
    ||K^L||^2 = (8/3) sum_l lambda_l(lambda_l-3/a^2) |w_l|^2,
    kappa a ||j||/sqrt(30) <= ||K^L||
        <= kappa a sqrt(3/10) ||j||.

The additive constant in J or w has no effect. The standing wave has no
transverse-vector forcing; this does not remove freely specified vector
data or solve generic momentum constraints. Check compatibility with all
six Killing vectors and four gradient dipole conformal Killing vectors.
A control with p proportional to D_i q has nonzero Killing charge
-p^T D_i q and must not be silently projected into a solved constraint.

### What a comparison can establish

The leading scalar-induced tensor with zero initial tensor data is exact:

    beta_ind(t) = tensor(F(q0))/C * B(t), C=Vol/kappa,
    B(t) = (1-cos(Omega t))/(2 Omega^2)
        + (cos(2 omega t)-cos(Omega t))/(2(Omega^2-4 omega^2)),
    Omega^2 = 8/a^2.

Its all-time norm is bounded by ||F(q0)||/C times
`1/Omega^2 + 1/abs(Omega^2-4 omega^2)`. Compare RMS spatial metric
norms `4 sqrt(3)||u||_RMS`, also with its mean removed, against
`2||beta_ind||_F`, and separately the previously published primary TT
history. Report upper bounds and actual values without equating them.

Because u_l is affine in cos^2(omega t), its exact maximum RMS over a
full scalar cycle occurs at one of the two endpoints. Compute its minimum
from the quadratic norm on [0,1]; K^L is proportional to sin(2 omega t).
Record continuous-time envelopes, rather than calling sampled maxima bounds.
Use s=0.2,0.1,0.05 to check the predicted s^2 scaling of both constraint
responses and beta_ind. This does not prove equality of their effects on phi:
both induced metric sectors first feed back at order s^3, but a coefficient
or cancellation in the full scalar equation requires an evolution completion.
The independent initial TT field in the original primary run need not scale
as s^2; distinguish that run from the scalar-induced comparison.

Constraints determine neither a lapse nor the supporting pressure/stress
response. Consequently a bound on u or K^L alone is **not** a bound on the
complete omitted scalar acceleration/history. State whether an additional
closure is found in the cited machinery; do not manufacture one. Provide
the conditional metric bounds even if the evolution question remains open.

## Frozen checks and additional review follow-ups

Primary a=kappa=1, t in [0,4] at 401 times. Project Y^2 on full l=0,2,4,6
multiplets using the existing polynomial harmonic machinery, with S3 rules
(8,16) and (12,24). Check an independent rotated mode and random degree-3
q,p energy data, without treating generic momentum as the standing wave.
Keep mode coefficients in a compact archive, not all quadrature fields.

Required gates (absolute residual unless called scaled):

1. Reconstruct Y^2 and the random-data rho to relative L2 error below 1e-9;
   dipole overlaps below 1e-10; two quadrature results agree to scaled 1e-9.
2. rho and j formulas agree with inherited stress to scaled 1e-10. Check
   nonzero momentum times, not only the p=0 initial slice.
3. Pointwise Hamiltonian residual and momentum divergence residual below
   scaled 1e-9, using differentiated reconstructed fields. Derive the latter
   through invariant-frame derivatives and connection terms independently
   of the spectral denominator. Check norm identities to scaled 1e-9.
4. All spectral bounds hold to tolerance 1e-10; verify sharp upper constants
   on individual l=0 and l=2 sources. Dipole and charged momentum controls
   must be rejected, not silently removed. Amplitude-halving norm ratios
   equal four to 1e-8 when nonzero.
5. Independent forced-oscillator integration agrees with beta_ind to 1e-9.
   Continuous envelopes contain the sampled histories to 1e-10.
6. Replace the two literal TT constraint zeros in the original probe with
   computed certificates: existing ADM first curvature variation for the
   Hamiltonian, explicit connection contraction with a general symmetric
   STF rate tensor for momentum. Require exact symbolic zeros.
7. Explain n=1's quaternion anticommutator and n>1's nonzero coherent
   coupling; the constant n=0 mode is also dark. A possible mouth-position
   interpretation remains conditional on a field map.
8. Measure the full coupled history's source norm and distance to the
   uniaxial cone, including t=0 and t=2, without making nonzero leakage a
   success criterion or using that history to improve the leading order.
9. A failed gate produces UNRESOLVED numerical verdicts, named failures,
   valid JSON, overwritten stale reports, and nonzero exit status.

Separate verdicts: conditional initial-data constraints, transverse forcing
for the leading standing wave, size comparison, evolution/backreaction bound,
triangle map and readout. No full Einstein solution is an allowed success
verdict. Any correction to a frozen prediction is recorded explicitly.

Deliver a reusable constraint module, independent probe/tests, compact
archive, and derivation; update #289. Standard ADM/York conventions can be
cross-checked against Gourgoulhon, *Construction of initial data for 3+1
numerical relativity*, equations (1)-(2), (15), (18), (27), (34):
https://arxiv.org/abs/0704.0149. The bounds and the specialization above
are to be derived here, not attributed to that reference.
