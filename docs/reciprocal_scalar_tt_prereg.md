# Pre-registration: reciprocal scalar–TT dynamics before a source readout

Baseline: `22f77a373a67fcda91078ad0284be6a37c7ca20b` (main, including #287).
This branch is independent of the revisions to #288. This document is to be
committed and pushed before implementation or numerical results.

## Question and scope

Can the existing conformally coupled scalar and the five-component homogeneous
TT mode be evolved reciprocally from one action, and what does that establish
about a physical realization of the triangle variable?

The starting points are `waves/two_wave.py`, `waves/backreaction.py`, and
`bulk/tt_triangle_rotor.py`. The first two currently solve the scalar on a
fixed 3+1 Einstein static universe and use its stress as a prescribed tensor
drive. The third derives the tensor action and finds that a rotating uniaxial
restriction fails the omitted tensor equations. Keep all five tensor
components here. This is the 3+1 ESU channel, not the 4+1 throat solver.

The first closed calculation uses smooth scalar harmonics on all of S3,
without point emitters, excised mouths, a throat, or an external drive.
Consequently it is not already a localized mouth history. The smooth compact
domain also makes integration by parts legitimate; that step cannot silently
be transferred to punctured domains with unresolved boundary terms.

No rotor constraint, manufactured restoring stress, extra field, fitted
holonomy potential, counting function, ensemble, or future-setting boundary
condition is to be added. An initial value history is not a solution of the
two-boundary measurement problem merely because both of its endpoints exist.

## Prior analytic predictions, before measurement

Write the spatial metric in the inherited invariant frame as
`g_ij=a^2(exp(2 beta))_ij`, with `beta` symmetric and trace-free, and let
`C=Vol(S3)/kappa=2 pi^2 a^3/kappa`. The inherited tensor action is

    L_TT = C/2 [tr(beta_dot^2) - (8/a^2) tr(beta^2)].

Use the existing massless scalar action with conformal coupling `xi=1/6`:

    S_phi = -1/2 int sqrt(-g) [g^{mu nu} d_mu phi d_nu phi + xi R phi^2].

At first order in homogeneous STF beta, the volume variation and scalar
curvature variation vanish. The predicted interaction is therefore

    L_int = int beta_ij (grad_i phi)(grad_j phi) dV.

The improved stress terms must be checked: their integrated homogeneous TT
projection is predicted to vanish by transversality, compactness, and the
isotropy of the background Einstein tensor. Verify this against the existing
`backreaction.stress_series` implementation, not just a minimal-stress formula.

Expand the scalar in a complete real degree-n harmonic multiplet, normalized
by `int Y_alpha Y_beta dV=delta_alpha_beta`. Denote its canonical coefficients
by `q`, and the real matrices of the unit-sphere invariant derivatives by
`D_i`. They are antisymmetric, preserve harmonic degree, and obey
`-sum D_i^2=n(n+2) I`. The scalar frequency is `(n+1)/a`.

For an orthonormal STF basis `E_A`, define

    F_A = -sum_ij E_Aij (D_i D_j + D_j D_i)/(2 a^2).
    beta = sum_A b_A E_A.

The candidate reciprocal action and Hamiltonian are

    L = C/2 (b_dot^2 - omega_T^2 b^2)
        + 1/2 (q_dot^2 - omega_n^2 q^2) + sum_A b_A q^T F_A q,
    H = P^2/(2C) + C omega_T^2 b^2/2
        + p^2/2 + omega_n^2 q^2/2 - sum_A b_A q^T F_A q.

They predict BOTH equations, including the relative factor two:

    C (b_ddot_A + omega_T^2 b_A) = q^T F_A q,
    q_ddot + omega_n^2 q - 2 sum_A b_A F_A q = 0.

These formulas are design predictions to derive and verify, not discoveries
to claim after running. The Hamiltonian is a small-field truncation and need
not be globally bounded below. Do not infer global stability or equilibrium.

An anticipated selection rule is that the degree-1 scalar multiplet has
`F_A=0`: its symmetric generator product is isotropic. Degree 3 is the first
candidate in the odd-parity scalar sector that can source this TT mode.
The complete degree-3 multiplet has 16 real modes. Antipodal odd parity is a
stated sector choice, not a new proof of the repository's topology.

For the normalized real harmonic `Re[(x0+i m.x_vec)^n]`, with unit m, the
predicted integrated gradient tensor is
`[n I+n(n-1) m m^T]/a^2`. Its STF coupling is
`n(n-1) m^T beta m/a^2` per squared modal amplitude. Test this prediction
for n=1 and n=3. It would identify a field-derived quadrupolar interaction,
not a source-local pointer or the triangle-holonomy energy.

## Approximation and constraints: a separate gate

Retain the quadratic tensor and scalar actions and the interaction of order
`beta phi^2`. Omit `beta^3`, `beta^2 phi^2`, and the metric sectors outside
the homogeneous TT channel. The tensor force is accurate to the retained
order, and the reciprocal scalar correction is first order in beta.

Left-invariant derivatives preserve each complete scalar harmonic multiplet.
Thus a finite direct sum of complete multiplets is expected to be invariant
under this projected scalar equation. This is an algebraic statement about
the specified homogeneous metric, not a consistent truncation theorem for
Einstein's equations.

In particular, linear homogeneous TT perturbations have `delta G_00=0` and
`delta G_0i=0`, whereas a generic scalar history has nonzero energy density
and momentum density at order `phi^2`. The unknown ESU supporting matter,
lapse/shift, scalar/vector metric responses, higher tensor harmonics, and
nonlinear gravitational terms cannot simply be discarded while calling the
result a complete Einstein–matter solution. Measure these omitted sources,
including inhomogeneous energy after subtraction of its spatial mean.
Their absence from the retained equations must be reported explicitly.

The predicted modal outcome is: a reciprocal variational TT–scalar model
exists; a full gravitational constraint completion and a triangle-history
map are not supplied. A nonzero omitted constraint source is a scoped
obstruction to treating this TT-only calculation as a full solution, not a
no-go for the model with the missing metric/matter response included.

## Frozen implementation and numerical checks

Seed `2026090711`. Dimensionless primary run: `a=kappa=1`.
Use all five STF components and the full n=3 scalar multiplet. Construct
real harmonic polynomials and their invariant derivatives independently of
the final interaction matrices, so the representation and stress can be
checked against pointwise fields.

Primary initial tensor data: `beta=0.01(nn^T-I/3)`, `n=e_z`,
`beta_dot=0.01(v n^T+n v^T)`, `v=0.4 e_x`. Scalar data:
`q` is 0.2 times the normalized real n=3 harmonic along
`m=(1,2,3)/sqrt(14)`; `p=0`. Integrate `t in [0,4]`, recording 401 times,
with DOP853 at `(rtol,atol)=(1e-10,1e-12)` and `(1e-12,1e-14)`.
Compare with a one-way control in which the scalar follows its free history
and drives the full tensor; label the missing reciprocal term in that control.

Required checks:

1. Harmonic dimensions `(n+1)^2`, orthonormality, antisymmetric derivative
   matrices, Casimir identity, and preservation of the harmonic subspace:
   scaled residual below `1e-10` for n=1 and n=3. Check the n=5 multiplet's
   algebra as a higher-degree control; do not claim a tensor-tower refinement.
2. Integrated improved TT stress agrees with the action source within `1e-9`
   using two deterministic S3 quadrature orders. The n=1 null control and the
   n=3 coherent-harmonic quadrupole agree within `1e-10` in scaled norm.
3. Compare the first-order matter action with the scalar action on a static
   anisotropic metric at beta scales `0.02,0.01,0.005`. The omitted remainder
   must fall quadratically (successive ratios between `3.5` and `4.5`).
   Track the conformal curvature term, not just the gradient energy.
4. Independent finite differences of H give both equations to scaled `1e-7`;
   the mixed force derivatives satisfy reciprocity within `1e-10`. A control
   dropping the scalar reaction must fail the reciprocity identity for a
   nonzero coupling. Energy conservation is required for the reciprocal model.
5. The two ODE tolerances agree to scaled `1e-8`; relative energy drift is
   below `1e-8`. Require `max ||beta||_F < 0.05` on the primary interval.
   Record, without a preselected nonzero lower gate, the scalar-history
   difference from the one-way control, tensor shape, and Q_m(t).
6. Free scalar (n=1 null coupling) and scalar-zero free TT controls agree
   with their independent harmonic solutions within `1e-9`.
7. Common SO(3) rotations of the tensor, scalar harmonic, and apparatus axis
   preserve the action and transform the source covariantly within `1e-9`.
   The inherited identification n ~ -n makes Q_m even; it does not establish
   equality of early-record laws under future interventions.
8. Evaluate the scalar's energy and momentum sources with the inherited
   improved stress, on two quadrature grids and at initial and later primary
   times. Record the projection scope and nonzero omitted constraint sources.
   Numerical zeros are not a universal constraint-completion theorem.
9. A deliberately failed required check must render and overwrite an old
   report with UNRESOLVED verdicts, valid JSON, named failures, and exit code 1.

All residual definitions and normalizations belong in the code/report. A
failure is not permission to retune these gates silently. Preserve this freeze
and record any correction or changed test statement separately. Algebraic
certificates take precedence over apparent numerical zeros or convergence.

## What the field histories do and do not identify

At one uniaxial instant `beta=A(nn^T-I/3)`, Q_m equals
`A[(m.n)^2-1/3]`. Check that identity and the full evolving shape. A nearest
eigenline of a biaxial tensor is not a recovered autonomous triangle rotor.

Report separate verdicts for reciprocal action/evolution, scalar modal
closure, full Einstein–matter constraints, triangle-history map, physical
source localization/readout, and probability selection. No full source
readability test, operational retrocausal signal, non-readability theorem,
canonical ensemble, Phi, Born law, or tensor product is inferred here.

Successful closure of this round means a verified reciprocal model at its
stated order plus an honest account of the constraint and history-map gaps,
or a certified obstruction to that model. Numerical failure remains unresolved.

## Deliverables

One reusable `waves/reciprocal_scalar_tt.py` module, an executable probe,
independent tests, an archived report, and a derivation with assumption and
constraint ledgers. Link from the existing TT write-up. Keep existing solvers,
prior freezes, and #288's files unchanged. Open a draft PR for review.
