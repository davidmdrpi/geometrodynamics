# Reciprocal scalar–TT dynamics: a field interaction, with a constraint boundary

The existing conformal scalar and all five homogeneous TT components can be
evolved from one action. A complete degree-3 scalar multiplet supplies a
quadrupolar tensor interaction without a fitted coupling. Its reciprocal
reaction changes the scalar history and conserves the projected Hamiltonian.
**This is not yet a complete Einstein–matter history:** the scalar also has
nonuniform energy and momentum sources for metric equations omitted by the
TT-only projection. An explicit field configuration certifies that gap.

The [constraint-response follow-up](scalar_tt_constraints.md), frozen after
the review of `67940cb`, now bounds the linear scalar and longitudinal
responses on independent CMC slices with zero support perturbations. In that
problem the inhomogeneous scalar metric response exceeds an all-time upper
bound on the scalar-induced TT response. This does not yet bound the full
scalar backreaction or complete the Einstein–matter evolution.

The [public freeze](reciprocal_scalar_tt_prereg.md) was pushed at
[`d8dc90d`](https://github.com/davidmdrpi/geometrodynamics/commit/d8dc90d6d66e14824c337f96ee93a512dc9ed84f)
before implementation and numerical results. The baseline is main `22f77a3`,
including #287; this branch does not depend on #288. The interaction,
degree-selection rule, and constraint issue were analytic predictions in
that freeze. This implementation derives and verifies them.

## 1. Fields and approximation

Use the massless conformal scalar (`xi=1/6`) in
[`waves/two_wave.py`](../geometrodynamics/waves/two_wave.py) and the homogeneous
3+1 ESU tensor in
[`waves/backreaction.py`](../geometrodynamics/waves/backreaction.py).
The latter previously used the scalar stress as a prescribed drive. The
[TT reduction](tt_triangle_rotor.md) derived its action and showed that
discarding the two biaxial shape equations does not produce a free rotor.

Here the tensor is a general STF matrix in the invariant frame on S3:

\[
 g_{ij}=a^2(e^{2\beta})_{ij},\qquad
 \beta=\sum_{A=1}^5 b_A E_A,\qquad
 \operatorname{tr}(E_AE_B)=\delta_{AB},\qquad
 C=\frac{2\pi^2a^3}{\kappa}.
\]

The domain is the whole smooth S3. There are no punctures, point emitters,
mouth boundary terms, or throat in this calculation. In particular a global
harmonic is not a source-localized mouth field.

Retain `beta^2`, `phi^2`, and the mixed term `beta phi^2` in the action;
omit `beta^3`, `beta^2 phi^2`, and metric sectors outside homogeneous TT.
The resulting finite Hamiltonian is evolved exactly up to ODE error, but
that does not improve its physical truncation order. It is not asserted to
be bounded below globally, to equilibrate, or to define a probability law.

## 2. One action fixes both forces

Begin with the inherited scalar action

\[
 S_\phi=-\frac12\int\sqrt{-g}
 \left(g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi+\frac16R\phi^2\right)d^4x.
\]

For homogeneous STF beta on the round ESU, the first variations of volume
and scalar curvature vanish. The inverse spatial metric varies by
`-2 beta/a^2`. Thus the first-order interaction is

\[
 L_{\rm int}=\int_{S^3}\beta_{ij}\nabla_i\phi\nabla_j\phi\,dV.
\]

This also follows from the **improved** stress used in the repository.
The isotropic stress pieces vanish on STF contraction. For the Hessian
improvement term, transversality gives
`int beta_ij nabla_i nabla_j(phi^2) dV=0` by integration by parts on compact
S3; the background Einstein tensor is isotropic. This identity relies on
the smooth domain and does not justify dropping unresolved mouth terms.

Expand phi in an orthonormal real scalar harmonic multiplet of polynomial
degree n, `int Y_alpha Y_beta dV=delta_alpha_beta`. Its real dimension is
`(n+1)^2`. The unit-sphere invariant derivatives have antisymmetric matrices
`D_i` with `-sum_i D_i^2=n(n+2)I`. Scalar n here denotes polynomial degree;
the degree-3 scalar frequency is `4/a`, distinct from the homogeneous tensor
frequency `sqrt(8)/a`.

With scalar amplitudes q and

\[
 F_A=-\frac{1}{2a^2}\sum_{ij}(E_A)_{ij}(D_iD_j+D_jD_i),
 \quad\omega_T^2=8/a^2,\quad\omega_n^2=(n+1)^2/a^2,
\]

the projected action is

\[
 L=\frac C2(\dot b^2-\omega_T^2b^2)
   +\frac12(\dot q^2-\omega_n^2q^2)+\sum_A b_Aq^TF_Aq.
\]

Independent variations give

\[
 C(\ddot b_A+\omega_T^2b_A)=q^TF_Aq,\qquad
 \ddot q+\omega_n^2q-2\sum_A b_AF_Aq=0.
\]

The relative factor two and its sign are fixed by the action. There is no
new interaction strength: radius, gravitational coupling, mode normalization,
and the scalar initial data determine the response. With `P=C b_dot` and
`p=q_dot`, the Hamiltonian includes the interaction energy:

\[
 H=\frac{P^2}{2C}+\frac{C\omega_T^2b^2}{2}
  +\frac{p^2}{2}+\frac{\omega_n^2q^2}{2}-\sum_A b_Aq^TF_Aq.
\]

It is the energy of this projected system, not the full gravitational
Hamiltonian constraint. Omitting the scalar reaction while retaining its
tensor source fails the mixed-derivative reciprocity test.

## 3. Which scalar modes can couple?

For the complete n=1 multiplet, the symmetric derivative product is isotropic:
every `F_A` vanishes. This is an operator identity, not a zero found for a
particular scalar trajectory. The first active **odd-degree** candidate is
n=3, with 16 real scalar modes. Odd parity is an explicit choice of sector;
it is not a new derivation of the antipodal identification.

A concrete normalized mode is

\[
 Y_{n,m}(x)\propto\operatorname{Re}(x_0+i\,m\cdot\mathbf x)^n,
 \qquad |m|=1.
\]

Its gradient energy along m is `n^2/a^2`, since the invariant derivative
along m rotates its phase at rate n. Axial symmetry makes the other two
entries equal; the Casimir fixes their sum, giving each `n/a^2`. Therefore

\[
 \int\nabla_iY_{n,m}\nabla_jY_{n,m}\,dV
 =\frac{n\delta_{ij}+n(n-1)m_im_j}{a^2}.
\]

For a scalar amplitude s in this mode, at that instant,

\[
 L_{\rm int}=\frac{n(n-1)s^2}{a^2}\,Q_m,
 \qquad Q_m=m^T\beta m.
\]

At n=3 the coefficient is `6s^2/a^2`. This realizes the quadrupolar source
term with an existing scalar field. The full 16-mode multiplet evolves; the
coherent mode is initial data and is not kept fixed by an external actuator.
This field interaction does not identify m with an operational analyzer, or
Q_m with a completed pointer record.

The polynomial construction tests n=1,3,5 multiplets of dimensions 4,16,36.
Invariant differentiation commutes with the round Casimir and preserves
polynomial degree and harmonicity. Each complete multiplet is consequently
invariant under the **projected scalar equation**. The n=5 algebra control
does not constitute a refinement of the gravitational tensor tower.

## 4. Coupled histories and independent checks

The frozen primary data use `a=kappa=1`,
`beta(0)=0.01(e_z e_z^T-I/3)`, initial tangent speed `0.4 e_x`, and scalar
amplitude `0.2 Y_{3,m}` with `m=(1,2,3)/sqrt(14)` and zero scalar momentum.
Both runs record 401 times on `[0,4]`. The one-way control evolves the free
scalar and lets it drive the full tensor; it deliberately drops reciprocal
scalar back-action.

| Quantity | Measured value |
|---|---:|
| Fine relative Hamiltonian drift | `6.74e-12` |
| Maximum state difference between the two ODE tolerances | `3.70e-11` |
| Maximum tensor Frobenius norm | `0.008188` |
| Maximum scalar coefficient distance from one-way history | `0.001081` |
| Maximum tensor coefficient distance from one-way history | `1.24e-05` |
| One-way history's relative defect in the reciprocal Hamiltonian | `0.002226` |
| Maximum distance to the full uniaxial cone | `0.001157` |

The trajectory is not confined to the uniaxial cone. All five equations
remain active; no nearest-axis projection is used as an evolution rule.
The scalar-history distance is invariant under a common orthogonal change
of modal basis. Individual archived scalar coefficients depend on the
orthonormal basis returned by the polynomial construction.

The action source agrees with a pointwise polynomial evaluation passed
through the repository's improved stress and frame projection to `1.5e-15`,
on deterministic grids of 2048 and 6912 points. The modal scalar equation
agrees with the independently differentiated pointwise field to `8.1e-16`.
The exact static anisotropic scalar action, **including conformal curvature**,
differs from the first-order action by a quadratic remainder; successive
halving ratios are `3.986` and `3.993`.

Hamilton equations are checked by differences of H, and reciprocity by
differences of both force implementations, including the one-way control.
Those checks use a coordinate step `1e-4`: each expression has degree at most
two in an individual coordinate, so centered differences have no Taylor
truncation error in real arithmetic. The step avoids unnecessary floating
point cancellation without changing any frozen acceptance gate.
Free TT, n=1 decoupling, common SO(3) covariance, and stale-report replacement
have separate tests. The probe has 13 numerical gates; CLI failure archival
is an additional end-to-end test.

## 5. Why this is not a complete gravitational history

The inherited linear homogeneous TT mode has
`delta G_00=delta G_0i=0`. Generic scalar matter does not have zero sources
for those equations. Removing its mean energy is insufficient.

There is a simple certificate for the frozen initial scalar. For n=3,
`Y=sqrt(8/Vol) Re(x0+i m.x_vec)^3`. At `p=0`, the improved energy density is

\[
 T_{00}=\frac12|\nabla\phi|^2+
 \frac16\left(-\Delta\phi^2+\frac3{a^2}\phi^2\right).
\]

At unit radius and `x0=1`, phi has zero spatial gradient and
`Delta phi=-15 phi`; hence `T00=44 s^2/Vol`. At `x0=m.x_vec=0`, phi and
its first derivatives vanish and `T00=0`. For `s=0.2` the contrast is
`0.0891626416053`; the direct stress evaluation matches the certificate to
`4.2e-17`. A homogeneous adjustment to the supporting fluid cannot remove
this nonuniform part.

| Time | Mean scalar energy density | Inhomogeneous energy RMS | Momentum density RMS |
|---:|---:|---:|---:|
| 0 | `0.0162114` | `0.0185345` | `0` |
| 1 | `0.0161951` | `0.0221626` | `0.00553063` |
| 4 | `0.0161900` | `0.0189900` | `0.00306291` |

These are round-background improved stress sources at the retained order
`phi^2`, not a computation of the full metric-corrected stress. Grid
disagreement in the listed quantities is below `4.6e-17`.

This is a scoped obstruction to calling the TT-only history a solution of
all Einstein equations. Additional metric response and/or a specified
inhomogeneous supporting-matter response is needed. No impossibility of
such a completion is established. Higher tensor harmonics and the other
spatial equations are also outside this projection.

For a scalar-induced tensor with `beta=O(phi^2)`, reciprocal scalar feedback
starts at order `phi^3`. Scalar/vector metric perturbations sourced at order
`phi^2` can affect the scalar at that same order. Thus the feedback measured
here cannot be promoted to a controlled full-GR correction while those
responses are missing. The quadratic and cubic terms omitted from the
action also limit the physical accuracy of the tensor feedback.

## 6. Consequence for the triangle and the causality gate

At the initial uniaxial instant,
`Q_m=A[(m.n)^2-1/3]`, an even function of n. This is a real field contraction
with a derived scalar interaction. It supplies a stronger starting point
than declaring the rotor locally accessible, but supplies no preparation
law for closure-conditioned source ensembles.

After the tensor becomes biaxial, its continued eigenlines do not establish
an autonomous triangle rotor. No map of these IVP histories to analyzer
itineraries or to the triangle two-boundary problem has been derived. Nor
has the global scalar harmonic been turned into a localized apparatus.

| Item | Status |
|---|---|
| Reciprocal scalar–TT interaction and projected histories | Derived and checked |
| Complete scalar multiplet closure under the projected equation | Algebraically invariant |
| Conditional linear CMC constraint response | Bounded in the [follow-up](scalar_tt_constraints.md), with zero support perturbations |
| Full Einstein–matter evolution and constraint propagation | Missing |
| Physical triangle-history map | Not derived |
| Local apparatus and early-record law | Not derived |
| Phi, Born law, canonical preparation | Not derived |

The next completion needs a specified supporting-matter response and the
metric sectors demanded by the constraints, followed by a source/apparatus
history map. This round establishes neither an operational retrocausal
channel nor a dynamical non-readability theorem.

## Reproduction

```bash
python -m experiments.closure_ledger.reciprocal_scalar_tt_probe
python -m pytest -q tests/test_reciprocal_scalar_tt.py tests/test_tt_triangle_rotor.py
```

The [archived report](../experiments/closure_ledger/runs/20260907_reciprocal_scalar_tt/probe.md)
links its full JSON companion through the same directory. It records the
frozen data, checks, histories, and separate physical verdicts. The original
unchanged freeze is the parent of the original implementation; the separately
frozen constraint follow-up preserves its archive and acceptance thresholds.
This PR adds modules and evidence, and leaves prior field solvers and #288's
files unchanged.
