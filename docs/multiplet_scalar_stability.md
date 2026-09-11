# The four-component ESU has an admissible growing scalar mode

**The background is linearly unstable in its homogeneous scalar sector.**
The physical scale perturbation obeys `r''=2r` in background conformal time
`eta=t/a`. Its proper-time growth rate is `sqrt(2)/a`. A nearby exact,
constraint-satisfying Einstein–four-scalar family has this growing tangent;
the instability is not an unconstrained energy perturbation or a clock mode.

The scalar-type degree-1 block on the unrestricted S3 cover is different:
after constraints and gauge fixing it has a semisimple period map `-I`.
Its nonzero field perturbations are even, and therefore are excluded if the
inherited componentwise odd field condition is retained. Exclusion is not
stability. No nonlinear completion of that cover block is established here.

Freeze: [`8eb33d5`](multiplet_scalar_stability_prereg.md), published before
implementation and measurements, on #295 at `9a0bfbd`. The homogeneous
predictions required no correction. The dipole reduction and its neutral
answer were derived **after** the freeze; P2 did not predict their outcome.

## 1. The matter model, rather than a fluid substitution

Use exactly #294's action for four independent real conformal scalars:

\[
 I=\int\sqrt{-g}\left[\frac{R-2\Lambda}{2\kappa}
 -\frac12\sum_{I=0}^3\left((\nabla\phi_I)^2+\frac{R}{6}\phi_I^2\right)\right].
\]

The background is `phi_I=P x_I`, with
`P=sqrt(3/(4 kappa)) cos(2t/a+delta)` and `Lambda=3/(2a^2)`.
Each component is odd on the unit S3 embedding. No extra fluid, constitutive
law, signal scalar, damping or support adjustment is introduced.

Scalar degree here labels representations under simultaneous spatial and
internal SO(4), as in #295. It does not label the harmonic degree of each
field perturbation. These low scalar blocks do not exhaust arbitrary
four-field perturbations, and no vector block is calculated in this round.

## 2. Exact homogeneous Einstein–matter reduction

Allow a variable scale and use conformal time:

\[
 ds^2=A(\eta)^2(-d\eta^2+\gamma),\qquad \phi_I=\frac{q(\eta)}{A(\eta)}x_I.
\]

This is the general diagonal-SO(4)-invariant scalar field perturbation:
at fixed x the stabilizer SO(3) leaves only the radial internal vector x.
The independent spatial identities are `sum x_I^2=1`,
`sum dx_I tensor dx_I=gamma`, `Delta x_I=-3x_I`.
Evaluating the full off-shell improved stress, including its `G_mu_nu phi^2/6`
term, gives

\[
 \rho=\frac{q'^2+4q^2}{2A^4},\qquad
 p=\frac{q'^2-2qq''-4q^2}{6A^4},\qquad
 (\Box-R/6)\phi_I=-\frac{q''+4q}{A^3}x_I.
\]

Thus `q''+4q=0` implies `E=(q'^2+4q^2)/2` constant and `p=rho/3`.
This radiation law has been derived for the **whole evolving homogeneous
family**, not transferred from the ESU's background equation of state.
Einstein's Hamiltonian constraint and trace evolution are

\[
 C=A'^2+A^2-\frac{\Lambda A^4}{3}-\frac{\kappa E}{3}=0,
 \qquad A''+A-\frac{2\Lambda}{3}A^3=0.
\]

The momentum and spatial TF equations vanish pointwise. Direct differentiation
gives `C'=0` on the two evolution equations. The full coordinate-curvature
route verifies these expressions off shell as well as on exact solutions.

At the ESU, `A=a`, `E0=3a^2/(2 kappa)` and
`q0=a sqrt(3/(4 kappa)) cos(2eta+delta)`. Put
`A=a(1+epsilon r)`, `q=q0+epsilon z`. Then

\[
 r''=2r,\qquad z''+4z=0,\qquad
 \delta E=q_0'z'+4q_0z=0.
\]

The oscillator constraint removes its amplitude perturbation. The surviving
solution `z=c q0'` changes the phase/time origin. In the conformal lapse gauge,
a residual constant time translation generates precisely that solution.
With a fixed external clock it can instead be retained as a neutral phase
coordinate; it contributes no growing exponent.

In contrast, `delta A` is invariant under a linear time shift because
`A0'=0`. The proper spatial volume is `2pi^2 A^3`, with fractional variation
`3 epsilon r`. This observable detects the growing mode without reading a
lapse, field phase or coordinate-time artifact.

## 3. The degenerate constraint is checked beyond first order

The linear constraint fixes `delta E`, but neither r nor r' appears in it.
That degeneracy is not permission to regard every unconstrained evolution
mode as physical. Keep **the same** Lambda and E=E0. The exact constraint is

\[
 A'^2=\frac{(A^2-a^2)^2}{2a^2}.
\]

Both branches

\[
 A'=\sigma\frac{A^2-a^2}{\sqrt2 a},\quad
 w=\frac{A-a}{A+a}=w(0)e^{\sigma\sqrt2\eta},\quad
 A=a\frac{1+w}{1-w},\qquad \sigma=\pm1,
\]

solve the trace equation and all remaining field equations with unchanged q.
Taking `A(0)=a(1+epsilon)` gives `w(0)=epsilon/(2+epsilon)` and
`delta A/a=exp(sigma sqrt(2) eta)`. This establishes exact local continuations
of both eigenvectors, including the velocity needed by the quadratic
constraint. The field remains componentwise odd. Positivity of A and
`f=1-kappa q^2/(6A^2)` holds in a neighborhood of the background, where
`f>=7/8` at epsilon=0. The formula is used only before its scale chart fails;
it is not a global regularity claim.

The proper clock is `dt=A d eta`. At first order and fixed `t/a=tau`,

\[
 \eta=\tau-\epsilon\int_0^\tau r(s)ds,
 \quad \delta\phi_I=\left[-P r-P'\int_0^\tau r(s)ds\right]x_I
\]

for the growing branch with unchanged conformal q and synchronized initial
clocks. Omitting the second term gives the wrong matter response. The probe
independently evolves the same exact branch in proper time, checks A, the
field and volume, and verifies the O(epsilon^2) central-variation convergence.
The growth rate is `sqrt(2)/a` in proper time in either calculation.

## 4. Degree 1 requires its own gauge and constraint reduction

Take `Y=d.x`, `|d|=1`. At a point, the internal radial and tangential parts of
d give the two scalar-type equivariant vector structures. Thus the complete
pair in this representation is

\[
 \sqrt\kappa\,\delta\phi_I=u x_IY+v\nabla^a x_I\nabla_aY
 =(u-v)x_IY+v d_I.
\]

It combines individual field degrees 0 and 2, unlike the degree-1 background.
No third equivariant structure exists: the stabilizer SO(3) decomposes d
into its radial singlet and tangent vector. Vector-type representations and
other internal modes are not being silently absorbed into this pair.

Write the dimensionless scalar metric perturbation as

\[
 a^{-2}ds^2=-(1+2\alpha Y)d\eta^2
       +2B\nabla_iY\,d\eta dx^i+(1+2hY)\gamma_{ij}dx^idx^j.
\]

Since `Hess Y=-gamma Y`, there is no scalar TF spatial harmonic. Under minus
the Lie derivative by `xi=T Y partial_eta + L grad Y`,

    h -> h+L,  B -> B+T-L',  alpha -> alpha-T',
    u -> u-T P',  v -> v-L P,

where now `P=sqrt(3/4) cos(2eta+delta)` is dimensionless.
Choose `L=-h` and `T=-B-h'` to set h=B=0. No residual scalar-type degree-1
gauge transformation preserves both conditions. In particular a nonzero
matter perturbation with zero metric in this gauge is not a discarded
coordinate mode. A pure-gauge control with all lapse, shift and spatial terms
included gives zero full Einstein/KG response; omitting its metric does not.

In this gauge, direct geometric variation yields

    delta G^0_0=delta G^0_i=0,
    delta G^i_j=-2 alpha Y delta^i_j/a^2,
    delta R=6 alpha Y/a^2.

The on-shell conformal stress is traceless, so the Einstein trace forces
alpha=0. This uses no nonexistent TF equation and divides by neither P nor P'.
All four Klein–Gordon equations then reduce to

\[
 u''=-7u+6v,\qquad v''=2u-3v.
\]

The two remaining Einstein constraints are

\[
 H=P'u'+5Pu-3Pv=0,\qquad
 J=Pv'-Pu'/3+2P'u/3-P'v=0.
\]

In units `1/(kappa a^2)`, H is the density coefficient and J is the
covariant momentum coefficient. The spatial pressure coefficient reduces
to H/3 on KG; spatial TF is identically zero. Differentiating gives
`H'=-3J`, `J'=H/3`, proving constraint preservation. The constraint matrix
has rank two at **every** phase. For P nonzero its J row has a v' coefficient
while H does not; for P=0, P' is nonzero and the two rows have disjoint
nonzero supports. No field-zero singularity or phase patch is needed.

The unconstrained two-oscillator frequencies are 1 and 3. The constraints
relate their amplitudes; at delta=0 a useful exact parameterization is

    u=A3(cos 3eta - 3 cos eta) + B3(sin 3eta + 3 sin eta),
    v=-A3(cos 3eta/3 + 3 cos eta) + B3(-sin 3eta/3 + 3 sin eta).

All these perturbations have zero summed linear stress, while the individual
fields change. Their coefficients and derivatives are antiperiodic after pi.
The **matrix** map on the two-dimensional constraint subspace is `-I`; this
excludes a Jordan block as well as exponential growth. Each of the four
independent directions d has the same block. This is a constrained linear
cover result, not a nonlinear integrability theorem for these dipoles.

## 5. Periods and antipodal admissibility

The primary period is the full background field period `T_eta=pi`. The
fields reverse sign after pi/2, although their stress and #295's tensor
coefficients repeat. Quotienting by the internal sign symmetry would require
an explicit endpoint identification for a matter period map. No such
identification is used for the dipole verdict here.

For the homogeneous metric pair `(r,r')` the constant generator is
`[[0,1],[2,0]]`. Its exact map over any T is

\[
 M(T)=\begin{pmatrix}
 \cosh(\sqrt2T)&\sinh(\sqrt2T)/\sqrt2\\
 \sqrt2\sinh(\sqrt2T)&\cosh(\sqrt2T)
 \end{pmatrix}.
\]

| Block over pi | Map classification | Multipliers |
|---|---|---|
| Homogeneous physical scale | Hyperbolic | 85.019695223207, 0.011761980531 |
| Homogeneous phase, if retained | Neutral time-origin direction | 1 |
| Constrained degree-1 cover pair | Neutral, semisimple, antiperiodic | -1, -1 |

The homogeneous trace is `85.031457203739`, determinant 1. The archive also
contains its half-period map, whose square gives the full map. Positivity of
the supporting kinetic matrix has not been used as a stability criterion.

Under x -> -x, both `x_I Y` and `d_I` are even. Therefore all nonzero dipoles
above violate the componentwise odd field condition inherited from #294.
If the metric is also required to be invariant under this spatial antipodal
map, its degree-1 lapse/shift/trace variations are excluded too. That metric
restriction is a named preparation/quotient condition, not a field equation.
On an unrestricted S3 cover the neutral block exists; in the stated
restricted class there is no admissible physical dipole map to label elliptic.
The homogeneous instability survives either choice.

## 6. Verification and what changes downstream

The probe verifies symbolic off-shell FRW stress/KG identities, the
Hamiltonian and evolution linearizations, exact constraint propagation,
and the dipole pressure/constraint reduction. It also uses coordinate metric
jets and the existing off-shell improved-stress routine for an independent
route that never calls the reduced evolution equation to calculate curvature.

There are 18 off-shell homogeneous and 18 off-shell dipole variation cases,
18 pure-gauge cases, 96 constrained dipole cases, 10 finite off-shell FRW
cases and 72 exact finite FRW cases, including both field and velocity zeros,
three radii, two kappa values and generic spatial points. Exact FRW normalized
Einstein residuals are at most `2.73e-15`. Both frozen DOP853 tolerances agree
with the analytic maps, with proper-time and phase controls. The constrained
dipole map differs from `-I` by `4.27e-15` in the principal run.

Ten frozen gates pass. The validator recomputes map comparisons and
determinants from the matrices, and does not accept true gate flags without
period evidence. Missing identities, missing/failed gates, nonfinite or
malformed results and a wrong map produce UNRESOLVED, including through the
CLI. Separate homogeneous evidence is retained if a dipole check fails.

This refutes linear stability of this four-component background. It does
not overturn #294's existence or equal-stress construction, or #295's
elliptic tensor block. It shows that restricting to tensor-only initial data
omits a growing admissible direction. In the linear homogeneous plane,
removing future exponential growth imposes the condition `c_plus=0` on
`r=c_plus exp(sqrt(2) eta)+c_minus exp(-sqrt(2) eta)`; no measure on the
full preparation/history space has been supplied, so a global measure-zero
claim does not follow.

The next rotor calculation would concern an unstable supporting background
unless a stated preparation or boundary mechanism excludes or controls this
mode. The instability alone does not construct such a mechanism, establish
inevitable collapse, or rule out two-boundary histories. Tensor-driven
quadratic scalar/vector sources, vector stability, a coupled rotor,
O(s^4) persistence, Phi selection and the causality gate remain uncomputed.

Reproduce with:

```bash
python -m experiments.closure_ledger.multiplet_scalar_stability_probe \
  --output experiments/closure_ledger/runs/20260911_multiplet_scalar_stability
python -m pytest -q tests/test_multiplet_scalar_stability.py
```
