# The existing TT mode does not close on a free triangle rotor

The field map and its restricted action can be derived, but the proposed
rotating family **fails the full field equations**. In the repository's
homogeneous ESU tensor mode, an axisymmetric tensor with a moving director
immediately develops biaxial shape components. A round angular kinetic term
appears after substituting the family into the action; it does not make that
family an autonomous rotor.

The obstruction is exact within the existing linear TT model. It concerns
the candidate `beta=A(nn^T-I/3)`, with no added shape constraint or stress.
A manufactured source can sustain its rotation, and is retained as a control.
This is not a no-go for every driven solution, field mode or BAM completion.
Neither the history weight `Phi` nor an operational source-readout channel
is derived by the failed reduction.

The [public preregistration](tt_triangle_rotor_prereg.md), commit
[`0eb684b`](https://github.com/davidmdrpi/geometrodynamics/commit/0eb684b86a57be9cfbd2d41d89377cc11f6b76cf),
preceded implementation and all new numerical calculations. The normal
obstruction was an analytic prediction in that freeze. This round tests and
derives it; it is not described as an unexpected numerical discovery.

## 1. Which existing field is being reduced?

[`waves/backreaction.py`](../geometrodynamics/waves/backreaction.py) derives
the homogeneous, symmetric trace-free tensor equation on a fixed ESU with
spatial `S3`:

\[
 \ddot\beta+\omega^2\beta=S,\qquad \omega^2=8/a^2.
\]

There are five independent tensor components. The existing Cartan/Einstein
derivation checks the round spatial curvature, background Einstein tensor,
trace and linear momentum constraint. Its scalar-wave machinery computes
the anisotropic stress sourcing this mode.

This is a **linear response in a 3+1 ESU sector**. It is not the nonlinear
4+1 spherically symmetric Einstein-scalar solver, nor a resolved mouth
boundary problem. The scalar sources in the response calculation are not
evolved reciprocally with the metric. No nonlinear GR constraint completion
or self-consistent matter feedback is asserted here.

Choose the proposed field family

\[
 \beta=A Q(n),\qquad Q(n)=nn^\mathsf T-I/3,\quad |n|=1.
\]

`A` is the signed anisotropy amplitude, unrelated to analyzer labels. For
`A != 0`, its eigenvalues are `2A/3,-A/3,-A/3`. Thus `n` and `-n` are exactly
the same field. The angular space is `RP2`; the family is a three-dimensional
cone inside the five-dimensional STF space. At `A=0`, the director is
undefined. The two absent directions split the repeated eigenvalues.

TT gauge invariance at linear order does not provide an absolute physical
orientation. Tensor components and a source/apparatus frame must be rotated
together. Relational contractions such as `m^T beta m` survive that common
change of frame. Identifying this director with the triangle variable `x`,
its source event and its analyzer itinerary would still require a map of
the actual field histories.

## 2. The normalized action follows from ADM

Use lapse one, fixed radius `a`, and homogeneous spatial metric
`g=a^2 exp(2 beta)` in the unit-quaternion coframe. The gravitational ADM
action is

\[
 S_g=\frac1{2\kappa}\int dt\,d^3z\sqrt g
       ({}^{(3)}R+K_{ij}K^{ij}-K^2-2\Lambda).
\]

Since `tr beta=0`, the spatial volume is exactly `V3=2 pi^2 a^3`. For
diagonal principal lengths `l_i`, an independent three-dimensional Koszul
calculation, using `[e_i,e_j]=2 epsilon_ijk l_k/(l_i l_j)e_k`, gives

\[
 {}^{(3)}R=
 \frac{2[-\sum_i l_i^4+2\sum_{i<j}l_i^2l_j^2]}{l_1^2l_2^2l_3^2}.
\]

Set `l_i=a exp(beta_i)`, `sum beta_i=0`. The expansion is

\[
 {}^{(3)}R=6/a^2-(8/a^2)\operatorname{tr}\beta^2+O(\beta^3),\qquad
 K_{ij}K^{ij}-K^2=\operatorname{tr}\dot\beta^2+O(\beta^4),
\]

where the order counts the perturbation and its time derivative together.
Rotational invariance extends the quadratic result to all STF components.
A separate matrix-exponential/Frechet-derivative calculation checks the
kinetic expression for noncommuting `beta,beta_dot`.

With `C=V3/kappa`, and the matter variation
`delta S_m=integral dt integral dV T_TT:delta beta`, the quadratic action is

\[
 L=\frac C2[\operatorname{tr}\dot\beta^2-
             \omega^2\operatorname{tr}\beta^2]
      +C\operatorname{tr}(\beta S),\qquad
 S=\frac{\kappa}{V_3}\int T_{\rm TT}\,dV.
\]

The isotropic supporting stress and cosmological term do not contribute an
anisotropic source in this approximation. This treats the matter stress as
the specified drive, as in the inherited linear response calculation.

The normalization matters: `backreaction.shear_sources` returns an **integral**
of the stress, while `shear_response` solves its normalized oscillator
equation for whatever `S` is passed. A physical use with explicit `kappa,a`
requires the factor `kappa/V3`. The new module makes that conversion explicit
and leaves the existing solver and its archived convention unchanged.
At `a=kappa=1`, `C=2 pi^2`, and the factor is `0.05066059182`.

## 3. The restricted action looks like a rotor

For `v=dot n`, `n.v=0`,

\[
 \operatorname{tr}\dot\beta^2=\frac23\dot A^2+2A^2|v|^2,
 \qquad \operatorname{tr}\beta^2=\frac23A^2.
\]

Substitution gives

\[
 L_{\rm restricted}=C\left[\frac{\dot A^2}{3}
                 +A^2|\dot n|^2-\frac{\omega^2 A^2}{3}
                 +A n^\mathsf T S n\right].
\]

In coordinates `(A,theta,phi)`, its kinetic metric is

\[
 G=\operatorname{diag}\left(\frac{2C}{3},\ 2CA^2,
                                  \ 2CA^2\sin^2\theta\right).
\]

The formal angular inertia is `I(A)=2CA^2`. It is round at fixed amplitude,
but neither a fixed amplitude nor an invariant angular family has yet been
established. The induced configuration volume scales as `A^2 dA dOmega`
on an `RP2` chart. It is a geometric pullback, not a selected ensemble of
the unconstrained field: this cone has zero five-dimensional configuration
volume, and restricting to it requires an additional conditioning or
confinement prescription. An isotropic full tensor ensemble can have uniform
axis marginals without being an autonomous rotor ensemble.

## 4. The omitted field equations obstruct the reduction

Let `P=I-nn^T`. Varying the restricted action, including the amplitude,
gives

\[
 \ddot A=3A|v|^2-\omega^2 A+\frac32n^\mathsf T S n,
 \qquad
 \ddot n=-|v|^2n-2\frac{\dot A}{A}v+\frac{PSn}{A}.
\]

These solve the three tangent equations of the cone. They need not solve
the remaining two. For a symmetric tensor `M`, define the transverse
trace-free projector

\[
 \mathcal N_n(M)=PMP-\frac12\operatorname{tr}(PMP)P.
\]

Twice differentiating the **field embedding**, before discarding any
equation, gives

\[
 \mathcal N_n(\ddot\beta+\omega^2\beta-S)
 =2A\left(vv^\mathsf T-\frac{|v|^2}{2}P\right)-\mathcal N_n(S).
\]

The amplitude derivatives, radial acceleration and angular acceleration
cannot cancel this component. In an orthonormal transverse basis with its
first vector along nonzero `v`, the free residual has eigenvalues
`A|v|^2,-A|v|^2`, hence

\[
 \|\mathcal N_n(\ddot\beta+\omega^2\beta)\|_F
       =\sqrt2\,|A|\,|v|^2>0\quad(A\ne0,\ v\ne0).
\]

**No source-free rotating solution stays in this uniaxial family on a regular
interval.** The full field develops biaxial components. A constant director
with freely oscillating amplitude is the invariant control. Freezing `A`
would add a separate radial support requirement; even a stationary tensor
of nonzero fixed amplitude needs `S=omega^2 beta`.

This is not a higher-order gravitational correction that can be dropped
while retaining director motion. At fixed director speed the residual is
linear in `A`, the same field order as the retained TT equation. Taking
smaller amplitudes leaves the relative departure unchanged. Slow-motion
approximations are a different question: the normal force and the angular
kinetic dynamics both first enter at order `|v|^2`.

For a driven solution, the necessary normal condition is instead

\[
 \mathcal N_n(S)=2A(vv^\mathsf T-|v|^2P/2).
\]

We exhibit `A=0.01`, `n=(sin(0.4t),0,cos(0.4t))` and the manufactured source
`S=beta_ddot+omega^2 beta`. Independent full-tensor integration reproduces
that rotation. Its radial and normal source norms are `0.0614005` and
`0.00226274`, respectively. This is added forcing data, not stress derived
from the scalar source histories. Whether existing coupled matter supplies
the required force remains open.

## 5. Unconstrained evolution independently confirms the departure

Use the frozen data `A0=0.01`, `A_dot0=0`, `n0=e_z`, `v0=0.4 e_x`.
The tensor evolves exactly as
`beta(t)=beta0 cos(omega t)+beta_dot0 sin(omega t)/omega`.
No projection back to the cone is applied.

At each time we minimize the Frobenius distance over **all** directors and
signed amplitudes. For a fixed director the minimizing amplitude is
`3 n^T beta n/2`; globally the director is an eigenvector of largest absolute
eigenvalue. Equivalently the distance is the absolute difference of the
remaining two eigenvalues divided by `sqrt(2)`. Both evaluations are checked.

| time | distance to the full uniaxial cone | distance / time squared |
|---:|---:|---:|
| 0.04 | 1.81361262e-6 | 0.00113350789 |
| 0.02 | 4.52761062e-7 | 0.00113190265 |
| 0.01 | 1.13150365e-7 | 0.00113150365 |
| 0.005 | 2.82851010e-8 | 0.00113140404 |

The analytic coefficient is `|A0| |v0|^2/sqrt(2)=0.00113137085`.
Step-halving ratios are `4.00567,4.00141,4.00035`; the final coefficient
differs by `2.94e-5` relatively. Amplitudes `0.01,0.005,0.0025` give identical
distance/amplitude curves. The stationary-director control stays on the cone.

Additional checks: independent DOP853 free-flow error below `6.86e-13`;
manufactured driven-flow error below `4.80e-14`; matrix normal-residual
identity below `3.88e-17`; normalized finite-difference metric error below
`3.82e-10`; conserved free energy and rotation charge. The repository's
response kernel agrees with an independent normalized-source ODE to
`1.30e-11` on its refined time grid.

### Post-review: exact departure, returns and eigenline motion

The user review supplied a stronger closed-form check after the frozen
experiment had already passed locally. These results and their four new
probe checks are explicitly **post-review additions**; the original freeze
and original 22-check archive remain intact.

For the frozen `A_dot0=0` family, use the frame `(n0,vhat,n0 cross vhat)`
and write `c=cos(omega t)`, `s=|v0|sin(omega t)/omega`. Then

\[
 \beta=A_0\begin{pmatrix}2c/3&s&0\\s&-c/3&0\\0&0&-c/3\end{pmatrix},
 \quad \mu_3=-A_0c/3,\qquad
 \mu_\pm=A_0[c/6\pm\sqrt{c^2/4+s^2}].
\]

The minimum pairwise eigenvalue gap therefore gives the exact distance

\[
 d(t)=\frac{|A_0|}{\sqrt2}
       \left[\sqrt{c^2/4+s^2}-\frac{|c|}{2}\right].
\]

The implementation rationalizes the bracket as
`s^2/(sqrt(c^2/4+s^2)+|c|/2)` near returns. For nonzero `A0,v0`, the field
returns to the cone **exactly at** `omega t=k pi`, with
`beta(t)=(-1)^k beta0`, and is biaxial between those isolated returns.
The frozen quadratic departure law is the small-time expansion of this
formula. Return at discrete times does not make the family dynamically
invariant on an interval.

The review also proposed interpreting the fitted director's angle as
libration, and strengthening the result to exclude monotone angular
advance. That does not follow. Once the field is biaxial, its nearest
uniaxial approximation is not an actual uniaxial field director. At `c=0`
the two extreme eigenvalues have equal absolute value: two perpendicular
axes, with opposite fitted amplitude signs, are equally good minimizers.
The nearest fit switches between those branches as `c` changes sign.
Folding that switched axis into an acute angle does not describe a
continuous director trajectory.

For `A0>0`, a continuously continued eigenline in the same plane has

\[
 \alpha(t)=\tfrac12\operatorname{unwrap}\arg[c+2is],\qquad
 \dot\alpha=\frac{|v_0|}{c^2+4s^2}>0.
\]

At phases `omega t/pi=0,0.25,0.375,0.5,0.625,0.75,1`, the continued angle
is approximately `0,7.897,17.163,45,72.837,82.103,90` degrees. It passes
45 degrees rather than reversing. At a repeated eigenvalue its continuation
is a choice within an eigenspace; it must not be identified with the unique
axis of the returned uniaxial tensor. Thus the free field can have an
advancing eigenline while failing to be an autonomous uniaxial rotor.

The new controls compare exact eigenvalues and distances with full matrix
evolution over a period for all three speeds and both amplitude signs;
verify periodic returns; check the continuous eigenline equation; and
exhibit both distinct minimizers at the nearest-axis switch. The physical
obstruction remains the missing shape equations, not a prohibition of
eigenframe rotation.
Across those controls, the spectrum, distance, return and eigenline
residuals are all below `3.5e-18` in the stated amplitude units.

## 6. What this supplies for selection and causality

The field has an existing anisotropic source coupling. On the candidate
family its energy is `-CA n^T S n`. For an axisymmetric prescribed source
`S=s(mm^T-I/3)`, this is a `CA s sin^2(chi)` alignment energy plus a constant,
where `chi` is the angle between the director and source axis. This is a
field-derived quadrupolar dependence, but `chi` has not been identified
with the triangle holonomy angle of #283. No analyzer itinerary or its
frame comparison enters this source term automatically.

A directional tensor contraction is

\[
 m^\mathsf T\beta m=A[(m\cdot n)^2-1/3].
\]

It is even under the exact director identification. That alone gives no
readability protection: #284/#286 already show informative even functions.
This contraction identifies the sort of relational field quantity an
apparatus could couple to, not a completed local measurement protocol
compatible with the two-boundary history problem. The homogeneous tensor
also is not a source-localized mouth solution merely because it has a
director. There is no inherited map from the closure-conditioned ensembles
to these tensor histories.

| Question | Result |
|---|---|
| An explicit director-to-existing-field embedding | Yes; uniaxial TT tensors, with `n ~ -n` and a singular apex |
| Quadratic field action and its normalization | Derived in the inherited linear ESU approximation |
| Round angular kinetic term after substitution | Yes; inertia `2CA^2` |
| Autonomous free triangle rotor in that family | Obstructed by two nonzero omitted shape equations |
| Rotating uniaxial fields under every possible drive | Not excluded; manufactured forced control works |
| Existing anisotropic source coupling | Quadrupolar; not an identified triangle-holonomy energy |
| Canonical preparation, history weight `Phi`, Born law | Not derived |
| Operational BAM source readout or universal non-readability | Not derived |

The next field calculation must retain the extra shape dynamics or derive
a source/confinement mechanism that controls them. Substituting the desired
rotor into an action and ignoring its normal equations would miss exactly
the obstruction established here.

## Reproduction and provenance

```bash
python -m pytest -q tests/test_tt_triangle_rotor.py
python -m experiments.closure_ledger.tt_triangle_rotor_probe \
  --output-dir experiments/closure_ledger/runs/20260906_tt_triangle_rotor_reviewed
```

**27 tests, 22/22 frozen checks and four post-review checks pass.** The physical verdict is
`FREE_ROTATING_UNIAXIAL_TT_FAMILY_NOT_INVARIANT`. A failed required numerical
check returns `UNRESOLVED` and a nonzero probe exit; it cannot be promoted
to a no-go. See the [reviewed report](../experiments/closure_ledger/runs/20260906_tt_triangle_rotor_reviewed/probe.md)
and [full reviewed archive](../experiments/closure_ledger/runs/20260906_tt_triangle_rotor_reviewed/probe.json).
The [original frozen run](../experiments/closure_ledger/runs/20260906_tt_triangle_rotor_probe/probe.json)
is preserved separately with its original 22 checks.

The integrated focused suite passes **210 tests**, covering this module,
positive counting, pointer spread, source readout, the five preceding
closure/history modules and six relevant inherited tensor tests.

An integration issue was resolved before publication: #285 had been merged
into #284's branch after `main` received #284's earlier head. Consequently
#286's `cd7f300` did not yet include the pointer-spread files. The integration
commit brings in #285's merge commit `7ec9ee7` and its exact files, retaining
round 9's additions. No prior probe output or preregistration was changed.
After #286 merged, `dcf6716` was incorporated before publishing the TT
implementation. It has exactly `cd7f300`'s tree: round 9 was already in the
baseline, while the pointer-spread integration remains needed at that main
revision. The TT implementation and frozen experiment were already complete
locally when the user review arrived; the exact-period and eigenline work
followed that review.

No physical parameter, hypothesis or numerical tolerance changed after the
freeze. The first symbolic run encountered a truncated pre-existing SymPy
installation; a fresh dependency installation resolved that environment
failure before the symbolic checks ran. No implementation reconstruction or
post-result retuning was required. The numerical seed `20260907` is the
frozen seed, not the execution date; the archive was produced on 2026-09-06.
