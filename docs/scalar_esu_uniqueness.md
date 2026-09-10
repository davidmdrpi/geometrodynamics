# Post-freeze extension: the homogeneous scalar support is unique

This result was derived **after** freeze `de55f3f`. It closes the
unrestricted-history classification left open in that freeze and in the
first implementation of [the scalar-support audit](scalar_esu_support.md).
The original freeze and its nine gates are unchanged. This supplement has
its own exact certificates, numerical controls and fail-closed verdict.

Within the stated action and exact round ESU geometry, **every nontrivial
smooth support by one real conformal scalar alone is homogeneous**, with

\[
 \phi=\sqrt{3/\kappa}\cos(t/a+\delta),\qquad
 \Lambda=3/(2a^2).
\]

The imposed odd sector therefore contains no such support. Without that
condition, the even control is the full family, up to its time phase.
This does not classify other matter contents, geometries or averaged
models, or select any ensemble, history weight or readout.

## The three coefficient conditions

The spatial-isotropy argument already proves that any nontrivial slice is
nowhere zero and has \(h=1/\phi=A+B\cdot x\), with \(|A|>|B|\).
By compactness and smoothness, a nowhere-zero slice has a time neighborhood
on which the whole sphere remains nowhere zero. The unique coefficients
A(t) and B(t) are smooth there, so differentiating them is legitimate.

With unit ambient coordinate x, the scalar equation gives

\[
 h^3(\Box-a^{-2})(1/h)=C+L\cdot x+x^TMx,
\]
\[
\begin{split}
 C&=A\ddot A-2\dot A^2+(2|B|^2-A^2)/a^2,\\
 L&=A\ddot B+\ddot A B-4\dot A\dot B+AB/a^2,\\
 M&=\operatorname{sym}(B\ddot B^T)-2\dot B\dot B^T,
 \quad \operatorname{sym}(UV^T)=(UV^T+VU^T)/2.
\end{split}
\]

For this polynomial to vanish on the unit sphere, its odd part requires
L=0. Its remaining quadratic form is constant on every unit vector, so
M=mu I4 and C+mu=0. These are the three conditions; setting C=0 before
disposing of mu would omit a sphere-constraint term.

## The rank argument and its zero-vector exception

Put V=dot B and W=ddot B. The exact factorization

\[
 M=H\begin{pmatrix}0&0&1/2\\0&-2&0\\1/2&0&0\end{pmatrix}H^T,
 \qquad H=[B\ V\ W]\in\mathbb R^{4\times3},
\]

proves rank(M)<=3, including every degenerate case. Equivalently, a nonzero
vector orthogonal to span{B,V,W} is in its kernel. A nonzero mu I4 has rank
four, so **mu=0**, and consequently C=L=M=0. The ambient dimension four
is essential to this particular rank step.

For **B!=0**, test M=0 on any w perpendicular to B:

\[
 0=w^TMw=-2(w\cdot V)^2.
\]

The vectors and field are real, so V=lambda B. Now M=0 becomes
sym(BW^T)=2 lambda^2 BB^T. Applying it to w perpendicular to B gives
\(B(W\cdot w)/2=0\), hence W is parallel to B as well. Its component
along B fixes

\[
 \dot B=\lambda B,\qquad \ddot B=2\lambda^2B.
\]

This proves the positive-semidefinite rank-one inference without relying
on an unstated matrix lemma. It also covers V=0: when B!=0 it forces W=0.

**At B=0 the conclusion is weaker.** The matrix is -2VV^T, so V=0, but
W is unrestricted by this one pointwise equation. For example
B=V=0, W=(1,0,0,0) gives M=0. This is a counterexample to an unqualified
pointwise assertion that all three vectors must be proportional to B.
It is not an admissible full-history counterexample. The uniqueness proof
below works on B!=0 intervals and never divides by B at a zero.

## Regularity rules out every B!=0 interval

On such an interval write B=b e, with b=|B|>0. The alignment condition
implies e is constant and b ddot b=2 dot b squared. Define q=A/b. The
longitudinal component of L=0, divided by b squared, is

\[
 \ddot q-2(\dot b/b)\dot q+q/a^2=0.
\]

Subtracting q times that equation from C/b squared=0 gives

\[
 \dot q^2+(q^2-1)/a^2=0,
 \qquad q^2+a^2\dot q^2=1.
\]

Thus |q|<=1. But smoothness on the complete sphere requires
|A|>|B|, hence |q|>1. This is an immediate contradiction. It does not
require integrating an inverse-linear b(t), extending that representation
past a pole, or assuming that a zero of B is generic.

If B were nonzero at any point of a regular reciprocal time patch, it
would be nonzero on a neighborhood, which has just been excluded. Therefore
B=0 on the entire patch and phi is spatially homogeneous there.

## Independent route from the momentum Einstein equation

There is a shorter independent check for *support*, since the exact ESU
also requires T0i=0. Direct substitution of phi=1/h in the improved stress
gives, without the wave equation,

\[
 T_{0i}=\frac23\dot\phi D_i\phi-\frac13\phi D_i\dot\phi
       =\frac{D_i\dot h}{3h^3}
       =\frac{D_i(\dot B\cdot x)}{3h^3}.
\]

The gradient of a linear ambient function vanishes everywhere on S3 only
when its coefficient vector vanishes. Thus **dot B=0** throughout any
regular time patch. This does not use the rank or alignment argument.

If the constant B were nonzero, L=0 would force
ddot A+A/a squared=0. Write A=A0 cos(t/a+delta), A0>=0, and b=|B|.
The remaining scalar numerator is exactly

\[
 C=\frac{2}{a^2}(b^2-A_0^2).
\]

Hence A0=b, so |A(t)|<=b and the reciprocal has a pole somewhere on the
sphere at every time. This contradicts smoothness again. The independently
verified residual therefore corroborates the rank route, with the same
conclusion B=0.

## Identically zero slices and the global history

No reciprocal representation is imposed on a slice where phi vanishes
identically. Instead choose any time at which the nontrivial history is
not identically zero. The preceding arguments give an open time interval
on which both phi and dot phi are homogeneous. At a time within that
interval every nonconstant scalar harmonic has zero field coefficient and
zero velocity. Each coefficient obeys the linear ESU oscillator equation,
so uniqueness makes it zero throughout the connected history. This extends
homogeneity through all simultaneous zero slices without dividing by phi.

The only alternative, the zero history, cannot supply the Einstein
enthalpy 2/(kappa a squared). The homogeneous wave equation and the two
background Einstein equations then fix the amplitude and Lambda quoted
above. No antipodal condition was needed for this uniqueness theorem.

## What was checked

- Symbolic expansion into C, L and M, including the spatial sphere terms.
- Exact rectangular factorization of M, transverse-square and alignment
  equations in a B-adapted orthogonal frame.
- Exact q identity giving the regularity contradiction.
- Independent improved-stress momentum identity and the b squared minus
  A0 squared residual.
- Direct off-shell reciprocal field jets at 120 sphere points across three
  radii, compared with the polynomial numerator and full improved stress.
- Generic rank-three matrices, exact aligned solutions, transverse velocity
  and acceleration failures, and B=dot B=0 with nonzero ddot B.
- Separate supplementary verdict gates, with missing/failed and unknown-key
  controls. The topological/continuation steps remain written mathematics
  for review; these checks are not formal verification of the whole proof.

```bash
python -m experiments.closure_ledger.scalar_esu_uniqueness_probe --output-dir experiments/closure_ledger/runs/20260909_scalar_esu_support
python -m pytest -q tests/test_scalar_esu_uniqueness.py tests/test_scalar_esu_support.py tests/test_esu_support_response.py
```

The original probe archive is preserved. The extension is recorded separately
in [uniqueness.md](../experiments/closure_ledger/runs/20260909_scalar_esu_support/uniqueness.md)
and [uniqueness.json](../experiments/closure_ledger/runs/20260909_scalar_esu_support/uniqueness.json).
