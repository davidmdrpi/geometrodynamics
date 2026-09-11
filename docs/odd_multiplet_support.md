# Odd multiplet ESU support and its preparation kernel

A complete odd multiplet of independent real conformal scalars supports the
exact round Einstein static universe. The supplied background construction
is confirmed. The new preparation result is degree dependent: **the fixed-
trace Gram matrix is isolated at k=1 and k=3, but at k=5 there is an exact
84-dimensional family of changes invisible to the entire stress tensor.**
Equal coefficients are therefore not universally necessary in this class.
Neither an isolated preparation nor an equal-stress family establishes
attraction, dynamical stability, or a preparation probability.

The [freeze](odd_multiplet_support_prereg.md) was published at `a12cb2f`.
The review amendment `9c33867` added componentwise parity, twelve named
required checks, their total verdict mapping, and the two-phase correction
to the component-bound proof **before implementation or a P3 scan**.
P1/P2 were prior candidate results supplied by the user. The ranks and
susceptibilities below were first computed after that amendment. The
analytic explanation of equality of the two kernels was derived after
seeing those ranks and is labelled separately below.

## 1. Independent components and exact pointwise support

Use the action and normalization in the freeze:

\[
 I=\int\sqrt{-g}\left[\frac{R-2\Lambda}{2\kappa}
 -\frac12\sum_{I=1}^{m}\left((\nabla\phi_I)^2+\frac R6\phi_I^2\right)\right],
 \qquad ds^2=-dt^2+a^2d\Omega_3^2.
\]

Let \(Y_i\) be a real Haar-orthonormal basis of degree k harmonics and
\(N=(k+1)^2\). For odd k, every component obeys \(Y_i(-x)=-Y_i(x)\).
The addition theorem and its derivatives give, pointwise,

\[
 \sum_iY_i^2=N,\qquad \sum_iY_iD_aY_i=0,\qquad
 \sum_iD_aY_iD_bY_i=\frac{Nk(k+2)}{3a^2}\delta_{ab}.
\]

The first sum is invariant under SO(4) and hence constant on the sphere;
its Haar integral fixes N. The derivative tensor is isotropic under the
point stabilizer, and the harmonic equation fixes its trace. These are
pointwise statements. Spatial averaging is not the support mechanism.
The implementation additionally checks all their polynomial coefficients
over the rationals at k=1,3,5, independently of numerical orthonormalization.

Prepare independent fields
\(\phi_i=sY_i\cos\theta\), where \(\theta=(k+1)t/a+\delta\).
The improved stress, including \(G_{\mu\nu}\phi_i^2/6\), gives

\[
 \sum_iT_{0a}=0,\quad \sum_iT_{ab}^{\rm TF}=0,\quad
 \rho=\frac{Ns^2(k+1)^2}{2a^2},\quad p=\rho/3.
\]

The Einstein equations require

\[
 \Lambda=\frac{3}{2a^2},\qquad
 s^2=\frac{3}{\kappa N^2},\qquad
 \rho=\frac{3}{2\kappa a^2}.
\]

This is exactly the density and cosmological constant of
[#293's even homogeneous control](scalar_esu_support.md), with different
field content and spatial configuration. Equal background stress does not
establish equal perturbative response. If the harmonics instead use physical
normalization \(Y_i/\sqrt V\), the coefficient squared is \(Vs^2\),
where \(V=2\pi^2a^3\).

Each individual component has anisotropic stress. The coherent field
\(\phi=\sum_i\phi_i\) has cross terms and fails the pointwise background
conditions. That negative control is a regression of #293's theorem, not a
new exclusion. An even-degree basis is also tested: isotropic stress alone
cannot qualify a preparation as odd.

## 2. Kinetic regularity and the component bound

Writing \(f=\kappa F\),

\[
 f=1-\frac{\cos^2\theta}{2N},\qquad
 K_{IJ}=\frac{\delta_{IJ}}f+\frac{\kappa\phi_I\phi_J}{6f^2}.
\]

The full Einstein-frame kinetic matrix has eigenvalues \(1/f\) transverse
to the field vector and \(1/f^2\) along it. At simultaneous field zeros it
is the identity. Its positivity follows from \(f>0\), with minima

| k | N independent real components | minimum f |
|---|---|---|
| 1 | 4 | 0.875 |
| 3 | 16 | 0.96875 |
| 5 | 36 | 0.9861111111 |

This is a regularity/kinetic-sign test, not a stability theorem.

For the scoped lower bound, write \(\phi_I=f_I(x)\cos\theta\) with all
\(f_I\) in the same positive degree. At a field-zero phase, constant
positive density forces \(\sum_I f_I^2=r^2>0\) constant in space. Its
Hessian vanishes. At a phase with **nonzero cosine**, isotropic stress then
forces \(\sum_I D_af_ID_bf_I\) proportional to the metric; the eigenfunction
equation fixes it to \(k(k+2)r^2\delta_{ab}/(3a^2)\). This has rank three.
But \(\sum_I f_ID_af_I=0\), so three independent derivative vectors lie
in the (m-1)-dimensional subspace orthogonal to f. Thus \(m\ge4\), attained
by k=1. This does not classify arbitrary phases, mixed degrees, other
interactions, or the minimum component count at each higher k.

## 3. The complete preparation map

Use unit-radius derivatives temporarily. For a symmetric Gram perturbation
H in the Haar basis, define

\[
 Z=Y^THY,\qquad U_{ab}=(D_aY)^TH(D_bY),\qquad
 A_{ab}=\left(U_{ab}-\tfrac16D_aD_bZ\right)^{\rm TF}.
\]

For a general common-phase Gram matrix the same formulas give the full
stress, with the appropriate amplitude included in H:

\[
 \rho=\frac{(k+1)^2}{2}Z+\frac{\cos^2\theta}{12}\Delta Z,
 \qquad T_{0a}=-\frac{k+1}{6}\sin\theta\cos\theta\,D_aZ,
\]
\[
 T_{ab}=\tfrac13\rho\delta_{ab}+\cos^2\theta A_{ab}.
\]

Physical derivatives and frequency restore the radius factors explicitly in
code. These formulas are checked against independently factorized component
fields and the full improved stress, rather than used as their own oracle.

At the field-zero phase, zero density response forces Z=0. Thus the complete
stress kernel is exactly the kernel of the polynomial map \(H\mapsto(Z,A)\).
The anisotropy-only map is \(H\mapsto A\). A fixed-trace perturbation obeys
\(\langle Z\rangle=\operatorname{tr}H=0\). Although the rational harmonic
basis is not orthonormal, its exact moment metric imposes this same trace
condition before whitening.

All coefficient polynomials are homogeneous of degree 2k in four ambient
coordinates. A homogeneous polynomial vanishing on the unit sphere vanishes
everywhere by radial scaling; no quotient relation is omitted. Symmetrized
invariant derivatives give the covariant scalar Hessian. The integer
certificate uses the five independent components of 12A and the full Z.

| k | fixed-trace input dimension | rank A | nullity A | rank L | nullity L |
|---|---|---|---|---|---|
| 1 | 9 | 9 | 0 | 9 | 0 |
| 3 | 135 | 135 | 0 | 135 | 0 |
| 5 | 665 | 581 | 84 | 581 | 84 |

The archived integer kernel columns annihilate every exact coefficient row.
Their independence and a matching rank lower bound are verified modulo
1000003. Together these establish the rank over Q: neither modular rank
alone nor small numerical singular values establish it. Exact fraction-free
row reduction regenerates the certificates; subsequent runs independently
verify the stored witnesses against rebuilt operators.

For any certified nonzero H in the full kernel,

\[
 G(\epsilon)=s^2(I+\epsilon H),\qquad
 |\epsilon|\|H\|_{\rm op}<1
\]

is positive definite and has **exactly the same full stress for all x and t**.
This is an 84-dimensional affine family near the degree-5 preparation.
Every archived kernel basis direction is checked by factoring G into actual
independent fields at the frozen positive and negative amplitudes. Internal
field rotations leave G unchanged; these nonzero H are not that redundancy.
Since Z=0, the summed field square and f are unchanged along these families,
so kinetic positivity also persists. Their response away from the round
metric has not been derived.

## 4. Post-freeze explanation: the two kernels coincide

The numerical rank equality prompted a further analytic derivation; this
was not a prediction in the freeze. On the unit sphere,

\[
 D^a A_{ab}=-\frac1{36}D_b\left(\Delta Z+12(k+1)^2 Z\right).
\]

It follows from \(\operatorname{tr}U=k(k+2)Z+\Delta Z/2\), the scalar
Hessian divergence identity, and Ric=2g; equivalently it is the spatial
stress-conservation equation. The implementation checks the complete
polynomial identity, not merely selected kernel witnesses.

If A=0, the expression in parentheses is constant. Z has harmonic degrees
at most 2k, and

\[
 12(k+1)^2-\ell(\ell+2)>0\quad(0\le\ell\le2k).
\]

Therefore Z is constant. Its fixed zero mean forces Z=0. Consequently
\(\ker A=\ker L\) throughout this common-phase, single-degree class.
This explains the measured equality without asserting uncomputed ranks at
other degrees.

## 5. Finite susceptibility, not a dynamical selection law

Use the frozen input norm \(\|H\|_F/\sqrt N\) and output norm
\(\langle\delta\rho^2+2|\delta j|^2+\|\delta T_{ab}\|_F^2\rangle^{1/2}/\rho_0\).
Time averaging is used only to define this norm. The stress remains a
pointwise function in all equality tests. Exact polynomial moments and
exact Fourier time averages determine the norm matrices; only whitening
and positive singular values use floating-point linear algebra.

| k | smallest positive A singular value | largest A | smallest positive L | largest L |
|---|---|---|---|---|
| 1 | 0.3042903097 | 0.3042903097 | 0.6526300069 | 0.6526300069 |
| 3 | 0.1463850109 | 0.3750000000 | 0.1677050983 | 0.7318798728 |
| 5 | 0.0857240833 | 0.4157397096 | 0.0857240833 | 0.7474694391 |

Certified zero singular values are reported as exact zeros; residual
eigenvalues of the numerical moment Gram matrix are separately recorded.
The frozen relative singular reporting cutoff is 1e-10. It does not decide
rank or alter the exact witnesses. The trace/amplitude direction, excluded
from the table, has anisotropy zero and full-source norm sqrt(4/3).

For k=1 and k=3 these spectra give finite lower and upper bounds on every
fixed-trace Gram mismatch. At k=5 the lower bound applies to the component
orthogonal to the certified kernel. A diagonal amplitude scan explores
only one basis-dependent slice. The archive records amplitude mismatch,
induced Gram mismatch, source response, and both signs of the three frozen
step sizes. It also records all singular extreme directions, twenty seeded
random directions per degree, and every kernel witness.

P1 uses both (8,16) and (12,24) Hopf rules plus off-grid points. Squared
stress at k=5 has degree 20, so the freeze's illustrative angular order 16
is insufficient for its norm. Norm verification uses (12,24) and (16,32),
with exact polynomial moments as an independent route. No rank conclusion
uses a quadrature fit.

## 6. Reproduction and remaining scope

Run from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.odd_multiplet_support_probe \
  --output-dir experiments/closure_ledger/runs/20260910_odd_multiplet_support
pytest -q tests/test_odd_multiplet_support.py
```

Use `--recompute-certificates` to repeat exact row reduction instead of
verifying the archived integer witnesses. SymPy's optional python-flint
backend accelerates regeneration; it is not needed for verification.

The [archive](../experiments/closure_ledger/runs/20260910_odd_multiplet_support/)
contains JSON/Markdown results and separate exact kernel certificates.
Every required gate is exercised as missing and false, including CLI
failure and stale-report replacement. A computation exception also writes
an unresolved report and exits nonzero.

The construction extends the matter action by independent components;
it does not derive their count, common phase, or preparation. A static
support with positive kinetics does not establish coupled stability.
Background gradients can contribute tensor stress under perturbation, so
#292's bare tensor frequency is not automatically transferable. Using the
same components for support and signal creates stress cross terms; adding
signal fields changes the field-content assumption. The coupled support
response is required before an O(s^4) persistence calculation. Phi selection
and the source-local causality gate remain open.
