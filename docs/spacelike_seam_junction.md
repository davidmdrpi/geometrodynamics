# Deflected spacelike seam: leading junction sources and null limit

**The deflection matches the induced metric, but requires both surface stress
and a scalar surface source.** Its large unit-normal angular stress has a finite
coordinate-densitized null limit; that distinction prevents interpreting a
normalization singularity as a new general no-go theorem. No supporting action
or finite-thickness metric has been constructed.

This follows #335's 2D packet and leading Einstein response. Its original files
are unchanged. The [validation protocol](spacelike_seam_protocol.md) was written
after analytic exploration; this is not a blind preregistration. The calculation
keeps only the first nonzero order of each quantity. Near a null seam those
orders differ, so an apparent finite-amplitude Israel tensor obtained by mixing
known and unknown orders would not be justified.

## 1. Matched embeddings

In common parameter V write the incoming and outgoing embeddings as

\[
 (U_-,V_-)=(x(V),V),\qquad (U_+,V_+)=(g(V),y(V)),\quad g=-c/V.
\]

x,y are O(A²). Set K=2/(b e²), F_h=4/e²=2bK and

\[
 S={\kappa_5 b\over3}(\Phi_1')^2,\quad
 I=\int S\,ds,\quad B=\int sS\,ds,\quad C=\int s^2S\,ds,
 \quad L=\int s\log(V/s)S(s)\,ds.
\]

All integrals begin in the vacuum past V0. The two horizon radius responses are
rho_-=-VI+B and rho_+=-B+C/V. Matching the angular and longitudinal components
of the first fundamental form gives

\[
 \rho_- -\rho_+=K(Vx-gy),\qquad x'=g'y',
\]
\[
 x={-VI+B+L\over KV},\quad x'=-{L\over KV^2},\quad
 y={V\over c}\left({\rho_- -\rho_+\over K}-Vx\right).
\]

The common metric is h=H dV²+R² dOmega3² with

\[
 H=-F_h x'={2bL\over V^2}>0,\qquad
 R=b+\rho_- -KVx,\quad R'=-B/V.
\]

It is spacelike wherever L>0. No positive threshold or division-by-zero fix
extends the unit-normal description onto the initial null segment L=0.
The archived CSV uses NaN there for quantities that are undefined.

The deflection changes the retained domains. With the computed signs x<0,y<0,
the outgoing point has g<0,y<0, on the **white-hole side** of the unperturbed
horizon. Thus this is not a gluing confined to the two original exterior
rectangles. The leading calculation uses regular horizon Taylor data from those
rectangles; no full interior evolution or global extension is claimed.

## 2. Future timelike normals and second fundamental form

The normals, normalized by n.n=-1 and pointing from the retained past to future,
are, to the required leading orders,

\[
 n_-={(-x',1)\over\sqrt H},\qquad
 n_+={(g',-y')\over\sqrt H}.
\]

Use K_ab=e_a^mu e_b^nu nabla_mu n_nu. Let k_s be its eigenvalue in the unit
longitudinal tangent direction, and k_a the eigenvalue in each angular direction.
The embedding-acceleration formula gives

\[
 k_s^-={1-B/(2L)\over V\sqrt H}+O(A),\qquad
 k_s^+={-1-B/(2L)\over V\sqrt H}+O(A).
\]

The angular calculation must include the O(A²) bulk radius response; evaluating
only the unperturbed metric would miss contributions of the same order:

\[
 k_a^-=-{B+2L\over bV\sqrt H}+O(A^3),\qquad
 k_a^+={-B+2L\over bV\sqrt H}+O(A^3).
\]

Consequently

\[
 [k_s]=-{2\over V\sqrt H}+O(A),\qquad
 [k_a]={2V\sqrt H\over b^2}+O(A^3).
\]

These expressions are tested independently using -n.acceleration/H and n(r)/r.
Metric-response corrections to the longitudinal curvature occur at its omitted
O(A) order; we do not infer that coefficient from the leading embeddings.

## 3. Required spacelike Israel source

For this convention, a spacelike interface in Gaussian normal coordinates has
ds²=-d tau²+h_ab dx^a dx^b. Its distributional Einstein tensor gives

\[
 \kappa_5 S^a{}_b=[K^a{}_b]-\delta^a{}_b[K].
\]

The four surface eigenvectors are spatial. S_s is **not a surface rest energy
density** of a timelike material shell. The leading required eigenvalues are

\[
 \kappa_5 S_s=-{6V\sqrt H\over b^2}+O(A^3),\qquad
 \kappa_5 S_a={2\over V\sqrt H}+O(A).
\]

In particular S_s<0: the distributional stress violates the ambient NEC in the
null direction n+e_s. The earlier null shell already required nonordinary stress;
the deflection does not remove the need for a specified support law. Do not
append the known -2[k_a] term to S_a and label the result accurate at O(A): the
unknown subleading [k_s] contributes at exactly that same order.

For b=c=kappa5=1, A=.01, V=4, N=512:

| Quantity | Leading value |
|---|---:|
| x | -7.17921210e-5 |
| y | -7.04904436e-4 |
| H | 2.13927982e-5 |
| Shared R-b | -1.71142386e-4 |
| [k_s] | -108.10261688 |
| [k_a] | +.03700187947 |
| S_s | -.11100563841 |
| S_a | +108.10261688 |

## 4. Scalar junction: common trace does not fix the normal derivative

Let P=partial_U Phi_- on the incoming horizon, Q=partial_Vout Phi_+ on the
outgoing horizon. For the spherical mode their values follow from integrating
u_UV=-Q_h u, Q_h=3 exp(-2)/b², with zero data on the second characteristic:

\[
 P={Q_h\over b^{3/2}}\left(Va-\int a\,ds\right),\qquad
 Q={Q_h\over b^{3/2}}\left(ga-\int g'a\,ds\right).
\]

Simply sampling the old fields on the moved surfaces gives an O(A³) scalar trace
discrepancy. We impose the necessary common trace at that order symbolically;
we do not mistake that discrepancy for an independent leading obstruction.
The normal derivatives can then be expressed in terms of that same trace:

\[
 (n\Phi)_-={\Phi_\Sigma'-2x'P\over\sqrt H},\quad
 (n\Phi)_+={\Phi_\Sigma'-2y'Q\over\sqrt H},
\]
\[
 [n\Phi]={2x'\over\sqrt H}(P-Q/g')+O(A^4).
\]

The common unknown trace correction cancels. This determines the leading O(A²)
source without pretending to have evolved the O(A³) bulk scalar. In the stated
normal convention Box Phi has delta-source J=-[n Phi]. No scalar surface
coupling means [n Phi]=0, which the computed packet does not satisfy.

At V=4, the characteristic-integral result is [n Phi]=-3.98246062e-5.
Independent one-sided derivatives of the 2D fields give -3.94317729e-5 at N=512.
Their relative discrepancy decreases through .31485, .11437, .03528, .009864
at N=64,128,256,512. The signed delta-source is therefore +3.98246062e-5 at
leading order. This rejects **source-free scalar continuation**, not a model
that supplies and varies an explicit surface coupling.

## 5. Conservation: what is checked and what remains

The O(A) Codazzi terms cancel exactly:

\[
 (S_s)' -3(R'/b)S_a=0+O(A^3).
\]

The omitted product (R'/b)S_s and the bulk flux jump [T_nV]=[n Phi]Phi_Sigma'
are O(A³). Verifying their balance requires the omitted stress and metric terms;
this calculation does not claim full conservation at that order.
The leading O(1) Gauss jump likewise cancels:

\[
 \overline{k_s}[k_a]+\overline{k_a}[k_s]=0.
\]

Both identities are independently checked symbolically. Remaining O(A²) normal
stress balance is uncomputed. They establish consistency of the retained orders,
not existence of a surface action or a nonlinear junction solution.

## 6. Near-null limit and necessary thickness budget

For positive amplitudes, H scales as A², S_s as A, S_a as 1/A, and J as A².
However the coefficient relevant to the surface measure stays finite:

\[
 \sqrt H\,S_a\longrightarrow {2\over\kappa_5 V}.
\]

For the longitudinal component, the unit tangent is e_s=t/sqrt(H). Its
contribution to the tensor density is therefore
sqrt(H) S_s e_s^mu e_s^nu=(S_s/sqrt(H))t^mu t^nu, with
S_s/sqrt(H)=-6V/(kappa5 b²). These are exactly the null shell density and
pressure computed independently in #334, in the common generator normalization.
Thus both independent stress coefficients have the correct tensor-density null
limit, not just the angular component.

These finite coordinate-densitized coefficients do not demonstrate that
every component of a regularized spacetime stress or curvature stays bounded.
It explains why the unit-normal divergence alone cannot rule out a null
distributional limit. Onset from the null prefix still requires a formalism
that remains valid through changing causal character.

For a specified proper-normal smearing with support [-L_perp,L_perp] and
integrated angular source S_a, a necessary budget is
max |T_a| >= |S_a|/(2 L_perp). At fixed L_perp this grows as 1/A. A uniformly
bounded target stress would require L_perp to grow at least as 1/A in this
prescription. Meanwhile the large longitudinal extrinsic curvature makes the
controlled thin-tube scale shrink as A. Thus a uniformly bounded-stress,
perturbatively thin proper-normal mollification is not obtained by this limit.
Other smoothing prescriptions or a fully changed focus geometry require their
own equations; this estimate is not a universal finite-thickness impossibility.

## Reproduction

```
python -m experiments.throat_dynamics.spacelike_seam_probe
pytest -q tests/test_spacelike_seam.py tests/test_reciprocal_characteristic_packet.py tests/test_null_horizon_continuation.py
```

[Result archive](../experiments/throat_dynamics/runs/20261014_spacelike_seam/result.json)
records source hashes, parent evidence, refinement, amplitudes, and width budgets.
`seam_N512.csv` contains the full sampled leading profiles. Amplitude scaling is
an asymptotic identity, not a nonlinear simulation. The missing deliverable is
an explicit support/scalar action with its next-order conservation and a regular
finite-thickness evolution; matching the first fundamental form does not supply it.
