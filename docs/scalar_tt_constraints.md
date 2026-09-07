# Bounding the omitted scalar–TT constraint response

The leading degree-3 standing scalar admits bounded particular solutions of
both omitted Einstein constraints. In the frozen spatial gauge, **the scalar
metric distortion, even with its mean removed, exceeds an all-time upper
bound on the scalar-induced homogeneous TT distortion**. Both scale as the
square of the scalar amplitude. Amplitude reduction cannot separate them.

This result assumes zero supporting-matter energy and momentum perturbations
on independently considered CMC slices. It is not an Einstein evolution.
**The complete omitted scalar backreaction remains unbounded by this work**:
its calculation needs the supporting stress and metric evolution. The result
rules out justifying this particular omission by smaller metric amplitude;
it does not prove that the scalar force is larger or cannot cancel.

The [follow-up freeze](scalar_tt_constraints_prereg.md) was published at
[`e1706b1`](https://github.com/davidmdrpi/geometrodynamics/commit/e1706b166b24b69afceda85932e88a43b7fb61f9)
before implementation or new measurements, on reviewed PR #289 head `67940cb`.
The [original freeze](reciprocal_scalar_tt_prereg.md) and archive are preserved.
The exact harmonic-power certificate below was derived after this freeze;
it strengthens the comparison without changing any frozen acceptance gate.

## 1. Assumptions and inherited equations

[`waves/initial_data.py`](../geometrodynamics/waves/initial_data.py) derives

\[
 \delta R^{(3)}=2\kappa\delta\rho=-8(\Delta+3/a^2)u,
 \qquad (\Delta+3/a^2)u=-\kappa\delta\rho/4.
\]

The scalar spatial perturbation is `h^S_ij=4u delta_ij` in the background
orthonormal frame, while `h^TT_ij=2 beta_ij`. Since `Kbar=0`, extrinsic
curvature enters the Hamiltonian constraint quadratically. Keep the linear
metric response to the order-`phi^2` round-background stress; source corrections
such as `u phi^2` are beyond this solve. Retain the spatial mean.

| Ingredient | Status |
|---|---|
| Conformal scalar stress and smooth S3 harmonics | Inherited |
| Hamiltonian constraint operator | Inherited from `initial_data.py` |
| Zero support `delta rho` and `delta j` | Chosen on individual slices |
| Spatial conformal gauge; zero dipole kernel coefficients | Chosen |
| Constant mean curvature `delta K=0` | Chosen |
| Minimal longitudinal extrinsic-curvature response | Chosen particular solution |
| Supporting pressure/stress law and lapse/shift evolution | Unspecified |
| Full Einstein evolution, constraint propagation, scalar reaction bound | Not established |

The inherited initial-data text connects rigidity to separate scalar stress
conservation. That conservation alone does **not** establish that a supporting
fluid stays unperturbed in a changing metric. Rigidity is initial-data input
here, not an evolution theorem.

Use `K_ij=-Lie_normal(g_ij)/2`, `j_i=-T_0i`, and the York operator
`(LX)_ij=nabla_i X_j+nabla_j X_i-(2/3)g_ij div X`. The CMC momentum equation
is `div K^L=kappa j`, and `div LX=Delta X+(1/3)grad div X+Ric(X)`.
These conventions are given in
[Gourgoulhon's initial-data lecture, equations (1)–(2), (15), (18), (27), (34)](https://arxiv.org/abs/0704.0149).
The specialization and bounds below are derived here.

## 2. A compatible finite Hamiltonian inverse

For physical orthonormal harmonics, `Delta Y_l=-lambda_l Y_l`, where
`lambda_l=l(l+2)/a^2`. The four degree-1 modes are a genuine kernel and impose
`int rho x^A dV=0`. A sourced kernel cannot be fixed by discarding its coefficient.

Degree-3 phi and phi_dot are antipodally odd, so their quadratic energy is
even with degrees `0,2,4,6`. Dipole compatibility therefore holds exactly.
The smooth source can be solved on whole S3; the singular-mouth problem
addressed in `initial_data.py` is absent here, not cured by this argument.

Set the unsourced kernel coefficients to zero. For the remaining coefficients,

\[
 u_l=-\frac{\kappa\rho_l}{4(3/a^2-\lambda_l)},\qquad
 \|u\|\le\frac{\kappa a^2}{12}\|\rho\|,
\]
\[
 \frac{\kappa a^2}{180}\|\rho_\perp\|
 \le\|u_\perp\|\le\frac{\kappa a^2}{20}\|\rho_\perp\|.
\]

These apply to physical L2 or volume-normalized RMS norms. Perpendicular
means the mean is removed. The inverse denominators for degrees `0,2,4,6`
are proportional to `3,-5,-21,-45`; the upper constants are sharp at degrees
0 and 2. The lower bound uses the degree-6 cutoff, whereas the upper
inhomogeneous bound extends to compatible sources of all degrees at least 2.

## 3. The standing-wave momentum is a gradient

Set `V=2 pi^2 a^3`, `omega=4/a`, and

\[
 \phi=A(t)Y(x),\quad A=s\cos\omega t,\quad
 Y=\sqrt{8/V}\operatorname{Re}(x_0+i\,m\cdot\mathbf x)^3.
\]

The improved stress reduces to

\[
 \rho=\tfrac12(\dot A^2+\omega^2A^2)Y^2+\tfrac1{12}A^2\Delta(Y^2),
 \qquad j=\nabla J,\quad J=-\tfrac16 A\dot A Y^2.
\]

These follow from
`rho=(phi_dot^2+|grad phi|^2)/2-Delta(phi^2)/6+phi^2/(2a^2)` and
`T_0i=(2/3)phi_dot grad_i phi-(1/3)phi grad_i phi_dot`. Both are checked
against the inherited full stress routine at nonzero momentum times.

For `Y^2=sum h_l Y_l`, the spectral sources and CMC solution are

\[
 \rho_l=\left(s^2\omega^2/2-\lambda_l A^2/12\right)h_l,
 \quad K^L=L\nabla w,\quad
 w_l=\frac{3\kappa J_l}{4(3/a^2-\lambda_l)},\quad l=2,4,6.
\]

Set the irrelevant constant in w to zero. Commuting derivatives using
`Ric=2g/a^2` gives `div L grad w=(4/3)grad(Delta+3/a^2)w`.
Integration by parts gives

\[
 \|K^L\|^2=\frac83\sum_l\lambda_l(\lambda_l-3/a^2)|w_l|^2,
 \qquad \frac{\kappa a}{\sqrt{30}}\|j\|
 \le\|K^L\|\le\kappa a\sqrt{3/10}\|j\|.
\]

For RMS norms divide the spectral sum by V. The lower bound uses the cutoff;
the upper constant is sharp at degree 2.

The **leading standing wave** has no transverse-vector forcing. Its momentum
is orthogonal to all six Killing vectors because it is a gradient, and to
all four gradient dipoles by parity. Generic q,p need not have this property:
the control `p=D_1 q` has Killing charge `-p^T D_1 q=-0.1371428571` and is
rejected. A trace K cannot absorb a Killing charge on compact S3 either,
because Killing vectors are divergence-free. No general momentum completion
is inferred from the special standing wave.

## 4. Exact metric size comparison

Isolate the scalar-induced tensor with zero initial tensor data. For
`Omega^2=8/a^2`, its exact forced solution is

\[
 \beta_{\rm ind}(t)=\frac{\operatorname{tensor}(F(q_0))}{C}
 \left[\frac{1-\cos\Omega t}{2\Omega^2}
 +\frac{\cos2\omega t-\cos\Omega t}{2(\Omega^2-4\omega^2)}\right],
 \quad C=V/\kappa.
\]

Since `||F(q0)||=2 sqrt(6) s^2/a^2`, the cosine bound implies

\[
 \|h^{TT}_{\rm ind}(t)\|_F\le\frac{4\sqrt6}{7}\frac{\kappa s^2}{V}.
\]

An independent exact decomposition bounds the scalar response. Under
normalized S3 volume, `R=x0^2+(m.x)^2` is uniform on `[0,1]`, with uniform
phase in that two-plane. Then

\[
 Y^2=\frac4V[R^3+\operatorname{Re}(x_0+i\,m\cdot\mathbf x)^6],
 \quad R^3=\tfrac14P_0+\tfrac9{20}P_1+\tfrac14P_2+\tfrac1{20}P_3,
\]

where `P_k=P_k(2R-1)` has S3 degree `2k`. The extra degree-6 harmonic has
zero phase average and is orthogonal to the zonal terms. Integration gives

\[
 V\|h_l\|^2=1,\quad27/25,\quad1/5,\quad201/175\quad(l=0,2,4,6).
\]

For `r=cos^2(omega t)`, `u_2=kappa s^2(12-r)h_2/30`. Keeping degree 2
alone proves, for nonzero s,

\[
 \|h^S_\perp(t)\|_{\rm RMS}=4\sqrt3\|u_\perp(t)\|_{\rm RMS}
 \ge\frac{66}{25}\frac{\kappa s^2}{V}
 >\frac{4\sqrt6}{7}\frac{\kappa s^2}{V}
 \ge\|h^{TT}_{\rm ind}(t)\|_F.
\]

This compares a scalar **lower bound** with a tensor **upper bound** at
every time, with the mean already removed. Their ratio is
`231/(50 sqrt(6)) > 1.88`, independent of s, a and kappa. The exact polynomial
decomposition, eigenvalues, harmonicity and positive squared gap are symbolic
checks; numerical harmonic projection independently checks the powers.

Both induced metric responses scale as `s^2`, with prospective scalar
reactions at `s^3`. Their norm comparison does not determine the coefficients
or cancellations in the full scalar equation. These are spatial metric norms
in the frozen gauge, not gravitational energies or observable-record bounds.

The original primary history also includes an independently specified tensor
of amplitude 0.01; its sampled total TT metric norm ranges from `0.003693`
to `0.016376`. The induced-field inequality must not be claimed for that
different total tensor.

## 5. Measured bounds and checks

For `a=kappa=1`, `s=0.2`, `m=(1,2,3)/sqrt(14)`:

| Quantity | Continuous extremum or bound |
|---|---:|
| Mean u, retained | `-0.001350949115` |
| u RMS | `0.001558149892` to `0.001597314609` |
| u RMS without mean | `0.000776381076` to `0.000852261959` |
| Scalar metric RMS without mean, exact degree-2 lower bound | `0.005349758496` |
| Induced TT metric norm, all-time upper bound | `0.002836402286` |
| Longitudinal extrinsic-curvature RMS, maximum | `0.001475656292` |

The u extrema are obtained from the norm of an affine function of
`cos^2(omega t)`; longitudinal curvature is proportional to `sin(2 omega t)`.
These are continuous extrema rather than sampled maxima.

The two spatial rules use 2048 and 6912 points, retaining all 84 harmonics
of degrees `0,2,4,6`. Reconstruction also passes for random q,p energy data
and a rotated coherent mode. Maximum scaled Hamiltonian and momentum
residuals are `2.27e-15` and `1.91e-15`. The momentum check differentiates
the reconstructed tensor with invariant-frame derivatives and explicit
connection terms, independently of the inverse formula. Quadrature changes
the compared quantities by at most `1.52e-15`. Independent forced-oscillator
integration agrees with the tensor closed form to `2.49e-14`.

All 14 gates pass. Tests also check sharp constants, radius normalization,
charged and arbitrarily small sourced-dipole rejection, and failed-gate
report overwrite behavior. The original archive is preserved.

## 6. What prevents a complete backreaction bound

With responsive support, the source is `rho_phi+delta rho_support`. For
compatible inhomogeneous support data, the triangle inequalities give

\[
 \max(0,\|u_{\phi,\perp}\|-\kappa a^2\|\delta\rho_{{\rm support},\perp}\|/20)
 \le\|u_{{\rm total},\perp}\|
 \le\|u_{\phi,\perp}\|+\kappa a^2\|\delta\rho_{{\rm support},\perp}\|/20.
\]

The full-source upper bound uses `kappa a^2/12`. A compatible gradient
support momentum has the analogous upper constant `kappa a sqrt(3/10)`;
arbitrary support momentum needs its transverse solve and compatibility
check. Sourced dipoles or Killing charges cannot be absorbed by these bounds.

No supporting-response norm bound or evolution law is supplied in this ESU
channel. As constraint data, proportional support can cancel or reinforce
the scalar source. Neither is a constructed matter history. Nor can an
arbitrarily chosen lapse establish physical instability: lapse is slicing
data, and an observable comparison needs a specified clock.

`waves/backreaction.py` explicitly leaves the supporting fluid and sound
speed unspecified. `initial_data.py` provides rigidity on a slice, not a
pressure law. `tangherlini/dynamics.py` evolves a different 4+1, spherical,
minimally coupled problem and does not complete this conformal 3+1 ESU
channel. The missing metric time derivatives, supporting stress, and
constraint propagation cannot be supplied by independently solving each slice.

The next derivation must specify or derive the ESU support's dynamical response
and evolve the scalar metric and matter equations with propagated constraints.
That completion can test whether the omitted reaction cancels, dominates, or
remains controlled beside the TT reaction. The source-local causality/readout
gate remains downstream of it.

## 7. The other review observations

The n=1 generators obey `D_i D_j+D_j D_i=-2 delta_ij I` exactly, so every
homogeneous TT coupling vanishes. For every integer `n>1`, the coherent
source has the nonzero term `n(n-1)(m m^T-I/3)/a^2`; these multiplets are
not identically dark. Thus n=1 is the unique **nonconstant** dark multiplet;
n=0 is also dark. If a mouth-position variable were represented solely by
the n=1 scalar sector, this TT channel could not read it at this order.
The repository has not supplied that field map.

In the original reciprocal history the source is initially uniaxial to
relative distance `1.46e-15`. At t=2 its norm is `0.00412219` and relative
biaxial distance `0.00956152`. Across 401 samples its norm ranges from
`2.47e-6` to `0.195959`, and its maximum relative distance is `0.27166`.
These are diagnostics, not a success criterion or a constraint evolution.
The contraction `Q_m=m^T beta m` remains meaningful; an autonomous source
director is not derived.

The original probe's two literal homogeneous TT constraint zeros now come
from the inherited ADM first curvature variation and the invariant connection
contracted with a general STF rate tensor. Their computed exact zeros and
limited linear homogeneous scope are recorded in the new report.

## Reproduce

```bash
python -m experiments.closure_ledger.scalar_tt_constraints_probe
python -m pytest -q tests/test_scalar_tt_constraints.py tests/test_reciprocal_scalar_tt.py tests/test_tt_triangle_rotor.py
```

See the [archived report](../experiments/closure_ledger/runs/20260907_scalar_tt_constraints/probe.md)
and [JSON data](../experiments/closure_ledger/runs/20260907_scalar_tt_constraints/probe.json).
