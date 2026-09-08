# ESU support response: the first scalar reaction does not cancel

For the preparation frozen in this round, the omitted scalar metric response
produces a nonzero initial cubic correction in the supporting fluid's proper
time:

\[
 \int_{S^3}Y\,F_{\rm proper}(0)\,dV
 =-\frac{7976}{875}\frac{\kappa s^3}{Va^2},\qquad V=2\pi^2a^3.
\]

The induced homogeneous TT force is initially zero. These responses therefore
cannot cancel identically for this preparation. This moves beyond #289's
metric-norm comparison to an actual scalar acceleration coefficient.

The coefficient does not depend on the support's sound speed, because its
initial density and pressure perturbations are zero. Later scalar-metric and
fluid responses do depend on the constitutive law. **That law is still not
derived from BAM's field content.** The result is conditional on the support
class, preparation and clock specified below; it is not a universal no-go for
BAM, a full corrected scalar evolution, or an operational readout result.

The [freeze](esu_support_response_prereg.md) was published at
[`c07891d`](https://github.com/davidmdrpi/geometrodynamics/commit/c07891d573ab719ba90e73f8f3309c8aad674350)
on main `080c1cc`, including #289, before implementation and measurements.
Both coefficient values were predictions in the freeze. The parity-necessity
question is separate and is not addressed here.

## 1. What fixes the ESU support, and what does not

The background Einstein equations give

\[
 \kappa\rho_0=3/a^2-\Lambda,\quad
 \kappa p_0=\Lambda-1/a^2,\quad
 H_f:=\rho_0+p_0=2/(\kappa a^2).
\]

They determine the enthalpy required to support the curvature. They determine
neither `dp/d rho` nor an entropy/anisotropic-stress response. The existing
[`waves/backreaction.py`](../geometrodynamics/waves/backreaction.py) names
this missing supporting fluid explicitly; the fixed ESU in `two_wave.py`
does not supply its evolution. `initial_data.py` holds its density fixed on
a slice. The 4+1 minimally coupled spherical solver in `tangherlini/dynamics.py`
is a different system and does not close this 3+1 conformal ESU problem.

Use the following explicit completion family: a perfect fluid with no
anisotropic stress, no exchange with the conformal scalar, and adiabatic
linear pressure response `delta p_f=c_s^2 delta rho_f`. Keep `c_s^2` as a
parameter. The demonstrated interval is `[0,1]`; the primary value `1/3`
is a chosen control, not a derived equation of state. The pressure derivative
does not need to equal the background ratio `p_0/rho_0`.

| Ingredient | Status |
|---|---|
| Existing conformal scalar and ESU geometry | Inherited |
| Required background fluid enthalpy | Derived from Einstein |
| Perfect-fluid, adiabatic perturbation class | Assumed completion family |
| Sound speed / microscopic BAM support | Not selected or derived |
| Initially unperturbed support and zero free tensor data | Chosen preparation |
| Linear scalar Einstein equations and fluid conservation | Derived and independently evolved |
| Initial proper-time cubic coefficient | Exact within these assumptions |
| Full corrected scalar field, apparatus, early-record law | Not derived |

Separate conservation does not make the support dynamically rigid. Its
density and velocity are allowed to respond here, and do so in the computed
histories.

## 2. Preparation and scalar Einstein equations

Use Newtonian gauge for degrees at least 2, with zero scalar shift and shear:

\[
 ds^2=-(1+2\alpha)dt^2+(1-2\psi)\bar g_{ij}dx^i dx^j.
\]

Set the homogeneous lapse `alpha_0=0`, fixing the background clock, and retain
the homogeneous spatial scale perturbation. The scalar is
`phi_1=A(t)Y`, `A=s cos(omega t)`, `omega=4/a`, with the same normalized
degree-3 coherent Y as #289. Initially the fluid density/pressure perturbations
and velocity vanish; `psi_dot=0`. All free tensor initial data vanish.

The initial spatial metric is not zero: the Hamiltonian constraint requires
`psi(0)=-2u(0)` from #289. The initial slice is CMC. Its subsequent Newtonian
gauge evolution is not the sequence of independently rigid CMC slices in
the previous round. The fluid's proper time is not the coordinate time t.

Let `lambda_l=l(l+2)/a^2`, `L_l=3/a^2-lambda_l`, and expand `Y^2` into
degrees `0,2,4,6`. Use `H_l=(Y^2)_l`, so a group amplitude multiplies a whole
harmonic projection, not one arbitrarily chosen basis vector. The exact powers
are `V int H_l^2 dV = 1,27/25,1/5,201/175`.

The improved stress has

\[
 \rho_\phi=\tfrac12(\dot A^2+\omega^2A^2)Y^2
             +\tfrac1{12}A^2\Delta(Y^2),\quad
 j_\phi=\nabla J_\phi,\quad J_\phi=-A\dot A Y^2/6,
 \quad p_\phi=\rho_\phi/3.
\]

Define the scalar anisotropic-stress potential by
`pi_phi^S=(Hess-gbar Delta/3)Pi`. For every nonconstant source mode,

\[
 \Pi_l=-\frac{3\dot J_{\phi,l}+\rho_{\phi,l}}{2L_l}
       =-\frac{A^2(\omega^2-\lambda_l/12)}{2L_l}(Y^2)_l.
\]

The first form follows from source momentum conservation and the divergence
of a scalar STF Hessian. It is independently checked by projecting the full
pointwise improved spatial stress onto all 83 nonconstant even scalar
harmonics. The remaining tensor anisotropic stress has not been set to zero.

A direct metric/connection variation on S3 gives

\[
 2(\Delta+3/a^2)\psi=\kappa(\rho_f+\rho_\phi),\quad
 -2\nabla\dot\psi=\kappa(j_f+j_\phi),\quad
 \psi_l-\alpha_l=\kappa\Pi_l\ (l\ge2),
\]
\[
 \ddot\psi-\psi/a^2-\tfrac13\Delta(\psi-\alpha)
       =\tfrac\kappa2(p_f+p_\phi).
\]

Substituting the pressure law reduces the scalar response to four forced
oscillators, retaining all source multiplets:

\[
 \ddot\psi_l=(a^{-2}+c_s^2L_l)\psi_l
       +\tfrac\kappa2(\tfrac13-c_s^2)\rho_{\phi,l}
       -\tfrac\kappa3\lambda_l\Pi_l.
\]

The last term is absent at degree 0. Recover the support by the constraints,
or independently evolve

\[
 \dot\rho_f=-\Delta J_f+3H_f\dot\psi,\qquad
 \dot J_f=-c_s^2\rho_f-H_f\alpha\quad(l\ge2).
\]

The independent run evolves fluid conservation and the spatial Einstein
equation. It does not reconstruct the fluid from the constraints, so its
Hamiltonian and momentum residuals are meaningful unused-equation checks.
The constant momentum potential is immaterial; its Euler equation is not
imposed. The homogeneous fluid energy equation remains active.

Momentum compatibility has its own reason here. The coherent standing wave
has `q=A q_dir`, `p=dot(A) q_dir` and `j_phi=grad J_phi`. The fluid starts at
rest, and its linear Euler forcing is a gradient, so `j_f=grad J_f` remains
true throughout this scalar response. For every one of the six round-S3
Killing fields K, compactness and `div K=0` give

\[
 \int K\cdot(j_\phi+j_f)\,dV
 =-\int(J_\phi+J_f)\operatorname{div}K\,dV=0.
\]

This is a gradient-current argument, independent of parity. General
parity-pure data can carry Killing charge: the inherited constraint tests
include and reject that negative control. The unused momentum residuals
above check local propagation separately. The homogeneous linear operator
and the leading standing-wave source keep this order-s^2 response within
degrees 0,2,4,6; hence no degree-1 Hamiltonian source is generated here.
Since `l(l+2)=3` only at `l=1`, additional source multiplets would require
their own dipole check and their own momentum accounting. Parity purity is
not asserted to be necessary for solvability. At later times this gauge is
not CMC, so the initial CMC conformal-Killing test is not imposed as a new
evolution condition.

The free-support scalar frequency is
`Omega_l^2=[c_s^2(l(l+2)-3)-1]/a^2`. This recovers the established degree-2
threshold `c_s^2=1/5`; the homogeneous scale mode remains unstable. This is a
cross-check against [Barrow et al., equation (15)](https://arxiv.org/abs/gr-qc/0302094),
not a new stability claim. No mean subtraction is used to remove it.

## 3. Coordinate acceleration and proper acceleration

Varying the conformal wave operator gives the coordinate-time cubic force

\[
 F_t=2\alpha\ddot\phi_1+(\dot\alpha+3\dot\psi)\dot\phi_1
       +2\psi\Delta\phi_1+\nabla(\alpha-\psi)\cdot\nabla\phi_1
       -\tfrac16\delta R^{(4)}\phi_1,
\]
\[
 \delta R^{(4)}=4\Delta\psi+12\psi/a^2-6\ddot\psi-2\Delta\alpha
              =\kappa(1-3c_s^2)\rho_f.
\]

The second equality is Einstein's trace equation: the leading conformal
scalar has zero trace, and Lambda is fixed. It is checked against the
geometric curvature expression, not substituted to make the check vanish.

Write the fluid coordinate velocity as `v^i=dx^i/dt=O(s^2)`. Normalization
gives `U^0=1-alpha+O(s^4)`, `U^i=v^i+O(s^4)`. Applying the material
derivative twice, with all displayed spatial contractions in the round
metric, gives

\[
 D_U^2\phi=\partial_t^2\phi
 -2\alpha\ddot\phi_1-\dot\alpha\dot\phi_1
 +\dot v\cdot\nabla\phi_1+2v\cdot\nabla\dot\phi_1+O(s^5).
\]

Initially the fluid is at rest, its pressure gradient vanishes and its
four-acceleration is zero. The spatial geodesic equation still gives
`v_dot=-grad alpha`: the fluid moves away from fixed spatial coordinates.
Also `dot(phi_1)(0)=0` and `delta R^(4)(0)=0`, since the initial support
perturbations vanish. Thus the conversion of the *cubic correction* is

\[
 F_{\rm proper}(0)=F_t(0)-2\alpha(0)\ddot\phi_1(0)
                         -\nabla\alpha(0)\cdot\nabla\phi_1(0)
 =2\psi(0)\Delta\phi_1(0)-\nabla\psi(0)\cdot\nabla\phi_1(0).
\]

Both correction terms matter. In the same units as the modal coefficients,

| Contribution to the initial conversion | Coefficient |
|---|---:|
| Coordinate-time correction | `+55096/875` |
| Clock normalization, `-2 alpha ddot(phi_1)` | `-76928/875` |
| Fluid coordinate acceleration, `-grad(alpha).grad(phi_1)` | `+13856/875` |
| Proper-time correction | **`-7976/875`** |

Equivalently, on this initial slice,
`D_U^2 phi=Delta_g phi-R^(4)phi/6`. The conformal spatial Laplacian variation
is `delta Delta phi=2 psi Delta phi-grad psi.grad phi`, giving the same
result without choosing a spatially varying lapse.

The same-metric numerical control strengthens the already frozen clock
gate without changing the freeze or its predictions. In the nonlinear test
metric `N=exp(epsilon alpha)`, `g_ij=exp(-2 epsilon psi)gbar_ij`, use
`U^0=N^-1`, `U^i=0` initially and construct `Gamma^mu_00` from the metric
derivatives. The geodesic scalar Hessian then gives, at the same event,

\[
 D_U^2\phi=N^{-2}\bigl(\partial_t^2\phi
              -\epsilon\dot\alpha\dot\phi\bigr)
       -\epsilon e^{2\epsilon\psi}\nabla\alpha\cdot\nabla\phi.
\]

Insert the coordinate acceleration obtained from that metric's full
conformal wave operator, including its scalar curvature. Differentiate both
clock results with respect to epsilon at zero, using centered differences
and Richardson extrapolation at the frozen metric-control steps. This route
does not use the first-order force or its manual clock conversion. It checks
both signed coefficients and both pointwise force fields on both spatial
rules and the radius controls. Its tolerance is the nonlinear metric
control's `1e-7`; the original first-order conversion still must pass `1e-9`.
A negative control substitutes coordinate acceleration for the proper
derivative: the independent clock gate fails and the verdict is UNRESOLVED.
The finite-epsilon metric is a variation control, not a nonlinear
Einstein-fluid solution, and this geodesic construction applies to the
initial fluid clock only.

The projection onto the fixed initial Y is a diagnostic of this local scalar
correction. It is not an operational detector record, and the preparation
does not include an order-s^3 redefinition of the initial field profile.
Changing the whole preparation or its normalization at that order is a
different comparison.

## 4. Exact coefficient and cancellation scope

At the initial instant, integration by parts gives

\[
 \int Y F_{\rm proper}(0)dV
 =-s\sum_l(30/a^2+\lambda_l/2)
                 \int\psi_l(0)(Y^2)_l\,dV.
\]

Insert the constraint solution and exact harmonic powers. The contributions
in units `kappa s^3/(V a^2)` are:

| Degree | Proper-time contribution |
|---:|---:|
| 0 | `-40` |
| 2 | `3366/125` |
| 4 | `6/5` |
| 6 | `2412/875` |
| Sum | **`-7976/875`** |

The mean is consequential: deleting degree 0 would change the result and
its sign. It is part of the prepared constraint solution, not an optional
background adjustment.

The same calculation in the fixed Newtonian coordinate time gives
`+55096/875`. The opposite signs do not describe one physical observable
with two predictions; they distinguish different time derivatives. The
proper-fluid-clock conversion is independently checked pointwise.

The result is independent of `c_s^2` within the frozen preparation because
the support has no initial density or pressure perturbation. Its initial
anisotropic stress is zero, as required by the support class. These conditions
are not a derivation of a universal BAM preparation.

With zero free tensor data, the induced homogeneous TT force vanishes at
t=0. This exact zero follows from the frozen preparation. The gate checks
that the implementation respects these initial data; its independent
numerical evidence is the comparison of the induced response over the
interval with the inherited TT solver. Thus the TT force cannot identically
cancel this nonzero scalar force. All free
tensor data were set to zero; unevolved induced tensor modes likewise have
zero initial metric perturbation and do not contribute to this initial
scalar wave-operator variation. Their later contribution is not computed.

This excludes identical cancellation for the specified data. It does not
exclude later zeros, cancellation under other free gravitational data,
different pressure/entropy preparations, or anisotropic support. The higher
odd scalar degrees sourced by the cubic force are not evolved, so the full
corrected scalar field has not been confined to the original 16 modes.

## 5. Verification and the remaining derivation

The primary run uses `a=kappa=1`, `s=0.02`, `c_s^2=1/3`, over `[0,2]`.
Controls use sound-speed squares `0,1/5,1`, two lower amplitudes, and radii
0.7 and 2. All 15 frozen gates pass:

- Generic scalar metric/connection variation yields exact symbolic zeros
  for Hamiltonian, momentum, spatial Einstein and curvature residuals;
  direct divergence-form wave-operator variation also agrees exactly.
- Full improved stress fixes the scalar anisotropic projection on two
  spatial rules. The initial coefficients agree to `2.7e-13` or better,
  including radius controls; proper-clock conversion agrees to `5.1e-16`.
- The additional same-metric geodesic check recovers both clock signs.
  Across both rules and the radius controls, its maximum coefficient error
  is `3.1e-8` and its scaled pointwise error is `1.6e-8`, within the `1e-7`
  nonlinear metric gate. Substituting the coordinate clock makes this gate
  fail even while the original first-order conversion still passes.
- Independently evolved fluid and reduced metric responses agree to
  `5.2e-12`. Maximum unused Hamiltonian and momentum residuals are
  `6.3e-13` and `6.6e-15`; tolerance refinement is below `6.6e-11`.
- The exponential-metric wave-operator check, including full ADM scalar
  curvature, agrees after Richardson extrapolation to `1.2e-10`.
  These small-metric differences are near the floating-point error floor;
  the result is not advertised as an observed step-halving convergence rate.
- Force and metric amplitude-halving ratios are eight and four. An
  independent run of the inherited TT code agrees with the induced TT
  force to `1.5e-11` after division by the cubic scale.

The finite-time demonstrations remain small by a continuous bound, not just
a sampled maximum. Decompose each oscillator forcing into a constant and
`cos(2 omega t)`, bound its hyperbolic/oscillatory response, and use the S3
addition theorem `||H_l||_infinity <= (l+1)sqrt(power_l)/V`. Across all
controls the resulting all-space/all-time potential bounds stay below
`0.00561`, comfortably inside the frozen `0.05` gate. This does not make
the homogeneous ESU mode stable at arbitrarily late times.

Missing or failed required gates produce UNRESOLVED verdicts, named failures,
overwritten JSON/Markdown reports and exit code 1. Passing a subset of the
gates is insufficient for a physical verdict.

The support now has a parameterized linear response with propagated scalar
constraints, and the initial proper-time coefficient needs no choice of
sound speed. What is missing is a BAM derivation selecting that support
class and its constitutive law, plus the remaining metric/scalar evolution
and a physical history/readout map. The causality gate remains downstream.

## Reproduce

```bash
python -m experiments.closure_ledger.esu_support_response_probe
python -m pytest -q tests/test_esu_support_response.py tests/test_scalar_tt_constraints.py tests/test_reciprocal_scalar_tt.py tests/test_tt_triangle_rotor.py
```

See the [archived report](../experiments/closure_ledger/runs/20260908_esu_support_response/probe.md)
and its [JSON data](../experiments/closure_ledger/runs/20260908_esu_support_response/probe.json).
