# One real conformal scalar cannot supply exact odd-sector ESU support

**A smooth, nonzero real scalar with pointwise isotropic improved spatial
stress on a complete round S3 cannot be antipodally odd.** Consequently the
existing single conformal scalar cannot, within that imposed odd sector,
replace the unspecified supporting fluid of the exact ESU.

There is an exact homogeneous scalar support if the odd condition is
relaxed. Its gravitational and scalar kinetic coefficients stay positive.
But its generic scalar perturbations have anisotropic stress, so its
radiation-like background does not derive the perfect-fluid response used
in [the ESU response round](esu_support_response.md).

The [freeze](scalar_esu_support_prereg.md) was published at `de55f3f` before
implementation or numerical checks. These are its predicted outcomes;
no analytic prediction required correction. The global exclusion below is
an analytic argument, not an inference from failed harmonic searches.

## The question is conditional, and local support is not a history measure

Use precisely one real scalar and the inherited 3+1 action

\[
 I=\int\sqrt{-g}\left[\frac{R-2\Lambda}{2\kappa}
   -\frac12(\nabla\phi)^2-\frac1{12}R\phi^2\right]d^4x,
 \qquad ds^2=-dt^2+a^2d\Omega_3^2.
\]

The improved stress is the one in `waves/two_wave.py` and
`waves/backreaction.py`, including its \(G_{\mu\nu}\phi^2/6\) term.
No new matter, potential or statistical average is inserted. The odd
condition \(\phi(t,-x)=-\phi(t,x)\) is imposed, not established by geometry
or constraint solvability; [round 11](parity_solvability.md) remains intact.

The theorem covers all smooth real configurations on the sphere, including
arbitrary harmonic mixtures, under exact pointwise isotropy. It does not
exclude complex or multiple fields, extra stresses which cancel anisotropy,
singular configurations, averaged backgrounds, non-round geometries or
different actions. This 3+1 ESU question is distinct from the 5D throat
support audited in `bulk/source_audit.py` and `shells/junction.py`.

Neither the negative theorem nor the even control selects a preparation
law, a triangle-to-field map, a history weight \(\Phi\), or a physical
early readout. The causality question remains open.

## Pointwise isotropy

On any fixed round spatial slice, let \(D\) be its covariant derivative and
TF denote the spatial trace-free part. The metric is static, so the spatial
Hessian of \(\phi^2\) has no time-connection contribution. Since the spatial
Einstein tensor is isotropic,

\[
 T_{ij}^{\rm TF}
 =\frac23(D_i\phi D_j\phi)^{\rm TF}
  -\frac13\phi(D_iD_j\phi)^{\rm TF}.
\]

This identity does not invoke the wave equation. On a connected nonzero
component \(U\), write \(h=1/\phi\). Direct differentiation gives

\[
 T_{ij}^{\rm TF}=\frac{\phi^3}{3}(D_iD_jh)^{\rm TF}.
\]

Thus isotropy implies \(D_iD_jh=\lambda g_{ij}\) on \(U\). Taking the
divergence and commuting derivatives on the radius-a sphere yields

\[
 D_j\lambda=D_j\Delta h+R_j{}^kD_kh
 =3D_j\lambda+\frac2{a^2}D_jh,
 \qquad D\lambda=-\frac1{a^2}Dh.
\]

Hence \(A=h+a^2\lambda\) is constant on the connected component. Put
\(f=h-A\); then \(D_iD_jf=-fg_{ij}/a^2\). This fixes the local solution as

\[
 h=A+B\cdot x,
\]

where x is the unit ambient coordinate and B is a constant vector in R4.
One can see that B is constant, rather than simply quote the harmonic
spectrum: on the unit sphere define
\(B=fx+\operatorname{grad}_{S^3}f\). Its ambient derivative in a tangent
direction v is
\((vf)x+fv-fv-(vf)x=0\), by the Hessian equation and the sphere's second
fundamental form. Thus A and B are the same throughout U, even if U is
not the whole sphere.

## The global step

Division by \(\phi\) has only been used on U. To justify extending the
conclusion across a putative zero set, suppose U has a boundary point p
in the sphere. For any sequence \(x_n\in U\) approaching p,

\[
 |h(x_n)|\le |A|+|B|,
 \qquad |\phi(x_n)|\ge\frac1{|A|+|B|}>0.
\]

The denominator is finite and nonzero since h exists on U. Continuity
therefore makes \(\phi(p)\ne0\). A connected neighborhood of p is then
nonzero and meets U, so p belongs to the same component, contradicting
that it is a boundary point. U is both open and closed. Connectedness
of S3 makes U all of S3.

Smoothness also forbids \(A+B\cdot x=0\); on the complete sphere its range
is \([A-|B|,A+|B|]\), so \(|A|>|B|\). Every nontrivial smooth isotropic
slice is therefore nowhere zero and has one sign. A continuous real odd
function necessarily has a zero, by following a path from x to -x.
An odd isotropic slice must be identically zero.

Apply this at **every** time: a smooth odd isotropic history is the zero
history, including its time derivative. An isolated identically zero
slice of an even oscillating field is not an exception to this all-time
argument. The zero history cannot supply
\(\rho+p=2/(\kappa a^2)\), regardless of Lambda. This finishes the
odd-sector exclusion.

The reciprocal form alone is a necessary spatial classification. The
initial implementation left its time, momentum and density equations open.
A [separate post-freeze extension](scalar_esu_uniqueness.md) now proves
that the homogeneous control below is the full smooth single-real-scalar
support family, even without imposing odd parity. It checks the rank
argument and its zero-vector exception, and supplies an independent route
through the zero-momentum Einstein equation. The freeze and original probe
archive remain unchanged.

## An exact even control

For a spatially constant scalar, the wave equation and stress become

\[
 \ddot\phi_0+\phi_0/a^2=0,\qquad
 \rho=\frac12(\dot\phi_0^2+\phi_0^2/a^2),\qquad p=\rho/3.
\]

The Einstein equations require
\(\kappa\rho=3/a^2-\Lambda\) and
\(\kappa p=\Lambda-1/a^2\). They are solved by

\[
 \phi_0=\sqrt{3/\kappa}\cos(t/a+\delta),\quad
 \Lambda=\frac3{2a^2},\quad
 \rho=\frac3{2\kappa a^2},\quad p=\frac1{2\kappa a^2}.
\]

All momentum components vanish. These are pointwise, all-time identities,
including the times \(\phi_0=0\) and \(\dot\phi_0=0\). Lambda and amplitude
are required background relations; the phase is free initial data. The
control is even and is excluded from the primary odd sector.

Write the gravitational term as \(F(\phi)R/2\), where
\(F=1/\kappa-\phi^2/6\). The regular conformal change
\(g_E=(\kappa F)g\) gives the scalar kinetic coefficient

\[
 K_E=\frac1{\kappa F}+\frac3{2\kappa}\left(\frac{F_\phi}{F}\right)^2
     =\frac1{(1-\kappa\phi^2/6)^2}.
\]

For the control, \(F\ge1/(2\kappa)\) and \(1\le K_E\le4\).
There is no gravitational-sign crossing or scalar ghost in this check.
This establishes neither nonlinear stability nor stability of every
perturbation multiplet.

## Its response is not the assumed fluid

Use longitudinal gauge and physical normal-frame density and pressure:

\[
 ds^2=-(1+2\alpha)dt^2+a^2(1-2\psi)d\Omega_3^2,
 \quad \phi=P(t)+\chi(t,x),\quad P=\phi_0.
\]

Let \(j_i=-T_{ni}=D_iJ\) and
\(\delta T_{ij}^{\rm TF}=(D_iD_j\Pi)^{\rm TF}\). Direct variation of the
action's stress, retaining the metric dependence, gives

\[
\begin{split}
 \delta\rho={}&\dot P\dot\chi+P\chi/a^2-P\Delta\chi/3
 -\alpha\dot P^2+P^2\Delta\psi/3+P^2\psi/a^2-P\dot P\dot\psi,\\
 J={}&(P\dot\chi-2\dot P\chi-P^2\dot\psi-P\dot P\alpha)/3,\\
 \delta p={}&\dot P\dot\chi/3-P\ddot\chi/3-\alpha\dot P^2/3
 -2\alpha P^2/(3a^2)+P\dot P\dot\alpha/3+2P\dot P\dot\psi/3\\
 &+2P\Delta\chi/9+P^2\ddot\psi/3-P^2\psi/(3a^2)
 -P^2\Delta(\psi-\alpha)/9,\\
 \Pi={}&-P\chi/3+P^2(\psi-\alpha)/6.
\end{split}
\]

These expressions are off shell for the perturbations. The perturbed
curvature and scalar equation are

\[
 \delta R=-6\ddot\psi+12\psi/a^2+4\Delta\psi-2\Delta\alpha,
\]
\[
 \mathcal K=-\ddot\chi+\Delta\chi-\chi/a^2-2\alpha P/a^2
   +(\dot\alpha+3\dot\psi)\dot P-P\delta R/6=0.
\]

The off-shell trace identity is
\(-\delta\rho+3\delta p=P\mathcal K\). On shell the pressure is
radiation-like, but the anisotropic stress generally is not zero. The
trace-free Einstein equation for degrees l>=2 requires

\[
 (1-\kappa P^2/6)(\psi-\alpha)=-\kappa P\chi/3.
\]

Thus generic physical scalar perturbations have a nonzero potential
difference and anisotropic stress. Setting Pi to zero for these same
histories violates the trace-free Einstein equation. No claim is made
that every special perturbation has nonzero Pi.

For evolution, solve this equation algebraically for alpha. Evolve
\(\mathcal K=0\) and the spatial Einstein equation

\[
 2\ddot\psi-2\psi/a^2-\tfrac23\Delta(\psi-\alpha)=\kappa\delta p.
\]

The coefficient matrix for \((\ddot\chi,\ddot\psi)\) is
\(\begin{pmatrix}-1&P\\\kappa P/3&2(1-\kappa P^2/6)\end{pmatrix}\),
with determinant **-2**. It never divides by P or its derivative.
The initial Hamiltonian and momentum constraints are

\[
 2(\Delta+3/a^2)\psi=\kappa\delta\rho,
 \qquad -2\dot\psi=\kappa J.
\]

Their coefficient determinant when solving for \((\psi,\dot\psi)\),
multiplied by \(a^2\), is
\((L-3)+(2L-12)u+(L-3)u^2\), with \(L=l(l+2)\) and
\(u=\sin^2(t/a+\delta)\). It is strictly positive for l>=2.
The checks evolve the scalar and spatial equations without enforcing
these two constraints again, so their later residuals test propagation.
The determinant observations strengthen the frozen turning-point gate;
they do not change the predictions.

All responses are coefficients of a small perturbation. The archive uses
unit response coefficients and quotes a physical amplitude of 1e-4;
it does not treat order-unity metric coefficients as a nonlinear solution.
Only degrees 2 and 3, a full scalar period and the stated radius control
are evolved. Homogeneous/degree-1 stability and nonlinear continuation are
not established.

## Why this does not fill the existing support slot

| Ingredient | Status |
|---|---|
| One existing real conformal scalar and its improved stress | Inherited |
| Antipodally odd field condition | Imposed, not derived |
| Exact odd-sector ESU support by that scalar alone | Excluded by the global isotropy argument |
| Homogeneous even support | Exact control, outside the primary sector |
| Positive kinetic coefficients for the control | Derived |
| Zero-anisotropic-stress fluid response | Not supplied by the generic control response |
| Lambda/amplitude relation | Required for this control |
| Phase and perturbation data | Chosen |
| BAM support selection, triangle map, history measure, physical readout | Not derived |

There is also a change of perturbative hierarchy: for support and signal
in the same field, \(T[P+\epsilon\chi]\) contains a term linear in epsilon.
The response explicitly has a nonzero such term. The zero-background
\(\phi=O(s)\), stress=O(s^2) hierarchy of #289--#292 therefore cannot be
imported unchanged into this scalar-supported control. A separate support
species would be another assumption and is outside this round.

## Reproduce and inspect

All nine frozen gates pass. Independent checks give:

| Check | Result |
|---|---|
| Isotropy identity versus inherited improved stress | 1.4e-16 normalized |
| Nonconstant reciprocal spatial-isotropy controls | 1.4e-17 absolute |
| Homogeneous full Einstein residual, all radius/coupling controls | 1.9e-15 normalized |
| Stress response versus direct metric/Ricci central differences | 7.9e-8 normalized at the coarser step |
| Observed error ratio when halving that step | 3.998 to 4.009 |
| Fine evolution Hamiltonian/momentum residual, all controls | 2.7e-11 absolute |
| Coarse versus fine response trajectories | 1.7e-9 absolute |

Every failed or missing gate is tested through the CLI, including
overwriting stale successful verdicts. The zero history, reversed Lambda,
and a 10% amplitude change fail the required background equations.

```bash
python -m experiments.closure_ledger.scalar_esu_support_probe --output-dir experiments/closure_ledger/runs/20260909_scalar_esu_support
python -m pytest -q tests/test_scalar_esu_support.py tests/test_esu_support_response.py
```

The [archived report](../experiments/closure_ledger/runs/20260909_scalar_esu_support/probe.md)
and [JSON](../experiments/closure_ledger/runs/20260909_scalar_esu_support/probe.json)
record exact local identities, inherited stress checks, direct metric/Ricci
variation, turning-point coverage, constraint propagation, and refinement.
The written global proof requires mathematical review; checking that its
file exists is not formal verification. Every mandatory gate must pass,
including the independent geometric response and the missing/failed-gate
controls. The CLI recomputes verdicts before overwriting reports, and exits
nonzero on failure.
