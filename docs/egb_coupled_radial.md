# Coupled metric–scalar radial problem and physical kinetic norm

**Result: `COUPLED_RADIAL_ADMISSION_FAILED`.** For the precise canonical
bulk/surface action of #331, the linear spherical problem is now reduced,
including the metric response, both shell junctions, scalar current and total
charge. The physical kinetic form is computed from that action, with its
Gauss–Bonnet boundary term. The oscillatory shell mode near #330's frequency
has **negative kinetic norm**. There are also growing modes in both reflection
parities. The earlier positive prescribed-exterior stability window does not
provide a positive-energy, stable radial mode of this material completion.

**Review correction of attribution:** the shell-free canonical-supported S4
already has one conformal negative kinetic direction and an Eddington growing
mode. The shell adds one negative direction, the negative-energy oscillation,
and two growing modes at the benchmark. Focus admission must be measured
relative to that background. See the [no-shell and small-mouth controls](egb_coupled_background_review.md);
the original numerical label and results below are preserved.

At R=1, a=.3, alpha=-.06, eta=.1:

| Parity / mode | omega^2 R^2 | Physical kinetic coefficient | Consequence |
|---|---:|---:|---|
| Even, lowest | -6.37353785051 | -623.939808035 | Growing, negative norm |
| Even, next | -2.48065462824 | +14.0419695103 | Growing, positive norm |
| Even, shell oscillation | +1.96815015124 | **-0.199921359973** | Real frequency, negative norm |
| Even, next oscillation | +6.42283192043 | +37.2757908273 | Positive norm |
| Odd, lowest physical | -1.41649415649 | +5.54078664666 | Growing, positive norm |

Even modes use fractional areal shell displacement z=delta r_shell/r0=1.
Odd modes have z=0 and unit integral of sin(u)^3 psi^2 over the half-bulk.
The displayed coefficients omit the common positive Omega3/kappa5 and scale
as R^4. Their magnitudes therefore depend on the declared normalization;
their signs do not. The first two even growth rates are 2.52458667/R and
1.57500941/R; the odd growth rate is 1.19016560/R.

This completes the **linear spherical** constraint and kinetic calculation for
#331's nominated action. It does not compute nonlinear evolution, nonspherical
stability, smooth finite thickness, or a quantum event. It is not a no-go for
all regular focus laws. No bulk/surface material parameter has been refitted.

## Background, gauge and retained degrees of freedom

Retain the same supported ultrastatic S4 interval u in [a,pi-a], identifying
its ends as one shell with two incident faces. R is the curvature radius,
r0=R sin(a), C=1+4 alpha/R^2>0, q^2=3C/R^2, and chi0=q t. The bulk action is
P_b=X-V (lambda=0); the surface action is exactly the power law plus constant
of [#331](egb_shell_kinetic_admission.md). Its trace scalar is continuous and
its linear EOS is delta p=eta delta s, with W=s0+p0>0 and 0<eta<=1.
Stresses in this document are multiplied by kappa5.

Use a scalar longitudinal gauge in the smooth bulk:

\[
ds^2=-(1+2n)dt^2+R^2(1+2\psi)(du^2+\sin^2u\,d\Omega_3^2),
\qquad \chi=q(t+v).
\]

The lower shell is at u=a+xi(t); its invariant fractional areal displacement is
z=psi(a)+cot(a) xi. A moving coordinate boundary is necessary in this gauge.
For the reflection-even sector the upper face is its mirror. Lapse and shift
constraints are imposed before assigning a physical norm; setting the shift
to zero does not discard its equation.

For reference, before the conformal spatial gauge, write the radial spatial
metric perturbation as ell and the angular one as psi. The linear momentum
and density equations are

\[
3C\,\partial_t[\cot u(\ell-\psi)-\psi']=q^2v',
\]
\[
H=v_t-n,
\quad H=-\psi''-4\cot u\,\psi'+\cot u\,\ell'
 +2(\cot^2u-1)\ell-2\csc^2u\,\psi.
\]

The second equation follows from (C/2) delta R4=delta rho and
q^2=3C/R^2. In the conformal gauge ell=psi it becomes
H=-(Delta_u+4)psi, where Delta_u=partial_u^2+3cot(u)partial_u.
The spatial traceless equation is
(2psi+Cn)''-cot(u)(2psi+Cn)'=0. Removing the residual clock and dipole
coordinate functions, the remaining equations give

\[
 n=-2\psi/C,\qquad v=-R^2\dot\psi,
\qquad R^2\ddot\psi=(\Delta_u+4+2/C)\psi.                 \tag{1}
\]

An additive constant in v is the constant scalar shift, quotiented at fixed
linear charge. The bulk scalar equation and remaining spatial metric equation
then follow. These equations use the canonical bulk law; changing P_b,XX
changes the reduction and is not covered by this calculation.

Independent tests reconstruct the full five-dimensional orthonormal Riemann
and contract G_AB+alpha H_AB, rather than inserting (1) back into the same
reduced equation. Density and all four spatial equations match the scalar
stress perturbations. The bulk equation has real local characteristic speed;
the instabilities below are global low modes, not an elliptic bulk principal
symbol.

## Coupled boundary equations

Write a mode as psi=f(u)Q(t), with Qddot=-omega^2 Q and Lambda=omega^2 R^2.
Define

\[
e=2/C,\quad h=\Lambda+e,\quad x=\sin^2a,\quad
D_0=\frac{2C}{\cos a}+\frac{16\alpha\cos a}{r_0^2},
\quad d=\frac4{\cos a\,D_0},
\]
\[
b=\frac{RW}{6C\eta},\qquad j=\frac{RD_0}{2Cr_0},\qquad
w(u)=\sin^3u.
\]

The bulk boundary-value equation is

\[
-(wf')'=(\Lambda+4+e)wf.                              \tag{2}
\]

At the shell, H=h f. The scalar material gives

\[
\delta s=(W/\eta)h f(a),\qquad\delta p=W h f(a).
\]

The geometrical variations, with inward normal on each retained face, are

\[
\delta Q_n=\frac1R[f'-\cot a\,f-\csc^2a\,\xi],\qquad
\delta k_0=-\frac{\Lambda}{R}\xi-\frac{2f'}{CR},
\]
\[
\delta s=-6C\delta Q_n+48\alpha Q_{n0}z/r_0^2,
\quad\delta p=4\delta Q_n+2C\delta k_0
             -16\alpha Q_{n0}\Lambda z/R^2,
\quad Q_{n0}=\cot a/R.
\]

Using W=2C tan(a)/R-D0/r0, both Davis junctions reduce to

\[
f'(a)+\tan a\,f(a)+b h f(a)-jz=0,
\qquad (x\Lambda+d)z-h f(a)=0.                         \tag{3}
\]

The full Davis tensor contraction independently reproduces these conditions.
The integrated scalar-current condition is

\[
\delta s+3Wz=\frac{6C}{R}[-f'(a)+\xi],                 \tag{4}
\]

and agrees with the density junction. It has not been imposed as a separate
external energy reservoir. Perturbed bulk scalar gradients carry the exchange.
For even modes, f'(pi/2)=0. Root finding uses the **pole-free** determinant

\[
(x\Lambda+d)[f'(a)+\tan a f(a)+b h f(a)]-j h f(a)=0.    \tag{5}
\]

Dividing by x Lambda+d before root finding would create a false singularity.
The exceptional value h=0 is treated separately, not removed simply because
its squared frequency is negative.

For odd physical modes, continuity of the induced lapse and scalar trace on
the identified shell requires H(a)=H(pi-a). For h != 0 this gives f(a)=0;
odd parity gives f(pi/2)=0 and z=0. The first such mode is growing despite its
positive norm. The exceptional odd solution f=cos(u), h=0, z=0 is generated
by a radial coordinate displacement proportional to sin(u), accompanied by
the clock displacement that preserves zero shift. Its kinetic form is exactly
zero. It is gauge, and is not the odd instability reported in the table.

## Total charge and the zero mode

For a nonzero-frequency even mode, the change of the bulk shift charge on a
constant-t slice includes the displacement of both excised caps:

\[
\delta Q_b=2qR^4[(\Lambda+4+e)\int_a^{\pi/2}wf\,du-w(a)\xi].
\]

The single shell contributes

\[
\delta Q_\Sigma=\frac{r_0^3}{q}W[3z+h f(a)/\eta]
=2qR^4w(a)[-f'(a)+\xi].
\]

Integrating (2) with f'(pi/2)=0 gives delta Q_b+delta Q_Sigma=0.
Odd physical modes cancel between the two halves. Thus the admitted modes
remain in the same **linearized** charge sector, unlike #331's frozen trial.
The constant scalar shift carries no perturbation of stress and is removed.
A variation of the background action constants or a change of the conserved
charge is not counted as a radial eigenmode. Exact fixed-charge nonlinear
continuations would additionally have to satisfy second-order constraints;
no such continuation is claimed here.

## The physical kinetic form from the action

The auxiliary spectral mass is not the answer. Compute the original action's
symplectic form in a gauge with the shell fixed and the time slices orthogonal
to it. For general ell,psi the quadratic gravitational time-derivative terms
in the smooth bulk are

\[
S^{(2)}_{g,\mathrm{kin}}=\frac{\Omega_3}{2\kappa_5}
\int dt\int_a^{\pi-a}du\,R^4w C
 [-6\dot\ell\dot\psi-6\dot\psi^2]
+\frac{\Omega_3}{\kappa_5}\int dt\,
 24\alpha r_0^3Q_{n0}\dot z^2.                         \tag{6}
\]

**The last term is essential.** It is the combined bulk integration-by-parts
and two-face Myers contribution. It is negative for alpha<0. No time-boundary
term is silently dropped to change the symplectic form: equivalent endpoint
terms change the potential by an exact variation, while (6) retains the
remaining spatial boundary kinetic term.

One way to check its coefficient independently is the exact warped-product
reduction, for base metric -dt^2+L^2 du^2 and sphere radius r(t,u):

\[
\mathcal R_5=\mathcal R_2+6B-6\Box r/r,
\quad\mathcal L_{GB}=12\mathcal R_2 B
 +24[(\Box r)^2-(\nabla\nabla r)^2]/r^2-24B\Box r/r,
\quad B=[1-(\nabla r)^2]/r^2.
\]

On each inward face the Myers integrand is
-r^3[k0+3Qn+2alpha(6k0(B_Sigma-Qn^2)-2Qn^3+12A Qn+6B_Sigma Qn)].
For k0=0 its quadratic time terms, after time integration by parts, are
36alpha r^3Qn psidot^2+24alpha r^3/R psidot psidot'
-24alpha r^3Qn elldot psidot. The lower-end boundary from the bulk GB
reduction cancels the derivative and mixed terms, leaving
12alpha r^3Qn psidot^2 **per face**. The source check also evaluates the full
curvature action and both Myers terms on periodic trial fields before expansion;
its error decreases quadratically with amplitude, reaching 9.74e-8 relative.

Include the scalar momenta as well:

\[
\delta p_\chi=q\sqrt{h_0}(\ell+3\psi+v_t-n)/\kappa_5,
\quad\delta p_\Sigma=q r_0^3[3P_Y z+Z_\Sigma(v_t-n)]/\kappa_5.
\]

The bulk metric momenta follow from (6). To compare gauges, use a radial
coordinate displacement epsilon=-xi E(u), with E(a)=1, E(pi/2)=0, and clock
shift T=-R^2 xidot J(u), J'=E, J(pi/2)=0. The transformed amplitudes are
ell=f+xi E', psi=f+xi cot(u)E, and v=R^2[-f+xi J]Qdot.
Two choices E=[(pi/2-u)/(pi/2-a)] and its cube give identical total symplectic
norms. The bulk, shell-scalar and gravitational-boundary terms separately
change under this transformation; retaining all three is necessary.

After the constraints and (3), the physical kinetic bilinear on even profiles
and their shell displacements is

\[
\boxed{\quad
\mathcal G[(f,z),(g,\zeta)]
=6CR^4\left[\int_a^{\pi/2}w(f'g'-4fg)du
 -w(a)\tan a\,f(a)g(a)\right]
 +3D_0r_0^4 z\zeta .\quad}                           \tag{7}
\]

For one eigenmode its reduced quadratic action is

\[
S^{(2)}=\frac{\Omega_3}{2\kappa_5}\int dt\,
\mathcal G[(f,z),(f,z)](\dot Q^2-\omega^2 Q^2).
\]

This is a physical time-evolution symplectic coefficient, not a normal-to-shell
momentum. The corresponding energy is the quadratic generator for the
background time translation combined with its scalar shift (the H-qJ
convention of #331). It is not an ADM energy at spatial infinity; this universe
has no such boundary. The linear norm calculation is distinct from nonlinear
fixed-charge integrability or claims of quantum vacuum decay.

Equation (7) establishes an indefinite kinetic form throughout the WEC/C>0,
D0<0 branch of this canonical completion. The shell contribution is negative.
The bulk block has exactly one negative direction: the Sturm–Liouville operator
-Delta_u-4 with f'(a)=-tan(a)f(a), f'(pi/2)=0. The positive function cos(u)
is the ground eigenfunction at eigenvalue 4 when the equator endpoint is
Dirichlet instead. Separated-boundary eigenvalue interlacing places exactly one
Neumann-equator eigenvalue below 4 and none equal to it. Thus the even physical
kinetic form has **two negative directions**, not just one numerical negative
value. One belongs to the supporting closed-universe background; only the
shell block is the additional negative direction attributable to this focus
construction. The background-relative distinction is essential. The odd physical Dirichlet/Dirichlet block is positive; its instability
is a potential instability instead. This count concerns kinetic directions,
not the number of growing eigenmodes.

## A self-adjoint spectral formulation and complete linear action

For a reproducible matrix formulation introduce y=xz-f(a), set
m=4+e and mu=-w(a)j/(d-xe)<0. For U=(f,y), V=(g,v), define the auxiliary forms

\[
M(U,V)=\int wfg\,du+w(a)b f(a)g(a)+\mu yv,
\]
\[
K(U,V)=\int w(f'g'-mfg)du
+w(a)[-\tan a-be+j/x]f(a)g(a)
+\frac{w(a)j}{x}[f(a)v+yg(a)]-\mu(d/x)yv.
\]

Equations (2)–(3) are K U=Lambda M U. In particular M has one negative
boundary degree of freedom, and **M alone is not the physical norm**.
The physical forms in these variables are

\[
G=6CR^4(K+eM),\qquad
V_{\rm phys}=6CR^2(KM^{-1}K+eK).                       \tag{8}
\]

Here products denote the associated operators on their constrained domain,
or their conforming finite-element matrices. On the even sector G has no gauge
null direction. The quadratic action (Qdot^T G Qdot-Q^T V_phys Q)/2 gives the
same evolution R^2 M Qddot+K Q=0. Substitution y=xz-f(a) reduces G to (7),
which independently agrees with the original symplectic calculation (6).
This supplies the action normalization that cannot be inferred from the
spectral pencil by itself.

For h != 0, x Lambda+d != 0, eliminating z gives f'(a)=F(Lambda)f(a), with

\[
F=-\tan a-bh+\frac{jh}{x\Lambda+d}.
\]

A third norm formula is

\[
\mathcal G=6CR^4h\left[\int wf^2du-w(a)F'(\Lambda)f(a)^2\right],
\quad F'=-b+\frac{j(d-xe)}{(x\Lambda+d)^2}.             \tag{9}
\]

All three norm routes agree. Omitting the boundary derivative in (9), or
replacing G by M, would produce an incorrect energy-sign test.

## Numerical validation and limits

The [source-bound result](../experiments/closure_ledger/runs/20261013_egb_coupled_radial/result.json)
records eight even low modes, the odd roots in a specified scan interval, charge
and junction residuals, norm decompositions, and all validation gates. It is a
new calculation; the original #330/#331 protocols, producers and archives are
unchanged. The [validation protocol](egb_coupled_radial_protocol.md) explicitly
discloses prior exploratory results and is not blind preregistration.

Validation: **97 focused tests pass** with NumPy 2.3.5 and 2.2.6, SciPy
1.17 and Python 3.12.14. Nine new tests cover bulk density/spatial/momentum
tensor contractions, the Davis junction, charge, three norm routes, independent
finite elements, scaling, exceptional/gauge modes, the full action Hessian and
source provenance. These are separately implemented checks in this repository,
not a claim of an external independent review.

DOP853 shooting is checked at two tolerances. A separately assembled linear
finite-element pencil uses 80,160,320 intervals and five-point element
quadrature. The first four eigenvalue errors fall by about four per doubling;
the finest errors are 1.57e-5, 8.38e-6, 1.17e-6 and 3.57e-4. Each matrix has
two negative physical kinetic directions and one negative auxiliary mass
direction. No complex eigenvalues appear in these finite pencils. This finite
calculation is not a certified count of the entire infinite spectrum.

The shell oscillation has omega^2 R^2=1.96815015124 compared with #330's
prescribed-exterior value 1.94540699079. The small frequency shift does not
repair its negative kinetic norm. At unit shell displacement its direct
Myers kinetic contribution is -0.240283095; cancellations with bulk and
surface matter leave the reported -0.199921360. A frequency-only test would
have missed the decisive failure.

The general coupled equations and physical form are given above; the computed
low modes already falsify radial stability and positive-energy admission at
the nominated point. Equation (7) extends the indefinite-kinetic conclusion
across the stated branch without claiming the same eigenfrequencies throughout
it. A different material law, geometric background, or regularized shell needs
its own constraint reduction. Adding surface fields solely to rescue the old
frequency would not address this action-level failure automatically.

Reproduction:

```bash
python -m experiments.closure_ledger.egb_coupled_radial_probe
python -m pytest tests/test_egb_coupled_radial.py -q
```

Primary conventions: [Davis](https://arxiv.org/abs/hep-th/0208205) for the
EGB/Myers action and junction; [Julie and Berti](https://arxiv.org/abs/2004.00003),
constant GB-coupling limit, for the ADM constraint/action structure. Our scalar
normalization is that of #331. Neither source supplies this matter-supported
compact-bulk spectrum; equations (1)–(9) and its checks are derived here.
