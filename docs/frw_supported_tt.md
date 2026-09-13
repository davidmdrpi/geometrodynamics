# Fully supported TT propagation along the exact FRW family

**The homogeneous five-component TT sector closes at linear order on the
evolving four-scalar background.** Its supported equation is

\[
 (M\beta')'+K\beta=0,\qquad
 M=\frac{A^2}{\kappa}-\frac{q^2}{6},\qquad
 K=\frac{8A^2}{\kappa}+\frac{2q^2}{3}.
\]

This is derived from the action and independently from the full improved
stress. Coordinate-curvature checks verify all four scalar equations and
all linear Einstein equations, including the constraints. Conformal and
proper clocks give the same response and finite-interval transport.

Freeze: [`c15acd4`](frw_supported_tt_prereg.md), published before implementation
or numerical measurements, on #296 at `8ac4fd3`. No frozen analytic
prediction required correction. The result concerns the inherited invariant
homogeneous TT block; it is not a derivation for every tensor harmonic.

## 1. An exact evolving background, with a separate tensor perturbation

Keep the four independent real conformal scalars and the action of #294–#296:

\[
 I=\int\sqrt{-g}\left[\frac{R-2\Lambda}{2\kappa}
 -\frac12\sum_{I=0}^3\left((\nabla\phi_I)^2+\frac{R\phi_I^2}{6}\right)\right].
\]

In conformal time,

\[
 g=A^2(-d\eta^2+\gamma),\qquad \phi_I=P x_I,\qquad P=q/A,
 \quad q=a\sqrt{\frac{3}{4\kappa}}\cos(2\eta+\delta),
 \quad \Lambda=\frac{3}{2a^2}.
\]

The supporting fields obey `q''+4q=0`. Their exact density and pressure are
`rho=3a^2/(2 kappa A^4)` and `p=rho/3`. The fixed-energy branches are

\[
 A'=\sigma\frac{A^2-a^2}{\sqrt2 a},\qquad
 w=\frac{A-a}{A+a}=\frac{d}{2+d}e^{\sigma\sqrt2\eta},
 \qquad A=a\frac{1+w}{1-w},\qquad \sigma=\pm1.
\]

Here d is the initial fractional departure from a. It is **not** the tensor
variation parameter epsilon. A is not expanded around the ESU in this round.
The exact backgrounds include departures d=0.1,0.3,-0.1 and both signs of
the branch; d=0 is a separate static control.

Use the invariant coframe and STF basis of #295, with

    g_ij=A^2 [exp(2 epsilon beta)]_ij,
    beta=sum_B b_B E_B, tr beta=0, tr(E_B E_C)=delta_BC.

The background remains an unstable support in the homogeneous scalar
direction established by #296. Choosing initially zero scalar perturbations
defines the sector here; it does not establish attraction or remove that
other admissible direction.

## 2. The scalar equations close without fixing the stress

The spatial identities on the round unit S3 are

    Hess(x_I)=-gamma x_I,
    sum_I x_I^2=1, sum_I dx_I tensor dx_I=gamma.

A symmetric invariant STF tensor is transverse: its divergence contracts
the antisymmetric invariant-frame connection with a symmetric matrix.
Consequently the linear variation of the spatial Laplacian on every x_I is
zero. Its Hessian contraction vanishes by the zero trace, and its first-
derivative term vanishes by transversality.

The time part of Box changes only through the trace of the spatial volume
derivative. That trace is exactly independent of beta because
`det exp(2 epsilon beta)=1`. Thus the changing A and P add no time-dependent
scalar source. The spatial curvature has no linear STF variation; the trace
of extrinsic curvature and the linear variation of its squared norm give
`delta R=0` as well. Therefore

\[
 \delta(\Box-R/6)\phi_I=0\qquad\text{for all four components}.
\]

The zero initial field perturbation and its first derivative remain zero by
uniqueness for the linear KG equation on the regular background interval.
This argument uses neither `1/q` nor `1/q'` and holds through field zeros.
It does not imply zero stress response: the metric, Hessian and Einstein
terms in the improved stress must still be varied.

## 3. Quadratic action with expansion and curvature coupling retained

Define `Q=P^2=q^2/A^2`, `F=1/kappa-Q/6` and `h=A'/A`. The geometric inputs are

\[
 R_3=\frac{6-8\epsilon^2\operatorname{tr}(\beta^2)}{A^2}
       +O(\epsilon^3),\qquad
 \delta^{(2)}R_4=\frac{\operatorname{tr}(\beta'^2)
                  -8\operatorname{tr}(\beta^2)}{A^2},
\]

and `sum |grad phi_I|^2=(Q/A^2) tr exp(-2 epsilon beta)`.
The curvature formula for a diagonal beta extends to all five components
by the invariant spatial rotation symmetry; independent finite variations
also use noncommuting beta and beta'.

Expansion changes an important premise from #295: the extrinsic trace is
now `tr K_extrinsic=3A'/A^2`, rather than zero. It is nevertheless independent
of beta **to all orders**. Jacobi's log-determinant identity gives this even
when beta and beta' do not commute. The fixed spatial volume, lapse, normal
and F make the associated weighted ADM boundary expression independent of
beta too. Integrating by parts with time-dependent F therefore introduces
no missing tensor term. The remaining extrinsic norm contributes precisely
`tr(beta'^2)/A^2` at quadratic order.

The full determinant is `A^8 det gamma`, independent of beta, so the Lambda
term has zero quadratic tensor variation. This has an explicit symbolic
assertion and a negative control: replacing the exponential metric by
`I+2 epsilon beta` in the action produces the spurious coefficient
`+Lambda Vol(S3_unit) A^4 tr(beta^2)/kappa`. Such a linear continuation is
used below only for the **first** variation, where its tangent is correct.

Combining gravitational curvature, nonminimal curvature coupling and scalar
gradients gives

\[
 L_2=\frac{\operatorname{Vol}(S^3_{\rm unit})}{2}
 \left[\left(\frac{A^2}{\kappa}-\frac{q^2}{6}\right)
           \operatorname{tr}(\beta'^2)
 -\left(\frac{8A^2}{\kappa}+\frac{2q^2}{3}\right)
           \operatorname{tr}(\beta^2)\right].
\]

No perfect-fluid response law or manufactured tensor source is used.

## 4. Full improved stress gives the same equation

Use mixed spatial components to avoid confusing background pressure times
the metric variation with anisotropic stress. Direct variation gives

\[
 \delta G^i{}_j\big|_{\rm TF}
   =\frac{\beta''+2h\beta'+8\beta}{A^2},
\]
\[
 \delta T^i{}_j\big|_{\rm TF}
   =-\frac{2Q}{A^2}\beta+\frac{Q'}{6A^2}\beta'
       +\frac{Q}{6}\delta G^i{}_j\big|_{\rm TF},
 \qquad Q'=\frac{2qq'}{A^2}-2hQ.
\]

The `-2hQ` term is the expansion contribution to the support response. It
is present even when the conformal oscillator has q'=0. Substitution yields

\[
 \delta(G-\kappa T)^i{}_j\big|_{\rm TF}
 =\frac{\kappa}{A^4}\left[M\beta''+M'\beta'+K\beta\right].
\]

The Hamiltonian and momentum variations vanish, and the spatial trace
variation is zero. These follow from transversality, the isotropic summed
gradient tensor, and the beta-independent volume and Q. On the displayed
tensor equation every component of the linear Einstein residual vanishes;
the scalar equations have already closed separately.

The independent implementation constructs coordinate metric jets and then
evaluates curvature, every component's off-shell improved stress and every
KG residual. It never calls the reduced tensor equation to compute curvature
or stress. Arbitrary beta, beta' and beta'' are compared off shell, while
on-equation cases check all remaining components jointly.

## 5. Two clocks and three equivalent evolution variables

Normalize `m=kappa M/a^2` and `k=kappa K/a^2`. A constant multiple of the
physical canonical pair is `(b,p=m b')`, with generator

\[
 \frac{d}{d\eta}\binom b p=
 \begin{pmatrix}0&1/m\\-k&0\end{pmatrix}\binom b p.
\]

The Hamiltonian is explicitly time-dependent. Canonical area preservation
does not establish conservation of oscillator energy or an adiabatic action.
The normal coordinate `y=sqrt(m)b` obeys

\[
 y''+\left[\frac{k}{m}-\frac{m''}{2m}
                    +\frac{m'^2}{4m^2}\right]y=0.
\]

Both endpoint transformations are required when comparing maps in y with
maps in `(b,p)`. In proper time `dt=A d eta`, the original equation is

\[
 \frac{d}{dt}\left(A M\frac{db}{dt}\right)+\frac KA b=0.
\]

The independent proper-clock calculation evolves A and eta along with the
tensor pair. The full-geometry route also changes lapse from N=A to N=1
and transforms **all** second derivatives, for example

    beta_dot=beta'/A,
    beta_ddot=beta''/A^2-beta' A'/A^3,

and identically for the field P and scale A. Comparing the finite metric
continuations in the two clocks tests covariance before taking a derivative.

At a field zero q=0, both Q and Q' vanish and the instantaneous stress
variation is zero. That does not make the normal potential bare: m'' still
contains q'^2. At that instant its difference from `8-A''/A` is
`kappa q'^2/(6A^2)`. Also q'=0 does not mean P'=0 on an evolving background.
Two extra deterministic phases, with `tan(2eta+delta)=-h/2` at eta=0.37,
explicitly test P'=0 in addition to the frozen q'=0 and field-zero cases.

## 6. Regular domain, controls and measured agreement

The primary interval is eta in [0,0.8]. All 108 combinations of the frozen
departure, branch sign, phase, radius and kappa values are integrated. A
certified bound uses monotonicity of A on each exact branch and
`q^2<=3a^2/(4 kappa)`:

    m >= min(A(0),A(0.8))^2/a^2 - 1/8.

The smallest bound in the grid is `0.3926292437`, exceeding the frozen 0.2
gate. A scale pole, nonpositive A, or nonpositive kinetic coefficient is
rejected. No claim is made about crossing such a boundary. The analytic
equation and closure proof apply wherever the stipulated background is
regular and F>0; the numerical checks cover the stated finite intervals.

There are 55 independent geometry cases: 15 individual basis/time-jet
variations, 18 noncommuting cases, 18 on-equation cases and four field/velocity
zero cases. Three generic coordinate points, both branches, field parity,
and all four individual scalar equations are checked. The independent
background Einstein residual is at most `6.22e-15` in the archive.

| Check | Maximum measured error |
|---|---:|
| Finest central variation vs analytic response, normalized | 1.35e-6 |
| Full finite-metric clock comparison, normalized | 5.24e-15 |
| Agreement among canonical tolerances, normal form and proper clock | 7.49e-15 |
| Canonical symplectic residual | 6.60e-15 |

Qualified halving ratios lie between 3.9940 and 4.0133. Both frozen DOP853
tolerances pass; composition and inverse-map checks also pass. At A=a the
coefficients and pi/2 map reproduce #295, including its trace
`-0.096306540195`. For d=1e-2,1e-3,1e-4 on a fixed interval, transport
converges monotonically to the **supported** ESU map.

The controls fail independently in the full geometry:

| Deliberately incorrect equation/input | Residual or discrepancy |
|---|---:|
| Bare FRW acceleration | 0.0178454 |
| Reversed M' term | 0.342391 |
| Missing expansion term in Q' | 0.00202796 |
| 5% wrong acceleration | 0.0339878 |
| Bare vs supported normal potential at the control point | 0.306307 |
| Spurious Lambda coefficient from linear volume, unit parameters | 2 |

All ten frozen gates pass. The archive contains the actual four maps for
every background; the validator checks their shapes, finiteness, agreement,
determinants, symplectic residuals and clock/domain evidence. Missing or
failed gates, missing identities, nonfinite/malformed evidence and exceptions
yield UNRESOLVED. The CLI replaces stale success files with a failing report.

## 7. What this supplies for the next question

This closes the supported **linear equation and sector-closure** question
posed after #296. It provides a regular classical finite-interval transport
problem on departing or approaching exact FRW histories.

These backgrounds are nonperiodic, so an arbitrary endpoint matrix is not
a Floquet period map. Its eigenvalues do not diagnose the stability of the
whole history. The bare potential `8-A''/A` drops the supporting fields and
fails even the ESU limit; it cannot be substituted to obtain a convenient
scattering problem.

The next invariant or frequency-mixing calculation now has a specified
operator. It still needs a stated classical quantity and endpoint basis.
No adiabatic invariant, preferred vacuum, particle production, general tensor
harmonic tower, director history, preparation selection or Phi is inferred.
The homogeneous scalar instability from #296, tensor-driven quadratic
scalar/vector sources, nonlinear persistence and operational causality remain
separate questions.

Reproduce:

```bash
python -m experiments.closure_ledger.frw_supported_tt_probe \
  --output-dir experiments/closure_ledger/runs/20260912_frw_supported_tt
python -m pytest -q tests/test_frw_supported_tt.py
```
