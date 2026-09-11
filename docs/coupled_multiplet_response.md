# Coupled degree-1 multiplet–metric tensor response

**The four-component support changes the homogeneous TT equation.** Its
field perturbations can consistently remain zero in this linear tensor
sector, but its stress changes with the metric. The resulting five tensor
components obey a periodic Hamiltonian equation rather than #292's bare
oscillator. The independent full Einstein–matter check includes the
Hamiltonian, momentum and spatial-trace equations, not only a TT projection.

The [freeze](coupled_multiplet_response_prereg.md) was published at
`64d7f9c` on #294 head `fe3d423`, before implementation, finite variations or
period-map measurements. The action, scalar-closure claim, response equation
and bare-oscillator failure were frozen analytic predictions; no prediction
needed correction. The numerical period-map outcome was not presumed.

This PR treats the first coupled sector: all five homogeneous TT metric
components around the degree-1 support, with zero scalar perturbation initial
data. It does not derive general scalar/vector responses, higher tensor
harmonics, responses around the degree-3/5 supports, or nonlinear persistence.
It is stacked on #294 while that PR remains open.

## 1. Background and the scalar closure argument

Use precisely [#294's action and support](odd_multiplet_support.md):

\[
 I=\int\sqrt{-g}\left[\frac{R-2\Lambda}{2\kappa}
  -\frac12\sum_{I=0}^3\left((\nabla\phi_I)^2+\frac R6\phi_I^2\right)\right],
\]
\[
 \phi_I=P(t)x_I,\qquad P=\sqrt{\frac{3}{4\kappa}}\cos\theta,
 \quad\theta=\frac{2t}{a}+\delta,\quad
 \Lambda=\frac{3}{2a^2},\quad\rho=\frac{3}{2\kappa a^2}.
\]

The four x_I are unit-S3 embedding coordinates. Every field is odd under
x -> -x. There is no added signal field or fluid. Work in background proper
time, lapse one and shift zero, and write the homogeneous perturbation as

\[
 g_{ab}=a^2[\exp(2\epsilon\beta)]_{ab},\qquad
 \beta=\sum_{A=1}^5 b_A E_A,\quad \operatorname{tr}\beta=0,
 \quad\operatorname{tr}(E_AE_B)=\delta_{AB}.
\]

The invariant-frame modes are spatially transverse. The degree-1 identity
\(D_aD_bx_I=-g_{ab}x_I/a^2\) implies that contraction with a TT metric
variation vanishes. Transversality removes the connection-divergence term,
and trace-freeness removes the time-volume term. The linear scalar-curvature
variation is also zero. Consequently

\[
 \delta\!\left(\Box-\frac R6\right)\phi_I=0
\]

for each of the four background fields. With zero field perturbations and
velocities initially, uniqueness for the linear Klein–Gordon initial-value
problem preserves that zero response. Arbitrary independent matter
perturbations have not been asserted absent; their zero initial data define
this experiment. This is a derived property of the tensor sector, not an
assumption that the matter stress stays rigid.

## 2. Full action: background-field terms cannot be dropped

Put \(Q=P^2\), \(C=2\pi^2a^3/\kappa\). Fixed volume gives
\(\operatorname{tr}K=0\). To quadratic tensor order,

\[
 R_3=\frac6{a^2}-\frac8{a^2}\operatorname{tr}\beta^2+O(\beta^3),
 \quad R_4=R_3+\operatorname{tr}\dot\beta^2+O(\beta^3),
\]
\[
 \sum_I|\nabla\phi_I|^2=
 \frac Q{a^2}\operatorname{tr}(e^{-2\beta}).
\]

Here the remainder counts beta and its time derivatives at the same small
order. The exact unit-volume condition eliminates the trace-extrinsic-
curvature boundary term, including when its coefficient depends on time.
The spatial curvature formula is checked against the inherited invariant
metric expression, and the independent full-curvature calculation below
tests the resulting coefficients without using this reduced action.

Keeping the nonzero background-field curvature coupling gives

\[
 L_2=\frac C2\left[f\operatorname{tr}\dot\beta^2
                -\frac g{a^2}\operatorname{tr}\beta^2\right],
 \quad f=1-\frac{\kappa Q}{6},\quad
 g=8+\frac23\kappa Q.
\]

Terms of order beta^2 times the supporting field squared contribute at the
leading tensor-response order. The older cubic truncation about a vanishing
signal scalar cannot supply them. No parameter is fitted here.

## 3. Independent stress variation and the full linear equations

Sum the complete improved stress before taking its mixed spatial trace-free
variation. With field perturbations zero, the three surviving terms are

\[
 \delta T^a{}_b\big|_{\rm TF}=
 -\frac{2Q}{a^2}\beta
 +\frac{\dot Q}{6}\dot\beta
 +\frac Q6\delta G^a{}_b\big|_{\rm TF},
\]
\[
 \delta G^a{}_b\big|_{\rm TF}
       =\ddot\beta+\frac8{a^2}\beta.
\]

The first term is the metric response of the spatial gradients. The second
comes from the connection in the Hessian of the summed field square. The
third is the nonminimal curvature term. Using mixed components removes the
background pressure times the metric variation; applying an ordinary STF
projection to covariant components would retain that spurious contribution.

Einstein's equation therefore agrees with action variation:

\[
 \boxed{f\ddot\beta+\dot f\dot\beta+\frac g{a^2}\beta=0.}
\]

The same calculation gives zero linear density, momentum, spatial-trace
and scalar-curvature responses. Their Einstein equations and all four KG
equations hold when the displayed tensor equation holds. Thus this is a
closed sector of the **linearized** Einstein–matter system, not merely the
projection of an arbitrary sourced metric history.

The independent implementation constructs coordinate metric jets through
second order and embedding-coordinate jets through third order. It then
computes Christoffel symbols, Ricci curvature and all four improved stresses.
It uses the finite metric continuation I+2 epsilon beta, whose linear tangent
matches the exponential parameterization. Finite epsilon is not claimed to
solve Einstein. Forty-two off-shell/on-equation cases cover all five STF
components, pure beta/velocity/acceleration controls, noncommuting matrices,
three spatial points, field and velocity zeros, and radius/coupling changes.
Symmetric differences at .002, .001 and .0005 have measured halving ratios
3.9948--4.0030 wherever the errors exceed the frozen noise floor. The largest
finest-step absolute tensor Einstein-coefficient error is below 9.7e-7;
full stress below 4.4e-8; scalar KG response below 1.3e-7. These are controlled
finite-difference errors, not exact-zero claims about the numerical samples.
The algebraic identities separately vanish exactly.

## 4. What replaces the bare oscillator

In dimensionless proper time tau=t/a,

\[
 f=1-\frac18\cos^2(2\tau+\delta),\qquad
 f'=\frac14\sin(4\tau+2\delta),\qquad
 g=8+\frac12\cos^2(2\tau+\delta),
\]
\[
 (fb')'+gb=0.
\]

The coefficient period is pi/2. The kinetic coefficient stays between 7/8
and 1, so no division by a field zero occurs. At tau=delta=0 the bare
initial acceleration b''=-8b, with b=1 and b'=0, leaves residual **3/2**.
The bare equation therefore fails even in a phase where the first-derivative
coefficient vanishes. A field-zero instant restores the bare stiffness
instantaneously; it does not restore the bare equation over a history.

Canonical momentum is \(\pi=Cf\dot b\). In normalized coordinates
\((b,p=fb')\),

\[
 \frac{d}{d\tau}\begin{pmatrix}b\\p\end{pmatrix}=\begin{pmatrix}0&1/f\\-g&0\end{pmatrix}\begin{pmatrix}b\\p\end{pmatrix}.
\]

The first-derivative term is not a dissipative constitutive law: the canonical
map is symplectic. The Hamiltonian depends periodically on time through the
chosen background. No claim of conserved tensor energy, damping, attraction,
or collapse prevention follows.

## 5. Measured periodic response and its limits

DOP853 integrations use the two frozen tolerance pairs and step caps
pi/100 and pi/200 (both within the frozen maximum). Independent evolution
of \(y=\sqrt f\,b\) uses

\[
 y''+\left(\frac gf-\frac{f''}{2f}+\frac{(f')^2}{4f^2}\right)y=0.
\]

At phase zero, the one-period canonical map is approximately

\[
 M=\begin{pmatrix}
 -0.0481532701 & -0.3716572780\\
  2.6844120149 & -0.0481532701
 \end{pmatrix},
 \qquad\operatorname{tr}M=-0.0963065402,\quad\det M=1.
\]

The eigenvalues are approximately -0.0481532701 +/- 0.9988399584 i.
The trace is well inside the frozen numerical ellipticity interval.
The corresponding bare-oscillator trace is -0.5325106841. These invariant
traces differ; no choice of a Floquet-frequency branch conceals that result.
Quasifrequency is deliberately not assigned, because its tau-unit value
is defined modulo 4 and up to the conjugate-branch choice.

The archive records tolerance/step refinement, symplectic residuals,
normal-form comparison, all five tensor initial components, physical-time
integrations at a=.7,1,2 and kappa=.4,1, and twenty periods compared with
powers of M. Shifting the background phase changes the period map by a
change of time origin and preserves its trace; this is also checked.
The eta=0 mathematical control reproduces the bare oscillator. It does not
claim to describe another self-consistent ESU with support removed.

The verdict is **NUMERICALLY_ELLIPTIC**, not a rigorous interval-certified
Floquet theorem and not stability of the full Einstein–matter system. Scalar,
vector and nonlinear perturbations remain outside this calculation.

## 6. Reproduction and consequence for subsequent work

```bash
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.coupled_multiplet_response_probe \
  --output-dir experiments/closure_ledger/runs/20260911_coupled_multiplet_response
pytest -q tests/test_coupled_multiplet_response.py
```

The [archive](../experiments/closure_ledger/runs/20260911_coupled_multiplet_response/)
contains the exact identities, full-geometry controls and evolution data.
All ten frozen gates are required. Missing/false gates, malformed or
nonfinite period evidence and computation exceptions invalidate the physical
verdicts, overwrite stale reports and return a failing CLI status.

For this explicit support, #292's bare tensor frequency cannot be copied
into a coupled history. The old driven orbit remains valid in its original
stated projection; this result changes the equations to be solved when that
projection is completed with the degree-1 multiplet. No replacement driven
orbit has been constructed here, and no extra signal component has been
silently introduced. General scalar/vector support response and the choice
of how signal and support fields coexist are still needed before an O(s^4)
persistence calculation. Preparation selection, Phi selection and the
source-local causality gate remain open.
