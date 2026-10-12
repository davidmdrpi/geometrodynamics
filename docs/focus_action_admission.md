# First focus-action admission test

12 October 2026. Follow-up to #328. This is a candidate rejection/control
calculation, not the antipodal event, and not a selection of a final BAM law.

## Candidate and gates fixed before running the producer

Use signature (-++++) and a dimensionless scalar chi. The explicit candidate is

\[
 S=\frac1{2\kappa_5}\int\sqrt{-g}\,[R+(\nabla\chi)^2]d^5x,
 \quad G_{ab}=-\partial_a\chi\partial_b\chi
                 +\tfrac12 g_{ab}(\nabla\chi)^2,\quad\Box\chi=0.
\]

Include the usual metric boundary term when prescribing an induced boundary
metric. There is no scalar potential, no curvature cutoff, no corrective
kicks, and no time-dependent activation rule. Scalar boundary variation is
the normal scalar flux; it must match to an exterior solution or an explicit
surface source. A zero exterior scalar gradient is not free matching data.
This is a complete local Einstein–scalar theory on a fixed topology, not a
complete closed-bulk BAM theory. A localized activation/collar action is
absent and may not be supplied after observing a result.

The reversed scalar sign is chosen as an explicit negative control because
ordinary scalar kinetic energy cannot supply the required negative radial
null contraction. It adds a physical field, not a derivation from pure GR.
The static waist radius b is an integration scale, not a prediction.

Pre-stated admission requirements:

1. A claim to recover #328's original spatial metric and mouth data must
   satisfy **every** Einstein component and the scalar equation. A nonzero
   analytic stress identity or boundary obstruction rejects that claim.
2. If the action instead supports another regular geometry, identify it and
   its changed initial data. Do not reuse the old neck relation silently.
3. For a robust static control require no resolved negative radial eigenvalue
   on the stated two-ended background. A negative variational bound or a
   converged eigenvalue below -0.01/b² rejects static stability. This is not
   a theorem about a different finite closed-bulk boundary problem.
4. Require a bounded-below local kinetic Hamiltonian for admission as the
   final support sector. A ghost rejects this candidate under this additional
   declared viability gate, even though its classical principal equation is
   hyperbolic. Ghost sign, loss of hyperbolicity and classical exponential
   instability are different diagnoses.
5. Eigenvalue refinement: h/b=0.04, 0.02, 0.01 on [-12b,12b], compare
   [-8b,8b] and [-16b,16b] at h/b=0.02. Demand successive mesh errors
   decrease, Richardson estimates agree within 1e-5, and domain change
   is below 1e-5. Independent Gaussian Rayleigh bound must be negative.
   Failure of numerical gates is unresolved, never stability or a bounce.

The published higher-dimensional Ellis instability is known in advance.
This is an independently implemented reproduction plus BAM admission audit,
not a blind discovery. No parameter or boundary is fitted to remove a mode.

## A wider obstruction on the original neck

For any minimally coupled scalar Lagrangian P(X,chi), with
X=-(grad chi)^2/2 and a **static radial** scalar, the orthonormal stress obeys
rho+p_angular=0. This includes potentials, either kinetic sign, and sums of
static radial scalars with such stresses. It does not include angular fields,
nonminimal curvature couplings, extra anisotropic matter or evolving throats.

For a static metric -N²dt²+ds²+r²dOmega_3² the condition is

\[
 N''+2\frac{r'}rN'+\left[2\frac{1-r'^2}{r^2}-\frac{r''}r\right]N=0.
\]

On the original scalar-flat neck r²=s²+b² this reduces to

\[
 (r^2N')'+\frac{b^2}{r^2}N=0.
\]

Shell-free matching to #328's ultrastatic mouths requires N'=0 at both
ends. Integrating gives 0=integral(b²N/r²)ds, impossible for N>0.
The obstruction concerns the **total** support; adding an unspecified
second sector would change the test. In z=atan(s/b), the equation is
N_zz+N=0. Equal unit lapse at the two mouths forces N=r_m/r and therefore
N'(S)=-S/r_m², not zero. This is a direct failure of the original boundary
contract, irrespective of the scalar potential.

## The actual massless-ghost solution

Solving the displayed action, with N=1, reflection symmetry, r(0)=b,
r'(0)=0 and chi(0)=0, gives

\[
 r''=2b^4/r^5,\quad r'^2=1-b^4/r^4,\quad
 \chi'=\sqrt6 b^2/r^3.
\]

Both Einstein and scalar equations are checked, not just the NEC.
The stress is (rho,p_s,p_angular)=(-3,-3,+3)b^4/(kappa_5 r^6);
the waist NEC is -6/b² and the Kretschmann scalar is 60b^8/r^12.
This is the five-dimensional Ellis control [1], not #328's spatial slice.
In particular its initial density is negative rather than zero.
Matching just r and r' to the original round exterior at r_m=R sin(a)
would require b=R sin(a)^(3/2), instead of R sin(a)^2. Its scalar flux
r³chi'=sqrt(6)b² is nonzero at every finite mouth: turning it off outside
needs additional field/boundary dynamics, not merely a geometric collar.

## Stability calculation and interpretation

Use x=s/b, y=r/b. Solve y''=2/y^5 from y(0)=1,y'(0)=0, reflect across
zero, and discretize the regular radial master operator [1, eqs. 3.17–3.18]

\[
 H=-d^2/dx^2+W(x),\qquad W=-\frac1{4y^2}-\frac{27}{4y^6},
 \qquad H u=(b\omega)^2u.
\]

The perturbation variable allows the waist to move; it is not the singular
fixed-waist gauge. Dirichlet endpoints approximate decay at the two infinities.
We also test a Gaussian in the continuum Rayleigh quotient. A negative
omega² gives exp(sqrt(-omega²)t), not a quantum transition. The operator
is from the cited linear perturbation reduction; this PR verifies its
Darboux identity away from the waist but does not rederive the full
gauge-invariant reduction from scratch. The paper's table-I benchmark is
(b omega)^2=-2.98495893027790 for spacetime dimension five.

The scalar Hamiltonian density in a local inertial frame is
-(dot(chi)^2+|grad chi|²)/(2 kappa_5), unbounded below. Adding a potential
does not repair that high-frequency kinetic sign. Neither a ghost nor an
unstable regular solution guarantees singularity avoidance under evolution.
An unstable equilibrium can participate in a transient process, so this
rejection is of a robust final support law under the declared gates, not of
every possible transient focus mechanism or of BAM as a whole.

## Measured verdict

`CANDIDATE_REJECTED_STATIC_AND_BOUNDARY_GATES`. All numerical checks pass;
all five admission checks fail. In units of the solution's own waist radius:

| Diagnostic | Result |
| --- | --- |
| Richardson radial eigenvalue (b omega)² | -2.98495891482145 |
| Instability e-folding time / b | 0.578803063475 |
| Independent Gaussian Rayleigh quotient | -2.865902515752 |
| Original-mouth lapse derivative demanded by scalar identity | -S/r_m², nonzero |
| Local scalar kinetic Hessian | negative |

The original-neck obstruction rejects the requested old boundary contract.
The instability rejects the alternative two-ended static solution under the
declared stability gate. Its finite-mouth closure/activation law is also
absent. The result does not establish instability for every possible closed
network boundary condition or rule out a transient focus event.

## Next gate

Any successor must specify how it avoids the scalar stress/boundary
obstruction, supply bulk and boundary equations, and justify its kinetic
and characteristic structure before a nonlinear packet run. Extra angular
stress, nonminimal geometry or genuinely dynamical/trapped data change the
problem and need their own constraints; the -3/b² benchmark cannot simply
be imposed on them. A successful classical continuation would still need
an independently stated quantum-event discriminator.

## Source and reproduction

[1] Torii and Shinkai, *Wormholes in higher dimensional space-time: Exact
solutions and their linear stability analysis*, Phys. Rev. D 88, 064027
(2013), [arXiv:1309.2058](https://arxiv.org/abs/1309.2058), especially
eqs. 2.16–2.18, 3.17–3.18 and table I. Our notation b is their throat
radius and spacetime dimension is five.

    OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.focus_action_admission_probe
    python -m pytest tests/test_focus_action_admission.py -q

Results and source hashes are saved separately from #328's archives.
The action/audit tests plus #328's focused controls total 35 passing tests
under both NumPy 2.3.5 and 2.2.6. Validation includes coordinate curvature
against the action stress, scalar flux conservation, the original-neck lapse
ODE, the Darboux identity, grid/domain refinement and the Rayleigh bound.
