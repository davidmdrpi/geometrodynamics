# Null-horizon continuation: preflight of the submitted protocol

Result: **SUBMITTED_PROTOCOL_REQUIRES_CORRECTION**. This is an equation-level
audit and finite-energy horizon-trace quadrature, **not a finite-packet spacetime
evolution**. It does not establish the proposed physical trilemma as written.
The original [submitted protocol](protocols/null_horizon_continuation_submitted.txt)
is preserved byte-for-byte. No corrections below are represented as frozen
before the analytic exploration. No original #128–#333 evidence is overwritten.

The local metric here is the 5D vacuum Tangherlini solution, independent of the
canonical-supported compact background and EGB surface action of #333. Changing
cosmological coordinates cannot remove that model's physical instability, but
that particular support action is not an admission requirement for this local
calculation. Neither its topology nor a closed cosmological embedding is solved.

## 1. Coordinate and operator corrections

Write b=r_s, kappa=1/b, u=t-r*, v=t+r*, U=-b exp(-u/b), V=b exp(v/b).
The submitted product UV is correct. The metric follows directly from
`-f du dv = -F dU dV`:

\[
 F=-\frac{f}{\kappa^2 UV}
   =\frac{(r+b)^2}{r^2}e^{-2r/b},\qquad F(b)=4e^{-2}.
\]

F is dimensionless for these dimensionful U,V. The first submitted expression
for F is four times this outside the horizon; its second expression is 4b²
times this, so the two submitted expressions also disagree when b is not set
numerically to one. An absolute value of UV cannot be used for continuation
of that expression into the interior. The regular metric has g_UV=-F/2.

For Phi=r^(-3/2) u_l Y_l, direct reduction of Box Phi=0 gives

\[
 V_l=f\left[\frac{l(l+2)}{r^2}+\frac{3f}{4r^2}
              +\frac{3f'}{2r}\right],\qquad
 \partial_U\partial_V u_l+Q_l u_l=0,\quad Q_l=\frac{F V_l}{4f}.
\]

The expression for Q is evaluated by cancelling f, not dividing numerically
at r=b. The submitted potential was already superseded by the repository's
canonical operator in #271 (`tangherlini.dynamics.master_potential`). Both that
old potential and the extra f in the submitted Kruskal wave equation matter.

\[
 Q_l(b)=\frac{e^{-2}}{b^2}\,[l(l+2)+3]\ne0.
\]

Thus vanishing *tortoise* potential is not a free Kruskal wave equation on the
horizon. One horizon trace is permissible characteristic data, but by itself
does not determine a field on an adjacent two-dimensional domain. A second
intersecting characteristic segment and compatible corner data are needed.

## 2. The three frozen laws

| Law | Actual preflight finding | Status of full packet gates |
|---|---|---|
| A: J(U,V,Omega)=(-U,-V,antipode) | Fixed-point-free local isometric quotient, reversing time orientation. A future-directed lift across right H+ goes into the black-hole interior. J sends right H+ to the *left* past-horizon branch, not the outgoing right branch specified in B. | Fails the requested retarded horizon-return requirement. No return evolution run. |
| B: U_out=-alpha V_in | Negative generator derivative reverses ordering and has incompatible transverse-side orientation for the two retained exteriors. A consistently normalized formal jump is zero, **not** the preregistered negative-density jump. | Fails before a physical null-shell tensor can be assigned to that exterior-to-exterior junction. |
| C: r=b+epsilon | Odd Dirichlet is well-defined. The even vector d_V-d_U is not the fixed-r normal except at t=0 in this chart. Correcting the vector gives a timelike reflecting test-field boundary whose null-ray return delay diverges as epsilon tends to zero. | No gravitational wall action or stress is specified. #333's action-specific obstruction cannot automatically be imported. |

For A, this is a statement about the local lift and causal return requirement,
not a claim that all exterior backscattering vanishes. The quotient need not
admit arbitrary retarded exterior data; its global causal initial-value problem
is not established. Horizon regularity also does not regularize r=0.

For C the proper static normal derivative on the rescaled field u is proportional
to `U d_U + V d_V`: d_(r*)=kappa(U d_U+V d_V). A physical Neumann condition on
Phi itself would instead induce a Robin condition on u. We retain the submitted
choice of u when describing the corrected *test-field* wall.

At r_obs=2b, the earliest radial-null round trip is

\[
 \Delta t_{\min}=2[r_*(2b)-r_*(b+\epsilon)].
\]

This is a geometric lower bound on the reflected signal arrival, not a simulated
packet peak. It grows as -b log(epsilon/b)+O(b). Observer proper time multiplies
it by sqrt(f(2b)). Numerical rows are in the result archive. Reflection at fixed
epsilon and observation at arbitrarily late time must not be confused with
finite-time return in the epsilon=0 limit. A Gaussian has no sharp wavefront;
this bound refers to a localized disturbance or compactly supported data.

## 3. Compute the common null frame before the jump

Let V parameterize both faces, with e_- = d_V and e_+ = g'(V) d_U.
At the horizon, r_(UV)=-2/(b e²). With e dot N=-1:

\[
 N_-={2\over F_h}\partial_U,\qquad
 N_+={2\over F_h g'}\partial_V,\qquad
 (Nr)_-=-V/b,\quad (Nr)_+=-{g\over bg'}.
\]

Consequently

\[
 [Nr]={V-g/g'\over b},\quad
 (8\pi G_5)\mu=-{3\over b}[Nr],\quad
 (8\pi G_5)p=-{g''\over g'}.
\]

These use C_ab=-N_mu e_a^nu nabla_nu e_b^mu, sides minus to plus, with
inaffinity C_VV=g''/g' on the plus face. The pressure sign in the submission is
opposite to this convention. The angular relation is C_ij=b(Nr) gamma_ij
=b² theta_N gamma_ij/3 in five dimensions, not b theta_N gamma_ij/2.
Spherical symmetry makes the surface current zero.

For the submitted g=-alpha V, both Nr values equal -V/b and g''=0. But this
does **not** produce a stress-free bounce. To see the problem without a sign
convention argument, use incoming x=U_in<0 and put on the other side

\[
 U_{out}=g(V),\qquad V_{out}=h(V)x+O(x^2),\qquad x>0.
\]

Keeping the outgoing exterior requires h>0. Continuity of the cross metric
requires g'h=1. Thus g'>0 is necessary. The negative-slope law has no continuous,
nondegenerate Lorentzian gluing of these retained sides in a regular chart.
Choosing h<0 would continue into a different Kruskal region. Merely comparing
induced degenerate metrics or forcing a sign flip on N misses this condition.

For every **orientation-compatible** spherical map with V>0, g<0, g'>0,
[Nr]>0 and mu<0. This is a local, restricted result for joining two equal-radius
vacuum Tangherlini exteriors with this angular identification. It is not a
theorem excluding regular finite-thickness focus laws or modified equations.

As a separately labelled exploratory control, choose g=-c/V, c>0. Then

\[
 g'>0,\quad [Nr]=2V/b,\quad
 (8\pi G_5)\mu=-6V/b^2,\quad (8\pi G_5)p=2/V.
\]

This map preserves the horizon Killing generator and has a boost-invariant
parameter c. It passes local orientation and formal field-trace conservation,
but requires a nonzero negative-density distributional shell. It is not an
admitted BAM focus under the submitted smoothness/stress-free criterion, and
no claim of global chronology, stability, or backreacted packet evolution follows.

## 4. Coordinate covariance is not a dynamical gate as submitted

Under a passive boost V'=lambda V, U'=U/lambda, the same physical linear map
has alpha'=alpha/lambda². Reusing the same numeric alpha changes the law.
For the reciprocal control, c'=c. The coordinate t reconstructed from V/U
changes by b log(lambda); a fixed observer's clock remains the same only after
accounting for this origin shift. The submitted assertion that this coordinate
t remains unchanged is incorrect.

Likewise, integral T_VV dV scales as 1/lambda. It is not an observer-normalized
energy. With unit angular harmonic normalization and a(V)=u(0,V), the stationary
Killing energies of the traces are

\[
 E_{in}=\kappa\int_0^\infty V(a')^2dV,\qquad
 E_{out}=\kappa\int_{-\infty}^0(-U)(\partial_U a_{out})^2dU.
\]

Both transforms g=-alpha V and g=-c/V preserve these positive trace energies.
Only the reciprocal map preserves the *oriented* symplectic pairing with both
generators directed to the future: the negative-slope map reverses its sign.
Parity (-1)^l cancels in these quadratic quantities for l=0 and l=1.

The probe directly integrates the specified Gaussian and a second independent
variation q=(V-2b)a/b, using the transformed coordinates for each boost. This
avoids a vacuous zero Wronskian for a single real profile. It records raw flux,
Killing energy, and symplectic orientation separately. The Gaussian is finite
energy but not compactly supported, contrary to the protocol's wording.
These are **trace identities**, not G3 passes or an observer-waveform test.

The metric is fixed; a nonzero Misner-Sharp mass response has not been solved.
A coupled calculation must specify the gravitational equations, boundary data,
and the surface matter response. Merely equating trace energies cannot supply it.

## 5. Reproduction and next admitted calculation

```
python -m experiments.throat_dynamics.null_horizon_continuation_probe
pytest -q tests/test_null_horizon_continuation.py
```

The archive records source hashes, environment, all trace quadratures, numerical
checks, law-specific gate dispositions, and an explicit list of uncomputed
quantities. Tests independently reduce the scalar operator, pull back the metric,
check the gluing orientation, and exercise nonunit b and coordinate boosts.

A corrected evolution protocol must specify a second characteristic segment,
a clock normalization, a surface/finite-thickness action if using a genuine
return map, and a compact packet if sharp support is required. The reciprocal
map gives a precise stress target for that next focus-law calculation. Repeating
stationary parity tests or passing raw coordinate-flux checks cannot replace it.

References: [Poisson, null-shell formalism](https://arxiv.org/abs/gr-qc/0207101),
especially equations (11)–(16), for conventions and normalization; the angular
trace above is explicitly extended from two to three angular dimensions.
[Manzano and Mars, general matching](https://arxiv.org/abs/2309.14874) provides
the broader hypersurface-data framework. The retained-side compatibility
condition used here is derived directly above.
