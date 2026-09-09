# Existing scalar stress supplies a leading-order TT rotor

**An explicit real degree-3 scalar supplies the normal stress of a uniformly
rotating uniaxial homogeneous TT projection at order s^2.** Its Hamiltonian
dipoles and all ten momentum compatibility charges vanish. The construction
uses chosen scalar phases and nonzero free tensor initial data. It establishes
existence at the stated perturbative order, not selection of that preparation
or a complete Einstein-matter history.

The [pre-registration](driven_normal_stress_prereg.md) was published at
[`11e625f`](https://github.com/davidmdrpi/geometrodynamics/commit/11e625fa877de2bb0669af06384347ede69d5777)
on main `c08f46a` before implementation or measurements. The scalar pair,
constant source, negative tensor amplitude and exceptional speed were all
frozen analytic predictions; their verification is not a numerical search
discovery. No analytic prediction required correction.

## 1. Which equations this solves

The [free rotor obstruction](tt_triangle_rotor.md) concerns the embedding
`beta=A(nn^T-I/3)` inside the five-component tensor equation. A rotating free
trajectory generally leaves that cone. A source must supply

\[
 \mathcal N_n(S)=2A\left(\dot n\dot n^\mathsf T
                       -\frac{|\dot n|^2}{2}P\right),\qquad
 P=I-nn^\mathsf T,
\]

as well as satisfy the radial and angular equations. The normal projector
is `N_n(M)=PMP-tr(PMP)P/2`.

The [scalar-TT action](reciprocal_scalar_tt.md) already fixes

\[
 S(q)=\frac1C\sum_\mu(q^\mathsf TF_\mu q)E_\mu,\qquad
 C=\frac{V}{\kappa},\quad V=2\pi^2a^3,\quad
 \ddot\beta+\frac8{a^2}\beta=S(q).
\]

The factor 1/C converts an integrated stress/action source to an acceleration
source. Both its sign and its normalization are independently checked against
the full improved pointwise stress, including nonzero scalar momenta.

Write `phi=s phi_1+O(s^3)`, `beta=O(s^2)`. At order s, the scalar obeys the
free round-background equation. Its quadratic stress determines the order-s^2
metric response. The scalar feedback from TT, scalar and vector metric
sectors enters at order s^3, changing stress at order s^4. The
[ESU response](esu_support_response.md) showed why those sectors must be
included in a later completion. It does not prevent this leading test.
Other metric harmonics are allowed to be sourced at order s^2; the linear
metric equation does not mix them into the homogeneous TT projection.

## 2. A constant quadrupole from a time-dependent scalar

Use unit embedding coordinates on S3 and define Haar-normalized harmonics

\[
 f=\sqrt8\,(x_0^3-3x_0x_3^2),\qquad
 g=\sqrt8\,(x_1^3-3x_1x_2^2),
\]
\[
 \phi(t,x)=\frac{s}{\sqrt V}
        \left[f(x)\cos(4t/a)+g(x)\sin(4t/a)\right].
\]

These are two spatial modes of one real scalar. No time-dependent actuator
holds their phase: the displayed history solves the leading free equation.
Their relative phase and equal amplitudes are prepared initial data.

Exact monomial moments, independent of the numerical harmonic basis, give

\[
 \langle f^2\rangle=\langle g^2\rangle=1,\quad \langle fg\rangle=0,
 \quad\langle D_i f D_j f\rangle
       =\langle D_i g D_j g\rangle=\operatorname{diag}(3,3,9),
\]
\[
 \langle D_{(i}fD_{j)}g\rangle=0.
\]

Here D is the dimensionless inherited invariant derivative and the brackets
use normalized Haar. The trace-free integrated gradient tensor is therefore
constant: the two diagonal contributions add with `cos^2+sin^2=1`, and their
symmetric cross contribution vanishes. The inherited improved stress yields
the same tensor:

\[
 S=kQ_z,\qquad Q_z=e_ze_z^\mathsf T-I/3,\qquad
 k=\frac{6s^2}{Ca^2}.
\]

This fixes the source from scalar data, rather than defining the source by
differentiating a desired tensor orbit.

## 3. The exact leading orbit and the frequency exception

Take

\[
 \Omega=\frac{\sqrt2}{a},\qquad
 A=-\frac{3s^2}{2C},\qquad
 n=(\cos\Omega t,\sin\Omega t,0),\qquad \beta=AQ_n.
\]

The amplitude is negative, on the oblate branch of the existing signed
uniaxial cone. Exact substitution gives `beta_ddot+8 beta/a^2=k Q_z`.
The radial, angular and normal equations all hold. The normal force is
nonzero; the free obstruction is supplied by the actual scalar stress.

The mechanism is transparent in the Fourier decomposition:

\[
 \beta=-\frac A2 Q_z+\frac A2
 \begin{pmatrix}
 \cos2\Omega t&\sin2\Omega t&0\\
 \sin2\Omega t&-\cos2\Omega t&0\\
 0&0&0
 \end{pmatrix}.
\]

The first term is the static particular solution, since
`k=-A omega_T^2/2`. The oscillatory term is a free tensor oscillation with
`2 Omega=omega_T`. Selecting its two quadratures gives the circular orbit.
The source thus supplies a static offset; it does not dynamically select the
free oscillation or its phase.

For general degree-3 scalar data, the source has frequencies 0 and 8/a.
A nonzero constant-A uniform rotor requires source frequencies 0 and
2 Omega, with oscillatory coefficient `omega_T^2-4 Omega^2`. Frequency
independence on an open interval therefore requires either `Omega=4/a`
or that this coefficient vanish. The latter is the realized constant-source
exception `Omega=sqrt(2)/a`. The other allowed frequency is not classified
by this round. The original speed 0.4/a cannot work in this particular
uniform, constant-A, single-multiplet leading-order class; variable amplitude
and nonuniform motion remain outside that frequency statement.

## 4. Complete leading constraint compatibility

The order-s^2 scalar energy density is even, so its four Hamiltonian dipoles
vanish. Momentum compatibility requires an independent argument. The code
integrates `j_i=-T_0i` against all six ambient rotations and all four gradient
dipoles, and an exact polynomial certificate verifies their vanishing for
arbitrary scalar phase. It does not equate the three invariant generators
with all six SO(4) rotations.

On the prepared CMC slice these are the kernel conditions for the linear
constraint inverses. The scalar operator `Delta+3/a^2` has only the degree-1
kernel. For the longitudinal momentum operator `div L`, integration by parts
gives `int W.div(LW)=-int |LW|^2/2`, identifying its kernel as the conformal
Killing fields. Orthogonality to that complete kernel supplies the Fredholm
compatibility condition. The remaining scalar and longitudinal initial
responses are required; they have not been set to zero. This round certifies
compatibility, rather than assembling or evolving the full initial metric
and extrinsic curvature. Later charge checks are checks of the leading free
source, not a substituted test of constraint propagation.

Two negative controls protect the distinction. One has all three invariant
charges below `2.3e-16`, yet an omitted rotation charge of `0.441`. The other
has zero Hamiltonian and rotation charges to roundoff, yet a gradient charge
of `0.357`. Both are rejected by the complete compatibility gate. This uses
the inherited stress directly and does not depend on #291's implementation.

## 5. Evidence and assumption budget

All ten gates pass. The primary run has `s=0.02`, `a=kappa=1`,
`A=-3.0396355093e-5` and `k=1.2158542037e-4`, over one tensor period
`pi/Omega=2.2214414691`, at 401 times.

| Check | Maximum error |
|---|---:|
| Action versus full improved stress, normalized by k | `7.7e-15` |
| Hamiltonian and all ten momentum charges, divided by s^2 | `2.8e-15` |
| Full tensor equation, normalized by k | `2.2e-15` |
| Unrestricted tensor trajectory, divided by abs(A) | `4.4e-12` |
| Distance to the full signed uniaxial cone, divided by abs(A) | `3.9e-12` |
| ODE tolerance refinement, divided by abs(A) | `4.4e-10` |
| Director projector agreement | `6.0e-14` |

The ODE evolves all five tensor components and all sixteen free scalar
components with their momenta, using the inherited action source. It never
projects back onto the cone. The isolated eigenline advances by pi/2 at half
the tensor period, modulo n~-n, confirming director motion in the actual
tensor. Radius controls preserve the source and tensor normalization; both
scale by four under scalar amplitude halving. The largest tensor Frobenius
norm across the frozen controls is below `2.9e-5`. This is a bound on the
tested tensor projection, not on every unevolved metric component.

Removing the source, reversing its sign, doubling the director speed or
changing the second scalar quadrature amplitude by ten percent produces
normalized full-equation residuals `0.816`, `1.633`, `4.243`, `0.171`.
Any failed or missing required gate makes every physical verdict UNRESOLVED;
the CLI recomputes the verdict and overwrites stale success reports.

| Ingredient | Status |
|---|---|
| Scalar and homogeneous TT modes; interaction strength | Inherited action |
| Perfect-fluid background support with no anisotropic stress | Assumed completion class |
| Odd degree-3 scalar sector, mode pair, relative phase and amplitude | Chosen preparation |
| Constant quadrupolar stress | Derived from that scalar |
| Signed A and Omega for the proposed orbit | Fixed by compatibility, conditional on that preparation |
| Nonzero free tensor oscillation and its phase | Chosen initial data |
| Complete leading constraint compatibility | Exact certificate and independent integrals |
| Full Einstein-matter evolution and higher-order persistence | Not established |
| Autonomous rotor ensemble, triangle map and Phi selection | Not derived |

The fixed source selects an axis and the matched orbit lies in its transverse
plane. This is not an autonomous round rotor over the entire director space.
It also differs from #290's zero-tensor preparation, so that round's initial
no-cancellation coefficient cannot be transferred to it. At higher order,
the scalar/support/vector response and other tensor modes can modify the
stress and the orbit. Those corrections must be included in any persistence
claim.

Finally, the scalar and director frequencies have irrational ratio
`(4/a)/(sqrt(2)/a)=2 sqrt(2)`. A periodic tensor does not close the full joint
history. No identification with the triangle-family holonomy, canonical
ensemble, operational readout or counting function follows from this result.

## Reproduce

```bash
python -m experiments.closure_ledger.driven_normal_stress_probe
python -m pytest -q tests/test_driven_normal_stress.py tests/test_reciprocal_scalar_tt.py tests/test_scalar_tt_constraints.py tests/test_tt_triangle_rotor.py
```

See the [archived report](../experiments/closure_ledger/runs/20260908_driven_normal_stress/probe.md)
and its [JSON certificates and trajectory](../experiments/closure_ledger/runs/20260908_driven_normal_stress/probe.json).
