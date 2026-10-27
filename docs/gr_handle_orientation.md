# GR handle continuation: physical time, clock offset and wave admission

This calculation continues **PR #349**, alongside the [deformable benchmark](deformable_alice_boundary.md).
It tests an actual vacuum GR geometry and its electromagnetic wave operator.
It does not evolve two mouths, supply a closed FRW embedding, or insert any
continuous force into the strict antipodal router. The earlier scalar failure
of sustained approach and collapse is unchanged.

## A reversed clock label need not reverse physical time

Let b=2M. In dimensionless Kruskal coordinates the Schwarzschild metric is

\[
 ds^2=\frac{4b^3}{r}e^{-r/b}(-dT^2+dX^2)+r^2d\Omega^2,
 \qquad X^2-T^2=(r/b-1)e^{r/b}.
\]

The antipodal map A on S2 has no fixed points. Both maps
`J_plus(T,X,Omega)=(T,-X,A Omega)` and
`J_minus(T,X,Omega)=(-T,-X,A Omega)` are free involutive isometries
away from the curvature singularity. They have different causal meanings:

| Map | Physical time orientation | Exterior Killing coordinate t=2b atanh(T/X) |
| --- | --- | --- |
| J_plus, ordinary RP3 geon | Preserved | t becomes -t |
| J_minus, time-nonorientable quotient | Reversed | t is unchanged |

The ordinary geon construction and its time orientability are given explicitly
by [Louko and Marolf, Eq. 4.12a](https://arxiv.org/abs/gr-qc/9802068).
Here the second map is checked directly against the metric; it is not a claim
that the ordinary geon already implements an advanced throat signal. On its
time-nonorientable quotient there is no global choice of future direction.
Any wave history must be compatible on the cover before descending. A freely
chosen retarded source on one exterior is not automatically such a history.

We also compute offsets. Let B(delta) be a proper Lorentz boost in the T,X
plane, and set `J(epsilon,delta)=B(delta) diag(epsilon,-1)`, with A on S2.
Every such map preserves the local metric. But

\[
 J(+,\delta)^2=I,\qquad J(-,\delta)^2=B(2\delta).
\]

For this **end-exchanging Schwarzschild family**, a nonzero boost offset on
the genuine time-reversing map fails the Z2 closure condition. It cannot be
silently used as a two-sheet quotient. A different deck group would require a
separate global analysis. This restriction is not a theorem about every
two-mouth MTY handle: two independently glued boundaries need not be the
quotient of an entire Kruskal cover. The exterior coordinate relation is
`t_destination=-epsilon*t_source+2b*delta`; this separates clock offset from
physical time orientation. Eighteen exterior/interior controls retain both
valid and invalid involutions. Their largest metric error is 3.34e-16.

## The spatial bridge does not determine a freely chosen lapse

For a static spherical metric

\[
 ds^2=-N(s)^2dt^2+ds^2+r(s)^2d\Omega^2,
\]

Einstein's equation in an orthonormal frame gives

\[
8\pi\rho=\frac{1-r'^2-2rr''}{r^2},\qquad
8\pi p_s=\frac{r'^2-1}{r^2}+\frac{2r'N'}{rN},
\]
\[
8\pi p_\perp=\frac{r''}{r}+\frac{r'N'}{rN}+\frac{N''}{N}.
\]

An independent Christoffel/Ricci construction tests these reduced expressions
on a nonconstant lapse and an asymmetric radius profile.

The repository's vacuum spatial bridge has `r'^2=1-b/r`, `r''=b/(2r^2)`.
With the Schwarzschild lapse `N=sqrt(1-b/r)` each exterior has zero stress;
the static lapse vanishes at the neck. Kruskal coordinates cover that regular
horizon, but do not turn it into a traversable static tube.

If the **same spatial slice** is assigned `N=1`, its required stress instead is

\[
8\pi\rho=0,\qquad 8\pi p_s=-b/r^3,\qquad
8\pi p_\perp=b/(2r^3).
\]

Thus zero spatial scalar curvature does not make this ultrastatic spacetime
vacuum. At any smooth flaring static throat with finite positive lapse,

\[
8\pi(\rho+p_s)=-2r''/r<0.
\]

Changing time orientation or the clock offset cannot erase that local null
energy condition violation. The Ellis profile `r=sqrt(s^2+b^2), N=1` gives
`8pi(rho,p_s,p_perp)=(-1,-1,+1)b^2/r^4`, explicitly displaying its support.
Fifteen radius/support controls include the repository relation `b=sin(a)^3`
at a=0.16, 0.3 and 0.5. These are static four-dimensional GR controls, not an
exclusion of dynamical, higher-dimensional or trapped-mouth alternatives.

## Maxwell waves on the vacuum geometry

The two radiative Maxwell parities share the radial equation

\[
 \frac{d^2\Psi}{dr_*^2}+
 \left[\omega^2-\left(1-\frac b r\right)\frac{\ell(\ell+1)}{r^2}\right]\Psi=0,
 \qquad r_*=r+b\log(r/b-1),\quad \ell\ge1.
\]

See [Lei, Wang and Jing, Eq. 6](https://arxiv.org/abs/2108.04146) for the
Maxwell reduction. We use its differential operator, not that paper's mirror
boundary condition or its quasinormal-mode calculation. Our horizon data are
purely ingoing; at the exterior end we solve for the incoming and reflected
coefficients. A radiative l=0 Maxwell mode is explicitly rejected.

| omega b | Horizon energy fraction, l=1 | Horizon energy fraction, l=2 |
| ---: | ---: | ---: |
| 0.2 | 0.00197606 | 1.88268e-7 |
| 0.5 | 0.471024 | 0.000467210 |
| 1.0 | 0.999538 | 0.785514 |
| 2.0 | Approximately 1 | 0.9999999 |

These are **horizon fluxes, not transmission to the other exterior**. Twenty-four
integrations use outer endpoints r*/b=120, 240 and 480; the maximum final
endpoint change in a horizon fraction is 1.95e-5. The largest flux defect
`|R+T-1|` is 3.20e-9. Flux precision does not remove endpoint truncation error;
unclipped values slightly above one at high frequency are retained in the data.
The horizon endpoint is r*/b=-30, where the potential is exponentially small.
The independent Wronskian drift is also recorded, including cancellation in
the small-transmission cases.

The two exterior wedges on the cover satisfy `X_A>|T_A|` and
`X_B<-|T_B|`. Consequently `|T_B-T_A|<|X_B-X_A|`, whereas a causal curve in this
metric requires `|Delta T|>=|Delta X|`. Angular motion only strengthens that
restriction. Ordinary geon identification does not add a traversing
characteristic. The time-reversing quotient requires a global field gluing
problem; it does not license relabelling ingoing horizon flux as an outgoing
wave at an independently moving second mouth. Our scattering calculation is
on the Schwarzschild cover and does not solve that quotient history.

## Consequence for the candidate attraction mechanism

An isometry of a round cross-section commutes with its Laplacian (and with
the Hodge Laplacian for Maxwell modes). Rotating, reflecting or antipodally
identifying it can mix degeneracies and parity, but cannot change its transverse
eigenvalue. A static clock reversal also preserves frequency magnitude.
These operations alone therefore do not supply the unequal channel momenta
that produced attraction in the earlier deformable benchmark. This statement
concerns spectral conversion; it does not prohibit local recoil, exchange
between mouths, interference, or ordinary gravitational attraction.

A displacement of a mouth relative to the source changes capture geometry.
A genuine nonisometric deformation can change mode matching. Neither effect
has been derived from the spherical quotient controls above. The next physical
closure must specify the dynamical geometry/support, solve the directional
field and its global time matching, and obtain impulses from that same action
at the allowed antipodal encounters. The positive defect benchmark is evidence
for a mechanism to seek, not a substitute for those equations. **PR #349 remains
open; a self-consistent BAM attraction/collapse solution is not established.**

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.gr_handle_probe
python -m pytest -q tests/test_gr_handle_checks.py
python -c 'from experiments.closure_ledger.gr_handle_probe import manifest; manifest()'
```

Source, tests and all controls are bound by the SHA-256 manifest in
`experiments/closure_ledger/runs/20261024_gr_handle/manifest.json`.
