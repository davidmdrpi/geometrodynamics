# TT turning and complete asymptotic transport

The expanding four-component FRW support has **one simple normal-form turning
point for every support phase**, inside the frozen interval
`2.97 <= A/a <= 3.03`. Physical tensor amplitude subsequently freezes, but
its independent decaying coefficient survives in a complete real symplectic
input/output map. Neither instantaneous oscillator action stays invariant.
An exact, basis-dependent quadratic invariant still exists.

These are classical results for the homogeneous TT sector on the specified
support. They do not select Phi, derive quantization, construct a rotor,
establish nonlinear persistence, or answer the operational causality gate.

## Provenance and review

The untouched freeze is [tt_turning_transport_prereg.md](tt_turning_transport_prereg.md),
commit `2b928d1a550e24019130544e9a4ff326969df965`, based on the supported
FRW derivation `6099a0baba1a14536ebbaad50323468c3281836e`. The integration
commit already contains main `2d9e62d`; no earlier implementation is removed.
The seed, sample grids, cutoffs, series orders, tolerances and twelve gates
are those in the freeze. No frozen analytic prediction required correction.

The review's concerns about WKB at a simple zero and retaining both future
coefficients are confirmed. Its supplementary, post-freeze formulas need two
corrections. In dimensionless conformal time the bare potential is
`W0=9-R^2`, without an extra `1/a^2`. The supported constant Frobenius column
has quartic coefficient `cos(alpha)^2/16 - 28/3`, not `-8`. We derive the
recurrence below and exercise the coefficient independently in tests. The
freeze only specified the constant, quadratic and free cubic coefficients;
its predictions remain intact. The review's local-slope estimate alone does
not certify a uniform derivative bound or global uniqueness; the exact
quadratic-form argument below supplies both required global sign control
and monotonicity in the root strip.

## Normalization and phase reduction

Let `x=eta_star-eta>0`, `R=coth(x/sqrt(2))`, and `theta=alpha-2x`.
The normalized action and equation are

\[
 L=\tfrac12(m\beta'^2-k\beta^2),\qquad
 (m\beta')'+k\beta=0,\qquad p=m\beta',
\]
\[
 m=R^2-\cos^2\theta/8,\qquad k=8R^2+\cos^2\theta/2.
\]

A prime always denotes **eta**, hence `p=-m beta_x`. Physical canonical
momentum has the extra factor `Vol(S3_unit) a^2/kappa`. The normal coordinate
`y=sqrt(m) beta` obeys `y''+W y=0`, where

\[
 W=k/m-m''/(2m)+m'^2/(4m^2).
\]

Here `m>=7/8` throughout the expanding branch. Matching the departure and
support phase to hold alpha fixed reproduces #297's independently expressed
coefficients and potential. Independent proper-time integrations use
`dt=A d eta` and reproduce the same normalized canonical maps.

## Exact all-phase certificate

Put `z=R^2`, `c=cos(theta)`, `s=sin(theta)` and `N=4m^2 W`. Direct algebra gives

\[
 N=4z^2(9-z)+2z+(3z^2/4-7z+1/4)c^2
       +\sqrt2 R(z-1)cs.
\]

On `c^2+s^2=1`, this is the quadratic form with entries

\[
 H_{11}=4z^2(9-z)+3z^2/4-5z+1/4,\quad
 H_{22}=4z^2(9-z)+2z,\quad H_{12}=\sqrt2 R(z-1)/2,
\]
\[
 \det H=z^3(z-9)P(z),\qquad P(z)=16z^2-147z+12.
\]

For `1<=z<9`, `H22>0`. The convex polynomial P is negative throughout
`[1,9]`: its endpoint values are -119 and -15 and it lies below their
secant. Thus `det H>0`, and N is strictly positive for every phase.
At `R=3` it becomes exactly `(4c+3 sqrt(2) s)^2>=0`.

For `z>=zu=(303/100)^2`, P and P' are positive at zu and P' increases.
Also `H22=-2z(2z^2-18z-1)<0`, since the parenthesis and its derivative are
positive there and increase. Thus `det H>0`, H is negative definite, and
N is strictly negative for every phase throughout `R>=3.03`.

It remains to exclude repeated roots in the intervening strip. Along the
actual history, `R'=(R^2-1)/sqrt(2)` and `theta'=2` give

\[
 N'=\frac{R(z-1)}{\sqrt2}
       [24z(6-z)+3(z-2)c^2]+24zcs.
\]

For `3<=R<=3.03`, the bracket is at most
`24z(6-z)+3(z-2)<=-627`: the last expression is decreasing since its
derivative is `147-48z<0`. The positive prefactor exceeds 16, and
`cs<=1/2`. Hence

\[
 N'<-10032+12z_u<-9900.
\]

At a root, `W'=N'/(4m^2)<-9900/(4zu^2)<-20`. Existence, uniqueness and
simplicity follow for **every alpha**, without a phase scan. In fact the
proof excludes roots below R=3; the primary verdict retains the frozen
window [2.97,3.03]. R=3 itself is attainable for the phase solving
`4cos(theta)+3sqrt(2)sin(theta)=0` there.

The code verifies the polynomial identities exactly modulo `c^2+s^2-1` and
archives all rational sign bounds. A supplied `verified=true` flag or a root
scan cannot pass the certificate validator. This proof is specific to these
supported coefficients, not to arbitrary matter.

## Actions and the future basis

Writing `h=x^2 m`, `j=x^2 k`, and `beta=sum b_n x^n`, the recurrence is

\[
 2n(n-3)b_n+
 \sum_{l=1}^{n}(n-l)(n-3)h_l b_{n-l}
 +\sum_{l=0}^{n-2}j_l b_{n-l-2}=0.
\]

Since `h0=2`, `j0=16`, and `h1=j1=b1=0`, the indices are 0 and 3,
`b2=4b0`, and the resonant n=3 source vanishes identically. There is no forced
logarithm. Set the free cubic coefficient to zero in the constant column,
and to -1/6 in the decaying column. Using
`h2=2/3-cos(alpha)^2/8`, `j2=16/3+cos(alpha)^2/2` gives the corrected b4 above.
The two normalized columns satisfy

\[
 b_c=1+4x^2+O(x^4),\quad p_c=-16/x+O(x),\qquad
 b_d=-x^3/6+O(x^5),\quad p_d=1+O(x^2).
\]

Their eta Wronskian tends to 1, so it is exactly 1 for the exact solutions;
the x Wronskian has the opposite sign. Thus `beta=C b_c+D b_d` has limit C,
while D remains independent. Discarding D is a rank-one projection.
Divergence of `y~sqrt(2)C/x` does not imply growing physical beta.

For `Omega=sqrt(k/m)` and `J_beta=(p^2/(m Omega)+m Omega beta^2)/2`,

\[
 J_\beta'=\tfrac12(m\Omega)'[\beta^2-p^2/(m\Omega)^2].
\]

This is generically nonzero. Endpoint asymptotics distinguish both modes:
`J_beta~sqrt(8) C^2/x^2` for C nonzero, and
`J_beta~D^2 x^2/(4sqrt(8))` for the pure decaying solution.

At the simple zero, `W~|W'_turn|(x-x_turn)`. Therefore the WKB diagnostic
`|W'|/(2W^(3/2))` diverges as distance to the power -3/2. The normal action
`J_y` is only defined here for W>0. Generic turning data with `y'=1` give
`J_y sqrt(W)->1/2`; the tuned data `(y,y')=(1,0)` give
`J_y/sqrt(W)->1/2`, and J_y tends to zero. No universal divergence is claimed.

As a countercontrol, for any fixed positive G and exact fundamental matrix S,
`I=(S^-1 z)^T G(S^-1 z)/2` is an exact invariant. A constant oscillator also
conserves its instantaneous action. Failure of these two named instantaneous
actions is not a no-go for classical invariants or a derivation of discreteness.

## Complete input/output transport

The real supported-ESU reference is fixed to identity at eta=0. Its
monodromy is elliptic, and its bounded periodic Floquet factors are checked.
Phase shifts conjugate the ESU monodromy. The FRW coefficients approach that
reference exponentially in the past. In an ESU interaction frame the
coefficient difference is integrable, so the incoming basis has a limit.

At match x=1 the numerical input is
`B_in(L)=U_FRW(-1,-L) S_ESU(-L,0)` with L=8,12,16,20. The outgoing basis
propagates both normalized Frobenius columns from x0=.04,.02,.01 with orders
8,10,12. The complete map is `T=B_out^-1 B_in`, returning `(C,D)`.
Both exact bases have unit symplectic determinant, so their limiting map
does too. The archived finite-cutoff maps confirm convergence; this numerical
accuracy is reported for the 16 frozen phases, not as a uniform phase error bound.

All four matrix entries, raw input/output bases, absolute and relative
comparison errors, determinants and symplectic errors are archived. Seven
final comparisons cover input cutoff, output order, output start, solver
tolerance, and both alternate matches x=.8,1.2. Omitting the ESU input
reference fails to converge. Wrong output normalization or an eta/x sign
error fails the Wronskian gate. A rank-one frozen-amplitude projection fails
the complete-map validator even though it can retain a correct C row.

## Numerical results and reproduction

The 80 root samples span `3.000007884569..3.017469157194`; these extrema are
only samples. The global theorem above certifies the frozen interval.

| Check | Largest observed error |
|---|---:|
| normalized coefficient reduction | 8.97e-15 |
| normal potential, independent #297 route | 4.84e-13 |
| proper-time canonical map | 6.11e-15 |
| series residual, all 144 order/start cases | 7.89e-9 |
| final complete-map relative disagreement | 2.54e-10 |
| complete-map determinant | 2.58e-14 |
| independent action derivative | 4.32e-11 |
| exact pulled-back invariant | 2.68e-15 |
| endpoint action ratios at x=.00125 | 2.56e-5 |

All twelve frozen gates pass. The coarsest order-8/start-.04 series is the
largest residual; it is retained and passes the original 1e-8 gate. No
cutoff or series column was discarded. Tests cover the frozen dependency
table, each missing/false gate, malformed/nonfinite evidence per target,
independent second-order propagation, both asymptotic columns, rank-one
failure, and CLI exceptions overwriting stale success artifacts.

```sh
python -m experiments.closure_ledger.tt_turning_transport_probe \
  --output-dir experiments/closure_ledger/runs/20260912_tt_turning_transport
python -m pytest -q tests/test_tt_turning_transport.py \
  tests/test_frw_supported_tt.py tests/test_multiplet_scalar_stability.py
```

The targeted suite reports **153 passed** (50 new tests). The run archive is
[probe.json](../experiments/closure_ledger/runs/20260912_tt_turning_transport/probe.json),
with a compact [gate report](../experiments/closure_ledger/runs/20260912_tt_turning_transport/probe.md).
This advances classical transport on the specified expanding family. The
preparation, count-selection, and causality questions remain open.
