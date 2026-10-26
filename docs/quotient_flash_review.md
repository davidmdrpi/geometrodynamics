# Second review: ray closure, Maxwell phases and antipodal rim loops

This follows [review 6402606052](https://github.com/davidmdrpi/geometrodynamics/pull/349#issuecomment-6402606052)
at `38ea0bd`. The reviewer accepts the earlier corrections about homogeneous
field freedom, capture bands and the missed first encounter. The new throat
storage condition is useful. Two stronger claims still need qualification:
an exact endpoint ordering is not reliably implemented in the cited #345
commit, and a geometrically closed ray ledger is not a matched Maxwell field.

## The interior closure band is consistent

For the stated symmetric coasting ray model, let `T=(pi-2a)/u`,
`t_f=(2k+1)pi` and `L=sin(2a)`. Capture plus retention until contact requires

\[
t_f\le T\le t_f+a/u,\qquad T<t_f+L/\sqrt{1-u^2}.
\]

Independent roots of these inequalities reproduce the review's edges:

| k | Lower speed | Upper speed | Lower-edge convention |
| ---: | ---: | ---: | --- |
| 0 | 0.7135211024 | 0.8090140683 | Included by the aperture rule |
| 1 | 0.2539419665 | 0.2696713561 | Excluded: exit coincides with contact |
| 2 | 0.1561211732 | 0.1618028137 | Excluded: exit coincides with contact |
| 3 | 0.1126622004 | 0.1155734383 | Excluded: exit coincides with contact |

The upper endpoint is included **only if capture is explicitly ordered before
contact**. These are ray-model storage bands; no spectrum or Maxwell boundary
condition follows from the timing inequalities alone. The independent test
checks transit inequalities on both sides of the lower roots. We have not
rerun the review's entire multicycle bisection scan.

## P2: the claimed contact-endpoint rule is not robust in the cited code

The audit pins #345 to `25ce28cc8e33ce7541b7d4dfac034a5574f6dcff`, reads its
actual `pair_birth.py` and `antipodal_pair_birth.py`, and verifies that their
remaining wave-module dependencies match this checkout. It does not edit #345.
Birth kicks are solved so the flash-paid speed equals the target. Runs use
the upstream tests' cutoff 12 and footprint 0.2, at dt=0.02, 0.01 and 0.005.

- At `u_k-1e-9`, all six k=1,2 runs capture both flashes, leaving no front.
- At `u_k+1e-9`, all six release first, leaving two fronts, as expected.
- At nominal `u_k`, **five of six release before capture and leave two fronts**.
  Only k=1 at dt=0.005 captures; its contact occurs about 5e-12 after focus.
  At k=1 and dt=0.02/0.01, recorded contact and focus times are exactly equal,
  yet the event order is contact, reflection, reflection.

`run()` dispatches `_continuous_event()` and release before the focal-event
block; `_release()` does not first process coincident foci. The new upstream
resonance tests use speeds slightly below the endpoint. They therefore test
interior capture, not the claimed exact priority rule. This is not a failure
of the analytic interior band: it is an unresolved numerical/boundary ordering
at its top endpoint. The fix belongs in that event dispatcher, with explicit
tie priority and exact/just-inside/just-outside tests across step sizes.
`pinned_endpoint_order.json` retains all 18 outcomes and source hashes.

## Whole-passage integrals now reproduce

Using the newly supplied endpoints, source A alone, a fixed probe at its n-th
focus, and the local increasing-birth-circle tangent gives:

| n | Integral over [n pi-pi/2, n pi+pi/2] |
| ---: | ---: |
| 1 | 0.0105520170 |
| 2 | -0.0080860608 |
| 3 | -0.0054176817 |
| 4 | 0.0131748749 |

These agree with the comment. They remain specified fixed-probe integrals,
not a derived finite mouth impulse or a universal caustic sign law.

## The actual Maxwell flash is not an exact antipodal template

For source A's modal state define `z_lm=omega_l Q_lm+i P_lm`, with
`omega_l=sqrt(l(l+1))`. Free evolution multiplies it by `exp(-i omega_l T)`.
We compare it to the antipodal spatial template, allowing the **best arbitrary
complex amplitude**. The normalized overlap is

\[
C(T)=\frac{\sum_{lm}|z_{lm}|^2(-1)^l e^{-i\omega_l T}}
{\sum_{lm}|z_{lm}|^2}.
\]

The unrepresented energy fraction `1-|C|^2` is therefore a deliberately weak
matching test: it already permits a fitted phase and gain. It is not an
absorption probability or the energy outside a physical aperture.

| T at the geometric endpoint | Unrepresented energy fraction |
| --- | ---: |
| pi | 0.00747490 |
| 3 pi | 0.06580014 |
| 5 pi | 0.17484821 |
| 7 pi | 0.32060952 |

The conformal control `omega_l=l+1/2` gives a perfect template at these times,
up to its common caustic phase. For Maxwell, frequencies sqrt(2) and sqrt(6)
already have an irrational ratio; nonzero contributions from both cannot
have a common exact finite period. This statement concerns that free modal
template, not every possible driven throat history. We retain cutoff 24, 48
and 64 calculations at band interiors and endpoints. The 48/64 projection
fractions agree within 3.4e-16; the 24/48 difference reaches 1.16e-6.

## A quotient-compatible linear field can be constructed, with changed data

The actual prescribed opposite-charge currents satisfy
`C_lm(T-t)=S_l C_lm(t)`, with **S_l=(-1)^(l+1)** for this code's transverse
Maxwell coefficients. This parity is measured from `PairBirth.source`, not
inferred by copying a scalar or spinor sign. Its largest residual in 24
controls is 3.0e-15. Positions still coast; no continuous force is applied.

For the free canonical field let `K(Q,P)=(S Q,-S P)` and U(T) be exact
oscillator evolution. A history whose specified endpoint jumps are images
must satisfy `U(T)z=Kz`. Since `A=K U(T)` is an energy-orthogonal involution,

\[
 z_+=\tfrac12(I+A)z
\]

is the nearest invariant initial state in the modal energy norm. This gives
an explicit consistent free history at **any** tested duration, and retains
homogeneous freedom. It changes the birth data; it does not show that the
original `P=0` flash already closes. The correction energy divided by original
flash energy at endpoints k=0,1,2,3 is **59.23%, 24.27%, 86.86%, 9.39%**.
The nominally exact conformal template still needs a 50% correction for this
real Q/P jump convention: a fitted common caustic phase is not equality of
prescribed real phase-space data. The projected residual near 1e-16 is a
construction check, not independent physical validation.

We also solve the **forced reference equation**
`Q''+omega^2 Q=C(t)` for these prescribed currents. Duhamel quadrature gives
the zero-initial-data endpoint `D=(Q_D(T),P_D(T))`. The time-symmetric
particular solution has initial data

\[
 z_{\rm particular}=\tfrac12 K D.
\]

For the compatible source this solves the inhomogeneous endpoint condition;
any invariant homogeneous history can still be added. At cutoff 48 its initial
energy is 0.0028254026 at the k=1 endpoint and 0.0055725818 at the band midpoint.
Orders 64,128,256,512 are retained. The coarse quadrature gives a deceptively
small symmetry residual while its initial field is wrong by about 10%:
symmetry alone is not a convergence test. The independent ODE regression and
quadrature refinement check the construction separately.
At cutoff 48 the 256-to-512 initial-state changes are 6.62e-14 and 8.26e-14
in relative energy norm. Independently evaluated interior matching defects
are at most 2.11e-13. These are accuracies of this prescribed linear reference,
not of the unsolved throat boundary or Einstein dynamics.

This advances the linear field calculation but is **not the strict BAM field**:
the source currents are continuously prescribed on the round sphere, the
initial field is selected rather than derived at birth, and neither throat
matching, finite-aperture absorption nor mouth backreaction is solved. No
momentum exchange is inferred from this external reference.

## The release pulse is not established as the flash's advanced half

The advanced flash focus coincides geometrically with previous contact under
the stated cycle timing. That is a statement about support locations and times.
The implemented release includes the pair's rest and kinetic energy and held
throat energy. The flash contains only the much smaller radiative allocation.
Moreover the release uses a conformal deformation operator and the flash uses
the Maxwell operator. Equality of their fields would need a conversion map,
normalization, polarization and phase matching. The focal coincidence does
not provide these, so the review's identification cannot yet be used as a
field equality or as the required absorber condition.

The review's clarification about momentum throughput is accepted: recycled
radiation must be counted at each encounter. Bounds for isolated propagating
packets should not be applied to a static field's momentum transfer divided
by net work. The quoted large ratio for static attraction is therefore not
an energy-conservation contradiction, nor does it install that continuous
force into the strict antipodal model.

## Both handle types remain geometrically admissible

In sphere-oriented rim frames the antipodal derivative is
`A_rim=diag(1,-1)`. The existing maps are `H_Klein=-I` and
`H_torus=diag(1,-1)`. Thus the aligned handle-plus-identification loops are

| Handle | Rim loop | Determinant | +1 eigenspace dimension in this 2D representation |
| --- | --- | ---: | ---: |
| Klein | diag(-1,+1) | -1 | 1 |
| Torus | I | +1 | 2 |

Both square to identity. The loop distinguishes orientation and allowed
parity sectors; it does not exclude either handle by itself. Selecting a
particular Maxwell polarization or twisted bundle needs the full field
matching, including normal transport, not only this two-dimensional rim map.

The ray storage band survives this review. Exact endpoint ordering and the
field/source conversion do not yet close. PR #349 remains open on those
requirements and the subsequent dynamical GR recoil calculation.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.quotient_flash_probe
# Requires the pinned upstream commit to have been fetched:
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.quotient_endpoint_review
python -m pytest -q tests/test_quotient_flash.py
python -c 'from experiments.closure_ledger.quotient_flash_probe import manifest; manifest()'
```

All diagnostics, source and tests are bound in the new run manifest. Previous
reports and full-state archives remain unchanged.
