# Supported fields and the mouth mixing operation

Public freeze: `e528c847d172a4629dc1bcd737fc46979b259a2f` (#301), unchanged.
Audited baseline: `4cd86541d3b836b35561b0c4a3a54629d28851cd` (#300).
Experiment date: 2026-09-14. [Pre-registration](field_apparatus_prereg.md).

**The frozen controls pass, but the supported fields still have no specified
map to a physical mouth.** The exact quartet has no scalar-intensity angular
source. The earlier graph mixer reproduces, while the covariant free scalar
on an intrinsic circle has no mixing-induced splitting. These statements
concern different actions and are mutually consistent.

The result does not exclude an embedded mouth with transverse dynamics,
extrinsic geometry, material structure or additional boundary interactions.
It identifies what a field-to-apparatus derivation must supply, and qualifies
the older claim that changing graph hoppings already derives that geometry.

## 1. The quartet cancellation is exact, and narrowly scoped

The inherited quaternion matrices give

    B^T B=Q I, Q=q.q,
    sum_alpha phi_alpha(x)^2 = Q x.x/A^2 = Q/A^2.

The symbolic matrix residual has sixteen exact zero entries. The field-square
identity holds before restricting x to any particular circle. Thus its
nonconstant Fourier coefficients vanish on every parameterized closed curve
in S3, throughout the exact homogeneous ansatz. No numerical evolution or
special support phase is needed to establish this.

All 2,304 frozen quadrature cases confirm the identity: 64 q draws, two
scales, six coordinate great circles and three grids, with harmonics 1,2,3
checked in each case. The maximum normalized Fourier residual is
`3.2181e-16`. The component-selective control has coefficient `1/4` as
predicted: cancellation belongs to the complete quartet, not every field.

This rules out a source proportional only to the complete scalar intensity
within this ansatz. The metric is still allowed to be anisotropic, and the
full improved stress includes derivatives. A spatially homogeneous tensor
can carry directional information. We have not identified any of those
directions with the older lattice's throat winding fiber.

## 2. The graph mixer survives as a graph result

Independent incidence-matrix assembly, with wrap-around links retained,
reproduces the capstone's `fiber_H` on its N=8 lattice. The flat-measure
graph has the exact first derivative

    <+1|partial_epsilon H_graph|-1>
      = (2-sqrt(2)) exp(i varphi)         for m=2.

The first derivatives for m=1 and m=3 vanish. The 24 multipole/phase cases
each use all three frozen central-difference steps. Maximum finest-step
relative matrix error is `1.2500e-5`, below `1e-4`; all 48 consecutive error
ratios lie in `3.99984999..3.99998479`. The earlier graph number is intact.

This operator is the graph form `D^T diag(R_mid^-2) D` with an unchanged
site kinetic measure. It describes, for example, a specified variation of
link stiffness at fixed site inertia. A geometric interpretation requires
deriving both forms from the relevant field action, rather than identifying
the stiffness with a metric coefficient by notation alone.

### Review clarification, 2026-09-20: normalization and phase

For the same midpoint convention, let theta=2 pi/N. The incidence matrix
acts on a winding mode by multiplication by exp(i k theta)-1. Combining
these two factors with the m=2 Fourier coefficient of the link derivative
gives, for non-aliased N>=5,

    <+1|partial_epsilon H_graph|-1>
      = -exp(i(theta+varphi)) (1-exp(-i theta))^2
      = 4 sin^2(pi/N) exp(i varphi).

The reversed bra/ket ordering is its complex conjugate. The magnitude is
the unweighted graph's winding-one eigenvalue, lambda1=4 sin^2(pi/N).
This identifies exactly how the N=8 value generalizes; it does not alter
the original registered N=8 experiment or its archived results.

The unscaled number tends to zero as N^-2. This alone is not a physical
continuum no-go: for a circle of fixed circumference, the Laplacian includes
the inverse square of the grid spacing. Dividing by (2 pi/N)^2 makes this
coefficient tend to exp(i varphi) at radius one, and its ratio to lambda1
is already one for every N. A specified material-stiffness model can
therefore have a finite continuum coupling. That does not derive it from
the covariant geometric action of section 3.

Nor does equality of magnitudes with lambda1 identify the graph coefficient
literally with the omitted covariant lambda0 delta W: the graph uses an
R^-2 stiffness, whose first variation has a factor -2, while the covariant
stiffness uses R^-1, with factor -1. At radius one their respective
off-diagonal stiffness derivatives are exp(i varphi) and exp(i varphi)/2;
the covariant kinetic contribution exp(i varphi)/2 cancels the latter.
Both the kinetic form and the stiffness must follow from the same action.

## 3. What the covariant circle changes

For the frozen static intrinsic-circle control, a free scalar has action

    S=1/2 int dt dchi [R u_t^2 - u_chi^2/R].

The kinetic measure is `R dchi`; it is not flat. Its generalized eigenproblem
is

    -partial_chi(R^-1 partial_chi u)=lambda R u.

Set `R=R0[1+epsilon cos(m chi+varphi)]`, with abs(epsilon)<1 and m>=1.
Arclength

    s=R0[chi+epsilon(sin(m chi+varphi)-sin(varphi))/m]

has period `2 pi R0` and gives exact modes `exp(i k s/R0)` and eigenvalues
`k^2/R0^2`. The scalar equation and periodicity are satisfied exactly.
The two winding directions remain degenerate at every allowed static
amplitude, not just to first order.

The first-order calculation pinpoints the missing piece. In the fixed,
W0-normalized k=+1,-1 basis, at m=2,

    delta W_+- = exp(i varphi)/2,
    delta K_+- = exp(i varphi)/(2 R0^2),
    lambda0=1/R0^2,
    (delta K-lambda0 delta W)_+-=0.

The diagonal entries also vanish, so the entire first-order degenerate
block is zero. Keeping delta K while dropping delta W incorrectly predicts
a splitting. At R0=1, varphi=0, the stiffness-only block has eigenvalue
spread 1, while the correctly normalized block has spread zero.

All 270 frozen circle/grid cases pass. Each retains four mapped modes,
the pointwise left and right sides of the coordinate equation, both
first-variation matrices and the kinetic Gram matrix. Maximum pointwise
normalized residual is `4.3774e-16`; the finest-grid Gram error is
`7.9454e-16`; the numerical generalized-block residual is `5.8879e-16`.
The zero verdict rests on the exact identities, not these residuals.

The independent positive control is a prescribed potential on a uniform
circle: `V=v cos(2 chi+varphi)` has matrix element `v exp(i varphi)/2`.
All 144 frozen potential/grid cases agree with the Fourier selection rule
to `4.1745e-17` in the registered norm. The controls can detect physical
operator perturbations when those perturbations are supplied.

An embedded elliptic mouth is not merely an isolated intrinsic circle.
Transverse modes, extrinsic curvature, connections, material degrees of
freedom, port locations and interface conditions may retain a physical
mixing effect. None was derived by the earlier radius-to-hopping assignment.
A moving-coordinate test also requires transforming its lapse/shift and
time-dependent basis; the static result does not exclude a driven mouth.

## 4. The source audit reaches a missing interface, not a no-go

The archive records 22 inspected source files, baseline content hashes,
named function/class anchors and their hashes, and four exact baseline
search commands with their returned paths. This is a reproducible manual
source audit, not a machine proof that no extension could work.

| Candidate | What exists | What it does not supply |
|---|---|---|
| #300 supported four fields | exact 3+1 bulk action and full homogeneous constraints | a mouth, a throat-fiber map, or a surface action |
| physical/areal throat | spherical spatial gluing and constraint response | Lorentzian m=2 Einstein--scalar interface evolution |
| point, finite-tube and excised-neck wave models | conserving scalar operators and boundary matching on specified geometry | a metric shape degree of freedom sourced by the quartet |
| 5D traversable benchmark and derived network | required benchmark stress and scattering/transport | realization of that stress by the 4D quartet |
| Einstein--Israel spherical shells | radial junction equations with a chosen constitutive law | nonspherical response with the quartet's nonminimal F R coupling |
| Newtonian shell multipoles | a genuine ell=2 mutual stiffness | derived shear modulus, relativistic kinetics or throat-fiber identification |
| conditional matrix-tube parent | restoring energy from nine endpoint-constrained channels | identification of those added channels and endpoints with the four fields |
| effective throat-order field | prescribed GL potential and vortex profile | microscopic potential and its coupled Einstein realization |
| measurement graph/pointer | specified mixing and quantum apparatus statistics | classical preparation and event-frequency derivation |

The Newtonian multipole and matrix-tube entries are important positive
partial constructions. Their existence prevents an argument from saying
that every geometric coupling vanishes. Their missing constitutive and
interface data prevent promoting them to the requested coupled solution.

In particular, #300 uses `F=1/kappa-sum phi_alpha^2/6`. Adding a physical
boundary requires a compatible variational boundary term and scalar/metric
matching. The constant-coupling vacuum Israel formula cannot substitute for
that derivation. Neither matching radii nor identifying two symbols named
chi provides the missing dimension, fiber, stress and kinetic data.

No complete inherited map was located, so no coupled mouth evolution was
run and no response addendum was issued. An external shape drive or a
postulated surface modulus would define a conditional apparatus extension;
the present experiment does not choose one to manufacture a result.

## 5. Verdicts, failure checks and reproducibility

| Verdict | Result |
|---|---|
| bulk-to-mouth map | `BULK_MOUTH_MAP_UNSPECIFIED` |
| reference lattice | `REFERENCE_LATTICE_REPRODUCED` |
| intrinsic circle | `INTRINSIC_CIRCLE_MIXING_CANCELS` |
| physical mouth operator | `UNRESOLVED` |
| quartet intensity | `NO_ANGULAR_SOURCE_IN_EXACT_QUARTET` |
| field response | `BLOCKED_BY_UNSPECIFIED_MAP` |

Eight verification gates pass. The map/interface and reciprocal full-field
response gates remain false, as visibly recorded. This is not a ten-gate
affirmative field-response result. Missing physical matrix elements and
generated amplitudes are null, not zero.

Review clarification, 2026-09-20: the two unavailable map/interface and
reciprocal-response gates are explicitly hard-coded false, not evaluated
physical residuals. Their substantive justification is the separately
checked source inventory. The current verdict selector has no code path
for the freeze's licensed outcomes `INHERITED_MAP_DERIVED`,
`SCOPED_MAP_OBSTRUCTION_PROVED`, `PHYSICAL_MOUTH_MIXING_DERIVED`,
`PREDICTION_REFUTED`, `INHERITED_FIELD_RESPONSE_DERIVED`, `SOURCE_TERM_ONLY`,
`CONDITIONAL_APPARATUS_RESPONSE`, or `SCOPED_RESPONSE_OBSTRUCTION_PROVED`.
Those outcomes are unimplemented, not measured exclusions. A later
implementation would need its own tests before reporting any of them.

The separate metric/stress-channel diagnostic requested by the freeze was
not performed in this round. The source audit identifies that the improved
stress can vary with the metric, but that is not a quantitative channel
measurement. Work stopped at the missing mouth map after the exact
scalar-intensity and intrinsic-circle controls; a missing map does not
prevent examining bulk stress anisotropy separately. That required
diagnostic remains outstanding, and no zero metric/stress response is
claimed. This reporting correction does not retroactively complete it.

The probe recomputes the control checks from raw matrices and samples; saved
passing flags do not determine the verdict. Its six failure-control records
exercise missing kinetic weight, altered graph entries, nonfinite potential
data, damaged symbolic evidence, missing inventory entries, and attempted
promotion of prescribed graph success into a physical map. Unit tests also
exercise phase-case duplication, forged maps, missing gates, cache mutation,
unknown gates and stale CLI output. Nonfinite evidence withdraws its affected
result while retaining independently established results.

The final targeted suite passes **29 tests**. An independent rerun reproduces
both the compressed raw archive and Markdown report byte-for-byte. The
decompressed JSON SHA256 is
`01a41e7afd34be51fc720d10cee774bf2be947221a93fc713ea6d8868c1c9144`.

The published freeze is unchanged. No prediction, grid or threshold required
correction. Implementation validation repairs are distinct from the frozen
scientific controls. The required controls and exact identities were run
after the public freeze. The reported positive potential is explicitly an
added interaction, and no Born-state or carrier-probability calculation was
used to establish these classical results.

Run from the repository root:

```bash
python -m experiments.closure_ledger.field_apparatus_probe \
  --output-dir experiments/closure_ledger/runs/20260914_field_apparatus
python -m pytest -q tests/test_field_apparatus.py
```

The [raw archive](../experiments/closure_ledger/runs/20260914_field_apparatus/probe.json.gz)
is deterministic gzip-compressed JSON, with the [verdict report](../experiments/closure_ledger/runs/20260914_field_apparatus/probe.md)
alongside it. Exit 0 means the scoped audit and control verification passed;
it does not mean the missing field-to-mouth response was derived.

The next physical calculation needs a specified mouth/field interface and
its shape action, derived from the chosen theory or openly introduced as
apparatus structure. The present result does not select Phi, derive Born
frequencies, establish canonical commutators, or close the causality gate.
