# The antipodal refocus of a trace-free deformation

> A faithful **representation** of a spin-2 field in the ℝ³ embedding of a
> fixed `S²` — not backreaction, and not linearised GR on a spacetime.

`docs/embedded_tidal_wave.md` established the projection: one potential `W`
carries both the shape `ρ = −½ΔW` and the shear
`h_ab = [2∇₍ₐ∇_b₎W]^TF`, so the drawn surface *has* the solved metric
perturbation as its own geometry. Every point of the sphere then runs its own
principal-strain history, and on a compact `S²` those histories all reconverge
at the antipode at `t = π`.

This asks what a continuous trace-free deformation does where they refocus.

## The answer: on a ring, never on the point

Because `h = sin²d · q` vanishes at both poles for **every** `q`, so does `ḣ`,
and the effective density `ρ_E ∝ ḣ_ab ḣ^ab = 2(ḣ₊² + ḣ_ˣ²)` cannot pile onto
the focal *point*. Measured on the antipode itself:

| quantity | value |
|---|---:|
| density on the antipode | `2.9e-08` of peak |
| peak amplification | `2.28×` |
| focal time | `2.99` |
| ring radius | `0.166` |
| conserved invariant drift | `1.7e-13` |

This is the same fact that forbids a spin-2 point **source**, seen at the
other end of the trip. A scalar has an `ℓ = 0` piece and can sit anywhere; a
weight-2 field has neither `ℓ = 0` nor `ℓ = 1`, so its smallest feature is a
ring of radius set by its own length scale.

### The ring radius *is* the pulse width

Not a numerical floor — the ratio is stable across a factor of four in width:

| pulse width | ring radius | ratio | amplification |
|---:|---:|---:|---:|
| 0.24 | 0.2238 | 0.933 | 2.320× |
| 0.18 | 0.1662 | 0.924 | 2.279× |
| 0.12 | 0.1139 | 0.949 | 2.238× |
| 0.09 | 0.0877 | 0.974 | 2.224× |
| 0.06 | 0.0589 | 0.982 | 2.154× |

Mean `0.952`, spread `0.058`.

## What is *not* a spin-2 effect

The amplification is finite — about `2.3×` — and it is tempting to read that
as the spin-2 structure protecting itself from a singularity. **It is not.** A
scalar pulse launched from a pole and refocused on the antipode amplifies by
the same `O(1)` factor, and neither runs away as the pulse is narrowed:

| pulse width | tensor | scalar |
|---:|---:|---:|
| 0.24 | 2.320× | 2.067× |
| 0.18 | 2.269× | 2.059× |
| 0.12 | 2.252× | 2.034× |
| 0.09 | 2.188× | 1.951× |
| 0.06 | 2.176× | 1.863× |

Launch and focus are geometrically the same situation on a sphere, so whatever
happens at one happens at the other. What belongs to the spin is the **node**
and the **ring**, not the factor. This is kept as a test
(`test_the_amplification_is_not_a_spin_two_effect`) so it cannot quietly turn
back into a claim.

### And `∫ρ_E dA` is not the conservation check

It is the *kinetic* half of the energy, oscillating against the gradient half —
its swing is `0.99` over a transit. The solver's own conserved invariant is the
check, and it drifts by `1.7e-13`.

## The material patch

A `MaterialPatch` is a geodesic disc of **labelled particles** on the reference
sphere, drawn wherever the displacement puts them. It is the only object here
that can change shape and be checked for size at the same time.

Its shape is read from the eigenvalues of its own area-weighted second moment,
in the **patch's own plane** — the least-variance eigenvector of the full 3×3
moment. That matters: at a display gain the surface tilts away from radial, so
projecting against the radial direction instead leaves part of that tilt in the
answer and rotates the measured long axis out of the surface.

- **Shape without size.** At `ε = 1e-2`, where the first-order statement is the
  whole statement, the patch's area swings by `1.9e-07` across a full return
  while its aspect ratio at display gain reaches `1.041`.
- **It reports the eigenvector, in the limit of a point.** At mid-latitude the
  patch's long axis matches the local stretch eigenvector to `1.000` at any
  size. Near the focus the shear turns over on the scale of the pulse, so a
  patch wide enough to straddle the focal ring averages a sign change:

  | patch radius | aspect ratio | alignment |
  |---:|---:|---:|
  | 0.24 | 1.1824 | 0.9591 |
  | 0.18 | 1.2889 | 0.9874 |
  | 0.12 | 1.3944 | 0.9981 |
  | 0.08 | 1.4460 | 0.9997 |
  | 0.05 | 1.4710 | 1.0000 |

  The convergence is the check that the disagreement is the patch's **size**
  and not the construction.

## Where the linear description runs out

Trace-free preserves area at *first* order. The residual is second order in
`ε` **times the local gradient of the field** — and that is not uniform. Push
to the display gain and carry two identical patches through a full return:

| quantity | mid-latitude | focal ring |
|---|---:|---:|
| worst area change | `1.93%` | `25.95%` |
| max aspect ratio | `1.0467` | `1.3974` |

Thirteen times the area error, eight and a half times the distortion, same
patch and same gain. Away from the focus the gradient scale is the wavelength;
on the focal ring it is the pulse width. The fitted residual exponent is
`ε^2.00`, so this is the linear description running out rather than a numerical
defect — and both numbers collapse as `ε → 0`.

**The area law fails first, and fails hardest, exactly where the wave
reconverges.**

## Scope

No singularity forms here and none can: this is a linear field on a fixed round
background, with no backreaction and no bulk crossing rule. What this
establishes is the geometry such a rule would have to act on — an annulus of
finite radius with a node at its centre — and the amplitude at which the linear
description stops being trustworthy. It does not supply the rule.

The gain `ε` is a display choice. The shape at any gain is the solved field.

## Reproduce

```bash
python -m experiments.closure_ledger.focal_refocus_probe
# Verdict: THE_STRAINS_REFOCUS_ON_A_RING  (7/7)

python scripts/geometrodynamics_v45_focal_patches.py                # animate
python scripts/geometrodynamics_v45_focal_patches.py --still out.png

python -m pytest tests/test_viz_focal_refocus.py -q                 # 14 passed
```

The renderer shows three things and nothing else: the continuous deformed
embedding surface coloured by the signed shear, sparse principal-strain vectors
along the eigenvectors of `h_ab` with length `∝ |λ|`, and two advected
material patches — one at mid-latitude, one on the focal ring. The camera
follows the wave, so the refocus happens facing the viewer.
