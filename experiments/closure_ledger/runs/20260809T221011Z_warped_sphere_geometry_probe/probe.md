# The restored geometry: one continuous S², warped by the wave it carries

_2026-08-09T22:10:11+00:00_

> A **display** of a solved classical field as a radial displacement of a fixed closed surface. Not backreaction.

## The dolls

| quantity | value |
|---|---:|
| `R_inner` | 0.7400 |
| `R_mid` | 1.0000 |
| `R_outer` | 1.2600 |
| half-gap `Δ` | 0.2600 |
| display gain `g` | 1.60 |
| pulse width | 0.18 |
| mesh | 61 × 91 |

## One closed surface

| check | value |
|---|---:|
| seam mismatch | `2.4e-16` |
| pole spread | `0.0e+00` |
| carries the poles | True |

## Nested, never touching

| quantity | value |
|---|---:|
| `r_min` over the run | 0.7594 |
| `r_max` over the run | 1.2226 |
| clearance, inner | 0.0194 |
| clearance, outer | 0.0374 |

## The focus is measured

| quantity | value |
|---|---:|
| deepest warp at distance | 3.141593 |
| distance error from π | `0.0e+00` |
| at time | 2.9814 |
| against π | 3.1416 |
| shortfall | 0.1602 |

The shortfall is the pulse's own width, not solver error:

| pulse width | peak time | shortfall |
|---:|---:|---:|
| 0.24 | 2.9452 | 0.1963 |
| 0.18 | 2.9814 | 0.1602 |
| 0.12 | 3.0363 | 0.1052 |
| 0.08 | 3.0725 | 0.0691 |

## And it rings

| | fraction of the gap | at time |
|---|---:|---:|
| driven out toward `R_outer` | 85.6% | 2.988 |
| pulled in toward `R_inner` | 79.6% | 3.211 |

## Verdict

**6/6 checks passed.**

**THE_SURFACE_ITSELF_CARRIES_THE_WAVE.** THE SURFACE ITSELF CARRIES THE WAVE. The archive's geometry is back, and this time the warp is solved rather than drawn.

ONE CLOSED SURFACE. The mesh carries its own poles and its φ seam matches to 2.4e-16, so it is a single manifold with nothing cut out of it. That is what makes 'a pulse sweeps every point once and fills its own void' a statement about a closed surface rather than about a patch — and it is exactly why a ring, not a pulse, is what a throat needs.

NESTED, NEVER TOUCHING. Over a full return the radius stays within [0.7594, 1.2226], clearing the inner doll by 0.0194 and the outer by 0.0374. The bound is structural: tanh cannot leave the gap.

AN HONEST DISPLAY. The displacement preserves the sign of the field everywhere and the ordering of every pair of amplitudes (worst inversion 0.0e+00), and no gain can flip a sign. It does NOT preserve ratios — the distortion spans 0.005 across the surface — which is the price of keeping the picture inside the vacuole, and is stated rather than hidden.

THE FOCUS IS MEASURED, NOT KEYED TO THE CLOCK. The deepest deformation over the run sits at geodesic distance 3.141593 from the source — the antipode, to 0.0e+00 — at t = 2.9814 against π = 3.1416. The shortfall of 0.1602 is the pulse's own width, not solver error: narrowing the pulse from 0.24 to 0.08 shrinks it monotonically. This is the difference from the archive, which grew its mound on a growth function tied to the frame number.

AND IT RINGS. The arrival is not one mound. The surface is driven 85.6% of the way to the outer doll at t = 2.988, then inverts and is pulled 79.6% of the way to the inner doll at t = 3.211. The focus ends up pulling the geometry INWARD — toward the inner shell, which is the one the ring caustic lands on.

SCOPE. This is a display of a solved field as a radial displacement of a FIXED surface, not backreaction: the wave does not feel what it is deforming, so no throat forms here and nothing about nucleation follows. What it restores is the object the intuition was built on — and it now moves because a solver says so.