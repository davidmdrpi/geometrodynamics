# Projecting `h_ab` into a continuous embedded surface

_2026-08-10T03:32:14+00:00_

> A faithful **representation** of a spin-2 field in the ℝ³ embedding of a fixed `S²` — not backreaction.

## The construction

```
X = (R + ερ) n̂ + ε∇W      ρ = −½ΔW      h_ab = [2∇₍ₐ∇_b₎W]^TF
```

## Why a height field cannot do it

| `d` | trace | trace-free part |
|---:|---:|---:|
| 0.60 | +1.4498 | `+1.7e-04` |
| 1.10 | -2.3537 | `+1.3e-04` |
| 1.90 | -3.1636 | `+7.5e-05` |
| 2.50 | +1.1350 | `+1.8e-04` |

A radial deformation is **conformal** at first order.

## The theorem: the drawn surface has the solved metric

| `d` | `h₊` measured | `h₊` solved | trace |
|---:|---:|---:|---:|
| 0.40 | +0.000024 | +0.000023 | `-6.4e-07` |
| 0.90 | +0.000802 | +0.000800 | `+3.2e-06` |
| 1.50 | +0.002622 | +0.002625 | `-6.2e-06` |
| 2.20 | +0.000000 | +0.000000 | `-6.5e-07` |
| 2.80 | +0.000000 | +0.000000 | `-6.7e-07` |

Worst relative error `4.7e-04`, trace `1.1e-03`, off-diagonal `6.7e-08`, at gain `1e-04`.

## The quadrupole is the Legendre shape

| quantity | value |
|---|---:|
| `ρ` against `P₂(cos d)` | `1.0e-07` |
| amplitude ratio | 1.000000 |
| residual dipole | `-4.2e-17` |

## Area, at first order

| quantity | value |
|---|---:|
| drawn area | 12.56631057 |
| round area | 12.56637061 |
| relative change | `4.8e-06` |

## Reach into the bulk

| quantity | value |
|---|---:|
| toward `R_outer` | 74.8% |
| at time | 3.085 |
| toward `R_inner` | 44.3% |
| display gain | 7.87 |

## Verdict

**7/7 checks passed.**

**ONE_POTENTIAL_CARRIES_SHAPE_AND_SHEAR.** ONE POTENTIAL CARRIES SHAPE AND SHEAR. A spin-2 field can be drawn as a continuous deformation of the embedded sphere, and the projection is forced rather than chosen.

A HEIGHT FIELD ALONE CANNOT. For X = r n̂ the induced metric is g_ab = r²ĝ_ab + ∂_a r ∂_b r, whose gradient term is second order, so at first order the perturbation is purely conformal: measured trace-free part 1.8e-04 against a trace of 3.164. Shape carries the trace and nothing else — which is exactly why the tangential slide has to be there, and why discrete ellipses were the only option before it was.

THE THEOREM. Adding ξ = ∇W and demanding tracelessness fixes the radial part completely, ρ = −½ΔW, and then the induced metric perturbation of the drawn surface IS the solved h_ab: measured against solved to 4.7e-04 of the peak amplitude, with a trace of 1.1e-03 and off-diagonal 6.7e-08. The surface is not an illustration of the tensor; it has the tensor as its own geometry.

THE QUADRUPOLE IS THE TEXTBOOK SHAPE. Feeding the exact ℓ = 2 mode returns ρ = P₂(cos d) to 1.0e-07, amplitude ratio 1.000000. The prolate–oblate picture of a quadrupole gravitational wave is not assumed here; it comes out of the projection.

THE FREE CONSTANT IS A RIGID TRANSLATION. C is the whole kernel of the construction — an ℓ = 1 displacement moves every point of the sphere by one vector, verified to 1.1e-16. Removing it leaves a residual dipole of -2.4e-18, and ℓ = 0 cannot appear at all because ∫ΔW dA = 0: a spin-2 wave can never breathe the sphere's area.

AND THE AREA HOLDS. The drawn surface keeps the round area to 4.8e-06 at gain 1e-03 — second order in ε, which is what trace-free means, measured on the surface rather than inferred from the tensor.

IT REACHES THE BULK. With the gain fixed from the run's own peak, the surface travels 74.8% of the way to R_outer at t = 3.085 and 44.3% toward R_inner, without touching either. The tensor wave now has an extrinsic amplitude — which is the whole reason to do this.

SCOPE. The gain is a display choice and the theorem is a first-order statement, so its residuals shrink with ε. This is a faithful representation of h_ab in the embedding; it does not make the wave act on the sphere it is deforming.