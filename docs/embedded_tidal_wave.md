# Projecting `h_ab` into a continuous embedded surface

> **Framing.** A faithful *representation* of a spin-2 field in the ℝ³
> embedding of a **fixed** `S²` — not backreaction, not linearised GR on a
> spacetime.

## The problem

`spin2_tidal.py` samples the tensor as discrete tangent-space ellipses.
Faithful, but flat: it never touches the embedding, so it gives no intuition
for how tidal shear meets a bulk. Can `h_ab` be drawn as a **continuous**
surface deformation `r(θ, φ, t)` instead?

## Why a height field cannot do it

For `X = r n̂` the induced metric is `g_ab = r² ĝ_ab + ∂_a r ∂_b r`, and the
gradient term is **second** order. So at first order

```
δg_ab = 2ρ ĝ_ab      — purely conformal
```

Shape carries the trace and nothing else. Measured on a deliberately radial
surface:

| `d` | trace | trace-free part |
|---:|---:|---:|
| 0.60 | +3.1642 | `+1.8e-04` |
| 1.10 | +1.3230 | `-5.6e-05` |
| 1.90 | -1.6259 | `+7.1e-05` |
| 2.50 | -3.4171 | `-1.6e-04` |

That is exactly why discrete ellipses were the only option before the
tangential part was added.

## The displacement that does

```
X = (R + ερ) n̂ + ε ξ ,   ξ = ∇W ,   δg_ab = 2ρ ĝ_ab + 2∇₍ₐξ_b₎
```

Demanding tracelessness fixes the radial part completely:

```
ρ = −½ ΔW        h_ab = [2∇₍ₐ∇_b₎W]^TF
```

**One potential carries both** — the shear you measure is its trace-free
Hessian, the shape you see is minus half its Laplacian. The shear lives in how
the material *slides*; the shape carries what is left.

For an axisymmetric `h` the Hessian condition `W'' − cot d W' = h` is first
order in `ψ = W'`, with integrating factor `sin d`:

```
ψ(d) = sin d [ C + ∫₀^d h/sin ]        ρ(d) = −½h − cos d [ C + ∫₀^d h/sin ]
```

No derivative of the solved field is ever taken — one integral, whose
integrand is regular because `h ~ sin²d` at both poles.

## The theorem

The induced metric perturbation of the drawn surface **is** the solved `h_ab`:

| `d` | `h₊` measured | `h₊` solved | trace |
|---:|---:|---:|---:|
| 0.40 | +0.000024 | +0.000023 | `-6.4e-07` |
| 0.90 | +0.000802 | +0.000800 | `+3.2e-06` |
| 1.50 | +0.002622 | +0.002625 | `-6.2e-06` |

Worst relative error `4.7e-04`, trace `1.1e-03`, off-diagonal `6.7e-08`, at
gain `1e-04`. The surface is not an illustration of the tensor; it *has* the
tensor as its own geometry.

## What falls out

**The quadrupole is the textbook shape.** Feeding the exact `ℓ = 2` mode
returns `ρ = P₂(cos d)` to `1.0e-07`, amplitude ratio `1.000000`. The
prolate–oblate picture of a quadrupole gravitational wave is not assumed here —
it comes out of the projection.

**The free constant is a rigid translation.** `C` is the whole kernel: an
`ℓ = 1` displacement moves every point by one vector, verified to `1.1e-16`.
Removing it leaves a residual dipole of `2.4e-18`. And `ℓ = 0` cannot appear at
all, since `∫ΔW dA = 0` — a spin-2 wave can never breathe the sphere's area.

**Area holds at first order.** The drawn surface keeps the round area to
`4.8e-06` at gain `1e-03` — second order in `ε`, which is what trace-free
means, measured on the surface rather than inferred from the tensor.

**It reaches the bulk.** With the gain fixed from the run's own peak, the
surface travels `74.8%` of the way to `R_outer` at `t = 3.085` and `44.3%`
toward `R_inner`, touching neither.

## Honest scope

* **A representation, not backreaction.** The wave now has an extrinsic
  amplitude; it still does not act on the sphere it deforms.
* **The gain `ε` is a display choice.** The *shape* at any gain is the solved
  field, and the theorem is first order, so its residuals shrink with `ε`.
* **Still a fixed background.** 2+1 gravity has no propagating tensor modes;
  this is the spin-2 analogue of the scalar wave, drawn where it can be seen.
* **The tangential slide is invisible as shape.** It is drawn as the material
  lattice, because a tangential displacement is a reparametrisation — that is
  the same fact as "shape carries only the trace", seen from the other side.

## Reproduce

```bash
python -m experiments.closure_ledger.embedded_wave_probe
# Verdict: ONE_POTENTIAL_CARRIES_SHAPE_AND_SHEAR  (7/7)

python scripts/geometrodynamics_v44_embedded_wave.py             # animate
python scripts/geometrodynamics_v44_embedded_wave.py --still sheet.png
```

```python
from geometrodynamics.viz.embedded_wave import EmbeddedTidalSurface
s = EmbeddedTidalSurface()
s.advance_to(3.05)
s.profiles()                       # shear h, slide ξ, radial ρ — one grid
s.induced_metric_perturbation(1.5) # what the drawn surface actually is
s.mesh()                           # the continuous deformed surface
```

## Where this goes next

1. **A ring source, embedded.** The spin-2 field insists on a ring and the
   focal geometry says a ring is what folds — now both can be drawn on the
   same surface, in the bulk.
2. **The throat.** A catenoidal neck has no degenerate frame but strong
   curvature; the same projection applies and the shape is not a sphere.
3. **Backreaction**, still the only thing that would turn any of this from a
   representation into a dynamics.
