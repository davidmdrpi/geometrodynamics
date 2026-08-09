# The restored geometry: one continuous S², warped by the wave it carries

> **Framing.** A *display* of a solved classical field as a radial displacement
> of a **fixed** closed surface. Not backreaction, not dynamical geometry.

## Why this exists

The archive's picture — `archive/geometrodynamics_v39.py` — was **one
continuous surface** whose radius carried the field, nested between two fixed
shells like Russian dolls. That is the object the BAM intuition was built on,
and it is what you could actually watch.

#242 and #243 replaced it with equirectangular maps, neck strips and meridional
sections. Those are correct and they measure things the picture cannot. They
are also not the picture. This restores the geometry first; projections and
slices go back on top of it, not in place of it.

## The object

```
r(θ, φ, t) = R_mid + Δ · tanh( g · u(θ, φ, t) / u_ref )
```

| | value | |
|---|---:|---|
| `R_inner` | 0.74 | the inner doll |
| `R_mid` | 1.00 | the surface the wave lives on |
| `R_outer` | 1.26 | the outer doll |
| half-gap `Δ` | 0.26 | |
| display gain `g` | 1.6 | |
| `u_ref` | 1.168 | the run's own peak, measured in a calibration pass |

The radii are deliberately the vacuole of `radial_caustic.py`, so the shell the
ring caustic lands on is the same inner doll drawn here.

## What is different from v39

v39 displaced the radius with a prescribed `sech²` envelope plus a mound grown
on a function of the frame number. Here the displacement is
`BareSphereSim` — an actual wave solve — so the deformation at the antipode
appears **because the wave focuses there**, and the probe measures that rather
than asserting it.

## One closed surface

| check | value |
|---|---:|
| φ seam mismatch | `2.4e-16` |
| pole spread | `0.0e+00` |
| carries its own poles | yes |

Nothing is cut out of it. That is what makes *"a pulse sweeps every point once
and fills its own void"* a statement about a closed manifold rather than about
a patch — and therefore why a **ring**, not a pulse, is what a throat needs.

## Nested, never touching

| quantity | value |
|---|---:|
| `r_min` over a full return | 0.7594 |
| `r_max` over a full return | 1.2226 |
| clearance, inner doll | 0.0194 |
| clearance, outer doll | 0.0374 |

The bound is structural rather than lucky: `tanh` cannot leave the gap.

## An honest display

`tanh` is strictly increasing, so the displacement preserves the **sign** of
the field everywhere and the **ordering** of every pair of amplitudes (worst
inversion `0.0e+00`), and no gain can flip a sign. It does **not** preserve
ratios — that distortion is the price of keeping the picture inside the vacuole,
and it is reported rather than hidden.

## The focus is measured, not keyed to the clock

| quantity | value |
|---|---:|
| deepest warp at geodesic distance | 3.141593 |
| error from π | `0.0e+00` |
| at time | 2.9814 |
| against π | 3.1416 |
| shortfall | 0.1602 |

The shortfall is the pulse's own width, not solver error — it shrinks
monotonically as the pulse narrows:

| pulse width | peak time | shortfall |
|---:|---:|---:|
| 0.24 | 2.9452 | 0.1963 |
| 0.18 | 2.9814 | 0.1602 |
| 0.12 | 3.0363 | 0.1052 |
| 0.08 | 3.0725 | 0.0691 |

## And it rings

The arrival is not a single mound:

| | fraction of the gap | at time |
|---|---:|---:|
| driven out toward `R_outer` | 85.6% | 2.988 |
| pulled in toward `R_inner` | 79.6% | 3.211 |

The focus arrives as a mound, inverts, and ends up pulling the geometry
**inward** — toward the inner shell, which is the one the ring caustic lands on.
v39 drew an outward spike there because its mound was prescribed outward.

## Honest scope

* **Not backreaction.** The wave does not feel the surface it is deforming, so
  no throat forms here and nothing about nucleation follows.
* **The displacement is a rendering choice**, bounded on purpose. Sign and
  ordering survive it; ratios do not.
* **A point source, on a closed surface.** The ring case — the one that can
  fold — is `radial_caustic.py`, and putting a ring on *this* surface is the
  obvious next step now that the surface is back.
* **The dolls are fixed.** They mark the vacuole; they are not solved for.

## Reproduce

```bash
python -m experiments.closure_ledger.warped_sphere_geometry_probe
# Verdict: THE_SURFACE_ITSELF_CARRIES_THE_WAVE  (6/6)

python scripts/geometrodynamics_v41_warped_sphere.py             # animate
python scripts/geometrodynamics_v41_warped_sphere.py --still sheet.png
python scripts/geometrodynamics_v41_warped_sphere.py --save out.gif
```

```python
from geometrodynamics.viz.warped_sphere import WarpedSphere
s = WarpedSphere()
s.advance_to(3.05)
s.excursion()      # how far toward each doll, right now
s.mesh()           # X, Y, Z of the warped surface — closed
```

## Where this goes next

1. **The ring on this surface.** `radial_caustic` says a ring is what folds;
   this surface is where a ring can be watched folding.
2. **The throat, drawn on the same object.** #242's catenoidal neck currently
   lives on a capped sphere with mouths; joining it to this continuous surface
   would put the wave, the warp and the handle in one picture.
3. **Backreaction.** The only way the displacement stops being a display.
