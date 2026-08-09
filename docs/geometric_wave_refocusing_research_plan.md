# Ring wavefronts on a surface with a throat (PR #242)

> **Framing (to avoid a category error).** A linear scalar wave on a *fixed
> classical* surface — geometry → field, **not** quantum gravity, and no
> backreaction anywhere.

## The question

The recent arc reached the throat through algebra: operators, pairings,
CHSH ceilings, winding charges. This probe goes back the other way and asks
the question the pictures were originally asking:

> **If we just let a classical wave run on a closed surface that has a
> throat, what does the geometry itself do to it?**

Nothing is fitted and no quantum ingredient is inserted. Every number below
is read off a live wave solve, and every stage of it can be watched.

## Three surfaces, one clock

The comparison needs all three terms, not two. Each isolates one change:

| | surface | what it isolates |
|---|---|---|
| 1 | **bare S²** — uncut | the canonical picture: point → expanding ring → great circle → contracting ring → antipodal focus |
| 2 | **plugged** — both caps cut out, mouths sealed by a mirror | what *cutting holes* does |
| 3 | **throat** — the same cut sphere, mouths joined by a neck | what *opening a second route* does |

The bare sphere is not decoration. A point source makes the field a function
of geodesic distance alone, so it is solved **exactly** by the axisymmetric
solver already in `viz/antipodal_focusing.py` and mapped onto the same
`(θ, φ)` grid — no new physics, no polar-grid pathology, and it supplies the
only case in which the front provably never meets another front.

## The neck is a true catenoid

The neck is the **minimal surface of revolution**, `H = 0`,
`r = b·cosh(z/b)`. Matching the circumference *and* its arclength slope to
the sphere at the mouth,

```
r = sin a,  dr/ds = −cos a   ⟹   tanh(z₀/b) = cos a,   b·cosh(z₀/b) = sin a
```

has the closed-form solution

```
b = sin²a,        L = sin 2a,        z₀ = b·artanh(cos a)
```

Its Gauss curvature `K = −b²/r⁴` **varies** — which is the point. It is
exactly **−1** at each mouth (the sphere's `+1` with its sign flipped, no
jump in magnitude) and deepens to `−1/sin⁴a` at the waist:

| `a` | `b = sin²a` | `L = sin 2a` | `\|r − b cosh(z/b)\|` | `K` mouth | `K` waist | `χ` |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.0612 | 0.4794 | 1.7e-16 | −1.000 | −266.92 | 0 |
| 0.50 | 0.2298 | 0.8415 | 2.2e-16 | −1.000 | −18.93 | 0 |
| 0.75 | 0.4646 | 0.9975 | 2.2e-16 | −1.000 | −4.63 | 0 |
| 1.00 | 0.7081 | 0.9093 | 2.2e-16 | −1.000 | −1.99 | 0 |
| 1.30 | 0.9284 | 0.5155 | 2.2e-16 | −1.000 | −1.16 | 0 |

So the wave crosses a real minimal-surface neck of continuously changing
negative curvature, not a manufactured constant-curvature strip.

### What χ = 0 actually tests

For *any* surface of revolution `K = −r''/r` and `dA = 2πr ds`, so

```
∫K dA = −2π[r']
```

depends **only on the boundary slopes**. The `C¹` join pins those to
`∓cos a`, so `4π cos a − 4π cos a = 0` closes for the catenoid and for every
other `C¹`-matched profile alike. **The closure is a check on the join, not
evidence for a particular neck** — a point the first version of this write-up
overstated, and one the test suite now pins with a deliberately different
(cubic Hermite) profile.

The glued surface is a **torus** (azimuth-preserving gluing) or a **Klein
bottle** (azimuth-reversing): `ds∧dψ = −dθ∧dφ` at the north mouth against
`−ε·dθ∧dφ` at the south, so a global orientation exists iff `ε = +1`.

At `a = 0.75` the geometry sets four times on one clock:

| quantity | value |
|---|---:|
| neck waist `b` | 0.4646 |
| neck length `L` (inner route) | 0.9975 |
| `π − 2a` (outer route) | 1.6416 |
| **shortcut ratio** | **1.65×** |
| ring reaches the mouths | 0.821 |
| sealed echo returns `2(θ₀−a)` | 1.642 |
| bulk crossing lands `2(θ₀−a)+L` | 2.639 |
| antipodal focus `π` | 3.142 |

## How the pieces are coupled

Each mouth is **one finite-volume face shared** by a sphere cell and a neck
cell. Its flux is evaluated once,

```
F = r_mouth · (f_outer − f_inner) / h,      h = ½(dθ + ds)
```

and handed to both cells with opposite signs, so the discrete divergence
theorem holds across the mouth. (The first version gave each domain its own
interpolated ghost ring; those two ghosts disagreed at `O(dθ)` and leaked.)
A sealed mouth is the zero-flux face — the same statement as a perfect
mirror, and conservative for the same reason.

Two further corrections make the diagnostics mean what they say:

* **the conserved quantity** is the exact leapfrog invariant, whose gradient
  term is the *cross* product `⟨∇uⁿ, ∇uⁿ⁻¹⟩` — the velocity lives at the
  half step, so pairing it with a same-time gradient leaves an `O(dt²)`
  wobble that peaks whenever a sharp front crosses something;
* **the launch** is a purely outgoing ring, `u⁻ = f(d + c·dt)`. A cap
  released from rest is d'Alembert data: it splits into an outgoing *and* an
  ingoing front, which puts two fronts on the surface for reasons that have
  nothing to do with the geometry.

A third correction belongs with them: a closed surface has no boundary, so
`d²/dt² ∫u dA = ∫Δu dA = 0` and the field's mean is a free mode that ramps
**linearly** whenever `∫u_t dA ≠ 0`. A one-way launch has exactly that
defect — it is a monopole as well as a ring — and the ramp lifts the whole
surface off zero and swamps the wake it is meant to reveal. Subtracting the
area-weighted mean of `u_t` pins the zero mode and leaves the ring untouched.

With all three in place the energy drift is `~10⁻¹⁵` — round-off — and the
integrated mouth power closes against the neck's stored energy to 1.8%.

## The answer (measured)

### Does a front ever meet another front?

Counting connected components of a level set could not answer this: a hole
*cutting* one ring into arcs raises the component count without any front
meeting another, and that is exactly what both the plugged and the open
surface do the moment the ring reaches a mouth. So the diagnostic now counts
**arrivals per point** — a per-cell hysteresis trigger on the energy density
`u_t² + |∇u|²`, armed above 35% and re-armed below 12% of that cell's own
peak. (Plain local-maximum counting also fails: a wave in 2+1 dimensions
violates Huygens' principle, so every front drags a rippling wake whose
ripples are local maxima too.)

Over `t < π`:

| surface | max arrivals | area with ≥2 | of the source side |
|---|---:|---:|---:|
| **bare S²** | 2 | **0.015** | **0.000** |
| plugged | 4 | 0.312 | **0.099** |
| throat | 3 | 0.443 | **0.000** |

On a closed surface with no boundary the front is the geodesic circle of
radius `t`: it sweeps each point exactly once, so it cannot meet itself —
1.5% of the bare surface ever sees a second front, and that only at the
antipodal caustic where the front converges and re-expands.

Sealing the mouths puts a second front back toward the source over 9.9% of
that hemisphere. Opening the throat puts a second front over *more* of the
surface (44%) but **none of it back home** — the mouth transmits instead of
reflecting. That is the echo result, resolved in space rather than in time.

### The echo delay is the neck length

The two routes share every segment except the crossing:

| route | predicted | measured |
|---|---:|---:|
| sealed mirror echo `2(θ₀−a)` | 1.6416 | 1.5956 |
| open bulk return `2(θ₀−a)+L` | 2.6391 | 2.5980 |
| **delay = `L`** | **0.9975** | **1.0024** |

0.49%, and 0.09–0.49% across a 2.25× grid refinement. Both absolute times sit
one pulse half-width early — the peak of a finite pulse is not its geodesic
front — but the bias is common to both and cancels in the delay. **The wave
measures the throat.**

### The mouth budget, by integrated flux

The first version reported the peak *stored* energy fraction inside the neck
and called it "transmitted". That is not a transmission coefficient: it is a
snapshot, it depends on how long the neck is, and it ignores everything that
has already passed through. The throughput is now measured by integrating
power through two surfaces — `offered` across a reference circle just inside
the mouth, `through` across the mouth face itself:

| mouth `a` | offered | through | transmission | reflection | peak stored |
|---:|---:|---:|---:|---:|---:|
| 0.40 | 0.4554 | 0.3596 | 0.790 | 0.210 | 0.236 |
| 0.55 | 0.5891 | 0.5076 | 0.862 | 0.138 | 0.336 |
| 0.75 | 0.7645 | 0.7018 | **0.918** | 0.082 | 0.467 |
| 0.90 | 0.8968 | 0.8408 | 0.937 | 0.063 | 0.558 |

On a closed surface only part of the wave ever reaches the mouth, so the
total energy is the wrong denominator; `offered` is the right one.

> **Mirror suppression and transmission are different measurements.** The
> sealed echo's *amplitude* at one watched point and time is suppressed by
> 85.5% when the throat is opened; the *energy* transmission at the mouth is
> 91.8%. They are two views of the same fact and must not be quoted
> interchangeably — the earlier PR description conflated them.

### The twist aims the bulk arrival

With `τ = π` the bulk route ends on the antipode at a predicted `2.6391`,
measured `2.6094` — **0.468 ahead** of the geodesic focus, and **3.2×**
stronger than the same neck at `τ = 0`, which sends the same energy back to
the source instead. The geometry, not a coupling constant, decides where it
lands.

### The orientation is real but hidden

A point source is symmetric under the reflection `R: φ → −φ` through its own
meridian, and `R` survives the gluing `g: φ ↦ εφ + τ` iff

```
R∘g = g∘R  ⟺  −(εφ + τ) = ε(−φ) + τ  ⟺  τ ≡ −τ  ⟺  τ ∈ {0, π}
```

| `τ/π` | torus vs Klein difference | mirror broken? |
|---:|---:|---|
| 0.000 | 0.0000 | no |
| 0.250 | 0.2766 | yes |
| 0.500 | 0.3377 | yes |
| 0.750 | 0.3201 | yes |
| 1.000 | 0.0000 | no |

Measurement matches argument to machine precision. The asymmetry is present
at *every* twist — the inner route is short and negatively curved, the outer
long and positively curved, and reversing the neck's normal flips the
extrinsic sign — but **it takes a twist that breaks the source's mirror to
make the orientation observable**. Intrinsic spherical symmetry and extrinsic
orientation asymmetry coexist, and the observable consequences of the second
are gated by how the source is placed.

## Honest scope

* **Linear, no backreaction.** A focus can sharpen but cannot nucleate, so
  this probe says nothing about the #175 threshold. It establishes the
  kinematic stage on which that question would be asked.
* **The join is `C¹`, not `C²`.** Each mouth carries a curvature ring which
  scatters a little on its own; that scattering is inside the reported mouth
  budget rather than removed from it.
* **`χ = 0` tests the join, not the profile.**
* **Arrival-time bias.** Absolute times carry a common leading-edge bias of
  about the pulse's half width. Every load-bearing number is a *difference*
  of arrival times taken on one grid with one pulse.
* **A 2-surface section, not S³.** The dimensional descent is faithful
  (`ring ↔ shell`, `great circle ↔ maximal shell`, `antipode ↔ antipode`),
  but a genuine `S³` with an embedded throat is a separate solve.

## Reproduce

```bash
python -m experiments.closure_ledger.geometric_wave_refocusing_probe
# Verdict: GEOMETRY_ALONE_ROUTES_THE_WAVE  (8/8)
```

Watch it (any interactive matplotlib backend):

```python
from geometrodynamics.viz import (
    BareSphereSim, ThroatWaveSim, plot_wavefront_panel, run_throat_animation)

anim = run_throat_animation(ThroatWaveSim(mode="throat", mouth_angle=0.75,
                                          twist_steps=96))

sim = ThroatWaveSim(mode="throat", mouth_angle=0.75, twist_steps=96)
sim.advance_to(1.4)
plot_wavefront_panel(sim)      # both hemispheres, the neck interior, the map
```

## Where this could go

This is meant as the beginning of a visual experimental platform rather than
one result: instead of asking geometry to reproduce a known formula, let the
geometry run and see what it produces.

1. **Backreaction.** Let the wave's stress-energy move `b`, and ask whether
   the focus can deepen the neck — the point at which #175's threshold
   becomes a question about *this* geometry rather than a reduced ring.
2. **Two throats.** Two handles give competing bulk routes and a genuine
   interference of geometric histories at a common focus.
3. **Up to S³.** The same construction with `S³` and a 3-dimensional neck:
   expanding `S²` shells, an embedded shell inside the throat, and the
   antipodal caustic in its native dimension.
4. **A moving object.** Replace the point source with a persistent structure
   and measure the momentum a passing pulse imparts to it — the reflection
   column of the mouth budget is the first term of that ledger.
