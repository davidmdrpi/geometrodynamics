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

## The surface

A unit 2-sphere with both polar caps `θ < a` and `θ > π − a` removed, and
the two mouth circles joined by a **catenoidal neck**
`r(s) = b·cosh((s − L/2)/b)` in arclength gauge. Demanding a `C¹` join —
the circumference *and* its slope continuous at both mouths — fixes the
neck completely from the single parameter `a`:

```
r(0) = sin a,  r'(0) = −cos a   ⟹   b = sin a / √(1 + cos²a),
                                    L = 2b·asinh(cos a)
```

In arclength gauge `r'' = r/b²`, so the neck carries **constant** Gauss
curvature `K = −1/b²`, and the two pieces cancel exactly in Gauss–Bonnet:

| piece | `∫K dA` |
|---|---:|
| capped sphere | `+4π cos a` |
| catenoidal neck | `−4π cos a` |
| **total** | **`0`** ⟹ `χ = 0` |

so the glued surface is a **torus** (azimuth-preserving gluing, `ε = +1`)
or a **Klein bottle** (azimuth-reversing, `ε = −1`). Gluing the neck frame
`(∂_s, ∂_ψ)` to the sphere frame `(∂_θ, ∂_φ)` gives `ds∧dψ = −dθ∧dφ` at
the north mouth and `−ε·dθ∧dφ` at the south, so a global orientation exists
iff `ε = +1`.

This is the 2-surface **section** of the S³ picture — `point → ring →
great circle → antipode` is the section of `point → shell → maximal shell →
antipode` — which is exactly why the low-dimensional movie is not a cartoon
of the high-dimensional one but a faithful slice of it.

At `a = 0.5` the neck is short and the outer way round is long:

| quantity | value |
|---|---:|
| neck waist `b` | 0.3603 |
| neck length `L` (inner route) | 0.5709 |
| `π − 2a` (outer route) | 2.1416 |
| **shortcut ratio** | **3.75×** |
| neck curvature `K = −1/b²` | −7.7014 |
| `χ = ∫K dA / 2π` | 0.0 (exact) |

## The controlled baseline

The comparison run is the **same domain on the same grid** with both mouths
sealed by a perfect mirror (`mode="plugged"`). Only the handle differs
between the two runs — there is no other change to compare away. A sealed
mouth is a filled-in throat, which is the physically meaningful control; a
"sphere with no holes at all" would change the domain as well.

## The answer (measured)

| finding | result |
|---|---|
| **the ring does not cross itself** | the front is a single connected circle from launch until it first reaches a mouth (`t = 1.071`); measured `1.104` plugged, `1.255` throat |
| **the echo delay is the neck length** | sealed echo at `2(θ₀−a)`, open return at `2(θ₀−a)+L`; delay measured `0.579` against `L = 0.571` — **1.5%**, stable to 1.1–1.6% across a 3× grid refinement |
| **the open mouth barely reflects** | opening the throat suppresses the mirror echo by **95.7%** |
| **the twist aims the arrival** | with `τ = π` the bulk route lands on the antipode `0.391` **ahead** of the geodesic focus and **9.7×** stronger than with `τ = 0` |
| **the orientation is real but hidden** | torus and Klein bottle are identical to a point source exactly at `τ ∈ {0, π}` (difference `0.0000`) and differ by ~18% elsewhere |
| **energy** | drift `< 7×10⁻⁴` over the whole coupled run |

### The ring is one circle in free flight

On a closed surface the front of a point pulse is the geodesic circle of
radius `t`: it expands to the great circle and contracts to the antipode
**without ever meeting itself**. Only the handle can put a second front on
the same surface, and on a closed surface two fronts must cross. Measured
by counting connected components of the smoothed energy density
`u_t² + |∇u|²` (not `|u|`, which splits one physical ring into concentric
shards at the crest/trough zero), the count is exactly `1` past the
free-flight time in both topologies.

The launch transient is excluded: a cap released from rest is d'Alembert
data and splits into an outgoing *and* an ingoing ring, which collapses
through the source and re-expands within about one pulse width. That is a
property of releasing from rest, not a wavefront crossing.

### The echo delay is the neck length

This is the load-bearing measurement, because the two routes share every
segment except the crossing:

| route | predicted | measured |
|---|---:|---:|
| sealed mirror echo `2(θ₀−a)` | 2.1416 | 2.0640 |
| open bulk return `2(θ₀−a)+L` | 2.7125 | 2.6434 |
| **delay = `L`** | **0.5709** | **0.5794** |

Both absolute times sit ≈ one pulse half-width early — the peak of a finite
pulse is not its geodesic front — but the bias is common to both, so it
cancels in the delay. **The wave measures the throat**: a classical scalar
field on a fixed surface reports the neck's arclength to 1.5% with no
fitted parameter anywhere.

### The mouth transmits rather than reflects

| mouth `a` | mouth radius | transmitted | reflected |
|---:|---:|---:|---:|
| 0.30 | 0.296 | 0.186 | 0.814 |
| 0.50 | 0.479 | 0.286 | 0.714 |
| 0.70 | 0.644 | 0.417 | 0.583 |

and the mirror echo that the sealed run produces is suppressed by 95.7%
when the throat is opened. The wave goes *through* the hole rather than
bouncing off it, and what comes back comes back later, by the bulk.

### The twist aims the bulk arrival

The gluing offset `τ = twist_steps · dφ` rotates where the bulk route
re-emerges, at no energy cost:

* `τ = 0` — the route runs source → north mouth → neck → south mouth →
  **back to the source**, a closed geodesic of length `2(θ₀−a)+L = 2.712`
  that generates `π₁` and does not exist without the handle;
* `τ = π` — the same route ends on the **antipode** instead, delivering a
  precursor there at `2.644` (predicted `2.712`), **0.391 ahead** of the
  geodesic focus at `3.035` (biased `π`), and 9.7× stronger than the
  untwisted throat delivers to the same point.

The geometry, not a coupling constant, decides where the energy lands.

### The orientation is real but hidden

A point source is symmetric under the reflection `R: φ → −φ` through its
own meridian. `R` survives the gluing `g: φ ↦ εφ + τ` iff

```
R∘g = g∘R  ⟺  −(εφ + τ) = ε(−φ) + τ  ⟺  τ ≡ −τ  ⟺  τ ∈ {0, π}
```

At those two offsets the mirror carries the torus gluing into the
Klein-bottle gluing, so **no point source can tell them apart**:

| `τ/π` | torus vs Klein difference | mirror broken? |
|---:|---:|---|
| 0.000 | 0.0000 | no |
| 0.250 | 0.1907 | yes |
| 0.500 | 0.1623 | yes |
| 0.750 | 0.1834 | yes |
| 1.000 | 0.0000 | no |

The measurement matches the argument to machine precision. This sharpens
the inner/outer intuition: the asymmetry is present at *every* twist — the
inner route is short and negatively curved, the outer route long and
positively curved, and reversing the neck's normal flips the extrinsic
sign — but **it takes a twist that breaks the source's mirror to make the
orientation observable**. Intrinsic spherical symmetry and extrinsic
orientation asymmetry coexist, and the observable consequences of the
second are gated by how the source is placed.

## Honest scope

* **Linear, no backreaction.** A focus here can sharpen but cannot nucleate
  a new throat, so this probe says nothing about the #175 threshold. It
  establishes the kinematic stage on which that question would be asked.
* **The join is `C¹`, not `C²`.** Each mouth carries a curvature ring which
  scatters a little on its own; that scattering is inside the reported
  mouth budget rather than removed from it.
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
from geometrodynamics.viz import ThroatWaveSim, run_throat_animation, plot_wavefront_panel
anim = run_throat_animation(ThroatWaveSim(mode="throat", twist_steps=96))

sim = ThroatWaveSim(mode="throat", twist_steps=96); sim.advance_to(1.4)
plot_wavefront_panel(sim)      # both hemispheres, the neck interior, the map
```

## Where this could go

1. **Backreaction.** Let the wave's stress-energy move `b`, and ask whether
   the focus can deepen the neck — the point at which #175's threshold
   becomes a question about *this* geometry rather than a reduced ring.
2. **Two throats.** Two handles give competing bulk routes and a genuine
   interference of geometric histories at a common focus.
3. **Up to S³.** The same construction with `S³` and a 3-dimensional neck:
   expanding `S²` shells, an embedded shell inside the throat, and the
   antipodal caustic in its native dimension.
4. **A moving object.** Replace the point source with a persistent
   structure and measure the momentum a passing pulse imparts to it — the
   reflection budget of T6 is the first term of that ledger.
