# What a wavefront has to be like to fold (PR #243)

> **Framing.** Ray and wavefront geometry on a *fixed classical* vacuole with a
> **flat** bulk — geometry → field, **not** quantum gravity. Nothing here is
> dynamical.

## The question

`geometric_wave_refocusing_probe` runs waves on a surface that *already has* a
throat. This is the prior question:

> **What kind of wavefront can fold at all?**

It is answered with the differential geometry of wavefronts. The closed forms
are exact, and the topology is measured **independently of them**, from the
front's own area element.

## The vacuole

Two concentric spheres, `R_inner` and `R_outer`, with a flat bulk between them.
Along a straight ray the impact parameter about the centre, `b = r sin α`, is
conserved, so a ray launched inward reaches radius `b` and no deeper.

| quantity | value |
|---|---:|
| `R_inner` | 0.7400 |
| `R_outer` | 1.2600 |
| `ΔR` | 0.5200 |
| grazing `sin α = R_in/R_out` | 0.5873 |
| grazing angle | 35.97° |
| critical ring `θ₀` | 54.03° |
| ring radius `ρ` | 1.0198 |
| **pulse crosses at** | **0.5200** |
| **ring's first caustic at** | **1.0198** |

## A point source does not fold — *in this flat bulk*

The front of a **point** is the metric sphere `|x − P| = t`, whose signed area
element is `t² sin θ`: positive for every `t > 0`. It never folds, and the
region behind it is the filled ball — the pulse *fills its own void*.

> **This is a statement about flat Euclidean space, not about point sources.**
> On a closed manifold the same source folds: on `S²` or `S³` a point's front
> converges on the antipode at `t = πR`, which is exactly what
> `throat_wavefront.py` measures. Curvature is what gives a point source a focal
> set; a flat bulk is what denies it one.

## A circle folds *coherently* — that is what is special

Any curved extended source has a focal set, so "only a ring folds" would be
false. What singles out the circle is not *that* it folds but *how*: by symmetry
the **whole ring arrives at one point simultaneously**.

The signed area element of the offset tube is `t(ρ + t cos v)`, which vanishes
where `cos v = −ρ/t`:

| | behaviour |
|---|---|
| `t < ρ` | immersed — no fold anywhere |
| `t = ρ` | the **first caustic**, infinitely degenerate: both roots coincide at the ring's centre and the *whole ring* arrives at once |
| `t > ρ` | **still singular**, at two axis points `z = ±√(t² − ρ²)` which run outward as the front grows |

So the fold is not an isolated event: a maximally degenerate first focus,
followed by a persistent singular locus along the symmetry axis.

### Measured independently

The fold time is found by scanning the front's own area element for a **sign
change** and bisecting — no radius of curvature is consulted, so the comparison
with `ρ` is a test rather than a tautology.

| | folds? | detected | closed form | error |
|---|---|---:|---:|---:|
| point source (flat bulk) | **no**, scanned to 2.20 | — | — | — |
| ring source | **yes** | 1.019806 | 1.019804 | `2.0e-06` |

and on four rings unrelated to this shell:

| ring `θ₀` | `ρ` | detected | error |
|---:|---:|---:|---:|
| 0.30 | 0.372355 | 0.372356 | `7.4e-07` |
| 0.70 | 0.811714 | 0.811716 | `1.6e-06` |
| 1.10 | 1.122921 | 1.122924 | `2.2e-06` |
| 1.40 | 1.241667 | 1.241669 | `2.5e-06` |

Two details make the detector mean what it says. The test is **relative**
(`min J < −tol·max|J|`), because a parametrisation whose area element merely
*vanishes* — the direction sphere's own poles do that at every `t` — has no
caustic; a fold is a change of *sign*. And the orientation is **referenced** at
small `t`, because `(X_u × X_v)·N` carries the handedness of whatever `(u, v)`
ordering a source happens to use: for the point source it is negative
everywhere, which an unreferenced test would read as an instant fold.

At the first caustic the ring's points are equidistant from the centre to
`4.4e-16` — degenerate, the whole ring — against multiplicity `2` just off it.

## The core result: the two conditions coincide

Ask for the ring whose first caustic lands *on the inner sphere*. Its centre
must sit at radius `R_inner`, which fixes `cos θ₀ = R_inner/R_outer` — and that
same ring's rays leave at `sin α = R_inner/R_outer`, exactly the **grazing ray**
tangent to the inner sphere.

| quantity | value |
|---|---:|
| first caustic at radius | 0.740000 |
| `R_inner` | 0.740000 |
| error | **0.0e+00** |
| launch `sin α` | 0.587302 |
| grazing `sin α` | 0.587302 |
| error | **0.0e+00** |
| ray turning radius | 0.740000 |

**The ring that focuses on the throat and the ray that grazes it are the same
ring**, and it forms at `t = ρ = √(R_outer² − R_inner²)`.

## Acceptance is asymmetric; propagation is not

| direction | closed form | Monte-Carlo |
|---|---:|---:|
| outer → inner | 0.1906 | 0.1927 |
| inner → outer | 1.0000 | 1.0000 |

Only launch directions with `sin α ≤ R_in/R_out` reach the inner sphere — about
**19%** of the inward hemisphere — while **every** direction from the inner
sphere reaches the outer one. A **5.2×** difference.

> **This is an angular (solid-angle) acceptance asymmetry, not nonreciprocal
> propagation.** Every individual ray is exactly reversible: `b` is unchanged
> under reversal and the accepted inward rays all have `b ≤ R_in`, so each one
> climbs back out along its own reverse — the probe asserts this rather than
> assuming it. What differs is the *measure of launch directions that connect*,
> because a hemisphere at `R_out` and a hemisphere at `R_in` are different sets
> of directions. No symmetry of the sphere is broken; the ordering of the two
> radii is the whole of it.

## It scales with the ratio, not this shell

| `R_in` | `sin α_crit` | `θ₀` | `ρ` | caustic at | pulse at | inward accept |
|---:|---:|---:|---:|---:|---:|---:|
| 0.20 | 0.1587 | 80.9° | 1.2440 | 1.2440 | 1.0600 | 0.0127 |
| 0.45 | 0.3571 | 69.1° | 1.1769 | 1.1769 | 0.8100 | 0.0660 |
| 0.74 | 0.5873 | 54.0° | 1.0198 | 1.0198 | 0.5200 | 0.1906 |
| 0.95 | 0.7540 | 41.1° | 0.8277 | 0.8277 | 0.3100 | 0.3431 |
| 1.15 | 0.9127 | 24.1° | 0.5149 | 0.5149 | 0.1100 | 0.5914 |

The caustic lands on the inner sphere at every ratio (error `< 1e-12`), always
*later* than the pulse crosses, and the acceptance tightens as the inner sphere
shrinks.

## Honest scope

* **A flat bulk.** Exact, and independent of any wave solve. `ShellGeometry`
  accepts a metric factor `f` and works in `r/√f(r)`, but **only where that is
  monotone on the shell** — it is validated at construction and refused
  otherwise, because a photon sphere admits trapped orbits and two-sided turning
  points these closed forms do not describe.
* **The point-source claim is scoped to the flat bulk** and does not survive on
  a closed manifold, where the same source focuses at the antipode.
* **A ring is not the only source that folds** — any curved extended source
  does. Its distinction is coherence: the whole ring at one point.
* **Nothing here is dynamical.** It says which fronts *can* fold, not what
  happens when one does.
* **No throat is shown forming.** That needs backreaction, which this programme
  does not have. The renderer marks the computed wavefront caustic faithfully
  and labels the would-be throat as schematic.
* **A wavefront caustic is not yet a topological defect of the geometry.** The
  front stops being immersed; the manifold does not change. Turning the first
  into the second is exactly the open step.

## Reproduce

```bash
python -m experiments.closure_ledger.ring_caustic_defect_probe
# Verdict: WHOLE_RING_FOCUSES_AT_ONE_POINT  (8/8)

python scripts/geometrodynamics_v40_ring_defect.py --still ring_defect.png
python scripts/geometrodynamics_v40_ring_defect.py          # animate
```

```python
from geometrodynamics.viz.radial_caustic import ShellGeometry, plot_pulse_vs_ring
sh = ShellGeometry(0.74, 1.26)
plot_pulse_vs_ring(sh)          # the whole argument in one figure
```

## Where this could go

1. **Backreaction.** Let the focused energy move the shell and ask whether the
   wavefront caustic becomes a caustic *of the geometry*. That is the step this
   probe deliberately does not take.
2. **The ring's own origin.** A ring is what an antipodal focus on the outer
   surface produces — so the natural chain is *point → antipodal ring → bulk
   caustic*, and the middle link is what `throat_wavefront.py` already solves.
3. **Curved bulk, done properly.** A monotone `R_eff` already works. Handling a
   photon sphere means abandoning the one-sided turning-point argument and
   treating trapped orbits, which is a different piece of work.
