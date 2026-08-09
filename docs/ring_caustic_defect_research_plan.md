# Why a throat needs a ring: front topology across the bulk (PR #242)

> **Framing.** Ray and wavefront geometry on a *fixed classical* vacuole —
> geometry → field, **not** quantum gravity. Nothing here is dynamical.

## The question

`geometric_wave_refocusing_probe` runs waves on a surface that *already has* a
throat. This is the prior question:

> **What kind of wavefront can make one?**

It is answered with the differential geometry of wavefronts — focal sets —
rather than with a simulation. Every number is a closed form, checked against
an independent numerical count.

## The vacuole

Two concentric spheres, `R_inner` and `R_outer`, with the bulk between them.
Along a straight ray the impact parameter about the centre,

```
b = r sin α          (α = angle to the radial direction)
```

is conserved, so a ray launched inward reaches radius `b` and no deeper.

At the programme's own vacuole (`R_MID = 1`, `ΔR/2 = 0.26`):

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
| **ring defect at** | **1.0198** |

## A pulse cannot make a defect

The wavefront of a **point** source is the metric sphere `|x − P| = t`. A
sphere is embedded for every `t`: the front never touches itself, and the region
behind it is the filled ball. The pulse *fills its own void*.

The reason is exact and is a statement about focal sets: **the focal set of a
point is empty**, so the front has nothing to fold on. It crosses the bulk at
`t = ΔR` and that is all it does.

## A ring must

The wavefront of a **circle** of radius `ρ` is the offset (tube) surface of that
circle — and a curve *does* have a focal set: the locus of its centres of
curvature. For a circle every point shares the **same** centre, so the focal set
collapses to a single point which the **entire ring reaches simultaneously** at

```
t = ρ
```

Measured, the ring's points are equidistant from that centre to `4.4e-16`. That
is a degenerate caustic of infinite multiplicity — not a smooth focus but a
point where the front ceases to be embedded. A **codimension-2 defect of the
wavefront, made by geometry alone.**

| | focal set | self-intersects | when |
|---|---:|---|---:|
| point source | 0 points | **never** | — |
| ring source | 1 point | **yes** | `t = ρ = 1.0198` |

Multiplicity at the focus is degenerate (the whole ring); one step off the axis
it is 2; outside the front it is 0.

## The two conditions coincide

Ask for the ring whose defect lands *on the inner sphere*. Its centre must sit
at radius `R_inner`, which fixes

```
cos θ₀ = R_inner / R_outer
```

and that same ring's rays leave the outer sphere at `sin α = R_inner/R_outer` —
exactly the **grazing ray**, tangent to the inner sphere. Both errors are
identically zero:

| quantity | value |
|---|---:|
| defect forms at radius | 0.740000 |
| `R_inner` | 0.740000 |
| error | **0.0e+00** |
| launch `sin α` | 0.587302 |
| grazing `sin α` | 0.587302 |
| error | **0.0e+00** |
| ray turning radius | 0.740000 |

**The ring that focuses on the throat and the ray that grazes it are the same
ring.** The defect forms at `t = ρ = √(R_outer² − R_inner²)`.

## The bulk is not symmetric between its faces

Because `b = r sin α` is conserved and `r` decreases going in:

| direction | closed form | Monte-Carlo |
|---|---:|---:|
| outer → inner | 0.1906 | 0.1927 |
| inner → outer | 1.0000 | 1.0000 |

Only rays with `sin α ≤ R_in/R_out` reach the inner sphere — about **19%** of
the inward hemisphere — while **every** ray from the inner sphere escapes. A
**5.2×** asymmetry across the same bulk, in opposite directions.

This is the inner/outer asymmetry in its plainest form: **nothing is broken.**
The sphere is as symmetric as it ever was; the asymmetry is the *ordering of the
two radii*, and it appears the moment you ask a ray to go one way rather than
the other.

## It scales with the ratio, not this shell

| `R_in` | `sin α_crit` | `θ₀` | `ρ` | defect at | pulse at | inward accept |
|---:|---:|---:|---:|---:|---:|---:|
| 0.20 | 0.1587 | 80.9° | 1.2440 | 1.2440 | 1.0600 | 0.0127 |
| 0.45 | 0.3571 | 69.1° | 1.1769 | 1.1769 | 0.8100 | 0.0660 |
| 0.74 | 0.5873 | 54.0° | 1.0198 | 1.0198 | 0.5200 | 0.1906 |
| 0.95 | 0.7540 | 41.1° | 0.8277 | 0.8277 | 0.3100 | 0.3431 |
| 1.15 | 0.9127 | 24.1° | 0.5149 | 0.5149 | 0.1100 | 0.5914 |

The defect lands on the inner sphere at every ratio (error `< 1e-12`), always
*later* than the pulse crosses, and the acceptance tightens as the inner sphere
shrinks.

## Honest scope

* **Ray and front geometry in a flat bulk.** Exact, and independent of any wave
  solve. A curved bulk replaces `r` by the effective radius `r/√f(r)` and the
  structure carries over — `ShellGeometry` accepts an `f` — but the closed forms
  quoted here are the flat ones.
* **Nothing here is dynamical.** It says which fronts *can* fold, not what
  happens when one does.
* **No throat is shown forming.** That needs backreaction, which this programme
  does not have. The renderer marks the computed wavefront defect faithfully and
  labels the would-be throat as schematic; the shell deformation it draws is a
  presentation cue keyed to the computed defect time, not a solution.
* **A wavefront defect is not yet a topological defect of the geometry.** The
  front stops being embedded; the manifold does not change. Turning the first
  into the second is exactly the open step.

## Reproduce

```bash
python -m experiments.closure_ledger.ring_caustic_defect_probe
# Verdict: ONLY_A_RING_FOLDS  (8/8)

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
   wavefront defect becomes a defect *of the geometry*. That is the step this
   probe deliberately does not take.
2. **The ring's own origin.** A ring is what an antipodal focus on the outer
   surface produces — so the natural chain is *point → antipodal ring → bulk
   defect*, and the middle link is what `throat_wavefront.py` already solves.
3. **Curved bulk.** Feed the 5D Tangherlini `f(r)`; the effective radius
   develops a minimum (a photon sphere), which adds a second grazing condition
   the flat case does not have.
