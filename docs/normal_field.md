# Draw the wave as vectors, and they intersect

> Still a **representation** of a scalar field on a fixed round `S²`. What
> changes is which object is drawn — the vectors, not their tips.

## The abstraction that was missing

Four rounds of this series established that a height field cannot wind
(`circle_slice_bulk.md`, `seam_scale.md`), cannot self-intersect at any gap or
gain (`ring_and_fold.md`), and needed an invented tangential freedom before
anything would fold.

All four were about the same object: **the graph of the displacement's tips**,
`r = f(σ)`. A graph is embedded by construction.

Draw the **vectors themselves** and the obstruction is gone, for a reason that
is entirely classical: *neighbouring normals to a curve meet at its centre of
curvature*. The normal family has an envelope — the evolute — and a normal of
length `L` crosses its neighbours as soon as `L` exceeds the local radius of
curvature `ρ = 1/κ`. Nothing is added by hand.

## The same wave, drawn two ways

At one instant, at the focus:

| object | self-intersections |
|---|---:|
| graph of the tips | **0** |
| normal field, `L = 0.26` | 146 |
| normal field, `L = 0.35` | 304 |
| normal field, `L = 0.50` | 520 |

Nothing about the field changed. Only which object was drawn.

## The threshold is the radius of curvature — and the wave drives it down

| t | `ρ_min` | first drawn crossing |
|---:|---:|---:|
| 1.2 | 0.0755 | 0.452 |
| 2.0 | 0.0762 | 0.432 |
| 2.6 | 0.0573 | 0.319 |
| 3.0 | **0.0338** | **0.182** |

The converging ring **sharpens its own surface** by `2.23×`, and the crossing
threshold falls with it. This is the ring concentration finally showing up as
something visible: not as height — which barely changes, and never beats the
launch — but as curvature, which is what decides whether the normals meet.

The drawn crossing always lags `ρ_min`, as it must: a finite sampling stride can
only ever trail the continuous envelope, never precede it. That the two move
together is the check.

## The reset is a second, separable mechanism

A normal long enough to leave through `R_outer` re-enters at `R_inner` — at the
angle where it left, continuing in the same direction, which is the same crossing
rule the slice itself uses, applied to the vector rather than to the surface.
The re-entered stub starts deep inside the annulus and shoots outward across
vectors it could never otherwise have reached.

| L | wrapped | without reset | with reset | added |
|---:|---:|---:|---:|---:|
| 0.20 | 22 | 30 | 36 | **6** |
| 0.26 | 435 | 146 | 152 | **6** |
| 0.35 | 501 | 304 | 402 | **98** |
| 0.50 | 501 | 520 | 810 | **290** |

Both mechanisms are real and they are counted separately, so neither is
smuggled into the other.

## And the gap matters again

This is what the height representation had severed. `ring_and_fold.md` showed
the fold threshold was *independent* of the gap — spread `0.0` across a fivefold
range — so shrinking the vacuole could never buy an intersection.

For normals, the vector length **is** what spans the gap, so they are one knob.
Drawing each bundle at `L = δ`, the length that just reaches across:

| δ | `L = δ` | normals alone | with reset |
|---:|---:|---:|---:|
| 0.40 | 0.40 | 386 | 386 |
| 0.26 | 0.26 | 146 | 152 |
| 0.16 | 0.16 | 0 | **206** |
| 0.09 | 0.09 | 0 | **522** |

At the tightest gap the normals alone are too short to reach each other — but
almost all of them wrap, and the stubs re-entering at `R_inner` cross everything.
**Reducing the distance between the shells now produces intersections rather than
being unable to.**

## What this changes about the earlier rounds

Nothing in them was wrong; they were all statements about a graph, and they hold.
What changes is their reach:

| | graph of the tips | normal field |
|---|---|---|
| self-intersects | never | above `L = ρ` |
| gap is a knob on it | no | yes |
| needs an invented `λ` | yes | no |
| the focusing shows up as | height (`0.935×` the launch) | curvature (`2.23×` sharper) |

## Scope

The vector length `L` is a display choice, like every gain in this series. The
directions and the curvature are the surface's own. What is derived is that the
crossing threshold is the radius of curvature, that the focusing drives it down,
and that the reset adds an independent mechanism on top.

One bug worth recording: `CircleSlice.sigma` deliberately keeps both `σ = −π` and
`σ = +π` so the drawn curve closes, but they are the *same point*. Sampling both
puts two coincident normals in the bundle, and the orientation test scores that
as a crossing at every length — including `1e-04`, at which nothing can possibly
have met. The sample excludes the duplicate, and a test holds it there.

## Reproduce

```bash
python -m experiments.closure_ledger.normal_field_probe
# Verdict: THE_NORMALS_INTERSECT  (5/5)

python -m pytest tests/test_viz_normal_field.py -q       # 13 passed
```
