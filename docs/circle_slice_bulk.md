# A circle slice, and a bulk the wave can wrap through

> A **representation** of the scalar field on a fixed round `S²`, cut by one
> great circle and drawn inside a glued annulus — not a derived boundary
> condition, and not backreaction.

Everything up to here has been the whole sphere. This cuts it with the great
circle through the source and its antipode, and watches the cross-section.

## The slice

Parametrised by `σ ∈ [−π, π)`, the geodesic distance from the source along that
circle is simply `d = |σ|`. So **one circle carries both halves of the wave**:
two lobes leaving the source in opposite directions, running around either
side, and meeting head-on at `σ = ±π` — the antipode.

Nothing is re-solved. The field is the 2-D solve sampled at `d(σ)`:

| quantity | value |
|---|---:|
| against the full `(θ, φ)` route | `1.4e-14` |
| mirror asymmetry `σ → −σ` | `6.3e-15` |
| lobe arrival time | `2.9908` |
| lobe positions at arrival | `−2.6922` / `+2.6922` |

That identification matters. A 1-D wave living on a circle would let two
counter-propagating pulses simply superpose to `2×` at the meeting point. The
slice instead inherits the sphere's real caustic, because what it is drawing is
a collapsing ring seen edge-on.

## The bulk, and its crossing rule

The slice lives in the vacuole — the annulus between `R_inner = 0.74` and
`R_outer = 1.26`. The crossing rule is the obvious one:

> a radius that would pass `R_outer` re-enters at `R_inner`.

So the wave that reaches up into the bulk **comes back inside the circle**.
That glues the two boundaries, makes the radial coordinate periodic with period
`gap = R_outer − R_inner`, and turns the space the curve lives on into a torus
`S¹_σ × S¹_ρ`.

### The threshold is exactly the half-gap over the peak

| quantity | value |
|---|---:|
| run peak `max\|u\|` | `1.179568` |
| predicted `(R_outer − R_mid)/peak` | `0.220420` |
| found by bisection | `0.220420` |
| relative error | `3.8e-16` |

Below it the slice stays in the annulus; above it, it wraps — and it wraps at
the antipode first, because that is where the amplitude is.

## What the bulk gives, and what it does not

On a torus a closed curve has an integer winding number in each direction, and
integers are exactly the stable objects a crossing rule is supposed to produce.

**This rule does not produce one.** Driving the gain up buys crossings and
never buys charge:

| gain / threshold | unsigned | signed | winding | sheets |
|---:|---:|---:|---:|---:|
| 1.0 | 0 | `+0` | `+0` | 1 |
| 1.6 | 2 | `+0` | `+0` | 2 |
| 2.4 | 2 | `+0` | `+0` | 2 |
| 3.6 | 4 | `+0` | `+0` | 2 |
| 5.0 | 6 | `+0` | `+0` | 3 |

The reason is not subtle and is worth stating plainly. The drawn curve is a
**graph** `r = f(σ)` over the circle with `f` single-valued, so `f(π) = f(−π)`
and its degree as a map `S¹ → S¹` is identically zero. Every outward crossing of
the seam is paid for by an inward one. **A height field cannot wind.**

This is checked two independent ways — a signed crossing ledger, and a degree
computed from unwrapped increments around the loop — which agree at every gain
and every time, including at absurd gains where the curve visits many sheets.

## What separates different waves

Driven at one common gain, four pulses meet the same bulk:

| pulse width | crossings | arc on the far sheet | arc / width | winding |
|---:|---:|---:|---:|---:|
| 0.36 | 4 | `0.1553` | 2.711 | `+0` |
| 0.24 | 4 | `0.0999` | 2.614 | `+0` |
| 0.14 | 4 | `0.0583` | 2.614 | `+0` |
| 0.08 | 4 | `0.0333` | 2.614 | `+0` |

The crossing *count* is the same for all of them — the launch amplitude is `1`
for every pulse and their peaks barely differ, so they all cross in and out the
same number of times. What varies is **how much of the circle rides the far
sheet**, a `4.7×` spread. Divided by the pulse width that is `2.61` for all of
them: the far-sheet arc simply *is* the pulse.

Normalising each pulse by its own wrap threshold instead would have made them
identical by construction, which is why they are all driven at one fixed gain.

## What this rules out

A crossing rule of this kind cannot manufacture topological charge from the
amplitude of a scalar height — however hard it is driven, however many times the
seam is crossed, and whatever the pulse looks like. A stable topological object
has to come from somewhere else: from a curve that is **free to stop being a
graph**.

That is a statement about what the next construction needs, not a defect in this
one. It also generalises past this particular gluing: any representation that
draws the bulk excursion as a height over the slice inherits the same zero.

## Scope

The crossing rule is a representation choice, not a derived boundary condition —
nothing here makes the wave dynamically aware of the seam, and the field is a
linear scalar on a fixed round background. The gain is a display choice; the
shape at any gain is the solved field.

## Reproduce

```bash
python -m experiments.closure_ledger.circle_slice_probe
# Verdict: THE_SEAM_IS_CROSSED_IN_PAIRS  (7/7)

python scripts/geometrodynamics_v46_circle_slice.py                 # animate
python scripts/geometrodynamics_v46_circle_slice.py --still out.png
python scripts/geometrodynamics_v46_circle_slice.py --waves out.png  # four waves

python -m pytest tests/test_viz_circle_slice.py -q                  # 18 passed
```

The renderer shows the slice in its annulus, a zoom on the antipode where the
two arms arrive, the same wave driven past threshold so it wraps through the
seam and reappears inside the circle, and the unrolled torus `(σ, ρ)` where the
seam is a line and the crossings are countable.
