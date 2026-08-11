# A circle slice, and a bulk the wave can wrap through

_2026-08-11T06:09:49+00:00_

> A **representation** of the scalar field on a fixed round `S²`, cut by one great circle — not a derived boundary condition.

## The slice is the sphere

| quantity | value |
|---|---:|
| against the (θ, φ) route | `1.4e-14` |
| mirror asymmetry `σ → −σ` | `6.3e-15` |
| lobe arrival time | 2.9908 |
| lobe positions | -2.6922 / +2.6922 |

## The wrap threshold

| quantity | value |
|---|---:|
| run peak `max\|u\|` | 1.179568 |
| predicted `(R_outer − R_mid)/peak` | 0.220420 |
| found by bisection | 0.220420 |
| relative error | `3.8e-16` |

## Amplitude buys crossings, never charge

| gain / threshold | unsigned | signed | winding | sheets |
|---:|---:|---:|---:|---:|
| 1.0 | 0 | +0 | +0 | 1 |
| 1.6 | 2 | +0 | +0 | 2 |
| 2.4 | 2 | +0 | +0 | 2 |
| 3.6 | 4 | +0 | +0 | 2 |
| 5.0 | 6 | +0 | +0 | 3 |

The curve is a graph `r = f(σ)` with `f` single-valued, so every outward crossing is paid for by an inward one. **A height field cannot wind.**

## Different waves, one bulk, one gain

| pulse width | crossings | arc on the far sheet | arc / width | winding |
|---:|---:|---:|---:|---:|
| 0.36 | 4 | 0.1553 | 2.711 | +0 |
| 0.24 | 4 | 0.0999 | 2.614 | +0 |
| 0.14 | 4 | 0.0583 | 2.614 | +0 |
| 0.08 | 4 | 0.0333 | 2.614 | +0 |

A 4.7× spread in arc, 2.64 per unit pulse width for all of them — the far-sheet arc **is** the pulse.

## Verdict

**7/7 checks passed.**

**THE_SEAM_IS_CROSSED_IN_PAIRS.** THE SEAM IS CROSSED IN PAIRS. Gluing the vacuole's outer boundary to its inner one lets the wave that reaches into the bulk come back inside the circle, and makes the space the curve lives on a torus — but it does not give the curve anywhere to wind.

THE SLICE IS THE SPHERE, NOT A NEW PROBLEM. The field on the circle matches the sphere's field at geodesic distance |σ|, reached through the full (θ, φ) route, to 1.4e-14, and the two arms agree to 6.3e-15. So the slice inherits the real caustic rather than the 2× superposition a 1-D wave on a circle would give. The lobes leave in opposite directions and arrive together at t = 2.991, symmetric to 4.4e-16.

THE WRAP THRESHOLD IS EXACTLY THE HALF-GAP OVER THE PEAK. ε_crit = 0.220420 predicted, 0.220420 found by bisection — a relative error of 3.8e-16. Below it the slice stays in the annulus; above it, it wraps.

AND THEN NOTHING ACCUMULATES. Driving the gain up buys unsigned crossings — 6 of them at the hardest drive tested — while the SIGNED total stays at 0 and the winding number at 0, at every gain and every time, with the crossing ledger and an independently computed degree agreeing. The reason is not subtle and is worth stating plainly: the curve is a GRAPH r = f(σ) over the circle with f single-valued, so every outward crossing of the seam is paid for by an inward one. A height field cannot wind.

DIFFERENT WAVES DIFFER IN ARC, NOT IN TOPOLOGY. At one common gain, pulses from 0.36 down to 0.08 cross the seam the same number of times but put 0.155 down to 0.033 of the circle on the far sheet — a 4.7× spread. Divided by the pulse width that is 2.64 for all of them, spread 0.10: the far-sheet arc simply IS the pulse. The winding column is zero all the way down.

WHAT THIS RULES OUT. A crossing rule of this kind cannot manufacture topological charge from the amplitude of a scalar height, however hard it is driven and however many times the seam is crossed. A stable topological object would have to come from somewhere else — from a curve that is free to stop being a graph, which is a statement about what the next construction needs rather than a defect in this one.

SCOPE. The crossing rule is a representation choice, not a derived boundary condition: nothing here makes the wave dynamically aware of the seam, and the field is a linear scalar on a fixed round background.