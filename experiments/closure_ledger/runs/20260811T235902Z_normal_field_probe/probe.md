# Draw the wave as vectors, and they intersect

_2026-08-11T23:59:02+00:00_

## The same wave, drawn two ways

Graph of the tips: **0** self-intersections.

| vector length | normal crossings |
|---:|---:|
| 0.10 | 0 |
| 0.20 | 30 |
| 0.26 | 146 |
| 0.35 | 304 |
| 0.50 | 520 |

## The threshold is the radius of curvature

| t | `ρ_min` | first drawn crossing |
|---:|---:|---:|
| 1.2 | 0.0755 | 0.4518 |
| 2.0 | 0.0762 | 0.4316 |
| 2.6 | 0.0573 | 0.3192 |
| 3.0 | 0.0338 | 0.1819 |

The focus sharpens the surface by 2.23×.

## The reset adds more

| L | wrapped | without reset | with reset | added |
|---:|---:|---:|---:|---:|
| 0.20 | 21 | 30 | 32 | **2** |
| 0.26 | 434 | 146 | 148 | **2** |
| 0.35 | 500 | 304 | 398 | **94** |
| 0.50 | 500 | 520 | 804 | **284** |

## And the gap matters again

| δ | L = δ | normals alone | with reset |
|---:|---:|---:|---:|
| 0.40 | 0.40 | 386 | **386** |
| 0.26 | 0.26 | 146 | **148** |
| 0.16 | 0.16 | 0 | **178** |
| 0.09 | 0.09 | 0 | **470** |

## Verdict

**5/5 checks passed.**

**THE_NORMALS_INTERSECT.** THE NORMALS INTERSECT. Four rounds of this series established that a height field cannot wind and cannot cross itself. All four were about the same object — the graph of the displacement's tips — and drawing the displacement itself removes the obstruction entirely.

THE SAME WAVE, DRAWN TWO WAYS. At one instant the graph of the tips gives 0 self-intersections and the normal field gives 520. Nothing about the field changed and nothing was added by hand; only which object was drawn. The reason is classical: neighbouring normals to a curve meet at its centre of curvature, so the normal family has an envelope and a normal longer than the local radius of curvature has already crossed its neighbours.

THE THRESHOLD IS THE RADIUS OF CURVATURE, AND THE WAVE DRIVES IT DOWN. ρ_min falls from 0.0755 in mid-flight to 0.0338 at the focus, a factor of 2.23 — the converging ring sharpens its own surface — and the first drawn crossing falls with it, 0.452 to 0.182. The envelope bounds the drawn crossing from below at every time, as it must: a finite sampling stride can only lag the continuous envelope, never precede it.

AND THE RESET IS A SECOND MECHANISM ON TOP. A normal long enough to leave through R_outer re-enters at R_inner at the angle where it left, and that stub shoots outward from deep inside the annulus across vectors it could never otherwise have reached. At the focus the count goes from 304 to 398 at L = 0.35, and the addition grows with the stub, reaching 284 at the longest length tested. The two mechanisms are separable and both are real.

AND THE GAP MATTERS AGAIN. This is what the height representation had severed: slice_folding showed the fold threshold did not know the gap existed, so shrinking the vacuole could never buy an intersection. For normals the vector length IS what spans the gap, so they are one knob — and at the tightest gap the reset dominates completely: 0 crossings from the normals alone against 470 once they re-enter. Reducing the distance between the shells now produces intersections rather than being unable to.

SCOPE. The vector length is a display choice like every gain in this series; the directions and the curvature are the surface's own. What is derived is that the crossing threshold is the radius of curvature, that the focusing drives it down, and that the reset adds an independent mechanism.