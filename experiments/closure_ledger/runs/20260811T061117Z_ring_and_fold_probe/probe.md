# The ring reaches across, but only a fold crosses

_2026-08-11T06:11:17+00:00_

## The ring grows as it converges

| quantity | value |
|---|---:|
| equator height (d = 1.507) | 0.1563 |
| peak height | 0.9331 |
| growth | 5.97× |
| against `1/√(sin d)` | 1.034 ± 0.071 |

## Reaching across — and the lead it gets on the pulse

| δ | ε | dwell | spans from d | lead |
|---:|---:|---:|---:|---:|
| 0.26 | 0.40 | 0.062 | 3.142 | 0.209 |
| 0.26 | 0.80 | 0.158 | 2.893 | 0.366 |
| 0.16 | 0.20 | 0.029 | 3.142 | 0.182 |
| 0.16 | 0.40 | 0.108 | 2.990 | 0.287 |
| 0.16 | 0.80 | 0.371 | 2.437 | 0.811 |
| 0.09 | 0.10 | 0.021 | 3.142 | 0.156 |
| 0.09 | 0.20 | 0.096 | 3.037 | 0.261 |
| 0.09 | 0.40 | 0.296 | 2.597 | 0.653 |
| 0.09 | 0.80 | 1.000 | 1.832 | 1.413 |
| 0.05 | 0.10 | 0.096 | 3.037 | 0.261 |
| 0.05 | 0.20 | 0.250 | 2.701 | 0.550 |
| 0.05 | 0.40 | 1.000 | 1.832 | 1.413 |
| 0.05 | 0.80 | 1.000 | 1.832 | 1.413 |
| 0.03 | 0.10 | 0.179 | 2.838 | 0.418 |
| 0.03 | 0.20 | 0.771 | 1.832 | 1.413 |
| 0.03 | 0.40 | 1.000 | 1.832 | 1.413 |
| 0.03 | 0.80 | 1.000 | 1.832 | 1.413 |

## ...and never intersecting

| δ | ε | seam crossings | self-intersections |
|---:|---:|---:|---:|
| 0.26 | 0.5 | 2 | **0** |
| 0.26 | 2.0 | 8 | **0** |
| 0.26 | 10.0 | 38 | **0** |
| 0.09 | 0.5 | 6 | **0** |
| 0.09 | 2.0 | 22 | **0** |
| 0.09 | 10.0 | 114 | **0** |
| 0.03 | 0.5 | 16 | **0** |
| 0.03 | 2.0 | 72 | **0** |
| 0.03 | 10.0 | 346 | **0** |

## The fold threshold

| quantity | value |
|---|---:|
| predicted `λε` | 0.012098 |
| found by bisection | 0.012098 |
| relative error | `1.8e-12` |
| folds at, from the antipode | 0.0000 |
| convergence drift | `1.5e-04` |

## The two knobs are orthogonal

| δ | span threshold | fold threshold |
|---:|---:|---:|
| 0.26 | 0.260 | 0.012362 |
| 0.12 | 0.120 | 0.012362 |
| 0.05 | 0.050 | 0.012362 |

Fold-threshold spread `0.0e+00`, and `λε ≈ 0.374 w²` across the pulse widths (spread 0.036).

## Verdict

**8/8 checks passed.**

**THE_RING_REACHES_BUT_ONLY_A_FOLD_CROSSES.** THE RING REACHES ACROSS, BUT ONLY A FOLD CROSSES. The converging ring really does get there first, and the gap really is a knob — but it is not the knob that buys an intersection, and no setting of it ever will be.

THE RING IS REAL AND IT GROWS. It thins to 0.156 at the equator and rises to 0.933 at the focus, a factor of 5.97, following the 1/√(sin d) law for a closing ring to a mean ratio of 1.034 (spread 0.071). All of that is before the focal pulse.

AND IT REACHES ACROSS EASILY. The threshold is exactly δ/max|u|, so raising the energy and shrinking the gap both buy it. Shrinking the gap buys something extra — LEAD. At the tightest setting tested the ring spans from just past the equator for the whole converging leg, a lead of 1.41 on the pulse and a dwell of 1.00 of the run. The reaching stops being an instant at the focus and becomes a sustained state.

AND IT STILL NEVER INTERSECTS. Swept over gap and gain with a real segment-intersection test — validated against a limaçon, a lemniscate and a folded loop so it is not a broken counter — the answer is 0, at up to 346 seam crossings. A curve r = f(σ) with f single-valued puts exactly one radius at each direction, so it is embedded by construction and two of its points cannot occupy the same place. This is the winding obstruction seen from the side, and neither knob touches it.

WHAT DOES FOLD IS TANGENTIAL FREEDOM. Let each material point move sideways as well as outward and the map σ₀ ↦ σ can fold. The threshold λε = 0.01210 is predicted from the front's curvature alone and found by bisection to a relative 1.8e-12; it converges under joint refinement of the angular and radial grids (last step 1.5e-04). Past it the curve self-intersects, below it it does not — and it folds first at 0.0000 from the antipode, ON the converging ring, at the moment of tightest focus. The intuition was right about where to look.

FOLDING IS NECESSARY, NOT SUFFICIENT. A crossing always comes with a fold — crossing-without-fold is 0 at every drive tested, which is the embedded-graph theorem again — while fold-without-crossing happens 8 times, because the two branches of a folded map can stay radially apart and never meet in the plane.

AND THE TWO KNOBS ARE ORTHOGONAL. The fold threshold does not know the gap exists: spread 0.0e+00 across a fivefold range of δ, while the spanning threshold scales with δ directly. What the fold threshold does scale with is the pulse — λε ≈ 0.374 w², spread 0.036 — so narrow fronts fold sooner because folding is about how sharply the front turns, not how tall it is. Reducing the distance between the shells changes when the wave arrives at them; it changes nothing about whether it can cross itself.

SCOPE. λ is a modelling choice, not derived from the scalar equation, and λ = 0 is exactly the height field of the earlier probes. What is derived is the threshold given λ, its independence from the gap, and its w² scaling.