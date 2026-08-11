# The scaling at the seam is a choice

_2026-08-11T05:21:37+00:00_

> Both gluings are **representation** choices, not derived boundary conditions.

## The emerging feature

| rule | height in | height out | arc out | aspect distortion |
|---|---:|---:|---:|---:|
| `translate` | 0.1000 | 0.1000 | 0.0370 | **1.7027** |
| `conformal` | 0.1000 | 0.0587 | 0.0370 | **1.0000** |

The circumference ratio is 1.7027, and the translate rule pays it in full.

## The inward sheets

| rule | sheet −1 | −2 | −3 | −4 |
|---|---:|---:|---:|---:|
| `translate` | 0.220 | -0.300 | -0.820 | -1.340 |
| `conformal` | 0.435 | 0.255 | 0.150 | 0.088 |

Arithmetic sheets walk through `r = 0`; geometric ones accumulate at it.

## What a winding number would look like

| seam | magnifications from r₀ = 0.80, 1.00, 1.20 | spread |
|---|---|---:|
| `translate` | 1.6500, 1.5200, 1.4333 | `2.2e-01` |
| `conformal` | 1.7027, 1.7027, 1.7027 | `2.2e-16` |

A scale is start-independent; a shift is not.

## ...and the winding is still zero

| seam | radial law | unsigned | signed | winding |
|---|---|---:|---:|---:|
| `translate` | `additive` | 274 | +0 | +0 |
| `conformal` | `multiplicative` | 268 | +0 | +0 |
| `conformal` | `additive` | 36 | +0 | +0 |

## Verdict

**6/6 checks passed.**

**THE_GLUING_SETS_THE_SCALE_BUT_NOT_THE_TOPOLOGY.** THE GLUING SETS THE SCALE, BUT NOT THE TOPOLOGY. How the seam is glued was a free choice, and it decides what the emerging wave looks like — but not whether it can carry a winding number.

THE TRANSLATE RULE RETURNS A CARICATURE. Carrying a radial offset across unchanged lands the feature on an arc shorter by 1.7027, so its aspect ratio is multiplied by 1.7027 — the emerging wave is not the same wave. The conformal rule scales the offset with the boundary it crosses, so height and arc length shrink together and the distortion is 1.0000: a faithful scaled copy.

AND THE TRANSLATE RULE RUNS OFF THE ORIGIN. Its inward sheets are arithmetic, so they march 0.22 → -0.30 → -0.82: straight through r = 0 into negative radius, because subtracting a fixed gap has nothing to stop it. Conformal sheets are geometric — 0.435 → 0.255 → 0.150 — accumulating at the origin and never arriving, and they pair with a multiplicative radial law that is positive by construction.

THE CHOICE ALSO DECIDES WHAT WINDING WOULD LOOK LIKE. A curve that genuinely winds is a logarithmic spiral on the conformal seam: it returns to the same point of the quotient magnified by 1.7027, and by the same factor from every starting radius (spread 2.2e-16). On the translate seam the same curve returns displaced instead, with a ratio that drifts by 0.217 depending on where it started — not a scale at all. A conformal gluing makes topological charge observable as a magnification; a translate gluing hides it.

BUT THE WINDING NUMBER IS STILL ZERO. Rebuilt on a conformal seam with a multiplicative radial law — a different identification, a different sheet structure, a different notion of size — the winding is identically zero at every gain tested, driving up to 274 unsigned crossings. ρ(σ) comes from a single-valued function on the circle whichever coordinate the seam translates in, so its degree is zero either way. The earlier negative result was not an artefact of an arbitrary scaling.

THE TWO RULES ARE THE SAME RULE UNTIL THE WAVE REACHES THE SEAM. Their difference is proportional to the excursion, with slope 0.4127 — so the earlier pictures were not wrong where the wave was small. They part company exactly where it matters.

SCOPE. The conformal rule is preferred on grounds of consistency — it returns a scaled copy and cannot produce a negative radius — not because any dynamics picks it out. Neither rule makes the seam something the wave can feel, and neither gives a height field anywhere to wind.