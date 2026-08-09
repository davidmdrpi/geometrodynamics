# Why a throat needs a ring: front topology across the bulk (PR #242)

**Run:** 2026-08-09T05:45:27+00:00

**What kind of wavefront can make a throat?** Answered with the focal-set geometry of wavefronts, in closed form, checked against an independent numerical count — no wave solve. *(Fixed classical vacuole; nothing here is dynamical.)*

- **The pulse**: front is an embedded sphere; fills its own void
- **The ring**: front folds onto its centre at t = ρ — a degenerate caustic
- **The coincidence**: the ring that focuses on the inner sphere grazes it
- **The asymmetry**: outer→inner ~19% accepted, inner→outer 100%

## The vacuole

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

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the vacuole and the two sources | **PASS** |
| T2 | `T2_point_never_folds` | a point's focal set is empty — never folds | **PASS** |
| T3 | `T3_ring_must_fold` | a circle's focal set is one point — folds at ρ | **PASS** |
| T4 | `T4_defect_is_degenerate` | the defect is degenerate: the whole ring at once | **PASS** |
| T5 | `T5_critical_ring_grazes` | the critical ring grazes the inner sphere | **PASS** |
| T6 | `T6_acceptance_asymmetry` | acceptance asymmetry, closed form vs Monte-Carlo | **PASS** |
| T7 | `T7_scales_with_the_ratio` | a statement about the radii, not this shell | **PASS** |
| T8 | `T8_assessment` | ONLY_A_RING_FOLDS | **PASS** |

## A pulse cannot; a ring must

| | focal set | self-intersects | when |
|---|---:|---|---:|
| point source | 0 points | **never** | — |
| ring source | 1 point | **yes** | `t = ρ = 1.0198` |

At the focus every point of the ring is equidistant to `4.4e-16` — the whole ring arrives at once (degenerate (whole ring)), against multiplicity 2 just off it. That is the defect: a codimension-2 point where the front stops being embedded, not a smooth focus.

## The critical ring grazes the inner sphere

| quantity | value |
|---|---:|
| ring `θ₀` | 54.0342° |
| defect forms at radius | 0.740000 |
| `R_inner` | 0.740000 |
| error | 0.0e+00 |
| launch `sin α` | 0.587302 |
| grazing `sin α` | 0.587302 |
| error | 0.0e+00 |
| ray turning radius | 0.740000 |

## The bulk is not symmetric between its faces

| direction | closed form | Monte-Carlo |
|---|---:|---:|
| outer → inner | 0.1906 | 0.1927 |
| inner → outer | 1.0000 | 1.0000 |

A **5.2×** asymmetry. Not a broken symmetry of the sphere — the ordering of the two radii is the asymmetry.

## It scales with the ratio, not the shell

| `R_in` | `sin α_crit` | `θ₀` | `ρ` | defect at | pulse at | inward accept |
|---:|---:|---:|---:|---:|---:|---:|
| 0.20 | 0.1587 | 80.9° | 1.2440 | 1.2440 | 1.0600 | 0.0127 |
| 0.45 | 0.3571 | 69.1° | 1.1769 | 1.1769 | 0.8100 | 0.0660 |
| 0.74 | 0.5873 | 54.0° | 1.0198 | 1.0198 | 0.5200 | 0.1906 |
| 0.95 | 0.7540 | 41.1° | 0.8277 | 0.8277 | 0.3100 | 0.3431 |
| 1.15 | 0.9127 | 24.1° | 0.5149 | 0.5149 | 0.1100 | 0.5914 |

## Verdict

**ONLY_A_RING_FOLDS.** ONLY A RING FOLDS. The difference between a pulse and a ring is not a matter of degree, it is the focal set. A point has none, so its front is the metric sphere — embedded at every t, never touching itself, sweeping the filled ball behind it. It crosses the bulk at t = ΔR = 0.520 and does nothing else; it fills its own void. A circle's focal set is its centres of curvature, and because every point of a circle shares the same centre that set collapses to a single point which the entire ring reaches simultaneously at t = ρ. Measured, the ring's points are equidistant from it to 4.4e-16: a degenerate caustic of infinite multiplicity, where the front stops being embedded. That is a codimension-2 defect made by geometry alone.

The two conditions coincide. Requiring the defect to land on the inner sphere gives θ₀ = 54.03°, and that same ring launches at sin α = 0.587302 = R_in/R_out — the grazing ray, tangent to the inner sphere (turning radius 0.7400 against R_in = 0.7400). The ring that focuses on the throat and the ray that grazes it are the same ring.

And the bulk is not symmetric between its two faces. Outer to inner, only 19.1% of the hemisphere arrives (Monte-Carlo 19.3%); inner to outer, all of it does — a 5.2× asymmetry from the ordering of two radii, with no symmetry broken anywhere.

SCOPE. Ray and front geometry in a flat bulk: exact, and independent of any wave solve. It says which fronts CAN fold, not what happens when one does — showing a throat actually form needs backreaction, which this programme does not yet have.
