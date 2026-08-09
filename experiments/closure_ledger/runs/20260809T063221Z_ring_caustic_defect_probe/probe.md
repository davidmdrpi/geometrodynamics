# What a wavefront has to be like to fold (PR #243)

**Run:** 2026-08-09T06:32:21+00:00

**What kind of wavefront can fold at all?** Answered with the focal geometry of wavefronts. The closed forms are exact; the topology is measured *independently* of them, from the front's own area element. *(Flat bulk, fixed classical vacuole; nothing here is dynamical.)*

- **The pulse**: front stays immersed in a flat bulk; fills its own void
- **The ring**: whole ring arrives at its centre at t = ρ, then stays singular on the axis
- **The coincidence**: the ring that focuses on the inner sphere grazes it (core result)
- **The asymmetry**: solid-angle acceptance ~19% vs 100%; rays remain reversible

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
| **ring first caustic at** | **1.0198** |

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the vacuole and the two sources | **PASS** |
| T2 | `T2_point_does_not_fold_in_a_flat_bulk` | a point does not fold — in this flat bulk | **PASS** |
| T3 | `T3_ring_folds_detected_independently` | the ring's fold, detected from the area element alone | **PASS** |
| T4 | `T4_first_caustic_degenerate_then_persists` | degenerate first caustic, singular afterwards | **PASS** |
| T5 | `T5_critical_ring_grazes` | the critical ring grazes the inner sphere | **PASS** |
| T6 | `T6_solid_angle_acceptance_asymmetry` | solid-angle acceptance; rays stay reversible | **PASS** |
| T7 | `T7_scales_with_the_ratio` | a statement about the radii, not this shell | **PASS** |
| T8 | `T8_assessment` | WHOLE_RING_FOCUSES_AT_ONE_POINT | **PASS** |

## A point does not fold here; a circle folds coherently

| | folds? | detected | closed form | error |
|---|---|---:|---:|---:|
| point source (flat bulk) | **no**, scanned to 2.20 | — | — | — |
| ring source | **yes** | 1.019806 | 1.019804 | 2.0e-06 |

> this is a property of the FLAT bulk, not of point sources.  On a closed manifold the same source folds: on S²/S³ a point's front converges on the antipode at t = πR, which throat_wavefront.py measures directly.  Curvature gives a point a focal set; a flat bulk denies it one.

the fold time is found by scanning the front's own area element for a sign change and bisecting — no radius of curvature is consulted, so the comparison with ρ is a test and not a tautology.  The test is relative (min J < −tol·max|J|) because a parametrisation whose area element merely vanishes, as the direction sphere's poles do, has no caustic; and the orientation is referenced at small t because (X_u × X_v)·N carries the handedness of whatever (u,v) ordering a source happens to use.

| ring `θ₀` | `ρ` | detected | error |
|---:|---:|---:|---:|
| 0.30 | 0.372355 | 0.372356 | 7.4e-07 |
| 0.70 | 0.811714 | 0.811716 | 1.6e-06 |
| 1.10 | 1.122921 | 1.122924 | 2.2e-06 |
| 1.40 | 1.241667 | 1.241669 | 2.5e-06 |

At the first caustic the ring's points are equidistant to `4.4e-16` — the whole ring arrives at once (degenerate (whole ring)), against multiplicity 2 just off it. And it does not end there: for `t > ρ` the tube stays singular at two axis points, separating as `√(t²−ρ²)` — measured `1.998399` against `1.998399` at `t = 1.4ρ`.

## The critical ring grazes the inner sphere

| quantity | value |
|---|---:|
| ring `θ₀` | 54.0342° |
| first caustic at radius | 0.740000 |
| `R_inner` | 0.740000 |
| error | 0.0e+00 |
| launch `sin α` | 0.587302 |
| grazing `sin α` | 0.587302 |
| error | 0.0e+00 |
| ray turning radius | 0.740000 |

## Solid-angle acceptance is asymmetric; propagation is not

| direction | closed form | Monte-Carlo |
|---|---:|---:|
| outer → inner | 0.1906 | 0.1927 |
| inner → outer | 1.0000 | 1.0000 |

A **5.2×** difference in accepted solid angle. Rays reversible: **True**.

> this is an ANGULAR (solid-angle) acceptance asymmetry, NOT nonreciprocal propagation.  Every individual ray is exactly reversible — b is unchanged under reversal, and the accepted inward rays all have b ≤ R_in, so each one climbs back out along its own reverse.  What differs is the measure of launch directions that connect, because a hemisphere at R_out and a hemisphere at R_in are different sets of directions.  No symmetry of the sphere is broken; the ordering of the two radii is the whole of it.

## It scales with the ratio, not the shell

| `R_in` | `sin α_crit` | `θ₀` | `ρ` | caustic at | pulse at | inward accept |
|---:|---:|---:|---:|---:|---:|---:|
| 0.20 | 0.1587 | 80.9° | 1.2440 | 1.2440 | 1.0600 | 0.0127 |
| 0.45 | 0.3571 | 69.1° | 1.1769 | 1.1769 | 0.8100 | 0.0660 |
| 0.74 | 0.5873 | 54.0° | 1.0198 | 1.0198 | 0.5200 | 0.1906 |
| 0.95 | 0.7540 | 41.1° | 0.8277 | 0.8277 | 0.3100 | 0.3431 |
| 1.15 | 0.9127 | 24.1° | 0.5149 | 0.5149 | 0.1100 | 0.5914 |

## Verdict

**WHOLE_RING_FOCUSES_AT_ONE_POINT.** THE WHOLE RING FOCUSES AT ONE POINT. Not every source folds the same way, and the circle's distinction is coherence rather than the mere existence of a caustic — any curved extended source has a focal set. What the circle does is arrive all at once.

In this FLAT bulk a point source does not fold at all: its front is the metric sphere, whose signed area element stays positive, so it crosses at t = ΔR = 0.520 and sweeps the filled ball behind it. That is a property of the flat bulk and not of point sources — on a closed manifold the same source converges on the antipode at t = πR, which the companion wave study measures directly.

A circle's offset tube has area element t(ρ + t cos v), and scanning that element for a sign change — with no radius of curvature consulted — locates the fold at 1.019806 against ρ = 1.019804, an error of 2.0e-06, and reproduces ρ on four unrelated rings to the same precision. At that first caustic the ring's points are equidistant from the centre to 4.4e-16: infinitely degenerate, the whole ring at once. And it does not end there — for t > ρ the tube stays singular at two axis points that separate as √(t²−ρ²), measured 1.998399 against 1.998399 at t = 1.4ρ.

THE CORE RESULT. Requiring that first caustic to land on the inner sphere gives θ₀ = 54.03°, and that same ring launches at sin α = 0.587302 = R_in/R_out — the grazing ray, tangent to the inner sphere (turning radius 0.7400 against R_in = 0.7400). The ring that focuses on the throat and the ray that grazes it are the same ring, and it forms at t = √(R_out² − R_in²) = 1.0198.

ACCEPTANCE, NOT NONRECIPROCITY. Outer to inner, 19.1% of the hemisphere arrives (Monte-Carlo 19.3%); inner to outer, all of it. That is a 5.2× difference in the *measure of launch directions that connect*, not in propagation: every ray is reversible and the probe verifies it (True). No symmetry is broken.

SCOPE. Ray and front geometry in a flat bulk: exact, and independent of any wave solve. It says which fronts CAN fold, not what happens when one does — showing a throat actually form needs backreaction, which this programme does not have.
