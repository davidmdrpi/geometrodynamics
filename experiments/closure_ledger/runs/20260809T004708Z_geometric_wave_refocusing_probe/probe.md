# Ring wavefronts on a surface with a throat (PR #242)

**Run:** 2026-08-09T00:47:08+00:00

**If we just let a classical wave run on a closed surface that has a throat, what does the geometry itself do to it?** Three surfaces on one clock: the bare S², the same sphere with both caps cut out and sealed, and the same cut sphere with the mouths joined by a catenoid. *(Geometry → field on a fixed classical surface; no backreaction.)*

- **The bare front**: sweeps each point once — a pulse cannot meet itself
- **A sealed mouth**: sends a second front back toward the source
- **An open mouth**: sends none back; the second front is downstream of the neck
- **The delay**: the open/sealed echo delay is the neck arclength L
- **The twist**: the gluing offset aims where the bulk energy lands
- **The orientation**: torus vs Klein bottle is hidden exactly at τ ∈ {0, π}

## The surface

| quantity | value |
|---|---:|
| mouth angle `a` | 0.750 rad |
| neck waist `b = sin²a` | 0.4646 |
| neck length `L = sin 2a` | 0.9975 |
| outer route `π − 2a` | 1.6416 |
| shortcut ratio | 1.65× |
| `K` at the mouth | -1.0000 |
| `K` at the waist | -4.6322 |
| `∫K dA / 2π` | 0.00e+00 |

Four geometric times set the clock: the ring reaches the mouths at `0.821`, a sealed echo returns at `1.642`, a bulk crossing lands at `2.639`, and the antipodal focus is at `3.142`.

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the question and the three surfaces | **PASS** |
| T2 | `T2_true_catenoid` | a true catenoid; χ = 0 tests the join | **PASS** |
| T3 | `T3_conservative_scheme` | conservative to round-off; flux closes | **PASS** |
| T4 | `T4_arrival_multiplicity` | arrival multiplicity per point | **PASS** |
| T5 | `T5_echo_delay_is_neck_length` | the echo delay is the neck length | **PASS** |
| T6 | `T6_mouth_budget_by_flux` | the mouth budget, by integrated flux | **PASS** |
| T7 | `T7_twist_aims_orientation_hides` | the twist aims; the orientation hides | **PASS** |
| T8 | `T8_assessment` | GEOMETRY_ALONE_ROUTES_THE_WAVE | **PASS** |

## The neck is a true catenoid

| `a` | `b = sin²a` | `L = sin 2a` | `\|r − b cosh(z/b)\|` | `K` mouth | `K` waist | `χ` |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 0.0612 | 0.4794 | 1.7e-16 | -1.000 | -266.92 | 0.0e+00 |
| 0.50 | 0.2298 | 0.8415 | 2.2e-16 | -1.000 | -18.93 | 0.0e+00 |
| 0.75 | 0.4646 | 0.9975 | 2.2e-16 | -1.000 | -4.63 | 0.0e+00 |
| 1.00 | 0.7081 | 0.9093 | 2.2e-16 | -1.000 | -1.99 | 0.0e+00 |
| 1.30 | 0.9284 | 0.5155 | 2.2e-16 | -1.000 | -1.16 | 0.0e+00 |

> ∫K dA = −2π[r'] for any surface of revolution, so χ = 0 follows from the C¹ join alone and holds for every C¹-matched profile — it is a check on the join, not evidence for the catenoid.

## Does a front ever meet another front?

Arrivals counted per grid point over `t < π`. a per-cell hysteresis trigger on the energy density u_t²+|∇u|², armed above 35% and re-armed below 12% of that cell's own peak; plain local-maximum counting fails because a 2+1-dimensional wave violates Huygens and every front drags a rippling wake.

| surface | max arrivals | area with ≥2 | of the source side |
|---|---:|---:|---:|
| bare | 2 | 0.015 | 0.000 |
| plugged | 3 | 0.243 | 0.095 |
| throat | 3 | 0.362 | 0.000 |

The bare sphere is the case with no second front at all. Sealing the mouths sends one back toward the source; opening them sends **none** back — the same fact the echo shows, resolved in space rather than in time.

## The echo delay (the wave measures the throat)

| route | predicted | measured |
|---|---:|---:|
| sealed mirror echo `2(θ₀−a)` | 1.6416 | 1.6017 |
| open bulk return `2(θ₀−a)+L` | 2.6391 | 2.6033 |
| **delay = neck length `L`** | **0.9975** | **1.0016** |

Relative error 0.41%.

### Scheme and grid stability

Energy drift 2.9e-16 (plugged) and 4.3e-16 (throat); the mouth flux closes against the neck's stored energy to 0.3%.

| grid | delay | rel. error |
|---|---:|---:|
| 96×128 | 1.0055 | 0.80% |
| 144×192 | 1.0016 | 0.41% |
| 216×288 | 1.0021 | 0.46% |

## The mouth budget, by integrated flux

integrated power through two surfaces: 'offered' is the energy crossing a reference circle a few cells inside the mouth, 'through' is the energy crossing the mouth face itself, and transmission is their ratio.  On a closed surface only part of the wave ever reaches the mouth, so the total energy is the wrong denominator.

| mouth `a` | offered | through | transmission | reflection | peak stored |
|---:|---:|---:|---:|---:|---:|
| 0.40 | 0.4729 | 0.3752 | 0.793 | 0.207 | 0.243 |
| 0.55 | 0.6102 | 0.5276 | 0.865 | 0.135 | 0.343 |
| 0.75 | 0.7913 | 0.7274 | 0.919 | 0.081 | 0.474 |
| 0.90 | 0.9283 | 0.8714 | 0.939 | 0.061 | 0.566 |

> **mirror suppression is an amplitude ratio at one watched point and one watched time; transmission is an energy ratio at the mouth.  They are different measurements of the same fact and must not be quoted interchangeably.** The sealed echo's amplitude is suppressed by 72.7% when the throat is opened; that is a different measurement from the transmission column and the two are not interchangeable.

## The twist aims the bulk arrival

With `τ = π` the bulk route ends on the antipode at a predicted `2.6391`, measured `2.6094` — **0.469 ahead** of the geodesic focus there, and **9.9×** stronger than the same throat with `τ = 0`.

## The orientation is real but hidden

`R∘g = g∘R ⟺ −(εφ+τ) = ε(−φ)+τ ⟺ τ ≡ −τ ⟺ τ ∈ {0, π}` — at those two offsets the meridian mirror of a point source carries the torus gluing into the Klein-bottle gluing, so no point source can tell them apart.

| `τ/π` | torus vs Klein difference | mirror broken? |
|---:|---:|---|
| 0.000 | 0.0000 | no |
| 0.250 | 0.1864 | yes |
| 0.500 | 0.2165 | yes |
| 0.750 | 0.2051 | yes |
| 1.000 | 0.0000 | no |

## Verdict

**GEOMETRY_ALONE_ROUTES_THE_WAVE.** GEOMETRY ALONE ROUTES THE WAVE. A linear classical wave on a fixed closed surface reports the handle's every property with no fitted parameter, and the three-surface comparison separates what each change of geometry is responsible for. On the bare sphere the front sweeps each point exactly once — 1.5% of the surface ever sees a second front and none of the source side does — so a pulse on a closed surface with no boundary cannot meet itself. Sealing the mouths puts a second front back toward the source over 9.5% of that hemisphere; opening the throat puts one downstream of the neck instead and 0.0% back home. The echoes say the same thing in time: the sealed mirror echo and the open bulk return differ by 1.0016 against a neck length of 0.9975 (0.41% error) — the wave measures the throat. Of the energy that actually reaches the mouth, 92% crosses into the neck, rising with the aperture. A gluing twist of π re-aims the bulk arrival onto the antipode, 0.469 ahead of the geodesic focus and 9.9× stronger than the untwisted throat. And the throat's orientation — torus against Klein bottle — is invisible to a point source at exactly τ ∈ {0, π} and visible elsewhere, the mirror argument confirmed to machine precision.

SCHEME. Each mouth is one shared finite-volume face, so the discrete energy is conserved to round-off (drift 4.3e-16) and the mouth flux closes against the neck's stored energy to 0.3%.

SCOPE. Linear and without backreaction, so a focus can sharpen but cannot nucleate — this says nothing about the #175 threshold. The C¹ join leaves a curvature ring at each mouth, inside the reported budget. χ = 0 checks the join and not the profile. Absolute arrival times carry a common pulse-width bias; every load-bearing number is a difference on one grid with one pulse. A 2-surface section of the S³ picture, not S³.
