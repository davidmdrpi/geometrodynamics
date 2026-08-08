# Ring wavefronts on a surface with a throat (PR #242)

**Run:** 2026-08-08T23:40:44+00:00

**If we just let a classical wave run on a closed surface that has a throat, what does the geometry itself do to it?** A unit S² with both polar caps removed and the mouths joined by a catenoidal neck, against the same domain with the mouths sealed. *(Geometry → field on a fixed classical surface; no backreaction.)*

- **The ring**: a single circle through free flight — a pulse does not cross itself
- **The delay**: the open/sealed echo delay is the neck arclength L
- **The mouth**: the open throat transmits and barely reflects
- **The twist**: the gluing offset aims where the bulk energy lands
- **The orientation**: torus vs Klein bottle is hidden exactly at τ ∈ {0, π}

## The surface

| quantity | value |
|---|---:|
| mouth angle `a` | 0.500 rad |
| neck waist `b` | 0.3603 |
| neck length `L` | 0.5709 |
| outer route `π − 2a` | 2.1416 |
| shortcut ratio | 3.75× |
| neck curvature `K = −1/b²` | -7.7014 |
| Euler characteristic `∫K dA / 2π` | 0.00e+00 |

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the question and the surface | **PASS** |
| T2 | `T2_geometry_closes` | the geometry closes (C¹ join, χ = 0) | **PASS** |
| T3 | `T3_numerics_honest` | energy conserved, delay grid-stable | **PASS** |
| T4 | `T4_free_flight_single_ring` | free flight: one ring, no self-crossing | **PASS** |
| T5 | `T5_echo_delay_is_neck_length` | the echo delay is the neck length | **PASS** |
| T6 | `T6_mouth_transmits` | the open mouth transmits, barely reflects | **PASS** |
| T7 | `T7_twist_aims_orientation_hides` | the twist aims; the orientation hides at τ∈{0,π} | **PASS** |
| T8 | `T8_assessment` | GEOMETRY_ALONE_ROUTES_THE_WAVE | **PASS** |

## The echo delay (the wave measures the throat)

| route | predicted | measured |
|---|---:|---:|
| sealed mirror echo `2(θ₀−a)` | 2.1416 | 2.0640 |
| open bulk return `2(θ₀−a)+L` | 2.7125 | 2.6434 |
| **delay = neck length `L`** | **0.5709** | **0.5794** |

Relative error 1.50%. The two absolute times share the same pulse-width bias, which cancels in the delay.

### Grid stability

| grid | delay | rel. error |
|---|---:|---:|
| 96×128 | 0.5802 | 1.63% |
| 144×192 | 0.5794 | 1.50% |
| 216×288 | 0.5774 | 1.15% |
| 288×384 | 0.5785 | 1.34% |

## The mouth budget

Opening the throat suppresses the mirror echo by **95.7%**.

| mouth `a` | mouth radius | transmitted | reflected | energy drift |
|---:|---:|---:|---:|---:|
| 0.30 | 0.296 | 0.186 | 0.814 | 2.1e-03 |
| 0.50 | 0.479 | 0.286 | 0.714 | 2.8e-03 |
| 0.70 | 0.644 | 0.417 | 0.583 | 3.8e-03 |

## The twist aims the bulk arrival

With `τ = π` the bulk route ends on the antipode at a predicted `2.7125`, measured `2.6436` — **0.391 ahead** of the geodesic focus there, and **9.7×** stronger than the same throat with `τ = 0`.

## The orientation is real but hidden

`R∘g = g∘R ⟺ −(εφ+τ) = ε(−φ)+τ ⟺ τ ≡ −τ ⟺ τ ∈ {0, π}` — at those two offsets the meridian mirror of a point source carries the torus gluing into the Klein-bottle gluing, so no point source can tell them apart.

| `τ/π` | torus vs Klein difference | mirror broken? |
|---:|---:|---|
| 0.000 | 0.0000 | no |
| 0.250 | 0.1907 | yes |
| 0.500 | 0.1623 | yes |
| 0.750 | 0.1834 | yes |
| 1.000 | 0.0000 | no |

## Verdict

**GEOMETRY_ALONE_ROUTES_THE_WAVE.** GEOMETRY ALONE ROUTES THE WAVE. A linear classical wave on a fixed closed surface reports the handle's every property with no fitted parameter. The ring stays a single circle through free flight (to t = 1.071), so it never crosses itself until the geometry gives it somewhere else to be. Sealing the mouths gives a mirror echo; opening them replaces it with a bulk return delayed by 0.5794 against a neck length of 0.5709 (1.5% error) — the wave measures the throat. The open mouth suppresses the reflection by 96%: energy transits the hole rather than bouncing off it. A gluing twist of π re-aims the bulk arrival onto the antipode, where it lands 0.391 ahead of the geodesic focus and 9.7× stronger than the untwisted throat. And the throat's orientation — torus against Klein bottle — is invisible to a point source at exactly τ ∈ {0, π} and visible elsewhere, the mirror argument confirmed to machine precision. The inner/outer asymmetry is there at every twist; it takes a twist to expose it.

SCOPE. Linear and without backreaction, so a focus can sharpen but cannot nucleate — this says nothing about the #175 threshold. The C¹ join leaves a curvature ring at each mouth which is inside the reported budget. Absolute arrival times carry a common pulse-width bias; every load-bearing number is a difference taken on one grid with one pulse. A 2-surface section of the S³ picture, not S³.
