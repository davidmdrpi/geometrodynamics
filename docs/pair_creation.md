# Pair creation belongs to the collision topology, not to focusing

`geometrodynamics/waves/pair_creation.py` · renderer
`scripts/geometrodynamics_v52_pair_creation.py` · probe
`experiments/closure_ledger/pair_creation_probe.py`

## 0. What this round corrects

Every earlier wave round in this arc — from PR #242 on — drew **one** wavefront
converging on its antipode and treated the caustic as the interesting event.

That was the wrong object, and the error is not that the drawing was inaccurate.
**A caustic is where the amplitude gets large; pair creation is a threshold on an
invariant.** Those are different quantities. Breit–Wheeler is

```
γ γ → e⁺ e⁻      s = 2 E₁E₂ (1 − cos θ) ≥ (2m)²
```

and the geometry lives in `θ`, not in `E`. Focusing raises `E`. It does not
create a `θ`.

So a caustic is a **venue** — the geometrically natural place to arrange a
collision, which is exactly how Breit–Wheeler has actually been observed — and
never the mechanism on its own. This round puts back what the original `S²`
picture had and the later ones dropped: **two** waveforms, one up and one down,
from two nearby sources.

## 1. Scope, because the headline depends on it

Kinematics on a fixed round sphere, `c = 1`, no metric evolution and no
backreaction.

* the Breit–Wheeler **threshold and cross-section are imported QED**, not
  derived — the cross-section is the textbook closed form, and is checked
  against its known peak rather than trusted;
* treating a wavefront's null rays as photons is a **correspondence**, stated
  rather than justified;
* **no rate is computed.** A rate needs a photon number density, which a
  classical amplitude does not supply;
* calling the crossing chord a **throat through the bulk** is this program's
  reading, and `shells/junction.py` (PR #249) priced that throat — the bill is
  inherited, not paid;
* the orientation labels `η = ±1` are **carried and checked, not produced**.
  Nothing here derives charge from geometry.

## 2. Focusing is neither sufficient nor necessary

Both halves, because the convenient half alone would be a slogan.

| construction | `s` | reaches `(2m)²`? |
| --- | ---: | --- |
| collinear, amplified `×1` | `0` | no |
| collinear, amplified `×10¹²` | `0` | **still no** |
| crossed head-on, `E = 0.5m`, no focus | `1.0` | no |
| crossed head-on, `E = 1.5m`, no focus | `9.0` | **yes** |

**Not sufficient:** collinear momenta have `s = 0` identically. Amplify by
`10¹²` and it is still exactly zero, so no caustic reaches threshold by being
bright. **Not necessary:** two beams crossed head-on with no focusing anywhere
clear threshold as soon as `E ≥ m`.

**The honest complication, recorded rather than buried.** A spherically
converging front *does* contain diametrically opposed rays, so its
self-invariant is not zero — a converging shell is in that sense self-colliding.
The real distinction is therefore **independence of the sources**, not
brightness, and claiming a single front has `s = 0` everywhere would be the easy
overclaim. What a single front cannot supply is two *independently propagated*
waves.

## 3. The invariant is an identity of geodesic triangles

Put two sources a geodesic distance `δ` apart and let both fire. Their fronts
are geodesic spheres of radius `t`, which intersect for `δ/2 ≤ t ≤ π − δ/2`, and
at every crossing point

```
1 − cos θ  =  (1 − cos δ) / sin²t        ⟹        s(t) = 4 E₁E₂ sin²(δ/2) / sin²t
```

Verified to `2.0e-14` against a control that **never uses the law of cosines**:
the crossing point is solved as a linear system, and the two propagation
directions are built as great-circle tangents in the embedding.

And it is checked on `S²` **and** `S³`, to the same precision, because that is
the identity's own claim — a geodesic triangle lies in a great 2-sphere whatever
it is embedded in, so the dimension cannot matter.

## 4. Which makes the collision head-on *twice*

`s(t)` is **U-shaped**. It is maximal — `4E₁E₂`, head-on — at *both* ends of the
crossing window, and minimal at the equator:

| `δ` | `θ` at `t = δ/2` | `θ` at equator | `θ` at `t = π − δ/2` | `s` equator / `s` head-on |
| ---: | ---: | ---: | ---: | ---: |
| 0.15 | `180°` | `0.15` rad | `180°` | `0.0056` |
| 0.42 | `180°` | `0.42` rad | `180°` | `0.043` |
| 1.00 | `180°` | `1.00` rad | `180°` | `0.23` |
| 2.00 | `180°` | `2.00` rad | `180°` | `0.71` |

The equator angle is **exactly** the source separation `δ`. So the moment the
wavefronts are largest is the moment the invariant is *smallest* — the opposite
of what the single-caustic pictures encouraged.

The head-on ends are tested on the **invariant**, not the angle: `arccos` near
`−1` has a square-root singularity that costs eight digits and turns an exact
statement into a tolerance argument.

## 5. So the threshold opens two windows, and never one

`s ≥ 4m²` ⟺ `sin t ≤ (E/m) sin(δ/2)`, which holds near *both* ends and fails in
the middle. Confirmed against a 40,000-point scan of `s(t)`:

| `E/m` | windows | where |
| ---: | ---: | --- |
| 0.6 | **0** | even head-on, `s_max = 4E² < 4m²` |
| 1.0 | **0** | zero width — threshold touched only at the two head-on instants |
| 1.4 | **2** | `[0.210, 0.296]` and `[2.845, 2.932]` |
| 3.0 | **2** | `[0.210, 0.676]` and `[2.466, 2.932]` |
| 6.0 | **1** | merged, above `E = m/sin(δ/2) = 4.797 m` |

## 6. And only the far window is a collision of independent waves

This is the answer to *why the second interaction has to be antipodal*, and it
is measured as a path length rather than asserted.

| | fronts have travelled | against a separation of |
| --- | ---: | ---: |
| near window closes at | `0.296` | `δ = 0.42` |
| far window opens at | `2.845` | — |

A factor of **9.6**. The near window sits on top of the sources: the fronts have
travelled less than the separation itself, so nothing there has propagated
independently of anything — it is the emission region. The far window is reached
only after each front has crossed a **half-circumference**.

That is the geometric content of the correction. Pair creation needs two
independently propagating waves to collide; on a round space the only place two
fronts from nearby sources meet head-on again, having genuinely separated
first, is **at the antipodes**.

## 7. One further trap, caught

The momenta are exact — perpendicular to their own wavefront to `1e-15`, with
the angle between them matching the closed form to `2e-13`. But a **figure**
shows their *projection*, and projection does not preserve angles.

| | worst error |
| --- | ---: |
| projected opening angle vs true | **67.4°** |
| disagreement *between the two crossing points* | **56.4°** |

The two crossing points have, exactly, the same opening angle. A reader
measuring the drawn arrows would conclude they differ by up to 56°. So the
renderer draws the opening angle in **the plane the two momenta actually span**,
where it is undistorted, and the arrows on the sphere are labelled as
projections.

This is the same trap as the previous rounds, and it keeps arriving in a new
costume: a plotting artefact that reads as physics.

## 8. What the cross-section adds, and what it does not

The textbook Breit–Wheeler `σ(s)`, verified against its known peak:

| | value | textbook |
| --- | ---: | ---: |
| `β` at peak | `0.7013` | `0.701` |
| `√s / 2m` at peak | `1.4028` | `1.40` |
| `σ_peak / σ_T` | `0.2556` | `0.256` |

It vanishes at threshold and **falls** at large `s`, which matters for the
picture: the most violent part of the crossing is not the most productive. But
`σ` is a cross-section, not a rate — turning it into events needs a photon
number density this round does not have.

## 9. What this closes, and what it does not

**Closes:** that a caustic is not a creation event; that the threshold is a
condition on the opening angle, which focusing cannot supply; that two sources
give a U-shaped invariant with two disjoint threshold windows; and that only the
antipodal window is a collision of independently propagated waves.

**Does not close:** whether the pair is *produced*. No rate, no field
quantisation, no backreaction, and no derivation of charge from orientation. The
throat interpretation of the crossing chord remains this program's reading, with
its exotic-matter bill inherited from PR #249.

## 10. The lesson

The arc's recurring finding, one level up. Previous rounds caught captions that
disagreed with the geometry, and objects drawn correctly that were the wrong
object. This one is neither: the drawings were accurate and the object was
right. **The quantity was wrong.** Amplitude is not an invariant, and no amount
of care about *how* to draw a caustic would have surfaced that — only asking
what the threshold is a threshold *on*.
