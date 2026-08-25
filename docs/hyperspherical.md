# What higher dimension does to the bulk picture

**Module:** `geometrodynamics/viz/hyperspherical.py`
**Probe:** `python -m experiments.closure_ledger.hyperspherical_probe` (11/11)
**Tests:** `tests/test_viz_hyperspherical.py` (58)
**Renderer:** `scripts/geometrodynamics_v70_hyperspherical.py`

---

Three-dimensional intuition models an extra dimension as *the same sphere with
another direction available*. That is wrong in ways that matter for this arc,
and this round measures them rather than asserting them.

The one that matters most is [§5](#5-the-collapse-is-fⁿ-not-f): a packet can
keep **exactly the same angular footprint** while its physical overlap shrinks
as `f₀ⁿ`. PR #268's review had already forced the overlap *size* law to be
scoped by dimension; this follows that scoping to its conclusion, and it
changes what the finite-bearing picture is a picture *of*.

## 0 · a correction first, because it is the usual one

The familiar "spheres peak at 5D" is about the **volume of the unit ball**,

    V_d = π^{d/2}/Γ(d/2+1) ,        peaking at  d = 5 ,

not the surface measure of the unit sphere,

    A_{d−1} = 2π^{d/2}/Γ(d/2) ,     peaking at  d = 7 — the sphere being S⁶ ⊂ ℝ⁷ .

| `d` | `V_d` (ball) | `A_{d−1}` (sphere) | |
|--|--|--|--|
| 4 | `4.93480` | `19.73921` | |
| **5** | **`5.26379`** | `26.31895` | ← ball peaks |
| 6 | `5.16771` | `31.00628` | |
| **7** | `4.72477` | **`33.07336`** | ← sphere peaks, `S⁶` |
| 8 | `4.05871` | `32.46970` | |

**And neither is a fact about dimension alone.** `V_d` and `A_{d−1}` carry
different physical units at different `d`, so comparing them across `d`
implicitly picks a length scale — and the peak follows it:

| `R` | ball peaks at | sphere peaks at |
|--|--|--|
| `0.5` | `d = 1` | `d = 2` *(tied with 3)* |
| `1.0` | `d = 5` | `d = 7` |
| `2.0` | `d = 24` | `d = 26` |
| `4.0` | `d = 100` | `d = 102` |

So "the peak at 5" is a statement about the **unit** ball and nothing more. Two
notes on the scan itself, both from mistakes it caught: the first version
searched `d < 80` and reported `79` for `R = 4`, which is a search range and not
a peak; and at `R = ½` the sphere measure is *exactly* tied between `d = 2` and
`d = 3` (`2πR = 4πR²` there), which a bare `max` resolves by floating-point
accident. Ties are now reported, not resolved.

The interesting geometry is elsewhere.

## 1 · almost all of a sphere is at the equator of any chosen point

The shell at geodesic angle `χ` from a pole has measure `∝ sin^{n−1}χ dχ`.
Writing `χ = π/2 + δ` gives `sin^{n−1}χ ≈ e^{−(n−1)δ²/2}`, so the occupied band
has width `Δχ ~ 1/√n`. Measured against the exact integral:

| `n` | `std(χ)` | `std·√n` | mass within `1/√n` of the equator |
|--|--|--|--|
| 2 | `0.683667` | `0.966852` | `0.6496` |
| 3 | `0.567862` | `0.983566` | `0.6587` |
| 15 | `0.258009` | `0.999264` | `0.6774` |
| 200 | `0.070710` | `0.999996` | `0.6823` |
| 1000 | `0.031623` | `1.000000` | `0.6826` |

> **In high dimension almost every point is nearly 90° from any chosen point** —
> not "somewhere randomly around the sphere", but concentrated in a band that
> shrinks like `1/√n`.

## 2 · so the antipode is an extremely non-generic relation

For random unit vectors `x·y` has mean `0` and width `1/√n`, so
`α = arccos(x·y) → π/2`. Random points do **not** tend toward `α = π`:

| ambient `n` | `std(x·y)·√n` | mean `α/π` | fraction with `α > 0.99π` |
|--|--|--|--|
| 3 | `0.9987` | `0.5000` | `3.2e-04` |
| 4 | `0.9991` | `0.5000` | `5.0e-06` |
| 10 | `1.0034` | `0.5000` | `0` |
| 500 | `1.0014` | `0.5000` | `0` |

That sharpens what this repository's antipodal identification is claiming.
Selecting `x ↔ −x` is not "pairing a point with the far one" — the far ones are
a vanishing fraction of an enormous nearly-orthogonal majority, and the fraction
*shrinks as dimension rises*. The identification picks a measure-zero relation
out of an overwhelmingly generic alternative.

**What that does not do** is make the identification correct. It removes one
reading of it — the bland one — and nothing more.

## 3 · the wavefront is a point, then an equatorial shell, then a point

A wave launched at a pole of `S^n_R` has transverse measure

    𝒜(χ) = |S^{n−1}| · (R sin χ)^{n−1} ,

zero at the pole, maximal at `χ = π/2`, zero again at the antipode. Raising `n`
makes the middle dominate more — the share of the path-integrated measure
within `0.1` rad of the equator runs `10.0%` on `S²` to `20.2%` on `S⁷` — while
the poles get relatively tinier.

So the dimensional-descent picture **strengthens** with dimension, and with it
the contrast the model turns on: almost all of the propagation measure sits
around the equator, yet the whole radial congruence still reconverges at one
antipodal locus. For `S³`, the repo's case, the front goes as `sin²χ`.

## 4 · the embedding centre is a direction space collapsed to a point

An intrinsic `S^n` has no centre. Embedded in `ℝ^{n+1}` it has exactly one, and
that point is the unique fixed point of the whole rotation group `SO(n+1)`.
Every direction `n̂ ∈ S^n` runs a radial line through it, so the origin is where
an entire `S^n`'s worth of directions coincide.

**Which is a reading of what PR #268 did.** Replacing `f = 0` with `f₀ > 0` does
not merely avoid a singularity — it restores, at finite size, the direction
space the origin had crushed. The bearing *is* the blown-up direction space, and
its measure is `f₀ⁿ·A_n` where `A_n = |S^n|`:

| `n` | `f₀` | `A_n` | bearing measure `f₀ⁿ·A_n` |
|--|--|--|--|
| 1 | `1e-03` | `6.2832` | `6.283e-03` |
| 2 | `1e-03` | `12.5664` | `1.257e-05` |
| 3 | `1e-03` | `19.7392` | `1.974e-08` |

Nonzero in every case, which is the point — and vanishingly small in a way the
2-D drawing cannot suggest.

### so the centre is a routing manifold, not a hub

That reframes what the singular origin *is*. It is **not** obtained because
every direction becomes equivalent there. It is obtained because an entire
finite direction space is compressed to **zero proper measure while its angular
structure stays intact** — a different statement, and the one the measurements
support.

Directional capacity — how many distinguishable directions a bearing carries at
angular resolution `ε` — is **dimensionless**. `f₀` does not enter it at all:

| bearing `S^n` | capacity at `60°` | capacity at `20°` | proper measure at `f₀ = 1e-03` |
|--|--|--|--|
| `S¹` | `3.0` | `9.0` | `6.283e-03` |
| `S²` | `4.0` | `33.2` | `1.257e-05` |
| **`S³`** | **`5.1`** | **`113.5`** | `1.974e-08` |
| `S⁴` | `6.4` | `374.1` | `2.632e-11` |
| `S¹⁰` | `20.4` | `3.524e+05` | `2.073e-29` |
| `S²⁰` | `112.3` | `2.236e+10` | `2.929e-61` |

Scanning `f₀` at fixed `n = 3` makes the independence explicit — the capacity is
the *same float* while the proper measure runs off the bottom:

| `f₀` | proper measure `f₀³·A₃` | capacity at `20°` |
|--|--|--|
| `1e-01` | `1.974e-02` | `113.5295` |
| `1e-04` | `1.974e-11` | `113.5295` |
| `1e-07` | `1.974e-20` | `113.5295` |

And the routing stays usable in high dimension for the reason §2 gave: in `ℝ¹⁰⁰⁰`,
`1024` random directions are *all* pairwise within `0.145` of orthogonal, so a
compressed direction space is not a crowded one.

### what the `f₀ → 0` limit actually separates

Three things were being run together. Taking two sectors at fixed angular
separation `0.3` rad and shrinking the neck:

| `f₀` | they meet? | overlap angle | proper overlap `n=2` | proper overlap `n=3` |
|--|--|--|--|--|
| `1e-01` | yes | `0.225` | `5.06e-04` | `1.14e-05` |
| `1e-03` | yes | `0.225` | `5.06e-08` | `1.14e-11` |
| `1e-06` | yes | `0.225` | `5.06e-14` | `1.14e-20` |
| `1e-09` | yes | `0.225` | `5.06e-20` | `1.14e-29` |

- **angular incidence** — survives, unchanged;
- **the overlap verdict** (do they meet at all) — survives, unchanged;
- **the proper interaction measure** — collapses as `f₀ⁿ`.

The first reading treated the limit as *merging the routes*. It merges their
sizes, not their labels.

*§4 is an interpretation, not a result.* PR #268's bearing stands on its own
geometry, and nothing in `regularized_center` depends on this section.

## 5 · the collapse is `fⁿ`, not `f`

For `dℓ² = ds² + f(s)²dΩ_n²` the transverse measure of an angular patch is
`fⁿ dΩ_n`. So squeezing from `F` to `f₀` costs `(f₀/F)ⁿ`:

| squeeze `f₀/F` | `S¹` | `S²` | `S³` | `S⁴` |
|--|--|--|--|--|
| `1e-01` | `1.0e-01` | `1.0e-02` | `1.0e-03` | `1.0e-04` |
| `1e-02` | `1.0e-02` | `1.0e-04` | `1.0e-06` | `1.0e-08` |
| `1e-03` | `1.0e-03` | `1.0e-06` | `1.0e-09` | `1.0e-12` |

> **The angular overlap can stay finite while the physical overlap collapses as
> `f₀ⁿ`.**

### which `n` is physical depends on which object is drawn

This has to stay explicit, or the exponent will migrate between objects that do
not share it. `n` is the dimension of the object's own **transverse sphere**,
which is a fact about the object, not a modelling choice:

| object | transverse sphere | `n` | patch at `f₀/F = 1e-03` | vs the 2-D drawing |
|--|--|--|--|--|
| the drawn 2-D cross-section | `S¹` | 1 | `1.0e-03` | — |
| **PR #265's spatial throat** | **`S²`** | **2** | `1.0e-06` | **`1e+03`** |
| a bearing in a 4-spatial-`d` embedding | `S³` | 3 | `1.0e-09` | `1e+06` |
| a bearing in a 5-spatial-`d` embedding | `S⁴` | 4 | `1.0e-12` | `1e+09` |

So the millionfold figure belongs to the `S³` **bearing**, not to every throat
quantity in the repository. The `#265` throat is `n = 2`: its own understatement
against the drawing is a **thousand**, and `physical_throat`'s neck area is
`4π f₀²` — measured `1.958592e-07` at `f₀ = 1.248e-04`, an `f²` law — which is
what an `S²` cross-section is required to give. `measure_which_n_is_physical_for_which_object`
pins each row.

That changes what the finite-bearing picture is a picture *of*. Not two thick
ribbons physically squeezing until they touch — the faithful reading is

    large angular structure → tiny proper measure → large angular structure

with the angular labels intact throughout. Much less like ribbons squeezing;
much more like a **vast angular configuration space packed into an extremely
small proper region**.

The consequence for the drawing is practical: two fronts with finite angular
overlap `ΔΩ` keep that overlap constant all the way down, while
`V_overlap ∝ f(s)ⁿ`. So the picture does not have to force a dramatic
macroscopic crossing to be honest. It can show an unmistakable intersection in
*angular* coordinates while representing an interaction region that is
extraordinarily small — the intersection is real in the topology of the
coordinate mapping, and tiny in proper measure. Two angular sectors can
become coincident on a greatly compressed direction space at finite `f₀`, then
re-expand into different radial sectors — which a 2-D cross-section can show and
3-D intuition would not suggest.

The yes/no overlap criterion is untouched. It was always angular, and it stays
angular; it is only the *size* of the meeting that carries the exponent. This is
the same `f^q` weighting `regularized_center.monopole_resistance` carries,
followed through to the overlap region instead of to the resistance.

## 6 · orientability flips with parity — and this repo uses two quotients

`ℝP^n` is orientable **iff `n` is odd**: the antipodal map is the restriction of
`−I` on `ℝ^{n+1}`, so the orientation factor is `det(−I_{n+1}) = (−1)^{n+1}`.
`ℝP²` non-orientable, `ℝP³` orientable, `ℝP⁴` non-orientable. A **parity**
effect, not a size effect — nothing "the next sphere up" would suggest.

The consequence here is specific, because this repository carries **two**
antipodal quotients, at *consecutive* `n`:

| spatial `d` | spatial `S^d/± = ℝP^d` | exchange `(ℝ^d∖0)/± ≃ ℝP^{d−1}` | `π₁`(exchange) |
|--|--|--|--|
| 2 | `ℝP²` **non-orientable** | `ℝP¹` orientable | `ℤ` (braid) |
| **3** | **`ℝP³` orientable** | **`ℝP²` non-orientable** | `ℤ₂` |
| 4 | `ℝP⁴` **non-orientable** | `ℝP³` orientable | `ℤ₂` |
| 5 | `ℝP⁵` orientable | `ℝP⁴` **non-orientable** | `ℤ₂` |

At `d = 3` the spatial quotient is orientable and the exchange space is not —
and it is that `ℝP²` whose Pin⁻ structure makes the throat a spinor (#170,
#188). Being one apart, the two are **always of opposite parity**, so raising
the spatial dimension by one *swaps which of them is non-orientable*. At `d = 4`
the exchange space would be the orientable `ℝP³`, and the spin-statistics
mechanism would have to be re-derived rather than carried across. `d = 2` fails
differently again: `ℝP¹ ≃ S¹` has `π₁ = ℤ`, the braid group, and the exchange
is not a `ℤ₂` at all.

This is stated as a consequence **if** the dimension moved. It is not an
argument that it should.

## 7 · and `S³` is not a generic sphere

`S³ ≅ SU(2)`, the unit quaternions. It is **parallelizable** — the frame
`(q·i, q·j, q·k)` is tangent, orthonormal and nowhere zero, verified to
`6.7e-16` across 20 000 sampled points — which only `S¹`, `S³` and `S⁷` are;
`S²` admits no nowhere-zero tangent field at all. And it carries the Hopf
fibration `S¹ → S³ → S²`, checked by confirming the image lands on `S²` to
`7.8e-16` and that a common phase on both complex components moves the point
(median `1.41`) but not its image (`7.8e-16`).

So dimensions do not simply add room. Some are algebraically special, and the
one this repository is built on is one of them — which is a standing argument
against reading any of §1–§6 as a smooth trend to be extrapolated.

*A note on the Hopf check.* The first draft wrote the map as a quaternion
sandwich and paired the wrong components; the image missed the sphere by `0.96`.
It is written in complex coordinates now, where the fibre is visibly the common
phase.

## Scope

Everything here is about **measure, orientation and frames on round spheres**.

- **No field equation is solved, and nothing evolves.** No throat profile
  appears except through the `fⁿ` weighting §5 shares with PR #268.
- **§4 is an interpretation**, explicitly marked as one.
- **§6 is conditional.** It says what would have to be redone if the dimension
  moved; it does not argue that it should, and it does not touch any result
  currently derived at `d = 3`.
- **The Monte-Carlo sections are asserted as scalings**, never on a single
  sampled number — `std·√n → 1` is the claim, not any one row.
