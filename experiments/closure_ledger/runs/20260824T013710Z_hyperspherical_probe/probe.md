# Hyperspherical probe — what higher dimension does to the picture

**Question.** what does higher dimension actually do to the bulk and intersection picture -- as opposed to what 3-D intuition suggests it does?

**Answer.** a packet can keep exactly the same angular footprint while its physical overlap collapses as f0^n; almost every point is at the equator of any other, so the antipodal relation is vanishingly non-generic; and orientability of the antipodal quotient flips with dimension parity -- with the exponent pinned to the object it belongs to, since n is the dimension of that object's own transverse sphere

**11/11 checks pass.**

| id | check | result |
|----|-------|--------|
| T0 | the peak is the ball's, and it needs a unit | PASS |
| T1 | *** almost all of a sphere is at the equator *** | PASS |
| T2 | *** the antipode is a vanishing-measure relation *** | PASS |
| T3 | point -> equatorial shell -> point | PASS |
| T4 | the bearing as the blown-up direction space | PASS |
| T5 | *** the collapse is f^n, not f *** | PASS |
| T6 | *** orientability flips with parity; two quotients *** | PASS |
| T7 | S^3 is not a generic sphere | PASS |
| T8 | *** which n is physical depends on which object *** | PASS |
| T9 | the finite centre is a routing manifold, not a hub | PASS |
| T10 | *** the limit separates incidence from measure *** | PASS |

## T0 — the peak is the ball's, and it needs a unit

The unit **ball** peaks at `d = 5`; the unit **sphere**'s surface peaks at `d = 7` — that is `S^6 in R^7`.

| `d` | `V_d` (ball) | `A_{d−1}` (sphere) | sphere |
|--|--|--|--|
| 1 | `2.00000` | `2.00000` | `S^0` |
| 2 | `3.14159` | `6.28319` | `S^1` |
| 3 | `4.18879` | `12.56637` | `S^2` |
| 4 | `4.93480` | `19.73921` | `S^3` |
| 5 | `5.26379` | `26.31895` | `S^4` |
| 6 | `5.16771` | `31.00628` | `S^5` |
| 7 | `4.72477` | `33.07336` | `S^6` |
| 8 | `4.05871` | `32.46970` | `S^7` |
| 9 | `3.29851` | `29.68658` | `S^8` |
| 10 | `2.55016` | `25.50164` | `S^9` |
| 11 | `1.88410` | `20.72514` | `S^10` |
| 12 | `1.33526` | `16.02315` | `S^11` |

And the peak is not a fact about dimension alone:

| `R` | ball peaks at | sphere peaks at |
|--|--|--|
| `0.5` | `d = 1` | `d = 2` |
| `1.0` | `d = 5` | `d = 7` |
| `2.0` | `d = 24` | `d = 26` |
| `4.0` | `d = 100` | `d = 102` |

**Why.** V_d and A_{d-1} carry different units at different d, so comparing them across d implicitly chooses a length scale — so the familiar claim is about the UNIT ball, and is not a dimension-invariant fact.

## T1 — the band narrows as `1/√n`

| `n` | `std(χ)` | `std·√n` | mass within `1/√n` of the equator |
|--|--|--|--|
| 2 | `0.683667` | `0.966852` | `0.6496` |
| 3 | `0.567862` | `0.983566` | `0.6587` |
| 4 | `0.495155` | `0.990311` | `0.6640` |
| 7 | `0.376711` | `0.996685` | `0.6716` |
| 15 | `0.258009` | `0.999264` | `0.6774` |
| 50 | `0.141412` | `0.999933` | `0.6811` |
| 200 | `0.070710` | `0.999996` | `0.6823` |
| 1000 | `0.031623` | `1.000000` | `0.6826` |

**Almost every point is nearly 90 degrees from any chosen point -- concentrated, not scattered.**

## T2 — so the antipode is non-generic

| ambient `n` | `std(x·y)·√n` | mean `α/π` | fraction with `α > 0.99π` |
|--|--|--|--|
| 3 | `0.998702` | `0.4996` | `3.20e-04` |
| 4 | `0.999082` | `0.4999` | `5.00e-06` |
| 10 | `1.003369` | `0.5004` | `0.00e+00` |
| 50 | `1.002857` | `0.5001` | `0.00e+00` |
| 500 | `1.001424` | `0.5000` | `0.00e+00` |

It selects a vanishing-measure relation out of an enormous set of nearly orthogonal alternatives, and gets more non-generic as dimension rises.

**What is not claimed:** that this makes the identification correct -- only that it is not the bland 'pair it with the far point' it can sound like.

## T5 — the collapse is `fⁿ`, not `f`

`transverse measure of an angular patch = f^n dOmega_n, so a squeeze from F to f0 costs (f0/F)^n`

| squeeze `f₀/F` | `S¹` | `S²` | `S³` | `S⁴` |
|--|--|--|--|--|
| `1e-01` | `1.0e-01` | `1.0e-02` | `1.0e-03` | `1.0e-04` |
| `1e-02` | `1.0e-02` | `1.0e-04` | `1.0e-06` | `1.0e-08` |
| `1e-03` | `1.0e-03` | `1.0e-06` | `1.0e-09` | `1.0e-12` |

**So the picture is of** not two ribbons squeezing together, but a vast angular configuration space packed into an extremely small proper region — and whether two fronts meet was always an angular question, and stays one.

## T8 — which `n` is physical depends on which object

N is the dimension of the angular cross-section of the object being represented -- not a single number shared by every throat quantity in the repository.

| object | cross-section | `n` | patch at `f₀/F = 1e-03` | vs the 2-D drawing |
|--|--|--|--|--|
| the drawn 2-D cross-section | `S^1` | 1 | `1.0e-03` | `1e+00` |
| PR #265 spatial throat (3-D slice) | `S^2` | 2 | `1.0e-06` | `1e+03` |
| bearing in a 4-spatial-D embedding | `S^3` | 3 | `1.0e-09` | `1e+06` |
| bearing in a 5-spatial-D embedding | `S^4` | 4 | `1.0e-12` | `1e+09` |

`physical_throat`'s neck area is `1.958592e-07`, which is `4π f₀²` to the last digit — an `f²` law, so the `#265` throat's own understatement against the drawing is `1e+03`, **not** the million. The millionfold figure belongs to *bearing in a 4-spatial-D embedding*.

**Why it matters.** Without this the same f^n law migrates between objects that do not share an n, and a figure derived for one gets quoted for another.

## T9 — the finite centre is a routing manifold, not a hub

| bearing `S^n` | capacity at `60°` | capacity at `20°` | proper measure at `f₀ = 1e-03` |
|--|--|--|--|
| `S^1` | `3.0` | `9` | `6.283e-03` |
| `S^2` | `4.0` | `33.16` | `1.257e-05` |
| `S^3` | `5.1` | `113.5` | `1.974e-08` |
| `S^4` | `6.4` | `374.1` | `2.632e-11` |
| `S^6` | `9.7` | `3818` | `3.307e-17` |
| `S^10` | `20.4` | `3.524e+05` | `2.073e-29` |
| `S^20` | `112.3` | `2.236e+10` | `2.929e-61` |

Holding `n = 3` and varying the neck makes the independence explicit — the capacity is the *same float*:

| `f₀` | proper measure | capacity at `20°` |
|--|--|--|
| `1e-01` | `1.974e-02` | `113.5295` |
| `1e-04` | `1.974e-11` | `113.5295` |
| `1e-07` | `1.974e-20` | `113.5295` |

**The reading:** a compressed routing manifold with its angular structure intact -- not a hub where things crowd together; the routing capacity is the same at any f0.

## T10 — what the `f₀ → 0` limit separates

| `f₀` | they meet? | overlap angle | proper overlap `n=2` | proper overlap `n=3` |
|--|--|--|--|--|
| `1e-01` | yes | `0.225` | `5.06e-04` | `1.14e-05` |
| `1e-02` | yes | `0.225` | `5.06e-06` | `1.14e-08` |
| `1e-03` | yes | `0.225` | `5.06e-08` | `1.14e-11` |
| `1e-06` | yes | `0.225` | `5.06e-14` | `1.14e-20` |
| `1e-09` | yes | `0.225` | `5.06e-20` | `1.14e-29` |

- angular incidence: survives
- overlap verdict: survives
- proper interaction measure: -> 0 as f0^n

**So the origin is** an entire finite direction space compressed to ZERO PROPER MEASURE with its angular structure intact -- not a place where every direction becomes equivalent.

*What the first reading got wrong:* treating the limit as merging the routes; it merges their sizes, not their labels.

## T6 — parity, and two quotients that are always opposite

| spatial `d` | spatial `ℝP^d` | exchange `ℝP^{d−1}` | `π₁`(exchange) | opposite? |
|--|--|--|--|--|
| 2 | `RP^2` **non-orientable** | `RP^1` orientable | `Z (braid)` | yes |
| 3 | `RP^3` orientable | `RP^2` **non-orientable** | `Z_2` | yes |
| 4 | `RP^4` **non-orientable** | `RP^3` orientable | `Z_2` | yes |
| 5 | `RP^5` orientable | `RP^4` **non-orientable** | `Z_2` | yes |
| 6 | `RP^6` **non-orientable** | `RP^5` orientable | `Z_2` | yes |

At this repository's dimension: spatial **RP^3 orientable**, exchange **RP^2 NON-orientable -- this is where the Pin- structure lives**.

**What is not claimed:** that the model should move dimension -- this says what would have to be re-derived if it did, nothing more.

## T4 — the bearing as the blown-up direction space

An entire s^n of directions collapses to one invariant point; the bearing restores that direction space at finite size, with measure f0^n |S^n|.

| `n` | `f₀` | `A_n` | bearing measure `f₀ⁿ·A_n` |
|--|--|--|--|
| 1 | `1e-02` | `6.2832` | `6.283e-02` |
| 1 | `1e-03` | `6.2832` | `6.283e-03` |
| 2 | `1e-02` | `12.5664` | `1.257e-03` |
| 2 | `1e-03` | `12.5664` | `1.257e-05` |
| 3 | `1e-02` | `19.7392` | `1.974e-05` |
| 3 | `1e-03` | `19.7392` | `1.974e-08` |
| 4 | `1e-02` | `26.3189` | `2.632e-07` |
| 4 | `1e-03` | `26.3189` | `2.632e-11` |

*Pr #268's bearing stands on its own geometry; this section only says what the construction can be read as doing, and nothing in regularized_center depends on it.*
