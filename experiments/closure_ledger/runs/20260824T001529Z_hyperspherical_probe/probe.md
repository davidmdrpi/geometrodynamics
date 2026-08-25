# Hyperspherical probe — what higher dimension does to the picture

**Question.** what does higher dimension actually do to the bulk and intersection picture -- as opposed to what 3-D intuition suggests it does?

**Answer.** a packet can keep exactly the same angular footprint while its physical overlap collapses as f0^n; almost every point is at the equator of any other, so the antipodal relation is vanishingly non-generic; and orientability of the antipodal quotient flips with dimension parity

**8/8 checks pass.**

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

The 2-D drawing understates the `S³` case by `1.0e+06`.

**So the picture is of** not two ribbons squeezing together, but a vast angular configuration space packed into an extremely small proper region — and whether two fronts meet was always an angular question, and stays one.

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

An entire s^n of directions collapses to one invariant point; the bearing that direction space at finite size, with measure f0^n |S^n|.

| `n` | `f₀` | `|S^n|` | bearing measure `f₀ⁿ|S^n|` |
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
