# Regularized-center probe — the middle as a finite bearing

**Question.** can the centre keep its special role -- the clock-hand hinge whose cost does not care where the mouths are -- without being a point?

**Answer.** yes, and the hinge is cheaper than the proposal assumed: the geodesic turn cost is quadratic, T(alpha) = alpha^2/(2I), not the bearing's arc f0*alpha

| the bearing | |
|--|--|
| neck `f₀` | `1.000e-03` |
| outer scale / arm | `1` / `1.003647` |
| inner scale / arm | `0.35` / `0.353121` |
| arm ratio `L_o/L_i` | `2.8422` — **not one, and need not be** |
| resistance `I` | `3.996141e+03` |
| conductance `4π/I` | `3.144627e-03` |

**The law.** `T(alpha) = alpha^2 / (2I) = alpha^2 * conductance / (8 pi);  I -> 4/f0 gives f0 alpha^2/8`

At a half turn (`α = π`): geodesic hinge cost `1.1420e-03` against the linear guess `3.1416e-03` — and `8.42e-04` of the arms' own `1.356768`.

**12/12 checks pass.**

| id | check | result |
|----|-------|--------|
| T1 | the arms are the repo's geometry, symmetry dropped | PASS |
| T2 | *** the two arms are independent *** | PASS |
| T3 | scale transport is explicit | PASS |
| T4 | *** the turn cost is quadratic, not linear *** | PASS |
| T5 | the linear guess is an upper bound, and loose | PASS |
| T6 | *** the hinge is never the expensive part *** | PASS |
| T7 | the point model is the limit, not a rival | PASS |
| T8 | *** intersection becomes overlap on the bearing *** | PASS |
| T9 | the drawn circle is honest; the identification bites | PASS |
| T10 | the law is about necks; the large-angle deficit is not | PASS |
| T11 | *** I2 is the universal hinge; I4 first remembers the shape *** | PASS |
| T12 | *** monopole flux and infinitesimal rotation are ONE form *** | PASS |

## The turn cost, against the integrated geodesic

`T(alpha) = alpha^2 / (2I),  I = int ds/f^2` — and `T(alpha) = alpha^2 * conductance / (8 pi)`, so the geometric hinge and the monopole channel are one integral.  `I -> 4/f0, so T -> f0 alpha^2 / 8`

| `f₀` | `α` | `T(α)` | `α²/(2I)` | ratio | `f₀α` | `T/(f₀α)` |
|--|--|--|--|--|--|--|
| `1e-03` | `0.0200` | `5.00481e-08` | `5.00483e-08` | `0.999997` | `2.00e-05` | `0.00250` |
| `1e-03` | `0.0500` | `3.12795e-07` | `3.12802e-07` | `0.999979` | `5.00e-05` | `0.00626` |
| `1e-03` | `0.1000` | `1.25110e-06` | `1.25121e-06` | `0.999916` | `1.00e-04` | `0.01251` |
| `1e-03` | `0.3000` | `1.12524e-05` | `1.12609e-05` | `0.999248` | `3.00e-04` | `0.03751` |
| `1e-03` | `0.6000` | `4.49084e-05` | `4.50435e-05` | `0.997002` | `6.00e-04` | `0.07485` |
| `1e-03` | `1.0000` | `1.24085e-04` | `1.25121e-04` | `0.991722` | `1.00e-03` | `0.12408` |
| `1e-03` | `1.5708` | `3.02505e-04` | `3.08723e-04` | `0.979857` | `1.57e-03` | `0.19258` |
| `1e-03` | `3.1416` | `1.14201e-03` | `1.23489e-03` | `0.924785` | `3.14e-03` | `0.36351` |
| `1e-05` | `0.0200` | `5.00003e-10` | `5.00005e-10` | `0.999997` | `2.00e-07` | `0.00250` |
| `1e-05` | `0.0500` | `3.12497e-09` | `3.12503e-09` | `0.999979` | `5.00e-07` | `0.00625` |
| `1e-05` | `0.1000` | `1.24991e-08` | `1.25001e-08` | `0.999917` | `1.00e-06` | `0.01250` |
| `1e-05` | `0.3000` | `1.12417e-07` | `1.12501e-07` | `0.999251` | `3.00e-06` | `0.03747` |
| `1e-05` | `0.6000` | `4.48659e-07` | `4.50004e-07` | `0.997010` | `6.00e-06` | `0.07478` |
| `1e-05` | `1.0000` | `1.23969e-06` | `1.25001e-06` | `0.991745` | `1.00e-05` | `0.12397` |
| `1e-05` | `1.5708` | `3.02233e-06` | `3.08428e-06` | `0.979913` | `1.57e-05` | `0.19241` |
| `1e-05` | `3.1416` | `1.14114e-05` | `1.23371e-05` | `0.924968` | `3.14e-05` | `0.36324` |

**Proposed:** L_turn ~ f0 alpha, the bearing's arc length.

**Measured:** T(alpha) = alpha^2/(2I), quadratic in the angle and much smaller -- the geodesic cuts the corner rather than walking the arc.

The resistance is not a new quantity: at `a = 0.05` this module gives `3.1999986656e+04` and `physical_throat.resistance()` gives `3.1999986656e+04` — difference `0.0e+00`.

## The arms do not have to match

| `f_i` | `L_o` | `L_i` | `L_o/L_i` | `T(1 rad)` | `α²/(2I)` |
|--|--|--|--|--|--|
| `0.002` | `1.003647` | `0.002296` | `437.21` | `1.4461e-04` | `1.4649e-04` |
| `0.01` | `1.003647` | `0.011305` | `88.78` | `1.2718e-04` | `1.2832e-04` |
| `0.1` | `1.003647` | `0.102492` | `9.79` | `1.2430e-04` | `1.2535e-04` |
| `0.5` | `1.003647` | `0.503300` | `1.99` | `1.2406e-04` | `1.2509e-04` |
| `1` | `1.003647` | `1.003647` | `1.00` | `1.2403e-04` | `1.2506e-04` |

The old vacuole picture had to give the inner and outer boundaries one shared arbitrary radial gap; a finite bearing with two arms has no such constraint, and the r_inner/r_outer ratio stops carrying any physical significance.

## What intersection becomes

**Criterion.** they meet iff the angular extents overlap: |separation| < (w_a + w_b)/2 — f0 does not appear in the criterion.

| `f₀` | separation | `w_a` | `w_b` | meet? | overlap on the bearing | gap |
|--|--|--|--|--|--|--|
| `1e-02` | `0.05` | `0.20` | `0.20` | **yes** | `1.50e-03` | `0.00e+00` |
| `1e-02` | `0.30` | `0.20` | `0.20` | no | `0.00e+00` | `1.00e-03` |
| `1e-02` | `0.30` | `0.40` | `0.30` | **yes** | `5.00e-04` | `0.00e+00` |
| `1e-02` | `1.20` | `0.20` | `0.20` | no | `0.00e+00` | `1.00e-02` |
| `1e-03` | `0.05` | `0.20` | `0.20` | **yes** | `1.50e-04` | `0.00e+00` |
| `1e-03` | `0.30` | `0.20` | `0.20` | no | `0.00e+00` | `1.00e-04` |
| `1e-03` | `0.30` | `0.40` | `0.30` | **yes** | `5.00e-05` | `0.00e+00` |
| `1e-03` | `1.20` | `0.20` | `0.20` | no | `0.00e+00` | `1.00e-03` |

**The point limit, correctly stated.** as f0 -> 0 the ANGULAR INCIDENCE SURVIVES and the PHYSICAL INTERACTION REGION COLLAPSES: which directions the fronts come in at, and whether their extents overlap, are untouched by f0; the region in which they actually share space is f0 x (overlap angle) and goes to zero. So f0 -> 0 does NOT make every route meet -- it shrinks the overlap AND the gap together, and the distinction survives as a yes/no while disappearing as a length.
