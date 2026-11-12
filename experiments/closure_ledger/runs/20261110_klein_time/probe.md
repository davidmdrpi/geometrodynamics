# Klein time attachment of a classical field

One linear zonal conformal scalar on the Einstein static universe, global quotients with caps placed by hand. Classical throughout; no hbar, probability rule or counting function is derived.

Public freeze: `bb6ff2c`. Lattice `N = 24`. Seed `2026111001`.

| verdict field | value |
|---|---|
| phase_space | **ADMISSIBLE_FIELDS_ARE_A_LAGRANGIAN_HALF** |
| complex_structure | **TIME_DECK_IS_A_CONJUGATION_OF_THE_MODE_SPACE** |
| phase_closure | **TWO_CAPS_FORCE_PHASE_0_OR_PI_FIXED_BY_TWIST** |
| field_freedom | **GENERIC_CLOSED_HISTORY_HAS_NO_FREE_FIELD** |
| hbar | **NOT_DERIVED** |
| born_rule | **NOT_DERIVED** |
| counting_function | **NOT_SELECTED** |
| lepton_index | **NOT_IDENTIFIED** |

## P1 Phase space

| cap | eta | dim + | dim - | max Omega on + | pairing min s.v. |
|---:|---:|---:|---:|---:|---:|
| 0 | +1 | 23 | 23 | 5.0e-16 | 0.7101 |
| 0 | -1 | 23 | 23 | 7.0e-16 | 0.7101 |
| 1 | +1 | 23 | 23 | 4.3e-17 | 1.0000 |
| 1 | -1 | 23 | 23 | 4.3e-17 | 1.0000 |
| 7 | +1 | 23 | 23 | 5.5e-15 | 0.1830 |
| 7 | -1 | 23 | 23 | 1.2e-15 | 0.1830 |
| 24 | +1 | 23 | 23 | 0.0e+00 | 0.7101 |
| 24 | -1 | 23 | 23 | 4.5e-15 | 0.7101 |

Orientable control: fixed dimension `24`, restricted Wronskian rank `24`.

## P3 Two caps

| m | tau/pi | untwisted modes | twisted modes |
|---:|---:|---|---|
| 48 | 1.0000 | `[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]` | `[]` |
| 24 | 0.5000 | `[2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22]` | `[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23]` |
| 16 | 0.3333 | `[3, 6, 9, 12, 15, 18, 21]` | `[]` |
| 12 | 0.2500 | `[4, 8, 12, 16, 20]` | `[2, 6, 10, 14, 18, 22]` |
| 7 | 0.1458 | `[]` | `[]` |

## P4 Closed history

Generic `m = 7`: smallest singular value `0.0171`, Green symmetry error `4.3e-15`, before/after ratio `1.000000` on the cap and `1.000000` off it.

Resonant `m = 24`: nullity `22`, residual `4.7e-13` for an orthogonal source and `4.1230` for a generic one (kernel component `4.1230`).

| Required check | Pass |
|---|---|
| lattice propagation is exactly periodic | True |
| deck map is an involution and reverses the Wronskian | True |
| each eigenspace has half the dimension | True |
| Wronskian vanishes on each eigenspace | True |
| the two eigenspaces pair nondegenerately | True |
| orientable control keeps a symplectic subspace | True |
| explicit marching confirms admissibility | True |
| j is a complex structure commuting with evolution | True |
| deck map anticommutes with j | True |
| j maps admissible onto twisted | True |
| metric is symmetric definite and h is real on admissible | True |
| frozen table equals the closure formula | True |
| untwisted two-cap modes are exactly the predicted set | True |
| twisted two-cap modes are exactly the predicted set | True |
| only the product of twists matters | True |
| periodic operator reproduces the two-cap count | True |
| generic closed history is uniquely solvable | True |
| its Green function is reciprocal and deck invariant | True |
| its response is not retarded | True |
| resonant history obeys the Fredholm alternative | True |

Passed 20/20 required checks.

These are classical kinematic structures. No quantum rule, lepton index or counting function follows from them.
