# Klein time attachment of a classical field

One linear zonal conformal scalar on the Einstein static universe, global quotients with caps placed by hand. Classical throughout; no hbar, probability rule or counting function is derived. Corrected after review: the frozen 'no free field' reading holds only in the frozen truncation.

Public freeze: `bb6ff2c`. Corrected after review of `5b79975`. Lattice `N = 24`. Seed `2026111001`.

| verdict field | value |
|---|---|
| phase_space | **ADMISSIBLE_FIELDS_ARE_A_LAGRANGIAN_HALF** |
| complex_structure | **TIME_DECK_IS_A_CONJUGATION_OF_THE_COVER_MODE_SPACE** |
| phase_closure | **TWO_CAPS_FORCE_PHASE_0_OR_PI_FIXED_BY_TWIST** |
| field_freedom | **UNIQUE_SOURCED_SOLUTION_IN_THE_FROZEN_TRUNCATION** |
| continuum_free_field | **COUNTEREXAMPLE** |
| irrational_spacing | **OPEN_SMALL_DIVISORS** |
| attachment_contribution | **QUADRATURE_SELECTION_NOT_THE_TIME_SYMMETRIC_RESPONSE** |
| sector_pairing | **NOT_SUPPLIED** |
| hbar | **NOT_DERIVED** |
| born_rule | **NOT_DERIVED** |
| counting_function | **NOT_SELECTED** |
| lepton_index | **NOT_IDENTIFIED** |

**Counterexample (continuum_free_field).** u = sin(48 chi) sin(48 t) satisfies both untwisted caps at tau = 7 pi / 48; every rational tau / pi = p / q admits the untwisted mode k = q. This refutes the frozen label GENERIC_CLOSED_HISTORY_HAS_NO_FREE_FIELD as a continuum statement.

Frozen labels narrowed after review: `TIME_DECK_IS_A_CONJUGATION_OF_THE_MODE_SPACE` -> `TIME_DECK_IS_A_CONJUGATION_OF_THE_COVER_MODE_SPACE`; `GENERIC_CLOSED_HISTORY_HAS_NO_FREE_FIELD` -> `UNIQUE_SOURCED_SOLUTION_IN_THE_FROZEN_TRUNCATION`.

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

## P3 Two caps against orientable closed time

| m | tau/pi | harmonics (untwisted / twisted) | orientable dim | Klein dim |
|---:|---:|---|---|---|
| 48 | 1.0000 | `[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]` / `[]` | 46 / 0 | 23 / 0 |
| 24 | 0.5000 | `[2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22]` / `[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23]` | 22 / 24 | 11 / 12 |
| 16 | 0.3333 | `[3, 6, 9, 12, 15, 18, 21]` / `[]` | 14 / 0 | 7 / 0 |
| 12 | 0.2500 | `[4, 8, 12, 16, 20]` / `[2, 6, 10, 14, 18, 22]` | 10 / 12 | 5 / 6 |
| 7 | 0.1458 | `[]` / `[]` | 0 / 0 | 0 / 0 |

## Fixed physical cap separation `tau = 7 pi / 48`

| N | m | untwisted survivors | twisted survivors |
|---:|---:|---|---|
| 24 | 7 | `[]` | `[]` |
| 48 | 14 | `[]` | `[24]` |
| 72 | 21 | `[48]` | `[24]` |
| 96 | 28 | `[48]` | `[24, 72]` |

## P4 Closed history, frozen truncation only

`m = 7`, `N = 24`: smallest singular value `0.0171`, Green symmetry error `4.3e-15`.

Without the deck constraint, a source that is not deck invariant already has a unique response, time-symmetric to `1.8e-15`. The attachment admits `81` of `161` source dimensions and reduces the resonant solvability conditions from `22` to `11`.

Resonant `m = 24`: nullity `22`, residual `4.7e-13` for an orthogonal source and `4.1230` for a generic one.

## Continuum divisors (diagnostic, not a required check)

| tau/pi | min gap, k <= 1e5 | at k |
|---|---:|---:|
| 7/48 | 0.000000 | 48 |
| golden | 1.000000 | 1 |
| sqrt2 - 1 | 1.000000 | 1 |
| pi - 3 | 0.048114 | 113 |
| 0.110001 | 0.181816 | 100 |

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
| j preserves no admissible field of a fixed twist | True |
| frozen table equals the closure formula | True |
| untwisted two-cap modes are exactly the predicted set | True |
| twisted two-cap modes are exactly the predicted set | True |
| only the product of twists matters | True |
| periodic operator reproduces the two-cap count | True |
| frozen truncation at m = 7 is uniquely solvable | True |
| its Green function is reciprocal and deck invariant | True |
| its response is not retarded | True |
| resonant history obeys the Fredholm alternative | True |
| fixed-geometry refinement matches the continuum survivors | True |
| explicit continuum modes satisfy both caps at tau = 7 pi / 48 | True |
| every sampled rational spacing admits the mode k = q | True |
| orientable closed time selects the same harmonics, both quadratures | True |
| reciprocal time-symmetric response exists without the deck | True |
| attachment restricts sources to the deck-invariant half | True |
| attachment halves the solvability conditions at resonance | True |

Passed 28/28 required checks.

These are classical kinematic structures. No quantum rule, lepton index, counting function or sector pairing follows from them.
