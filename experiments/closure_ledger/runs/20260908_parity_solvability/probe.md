# Does constraint solvability force antipodal parity?

Linearized Hamiltonian constraint at order phi^2 on the round S^3, zero supporting-matter perturbation, linear subspaces of pure-multiplet data. Not the nonlinear constraint, not evolution, not a variety-level characterization.

Public freeze: `495f1f185411a80cfd688b0b33291d1986f9702a`. Seed `2026090713`.

| verdict field | value |
|---|---|
| dipole_obstruction_structure | **BILINEAR_CROSS_PARITY_ADJACENT_DEGREE_ONLY** |
| antipodal_parity_status | **SUFFICIENT_NOT_NECESSARY** |
| f6_consequence | **CONSTRAINT_SOLVABILITY_DOES_NOT_DERIVE_THE_ANTIPODAL_CONDITION** |
| momentum_sector | **KILLING_DIAGONAL_IN_DEGREE_AND_GRADIENT_CKV_ADJACENT; NOT_PREDICTED_IN_ADVANCE** |
| triangle_map | **NOT_DERIVED** |
| readout | **NOT_DERIVED** |

## Selection rule

Minimum adjacent overlap norm `3.4641`; maximum non-adjacent `1.387e-14`.

| degrees | mixed parity | adjacent pair | max \|P^A\| |
|---|---|---|---:|
| `[1, 4]` | True | False | 5.618e-15 |
| `[2, 5]` | True | False | 1.523e-14 |
| `[3, 6]` | True | False | 2.862e-14 |
| `[1, 2]` | True | True | 5.076e+00 |
| `[2, 3]` | True | True | 6.933e+00 |
| `[3, 4]` | True | True | 9.032e+00 |
| `[1, 3]` | False | False | 0.000e+00 |
| `[2, 4]` | False | False | 0.000e+00 |
| `[1, 3, 5]` | False | False | 0.000e+00 |
| `[2, 4, 6]` | False | False | 0.000e+00 |

## Momentum sector (no prediction was frozen)

| degree | parity | Hamiltonian dipole | Killing charge |
|---:|---|---:|---:|
| 1 | odd | 0.000e+00 | 1.000000 |
| 2 | even | 0.000e+00 | 1.071142 |
| 3 | odd | 0.000e+00 | 3.745197 |
| 4 | even | 0.000e+00 | 6.992044 |

| Required check | Pass |
|---|---|
| selection rule holds for every pair up to degree 6 | True |
| predicted free mixed-parity sets have no obstruction | True |
| predicted obstructed sets do obstruct | True |
| parity-pure sets have no obstruction | True |
| overlap reduction matches the inherited improved stress | True |
| the obstruction is bilinear, measured independently | True |
| no nonzero subspace evades an adjacent partner | True |
| no nonzero graph subspace evades an adjacent pair | True |
| a small-projection mixed-parity subspace evades an adjacent pair | True |
| the complete six Killing and four gradient charges are audited | True |
| momentum charges are not a parity condition | True |

Passed 11/11 required checks.

### Implementation corrections

- C1: the overlap route summed ordered pairs (n,m) and (m,n), each with the full frozen P2 coefficient, double counting the total cross term. Found by the frozen improved-stress cross-check, which disagreed by exactly 2.000. One overall factor of 1/2 fixes it. No frozen prediction changed: a global factor cannot move a zero, so the selection rule and both counterexample classes are unaffected; only reported magnitudes halve.

Sufficiency of antipodal parity is not a derivation of it. No counting function, Born law, triangle map or readout follows.
