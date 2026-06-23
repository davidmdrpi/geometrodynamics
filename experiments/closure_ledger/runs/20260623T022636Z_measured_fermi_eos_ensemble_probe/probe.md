# The measured Fermi equation of state: a many-throat ensemble (PR #172)

**Run:** 2026-06-23T02:26:36+00:00

The second closing option, on the same branch as #171 for comparison: SIMULATE a many-throat ensemble and MEASURE the equation of state, rather than assuming antisymmetry and reading off 5/3. *(QFT on the classical throat, not quantum gravity.)*

- **Ensemble**: free fermions in a box; −1 exchange sign = Pauli single-occupancy
- **Measured P/u**: 2/3 (NR), 1/3 (UR) from P = −dE/dV
- **Measured Γ**: 5/3 (NR), 4/3 (UR) from extrapolated level-filling
- **Control**: Bose Γ = 1, degeneracy pressure vanishes
- **Comparison**: matches #170 (assumed) and confirms #171 (topological sign)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | goal: measure the EoS from the ensemble (vs assume, #170) | **PASS** |
| T2 | `T2_many_throat_ensemble` | the ensemble: free fermions; −1 sign = Pauli filling | **PASS** |
| T3 | `T3_measure_pressure_over_u` | measure P/u from the volume derivative (virial) | **PASS** |
| T4 | `T4_measure_polytropic_index` | measure Γ → 5/3, 4/3 (finite-size-extrapolated) | **PASS** |
| T5 | `T5_bose_control` | Bose control: Γ=1, degeneracy pressure vanishes | **PASS** |
| T6 | `T6_three_route_comparison` | three routes agree (#170 assumed, #171 topological, measured) | **PASS** |
| T7 | `T7_honesty_and_scope` | honest scope (measured vs input vs idealizations) | **PASS** |
| T8 | `T8_assessment` | MEASURED_FERMI_EOS_GAMMA_5_3_AND_4_3 | **PASS** |

## The three routes to the index

| regime | assumed (#170) | topological (#171) | measured (#172) |
|---|---|---|---:|
| non-relativistic | 5/3 | −1 ⇒ Fermi ⇒ 5/3 | 1.6665 |
| ultra-relativistic | 4/3 | −1 ⇒ Fermi ⇒ 4/3 | 1.3332 |

(Bose control: Γ = 1, degeneracy pressure vanishes — the index tracks the exchange sign)

## Verdict

**MEASURED_FERMI_EOS_FROM_PAULI_FILLING_GAMMA_5_3_AND_4_3.** MEASURED, NOT ASSUMED. The Fermi equation of state is extracted from a simulated many-throat ensemble, and it matches the two earlier routes.

THE ENSEMBLE. N identical throats are free fermions in a cubic box; the −1 exchange sign (Pin⁻/geon, #170/#171) is realised as Pauli single-occupancy of the box modes, and the ground state is the filled Fermi sea built by level-filling — no occupation distribution assumed.

THE MEASUREMENTS. From the volume derivative P = −dE/dV: P/u = 0.6667 = 2/3 (non-relativistic) and 0.3333 = 1/3 (ultra-relativistic). From the finite-size-extrapolated level-filling: Γ = 1.6665 ≈ 5/3 (NR) and 1.3332 ≈ 4/3 (UR) — the indices are outputs of the simulation, not read off a formula.

THE CONTROL. With Bose statistics the index is Γ = 1.0000 = 1 and the T = 0 degeneracy pressure vanishes — so the 5/3 stiffening is a measured consequence of the −1 exchange sign, not a universal property of the box.

THE COMPARISON. The measured indices reproduce #170's assumed-analytic values and confirm the equation of state implied by #171's topological −1 sign. Assumed, topological, and measured — three routes, one answer. (Input not re-derived: the exchange sign itself; idealizations: free, T=0, box.)
