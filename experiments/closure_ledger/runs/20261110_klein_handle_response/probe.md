# A localized, responding Klein handle

A lattice ring with one handle whose seam may reverse time and fibre; one waist coordinate per tube solved from the same Hamiltonian as the field; effective oscillator support; static prescribed potential. Classical throughout. No probability rule, spinor matching or statement about the cell BAM occupies.

Public freeze: `5d15963`. Seed `2026111002`. Cover sites: `32` and `52`; every required check must hold at both.

| verdict field | value |
|---|---|
| numerical_validity | **PASS** |
| support_dynamics | **EFFECTIVE_SUPPORT_ASSUMED** |
| electric_response | **EXTERNAL_POTENTIAL_PRESCRIBED** |
| localized_structure | **LAGRANGIAN_HALF_SURVIVES_A_LOCALIZED_RESPONDING_HANDLE** |
| local_bracket | **SHEET_DATA_FREE_AT_MIRROR_WRONSKIAN_NOT_CONSERVED** |
| neutral_sector_pairing | **NOT_SUPPLIED_BY_MECHANICAL_RESPONSE** |
| charged_sector_pairing | **CHARGE_PAIRS_THE_SECTORS_ONLY_WITH_BOTH_REVERSALS** |
| positive_response | **NOT_DERIVED** |
| spinor_matching | **NOT_DERIVED** |
| event_law | **NOT_DERIVED** |
| bam_cell | **NOT_ESTABLISHED** |

## Four seam types, charged field (primary size)

| s_t | s_f | fixed dim | Omega on fixed | charge on fixed | signature | Fix meet i Fix |
|---:|---:|---:|---|---|---|---:|
| +1 | +1 | 66 | rank 66 | nonzero | (32, 32) | 64 |
| +1 | -1 | 66 | rank 66 | zero | (0, 0) | 0 |
| -1 | +1 | 66 | zero | zero | (0, 0) | 64 |
| -1 | -1 | 66 | zero | nonzero | (32, 32) | 0 |

## Response (primary size)

| amplitude | seam defect | max mouth excursion | energy drift |
|---:|---:|---:|---:|
| 0.015625 | 2.5e-16 | 6.255e-06 | 2.6e-10 |
| 0.03125 | 3.3e-16 | 2.502e-05 | 1.0e-09 |
| 0.0625 | 3.1e-16 | 1.001e-04 | 4.1e-09 |
| 0.125 | 3.9e-16 | 4.002e-04 | 1.6e-08 |
| 0.5 | 2.5e-16 | 6.389e-03 | 2.6e-07 |

Supplementary, not frozen:

- amplitude `2.0`: seam defect `3.3e-16`, mouth excursion `0.115`, barrier factor `0.816`, energy drift `2.2e-06`
- amplitude `4.0`: seam defect `4.8e-16`, mouth excursion `0.446`, barrier factor `0.454`, energy drift `9.3e-06`

Order counting: mouth slopes `[2.0, 2.0, 1.9999]`, field-correction slopes `[3.0, 2.9999, 2.9998]`.

Sector coupling: imaginary part from mechanical response `0.0e+00`; from the potential, slopes `[1.0, 1.0]`. On the ordinary handle the potential gives slopes `[1.0, 1.0]` (supplementary).

| Required check (both sizes) | Pass |
|---|---|
| H1 deck map is an antisymplectic involution | True |
| H1 fixed space is a Lagrangian half | True |
| H1 frozen histories satisfy the seam relation | True |
| H1 even and odd spectra are the ordinary and twisted handle | True |
| H1 time-preserving seam keeps a symplectic subspace | True |
| H3 scheme is reversed by the deck map | True |
| H3 responding histories satisfy the seam relation at every amplitude | True |
| H3 tangent space is a Lagrangian half | True |
| H3 linearized tangent agrees with finite differences | True |
| H3 energy closes to 1e-5 for amplitudes up to 1/8 | True |
| H4 mouth response is second order in amplitude | True |
| H4 field correction is third order in amplitude | True |
| H4 the two mouths respond differently and mirror each other | True |
| H2 sheet data are free at the mirror time | True |
| H2 sheet Wronskian is symmetric about the mirror time | True |
| H2 sheet Wronskian is not conserved | True |
| H7 mechanical response generates no imaginary part | True |
| H7 the mouths were responding in those runs | True |
| H5 each deck map is an involution with the predicted symplectic sign | True |
| H5 charge density has parity s_t s_f | True |
| H5 fixed space is half, Lagrangian iff the seam reverses time | True |
| H5 charge vanishes identically in the two mixed cells | True |
| H5 twist is a quarter rotation in the antilinear cells only | True |
| H6 matching potential parity admits the seam | True |
| H6 mismatched potential parity breaks it | True |
| H7 charge is nondegenerate with balanced signature | True |
| H7 charge is the Wronskian pairing of the two sectors | True |
| H7 charge is conserved with potential and moving mouths | True |
| H7 sheet charge obeys Q1(t) + Q1(-t) = Q | True |
| H7 potential couples the sectors at first order | True |
| H8 charge is indefinite | True |
| H8 multiplication by i preserves no admissible field | True |

Passed 32/32 required checks.

An indefinite conserved pairing is not a probability, a norm or a spinor pairing. Which seam type BAM has is not established here.
