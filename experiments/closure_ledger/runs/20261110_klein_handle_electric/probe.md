# A responding electric field on the Klein handle

Lattice scalar electrodynamics on the cover network of a ring with one handle: charged field, link potentials, electric fields and two waist coordinates solved from one Hamiltonian. Effective oscillator support. Trapped flux enters as initial data. Classical throughout; no probability rule, norm, spinor matching or statement about the seam BAM has.

Public freeze: `4ab8bcd`. Seed `2026111003`. Phase-space dimensions `200` and `320`; every required check must hold at both sizes.

| verdict field | value |
|---|---|
| numerical_validity | **PASS** |
| support_dynamics | **EFFECTIVE_SUPPORT_ASSUMED** |
| electric_dynamics | **LATTICE_GAUGE_FIELD_SOLVED_FROM_THE_ACTION** |
| structure | **LAGRANGIAN_HALF_SURVIVES_ELECTRIC_RESPONSE** |
| flux_admissibility | **HANDLE_FLUX_ADMISSIBLE_IFF_CHARGE_UNTWISTED** |
| spontaneous_pairing | **NOT_GENERATED_BY_RESPONSE_WITHOUT_FLUX** |
| flux_pairing | **TRAPPED_FLUX_COUPLES_SECTORS_AND_SETS_ORIENTATION** |
| flux_dynamics | **FLUX_RESPONDS_AT_SECOND_ORDER** |
| pairing_carrier | **SHEET_PAIRING_IS_THE_SEAM_FLUX_BY_GAUSS** |
| seam_specificity | **FLUX_COUPLING_ALSO_ON_THE_ORDINARY_HANDLE** |
| flux_origin | **INITIAL_DATA_NOT_DERIVED** |
| positive_response | **NOT_DERIVED** |
| spinor_matching | **NOT_DERIVED** |
| event_law | **NOT_DERIVED** |
| bam_cell | **NOT_ESTABLISHED** |

## Admissible trapped flux (primary size)

| s_t | s_f | admissible flux dim | max on tube edges |
|---:|---:|---:|---:|
| +1 | +1 | 2 | 2.17e-01 |
| +1 | -1 | 1 | 4.46e-16 |
| -1 | +1 | 1 | 4.46e-16 |
| -1 | -1 | 2 | 2.17e-01 |

## Trapped flux, real field

| size | seam | Im slopes in flux | odd error | flux slopes in amplitude | max sheet charge |
|---|---|---|---:|---|---:|
| primary | klein | `[1.0, 0.9999]` | 0.0e+00 | `[2.0, 2.0001, 2.0005]` | 6.74e-07 |
| primary | ordinary | `[1.0, 0.9999]` | 0.0e+00 | `[2.0, 1.9998, 1.9994]` | 2.81e-19 |
| replicate | klein | `[1.0, 0.9999]` | 0.0e+00 | `[2.0, 1.9998, 1.9994]` | 5.20e-06 |
| replicate | ordinary | `[1.0, 0.9999]` | 0.0e+00 | `[2.0, 1.9999, 1.9998]` | 1.80e-19 |

## Supplementary, not frozen

| size | seam | screening at peak | max sheet charge | smallest per-field max sheet charge |
|---|---|---:|---:|---:|
| primary | klein | 4 of 12 | 1.50e-05 | 1.41e-06 |
| primary | ordinary | 3 of 12 | 4.63e-19 | 2.54e-20 |
| replicate | klein | 5 of 12 | 5.41e-06 | 5.19e-07 |
| replicate | ordinary | 7 of 12 | 3.15e-19 | 2.77e-20 |

| Required check (both sizes) | Pass |
|---|---|
| E1 each deck map is an involution with the predicted symplectic sign | True |
| E1 every deck map preserves H and maps the Gauss function covariantly | True |
| E1 fixed space is half, Lagrangian iff the seam reverses time | True |
| E1 scheme is reversed by the deck map | True |
| E1 mirror data are deck fixed and on the Gauss surface | True |
| E1 responding histories satisfy the seam relation | True |
| E1 time-preserving fixed set is invariant under the flow | True |
| E1 tangent space is a Lagrangian half | True |
| E1 Gauss law holds and total charge vanishes | True |
| E1 energy closes to 1e-5 | True |
| E1 the histories carried charge, field strength and mouth motion | True |
| E2 cover has three independent cycles | True |
| E2 handle loop is source free and swap even | True |
| E2 admissible flux dimension is 2 in matched cells, 1 in mixed | True |
| E2 flux through the tubes is admissible iff the cell is matched | True |
| E3 no imaginary part or gauge field from real data | True |
| E3 the mouths were responding in those runs | True |
| E4 imaginary part is first order in the trapped flux | True |
| E4 it is exactly odd in the flux, the real part even | True |
| E4 the flux histories are admissible | True |
| E5 flux response is second order in field amplitude | True |
| E5 flux response is first order in the flux | True |
| E6 sheet charge equals the flux through the seam | True |
| E6 total charge vanishes and the sheet charge is odd in time | True |
| E6 charge separates across the seam | True |
| E7 ordinary handle: flux couples the sectors at first order | True |
| E7 ordinary handle: exactly odd in the flux | True |

Passed 27/27 required checks.

Trapped flux is initial data here. An orientation is not a complex structure, and an indefinite pairing is not a norm.
