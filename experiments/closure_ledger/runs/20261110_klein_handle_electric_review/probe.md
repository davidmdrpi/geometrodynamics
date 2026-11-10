# A responding electric field on the Klein handle

Lattice scalar electrodynamics on the cover network of a ring with one handle: charged field, link potentials, electric fields and two waist coordinates solved from one Hamiltonian. Effective oscillator support. Classical throughout; no probability rule, norm, spinor matching or statement about the seam BAM has. Corrected after review: trapped flux is sufficient but not necessary for sheet pairing, the frozen tangent test was ambient, and the imaginary part of the field is gauge dependent.

Public freeze: `4ab8bcd`. Seed `2026111003`. Phase-space dimensions `200` and `320`; every required check must hold at both sizes.

| verdict field | value |
|---|---|
| numerical_validity | **PASS** |
| support_dynamics | **EFFECTIVE_SUPPORT_ASSUMED** |
| electric_dynamics | **LATTICE_GAUGE_FIELD_SOLVED_FROM_THE_ACTION** |
| structure | **AMBIENT_LAGRANGIAN_HALF_SURVIVES_ELECTRIC_RESPONSE** |
| constrained_structure | **LAGRANGIAN_IN_THE_GAUSS_REDUCED_SPACE** |
| flux_admissibility | **HANDLE_FLUX_ADMISSIBLE_IFF_CHARGE_UNTWISTED** |
| spontaneous_pairing | **NOT_GENERATED_FROM_REAL_DATA_WITHOUT_FLUX** |
| matter_phase_route | **MATTER_CURRENT_GENERATES_SHEET_PAIRING_WITHOUT_FLUX** |
| flux_pairing | **TRAPPED_FLUX_DRIVES_GAUGE_INVARIANT_CURRENT_ODD_IN_FLUX** |
| gauge_dependence | **IMAGINARY_PART_IS_GAUGE_DEPENDENT** |
| flux_dynamics | **FLUX_RESPONDS_AT_SECOND_ORDER** |
| pairing_carrier | **SHEET_PAIRING_IS_THE_SEAM_FLUX_BY_GAUSS** |
| seam_specificity | **FLUX_COUPLING_ALSO_ON_THE_ORDINARY_HANDLE** |
| orienting_datum | **CONJUGATION_ODD_INITIAL_DATA_FLUX_OR_MATTER_CURRENT_NOT_DERIVED** |
| positive_response | **NOT_DERIVED** |
| spinor_matching | **NOT_DERIVED** |
| event_law | **NOT_DERIVED** |
| bam_cell | **NOT_ESTABLISHED** |

Corrected after review of `562af71`. Withdrawn:

- *response becomes dynamically active only through trapped handle flux.* Admissible Klein data with x swap even, y swap odd and Pi = a = E = 0 have no charge, no field and no flux, and generate sheet charge with the handle-loop flux zero throughout. Trapped flux is sufficient, not necessary.
- *the tangent test was done with Gauss's law imposed.* The frozen test perturbed along every deck-fixed direction; those violate linearized Gauss. It establishes an ambient property only.
- *the generated imaginary part is a physical orientation.* A gauge transformation preserving the deck condition changes it. The gauge-invariant statements are about link currents, sheet charge, loop holonomy and the electric field.

## Review: sheet pairing without trapped flux

| size | amplitude | max sheet charge | max handle flux | max field |
|---|---:|---:|---:|---:|
| primary | 0.015625 | 2.639e-06 | 2.1e-22 | 1.603e-05 |
| primary | 0.03125 | 1.056e-05 | 8.5e-22 | 6.411e-05 |
| primary | 0.0625 | 4.223e-05 | 2.8e-21 | 2.565e-04 |
| primary | 0.125 | 1.690e-04 | 1.1e-20 | 1.026e-03 |
| replicate | 0.015625 | 4.486e-06 | 6.4e-23 | 9.271e-06 |
| replicate | 0.03125 | 1.794e-05 | 3.4e-22 | 3.709e-05 |
| replicate | 0.0625 | 7.177e-05 | 1.5e-21 | 1.483e-04 |
| replicate | 0.125 | 2.871e-04 | 5.4e-21 | 5.934e-04 |

## Review: the Gauss surface

| size | ambient fixed | violation | constrained fixed | gauge in it | reduced dim | reduced fixed |
|---|---:|---:|---:|---:|---:|---:|
| primary | 100 | 0.7071 | 84 | 16 | 136 | 68 |
| replicate | 160 | 0.7071 | 134 | 26 | 216 | 108 |

## Review: gauge dependence (primary size)

Imaginary norm `0.0e+00` before a deck-preserving gauge transformation, `0.0326` after. Largest change in a gauge-invariant observable: `1.4e-15`. Gauge-invariant slopes in the flux: current `[0.9999, 0.9998]`, sheet charge `[1.0, 1.0]`, holonomy `[1.0, 1.0]`.

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
| R2 frozen tangent directions leave the Gauss surface | True |
| R2 the presymplectic kernel on the Gauss surface is the gauge directions | True |
| R2 constrained fixed space is isotropic | True |
| R2 it is half the reduced space after removing its gauge directions | True |
| R2 constrained tangents stay on the surface and isotropic along a responding history | True |
| R1 the data are admissible with no charge, field, momentum or link potential | True |
| R1 the gauge-invariant input is a nonzero link current | True |
| R1 sheet charge is generated at second order in amplitude | True |
| R1 handle-loop flux and holonomy stay zero throughout | True |
| R1 the histories are admissible and obey Gauss's law | True |
| R1 charge conjugation reverses the sheet charge | True |
| R1 the sheet charge is already there with the coupling off | True |
| R1 an independent integrator and a halved step agree | True |
| R3 a deck-preserving gauge transformation changes the imaginary part | True |
| R3 it leaves field, currents, charge, holonomy and energy unchanged | True |
| E2 cover has three independent cycles | True |
| E2 handle loop is source free and swap even | True |
| E2 admissible flux dimension is 2 in matched cells, 1 in mixed | True |
| E2 flux through the tubes is admissible iff the cell is matched | True |
| E3 no imaginary part or gauge field from real data | True |
| E3 the mouths were responding in those runs | True |
| R3 gauge-invariantly: no current, charge or field from real data | True |
| E4 imaginary part is first order in the trapped flux | True |
| E4 it is exactly odd in the flux, the real part even | True |
| E4 the flux histories are admissible | True |
| R3 link current, sheet charge and loop holonomy are first order in the flux | True |
| R3 they are exactly odd in the flux | True |
| E5 flux response is second order in field amplitude | True |
| E5 flux response is first order in the flux | True |
| E6 sheet charge equals the flux through the seam | True |
| E6 total charge vanishes and the sheet charge is odd in time | True |
| E6 charge separates across the seam | True |
| E7 ordinary handle: flux couples the sectors at first order | True |
| E7 ordinary handle: exactly odd in the flux | True |

Passed 45/45 required checks.

The orienting datum is conjugation-odd initial data, trapped flux or a matter current; neither is derived. An indefinite pairing is not a norm.
