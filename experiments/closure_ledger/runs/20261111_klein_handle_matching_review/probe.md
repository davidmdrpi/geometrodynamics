# The matching map

Part I: the repository's recorded maps on the bulk mouth S^3 acting on fibre-charged functions. Part II: a lattice multiplet on the handle cover with a sign around the doubled loop. Part III: the Pin sign against the trapped-flux sign in round 14's model. Classical throughout. Which field BAM's matter is, and P against PT, are not decided here.

Public freeze: `3d3ee38`. Seed `2026111004`. Lattice checks must hold at both sizes.

| verdict field | value |
|---|---|
| numerical_validity | **PASS** |
| antipodal_map | **FIBRE_TRANSLATION_FIXING_THE_HOPF_BASE** |
| two_spheres | **SPATIAL_BRANE_MOUTH_IS_NOT_THE_HOPF_BASE** |
| transport_on_fibre_modes | **ANTILINEAR_SQUARE_MINUS_ONE_ON_ODD_CHARGE_J_ON_THE_DOUBLET** |
| pr129_cell | **FIBRE_PRESERVING_TOTAL_CHARGE_ZERO_SHEET_EXCHANGE_EVEN** |
| spinor_seam_consistency | **TRANSPORT_SEAM_NEEDS_ANTIPERIODIC_DOUBLED_LOOP** |
| transport_seam_charges | **CHARGE_SURVIVES_SPIN_VANISHES_WITH_TIME_REVERSAL** |
| flux_sign_vs_pin_label | **INDEPENDENT_LABELS** |
| bam_cell | **DEPENDS_ON_FIELD_TYPE_AND_TIME_EXTENSION** |
| support_dynamics | **EFFECTIVE_SUPPORT_ASSUMED** |
| positive_response | **NOT_DERIVED** |
| event_law | **NOT_DERIVED** |

## Part I: the two recorded maps

| | antipode `-I_4` | transport `sigma` |
|---|---|---|
| on the fibre | translation by pi | reversal |
| on the Hopf base | identity | antipode |
| on the spatial brane mouth | antipode | leaves it (`0.998`) |
| on charge k | scalar `(-1)^k` | antilinear, square `(-1)^k` |

Doublet matrix of `A` in the basis `(z1, z2)`: `[[-0.0, -1.0], [1.0, 0.0]]`.

## Part II: when the transport seam admits a field (primary size)

| field | B | loop sign W | D squared | fixed dim |
|---|---|---:|---:|---:|
| singlet | conjugate | +1 | +1 | 64 |
| singlet | conjugate | -1 | -1 | 0 |
| doublet | J | +1 | -1 | 0 |
| doublet | J | -1 | +1 | 128 |

| seam | time | charge signature | spin signature |
|---|---|---|---|
| J, W = -1 | reversed | (64, 64) | (0, 0) |
| J, W = -1 | kept | (0, 0) | (64, 64) |
| identity, W = +1 | reversed | (0, 0) | (0, 0) |
| identity, W = +1 | kept | (64, 64) | (64, 64) |

## Part III: Pin sign against flux sign

Quarter rotation: flips the Pin sign, flux unchanged (`0.0e+00`). Conjugation: keeps the Pin sign, flux reversed (`0.0e+00`). A rotated history tested against the original Pin sign misses by `0.300`.

| Required check | Pass |
|---|---|
| M1 antipode is the fibre translation by pi | True |
| M1 it fixes the Hopf base | True |
| M1 it acts on charge k as the scalar (-1)^k | True |
| M2 antipode preserves the spatial brane mouth | True |
| M2 Hopf map identifies antipodes of the spatial brane mouth | True |
| review: it is two-to-one where p.i is nonzero | True |
| review: the whole circle p.i = 0 maps to -i | True |
| M2 transport does not preserve the spatial brane mouth | True |
| M2 transport covers the Hopf-base antipode and reverses the fibre | True |
| M2 transport is the recorded rotation with square -I | True |
| M3 pullback sends charge k to charge -k | True |
| M3 A squares to (-1)^k | True |
| M3 on the doublet A is J up to the Pin sign | True |
| review: neutral sheets still exchange charge, evenly in time, with zero total | True |
| M4 total charge is deck odd and vanishes for both twists | True |
| M5 D squares to B^2 W in all four cases | True |
| M5 H is invariant under D in all four cases | True |
| M5 singlet is admissible only on the periodic loop | True |
| M5 doublet with J is admissible only on the antiperiodic loop | True |
| M6 every seam is an involution with the predicted symplectic sign | True |
| M6 J seam, time reversed: charge survives, spin vanishes | True |
| M6 J seam, time kept: charge vanishes, spin survives | True |
| M6 linear seam, time reversed: both vanish | True |
| M6 linear seam, time kept: both survive | True |
| M6 responding histories satisfy the seam relation | True |
| M6 fixed space is a Lagrangian half | True |
| M6 charge and spin are conserved and total spin is zero | True |
| M6 sheet charge obeys Q1(t) + Q1(-t) = Q | True |
| M7 quarter rotation is a symmetry that leaves the flux alone | True |
| M7 quarter rotation flips the Pin sign | True |
| M7 conjugation is a symmetry that flips the flux and keeps the Pin sign | True |
| M7 both map admissible histories to admissible histories | True |
| M7 the rotated history is not admissible for the original Pin sign | True |

Passed 33/33 required checks.

Which cell BAM occupies depends on the field type and on the time extension; neither is decided here.
