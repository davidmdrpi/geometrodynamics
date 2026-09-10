# Can one real conformal scalar supply the ESU support?

The odd-sector result is a pointwise global obstruction in the stated class.
The positive homogeneous control is even and has a generically anisotropic response.

Public freeze: `de55f3f3175adafcf2cb760e0aef5767ca3e5016`.

| Verdict | Result |
|---|---|
| odd_sector_exact_ESU_support | EXCLUDED_IN_STATED_SINGLE_REAL_SCALAR_CLASS |
| homogeneous_even_control | EXACT_BUT_OUTSIDE_ODD_SECTOR |
| control_kinetic_regularity | REGULAR_POSITIVE_KINETIC_COEFFICIENTS |
| control_response | GENERIC_SCALAR_ANISOTROPIC_STRESS |
| control_constraints | PROPAGATED_FOR_TESTED_L_GE_2_RESPONSES |
| BAM_support_selection | NOT_DERIVED |
| triangle_map | NOT_DERIVED |
| Phi_selection | NOT_DERIVED |
| causality_gate | OPEN |

| Required check | Pass |
|---|---|
| stress_identity | True |
| global_obstruction | True |
| homogeneous_background | True |
| admissibility | True |
| response_variation | True |
| constraint_propagation | True |
| nonfluid_control | True |
| scope_and_order | True |
| fail_closed | True |

The numerical controls do not establish the global exclusion by search.
See docs/scalar_esu_support.md for the component-boundary proof.

| a | kappa | degree | Fine constraint max (absolute) | Refinement (absolute) | max abs(Pi) |
|---|---|---|---|---|---|
| 1.0 | 1.0 | 2 | 2.027e-11 | 1.065e-09 | 1.15667 |
| 1.0 | 1.0 | 3 | 2.467e-11 | 1.610e-09 | 1.15548 |
| 0.7 | 0.4 | 2 | 2.684e-11 | 1.467e-09 | 1.82714 |

The even control is not a BAM support selection. The triangle map, history
preparation/weight and physical readout remain unconstructed.
