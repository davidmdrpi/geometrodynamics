# Round 12: odd multiplet support and preparation sensitivity

Prior candidate P1/P2 independently re-derived. P3 follows the published freeze amendment.

| Required check | Passed |
|---|---|
| harmonic_parity | True |
| addition_identities | True |
| full_stress_background | True |
| einstein_normalization | True |
| kinetic_matrix | True |
| component_bound | True |
| independent_field_controls | True |
| exact_preparation_kernel | True |
| normalized_sensitivity | True |
| gram_reconstruction | True |
| scope_and_order | True |
| failure_paths | True |

| Verdict | Result |
|---|---|
| background_existence | ODD_MULTIPLET_EXACT_ESU_SUPPORT |
| kinetic_regularity | POSITIVE_FULL_KINETIC_MATRIX |
| component_count_bound | FOUR_MINIMAL_IN_COMMON_PHASE_SINGLE_DEGREE_CLASS |
| diagonal_sensitivity | FINITE_FIXED_BACKGROUND_SUSCEPTIBILITY |
| full_preparation_kernel | {'1': 'ISOLATED_FIXED_TRACE_GRAM', '3': 'ISOLATED_FIXED_TRACE_GRAM', '5': 'EXACT_EQUAL_STRESS_FAMILY'} |
| coupled_dynamical_stability | NOT_ESTABLISHED |
| field_content_preparation_selection | NOT_DERIVED |
| TT_frequency_transfer | NOT_ESTABLISHED |
| coupled_support_response | NOT_DERIVED |
| Phi_selection | NOT_DERIVED |
| causality_gate | OPEN |
| failed_checks | [] |

| k | dim fixed trace | rank A | null A | rank L | null L |
|---|---|---|---|---|---|
| 1 | 9 | 9 | 0 | 9 | 0 |
| 3 | 135 | 135 | 0 | 135 | 0 |
| 5 | 665 | 581 | 84 | 581 | 84 |

No preparation selection, dynamical stability, tensor-frequency transfer, Phi selection or causality resolution is inferred.
