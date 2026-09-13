# Nonlinear supported tensor completion

Freeze: `a06778213efc0f92cca170d66e579a981617bbe3`.

| Gate | Pass |
|---|---|
| conventions_and_units | True |
| exact_field_closure | True |
| action_full_field_agreement | True |
| constraint_derivation | True |
| constraint_completion | True |
| constraint_propagation | True |
| linear_recovery | False |
| quadratic_response | False |
| finite_time_evolution | True |
| future_continuation_bound | True |
| negative_controls | True |
| scope | True |
| failure_paths | True |

```json
{
  "exact_reduction": "EXACT_HOMOGENEOUS_REDUCTION_VERIFIED",
  "constraint_completed_families": "LOCAL_CONSTRAINT_COMPLETED_FAMILIES_VERIFIED",
  "finite_amplitude_response": "UNRESOLVED",
  "future_persistence": "UNRESOLVED",
  "Phi_selection": "NOT_DERIVED",
  "quantization": "NOT_DERIVED",
  "causality_gate": "OPEN",
  "nonlinear_rotor": "NOT_DERIVED",
  "inhomogeneous_stability": "NOT_ESTABLISHED",
  "preparation_selection": "NOT_DERIVED",
  "tensor_only_symplectic_map": "NOT_ASSERTED",
  "explicit_initial_epsilon_radius": "NOT_COMPUTED",
  "failed_checks": {
    "D": [],
    "C": [],
    "N": [
      "linear_recovery",
      "quadratic_response",
      "evidence_linear_recovery",
      "evidence_quadratic_response"
    ],
    "F": [
      "linear_recovery",
      "evidence_linear_recovery"
    ]
  }
}
```
