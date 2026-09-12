# Low scalar modes of the four-field ESU

Freeze: `8eb33d5348138cafd73c77f0a4ba66a285e99a73`.

| Quantity | Value |
|---|---|
| Homogeneous full-period trace | 85.031457203739 |
| Homogeneous multipliers | 85.019695223207, 0.011761980531 |
| Homogeneous determinant | 1.000000000000 |
| Constrained cover dipole: norm(M + I) | 4.269e-15 |
| Exact FRW: maximum Einstein residual | 2.721e-15 |

| Gate | Pass |
|---|---|
| field_reduction | True |
| homogeneous_constraints | True |
| exact_continuation | True |
| clock_and_gauge | True |
| dipole_reduction | True |
| dipole_parity | True |
| independent_geometry | True |
| period_maps | True |
| negative_controls | True |
| failure_paths | True |

```json
{
  "homogeneous_physical_block": "HYPERBOLIC_GROWING_MODE",
  "homogeneous_constraint_completion": "EXACT_FRW_CONTINUATION",
  "dipole_cover_block": "CONSTRAINED_LINEAR_NEUTRAL_SEMISIMPLE",
  "dipole_restricted_admissibility": "EXCLUDED_UNDER_STATED_ANTIPODAL_RESTRICTIONS",
  "vector_response": "NOT_DERIVED",
  "tensor_quadratic_sources": "NOT_DERIVED",
  "dipole_nonlinear_completion": "NOT_ESTABLISHED",
  "coupled_rotor": "NOT_DERIVED",
  "preparation_measure": "NOT_SPECIFIED",
  "Phi_selection": "NOT_DERIVED",
  "causality_gate": "OPEN",
  "failed_checks": []
}
```

The cover dipole is excluded only under the stated antipodal restrictions.
This linear instability does not specify a measure on viable preparations.
