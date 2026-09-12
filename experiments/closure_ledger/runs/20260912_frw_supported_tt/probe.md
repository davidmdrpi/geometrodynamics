# Supported FRW tensor response

Freeze: `c15acd48e61d71378a1af2284767d52bfee94810`.

`(M beta')' + K beta = 0`, `M=A^2/kappa-q^2/6`, `K=8A^2/kappa+2q^2/3`.

Independent full-geometry cases: 55; transport backgrounds: 108.

| Gate | Pass |
|---|---|
| exact_background | True |
| action_derivation | True |
| full_stress_response | True |
| scalar_constraint_closure | True |
| independent_geometry | True |
| two_clocks | True |
| static_limit | True |
| canonical_transport | True |
| negative_controls | True |
| failure_paths | True |

```json
{
  "linear_tensor_sector": "CLOSED_FOR_ZERO_SCALAR_PERTURBATION_DATA",
  "supported_operator": "FULLY_SUPPORTED_FRW_EQUATION_VERIFIED",
  "clock_agreement": "CONFORMAL_AND_PROPER_AGREE",
  "finite_interval_transport": "CANONICAL_TRANSPORT_VERIFIED",
  "general_tensor_harmonics": "NOT_DERIVED",
  "nonlinear_stability": "NOT_ESTABLISHED",
  "tensor_quadratic_scalar_vector_sources": "NOT_DERIVED",
  "coupled_rotor": "NOT_DERIVED",
  "adiabatic_invariant": "NOT_TESTED",
  "Bogoliubov_bases": "NOT_SPECIFIED",
  "preparation_selection": "NOT_DERIVED",
  "Phi_selection": "NOT_DERIVED",
  "causality_gate": "OPEN",
  "failed_checks": []
}
```

These are finite-interval classical maps on a nonperiodic background, not Floquet or quantum verdicts.
