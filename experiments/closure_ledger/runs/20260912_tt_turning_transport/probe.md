# TT turning and complete asymptotic transport

Freeze: `2b928d1a550e24019130544e9a4ff326969df965`.

The global turning bound is certified by exact quadratic-form identities and rational sign bounds.
Sampled R_turn range: 3.000007884569 to 3.017469157194.
The sampled extrema do not replace the frozen all-phase window [2.97,3.03].

| Gate | Pass |
|---|---|
| normalization | True |
| phase_reduction | True |
| turning_location | True |
| adiabatic_diagnostic | True |
| instantaneous_actions | True |
| frobenius_basis | True |
| past_basis | True |
| transport_convergence | True |
| symplectic_completion | True |
| negative_controls | True |
| scope | True |
| failure_paths | True |

```json
{
  "all_phase_turning_bound": "CERTIFIED_UNIQUE_SIMPLE_ROOT_IN_FROZEN_WINDOW",
  "named_instantaneous_actions": "NAMED_ACTION_PREDICTIONS_VERIFIED_NOT_A_GENERAL_NO_GO",
  "future_basis": "TWO_COLUMN_FROBENIUS_BASIS_VERIFIED",
  "complete_transport": "REAL_SYMPLECTIC_ASYMPTOTIC_MAP_VERIFIED",
  "all_classical_invariants_excluded": false,
  "quantization": "NOT_DERIVED",
  "preferred_complex_structure": "NOT_SPECIFIED",
  "coupled_rotor": "NOT_DERIVED",
  "general_support_independence": "NOT_ESTABLISHED",
  "Phi_selection": "NOT_DERIVED",
  "causality_gate": "OPEN",
  "failed_checks": {
    "B": [],
    "A": [],
    "F": [],
    "T": []
  }
}
```

The named instantaneous actions fail to stay invariant; chosen exact quadratic invariants still exist.
The complete map retains both frozen and decaying coefficients. No Phi or quantization is derived.
