# Reciprocal scalar–TT projection

Reciprocal variational ESU TT–scalar projection; full Einstein constraints, localized source, and two-boundary triangle map remain open.

Public freeze: `d8dc90d6d66e14824c337f96ee93a512dc9ed84f`. Seed: `2026090711`.

| Question | Verdict |
|---|---|
| reciprocal_dynamics | VARIATIONAL_TT_SCALAR_PROJECTION_VERIFIED |
| scalar_modal_closure | COMPLETE_MULTIPLETS_INVARIANT_IN_PROJECTED_SCALAR_EQUATION |
| einstein_constraints | OMITTED_METRIC_AND_SUPPORT_RESPONSE_REQUIRED |
| triangle_history_map | NOT_DERIVED |
| source_local_readout | NOT_DERIVED |
| probability_selection | NOT_DERIVED |

| History diagnostic | Value |
|---|---:|
| fine_relative_energy_drift | 6.73902256049e-12 |
| max_absolute_state_difference_between_tolerances | 3.70466990418e-11 |
| max_beta_frobenius | 0.00818775623292 |
| max_scalar_difference_from_one_way | 0.00108059045007 |
| max_tensor_difference_from_one_way | 1.24281305775e-05 |
| one_way_relative_defect_in_reciprocal_H | 0.00222565537358 |

| Time | Mean energy density | Inhomogeneous energy RMS | Momentum RMS |
|---:|---:|---:|---:|
| 0 | 0.0162113893828 | 0.0185345379662 | 0 |
| 1 | 0.0161950883862 | 0.0221626288451 | 0.0055306342879 |
| 4 | 0.016189952316 | 0.018990014622 | 0.00306290525525 |

| Required check | Pass |
|---|---|
| complete harmonic multiplets and invariant derivative algebra | True |
| n1 null coupling and n3 field quadrupole | True |
| action source matches inherited improved stress on both grids | True |
| pointwise scalar equation matches the modal equation | True |
| static matter-action remainder is quadratic | True |
| Hamilton equations from independent energy differences | True |
| reciprocity holds and fails in the one-way control | True |
| ODE refinement and Hamiltonian conservation | True |
| tensor remains within the frozen small-field range | True |
| free controls agree with independent harmonic solutions | True |
| action and source are covariant under common SO3 rotations | True |
| omitted constraint sources converge and match an analytic certificate | True |
| Q_m has the uniaxial initial identity | True |

Passed 13/13 numerical checks.

The CLI failure/overwrite path is verified separately by an end-to-end test.

No counting function, Born law, physical source-local readout, or retrocausal channel is inferred.
