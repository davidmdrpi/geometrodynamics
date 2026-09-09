# Field-derived driven normal stress

O(s^2) homogeneous TT projection with chosen constraint-compatible initial data; not a full Einstein-matter or triangle history

The nonzero tensor initial data are chosen; existence does not select this preparation.

Public freeze: `11e625fa877de2bb0669af06384347ede69d5777`

| Verdict | Result |
|---|---|
| normal_balance | LEADING_ORDER_FIELD_SUPPORTED_ROTOR |
| constraint_compatibility | CHOSEN_CONSTRAINT_COMPATIBLE_PREPARATION |
| homogeneous_tensor_evolution | VERIFIED_AT_ORDER_S2 |
| full_Einstein_matter_evolution | NOT_ESTABLISHED |
| preparation_selection | NOT_DERIVED |
| triangle_map | NOT_DERIVED |
| Phi_selection | NOT_DERIVED |

The frozen analytic candidate has S = 6 s^2 Q_z/(C a^2), A = -3 s^2/(2 C), Omega = sqrt(2)/a.
Its scalar and tensor frequencies have irrational ratio; no periodic joint history is inferred.

| Required gate | Pass |
|---|---|
| exact_polynomial_certificate | True |
| action_and_improved_stress | True |
| complete_constraint_compatibility | True |
| omitted_charge_controls | True |
| all_tensor_equations | True |
| unrestricted_tensor_evolution | True |
| necessity_controls | True |
| amplitude_radius_and_smallness | True |
| frequency_restriction | True |
| fail_closed_verdicts | True |

Passed 10/10 required gates.

| Normalized residual | Maximum |
|---|---:|
| constant_source balance | 1.9973e-15 |
| radial balance | 1.60177e-15 |
| angular balance | 5.04373e-16 |
| normal balance | 1.89855e-15 |
| full balance | 2.13248e-15 |
| trajectory_error | 4.31178e-12 |
| refinement_error | 4.37983e-10 |
| max_cone_distance | 3.82559e-12 |
| director_projector_error | 5.92226e-14 |
