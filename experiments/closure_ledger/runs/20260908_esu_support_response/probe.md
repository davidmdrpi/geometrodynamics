# ESU support response and first cubic reaction

The support class is assumed; its sound speed is a parameter. The proper-time coefficient is preparation-scoped.
The continuous metric bounds apply to t in [0,2], with a=kappa=1 and s=0.02.

- support_class: **PERFECT_FLUID_ADIABATIC_RESPONSE_ASSUMED_WITH_PARAMETER_CS2**
- scalar_constraints: **PROPAGATED_AT_ORDER_S2_IN_SCALAR_SECTOR**
- initial_proper_coefficient: **MINUS_7976_OVER_875_TIMES_KAPPA_S3_OVER_V_A2**
- cancellation: **IDENTICAL_CANCELLATION_EXCLUDED_FOR_FROZEN_PREPARATION**
- complete_evolution: **NOT_ESTABLISHED**
- BAM_support_selection: **NOT_DERIVED**
- readout: **NOT_DERIVED**

| Initial projection, units kappa s^3/(V a^2) | Exact coefficient |
|---|---:|
| Fluid proper-time scalar force | -7976/875 |
| Newtonian coordinate-time scalar force | 55096/875 |
| Induced homogeneous TT force | 0 |

The initial TT zero follows from the frozen preparation; the independent TT check compares the induced response over the interval.

Both clocks are also evaluated in the same exponential test metric, using its connection for the initial proper derivative.
Maximum two-clock Richardson field/coefficient error: 3.1e-08.

| c_s^2 | Unused Hamiltonian residual | Unused momentum residual | All-space/all-time psi upper bound |
|---:|---:|---:|---:|
| 0.333333 | 5.05e-13 | 5.08e-15 | 0.00100504 |
| 0 | 6.22e-13 | 6.25e-15 | 0.00105686 |
| 0.2 | 5.3e-13 | 5.32e-15 | 0.000854372 |
| 1 | 2.71e-13 | 6.59e-15 | 0.0056069 |

| Required check | Pass |
|---|---|
| metric_Einstein_derivation | True |
| improved_stress_and_scalar_anisotropy | True |
| two_spatial_rules | True |
| independent_fluid_evolution | True |
| propagated_constraints | True |
| ODE_refinement | True |
| trace_curvature | True |
| fluid_only_frequency | True |
| nonlinear_metric_variation | True |
| initial_proper_and_coordinate_coefficients | True |
| clock_conversion | True |
| amplitude_and_radius_scaling | True |
| induced_TT_force | True |
| small_metric_regime | True |
| exact_preparation_obstruction | True |

Passed 15/15 required checks.

No full corrected scalar evolution, microscopic BAM support, or source-local readout is derived.
