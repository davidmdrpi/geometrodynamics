# Conditional scalar–TT constraint response

Each time is an independent linear CMC slice with zero support perturbations.
These are metric bounds; the complete scalar backreaction is not bounded here.

- conditional_constraints: **LINEAR_CMC_PARTICULAR_RESPONSE_BOUNDED_WITH_ZERO_SUPPORT_PERTURBATIONS**
- standing_wave_transverse_forcing: **ZERO_FOR_LEADING_STANDING_WAVE_ONLY**
- size_comparison: **CONTINUOUS_METRIC_BOUNDS_REPORTED**
- complete_scalar_backreaction_bound: **NOT_ESTABLISHED_WITHOUT_SUPPORT_AND_EVOLUTION_CLOSURE**
- triangle_map: **NOT_DERIVED**
- readout: **NOT_DERIVED**

| Continuous quantity | Value |
|---|---:|
| u minimum_rms | 0.00155814989228 |
| u maximum_rms | 0.00159731460862 |
| u_inhomogeneous minimum_rms | 0.000776381075802 |
| u_inhomogeneous maximum_rms | 0.000852261959117 |
| K_longitudinal_maximum_rms | 0.00147565629234 |
| j_maximum_rms | 0.00559346074107 |
| scalar_metric_inhomogeneous_all_time_lower | 0.00534975849632 |
| induced_TT_metric_all_time_upper | 0.00283640228638 |

| Time | u mean | u RMS without mean | K longitudinal RMS | scalar metric RMS without mean | induced TT metric norm |
|---:|---:|---:|---:|---:|---:|
| 0 | -0.001350949115 | 0.0007763810758 | 0 | 0.005378925877 | 0 |
| 1 | -0.001350949115 | 0.0008195486219 | 0.001459952722 | 0.00567799941 | 0.002278637747 |
| 2 | -0.001350949115 | 0.0008506316134 | 0.0004248463408 | 0.005893348692 | 0.0005489427374 |
| 3 | -0.001350949115 | 0.0007979756601 | 0.001336322408 | 0.005528537546 | 0.001793490477 |
| 4 | -0.001350949115 | 0.0007825753082 | 0.0008137162519 | 0.005421840778 | 0.0007603343271 |

| Required check | Pass |
|---|---|
| full_even_reconstruction | True |
| independent_improved_stress | True |
| Hamiltonian_and_momentum_residuals | True |
| dipole_and_CKV_compatibility | True |
| K_norm_identity | True |
| two_quadrature_rules | True |
| spectral_bounds_and_sharp_constants | True |
| incompatible_sources_rejected | True |
| quadratic_amplitude_scaling | True |
| independent_forced_TT_solution | True |
| continuous_envelopes | True |
| computed_homogeneous_TT_zeros | True |
| degree_one_anticommutator | True |
| exact_coherent_power_and_size_certificate | True |

Passed 14/14 checks.

A nonzero constraint response establishes neither an operational readout nor a completed Einstein history.
