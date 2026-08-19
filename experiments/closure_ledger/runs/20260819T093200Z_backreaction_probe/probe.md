# probe — backreaction

_2026-08-19T09:32:00.280520+00:00_

## T1_goal — PASS

- **question**: PR #260's gate is closed by #261 and #262, so ask the roadmap's first GR question: does A + B produce a metric response that rescaling A or B alone cannot reproduce? Everything reachable by rescaling is the two-parameter family {c^2 beta_A + d^2 beta_B}, so the question is a projection residual.
- **scope**: a linear conformally coupled scalar on a FIXED ESU with PR #257's POINT throat; the metric response to FIRST ORDER in the lowest (n=3, homogeneous) transverse-traceless harmonic. Linear response, not a solved coupled system.

## T2_the_response_is_derived — PASS

- **omega_squared_times_a_squared**: 8
- **tensor_mode_frequency**: 2.82843
- **s3_ricci**: 2/a**2
- **s3_scalar**: 6/a**2
- **esu_einstein_diagonal**:
    - 3/a**2
    - -1/a**2
    - -1/a**2
    - -1/a**2
- **frw_g00**: 3*(Derivative(a(t), t)**2 + 1)/a(t)**2
- **delta_g_tt**: Derivative(b1(t), (t, 2)) + 8*b1(t)/a**2
- **the_validations_pass**: True
- **the_first_order_piece_is_traceless**: True
- **the_momentum_constraint_holds**: True
- **the_frequency_is_eight_over_a_squared**: True
- **the_sector_is_stable**: True
- **why_this_channel**: a perfect fluid has T_ab = diag(rho,p,p,p) in an orthonormal frame whatever the anisotropy, so it and Λ contribute nothing traceless — the only channel whose answer does not depend on the matter this arc never specified

## T3_the_split_is_bilinear — PASS

- **rows**:
    - scale=0.5, self_norm=0.0308866, cross_norm=0.069958, self_over_c_squared=1, cross_over_c=1
    - scale=1, self_norm=0.123546, cross_norm=0.139916, self_over_c_squared=1, cross_over_c=1
    - scale=2, self_norm=0.494185, cross_norm=0.279832, self_over_c_squared=1, cross_over_c=1
    - scale=3, self_norm=1.11192, cross_norm=0.419748, self_over_c_squared=1, cross_over_c=1
- **cross_with_one_source_off**: 0
- **the_self_term_is_quadratic**: True
- **the_cross_term_is_linear**: True
- **it_vanishes_with_one_source_off**: True
- **why_it_matters**: beta_A scales as c^2 and beta_cross as c, so the set reachable by rescaling is the two-parameter cone {c^2 beta_A + d^2 beta_B} and the question is whether beta_cross is in it

## T4_the_quadrature_control — PASS

- **points**:
    - 7802
    - 15990
- **volume_errors**:
    - 0.000878516
    - 0.000630816
- **rows**:
    - component=A, from=7802, to=15990, correlation=0.994222, magnitude_ratio=0.992684
    - component=B, from=7802, to=15990, correlation=0.970415, magnitude_ratio=1.00944
    - component=cross, from=7802, to=15990, correlation=0.999107, magnitude_ratio=0.993737
- **residual_by_level**:
    - 0.929072
    - 0.921509
- **worst_correlation**: 0.970415
- **worst_magnitude_drift**: 0.00944464
- **residual_drift**: 0.00756281
- **every_component_is_converged**: True
- **the_residual_is_stable**: True
- **what_it_caught**: a first draft's 0.982 was pure quadrature noise — independent rules for the same quantity correlated at -0.04

## T5_the_interference_is_unreachable — PASS

- **rows**:
    - window=[4, 12], unreachable_off_the_span=0.996405, cross_over_single=2.19943, norm_a=0.0343934, norm_b=0.0348924, norm_cross=0.075646, best_fit=[0.190296, -0.00999828], cos_between_the_singles=0.412413
    - window=[4, 20], unreachable_off_the_span=0.963523, cross_over_single=1.28094, norm_a=0.0775475, norm_b=0.0737621, norm_cross=0.0993338, best_fit=[0.149482, 0.279348], cos_between_the_singles=0.309324
    - window=[4, 30], unreachable_off_the_span=0.921509, cross_over_single=1.02147, norm_a=0.118052, norm_b=0.10307, norm_cross=0.120587, best_fit=[0.214025, 0.342599], cos_between_the_singles=0.172523
    - window=[8, 45], unreachable_off_the_span=0.87874, cross_over_single=0.963522, norm_a=0.151692, norm_b=0.144831, norm_cross=0.146159, best_fit=[0.295656, 0.328851], cos_between_the_singles=0.137374
- **points**: 15990
- **volume_error**: 0.000630816
- **unreachable_fraction**: 0.921509
- **cross_over_single**: 1.02147
- **cos_between_the_singles**: 0.172523
- **spread_over_windows**: 0.117664
- **most_of_it_is_unreachable**: True
- **the_two_singles_are_independent**: True
- **the_answer**: yes — the interference response is comparable in size to the single-wave ones and almost orthogonal to both, so no rescaling reproduces it

## T6_never_on_resonance — FAIL

- **rows**:
    - carrier=0.7, spectral_peak=5.96903, peak_offset_from_integer=0.030974, worst_peak_offset=0.141593, nearest_peak_to_the_mode=0.313166, power_at_the_mode=0.801585, unreachable=0.935588, cross_over_single=2.87854
    - carrier=1.41421, spectral_peak=5.96903, peak_offset_from_integer=0.030974, worst_peak_offset=0.141593, nearest_peak_to_the_mode=0.313166, power_at_the_mode=0.725345, unreachable=0.976452, cross_over_single=2.29169
    - carrier=2.82843, spectral_peak=5.96903, peak_offset_from_integer=0.030974, worst_peak_offset=0.261063, nearest_peak_to_the_mode=0.734032, power_at_the_mode=1.59106, unreachable=0.929072, cross_over_single=1.02039
    - carrier=4, spectral_peak=0.314159, peak_offset_from_integer=0.314159, worst_peak_offset=0.271387, nearest_peak_to_the_mode=3.1406, power_at_the_mode=2.87732, unreachable=0.564932, cross_over_single=1.59444
- **tensor_mode_frequency**: 2.82843
- **distance_to_the_nearest_integer**: 0.171573
- **grid_frequency_spacing**: 0.10472
- **spectral_peak_spread**: 5.65487
- **unreachable_range**:
    - 0.564932
    - 0.976452
- **closest_any_peak_gets_to_the_mode**: 0.313166
- **the_source_rings_on_integers**: True
- **no_peak_lands_on_the_tensor_mode**: False
- **the_tensor_mode_is_irrational**: True
- **it_is_unreachable_at_every_carrier**: True
- **the_fraction_is_not_a_universal_constant**: True
- **what_the_first_draft_got_wrong**: that T being quadratic puts the source's power at 2*w0, so the carrier should be chosen to match; the measured peak is the same for every carrier because the ringing is the ESU's own integer spectrum, not the pulse's

## T7_which_branches — PASS

- **with_throat**: unreachable=0.929072, cross_over_single=1.02039, norm_cross=0.121347, points=7802
- **no_throat**: unreachable=0.986019, cross_over_single=9.45856, norm_cross=0.102563, points=7802
- **unreachable_both_ways**: True
- **the_throat_changes_the_size**: 0.183141
- **what_it_scopes**: the conclusion survives switching the throat off, so it is a statement about two waves rather than about the throat; the throat changes the amplitude, and PR #257's rule that an integrated quantity must name its branches is why this control exists

## T8_assessment — FAIL

- **n_passed**: 6
- **n_total**: 7

## verdict — INCONCLUSIVE

INCONCLUSIVE. Failed checks: T6_never_on_resonance, T8_assessment
