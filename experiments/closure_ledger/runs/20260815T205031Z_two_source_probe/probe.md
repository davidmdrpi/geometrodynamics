# probe — two_source

_2026-08-15T20:50:31.010310+00:00_

## T1_goal — PASS

- **question**: PR #253 ended rank counting by naming what it could not supply: a quantity that vanishes when a source is removed rather than becoming underdetermined. Build it as a field quantity, and ask whether it distinguishes a throat from two ordinary scatterers.
- **scope**: a linear scalar field on a fixed Einstein static universe. The throat is point-supported: no interior, no proper length, no delay. Everything is evaluated at a boundary matrix strictly inside PR #257's cone with its Loewner margin quoted, and the exact antipodal endpoint is tested separately. No backreaction, no stress tensor, no topology change, no rate.

## T2_it_vanishes_when_a_source_is_removed — PASS

- **n_pairs**: 40
- **boundary**:
    - 0.3
    - 0.35
    - 0.06
    - 0
- **loewner_margin**: 0.322775
- **inside_the_cone**: True
- **largest_value_with_a_source_removed**: 0
- **smallest_value_with_both_present**: 0.0366171
- **it_vanishes_exactly**: True
- **it_is_not_vacuous**: True
- **is_bilinear**: True
- **the_contrast**: a deleted equation costs a dimension whatever was deleted; a deleted source costs the whole cross term

## T3_the_throat_channel_is_rank_two — PASS

- **n_sources**: 12
- **loewner_margin**: 0.322775
- **chart_n_sources**: 12
- **chart_direct_rank**: 12
- **chart_throat_rank**: 2
- **chart_response_rank**: 2
- **chart_throat_norm**: 0.0386778
- **chart_direct_norm**: 0.263341
- **complex_strata**:
    - det_B=1.43974e-16, rank_B=1, rank_R=1, hermitian_defect=2.36248e-16, table_rank=2
    - det_B=3.37031e-17, rank_B=1, rank_R=1, hermitian_defect=1.24127e-16, table_rank=2
    - det_B=1.83212e-16, rank_B=1, rank_R=1, hermitian_defect=1.73404e-16, table_rank=2
    - det_B=2.49165e-16, rank_B=1, rank_R=1, hermitian_defect=2.73947e-16, table_rank=2
- **real_strata**:
    - det_B=3.8257e-17, rank_B=1, rank_R=1, hermitian_defect=1.8754e-16, table_rank=1
    - det_B=9.57335e-17, rank_B=1, rank_R=1, hermitian_defect=1.04256e-16, table_rank=1
    - det_B=7.28141e-16, rank_B=1, rank_R=1, hermitian_defect=8.75772e-15, table_rank=1
    - det_B=2.67746e-17, rank_B=1, rank_R=1, hermitian_defect=1.16759e-16, table_rank=1
- **the_throat_table_is_rank_two_in_the_chart**: True
- **the_direct_table_has_full_rank**: True
- **the_complex_response_has_the_rank_of_B**: True
- **a_complex_dirichlet_direction_still_fills_both_channels**: True
- **a_real_dirichlet_direction_drops_the_table_to_one**: True
- **the_response_is_hermitian_off_the_chart**: True
- **the_bound_is_two_at_every_source_count**: True

## T4_anisotropy_is_not_the_signature — PASS

- **loewner_margin**: 0.322775
- **chi**: 1
- **n_points**: 240
- **free_value**: 0.0644671
- **free_spread**: 8.32667e-17
- **throat_spread**: 0.0600309
- **throat_relative_spread**: 0.658594
- **disconnected_spread**: 0.0663908
- **disconnected_relative_spread**: 0.692051
- **the_free_field_is_isotropic**: True
- **both_break_it**: True
- **disconnected_over_connected**: 1.0508
- **anisotropy_does_not_discriminate**: True
- **what_it_does_show**: the interaction depends on more than the geodesic separation, which no free field on this background can do

## T5_the_disconnected_family_is_a_surface — PASS

- **n_draws**: 200
- **worst_defect_on_the_disconnected_family**: 1.38778e-16
- **smallest_defect_off_it**: 0.00500912
- **n_connected_draws_detected**: 199
- **the_disconnected_family_is_a_surface**: True
- **connected_throats_are_off_it**: True
- **the_equation**: S₁₂ = G₀ · det S,  two knobs and three observables

## T6_the_defect_is_the_coupling — PASS

- **n_draws**: 120
- **worst_error_in_W_plus_beta**: 4.996e-16
- **W_is_minus_beta**: True
- **margin_rows**:
    - loewner_margin=0.4, invariant=0.050965, defect=-0.06
    - loewner_margin=0.1, invariant=0.0791069, defect=-0.06
    - loewner_margin=0.02, invariant=0.150312, defect=-0.06
    - loewner_margin=0.004, invariant=0.196002, defect=-0.06
- **invariant_growth_toward_the_boundary**: 3.84582
- **worst_defect_drift_across_margins**: 2.08167e-17
- **every_row_is_inside_the_cone**: True
- **the_discriminator_does_not_see_the_resonance**: True
- **the_caution_it_answers**: PR #255: a test built from a resummed field can measure the pole instead of the source

## T7_it_is_a_protocol — PASS

- **n_observations**: 24
- **loewner_margin**: 0.322775
- **condition_number**: 9.85319
- **fit_residual**: 6.93889e-18
- **worst_entry_error**: 1.33227e-15
- **W_from_the_boundary_data**: -0.06
- **W_from_the_observations**: -0.06
- **W_error**: 1.11022e-16
- **minus_beta**: -0.06
- **the_protocol_closes**: True

## T8_the_blind_spot — PASS

- **separation**: 1.3
- **G_between_mouths**: 0.0484123
- **rows**:
    - alpha=[0.3, 0.35], re_beta=-0.05, branch=Re β < 0, im_beta=0.238993, abs_beta=0.244167, W=-2.08167e-17, loewner_margin=0.0906623, inside_the_cone=True, W_at_lambda_minus_one=0.0239973
    - alpha=[0.5, 0.4], re_beta=-0.02, branch=Re β < 0, im_beta=0.252889, abs_beta=0.253679, W=0, loewner_margin=0.208622, inside_the_cone=True, W_at_lambda_minus_one=0.0102898
    - alpha=[0.25, 0.25], re_beta=-0.1, branch=Re β < 0, im_beta=0.190361, abs_beta=0.215029, W=6.93889e-18, loewner_margin=0.033952, inside_the_cone=True, W_at_lambda_minus_one=0.0529442
    - alpha=[0.3, 0.35], re_beta=0.0684123, branch=Re β > G_d, im_beta=0.645221, abs_beta=0.648838, W=1.38778e-17, loewner_margin=-0.295685, inside_the_cone=False, W_at_lambda_minus_one=0.0515376
    - alpha=[0.5, 0.4], re_beta=0.0984123, branch=Re β > G_d, im_beta=0.659441, abs_beta=0.666744, W=6.93889e-18, loewner_margin=-0.187891, inside_the_cone=False, W_at_lambda_minus_one=0.0235888
- **the_blind_family_is_not_empty**: True
- **the_upper_branch_is_excluded_by_the_stability_gate**: True
- **the_lower_branch_survives_it**: True
- **they_are_invisible_at_lambda_zero**: True
- **they_are_visible_at_another_frequency**: True
- **largest_stable_invisible_coupling**: 0.253679
- **smallest_stable_invisible_margin**: 0.033952
- **the_cost**: a one-frequency two-source test cannot falsify a connected throat; the stability gate excludes only the Re β > G_d half of the blind family

## T9_two_frequencies_reconstruct_the_throat — PASS

- **n_draws**: 6
- **rows**:
    - true=[0.20077, 0.502979, 0.0446369, 0.125463], recovered=[0.20077, 0.502979, 0.0446369, 0.125463], max_parameter_error=1.66533e-16, residual=8.88178e-16, loewner_margin=0.180767
    - true=[0.400347, 0.423258, 0.194381, 0.148963], recovered=[0.400347, 0.423258, 0.194381, 0.148963], max_parameter_error=1.08247e-15, residual=4.44089e-16, loewner_margin=0.228259
    - true=[0.419509, 0.356846, 0.0438168, 0.168694], recovered=[0.419509, 0.356846, 0.0438168, 0.168694], max_parameter_error=8.32667e-17, residual=4.44089e-16, loewner_margin=0.241868
    - true=[0.293654, 0.465747, -0.149274, 0.126553], recovered=[0.293654, 0.465747, -0.149274, 0.126553], max_parameter_error=1.22125e-15, residual=1.55431e-15, loewner_margin=0.155032
    - true=[0.518626, 0.41012, -0.00415938, 0.284736], recovered=[0.518626, 0.41012, -0.00415938, 0.284736], max_parameter_error=6.10623e-16, residual=4.44089e-16, loewner_margin=0.195115
    - true=[0.358649, 0.398272, 0.0670723, 0.282549], recovered=[0.358649, 0.398272, 0.0670723, 0.282549], max_parameter_error=3.88578e-16, residual=1.77636e-15, loewner_margin=0.119935
- **worst_parameter_error**: 1.22125e-15
- **worst_residual**: 1.77636e-15
- **the_boundary_matrix_is_reconstructed**: True
- **blind_family_member**: lambdas=[0, -1], true=[0.3, 0.35, -0.05, 0.238993], recovered=[0.3, 0.35, -0.05, 0.238993], max_parameter_error=3.88578e-16, residual=8.88178e-16, sign_of_im_beta_is_not_observable=True
- **even_the_blind_family_is_reconstructed**: True
- **what_is_still_not_observable**: the sign of Im β — PR #256's time reversal, not a gap in the measurement

## T10_the_antipodal_endpoint — PASS

- **separation**: 3.14159
- **G_between_mouths_at_zero**: 0.0253303
- **minus_g0**: 0.0253303
- **the_antipodal_value_is_minus_g0**: True
- **rows**:
    - epsilon=0.01, loewner_margin=0.01, invariant=0.3913, W=3.46945e-18
    - epsilon=0.001, loewner_margin=0.001, invariant=3.42818, W=-3.46945e-17
    - epsilon=0.0001, loewner_margin=0.0001, invariant=33.7962, W=2.5327e-16
    - epsilon=1e-05, loewner_margin=1e-05, invariant=337.476, W=3.59088e-15
- **invariant_times_epsilon**:
    - 0.003913
    - 0.00342818
    - 0.00337962
    - 0.00337476
- **residue_of_the_divergence**: 0.00337476
- **the_invariant_diverges_like_one_over_epsilon**: True
- **the_defect_stays_zero**: True
- **W_of_a_connected_antipodal_throat**: -0.06
- **minus_beta**: -0.06
- **the_identity_survives_the_endpoint**: True
- **the_lesson**: at the marginal endpoint the signal is unbounded and the discriminator is exactly zero — size is not evidence

## T11_assessment — PASS

- **n_passed**: 10
- **n_total**: 10

## verdict — THE_TWO_SOURCE_INVARIANT_MEASURES_THE_MOUTH_MIXING

THE INVARIANT IS THE CROSS TERM OF A QUADRATIC FUNCTIONAL, AND ITS DISCONNECTION DEFECT IS MINUS THE MOUTH-MIXING AMPLITUDE. PR #253 ended rank counting by naming what it could not give: something that VANISHES when a source is removed rather than becoming underdetermined. Superposition makes every linear functional additive, so the object has to be quadratic, and its cross term is the throat's Green function between the two source points -- bilinear in the sources, exactly zero when either is switched off, and not vacuous (0.0366 at its smallest with both present). Written out it is PR #255's requested index, A MATRIX IN A PAIR OF BRANCHES: which mouth the field entered, which it left, plus the channel that used neither. THE THROAT CHANNEL IS RANK TWO AT ANY SOURCE COUNT -- the N x N table of throat-mediated cross terms is V^T S V with V of shape 2 x N, rank 2 against a direct table of rank 12 for 12 sources -- and off the chart rank R = rank B, though static sources see only Re R, whose rank is two even for a COMPLEX rank-one boundary condition and one for a real one. TWO THINGS THAT LOOK LIKE THE SIGNATURE AND ARE NOT. The cross term being nonzero is interference. And the interaction being ANISOTROPIC -- depending on more than the geodesic separation, which the free field on this background cannot do at all (8.3e-17 spread) -- is real, 0.66 of the mean, and two DISCONNECTED scatterers produce 0.69. It detects structure at the mouths, not a connection between them. WHAT DOES DISCRIMINATE IS A PARAMETER COUNT. The static invariant determines three numbers, the entries of S = Re R; two independent scatterers have two knobs, so their image is a SURFACE with the exact equation S12 = G0 det S, satisfied to 1.4e-16 on 200 draws. The defect W = S12/det S - G0 is therefore the discriminator -- and on real beta it is not merely nonzero but EQUAL TO THE COUPLING: W = -beta to 5.0e-16, independent of the self-energies, the separation, and the Loewner margin. That last independence answers PR #255's caution that a resummed field measures the pole instead of the source: driven from margin 0.4 to 0.004 the invariant grows 3.8x while W drifts 2.1e-17. AND IT IS A PROTOCOL: an observer who measures interaction energies and knows the background and the mouth positions, but is not told the boundary data, recovers S by least squares from source placements and gets W to 1.1e-16 at 24 observations, condition number 9.9. AGAINST THE ROUND: A ONE-FREQUENCY TEST IS BLIND ON A ONE-PARAMETER FAMILY. For complex beta, W = -Re beta - (G_d - Re beta)(Im beta)^2/P vanishes away from beta = 0 on two branches. PR #257's gate excludes one of them -- Re beta > G_d has det(A - Gamma) < 0 -- and leaves the other: connected throats with |beta| up to 0.254, strictly inside the cone at margin 0.034, that the static invariant cannot see at all. Not fine-tuned and not unstable. THE REPAIR IS MEASURED RATHER THAN ASSERTED: Gamma depends on lambda so the blind surface moves, two frequencies give six equations for four parameters, and the boundary matrix comes back to 1.2e-15 -- the blind family included. What is still not observable is the SIGN of Im beta, which PR #256 established is a time reversal. FINALLY, THE ANTIPODAL ENDPOINT ON ITS OWN, because PR #257 showed it is a different statement rather than a limit: at d = pi, Gamma(0) is negative semidefinite with a zero in the symmetric channel, so the static response is singular as A -> 0 and the invariant DIVERGES like 1/eps (residue 0.00337, flat across four decades). W stays exactly zero through all of it, and W = -beta still holds for a connected antipodal throat. The loudest available two-source signal carries no information about whether the mouths are connected -- size is not evidence. WHAT IS STILL PUT IN: the background, the mouth positions, and the boundary data itself, four real numbers chosen and not derived, with PR #249 still the thing that would fix them from matter. The throat is point-supported -- no interior, no proper length, no delay. No backreaction, no stress tensor, no topology change, no rate. What this round hands the next one is a quantity that is zero without a second source, equal to the non-local part of the operator when the throat is time-reversal invariant, and a stated blind spot with a stated repair.
