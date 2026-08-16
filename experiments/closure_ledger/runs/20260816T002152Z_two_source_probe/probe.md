# probe — two_source

_2026-08-16T00:21:52.247904+00:00_

## T1_goal — PASS

- **question**: PR #253 ended rank counting by naming what it could not supply: a quantity that vanishes when a source is removed rather than becoming underdetermined. Build it as a field quantity, and ask whether it distinguishes a throat from two ordinary scatterers.
- **what_this_is_not**: NOT the roadmap's two-wave collision invariant. This is a static source-interaction kernel with no local null momenta, so it cannot separate collinear from counterpropagating waves. The index (i,j) labels mouth channels, not the geodesic/winding branches of PRs #253-#255.
- **scope**: a linear scalar field on a fixed Einstein static universe. The throat is point-supported: no interior, no proper length, no delay. Everything is evaluated at a boundary matrix strictly inside PR #257's cone with its Loewner margin quoted, and the exact antipodal endpoint is tested separately. No backreaction, no stress tensor, no topology change, no rate.

## T2_it_vanishes_when_a_source_is_removed — PASS

- **n_pairs**: 40
- **source_strengths**:
    - 1.7
    - -0.9
- **boundary**:
    - 0.3
    - 0.35
    - 0.06
    - 0
- **loewner_margin**: 0.322775
- **inside_the_cone**: True
- **a_self_energy_term_for_scale**: -0.0231977
- **worst_error_in_Q_minus_Q_minus_Q**: 2.77556e-17
- **largest_value_with_a_source_removed**: 0
- **smallest_value_with_both_present**: 0.0560242
- **the_cross_term_is_the_invariant**: True
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
- **what_is_actually_detected**: off-diagonal mouth-boundary mixing, relative to the diagonal two-scatterer null model

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
    - alpha=[0.3, 0.35], re_beta=-0.05, branch=Re β < 0, im_beta=0.238993, abs_beta=0.244167, smaller_than_the_self_energies=True, real_field_compatible=False, W=-2.08167e-17, loewner_margin=0.0906623, inside_the_cone=True, W_at_a_second_spectral_parameter=0.270761
    - alpha=[0.5, 0.4], re_beta=-0.02, branch=Re β < 0, im_beta=0.252889, abs_beta=0.253679, smaller_than_the_self_energies=True, real_field_compatible=False, W=0, loewner_margin=0.208622, inside_the_cone=True, W_at_a_second_spectral_parameter=1.00254
    - alpha=[0.25, 0.25], re_beta=-0.1, branch=Re β < 0, im_beta=0.190361, abs_beta=0.215029, smaller_than_the_self_energies=True, real_field_compatible=False, W=6.93889e-18, loewner_margin=0.033952, inside_the_cone=True, W_at_a_second_spectral_parameter=0.204258
    - alpha=[0.3, 0.35], re_beta=0.0684123, branch=Re β > G_d, im_beta=0.645221, abs_beta=0.648838, smaller_than_the_self_energies=False, real_field_compatible=False, W=1.38778e-17, loewner_margin=-0.295685, inside_the_cone=False, W_at_a_second_spectral_parameter=3.64956
    - alpha=[0.5, 0.4], re_beta=0.0984123, branch=Re β > G_d, im_beta=0.659441, abs_beta=0.666744, smaller_than_the_self_energies=False, real_field_compatible=False, W=6.93889e-18, loewner_margin=-0.187891, inside_the_cone=False, W_at_a_second_spectral_parameter=-2.14378
- **the_blind_family_is_not_empty**: True
- **the_upper_branch_is_excluded_by_the_stability_gate**: True
- **the_lower_branch_survives_the_stability_gate**: True
- **but_no_blind_point_is_real_field_compatible**: True
- **they_are_invisible_at_lambda_zero**: True
- **they_are_visible_at_a_second_spectral_parameter**: True
- **largest_stable_invisible_coupling**: 0.253679
- **smallest_stable_invisible_margin**: 0.033952
- **every_stable_coupling_is_smaller_than_its_self_energies**: True
- **the_scope**: the real-static-source protocol at a single λ is blind on this family — and the family exists only in a complex-scalar, time-reversal-breaking extension, not in PR #254's real field

## T9_a_real_field_forces_beta_real — PASS

- **separation**: 1.3
- **rows**:
    - beta=0.06+0j, real_field_compatible=True, max_imaginary_part_of_R=0, imaginary_part_of_the_field=0, has_an_invisible_partner=True
    - beta=0.06+0.2j, real_field_compatible=False, max_imaginary_part_of_R=2.43986, imaginary_part_of_the_field=0.00038531, has_an_invisible_partner=True
- **a_real_beta_gives_a_real_field**: True
- **a_complex_beta_does_not**: True
- **the_condition**: A = A* ⟺ the domain is conjugation-invariant
- **on_the_real_slice_W_is_minus_beta**: -0.06
- **the_blind_family_needs_a_complex_scalar**: True
- **so_for_PR254s_field_there_is_no_blind_family**: True

## T10_one_spectral_parameter_suffices_with_phase — PASS

- **separation**: 1.3
- **n_placements**: 8
- **boundary**:
    - 0.3
    - 0.35
    - -0.05
    - 0.24
- **the_quadratures_give_the_kernel**: True
- **in_phase**: 0.0569442
- **quadrature**: 0.00177262
- **condition_number**: 43.1711
- **worst_response_error**: 1.77721e-14
- **worst_boundary_error**: 3.88583e-15
- **one_spectral_parameter_suffices**: True
- **the_identity**: A = Γ(λ) + R(λ)⁻¹
- **what_this_corrects**: the two-λ requirement is a restriction of the real-static-source protocol, not of the operator

## T11_two_spectral_parameters_reconstruct_the_throat — PASS

- **n_draws**: 6
- **rows**:
    - true=[0.20077, 0.502979, 0.0446369, 0.125463], recovered=[0.20077, 0.502979, 0.0446369, 0.125463], max_parameter_error=2.77556e-17, residual=8.88178e-16, loewner_margin=0.180767
    - true=[0.400347, 0.423258, 0.194381, 0.148963], recovered=[0.400347, 0.423258, 0.194381, 0.148963], max_parameter_error=5.55112e-17, residual=3.55271e-15, loewner_margin=0.228259
    - true=[0.419509, 0.356846, 0.0438168, 0.168694], recovered=[0.419509, 0.356846, 0.0438168, 0.168694], max_parameter_error=5.55112e-17, residual=8.88178e-16, loewner_margin=0.241868
    - true=[0.293654, 0.465747, -0.149274, 0.126553], recovered=[0.293654, 0.465747, -0.149274, 0.126553], max_parameter_error=1.11022e-16, residual=8.88178e-16, loewner_margin=0.155032
    - true=[0.518626, 0.41012, -0.00415938, 0.284736], recovered=[0.518626, 0.41012, -0.00415938, 0.284736], max_parameter_error=1.11022e-16, residual=1.77636e-15, loewner_margin=0.195115
    - true=[0.358649, 0.398272, 0.0670723, 0.282549], recovered=[0.358649, 0.398272, 0.0670723, 0.282549], max_parameter_error=5.55112e-17, residual=8.88178e-16, loewner_margin=0.119935
- **worst_parameter_error**: 1.11022e-16
- **worst_residual**: 3.55271e-15
- **the_boundary_matrix_is_reconstructed**: True
- **blind_family_member**: lambdas=[0.3, 0.8], n_starts_tried=4, true=[0.3, 0.35, -0.05, 0.238993], recovered=[0.3, 0.35, -0.05, 0.238993], max_parameter_error=1.45717e-16, residual=7.10543e-15, sign_of_im_beta_is_not_observable=True
- **even_the_blind_family_is_reconstructed**: True
- **what_is_still_not_observable**: the sign of Im β — PR #256's time reversal, not a gap in the measurement

## T12_the_antipodal_endpoint — PASS

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

## T13_assessment — PASS

- **n_passed**: 12
- **n_total**: 12

## verdict — STATIC_THROAT_TOMOGRAPHY_MEASURES_THE_MOUTH_MIXING

THIS IS STATIC TWO-SOURCE THROAT TOMOGRAPHY, AND ITS DISCONNECTION DEFECT IS MINUS THE MOUTH-MIXING AMPLITUDE. FIRST, WHAT IT IS NOT: not the roadmap's two-wave collision invariant. The object is a STATIC source-interaction kernel at a fixed spectral parameter, it carries NO local null momenta, and it therefore cannot distinguish equal-energy collinear from counterpropagating waves -- the load-bearing control behind I_A I_B (k_A.k_B)^2. The index (i,j) labels MOUTH CHANNELS, not the geodesic/winding branches of PRs #253-#255. The dynamical object built from T_A^{mu nu} T^B_{mu nu} is still owed and is the next PR; the roadmap entry stays open. WHAT IS BUILT: PR #253 ended rank counting by naming what it could not give, something that VANISHES when a source is removed rather than becoming underdetermined. Superposition makes every linear functional additive, so the object has to be quadratic, and its cross term is the throat's Green function between the two source points. Computed the honest way -- the quadratic functional carries its own self-energy terms and the cross term is Q[a,b] - Q[a,0] - Q[0,b] from three separate evaluations, matching a*b*C to 2.8e-17; removing a source is evaluating the same functional at b = 0, not multiplying the answer by zero. THE THROAT CHANNEL IS RANK TWO AT ANY SOURCE COUNT: rank 2 against a direct table of rank 12 at 12 sources, because the table is V^T S V with V of shape 2 x N. Off the chart rank R = rank B, but static sources see only Re R, whose rank is two even for a COMPLEX rank-one boundary condition and one for a real one. TWO THINGS THAT LOOK LIKE THE SIGNATURE AND ARE NOT. The cross term being nonzero is interference. And ANISOTROPY -- the interaction depending on more than the geodesic separation, which the free field on this background cannot do at all (8.3e-17) -- is real at 0.66 of the mean, and two DISCONNECTED scatterers give 0.69. It detects structure at the mouths, not a connection between them, and the off-diagonal response block is the same trap one level down: Gamma fills it for diagonal boundary data too, so it is a CROSS-MOUTH channel and not 'through the throat'. WHAT DOES DISCRIMINATE IS A PARAMETER COUNT. The static invariant determines three numbers, the entries of S = Re R; two independent scatterers have two knobs, so their image is a SURFACE with the exact equation S12 = G0 det S, satisfied to 1.4e-16 on 200 draws. The defect W = S12/det S - G0 is its defining function, and on real beta it EQUALS THE COUPLING: W = -beta to 5.0e-16, independent of the self-energies, the separation and the Loewner margin -- driven from margin 0.4 to 0.004 the invariant grows 3.8x while W drifts 2.1e-17, which answers PR #255's caution that a resummed field measures the pole instead of the source. STATED EXACTLY, and this is the claim the round proves: W DETECTS OFF-DIAGONAL MOUTH-BOUNDARY MIXING RELATIVE TO THE DIAGONAL TWO-SCATTERER NULL MODEL, inside this point-interaction model -- not topology, not a traversable interior. AND IT IS A PROTOCOL: an observer who measures interaction energies and knows the background and the mouth positions, but not the boundary data, recovers S by least squares and gets W to 1.1e-16 at 24 placements. THE BLIND FAMILY, AND THE TWO THINGS THAT REMOVE IT. For complex beta, W = -Re beta - (G_d - Re beta)(Im beta)^2/P vanishes away from beta = 0 on two branches. PR #257's gate removes one, Re beta > G_d, where det(A - Gamma) < 0. REALITY OF THE FIELD REMOVES THE REST: PR #254 solves a REAL scalar, a real solution needs the self-adjoint domain conjugation-invariant, and that is exactly A = A*, hence beta real -- measured rather than argued, since with complex beta a real unit static source produces a field with Im = 3.85e-04. Every blind point needs Im beta != 0, so the family belongs to a deliberately time-reversal-breaking COMPLEX-scalar extension and not to the arc's field. Its stable couplings are COMPARABLE TO AND SMALLER THAN the self-energies (largest 0.254 at margin 0.034); an earlier draft said larger, which is false. AND EVEN INSIDE THAT EXTENSION THE LIMITATION IS THE PROTOCOL, NOT THE OPERATOR: real static sources see only Re R, three numbers for four parameters, while PHASE-SENSITIVE COMPLEX SOURCES see both quadratures and give the full complex R, whence A = Gamma + R^-1 at ONE spectral parameter, reconstructed to 3.9e-15. The real-static-source reconstruction needs two spectral parameters -- both POSITIVE and below the free ground state, since lambda = omega^2 makes a negative lambda an imaginary frequency rather than a second driving frequency -- and returns the boundary matrix to 1.1e-16, blind family included, using several starts because a single start does land in a local minimum and the reported residual is what catches it. FINALLY, THE ANTIPODAL ENDPOINT ON ITS OWN, because PR #257 showed it is a different statement rather than a limit: at d = pi the static response is singular as A -> 0 and the invariant DIVERGES like 1/eps (residue 0.00337, flat across four decades) while W stays exactly zero throughout. The loudest available two-source signal carries no information about whether the mouths are connected -- size is not evidence. WHAT IS STILL PUT IN: the background, the mouth positions, and the boundary data itself, four real numbers chosen and not derived, with PR #249 still the thing that would fix them from matter. Everything is STATIC-SOURCE, so this is an interaction-energy statement and not a scattering one. The throat is point-supported -- no interior, no proper length, no delay. No backreaction, no stress tensor, no topology change, no rate, AND NOT THE TWO-WAVE INVARIANT: what this round hands the next one is a measured non-locality, W = -beta, and a stated reason the dynamical object still has to be built.
