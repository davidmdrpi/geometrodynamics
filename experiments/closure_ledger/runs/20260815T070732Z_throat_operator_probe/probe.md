# probe — throat_operator

_2026-08-15T07:07:32.384108+00:00_

## T1_goal — PASS

- **question**: what does making the throat a genuine self-adjoint extension buy -- and what does it not? in particular, does flux conservation imply stability?
- **scope**: a linear scalar field on a fixed Einstein static universe. The throat is POINT-SUPPORTED: no interior, no proper length, and no delay -- Delta is not a parameter of a point extension. The boundary data is a choice of four real parameters, not a derivation; PR #249 is what would fix it. No backreaction, no stress tensor, no topology change, no rate, no two-source invariant.

## T2_the_green_function_has_a_finite_part — PASS

- **rows**:
    - omega=0.37, chi=0.4, closed_form=0.189076, branch_series_limit=0.189076, abs_error=1.33643e-13
    - omega=0.37, chi=1.3, closed_form=0.0566811, branch_series_limit=0.0566811, abs_error=1.13382e-13
    - omega=0.37, chi=2.6, closed_form=0.0334809, branch_series_limit=0.0334809, abs_error=1.00815e-13
    - omega=1.63, chi=0.4, closed_form=0.216089, branch_series_limit=0.216089, abs_error=8.539e-13
    - omega=1.63, chi=1.3, closed_form=-0.0125391, branch_series_limit=-0.0125391, abs_error=1.59303e-13
    - omega=1.63, chi=2.6, closed_form=-0.12994, branch_series_limit=-0.12994, abs_error=9.40942e-13
    - omega=2.4, chi=0.4, closed_form=0.0628065, branch_series_limit=0.0628065, abs_error=4.35624e-13
    - omega=2.4, chi=1.3, closed_form=-0.0831472, branch_series_limit=-0.0831472, abs_error=3.0885e-13
    - omega=2.4, chi=2.6, closed_form=0.156391, branch_series_limit=0.156391, abs_error=8.87734e-13
    - omega=5.21, chi=0.4, closed_form=-0.329837, branch_series_limit=-0.329837, abs_error=6.33893e-12
    - omega=5.21, chi=1.3, closed_form=0.0227861, branch_series_limit=0.0227861, abs_error=5.4292e-13
    - omega=5.21, chi=2.6, closed_form=-0.0792028, branch_series_limit=-0.0792028, abs_error=2.20288e-12
- **worst_abs_error**: 6.33893e-12
- **the_closed_form_is_the_branch_series**: True
- **finite_part**:
    - omega=0.37, chi=0.01, finite_part=-0.0126634, g=-0.0127414, error=7.79758e-05
    - omega=0.37, chi=0.001, finite_part=-0.0127336, g=-0.0127414, error=7.814e-06
    - omega=0.37, chi=0.0001, finite_part=-0.0127406, g=-0.0127414, error=7.81565e-07
    - omega=1.63, chi=0.01, finite_part=0.055205, g=0.0561311, error=0.000926061
    - omega=1.63, chi=0.001, finite_part=0.0560386, g=0.0561311, error=9.24673e-05
    - omega=1.63, chi=0.0001, finite_part=0.0561218, g=0.0561311, error=9.24533e-06
    - omega=2.4, chi=0.01, finite_part=-0.0642093, g=-0.0620551, error=0.00215421
    - omega=2.4, chi=0.001, finite_part=-0.062271, g=-0.0620551, error=0.000215871
    - omega=2.4, chi=0.0001, finite_part=-0.0620767, g=-0.0620551, error=2.15915e-05
- **convergence_ratios_per_decade**:
    - 9.97898
    - 9.99789
    - 10.015
    - 10.0015
    - 9.97914
    - 9.99794
- **the_remainder_is_first_order_in_chi**: True
- **antipode**:
    - omega=0.7, limit=0.0688542, at_pi_minus_1e-5=0.0688542, relative_error=3.74647e-12
    - omega=2.3, limit=0.226235, at_pi_minus_1e-5=0.226235, relative_error=8.37464e-11
    - omega=4.4, limit=0.36816, at_pi_minus_1e-5=0.36816, relative_error=3.18247e-10
- **the_antipodal_focus_is_finite**: True
- **what_this_shows**: the coincidence divergence is the universal 1/(4πχ), so g(ω) exists and a point interaction is definable

## T3_hermiticity_is_flux_conservation — PASS

- **n_draws**: 200
- **worst_relative_net_flux**: 1.83305e-16
- **flux_is_conserved_identically**: True
- **worst_pairwise_imbalance**: 1.72379e-16
- **what_one_mouth_absorbs_the_other_emits**: True
- **control_median_net_flux**: 0.538982
- **the_control_does_not_conserve**: True
- **cayley**:
    - scale=0.05, diagonal=0.955331, off_diagonal=0.295538, unitarity_defect=2.22617e-16, column_norm=1
    - scale=0.1, diagonal=0.854761, off_diagonal=0.519021, unitarity_defect=2.22187e-16, column_norm=1
    - scale=0.2, diagonal=0.713089, off_diagonal=0.701074, unitarity_defect=2.22072e-16, column_norm=1
- **the_cayley_transform_is_unitary**: True
- **cayley_diagonal_spread_over_the_reference_scale**: 0.242242
- **the_cayley_entries_are_not_physical_amplitudes**: True
- **the_correction**: an earlier version of this module read the Cayley entries as reflection and transmission; their magnitudes depend on the arbitrary scale c, so they are boundary-mixing coefficients and the physical conservation result is the flux identity
- **the_identity**: Im(q† A q) = 0 for all q  ⟺  A = A†

## T4_self_adjointness_makes_lambda_real_not_omega — PASS

- **rows**:
    - alpha1=0.05, alpha2=0.05, beta=0.03, hermitian=True, worst_relative_imaginary_part_of_det=0, n_growing_modes=0, lambda_of_the_growing_modes=[], growth_rates=[], stable=True
    - alpha1=0.2, alpha2=-0.13, beta=0.15+0.07j, hermitian=True, worst_relative_imaginary_part_of_det=2.13718e-15, n_growing_modes=1, lambda_of_the_growing_modes=[-6.10353], growth_rates=[2.47053], stable=False
    - alpha1=-0.4, alpha2=0.07, beta=-0.09+0.31j, hermitian=True, worst_relative_imaginary_part_of_det=2.88658e-15, n_growing_modes=1, lambda_of_the_growing_modes=[-50.282], growth_rates=[7.09098], stable=False
- **worst_relative_imaginary_part**: 2.88658e-15
- **the_secular_function_is_real_in_lambda**: True
- **hermiticity_gives_real_lambda**: True
- **n_unstable_examples**: 2
- **hermiticity_does_not_give_positivity**: True
- **the_withdrawn_claim**: 'the spectrum is real for every coupling, so a conserving throat cannot ring up' — false: λ is real, ω need not be
- **why_it_was_missed**: the earlier complex-root search seeded only Re ω in [1.1, 6.9] and discarded roots leaving that window, so a root on the imaginary axis was outside its reach by construction
- **what_replaces_it**: a positivity audit: solve in λ, and map the region of boundary data with λ_min ≥ 0

## T5_the_stability_region_in_the_boundary_family — PASS

- **thresholds**: g_at_zero=-0.0253303, G_d_at_zero=0.0484123, symmetric_threshold=0.023082, antisymmetric_threshold=-0.0737426, the_rule=a negative-λ mode exists iff α+β < symmetric threshold or α−β < antisymmetric threshold
- **probes**:
    - alpha=0.05, beta=0.03, stable=True, closed_form_stable=True, agrees=True, n_growing_modes=0, worst_growth_rate=0
    - alpha=0, beta=0, stable=False, closed_form_stable=False, agrees=True, n_growing_modes=1, worst_growth_rate=0.442021
    - alpha=-0.05, beta=0, stable=False, closed_form_stable=False, agrees=True, n_growing_modes=1, worst_growth_rate=0.925107
    - alpha=0.05, beta=0.3, stable=False, closed_form_stable=False, agrees=True, n_growing_modes=1, worst_growth_rate=3.12371
    - alpha=0.2, beta=0, stable=True, closed_form_stable=True, agrees=True, n_growing_modes=0, worst_growth_rate=0
- **the_channel_functions_are_monotone_on_the_imaginary_axis**: True
- **channel_run_along_the_imaginary_axis**: symmetric=[0.0230819, -3.1831], antisymmetric=[-0.0737426, -3.1831]
- **the_closed_form_agrees_with_every_probe**: True
- **grid_points**: 221
- **grid_stable**: 56
- **grid_mismatches**: 0
- **the_closed_form_matches_everywhere**: True
- **both_signs_are_represented**: True
- **what_this_shows**: positivity is a condition on the boundary data, separate from self-adjointness, and it has a closed form for the exchange-symmetric family

## T6_the_mouth_active_sector_is_rank_two — PASS

- **sector_counts**: below the free ground state=1, interlacing=12
- **roots_per_interlacing_gap**: 1=2, 2=2, 3=2, 4=2, 5=2, 6=2
- **two_per_interlacing_gap**: True
- **n_below_the_free_ground_state**: 1
- **there_is_a_sector_below_the_ground_state**: True
- **level_bookkeeping**:
    - level=0, degeneracy=1, mouth_active=1, untouched=0
    - level=1, degeneracy=4, mouth_active=2, untouched=2
    - level=2, degeneracy=9, mouth_active=2, untouched=7
    - level=3, degeneracy=16, mouth_active=2, untouched=14
    - level=4, degeneracy=25, mouth_active=2, untouched=23
- **at_most_two_modes_move_per_level**: True
- **untouched_modes_at_level_4**: 23
- **channel_endpoints**:
    - gap=1, lambda_range=[1, 4], symmetric_at_lower=-1.01321e+08, symmetric_at_upper=6.42122e+07, antisymmetric_at_lower=-0.0382805, antisymmetric_at_upper=3.71089e+07, symmetric_spans_the_whole_line=True, antisymmetric_spans_the_whole_line=False
    - gap=2, lambda_range=[4, 9], symmetric_at_lower=-6.42122e+07, symmetric_at_upper=3.86071e+07, antisymmetric_at_lower=-3.71089e+07, antisymmetric_at_upper=6.27141e+07, symmetric_spans_the_whole_line=True, antisymmetric_spans_the_whole_line=True
    - gap=3, lambda_range=[9, 16], symmetric_at_lower=-3.86071e+07, symmetric_at_upper=3.90483e+07, antisymmetric_at_lower=-6.27141e+07, antisymmetric_at_upper=6.22728e+07, symmetric_spans_the_whole_line=True, antisymmetric_spans_the_whole_line=True
    - gap=4, lambda_range=[16, 25], symmetric_at_lower=-3.90483e+07, symmetric_at_upper=5.29226e+07, antisymmetric_at_lower=-6.22729e+07, antisymmetric_at_upper=4.83985e+07, symmetric_spans_the_whole_line=True, antisymmetric_spans_the_whole_line=True
    - gap=5, lambda_range=[25, 36], symmetric_at_lower=-5.29226e+07, symmetric_at_upper=5.94106e+07, antisymmetric_at_lower=-4.83985e+07, antisymmetric_at_upper=4.19106e+07, symmetric_spans_the_whole_line=True, antisymmetric_spans_the_whole_line=True
    - gap=6, lambda_range=[36, 49], symmetric_at_lower=-5.94106e+07, symmetric_at_upper=5.30573e+07, antisymmetric_at_lower=-4.19106e+07, antisymmetric_at_upper=4.82639e+07, symmetric_spans_the_whole_line=True, antisymmetric_spans_the_whole_line=True
- **the_first_gap_antisymmetric_endpoints_are_finite**: True
- **the_symmetric_channel_does_span_it**: True
- **first_gap_root_depends_on_alpha_minus_beta**:
    - alpha_minus_beta=0.02, antisymmetric_root_in_gap_1=True
    - alpha_minus_beta=-0.05, antisymmetric_root_in_gap_1=False
    - alpha_minus_beta=-0.09, antisymmetric_root_in_gap_1=False
- **existence_in_the_first_gap_is_conditional**: True
- **by_channel**:
    - gap=1, symmetric=2.82167, antisymmetric=2.10694, n_roots=2
    - gap=2, symmetric=7.4169, antisymmetric=5.75125, n_roots=2
    - gap=3, symmetric=13.0976, antisymmetric=12.2953, n_roots=2
    - gap=4, symmetric=20.2927, antisymmetric=20.9599, n_roots=2
    - gap=5, symmetric=30.445, antisymmetric=30.8861, n_roots=2
    - gap=6, symmetric=43.1737, antisymmetric=42.0209, n_roots=2
- **the_scope**: det(C − BΓ) = 0 is the rank-two mouth-active sector; (n+1)² − 2 modes per level are untouched and are not in it

## T7_the_pr255_embedding — PASS

- **rows**:
    - kappa=0.3, omega=1.3+0.2j, det_of_the_embedding=0.989722-0.0198847j, pr255_one_minus_L=0.989722-0.0198847j, embedding_error=3.46945e-18, old_wrong_control=-0.143797-0.0914521j, old_control_error=1.13578, maximal=True, self_adjointness_defect=0.366421
    - kappa=0.3, omega=2.77-0.4j, det_of_the_embedding=0.988445+0.00167728j, pr255_one_minus_L=0.988445+0.00167728j, embedding_error=0, old_wrong_control=0.128977-0.222736j, old_control_error=0.888283, maximal=True, self_adjointness_defect=0.201096
    - kappa=0.3, omega=4.11+0.05j, det_of_the_embedding=1.01472-0.0655695j, pr255_one_minus_L=1.01472-0.0655695j, embedding_error=0, old_wrong_control=-0.164383-0.891036j, old_control_error=1.43934, maximal=True, self_adjointness_defect=0.315381
    - kappa=1, omega=1.3+0.2j, det_of_the_embedding=0.965739-0.0662823j, pr255_one_minus_L=0.965739-0.0662823j, embedding_error=0, old_wrong_control=-0.0444365-0.0302421j, old_control_error=1.01082, maximal=True, self_adjointness_defect=1.2214
    - kappa=1, omega=2.77-0.4j, det_of_the_embedding=0.961484+0.00559092j, pr255_one_minus_L=0.961484+0.00559092j, embedding_error=8.67362e-19, old_wrong_control=0.00135066-0.0660239j, old_control_error=0.962801, maximal=True, self_adjointness_defect=0.67032
    - kappa=1, omega=4.11+0.05j, det_of_the_embedding=1.04907-0.218565j, pr255_one_minus_L=1.04907-0.218565j, embedding_error=0, old_wrong_control=0.229421-0.629134j, old_control_error=0.916733, maximal=True, self_adjointness_defect=1.05127
- **worst_embedding_error**: 3.46945e-18
- **the_embedding_is_exact**: True
- **worst_old_control_error**: 1.43934
- **the_old_control_was_a_different_model**: True
- **every_embedding_is_maximal**: True
- **none_is_self_adjoint**: True
- **it_needs_a_singular_B**: q₁ = 0 is not of the form φ^reg = A q for any finite Hermitian A
- **what_is_not_concluded**: that PR #255's off-axis poles were *caused* by non-self-adjointness — a self-adjoint throat can also have growing modes, as the stability measurement shows, so this is a classification, not a diagnosis

## T8_the_phase_of_beta — PASS

- **rows**:
    - arg_beta=0, det=-0.0651079, det_conjugate=-0.0651079, conjugation_defect=0
    - arg_beta=0.136591, det=-0.0648291, det_conjugate=-0.0648291, conjugation_defect=0
    - arg_beta=0.273182, det=-0.0639978, det_conjugate=-0.0639978, conjugation_defect=0
    - arg_beta=0.409773, det=-0.0626295, det_conjugate=-0.0626295, conjugation_defect=0
    - arg_beta=0.546364, det=-0.0607498, det_conjugate=-0.0607498, conjugation_defect=0
    - arg_beta=0.682955, det=-0.0583936, det_conjugate=-0.0583936, conjugation_defect=0
    - arg_beta=0.819546, det=-0.0556048, det_conjugate=-0.0556048, conjugation_defect=0
    - arg_beta=0.956137, det=-0.0524354, det_conjugate=-0.0524354, conjugation_defect=0
    - arg_beta=1.09273, det=-0.0489444, det_conjugate=-0.0489444, conjugation_defect=0
    - arg_beta=1.22932, det=-0.0451968, det_conjugate=-0.0451968, conjugation_defect=0
    - arg_beta=1.36591, det=-0.0412625, det_conjugate=-0.0412625, conjugation_defect=0
    - arg_beta=1.5025, det=-0.0372147, det_conjugate=-0.0372147, conjugation_defect=0
    - arg_beta=1.63909, det=-0.0331289, det_conjugate=-0.0331289, conjugation_defect=0
    - arg_beta=1.77568, det=-0.0290811, det_conjugate=-0.0290811, conjugation_defect=0
    - arg_beta=1.91227, det=-0.0251468, det_conjugate=-0.0251468, conjugation_defect=0
    - arg_beta=2.04886, det=-0.0213992, det_conjugate=-0.0213992, conjugation_defect=0
    - arg_beta=2.18546, det=-0.0179082, det_conjugate=-0.0179082, conjugation_defect=0
    - arg_beta=2.32205, det=-0.0147388, det_conjugate=-0.0147388, conjugation_defect=0
    - arg_beta=2.45864, det=-0.01195, det_conjugate=-0.01195, conjugation_defect=0
    - arg_beta=2.59523, det=-0.00959378, det_conjugate=-0.00959378, conjugation_defect=0
    - arg_beta=2.73182, det=-0.00771403, det_conjugate=-0.00771403, conjugation_defect=0
    - arg_beta=2.86841, det=-0.00634577, det_conjugate=-0.00634577, conjugation_defect=0
    - arg_beta=3.005, det=-0.00551449, det_conjugate=-0.00551449, conjugation_defect=0
    - arg_beta=3.14159, det=-0.00523566, det_conjugate=-0.00523566, conjugation_defect=0
- **modulus**: 0.1655
- **spread_over_the_phase**: 0.0598723
- **the_phase_of_beta_is_physical**: True
- **worst_conjugation_defect**: 0
- **the_spectrum_is_conjugation_symmetric**: True
- **why_the_phase_survives**: Γ is not diagonal — the mouths are joined through the bulk as well as through the throat, and that fixes the relative phase
- **the_withdrawn_claim**: that a complex β makes the throat non-reciprocal; that read the Cayley entries, which depend on an arbitrary scale, as physical amplitudes

## T9_assessment — PASS

- **n_passed**: 8
- **n_total**: 8

## verdict — SELF_ADJOINTNESS_IS_CONSERVATION_NOT_STABILITY

SELF-ADJOINTNESS BUYS CONSERVATION AND A REAL lambda; POSITIVITY IS A SEPARATE CONDITION, AND IT IS THE ONE THAT DECIDES STABILITY. A point-supported throat is a self-adjoint extension of the Laplacian on S3 minus two points, parametrized by U(2), and writing the boundary condition as the PAIR B phi^reg = C q makes the mouth-active spectrum det(C - B Gamma) = 0. FIRST, IT IS DEFINABLE AT ALL: G(chi,omega) = sin(omega(pi-chi))/(4 pi sin chi sin(pi omega)), real on the axis, poles exactly at omega = n+1, finite at the antipode because the numerator's zero cancels sin chi, and agreeing with PR #255's branch series to 6.3e-12; its short-distance split is 1/(4 pi chi) + g(omega) + O(chi), remainder first order in chi, so the subtraction a point interaction needs is forced rather than chosen. SECOND, HERMITICITY IS EXACTLY FLUX CONSERVATION: the current through a small sphere at mouth j is Im(q_j* phi_j^reg), so the total absorbed is Im(q^dag A q), zero for EVERY q when A = A^dag -- 1.8e-16 over 200 random draws, against a median 0.54 for a non-Hermitian control. THIRD -- AND THIS IS A CORRECTION -- THE CAYLEY ENTRIES ARE NOT AMPLITUDES. The transform is unitary for every reference scale c, which is a real fact about the parametrization, but the same A gives entry magnitudes spread by 0.242 across c, so they are boundary-mixing coefficients; a closed universe has no asymptotic region to normalize a scattering matrix against. FOURTH -- AND THIS IS THE MAIN CORRECTION -- SELF-ADJOINTNESS MAKES lambda = omega^2 REAL AND NOTHING MORE. Gamma is real symmetric for real lambda of either sign, so the secular function is real; but lambda can be NEGATIVE, and then omega = +/- i sqrt|lambda| with one member of the pair growing. 2 of the three boundary matrices this module previously advertised do exactly that -- sigma = 2.470532 and 7.090982 -- and they were missed because the earlier search seeded only Re omega in [1.1, 6.9] and discarded roots leaving that window, so a root on the imaginary axis was outside its reach by construction. The claim 'real spectrum for every coupling, so a conserving throat cannot ring up' is WITHDRAWN. FIFTH, WHAT REPLACES IT IS A STABILITY REGION WITH A CLOSED FORM: both channel functions fall monotonically along the imaginary axis from their lambda = 0 values, so a growing mode exists iff alpha + beta < +0.02308202 or alpha - beta < -0.07374262, verified against a negative-lambda scan at all 221 grid points with 0 mismatches and only 56 of them stable. SIXTH, SCOPE: det(C - B Gamma) = 0 is the RANK-TWO MOUTH-ACTIVE SECTOR, not the spectrum -- level n has degeneracy (n+1)^2 and only two combinations can move, so 23 of 25 modes at level 4 never leave the free eigenvalue. Inside the sector there is a mode below the free ground state that an omega-scan starting above 1 cannot see, then two per interlacing gap; and the convenient claim that both channels run -infinity to +infinity across every gap is false, since the m = 1 pole cancels in the antisymmetric channel and a first-gap root is conditional on alpha - beta. SEVENTH, PR #255's RELATION EMBEDS EXACTLY, as q_1 = 0 with q_2 = gain . phi_1^reg, i.e. B = [[0,0],[gain,0]], C = I, giving det(C - B Gamma) = 1 - gain . G_d to 3.5e-18 against that round's own expression. It is maximal but not self-adjoint, and needs the singular B -- no finite Hermitian A reproduces it. The previous version of this control was a different function entirely, off by 1.44, so nothing is concluded from it: this is a classification of PR #255's boundary condition, NOT a diagnosis of its off-axis poles, because a self-adjoint throat can have growing modes too. EIGHTH, the phase of beta is physical because Gamma is not diagonal -- the mouths are joined through the bulk as well as the throat -- while the spectrum is invariant under beta -> conj(beta), which is time reversal; the earlier 'non-reciprocal' reading is withdrawn. WHAT IS STILL PUT IN: the boundary data, four real numbers chosen and not derived, with PR #249 still the thing that would fix them. The throat is POINT-supported, so it has no interior and no proper length, and the delay Delta of PRs #251-#255 does not survive into a point extension -- a real loss of structure, stated rather than hidden. No backreaction, no stress tensor, no topology change, no rate, and no two-source invariant.
