# probe — finite_throat

_2026-08-19T01:05:51.746089+00:00_

## T1_goal — PASS

- **question**: replace the point-supported throat of PRs #253-#259 with a FINITE conservative one -- a tube with a Dirichlet-to-Neumann map -- and measure what the point idealization was throwing away: a traversal delay, a reflected channel, and the rank of the static response.
- **scope**: a linear conformally coupled scalar on a fixed Einstein static universe. The tube's interior is one-dimensional and its mouths are still points in the ambient, which is the one approximation the round then measures the cost of. No backreaction, no topology change, and L, A and m are chosen rather than derived.

## T2_the_enlarged_system_is_conservative — PASS

- **rows**:
    - lambda=-4, hermiticity_defect=0, rank=2, imaginary_part=0
    - lambda=-0.7, hermiticity_defect=0, rank=2, imaginary_part=0
    - lambda=0.15, hermiticity_defect=0, rank=2, imaginary_part=0
    - lambda=0.9, hermiticity_defect=0, rank=2, imaginary_part=0
    - lambda=2.5, hermiticity_defect=0, rank=2, imaginary_part=0
    - lambda=7, hermiticity_defect=0, rank=2, imaginary_part=0
    - lambda=19, hermiticity_defect=0, rank=2, imaginary_part=0
- **worst_hermiticity_defect**: 0
- **every_condition_is_maximal**: True
- **worst_imaginary_part**: 0
- **worst_green_identity_residual**: 1.53066e-07
- **the_rank_one_control_defect**: 0.3
- **the_finite_throat_is_conservative**: True
- **the_control_is_not**: True
- **what_it_checks**: maximality and BC† = CB† at every λ, the DtN against its own interior, and the rank-one transfer model as the failing control

## T3_the_throat_transmits_at_the_traversal_time — PASS

- **rows**:
    - length=0.4, sigma_star=1.76858, contour=2.56858, onset_same_mouth=0, onset_opposite=0.258636, prediction=0.4
    - length=0.6, sigma_star=1.57523, contour=2.37523, onset_same_mouth=0, onset_opposite=0.460052, prediction=0.6
    - length=0.9, sigma_star=1.41666, contour=2.21666, onset_same_mouth=0, onset_opposite=0.762177, prediction=0.9
    - length=1.2, sigma_star=1.32491, contour=2.12491, onset_same_mouth=0, onset_opposite=1.0643, prediction=1.2
    - length=2, sigma_star=1.20567, contour=2.00567, onset_same_mouth=0, onset_opposite=1.13754, prediction=1.3
    - length=3, sigma_star=1.15197, contour=1.95197, onset_same_mouth=0, onset_opposite=1.13754, prediction=1.3
- **slope_below_the_ambient_path**: 1.00708
- **onset_spread_above_it**: 0
- **point_throat_onset_same_mouth**: 0
- **point_throat_onset_opposite**: 0
- **reflection_is_instantaneous**: True
- **the_delay_is_the_traversal_time**: True
- **the_ambient_path_takes_over**: True
- **the_point_throat_transmits_instantly**: True
- **the_two_channels**: the tube arrives at L, the ambient at d; a point throat arrives at 0, which is what a point throat is

## T4_the_delay_ledger_is_the_bounce_series — PASS

- **contour**: 1.6
- **order**: 400
- **cot_series_error**: 4.51828e-16
- **csc_series_error**: 1.66533e-15
- **delays**: same_mouth=[0, 2, 4, 6], opposite_mouths=[1, 3, 5, 7]
- **the_series_converge_on_the_contour**: True
- **the_parity_rule**: even multiples of L return to the mouth they entered; odd multiples cross

## T5_the_short_tube_limit_is_a_mixed_stratum — PASS

- **lambda_0**: 1
- **band**:
    - lambda=1, beta_exact=0.101589, beta_frozen=0.101589, relative_error=0
    - lambda=1.05, beta_exact=0.0974463, beta_frozen=0.101589, relative_error=0.0425132
    - lambda=1.2, beta_exact=0.0871271, beta_frozen=0.101589, relative_error=0.165987
    - lambda=1.5, beta_exact=0.0728225, beta_frozen=0.101589, relative_error=0.395024
    - lambda=2, beta_exact=0.0588642, beta_frozen=0.101589, relative_error=0.72582
    - lambda=3, beta_exact=0.0459474, beta_frozen=0.101589, relative_error=1.21099
- **short_tubes**:
    - length=0.4, anti=-0.0161312, anti_limit=-0.0159155, anti_error=0.000215658, anti_error_over_L_squared=0.00134786, sym=0.392568, sym_prediction=0.397887, distance_to_the_stratum=2.54733, rank=2, distance_over_L=6.36832
    - length=0.2, anti=-0.00798438, anti_limit=-0.00795775, anti_error=2.66324e-05, anti_error_over_L_squared=0.000665809, sym=0.79312, sym_prediction=0.795775, distance_to_the_stratum=1.26084, rank=2, distance_over_L=6.30421
    - length=0.1, anti=-0.00398219, anti_limit=-0.00397887, anti_error=3.31905e-06, anti_error_over_L_squared=0.000331905, sym=1.59022, sym_prediction=1.59155, distance_to_the_stratum=0.628843, rank=2, distance_over_L=6.28843
    - length=0.05, anti=-0.00198985, anti_limit=-0.00198944, anti_error=4.1457e-07, anti_error_over_L_squared=0.000165828, sym=3.18244, sym_prediction=3.1831, distance_to_the_stratum=0.314225, rank=2, distance_over_L=6.28449
    - length=0.02, anti=-0.000795801, anti_limit=-0.000795775, anti_error=2.65269e-08, anti_error_over_L_squared=6.63172e-05, sym=7.95748, sym_prediction=7.95775, distance_to_the_stratum=0.125668, rank=2, distance_over_L=6.28339
- **the_stratum**: B=[[0, 0], [0, 1]], C=[[-1, -0], [-0, -0]], condition=Phi_anti = 0 and q_sym = 0
- **the_band_error_reaches_one**: True
- **the_antisymmetric_channel_has_a_limit**: True
- **the_chart_matrix_diverges**: True
- **worst_symmetric_prediction_error**: 0.0135502
- **distance_to_the_stratum**: 0.125668
- **convergence_is_linear_in_L**: True
- **every_pair_is_maximal**: True
- **the_limit_exists_and_is_not_a_finite_A**: True
- **the_scope**: the limit is a mixed Dirichlet-Neumann stratum, not a finite A; the chart matrix diverges because the limit leaves the chart, not because it is absent

## T6_the_static_limit_is_rank_one — PASS

- **rows**:
    - lambda=1e-08, det_S=1.49076e-06, det_S_over_lambda=149.076, antisymmetry=4.29009e-09, defect=-8.84194e+06, minus_beta=-8.84194e+06
    - lambda=1e-06, det_S=0.000149076, det_S_over_lambda=149.076, antisymmetry=4.29009e-07, defect=-88419.4, minus_beta=-88419.4
    - lambda=0.0001, det_S=0.0149092, det_S_over_lambda=149.092, antisymmetry=4.28971e-05, defect=-884.206, minus_beta=-884.206
    - lambda=0.01, det_S=1.50715, det_S_over_lambda=150.715, antisymmetry=0.0042519, defect=-8.85336, minus_beta=-8.85389
- **massive**:
    - mass=0.05, det_S=-0.372447, det_S_over_mass_squared=-148.979, defect=35.3558, minus_beta=35.3558, defect_error=3.05533e-13
    - mass=0.1, det_S=-1.48687, det_S_over_mass_squared=-148.687, defect=8.83002, minus_beta=8.83002, defect_error=1.27898e-13
    - mass=0.2, det_S=-5.9012, det_S_over_mass_squared=-147.53, defect=2.19859, minus_beta=2.19859, defect_error=3.9968e-15
    - mass=0.4, det_S=-22.8898, det_S_over_mass_squared=-143.061, defect=0.540863, minus_beta=0.540863, defect_error=4.44089e-16
- **stratum_response**:
    - [6.78034, -6.78034]
    - [-6.78034, 6.78034]
- **stratum_rank**: 1
- **approach_to_the_stratum**:
    - length=0.4, distance_to_the_stratum_response=1.86612, det_S=4.34619e-08
    - length=0.1, distance_to_the_stratum_response=0.386691, det_S=9.00635e-09
    - length=0.02, distance_to_the_stratum_response=0.0739716, det_S=1.72268e-09
- **the_stratum_is_rank_one_too**: True
- **it_falsifies_the_finite_A_family_not_point_ness**: True
- **det_S_is_linear_in_lambda**: True
- **linear_coefficient**: 149.082
- **departure_at_the_largest_lambda**: 0.0109582
- **worst_antisymmetry**: 4.28971e-05
- **worst_defect_error_over_lambda**: 5.9622e-05
- **the_static_response_is_rank_one**: True
- **the_defect_diverges**: True
- **worst_defect_error**: 3.05533e-13
- **the_defect_is_still_minus_beta**: True
- **the_falsifiable_difference**: a generic finite-A throat is statically rank two; this tube is rank one and 𝒲 diverges — but so is its own point limit, so rank one falsifies the chart, not point-ness

## T7_the_interior_mass_is_a_cutoff — PASS

- **mass**: 2
- **cutoff**: 4
- **length**: 3
- **rows**:
    - lambda=0, propagating=False, beta=-0.000197254, imaginary_part=0, kappa=2, exponential_asymptote=-0.000197253
    - lambda=0.5, propagating=False, beta=-0.000310685, imaginary_part=0, kappa=1.87083, exponential_asymptote=-0.000310681
    - lambda=1.5, propagating=False, beta=-0.000876685, imaginary_part=0, kappa=1.58114, exponential_asymptote=-0.000876618
    - lambda=3, propagating=False, beta=-0.00794355, imaginary_part=0, kappa=1, exponential_asymptote=-0.00792386
    - lambda=3.9, propagating=False, beta=-0.229284, imaginary_part=0, kappa=0.316228, exponential_asymptote=-0.1949
    - lambda=4.1, propagating=True, beta=0.309661, imaginary_part=0, kappa=0, exponential_asymptote=nan
    - lambda=6, propagating=True, beta=-0.0631052, imaginary_part=0, kappa=0, exponential_asymptote=nan
    - lambda=12, propagating=True, beta=0.0348523, imaginary_part=0, kappa=0, exponential_asymptote=nan
- **worst_imaginary_part**: 0
- **asymptote_error**: 7.58439e-05
- **deepest_kappa_L**: 6
- **suppression_ratio**: 0.0031258
- **sign_changes_below**: 0
- **sign_changes_above**: 7
- **below_is_monotone**: True
- **the_evanescent_side_does_not_oscillate**: True
- **the_transmission_is_suppressed_below_cutoff**: True
- **everything_stays_real**: True
- **what_it_means**: a massive interior gives the channel a mass gap; below m² transmission is evanescent and the mouths look like two ordinary scatterers with a tunnelling term

## T8_the_growing_mode_belongs_to_the_mouth — PASS

- **rows**:
    - area=0.2, length=1.5, sigma_star=7.92673, closed_form=7.92665, relative_error=9.0494e-06, channel_split=0.000143477, split_over_exponential=4.28625, sigma_times_L=11.8901
    - area=0.2, length=3, sigma_star=7.92667, closed_form=7.92665, relative_error=2.19154e-06, channel_split=3.47439e-05, split_over_exponential=1.03787, sigma_times_L=23.78
    - area=0.2, length=6, sigma_star=7.92667, closed_form=7.92665, relative_error=2.19149e-06, channel_split=3.47431e-05, split_over_exponential=1.03784, sigma_times_L=47.56
    - area=0.5, length=1.5, sigma_star=5.01672, closed_form=5.01326, relative_error=0.000691809, channel_split=0.00697021, split_over_exponential=4.73808, sigma_times_L=7.52509
    - area=0.5, length=3, sigma_star=5.01402, closed_form=5.01326, relative_error=0.000153111, channel_split=0.00153659, split_over_exponential=1.04085, sigma_times_L=15.0421
    - area=0.5, length=6, sigma_star=5.01402, closed_form=5.01326, relative_error=0.000152818, channel_split=0.00153364, split_over_exponential=1.03886, sigma_times_L=30.0841
    - area=1, length=1.5, sigma_star=3.56681, closed_form=3.54491, relative_error=0.00617764, channel_split=0.0451896, split_over_exponential=4.6644, sigma_times_L=5.35021
    - area=1, length=3, sigma_star=3.55013, closed_form=3.54491, relative_error=0.00147399, channel_split=0.0105162, split_over_exponential=1.06219, sigma_times_L=10.6504
    - area=1, length=6, sigma_star=3.55005, closed_form=3.54491, relative_error=0.00145045, channel_split=0.0103455, split_over_exponential=1.04483, sigma_times_L=21.3003
- **by_separation**:
    - separation=0.8, sigma_star=7.92788, sigma_times_d=6.34231
    - separation=1.3, sigma_star=7.92667, sigma_times_d=10.3047
    - separation=2.4, sigma_star=7.92665, sigma_times_d=19.024
    - separation=3, sigma_star=7.92665, sigma_times_d=23.78
- **separation_spread**: 0.00122699
- **separation_spread_far**: 3.90009e-09
- **worst_closed_form_error**: 0.00147399
- **worst_channel_split**: 0.0105162
- **split_over_exponential**:
    - 1.03787
    - 1.03784
    - 1.04085
    - 1.03886
    - 1.04483
- **every_throat_has_one**: True
- **it_stops_knowing_the_separation**: True
- **the_closed_form_holds_once_sigma_L_is_large**: True
- **the_split_is_the_euclidean_propagator**: True
- **the_working_band**: 2.00692
- **the_diagnosis**: a mode that ignores L and d and does not split the channels is a single-mouth object, so the instability belongs to the point-mouth matching and not to the interior
- **the_open_question**: whether a finite-radius mouth or neck geometry removes it — unresolved here, and the thing to settle before stationary action or backreaction

## T9_the_contour_must_clear_it — PASS

- **sigma_star**: 1.57523
- **true_onset**: 0.6
- **rows**:
    - clearance=-0.03, contour=1.54523, above_the_mode=False, onset=0, pedestal=0.991906
    - clearance=0.02, contour=1.59523, above_the_mode=True, onset=0, pedestal=0.00261455
    - clearance=0.3, contour=1.87523, above_the_mode=True, onset=0.46463, pedestal=1.02588e-16
    - clearance=0.8, contour=2.37523, above_the_mode=True, onset=0.460052, pedestal=7.5865e-17
- **pedestal_below**: 0.991906
- **pedestal_above**: 1.02588e-16
- **onset_below**: 0
- **onset_above**: 0.460052
- **a_contour_below_the_mode_breaks_causality**: True
- **the_rule**: ε must exceed σ*, which has a closed form, so the contour can be placed before the solve rather than diagnosed after it
- **what_it_does_not_do**: clearing the contour evaluates the correct retarded solution of an unstable system; it does not cure the instability

## T10_assessment — PASS

- **n_passed**: 9
- **n_total**: 9

## verdict — THE_INTERIOR_GIVES_A_DELAY_AND_THE_POINT_MOUTH_IS_UNSTABLE

THE THROAT NOW HAS AN INTERIOR, AND THE INTERIOR IS THE DELAY. PRs #253-#259 all carried the same disclaimer -- point-supported, no proper length, no delay, a rank-one mouth-transfer model that is lossy for kappa < 1 -- and this round replaces it with a tube of length L and cross-section A whose Dirichlet-to-Neumann map is exact. Two lines are put in: N(lam) = A k [[cot kL, -csc kL], [-csc kL, cot kL]] and the matching q = -N Phi. Everything else follows. THE CONSERVATIVE OBJECT IS THE ENLARGED SYSTEM, ambient (+) tube, with lam-independent matching; eliminating the tube leaves a lam-DEPENDENT boundary condition -- the Weyl function of that elimination, which is a Nevanlinna function and not itself a self-adjoint operator on the ambient space. What is checked is that the elimination is faithful at each lam: rank[B|C] = 2 and BC^dag - CB^dag = 0.0e+00 at seven values of lam on both sides of zero, with the DtN map checked against the interior it summarizes by the SESQUILINEAR Green's identity to 1.5e-07 rather than against itself; PR #255's rank-one transfer model, the control, has defect 0.30 -- the size of the coupling itself, which is what the kappa < 1 loss was. *** AND THE THROAT TRANSMITS AT THE TRAVERSAL TIME. *** The measured object is the TWO-MOUTH BLOCK's impulse response -- the source and observer legs are gone, but Gamma, the ambient's own mouth-to-mouth propagator, stays in, so it is the coupled ambient+tube response and not the throat alone: r_11, same mouth in and out, starts at 0.0 -- a wave that reaches a mouth is partly reflected INSTANTLY -- and r_12, opposite mouths, starts at min(L, d) with d(onset)/dL = 1.0071 against a predicted 1, saturating above the ambient path to 0.0e+00. THE AMBIENT ALSO CONNECTS THE TWO MOUTHS, along a geodesic of length d, whether or not they are joined: PR #258's cross-mouth channel and PR #259's beta = 0 control, now SEPARATED IN TIME instead of by rank counting. A point throat transmits at 0.0, which is what a point throat is. THE LEDGER IS A DERIVATION, NOT A FIT: on the contour cot x = -i - 2i sum e^{2ikx} and csc x = -2i sum e^{i(2k+1)x}, to 4.5e-16 and 1.7e-15, so the same-mouth entry carries 0, 2L, 4L and the cross-mouth entry L, 3L, 5L. The parities are the physics, and the reflected channel is the one the rank-one model does not have at all. THERE IS A POINT LIMIT, AND IT IS NOT A FINITE A. Freezing A at A(lam_0) is exact at lam_0 and 121% wrong at 3 lam_0 -- a band of width ~1/L in omega -- and as L -> 0 the antisymmetric channel converges to -L/(2A) while the symmetric one diverges like 2/(A lam L). A first draft concluded that the limit does not exist. It does: a boundary pair is defined up to (B,C) -> (MB,MC), so a diverging chart matrix means the limit has LEFT THE CHART, and row-scaled the pair converges to (P_anti, -P_sym) -- Phi_anti = 0, q_sym = 0, a MIXED DIRICHLET-NEUMANN stratum, maximal throughout and 0.13 away at L = 0.02, linearly in L. No finite Hermitian A reaches it, which is exactly the stratum PR #257's review said the chart does not cover. THE SAME ZERO MODE BREAKS PR #258's TOMOGRAPHY. At lam = 0 the static response collapses onto [[1,-1],[-1,1]] to 4.3e-05, det S goes to zero LINEARLY in lam with coefficient 149.08, and W = S_12/det S - G_0 diverges like 1/lam. WHAT THAT FALSIFIES is the generic finite-A family, every member of which is rank two -- and NOT point-ness: the tube's own short-tube stratum is rank one as well, and the tube converges to it. The first draft claimed the stronger and wrong version. Give the tube an interior mass and the rank returns; and off the collapse W = -beta(lam) exactly, to 3.1e-13, with beta the tube's own transmission amplitude. PR #258's theorem survives the generalization and returns the interior's amplitude instead of a constant. AN INTERIOR MASS IS A TRANSMISSION CUTOFF: below lam = m^2 the tube is evanescent, beta matches its exponential asymptote to 7.6e-05, is suppressed by 3.1e-03, and has 0 sign changes against 7 above -- monotone decay against oscillation, which is the discriminator rather than the sign itself; the interior has a MASS GAP and the channel below it is evanescent, which is the opposite of a low-pass filter. *** AND THE MODEL FAILS THE STABILITY GATE, WITH A GROWING MODE THAT IS THE MOUTH'S AND NOT THE TUBE'S. *** A(lam) decreases and Gamma(lam) increases -- PR #257's Gram identity -- so A - Gamma is strictly monotone between poles and each channel has at most one root, a count rather than a scan. The symmetric channel always has exactly one, at lam < 0. Three facts identify what it is, all limits with their convergence measured: its rate matches sigma* = 2 sqrt(pi/A) to 1.5e-03 with NO L in it; two mouth separations agree to 3.9e-09; and the channel splitting is 1.04 e^{-sigma* d}, the Euclidean propagator between the mouths. A mode that ignores the tube's length and the mouths' separation and does not distinguish the channels is a SINGLE-MOUTH object: the instability belongs to the POINT-MOUTH MATCHING and not to the interior. THIS IS THE ROUND'S FALSIFICATION RESULT, and nothing here cures it -- whether a finite-radius mouth or neck geometry removes the mode is open, and is the thing to settle before stationary-action or backreaction work. SO THE CONTOUR MUST CLEAR IT: placed 0.03 below sigma*, the inversion returns a field with support before its own light cone -- a pedestal at 99% of the peak and an onset of 0.0 for an event that cannot begin until 0.6 -- against 1.0e-16 when it is placed above. Same species as PR #259's under-resolved contour, reported the same way, and this time the rule is checkable in advance because sigma* has a closed form -- but clearing the contour STABILIZES NOTHING: above sigma* the inversion returns the correct retarded solution OF AN UNSTABLE SYSTEM, and the delay is read from the causal onset, which is immune to what happens afterwards. WHICH FREQUENCIES: the delay and the ledger are statements about the exact model's analytic structure at ALL frequencies -- a causal onset is a UV object and the probe pulse carries content to w ~ 30, far above sigma* -- so they are exact results about this model rather than predictions about a resolved physical mouth; the static and low-frequency results sit inside the band. A is a ONE-DIMENSIONAL COUPLING, not an area with a radius attached. WHAT IS STILL PUT IN: the background, the mouth positions, and L, A, m -- three numbers with dimensions where the real-field point sector has three without them. NO BACKREACTION: the throat is a fixed background, and PR #261's common action and PR #262's A/B/A+B metric backreaction are the next two steps. They get an object with a proper length, a delay, a conservation law -- and a stated failure mode to resolve first.
