# probe — two_wave

_2026-08-18T02:48:43.927292+00:00_

## T1_goal — PASS

- **question**: build the time-dependent two-wave invariant on the stable point-throat background and measure it against the known WKB collinear/head-on result -- not to re-derive that identity, but to quantify how the exact multipath throat-coupled field differs from it.
- **scope**: a linear conformally coupled scalar on a fixed Einstein static universe. The throat is point-supported and taken strictly inside PR #257's cone with its Loewner margin quoted. No backreaction: the stress tensor is computed from the field and never fed back. No topology change, no rate, and the boundary data is chosen rather than derived.

## T2_the_solver_is_the_image_sum — PASS

- **chi**: 1.7
- **rows**:
    - t=1.70288, solver=0.532946, image_sum=0.532946, difference=2.22045e-16, maslov_sign=1
    - t=4.58679, solver=-0.532598, image_sum=-0.532598, difference=-2.22045e-16, maslov_sign=-1
    - t=7.9834, solver=0.533558, image_sum=0.533558, difference=3.33067e-16, maslov_sign=1
    - t=10.8673, solver=-0.533496, image_sum=-0.533496, difference=0, maslov_sign=-1
- **worst_difference**: 3.33067e-16
- **peak_scale**: 0.533558
- **the_two_constructions_agree**: True
- **the_signs_alternate**: True
- **what_this_licenses**: the frequency-domain solve, on which every stress tensor below is built

## T3_the_wave_equation_holds — PASS

- **free**: worst_residual=4.26326e-14, scale=41.9358, relative=1.01661e-15
- **with_throat**: worst_residual=4.26326e-14, scale=41.9433, relative=1.01643e-15
- **loewner_margin**: 0.322775
- **the_equation_holds**: True
- **nothing_is_differenced**: True

## T4_the_stress_tensor_is_traceless — PASS

- **rows**:
    - throat=False, T00=1.65027, trace=-3.08087e-15, relative=1.86689e-15
    - throat=True, T00=1.65027, trace=-2.83107e-15, relative=1.71552e-15
- **worst_relative_trace**: 1.86689e-15
- **the_tensor_is_traceless**: True
- **why_it_is_not_vacuous**: □φ is taken from the solve, not substituted on shell; the trace then equals φ(□φ − φ)

## T5_the_wkb_limit_is_recovered — PASS

- **rows**:
    - carrier=6, collinear=0.000199441, collinear_wkb=1.97215e-31, collinear_dot=1, head_on=3.95291, head_on_wkb=4, head_on_dot=-1
    - carrier=8, collinear=0.000112649, collinear_wkb=1.97215e-31, collinear_dot=1, head_on=3.96252, head_on_wkb=4, head_on_dot=-1, collinear_slope=-1.98565
    - carrier=12, collinear=2.85956e-05, collinear_wkb=1.97215e-31, collinear_dot=1, head_on=3.98171, head_on_wkb=4, head_on_dot=-1, collinear_slope=-3.38137
    - carrier=16, collinear=5.62625e-06, collinear_wkb=1.97215e-31, collinear_dot=1, head_on=3.99083, head_on_wkb=4, head_on_dot=-1, collinear_slope=-5.65141
    - carrier=24, collinear=1.9048e-07, collinear_wkb=1.97215e-31, collinear_dot=1, head_on=3.99817, head_on_wkb=4, head_on_dot=-1, collinear_slope=-8.35005
    - carrier=32, collinear=1.38239e-08, collinear_wkb=1.97215e-31, collinear_dot=1, head_on=3.99956, head_on_wkb=4, head_on_dot=-1, collinear_slope=-9.11821
    - carrier=48, collinear=1.78896e-10, collinear_wkb=1.97215e-31, collinear_dot=1, head_on=3.99995, head_on_wkb=4, head_on_dot=-1, collinear_slope=-10.7219
- **width**: 0.1
- **head_on_at_the_largest_carrier**: 3.99995
- **head_on_target**: 4
- **head_on_error**: 5.02828e-05
- **collinear_at_the_largest_carrier**: 1.78896e-10
- **the_directions_are_exactly_parallel**: True
- **the_directions_are_exactly_antiparallel**: True
- **head_on_converges_to_four**: True
- **collinear_converges_to_zero**: True
- **the_collinear_exponent_steepens**: True
- **under_resolved_contour_value**: -0.000673358
- **converged_value_there**: 1.38239e-08
- **the_contour_needs_eps_above_the_spacing**: True

## T6_multipath_destroys_the_collinear_null — PASS

- **carrier**: 24
- **observer_reach**: 0.8
- **loewner_margin**: 0.322775
- **rows**:
    - branch_pair=A direct + B direct, invariant=1.9048e-07, geometric_prediction=1.97215e-31, energy_product=16.5274
    - branch_pair=A long-way image + B direct, invariant=3.99806, geometric_prediction=4, energy_product=16.8858
    - branch_pair=A direct + B via a mouth, invariant=0.56501, geometric_prediction=0.563669, energy_product=0.0124955
- **throat_delay**: 3.00087
- **throat_exit_mouth**: 2
- **collinear_floor**: 1.9048e-07
- **long_way_value**: 3.99806
- **through_the_throat_value**: 0.56501
- **through_the_throat_prediction**: 0.563669
- **throat_relative_error**: 0.00237922
- **free_control_energy_product**: 4.1489e-29
- **throat_energy_product**: 0.0124955
- **the_control_has_no_second_arrival**: True
- **the_direct_pair_is_null**: True
- **the_winding_image_reads_head_on**: True
- **the_throat_matches_its_geometry**: True
- **the_lesson**: the two-wave invariant is branch-resolved: the same sources at the same event give 0, 4 or the mouth's angle depending on which branch arrived
- **what_the_control_does_not_say**: that the mouths' CONNECTION supplies the branch — β = 0 gives the same invariant; see the cross-mouth audit

## T7_the_cross_mouth_audit — PASS

- **carrier**: 60
- **width**: 0.035
- **channels**:
    - exit_mouth=2, entry_mouth=2, delay=3.00087, predicted_invariant=0.563669
    - exit_mouth=2, entry_mouth=1, delay=3.20923, predicted_invariant=0.563669
    - exit_mouth=1, entry_mouth=2, delay=3.23688, predicted_invariant=0.651935
    - exit_mouth=1, entry_mouth=1, delay=3.44523, predicted_invariant=0.651935
- **distinct_predictions**:
    - 0.563669
    - 0.651935
- **rows**:
    - exit_mouth=2, entry_mouth=2, delay=3.00087, predicted_invariant=0.563669, measured_invariant=0.563951, control_beta_zero=0.563951, relative_error=0.00050065, beta_shift=3.31048e-08, weight=0.0915051, control_weight=0.0915346
    - exit_mouth=1, entry_mouth=1, delay=3.44523, predicted_invariant=0.651935, measured_invariant=0.651408, control_beta_zero=0.651332, relative_error=0.000806969, beta_shift=7.60463e-05, weight=0.0954324, control_weight=0.0977657
- **beta_sweep**:
    - beta=0, loewner_margin=0.295844, invariant=0.563951, weight_ratio=1
    - beta=0.06, loewner_margin=0.322775, invariant=0.563951, weight_ratio=0.999678
    - beta=0.12, loewner_margin=0.274503, invariant=0.563951, weight_ratio=0.998709
    - beta=0.2, loewner_margin=0.196695, invariant=0.563951, weight_ratio=0.996409
    - beta=0.26, loewner_margin=0.137271, invariant=0.563951, weight_ratio=0.993921
- **worst_relative_error**: 0.000806969
- **the_prediction_depends_only_on_the_exit_mouth**: True
- **the_field_picks_the_right_one**: True
- **worst_beta_shift**: 7.60463e-05
- **beta_sweep_spread**: 6.22376e-07
- **exit_mouth_separation**: 0.0882653
- **beta_spread_over_the_signal**: 7.0512e-06
- **the_invariant_is_beta_independent**: True
- **the_weight_moves_instead**: 0.0060785
- **every_sweep_point_is_inside_the_cone**: True
- **the_scope**: this observable sees structure at the mouths, not the connection between them; W = −β from the low-frequency limit is what sees the connection

## T8_the_interference_tensor — PASS

- **carrier**: 24
- **width**: 0.1
- **rows**:
    - configuration=collinear, invariant=1.9048e-07, delta_T00=8.12939, normalized_delta_T00=1.99966, max_component=8.12939, trace=1.77636e-15, with_a_source_removed=0
    - configuration=head_on, invariant=3.99817, delta_T00=4.03699, normalized_delta_T00=1.04402, max_component=4.03699, trace=-4.44089e-16, with_a_source_removed=0
- **collinear_invariant**: 1.9048e-07
- **collinear_interference**: 1.99966
- **head_on_invariant**: 3.99817
- **head_on_interference**: 1.04402
- **worst_trace**: 1.77636e-15
- **worst_value_with_a_source_removed**: 0
- **delta_T_is_traceless**: True
- **delta_T_vanishes_when_a_source_is_removed**: True
- **the_interference_is_maximal_where_the_invariant_is_null**: True
- **the_lesson**: ΔT and T_A:T_B are different diagnostics: the collinear configuration nulls the invariant and MAXIMISES the interference energy

## T9_the_arrivals_are_the_branch_ledger — PASS

- **free_arrivals**:
    - 2.4
    - 3.88319
    - 8.68319
- **two_leg_times**:
    - 2.84451
    - 2.97822
    - 3.08052
    - 3.21422
- **rows**:
    - branch=free, predicted_t=2.4, found_t=2.39868, offset=-0.00131836, value=0.783145
    - branch=free, predicted_t=3.88319, found_t=3.88184, offset=-0.00134937, value=-0.783136
    - branch=free, predicted_t=8.68319, found_t=8.68378, offset=0.000591548, value=0.783296
- **worst_free_offset**: 0.00134937
- **free_signs**:
    - 1
    - -1
    - 1
- **the_free_signs_alternate**: True
- **the_free_arrivals_are_sharp**: True
- **throat_onset_measured**: 2.88391
- **throat_onset_predicted**: 2.84451
- **onset_offset**: 0.0393977
- **the_throat_onset_is_causal**: True
- **the_throat_arrivals_are_new**: True
- **why_onsets_and_not_peaks**: R(ω) has poles, so a throat arrival rings up rather than pulsing; the onset is the sharp part

## T10_the_only_tail_is_the_throats — PASS

- **free**: peak=0.117775, between_arrivals=1.66007e-09, ratio=1.40953e-08
- **with_throat**: peak=0.117775, between_arrivals=0.00952193, ratio=0.0808483
- **free_ratio**: 1.40953e-08
- **throat_ratio**: 0.0808483
- **amplification**: 5.73585e+06
- **the_free_field_has_no_tail**: True
- **the_throat_has_one**: True
- **why**: S³ × R is conformally flat, so the conformal scalar obeys Huygens exactly; R(ω) has poles, so the mouths ring

## T11_the_caustic_is_where_wkb_stops — PASS

- **rows**:
    - carrier=8.5, omega_times_e=2, e=0.235294, exact=0.310384, wkb=0.341345, ratio=0.909297
    - carrier=8.5, omega_times_e=1, e=0.117647, exact=0.570493, wkb=0.677971, ratio=0.841471
    - carrier=8.5, omega_times_e=0.5, e=0.0588235, exact=0.648949, wkb=1.3536, ratio=0.479426
    - carrier=8.5, omega_times_e=0, e=0, exact=0.676409, wkb=inf, ratio=0
    - carrier=16.5, omega_times_e=2, e=0.121212, exact=0.598431, wkb=0.658125, ratio=0.909297
    - carrier=16.5, omega_times_e=1, e=0.0606061, exact=1.10555, wkb=1.31383, ratio=0.841471
    - carrier=16.5, omega_times_e=0.5, e=0.030303, exact=1.25919, wkb=2.62646, ratio=0.479426
    - carrier=16.5, omega_times_e=0, e=0, exact=1.31303, wkb=inf, ratio=0
    - carrier=32.5, omega_times_e=2, e=0.0615385, exact=1.17659, wkb=1.29395, ratio=0.909297
    - carrier=32.5, omega_times_e=1, e=0.0307692, exact=2.17661, wkb=2.58668, ratio=0.841471
    - carrier=32.5, omega_times_e=0.5, e=0.0153846, exact=2.47994, wkb=5.17274, ratio=0.479426
    - carrier=32.5, omega_times_e=0, e=0, exact=2.58627, wkb=inf, ratio=0
- **saturation_amplitudes**:
    - 0.676409
    - 1.31303
    - 2.58627
- **carrier_ratios**:
    - 1.94118
    - 1.9697
- **saturation_ratios**:
    - 1.94118
    - 1.9697
- **worst_collapse_spread**: 6.55032e-15
- **the_saturation_is_linear_in_omega**: True
- **the_ratio_collapses_in_omega_times_e**: True
- **the_caustic_scale**: e* ~ 1/ω
- **what_wkb_gets_wrong**: a divergence where the exact amplitude is finite and proportional to ω

## T12_the_low_frequency_limit_recovers_the_tomography — PASS

- **n_observations**: 12
- **rows**:
    - raw_eps_0p08=0.0564312, raw_eps_0p04=0.0568025, richardson=0.0569263, exact=0.0569273, error=-1.01923e-06
    - raw_eps_0p08=0.0598284, raw_eps_0p04=0.0601853, richardson=0.0603042, exact=0.0603052, error=-9.52126e-07
    - raw_eps_0p08=0.0395206, raw_eps_0p04=0.039847, richardson=0.0399558, exact=0.0399568, error=-9.57652e-07
    - raw_eps_0p08=0.0466189, raw_eps_0p04=0.0469672, richardson=0.0470833, exact=0.0470843, error=-9.91187e-07
- **worst_kernel_error**: 1.54185e-06
- **worst_response_error**: 0.000276675
- **W_from_the_time_dependent_solve**: -0.06001
- **minus_beta**: -0.06
- **W_error**: 1.00136e-05
- **the_bridge_closes**: True
- **what_it_checks**: the DC content of the contour integral, the least-squares recovery, and PR #258's defect, end to end

## T13_assessment — PASS

- **n_passed**: 12
- **n_total**: 12

## verdict — THE_TWO_WAVE_INVARIANT_IS_BRANCH_RESOLVED

THE EXACT TWO-WAVE INVARIANT IS BRANCH-RESOLVED, AND THAT IS WHERE IT DEPARTS FROM WKB. The known limit is not re-derived here; it is used as the control. Two sources and an observer on one great circle give arriving directions that are exactly parallel or exactly antiparallel -- to 1e-12, by construction rather than by fitting -- so WKB's N = (1 - cos theta)^2 predicts 0 and 4, and the exact field returns 3.999950 head-on (error 5.0e-05) and 1.8e-10 collinear. THE COLLINEAR NULL IS STRONGER THAN LEADING ORDER: on this geometry the two wavefronts share their normal exactly, so amplitude gradients cannot tilt either k, the residue falls faster than any fixed power, and the measured exponent steepens with omega. *** AND THAT IS WHAT MAKES THE MULTIPATH RESULT LARGE. *** Holding the sources and the observation point fixed and changing only WHICH BRANCH has arrived, the same pair gives 1.9e-07 on the direct branches, 3.9981 when A arrives on its LONG-WAY WINDING IMAGE -- that branch runs the other way round the sphere, so its arrival direction is reversed and a collinear pair reads head-on -- and 0.5650 when B arrives VIA A MOUTH, against 0.5637 predicted from the mouth's position alone, a 0.24% agreement with a number that was never fitted. The free-propagation control at the same instant has no second arrival at all -- energy product 4.1e-29 against 1.2e-02 -- so the mouths are CREATING the second arrival, not merely bending one. AND THE (i,j) AUDIT SCOPES THAT. All four two-leg paths are enumerated rather than minimised over -- j the mouth the source drives, i the mouth the signal leaves from -- and the predicted invariant depends ONLY ON i, giving two distinct predictions 0.563669 and 0.651935, both matched by the field to 8.1e-04 relative. Then the control PR #258's review taught this arc to run first: beta = 0, two DISCONNECTED mouths. The invariant does not move -- swept over beta in [0, 0.26], all inside the cone, N shifts by 6.2e-07, which is 7.1e-06 of the 0.0883 that separates the two exit mouths, while the channel's WEIGHT moves by 0.6%. So the honest scope is the dynamical version of #258's: THIS OBSERVABLE SEES STRUCTURE AT THE MOUTHS, NOT THE CONNECTION BETWEEN THEM. What sees the connection is W = -beta, from the low-frequency limit of the same solve. AND T[phi_A + phi_B] IS A DIFFERENT DIAGNOSTIC ENTIRELY. T is quadratic, so the two-wave content of the TOTAL stress tensor is the bilinear cross term dT = T[phi_A+phi_B] - T[phi_A] - T[phi_B], built from three evaluations of the same functional: traceless to 1.8e-15 and EXACTLY zero when either source is switched off, which is PR #253's missing property at tensor level. Normalised, dT^00/sqrt(T_A^00 T_B^00) reaches its MAXIMUM 2.000 in the COLLINEAR configuration -- two parallel waves add coherently -- precisely where T_A:T_B vanishes; head-on it is only 1.044. A backreaction estimate driven by C = T_A:T_B would look at the collinear case, see nothing, and be wrong about its own source by the size of the whole effect. THE CONCLUSION: the collinear null is not spoiled by curvature corrections, which are 1e-7 here; it is spoiled by MULTIPATH, at O(1). The invariant has to carry the branch index PR #255 named, and a single-branch WKB formula is not merely approximate on this background -- it is answering a different question. THE SOLVER EARNS THAT. Its free part reproduces PR #254's closed-form winding-image sum to 3.3e-16 including the alternating Maslov signs, two constructions sharing no code; the solved field satisfies the conformal wave equation to 4e-16 relative with and without the throat, with nothing finite-differenced; and the improved stress tensor is traceless to 1.9e-15 -- which is a real test because box phi is taken from the solve rather than substituted on shell, so the trace equals phi(box phi - phi) instead of vanishing algebraically. Convergence is measured too: with the contour offset at the frequency spacing the collinear value comes out -6.7e-04 instead of 1.4e-08, four orders wrong, so eps must sit well above 2pi/span. THE OTHER CORRECTIONS, QUANTIFIED. Arrivals: the free ones land at chi, 2pi-chi, 2pi+chi to 1.3e-03 with signs + - +, and the throat adds two-leg arrivals the free ledger does not contain, checked at the causal onset because R(w) has poles and a throat arrival rings up rather than pulsing. TAIL: S3 x R is conformally flat, so the conformal scalar obeys Huygens exactly -- between geometric arrivals the free field is 1.4e-08 of its peak against 8.1e-02 with the throat, a factor of 5.7e+06. Every bit of tail in this model is the throat's. CAUSTIC: geometric optics gives 1/(4 pi sin chi), divergent at the antipode, where the exact kernel is finite and LINEAR in omega; in between the exact/WKB ratio is |sin(omega e)|, a function of omega*e alone, identical across carriers to 6.6e-15, so the caustic is cut off at e* ~ 1/omega. AND THE ROUND CLOSES BACK ON PR #258: the DC content of the solved time series is exactly the static kernel that round did its tomography on, and running that protocol on numbers from the dynamic solver returns W = -0.060010 against -beta = -0.06, error 1.0e-05, through the whole contour integral with the O(eps^2) contour bias Richardson-extrapolated and both numbers reported. WHAT IS STILL PUT IN: the background, the mouth positions, and the boundary data -- four real numbers chosen and not derived, with PR #249 still the thing that would fix them from matter. NO BACKREACTION: the stress tensor is computed from the field and never fed back, which is the next step, and it now has a concrete object to feed.
