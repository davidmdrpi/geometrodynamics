# probe — wormhole_ledger

_2026-08-13T05:44:38.155503+00:00_

## T1_goal — PASS

- **question**: previous rounds visualised wormhole networks whose mouths supply antipodal wavefronts; can the emitted shell, the passing collapsing shell, the receiver's recoil and a time-delayed past mouth be shown to be ONE conserved wave rather than several that happen to balance?
- **scope**: kinematics and bookkeeping on a fixed round S3. The wormhole is an IDENTIFICATION MAP, not a solution; kappa and the delay are parameters; shells.junction (PR #249) priced the throat and this round inherits that bill without paying it. No Einstein equation, no backreaction, no cross-section.

## T2_the_destination_is_geometry_not_staging — PASS

- **rows**:
    - chi=0.4, sampled_volume=0.266923, closed_form_volume=0.259633, relative_error=0.0280779, standard_error=0.00360472, z_score=2.02233, area_at_chi=1.90565
    - chi=0.9, sampled_volume=2.60143, closed_form_volume=2.59543, relative_error=0.00231025, standard_error=0.0105573, z_score=0.567958, area_at_chi=7.71074
    - chi=1.5708, sampled_volume=9.87178, closed_form_volume=9.8696, relative_error=0.00022, standard_error=0.0156052, z_score=0.13914, area_at_chi=12.5664
    - chi=2.2, sampled_volume=16.8119, closed_form_volume=16.8126, relative_error=3.98288e-05, standard_error=0.0110921, z_score=0.0603694, area_at_chi=8.21421
    - chi=2.8, sampled_volume=19.5705, closed_form_volume=19.5761, relative_error=0.000284236, standard_error=0.00287271, z_score=1.93692, area_at_chi=1.41016
- **worst_z_score**: 2.02233
- **worst_relative_error**: 0.0280779
- **derivative_of_volume_is_the_area**: True
- **the_area_law_holds**: True
- **tested_against_the_sampling_error_not_a_fixed_percent**: True
- **total_volume_is_two_pi_squared**: True
- **area_vanishes_at_both_poles**: True
- **amplitude_diverges_at_the_antipode**: True
- **so_the_destination_is_geometry_not_staging**: True

## T3_the_arriving_shell_is_actually_collapsing — PASS

- **refocused_receiver**: arrival_chi=3.14159, area_at_arrival=1.88465e-31, amplitude_at_arrival=1e+09, dA_dchi_min=-12.5662, dA_dchi_max=-0.158703, collapsing_all_the_way_in=True, expanding_all_the_way_in=False
- **displaced_receiver**: arrival_chi=1.2, area_at_arrival=10.9164, amplitude_at_arrival=1.07292, dA_dchi_min=8.53272, dA_dchi_max=12.5664, collapsing_all_the_way_in=False, expanding_all_the_way_in=True
- **default_arrival_is_the_antipode**: True
- **displaced_arrival_chi**: 1.2
- **the_shell_collapses_only_at_the_refocus**: True
- **the_focused_arrival_has_vanishing_area**: True
- **the_displaced_arrival_does_not**: True
- **so_collapsing_is_a_sign_not_a_caption**: True

## T4_the_drawn_object_is_the_shell_and_its_size_is_the_area — PASS

- **worst_geodesic_distance_error**: 6.66134e-15
- **the_drawn_points_are_the_level_set**: True
- **worst_shadow_radius**: 1
- **the_shadow_never_leaves_the_unit_ball**: True
- **rows**:
    - chi=0.2, sin_chi=0.198669, screen_extent=0.365652, sqrt_area_over_4pi=0.198669, extent_over_sin_chi=1.84051
    - chi=0.6, sin_chi=0.564642, screen_extent=1.03923, sqrt_area_over_4pi=0.564642, extent_over_sin_chi=1.84051
    - chi=1.5708, sin_chi=1, screen_extent=1.84051, sqrt_area_over_4pi=1, extent_over_sin_chi=1.84051
    - chi=2.1, sin_chi=0.863209, screen_extent=1.58874, sqrt_area_over_4pi=0.863209, extent_over_sin_chi=1.84051
    - chi=2.9, sin_chi=0.239249, screen_extent=0.44034, sqrt_area_over_4pi=0.239249, extent_over_sin_chi=1.84051
- **extent_ratio_spread**: 3.61929e-16
- **the_drawn_size_is_sqrt_of_the_area**: True
- **stereographic_comparison**:
    - epsilon=0.01, stereographic_radius=199.998, radius_times_epsilon=1.99998, shadow_radius=0.562634
    - epsilon=0.001, stereographic_radius=2000, radius_times_epsilon=2, shadow_radius=0.560759
    - epsilon=0.0001, stereographic_radius=20000, radius_times_epsilon=2, shadow_radius=0.560577
    - epsilon=1e-05, stereographic_radius=200000, radius_times_epsilon=2, shadow_radius=0.560559
- **shell_reaches_the_pole_at_chi**: 1.02896
- **the_stereographic_radius_diverges_as_two_over_epsilon**: True
- **it_does_not_converge_under_refinement**: True
- **while_the_shadow_stays_in_the_unit_ball**: True
- **and_the_shell_reaches_every_pole**: True
- **so_no_stereographic_pole_is_safe_for_this_scene**: a shell launched from a point sweeps all of S³, so whatever pole is chosen the shell crosses it once; the first draft chose q₃=+1, which is the emitter's own position, and the emitter was a division by zero that never got drawn
- **and_the_shadow_is_two_to_one**: depth is returned, not hidden; a crossing in the image is not a crossing on S³, and no claim in this round rests on one

## T5_one_self_consistent_amplitude — PASS

- **rows**:
    - kappa=0.3, closed_form=(1.4285714285714286+0j), iterated=(1.4285714285714284+0j), abs_error=2.22045e-16, amplification=1.42857
    - kappa=-0.6, closed_form=(0.625+0j), iterated=(0.625+0j), abs_error=0, amplification=0.625
    - kappa=0.9, closed_form=(10.000000000000002+0j), iterated=(9.999999999999995+0j), abs_error=7.10543e-15, amplification=10
    - kappa=0.99, closed_form=(99.99999999999991+0j), iterated=(99.9999999999992+0j), abs_error=7.10543e-13, amplification=100
    - kappa=(0.3+0.4j), closed_form=(1.0769230769230769+0.6153846153846154j), iterated=(1.0769230769230769+0.6153846153846154j), abs_error=0, amplification=1.24035
- **worst_abs_error**: 7.10543e-13
- **the_fixed_point_is_what_the_loop_does**: True
- **no_tuning_was_required**: True
- **what_it_means**: a linear wave on a time-displaced loop has exactly one self-consistent amplitude; the loop does not get to choose and nothing has to be arranged

## T6_only_the_resonance_obstructs — PASS

- **rows**:
    - abs_kappa=0.2, iteration_converges=True, fixed_point_exists=True
    - abs_kappa=0.5, iteration_converges=True, fixed_point_exists=True
    - abs_kappa=0.9, iteration_converges=True, fixed_point_exists=True
    - abs_kappa=0.99, iteration_converges=True, fixed_point_exists=True
    - abs_kappa=1.01, iteration_converges=False, fixed_point_exists=True
    - abs_kappa=1.5, iteration_converges=False, fixed_point_exists=True
    - abs_kappa=3, iteration_converges=False, fixed_point_exists=True
- **failures**: 0
- **the_fixed_point_exists_everywhere_but_one_point**: True
- **the_resonance_is_refused**: True
- **it_exists_even_where_the_iteration_diverges**: True
- **so_divergence_of_a_sum_is_not_absence_of_a_solution**: True

## T7_nonlinearity_is_what_would_break_it — PASS

- **rows**:
    - source=0.05, discriminant=0.88, n_real_solutions=2, roots=[0.05159737336276171, 1.615069293303905], worst_residual=2.77556e-17
    - source=0.2, discriminant=0.52, n_real_solutions=2, roots=[0.23240812075600176, 1.434258545910665], worst_residual=2.77556e-17
    - source=0.416, discriminant=0.0016, n_real_solutions=2, roots=[0.7999999999999995, 0.8666666666666671], worst_residual=1.11022e-16
    - source=0.5, discriminant=-0.2, n_real_solutions=0, roots=[], worst_residual=None
    - source=1, discriminant=-1.4, n_real_solutions=0, roots=[], worst_residual=None
- **kappa**: 0.6
- **threshold_source**: 0.416667
- **two_solutions_below_threshold**: True
- **none_above_it**: True
- **so_uniqueness_is_a_linearity_result**: True
- **every_reported_root_solves_the_loop**: True

## T8_the_ledger_closes_and_that_is_the_assumption — PASS

- **rows**:
    - kappa=0.2, emitted_by_source=1, returned_through_throat=0.0625, total_at_emitter=1.5625, residual=0, closes=True
    - kappa=0.5, emitted_by_source=1, returned_through_throat=1, total_at_emitter=4, residual=0, closes=True
    - kappa=-0.4, emitted_by_source=1, returned_through_throat=0.0816327, total_at_emitter=0.510204, residual=1.11022e-16, closes=True
    - kappa=(0.3+0.5j), emitted_by_source=1, returned_through_throat=0.459459, total_at_emitter=1.35135, residual=0, closes=True
- **every_ledger_closes**: True
- **worst_residual**: 1.11022e-16
- **but_that_is_the_assumption_not_the_result**: flux conservation through the throat is put in when the mouths are identified; the residual checks the arithmetic

## T9_the_delay_decides_the_story_and_nothing_conserved — PASS

- **rows**:
    - delay=-2, strike_receiver=4.28319, receiver_before_emitter=False, amplitude=1.81818, ledger_closes=True
    - delay=-12, strike_receiver=-5.71681, receiver_before_emitter=True, amplitude=1.81818, ledger_closes=True
    - delay=-100, strike_receiver=-93.7168, receiver_before_emitter=True, amplitude=1.81818, ledger_closes=True
    - delay=-10000, strike_receiver=-9993.72, receiver_before_emitter=True, amplitude=1.81818, ledger_closes=True
- **amplitude_spread**: 0
- **the_delay_changes_no_conserved_quantity**: True
- **but_it_decides_the_ordering**: True
- **so_the_picture_depends_on_it_and_the_physics_does_not**: True

## T10_two_local_events_one_conserved_wave — PASS

- **geodesic_emitter_to_future_mouth**: 3.14159
- **it_is_exactly_the_antipode**: True
- **sweep_distance_past_the_emitter**: 1.35
- **times**: emit=0, reach_future_mouth=3.14159, leave_past_mouth=-8.85841, strike_receiver=-5.71681, sweep_past_emitter=-7.50841, receiver_before_emitter=True
- **receiver_struck_before_emission**: True
- **flux_out_of_emitter**: 3.30579
- **flux_into_receiver**: 0.669421
- **emitter_alone_conserves**: False
- **receiver_alone_conserves**: False
- **the_pair_closes_to**: 1.11022e-16
- **the_pair_conserves**: True
- **what_is_apparent**: that the receiver was struck by an independent particle; the momentum transfer is real, the independence is not

## T11_assessment — PASS

- **n_passed**: 10
- **n_total**: 10

## verdict — ONE_WAVE_AND_LINEARITY_IS_WHY

ONE CONSERVED WAVE, AND LINEARITY IS WHY IT COSTS NOTHING. The scene's staging is geometry rather than choice: a geodesic sphere at distance chi on S3 has area 4 pi sin^2 chi, so an energy-conserving shell refocuses EXACTLY at the antipode. Checked against uniform sampling of S3 through the enclosed volume — the area's integral — and scored against the estimator's own binomial standard error rather than a fixed percentage: worst z = 2.02. That fact is used TWICE. The future mouth sits at the emitter's antipode, and the receiver sits at the past mouth's, which is the only place the arriving shell is genuinely COLLAPSING (dA/dchi = 4 pi sin 2chi < 0 all the way in) rather than merely arriving; against a receiver displaced to chi = 1.2 the same wave is still expanding when it lands. THE CONTENT IS THE SELF-CONSISTENCY. A closed timelike loop carrying a LINEAR wave has exactly one amplitude, A = A_src/(1 - kappa), matching brute iteration of the loop — which solves nothing and just keeps sending the wave around — to 7.1e-13, including at kappa = 0.99 where the amplification is x100. The fixed point exists and is unique for EVERY kappa != 1, the single resonance being the only obstruction in the whole complex plane; it exists even where |kappa| > 1 and the iteration diverges, because divergence of a summation method is not absence of a solution. Nothing is tuned for consistency and no paradox is available. THAT IS A STATEMENT ABOUT LINEAR EQUATIONS, NOT ABOUT TIME TRAVEL, and the probe shows the difference rather than asserting it: the same loop with a quadratic return has two solutions below a source threshold of 0.4167 and NONE above it. What each observer sees is two unrelated events, and neither conserves anything alone — the emitter loses flux (3.3058), the receiver gains it (0.6694) — while the pair closes to 1.1e-16. That balance is the only thing making them one object rather than two. THE PICTURE MEASURES RATHER THAN ILLUSTRATES: every drawn point sits at geodesic distance chi from its centre to 6.7e-15, and the shadow's screen extent is proportional to sin chi with one constant to 3.6e-16 — that is sqrt(A/4pi), so the apparent size in the figure IS the area law. Getting there required changing the projection, which is this round's real correction. A stereographic chart is unbounded at its own pole and a shell launched from a point sweeps ALL of S3, so whatever pole is chosen the shell crosses it once: the radius grows as 2/epsilon across four decades and never converges. The first renderer projected from q3 = +1, which is the emitter's own position — the emitter was a division by zero and never got drawn. WHAT IS PUT IN, PLAINLY. The wormhole is an identification map, not a solution. Flux conservation through the throat is an ASSUMPTION made when the mouths are identified, so the ledger residual checks the arithmetic and is not evidence for the model. PR #249 priced this throat — a minimal surface has sigma < 0 identically, so a connected throat needs exotic matter — and this round inherits that bill without paying it. And the delay, which decides the entire story, changes NO conserved quantity: the amplitude spread across delays from -2 to -1e4 is 0.0e+00, while whether the receiver is struck before the emitter fires flips inside that range.
