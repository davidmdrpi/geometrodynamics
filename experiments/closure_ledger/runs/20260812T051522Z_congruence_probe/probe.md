# congruence probe

_generated 2026-08-12T05:15:22.165410+00:00_

## T1_goal — PASS
- **claim**: visualise geodesic/tidal focusing all the way to the pinch threshold, and let the equations say what the threshold produces
- **object**: a congruence with a deforming cross-section, not a point
- **quantity**: J = det F, the cross-sectional area of the bundle
- **the_three_outcomes**: ['caustic -> passage', 'caustic -> singular termination', 'caustic -> finite-radius reconnection']

## T2_raychaudhuri_is_exact — PASS
- **gain**: 25
- **worst_raychaudhuri_residual**: 6.66134e-15
- **worst_ricci_term**: 0
- **raychaudhuri_is_exact**: True
- **the_ricci_term_vanishes**: True
- **so_the_focusing_is_all_shear**: True

## T3_the_operator_form_matters — PASS
- **rows**:
    - gain=25, n_radial=600, operator_min_J=-13.0553, seeded_difference_min_J=0.627331
    - gain=25, n_radial=1200, operator_min_J=-12.9055, seeded_difference_min_J=0.627937
    - gain=25, n_radial=2400, operator_min_J=-12.8356, seeded_difference_min_J=0.62807
- **operator_spread_over_the_ladder**: 0.0168212
- **seeded_spread_over_the_ladder**: 0.00117763
- **both_forms_converge**: True
- **but_they_disagree**: True
- **refinement_cannot_reveal_it**: True

## T4_the_front_is_causal — PASS
- **rows**:
    - n_radial=800, gaussian_arrival=1.99923, compact_arrival=2.64758, gaussian_early_by=0.770669, compact_early_by=0.122323
    - n_radial=1600, gaussian_arrival=2.0008, compact_arrival=2.69765, gaussian_early_by=0.769098, compact_early_by=0.0722537
    - n_radial=3200, gaussian_arrival=2.0011, compact_arrival=2.72779, gaussian_early_by=0.768804, compact_early_by=0.042114
- **probe_distance**: 2.9499
- **light_arrival**: 2.9499
- **causal_bound**: 2.7699
- **gaussian_spread_over_the_ladder**: 0.00186532
- **the_gaussian_arrival_is_grid_converged**: True
- **the_gaussian_beats_the_bound**: True
- **the_compact_launch_respects_it**: True
- **compact_earliness_shrinks_with_refinement**: True

## T5_the_three_cases_separate — PASS
- **rows**:
    - gain=1, peak_strain=0.0010328, min_J_near_the_source=0.975901, min_J_on_the_converging_ring=0.999516, source_caustic=False, ring_caustic=False, case=ordinary focus
    - gain=4, peak_strain=0.00413121, min_J_near_the_source=0.614529, min_J_on_the_converging_ring=0.992272, source_caustic=False, ring_caustic=False, case=ordinary focus
    - gain=20, peak_strain=0.020656, min_J_near_the_source=-8.56095, min_J_on_the_converging_ring=0.816625, source_caustic=True, ring_caustic=False, case=caustic
    - gain=80, peak_strain=0.0826242, min_J_near_the_source=-133.155, min_J_on_the_converging_ring=-0.000849665, source_caustic=True, ring_caustic=True, case=caustic
    - gain=150, peak_strain=0.15492, min_J_near_the_source=-307.438, min_J_on_the_converging_ring=-0.00276657, source_caustic=True, ring_caustic=True, case=caustic
- **n_ordinary_focus**: 2
- **n_caustic**: 3
- **both_cases_appear**: True
- **a_weak_wave_only_dips**: True
- **the_source_ring_closes_first**: True
- **case_three_never_appears**: True

## T6_two_rings_two_thresholds — PASS
- **window_in_units_of_pi**: 1.2
- **unit_gain_peak_strain**: 0.00382739
- **source_ring_threshold_gain**: 6.6635
- **source_ring_threshold_strain**: 0.0255038
- **converging_ring_threshold_gain**: 64.5286
- **converging_ring_threshold_strain**: 0.246976
- **ratio**: 9.68389
- **the_source_ring_is_far_easier_to_close**: True
- **closing_the_ring_needs_an_enormous_strain**: True

## T7_the_neck_is_a_ring — PASS
- **rows**:
    - pulse_width=0.12, gain=116.975, min_J=-0.000358593, at_time_over_pi=1.21375, neck_ring_radius=0.0566838, radius_over_width=0.472365, cells_from_the_pole=13
    - pulse_width=0.18, gain=52.3561, min_J=0.00103571, at_time_over_pi=1.3, neck_ring_radius=0.0784853, radius_over_width=0.43603, cells_from_the_pole=18
    - pulse_width=0.24, gain=29.5595, min_J=0.0686602, at_time_over_pi=1.3, neck_ring_radius=0.104647, radius_over_width=0.43603, cells_from_the_pole=24
    - pulse_width=0.3, gain=18.9958, min_J=0.158585, at_time_over_pi=1.3, neck_ring_radius=0.130809, radius_over_width=0.43603, cells_from_the_pole=30
    - pulse_width=0.4, gain=10.7763, min_J=0.284989, at_time_over_pi=1.3, neck_ring_radius=0.174412, radius_over_width=0.43603, cells_from_the_pole=40
- **strain**: 0.2
- **mean_radius_over_width**: 0.443297
- **relative_spread**: 0.0819672
- **one_cell_in_the_same_units**: 0.0819672
- **spread_is_at_the_resolution_floor**: True
- **width_range**: 3.33333
- **grid_ladder**: [(721, 0.4357271364202214), (1441, 0.43602951472446627), (2881, 0.4482969963470332)]
- **grid_drift**: 0.0273646
- **one_cell_at_the_finest_grid**: 0.0273319
- **the_ratio_converges**: True
- **the_radius_scales_with_the_width**: True
- **the_neck_is_resolved_off_the_pole**: True
- **the_neck_is_never_at_the_antipode**: True
- **because_h_vanishes_at_the_pole**: True

## T8_the_caustic_is_a_passage — PASS
- **gain**: 150
- **peak_strain**: 0.1246
- **crossing_slope**: -17.8768
- **crossing_slope_at_half_the_timestep**: -17.8356
- **crossing_slope_drift**: 0.00231422
- **source_ring**: region=near, tracked_distance=0.0915027, n_zero_crossings=1, crossing_times=[0.11938052083641301], crossing_slopes=[-17.87683366101361], relative_slopes=[0.19000265125585908], depth_past_zero=-471.486, samples_below_zero=24992, final_J=-471.486, finite=True, peak_strain=0.1246, invariant_drift=2.45772e-14
- **source_ring_at_half_the_timestep**: region=near, tracked_distance=0.0915027, n_zero_crossings=1, crossing_times=[0.1195768703772584], crossing_slopes=[-17.835558333282346], relative_slopes=[0.19051018446322016], depth_past_zero=-469.155, samples_below_zero=49983, final_J=-469.155, finite=True, peak_strain=0.124297, invariant_drift=3.80777e-13
- **converging_ring**: region=far, tracked_distance=3.06316, n_zero_crossings=2, crossing_times=[3.2311280442181713, 3.234858685494311], crossing_slopes=[-0.13354492862767905, 0.12770128603769995], relative_slopes=[0.35538789172742646, 0.3398368719964471], depth_past_zero=-0.000121705, samples_below_zero=19, final_J=109.037, finite=True, peak_strain=0.1246, invariant_drift=2.45772e-14
- **source_depth_drift_under_refinement**: 0.00496879
- **converging_ring_depth_drift**: 3.06438
- **crossings_are_transversal**: True
- **the_source_excursion_is_resolved**: True
- **the_converging_ring_only_grazes**: True
- **sign_flipped**: True
- **everything_stayed_finite**: True
- **evolution_continued**: True
- **solver_invariant_drift**: 2.45772e-14
- **the_caustic_is_a_passage**: True
- **not_a_termination**: True
- **no_finite_radius_bounce**: True
- **reconnection_was_never_available**: each material point's F is driven only by the external h and never by its neighbours, so the congruence cannot act back on anything; a finite-radius reconnection needs exactly that feedback, so this program could not have produced one, and 'we did not see it' would have been the wrong thing to say

## T9_case_three_is_out_of_scope — PASS
- **gain**: 100
- **peak_strain**: 0.091417
- **background_curvature**: 1
- **background_curvature_range**: 0
- **worst_invariant_drift**: 1.23057e-13
- **the_field_stayed_finite**: True
- **evolution_terminated**: False
- **the_geometry_never_moves**: True
- **case_three_is_out_of_scope**: True
- **reason**: the background is a fixed round S² with curvature 1 at every time; there is no metric evolution, no Einstein equation and no backreaction, so nothing here can diverge or stop, and a caustic is the strongest thing available

## T10_every_threshold_carries_its_window — PASS
- **rows**:
    - window_in_units_of_pi=1.2, threshold_gain=66.4844, threshold_strain=0.255239
    - window_in_units_of_pi=2, threshold_gain=30.0859, threshold_strain=0.115502
    - window_in_units_of_pi=3, threshold_gain=21.3085, threshold_strain=0.0818049
    - window_in_units_of_pi=4, threshold_gain=16.325, threshold_strain=0.062673
- **threshold_falls_with_the_window**: True
- **ratio_first_to_last**: 4.07255
- **so_a_threshold_without_a_window_is_meaningless**: True
- **the_accumulation_is_a_hill_equation**: True

## T11_assessment — PASS
- **n_passed**: 10
- **n_total**: 10

## verdict — THE_CAUSTIC_IS_A_PASSAGE

THE CAUSTIC IS A PASSAGE. Focusing was carried to the pinch threshold without assuming what the threshold would produce, and the equations chose the first of the three offered outcomes. J crosses zero transversally — slope -17.877, converged to -17.836 under a halved timestep, where a tangency would have driven it to zero — plunges to -471, and the evolution continues with the solver's invariant unmoved at 2.5e-14. Neither of the other two outcomes was available, and for different reasons worth keeping apart. Singular termination needs the geometry to fail, and the background here is a fixed round S² with curvature 1 at every time. Finite-radius reconnection needs the congruence to act back on something, and each point's F is driven only by the external h and never by its neighbours. So this is not 'we looked and did not find them'; it is 'this program could not have produced them', which is a different and more useful statement. What the program CAN say is where the neck forms and what it costs: on a ring around the antipode, never at it, because spin weight forces h to vanish at the pole, at a radius of 0.44 w — the same ratio across a 3.3x range of pulse width, to within one grid cell; and at peak strain 0.247 for the converging ring against 0.026 for the source ring, a factor of ten apart. Even then the antipodal crossing only grazes zero and the depth of its excursion does not converge. A spin-2 focus has no centre, and closing its neck costs a strain nothing physical would reach.
