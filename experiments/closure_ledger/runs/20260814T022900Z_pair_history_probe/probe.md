# probe — pair_history

_2026-08-14T02:29:00.589031+00:00_

## T1_goal — PASS

- **question**: if two closed wormhole-wave histories must share their interaction event, is that event selected by the closure conditions, or is it still being inserted by hand?
- **scope**: a counting result on a kinematic skeleton. Fixed round S3; throats are identification maps with GIVEN mouths and delays; two of them, with PR #249's exotic-matter bill inherited and unpaid. No action principle, no field equations, no topology change, no dynamics, no rate, no worldline, and no claim of realisability. Conjugacy is a label, not a derivation of charge.

## T2_closure_is_a_geodesic_ellipsoid — PASS

- **rows**:
    - mouth_separation=1.69441, band=[1.6944147666652507, 4.588770540514336], sampled_min=1.69476, sampled_max=4.58871
    - mouth_separation=1.36208, band=[1.3620839518573977, 4.921101355322189], sampled_min=1.36209, sampled_max=4.92081
    - mouth_separation=1.90196, band=[1.9019634656870652, 4.381221841492521], sampled_min=1.90197, sampled_max=4.38115
    - mouth_separation=1.8633, band=[1.8632964678571309, 4.419888839322455], sampled_min=1.86365, sampled_max=4.41984
- **sampling_never_goes_below_the_band**: True
- **sampling_never_goes_above_the_band**: True
- **the_band_is_mouth_separation_to_two_pi_minus_it**: True
- **an_infeasible_delay_is_rejected**: True
- **so_closure_is_an_ellipsoid_condition**: the locus of points whose summed distance to the two mouths is |Δ| — a geodesic ellipsoid with the mouths as foci

## T3_the_event_is_selected_not_inserted — PASS

- **rows**:
    - n_events=2, ranks=[5, 5], worst_residual=2.63274e-11
    - n_events=2, ranks=[5, 5], worst_residual=3.75033e-11
    - n_events=0, ranks=[], worst_residual=None
    - n_events=0, ranks=[], worst_residual=None
    - n_events=2, ranks=[5, 5], worst_residual=4.02611e-12
    - n_events=2, ranks=[5, 5], worst_residual=9.76419e-12
    - n_events=2, ranks=[5, 5], worst_residual=9.23128e-12
    - n_events=0, ranks=[], worst_residual=None
    - n_events=0, ranks=[], worst_residual=None
    - n_events=2, ranks=[5, 5], worst_residual=2.05556e-11
    - n_events=0, ranks=[], worst_residual=None
    - n_events=0, ranks=[], worst_residual=None
- **configurations**: 12
- **configurations_with_a_selected_event**: 6
- **total_events_found**: 12
- **events_at_full_rank_5**: 12
- **every_event_is_nondegenerate**: True
- **the_counts_are_small_and_finite**: True
- **unknowns**: 5
- **equations**: 5
- **the_event_is_selected_not_inserted**: True
- **what_this_is_and_is_not**: a counting result on a kinematic skeleton: the closure conditions leave the interaction event isolated. No action principle, no field equations, no topology change, and no claim the configuration is realisable

## T4_removing_a_wave_removes_the_selection — PASS

- **rows**:
    - with_both_waves=0, ranks_both=[], with_one_wave=0, ranks_one=[], solution_spread_one_wave=0
    - with_both_waves=0, ranks_both=[], with_one_wave=0, ranks_one=[], solution_spread_one_wave=0
    - with_both_waves=2, ranks_both=[5, 5], with_one_wave=182, ranks_one=[4], solution_spread_one_wave=1.98883
    - with_both_waves=2, ranks_both=[5, 5], with_one_wave=147, ranks_one=[4], solution_spread_one_wave=1.72776
    - with_both_waves=2, ranks_both=[5, 5], with_one_wave=170, ranks_one=[4], solution_spread_one_wave=1.77061
    - with_both_waves=2, ranks_both=[5, 5], with_one_wave=159, ranks_one=[4], solution_spread_one_wave=1.71918
    - with_both_waves=0, ranks_both=[], with_one_wave=74, ranks_one=[4], solution_spread_one_wave=1.01293
    - with_both_waves=0, ranks_both=[], with_one_wave=0, ranks_one=[], solution_spread_one_wave=0
- **configurations**: 8
- **configurations_admitting_a_pair_history**: 4
- **two_waves_give_isolated_events**: True
- **one_wave_gives_a_one_parameter_family**: True
- **typical_family_size_with_one_wave**: 159
- **nullity_with_one_wave**: 1
- **the_solutions_do_not_vanish_they_stop_being_isolated**: True
- **dropping_a_constraint_can_even_create_solutions**: True
- **and_the_invariant_is_undefined_with_one_wave**: s = 2E₁E₂(1 − cos θ) needs two momenta; with a single front there is no second independent history to form an opening angle with
- **the_selection_requires_both_waves**: True

## T5_a_shared_throat_cannot_carry_the_pair — PASS

- **opposite_traversal_rows**:
    - delay_a=-3.23824, delay_b=3.23824, b_requires_negative_path_length=True, b_throat_is_feasible=False
    - delay_a=-1.58785, delay_b=1.58785, b_requires_negative_path_length=True, b_throat_is_feasible=False
    - delay_a=-2.2873, delay_b=2.2873, b_requires_negative_path_length=True, b_throat_is_feasible=False
    - delay_a=-1.96001, delay_b=1.96001, b_requires_negative_path_length=True, b_throat_is_feasible=False
    - delay_a=-4.60796, delay_b=4.60796, b_requires_negative_path_length=True, b_throat_is_feasible=False
    - delay_a=-5.26012, delay_b=5.26012, b_requires_negative_path_length=True, b_throat_is_feasible=False
- **same_traversal_rows**:
    - n_solutions=93, ranks=[4]
    - n_solutions=0, ranks=[]
    - n_solutions=106, ranks=[4]
    - n_solutions=167, ranks=[4]
    - n_solutions=116, ranks=[4]
    - n_solutions=105, ranks=[4]
- **opposite_traversal_is_infeasible**: True
- **configurations_admitting_a_solution**: 5
- **same_traversal_loses_a_rank**: True
- **same_traversal_gives_a_family_not_a_selection**: True
- **so_the_pair_needs_two_distinct_throats**: True
- **and_that_was_not_assumed**: the two-history picture does not require it; it follows from the rank of the system

## T6_the_delays_must_be_given_not_solved_for — PASS

- **events_with_delays_given**: 2
- **unknowns_with_delays_given**: 5
- **equations**: 5
- **unknowns_with_delays_free**: 7
- **nullity_with_delays_free**: 2
- **sampled_events_on_both_fronts**: 345
- **fraction_closable_by_choosing_a_delay**: 1
- **with_free_delays_almost_any_event_closes**: True
- **so_the_content_is_in_the_throat_being_given**: True
- **which_is_where_a_circular_version_of_this_would_hide**: solving for Δ after choosing the event would make every event a solution and the selection meaningless

## T7_the_threshold_is_a_separate_condition — PASS

- **configurations**: 14
- **configurations_with_a_selected_event**: 9
- **selected_events**: 18
- **rows**:
    - energy_over_mass=1, fraction_clearing_threshold=0
    - energy_over_mass=1.5, fraction_clearing_threshold=0.777778
    - energy_over_mass=2, fraction_clearing_threshold=0.944444
    - energy_over_mass=3, fraction_clearing_threshold=1
- **median_s_at_energy_equal_mass**: 2.4787
- **max_s_at_energy_equal_mass**: 3.97198
- **none_clear_threshold_at_energy_equal_mass**: True
- **most_clear_it_by_one_and_a_half_masses**: True
- **closure_selects_where_the_invariant_decides_whether**: True

## T8_the_conjugacy_is_carried_not_derived — PASS

- **orientations**:
    - 1
    - -1
- **net_orientation**: 0
- **the_labels_cancel**: True
- **a_same_sign_pair_is_refused**: True
- **but_nothing_here_derives_charge**: the kinematics is blind to the sign; the label is carried through the system and checked at the end, which is bookkeeping
- **and_the_throat_bill_is_still_inherited**: shells.junction priced a connected throat as necessarily exotic; this round adds two of them and pays for neither

## T9_assessment — PASS

- **n_passed**: 8
- **n_total**: 8

## verdict — THE_EVENT_IS_SELECTED_AND_LOSING_A_WAVE_LOSES_ISOLATION

THE INTERACTION EVENT IS SELECTED, NOT INSERTED — AND REMOVING A WAVE COSTS ISOLATION RATHER THAN EXISTENCE. Requiring two closed wormhole-wave histories to share one interaction event gives five equations in five unknowns: normalisation, the event lying on each of the two null fronts, and each history closing in coordinate time. Closure is a GEODESIC ELLIPSOID condition — the locus whose summed distance to the two mouths is |Delta| — feasible exactly on [d, 2pi - d], verified against 40,000 uniform samples of S3 per configuration rather than asserted. Solved BLIND from random starts, every event found sits at FULL JACOBIAN RANK 5 (12 of 12), so isolation is a property of the system and not of the solver stopping early. Existence is also RESTRICTIVE: only about half of random feasible configurations (6 of 12) admit a closed pair-history at all, so the closure conditions can forbid a configuration outright rather than merely locate an event inside it. THE FALSIFICATION LANDS DIFFERENTLY FROM THE EXPECTATION. Removing the second incoming wave does NOT delete the pair-history solution. The system becomes four equations in five unknowns, the Jacobian drops to RANK 4, and the solutions become a ONE-PARAMETER FAMILY — typically 159 distinct sampled points spanning the sphere. There is still a locus of events closing both histories; there is no longer a SELECTED one, and dropping the constraint can even CREATE solutions where two waves admitted none. So the Breit-Wheeler two-wave requirement appears here as loss of ISOLATION, which is a sharper and weaker statement than nonexistence, and it is the weaker one that is true. AND THE CONJUGATE PAIR NEEDS TWO DISTINCT THROATS, which the two-history picture did not assume — it falls out of the rank. A single shared throat fails in both senses: traversed oppositely, history B sees Delta_B = -Delta_A > 0 and its closure demands a sum of geodesic distances that is NEGATIVE, infeasible identically; traversed the same way, the two closure equations coincide, the rank drops to 4, and the event stops being selected. THE NON-CIRCULARITY CHECK IS THE ONE THAT MATTERS. If the delays were unknowns rather than given data the system would be five equations in seven unknowns with nullity 2, and every event on both fronts could be closed by choosing Delta afterwards — measured at 100% of 345 sampled events. The entire result therefore rests on the throat being GIVEN, and a version that solved for the delays would select nothing while looking identical from outside. FINALLY, CLOSURE SELECTS WHERE AND THE INVARIANT DECIDES WHETHER — kept strictly apart, because conflating amplitude with invariant is the error PR #252 unwound. NOT ONE selected event clears s >= 4m^2 at E = m (median 2.48, maximum 3.97, just under), while E = 1.5m clears 78% and E = 3m clears all of them. WHAT THIS IS NOT: a counting result on a kinematic skeleton. No action principle, no field equations, no topology change, no dynamics, no rate, no worldline, and no claim that such a configuration is realisable. The conjugacy Q_A + Q_B = 0 is a label carried and checked, never derived, and two throats now sit on the books with PR #249's exotic-matter bill unpaid for both.
