# probe — pair_history

_2026-08-14T03:33:38.371489+00:00_

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

## T3_the_results_are_scoped_to_the_principal_branch — PASS

- **principal_band_branch_counts**:
    - 1
    - 1
    - 1
    - 1
    - 1
    - 1
- **wide_delay_branch_counts**:
    - 4
    - 4
    - 4
    - 4
    - 4
    - 4
- **inside_the_band_only_the_principal_branch_is_feasible**: True
- **outside_it_more_branches_open**: True
- **locus_kinds_off_the_principal_branch**:
    - difference
    - sum
- **off_branch_loci_are_difference_type**: True
- **off_branch_rows**:
    - branches=[[0, 0, 1, 0], [1, 0, 0, 0]], n_events=2, ranks=[5]
    - branches=[[1, 0, 0, 0], [0, 0, 1, 0]], n_events=2, ranks=[5]
    - branches=[[0, 0, 1, 0], [0, 0, 1, 0]], n_events=0, ranks=[]
    - branches=[[1, 0, 0, 0], [0, 0, 1, 0]], n_events=2, ranks=[5]
    - branches=[[0, 0, 1, 0], [1, 0, 0, 0]], n_events=2, ranks=[5]
    - branches=[[0, 0, 1, 0], [1, 0, 0, 0]], n_events=0, ranks=[]
- **off_branch_rows_use_difference_type_branches**: True
- **discreteness_survives_on_a_fixed_off_branch**: True
- **so_the_other_results_are_principal_branch_scoped**: the prior draws |Δ| inside [D, 2π − D], where the principal branch is the only feasible one; branching changes the number of candidate events and the existence rate, not the local structure

## T4_the_events_are_discrete_and_locally_isolated — PASS

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
- **events_are_discrete_and_locally_isolated**: True
- **but_not_proved_exhaustive**: multi-start root-finding plus full rank shows each root found is locally isolated; it does not show all roots were found, nor that the event is unique
- **branch_scope**: principal (short-way, no winding)
- **what_this_is_and_is_not**: a counting result on a kinematic skeleton: on a fixed branch the closure conditions leave the allowed events discrete. No action principle, no field equations, no topology change, and no claim the configuration is realisable

## T5_removing_an_equation_is_a_dimensionality_control — PASS

- **rows**:
    - with_both_waves=0, ranks_both=[], with_one_wave=0, ranks_one=[], dropping_a_closure_instead=0, ranks_dropping_a_closure=[], solution_spread_one_wave=0
    - with_both_waves=0, ranks_both=[], with_one_wave=118, ranks_one=[4], dropping_a_closure_instead=0, ranks_dropping_a_closure=[], solution_spread_one_wave=0.754975
    - with_both_waves=2, ranks_both=[5, 5], with_one_wave=162, ranks_one=[4], dropping_a_closure_instead=187, ranks_dropping_a_closure=[4], solution_spread_one_wave=1.9745
    - with_both_waves=2, ranks_both=[5, 5], with_one_wave=56, ranks_one=[4], dropping_a_closure_instead=183, ranks_dropping_a_closure=[4], solution_spread_one_wave=0.46613
    - with_both_waves=2, ranks_both=[5, 5], with_one_wave=66, ranks_one=[4], dropping_a_closure_instead=179, ranks_dropping_a_closure=[4], solution_spread_one_wave=0.994533
    - with_both_waves=2, ranks_both=[5, 5], with_one_wave=133, ranks_one=[4], dropping_a_closure_instead=185, ranks_dropping_a_closure=[4], solution_spread_one_wave=1.85784
    - with_both_waves=0, ranks_both=[], with_one_wave=33, ranks_one=[4], dropping_a_closure_instead=212, ranks_dropping_a_closure=[4], solution_spread_one_wave=0.960875
    - with_both_waves=0, ranks_both=[], with_one_wave=0, ranks_one=[], dropping_a_closure_instead=211, ranks_dropping_a_closure=[4], solution_spread_one_wave=0
- **configurations**: 8
- **configurations_admitting_a_pair_history**: 4
- **two_waves_give_isolated_events**: True
- **one_wave_gives_a_one_parameter_family**: True
- **typical_family_size_with_one_wave**: 92
- **nullity_with_one_wave**: 1
- **the_solutions_do_not_vanish_they_stop_being_isolated**: True
- **deleting_a_closure_instead_drops_the_rank_the_same_way**: True
- **this_is_a_dimensionality_control_not_a_physics_result**: deleting ANY one scalar equation from a square nondegenerate system drops the rank by one and restores a continuous degree of freedom; nothing here singles out the wave constraint, and the two-photon content lives in the invariant s, measured separately
- **dropping_a_constraint_can_even_create_solutions**: True
- **and_the_invariant_is_undefined_with_one_wave**: s = 2E₁E₂(1 − cos θ) needs two momenta; with a single front there is no second independent history to form an opening angle with
- **the_square_system_behaves_nondegenerately**: True

## T6_a_shared_throat_cannot_carry_the_pair — PASS

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
- **branch_scan_rows**:
    - feasible_branches=3, distinct_branch_pairs=3
    - feasible_branches=3, distinct_branch_pairs=3
    - feasible_branches=3, distinct_branch_pairs=3
    - feasible_branches=3, distinct_branch_pairs=3
    - feasible_branches=3, distinct_branch_pairs=3
    - feasible_branches=3, distinct_branch_pairs=3
- **branch_pairs_that_restored_discreteness**: 0
- **no_branch_pair_rescues_a_shared_throat**: True
- **the_opposite_traversal_result_holds_on_every_branch**: leg lengths are non-negative on every branch, so a closure demanding a negative sum is infeasible regardless of winding
- **the_same_traversal_result_is_scoped**: shown in the minimal single-pass model and scanned to winding 1; a different gluing, or higher winding, is not excluded
- **so_in_this_model_the_pair_needs_two_distinct_throats**: True

## T7_the_delays_must_be_given_not_solved_for — PASS

- **events_with_delays_given**: 2
- **unknowns_with_delays_given**: 5
- **equations**: 5
- **unknowns_with_delays_free**: 7
- **measured_rank_of_the_five_by_seven_system**: 5
- **nullity_with_delays_free**: 2
- **the_nullity_is_measured_not_counted**: True
- **feasibility_checked_for_both_throats**: True
- **sampled_events_on_both_fronts**: 345
- **fraction_closable_by_choosing_a_delay**: 1
- **with_free_delays_almost_any_event_closes**: True
- **so_the_content_is_in_the_throat_being_given**: True
- **which_is_where_a_circular_version_of_this_would_hide**: solving for Δ after choosing the event would make every event a solution and the selection meaningless

## T8_the_threshold_is_a_separate_condition — PASS

- **analytic_bound**: s = 2E²(1 − cos θ) ≤ 4E², equality only at θ = π
- **zero_percent_at_E_equals_m_is_forced_not_measured**: True
- **fractions_are_prior_dependent_diagnostics**: conditioned on _random_system's arbitrary prior over mouths, delays and launch times; regression diagnostics, not predictions
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

## T9_the_conjugacy_is_carried_not_derived — PASS

- **orientations**:
    - 1
    - -1
- **net_orientation**: 0
- **the_labels_cancel**: True
- **a_same_sign_pair_is_refused**: True
- **but_nothing_here_derives_charge**: the kinematics is blind to the sign; the label is carried through the system and checked at the end, which is bookkeeping
- **and_the_throat_bill_is_still_inherited**: shells.junction priced a connected throat as necessarily exotic; this round adds two of them and pays for neither

## T10_assessment — PASS

- **n_passed**: 9
- **n_total**: 9

## verdict — DISCRETE_ON_A_FIXED_BRANCH_WITH_GIVEN_THROAT_DATA

WITH FIXED THROAT DATA AND ON A FIXED PROPAGATION BRANCH, TWO CLOSED HISTORIES SHARING AN EVENT LEAVE THE ALLOWED EVENTS DISCRETE. Requiring two closed wormhole-wave histories to share one interaction event gives five equations in five unknowns: normalisation, the event lying on each of the two null fronts, and each history closing in coordinate time. Solved BLIND from random starts, every root found sits at FULL JACOBIAN RANK 5 (12 of 12). Said precisely: multi-start root-finding plus full rank shows each root FOUND is locally isolated — not that all roots were found, and not that the event is unique. An earlier draft called this 'the event is selected'; that was an overstatement and is withdrawn. THE BRANCH SCOPE IS LOAD-BEARING AND IS STATED FIRST. d is the PRINCIPAL geodesic distance, and inside the principal delay band it is the only feasible branch — so every other measurement here is principal-branch BY CONSTRUCTION OF ITS PRIOR rather than by argument. Sampling the whole delay axis opens up to four branches, of both sum and DIFFERENCE type, and a mixed branch fixes the difference of the two distances: a hyperboloid, not an ellipsoid. What survives is that discreteness is a PER-BRANCH property — on any fixed branch the system is still 5x5 with isolated full-rank roots — so branching multiplies the candidate count and shifts the existence rate without touching the local structure. REMOVING A WAVE IS A DIMENSIONALITY CONTROL, NOT PHYSICS. Dropping wave B deletes one scalar equation from a square nondegenerate system, so rank 5 -> 4 and a one-parameter family is exactly what the implicit function theorem gives for ANY deleted equation. The probe therefore deletes a CLOSURE equation instead and gets the identical drop. This is evidence of nondegeneracy and NOT evidence that pair creation needs two photons; that content lives in the invariant s, which needs two independent momenta. What survives as interesting is only the direction of the surprise: the solutions do not vanish, they stop being isolated. IN THIS MODEL THE CONJUGATE PAIR NEEDS TWO DISTINCT THROATS. Traversed oppositely, history B's closure demands a NEGATIVE sum of leg lengths, infeasible on EVERY branch — the one conclusion here independent of branch scope. Traversed the same way on the same branch, the two closure equations coincide and the rank drops to 4; that half is scoped to the minimal single-pass model and SCANNED rather than argued, with every distinct branch pair at winding <= 1 either reducing to the identical constraint or jointly inconsistent, and no counterexample found. Not found is not impossible, and a different gluing is not excluded. THE NON-CIRCULARITY CHECK IS THE ONE THAT MATTERS: with the delays as unknowns the nullity is MEASURED on the actual 5x7 Jacobian as 2, and 100% of 345 sampled events on both fronts can then be closed by choosing Delta afterwards, with feasibility checked for BOTH throats. The entire result therefore rests on the throat being GIVEN. FINALLY, CLOSURE CONSTRAINS WHERE AND THE INVARIANT DECIDES WHETHER, and the numbers carry two warnings: that no event clears s >= 4m^2 at E = m is FORCED rather than measured, since s = 2E^2(1 - cos theta) <= 4E^2 with equality only at exactly head-on, a measure-zero set; and every fraction reported is conditioned on an arbitrary prior over mouths, delays and launch times, so they are regression diagnostics and NOT predictions. WHAT THIS IS NOT: a counting result on a kinematic skeleton. No action principle, no field equations, no topology change, no dynamics, no rate, no worldline, and no claim of realisability. Conjugacy is a label carried and checked, never derived, and two throats sit on the books with PR #249's exotic-matter bill unpaid for both.
