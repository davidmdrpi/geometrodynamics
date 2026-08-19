# probe — neck

_2026-08-19T06:13:09.934242+00:00_

## T1_goal — PASS

- **question**: PR #261 answered PR #260's gate with spheres in a FIXED ambient -- the balls were never removed, and only l = 0 coupled. Remove them: solve on S^3 minus two balls, glue the tube along the boundary spheres, and re-ask whether the negative mode survives.
- **scope**: a linear conformally coupled scalar on a FIXED round S^3 with two balls cut out. The tube has ONE transverse channel, so A is a coupling and the neck is a quantum-graph edge rather than a solved geometry. No backreaction.

## T2_the_exterior_map — PASS

- **ells**:
    - 0
    - 1
    - 2
    - 5
    - 10
- **rows**:
    - radius=0.02, sigma=0.05, by_ell=[0.252883, 0.578703, 0.847472, 1.61771, 2.88147], min_over_ell=0.252883
    - radius=0.02, sigma=1, by_ell=[0.256306, 0.578765, 0.847499, 1.61772, 2.88147], min_over_ell=0.256306
    - radius=0.02, sigma=5, by_ell=[0.27639, 0.580226, 0.848142, 1.61796, 2.88159], min_over_ell=0.27639
    - radius=0.02, sigma=30, by_ell=[0.402037, 0.622966, 0.870592, 1.62683, 2.886], min_over_ell=0.402037
    - radius=0.02, sigma=150, by_ell=[1.00514, 1.0982, 1.23921, 1.82586, 2.99222], min_over_ell=1.00514
    - radius=0.2, sigma=0.05, by_ell=[2.61661, 5.73803, 8.4007, 16.0484, 28.6015], min_over_ell=2.61661
    - radius=0.2, sigma=1, by_ell=[2.94555, 5.79121, 8.42623, 16.0583, 28.6064], min_over_ell=2.94555
    - radius=0.2, sigma=5, by_ell=[4.92673, 6.74024, 8.98769, 16.2935, 28.7244], min_over_ell=4.92673
    - radius=0.2, sigma=30, by_ell=[17.3264, 17.8636, 18.7274, 22.8638, 32.6779], min_over_ell=17.3264
    - radius=0.2, sigma=150, by_ell=[76.845, 76.9667, 77.169, 78.2526, 81.5726], min_over_ell=76.845
    - radius=0.5, sigma=0.05, by_ell=[6.38689, 13.6623, 20.0319, 38.4472, 68.7262], min_over_ell=6.38689
    - radius=0.5, sigma=1, by_ell=[8.20496, 14.2435, 20.3564, 38.5825, 68.7945], min_over_ell=8.20496
    - radius=0.5, sigma=5, by_ell=[19.7289, 22.4367, 26.3728, 41.652, 70.4137], min_over_ell=19.7289
    - radius=0.5, sigma=30, by_ell=[91.9381, 92.5297, 93.5077, 98.6285, 113.291], min_over_ell=91.9381
    - radius=0.5, sigma=150, by_ell=[438.542, 438.666, 438.873, 439.988, 443.481], min_over_ell=438.542
- **capacitance**:
    - radius=0.02, N0_at_zero=0.25287, four_pi_a=0.251327, ratio=1.00614
    - radius=0.01, N0_at_zero=0.126057, four_pi_a=0.125664, ratio=1.00313
    - radius=0.005, N0_at_zero=0.062931, four_pi_a=0.0628319, ratio=1.00158
- **worst_shooting_vs_closed_form**: 1.73198e-14
- **smallest_dtn_seen**: 0.252883
- **it_is_positive_everywhere**: True
- **it_increases_with_ell**: True
- **the_capacitance_normalization_holds**: True
- **what_it_replaces**: PR #261's f(a)G(a), which averaged over a sphere while keeping the ball in the ambient

## T3_the_quadratic_form_is_positive — PASS

- **trials**: 40
- **rows**:
    - quotient=49.1265, ambient_fraction=0.194923
    - quotient=0.359093, ambient_fraction=0.177375
    - quotient=8.64787, ambient_fraction=0.849048
    - quotient=1.57227, ambient_fraction=0.953339
    - quotient=0.778199, ambient_fraction=0.526272
    - quotient=18.4582, ambient_fraction=0.465587
    - quotient=5.83259, ambient_fraction=0.335603
    - quotient=0.444694, ambient_fraction=0.247168
- **smallest_quotient**: 0.359093
- **largest_quotient**: 49.1265
- **lowest_computed_mode**: 0.107516
- **interior_only_quotient**: 12.1847
- **poincare_floor**: 12.1847
- **every_quotient_is_positive**: True
- **no_trial_beats_the_lowest_mode**: True
- **the_interior_only_case_hits_the_poincare_floor**: True
- **the_theorem**: every term is non-negative and the degenerate case is floored by π²/L², so λ > 0 for all multipoles with no truncation and no sweep

## T4_the_negative_mode_does_not_survive — PASS

- **samples**: 1197
- **positives**: 0
- **worst_approach_to_zero**: -0.00163683
- **the_tube_side_is_always_negative**: True
- **the_exterior_side_is_always_positive**: True
- **there_is_no_growing_mode**: True
- **the_answer**: no, with the balls removed as well — and here the quadratic form settles it without any sweep

## T5_the_higher_multipoles_decouple — PASS

- **rows**:
    - ell=1, sigma=0.1, dtn=1.44613, monopole_dtn=0.637747, stiffer_than_monopole=True
    - ell=1, sigma=1, dtn=1.44708, monopole_dtn=0.658791, stiffer_than_monopole=True
    - ell=1, sigma=10, dtn=1.52614, monopole_dtn=0.941169, stiffer_than_monopole=True
    - ell=2, sigma=0.1, dtn=2.1177, monopole_dtn=0.637747, stiffer_than_monopole=True
    - ell=2, sigma=1, dtn=2.11811, monopole_dtn=0.658791, stiffer_than_monopole=True
    - ell=2, sigma=10, dtn=2.15823, monopole_dtn=0.941169, stiffer_than_monopole=True
    - ell=3, sigma=0.1, dtn=2.76635, monopole_dtn=0.637747, stiffer_than_monopole=True
    - ell=3, sigma=1, dtn=2.76661, monopole_dtn=0.658791, stiffer_than_monopole=True
    - ell=3, sigma=10, dtn=2.79303, monopole_dtn=0.941169, stiffer_than_monopole=True
    - ell=5, sigma=0.1, dtn=4.04256, monopole_dtn=0.637747, stiffer_than_monopole=True
    - ell=5, sigma=1, dtn=4.04271, monopole_dtn=0.658791, stiffer_than_monopole=True
    - ell=5, sigma=10, dtn=4.05839, monopole_dtn=0.941169, stiffer_than_monopole=True
- **every_channel_is_positive**: True
- **every_channel_is_stiffer_than_the_monopole**: True
- **smallest_ratio_to_the_monopole**: 1.62154
- **why_it_matters**: a one-channel tube drives only ℓ = 0, so the higher sectors are spectators; the (a/d)^ℓ screening bounds missed amplitude, not the stability answer

## T6_what_the_fixed_ambient_cost — PASS

- **rows**:
    - radius=0.02, lambda=0, fixed_ambient=3.95407, balls_removed=3.95459, ratio=0.999868, relative_error=0.000132488
    - radius=0.02, lambda=-4, fixed_ambient=3.82439, balls_removed=3.82684, ratio=0.999359, relative_error=0.000640646
    - radius=0.05, lambda=0, fixed_ambient=1.56752, balls_removed=1.56881, ratio=0.99918, relative_error=0.000820207
    - radius=0.05, lambda=-4, fixed_ambient=1.4437, balls_removed=1.44917, ratio=0.996225, relative_error=0.00377455
    - radius=0.15, lambda=0, fixed_ambient=0.508992, balls_removed=0.512659, ratio=0.992847, relative_error=0.00715264
    - radius=0.15, lambda=-4, fixed_ambient=0.401942, balls_removed=0.413552, ratio=0.971927, relative_error=0.0280735
    - radius=0.35, lambda=0, fixed_ambient=0.21049, balls_removed=0.218483, ratio=0.963416, relative_error=0.036584
    - radius=0.35, lambda=-4, fixed_ambient=0.127475, balls_removed=0.142798, ratio=0.892692, relative_error=0.107308
- **single_scattering**:
    - radius=0.02, separation=1.3, lambda=0, cross_over_self=0.0122437, over_a_over_d=0.795839, neglected_order=0.000149908
    - radius=0.02, separation=1.3, lambda=-4, cross_over_self=0.00160296, over_a_over_d=0.104192, neglected_order=2.56948e-06
    - radius=0.05, separation=1.3, lambda=0, cross_over_self=0.030885, over_a_over_d=0.803009, neglected_order=0.000953881
    - radius=0.05, separation=1.3, lambda=-4, cross_over_self=0.0042478, over_a_over_d=0.110443, neglected_order=1.80438e-05
    - radius=0.15, separation=1.3, lambda=0, cross_over_self=0.0951452, over_a_over_d=0.824592, neglected_order=0.00905261
    - radius=0.15, separation=1.3, lambda=-4, cross_over_self=0.0153884, over_a_over_d=0.133366, neglected_order=0.000236802
    - radius=0.35, separation=1.3, lambda=0, cross_over_self=0.230858, over_a_over_d=0.857474, neglected_order=0.0532956
    - radius=0.35, separation=1.3, lambda=-4, cross_over_self=0.0525254, over_a_over_d=0.195094, neglected_order=0.00275891
    - radius=0.35, separation=0.8, lambda=0, cross_over_self=0.394282, over_a_over_d=0.901215, neglected_order=0.155458
    - radius=0.35, separation=0.8, lambda=-4, cross_over_self=0.191886, over_a_over_d=0.438598, neglected_order=0.0368204
- **worst_error_at_the_working_radius**: 0.00377455
- **worst_error_overall**: 0.107308
- **the_error_grows_with_the_radius**: True
- **pr_261_was_right_where_it_looked**: True
- **cross_over_self_scales_like_a_over_d**: True
- **neglected_scattering_at_the_working_point**: 0.000953881
- **neglected_scattering_worst**: 0.155458
- **the_neglected_series_cannot_reach_the_sign**: True

## T7_the_soft_mode_is_forced — PASS

- **rows**:
    - length=0.3, f_sym_at_zero=5.30516e+08, f_sym_at_the_gap=-1.01323e+08, f_anti_at_zero=-1.5323, symmetric_roots=[0.316112], antisymmetric_roots=[], spurious_pole_crossings=0, tube_harmonics_in_the_gap=0, states=1
    - length=0.9, f_sym_at_zero=1.76839e+08, f_sym_at_the_gap=-1.01323e+08, f_anti_at_zero=-1.55617, symmetric_roots=[0.107516], antisymmetric_roots=[], spurious_pole_crossings=0, tube_harmonics_in_the_gap=0, states=1
    - length=2, f_sym_at_zero=7.95775e+07, f_sym_at_the_gap=-1.01323e+08, f_anti_at_zero=-1.59994, symmetric_roots=[0.0482144], antisymmetric_roots=[], spurious_pole_crossings=0, tube_harmonics_in_the_gap=0, states=1
    - length=3, f_sym_at_zero=5.30516e+07, f_sym_at_the_gap=-1.01323e+08, f_anti_at_zero=-1.63973, symmetric_roots=[0.0319287], antisymmetric_roots=[], spurious_pole_crossings=0, tube_harmonics_in_the_gap=0, states=1
    - length=4, f_sym_at_zero=3.97887e+07, f_sym_at_the_gap=-1.01323e+08, f_anti_at_zero=-1.67951, symmetric_roots=[0.0237708], antisymmetric_roots=[0.667445], spurious_pole_crossings=1, tube_harmonics_in_the_gap=1, states=2
    - length=6, f_sym_at_zero=2.65258e+07, f_sym_at_the_gap=-1.01323e+08, f_anti_at_zero=-1.75909, symmetric_roots=[0.0156063], antisymmetric_roots=[0.307777], spurious_pole_crossings=1, tube_harmonics_in_the_gap=1, states=2
    - length=8, f_sym_at_zero=1.98944e+07, f_sym_at_the_gap=-1.01323e+08, f_anti_at_zero=-1.83867, symmetric_roots=[0.0115248, 0.638527], antisymmetric_roots=[0.179227], spurious_pole_crossings=2, tube_harmonics_in_the_gap=2, states=3
- **gap_edge**:
    - lmbda=0.9, N0=0.486627, N0_over_one_minus_lambda=4.86627, closed_form_slope=19.7387, relative=0.753465
    - lmbda=0.99, N0=0.1511, N0_over_one_minus_lambda=15.11, closed_form_slope=19.7387, relative=0.234499
    - lmbda=0.999, N0=0.0191519, N0_over_one_minus_lambda=19.1519, closed_form_slope=19.7387, relative=0.0297296
    - lmbda=0.9999, N0=0.00196784, N0_over_one_minus_lambda=19.6784, closed_form_slope=19.7387, relative=0.00305477
    - lmbda=0.99999, N0=0.000197326, N0_over_one_minus_lambda=19.7326, closed_form_slope=19.7387, relative=0.00030632
    - lmbda=0.999999, N0=1.97381e-05, N0_over_one_minus_lambda=19.7381, closed_form_slope=19.7387, relative=3.06405e-05
- **the_symmetric_channel_starts_positive**: True
- **and_ends_negative**: True
- **so_a_symmetric_state_is_forced**: True
- **the_antisymmetric_channel_starts_negative**: True
- **gap_edge_slope_error**: 3.06405e-05
- **gap_edge_error_ratio_per_decade**: 0.100028
- **the_exterior_stiffness_vanishes_linearly_at_the_gap_edge**: True
- **short_tubes_hold_exactly_one_state**: True
- **long_tubes_hold_more**: True
- **every_extra_state_follows_a_tube_harmonic**: True
- **pole_crossings_that_would_have_been_miscounted**: 4
- **what_pr_261_got_wrong**: 'exactly one state below the gap' is a statement about L < π, not a structural one; it holds for every length that round used and fails above L = π

## T8_the_soft_mode_survives_the_removal — PASS

- **rows**:
    - radius=0.2, lambda_0_neck=0.349409, lambda_0_fixed_ambient=0.351588, closed_form=0.444444, shift_from_removal=0.00619955, ratio_to_closed_form=0.78617
    - radius=0.1, lambda_0_neck=0.203785, lambda_0_fixed_ambient=0.204256, closed_form=0.222222, shift_from_removal=0.00230597, ratio_to_closed_form=0.917035
    - radius=0.05, lambda_0_neck=0.107516, lambda_0_fixed_ambient=0.107591, closed_form=0.111111, shift_from_removal=0.000699234, ratio_to_closed_form=0.967644
    - radius=0.02, lambda_0_neck=0.0439781, lambda_0_fixed_ambient=0.0439836, closed_form=0.0444444, shift_from_removal=0.000124605, ratio_to_closed_form=0.989507
    - radius=0.01, lambda_0_neck=0.022115, lambda_0_fixed_ambient=0.0221157, closed_form=0.0222222, shift_from_removal=3.22376e-05, ratio_to_closed_form=0.995174
    - radius=0.005, lambda_0_neck=0.0110855, lambda_0_fixed_ambient=0.0110856, closed_form=0.0111111, shift_from_removal=8.19612e-06, ratio_to_closed_form=0.997693
- **every_one_is_positive_and_below_the_gap**: True
- **worst_shift_from_removal**: 0.00619955
- **shift_at_the_working_radius**: 0.000699234
- **worst_closed_form_error**: 0.010493
- **the_capacitance_formula_survives**: True
- **removing_the_balls_barely_moves_it**: True

## T9_assessment — PASS

- **n_passed**: 8
- **n_total**: 8

## verdict — THE_NEGATIVE_MODE_DOES_NOT_SURVIVE_THE_NECK

*** THE ANSWER IS STILL NO, AND IT IS NOW A THEOREM. *** PR #261 answered PR #260's gate with spheres in a FIXED ambient: the balls were never removed, and only l = 0 coupled. That was its own stated limitation and the strongest remaining doubt. This round removes them -- the ambient is S^3 minus two balls, the tube is glued along the boundary spheres -- and the answer does not move. WITH THE BALLS REMOVED THERE IS NO SUBTRACTION ANYWHERE, so the energy E = Int_Omega (|grad phi|^2 + phi^2) + A Int (|u'|^2 + m^2 |u|^2) is a sum of manifestly non-negative terms, and E = 0 forces phi = 0, whence matching gives u(0) = u(L) = 0 and Poincare gives A Int |u'|^2 >= (pi/L)^2 A Int |u|^2. Hence lam > 0 for EVERY configuration, ALL MULTIPOLES, NO TRUNCATION, NO SWEEP. That is a change of footing rather than a refinement: PR #261 had to establish a SIGN on a reduced 2x2, this round has POSITIVITY OF THE FORM ITSELF, and the renormalized g(lam) < 0 that PR #260's mode fed on has nowhere to enter because nothing is renormalized. The object the theorem is about is checked rather than asserted: the exterior DtN is computed by shooting the radial equation from the far pole and agrees with an independent closed form to 1.7e-14; it is positive in every channel and INCREASES with l, so the monopole is the softest and the higher channels cannot be the first to go soft; and N_0 -> 4 pi a, the capacitance of a sphere, which fixes the normalization as physical. Explicit trial configurations give Rayleigh quotients from 0.3591 up, all positive and all above the computed lowest mode 0.107516, and the degenerate purely-interior case lands on the Poincare floor pi^2/L^2 exactly. The sweep agrees with the theorem: 1197 samples give 0 positives, worst approach -1.6e-03, with the tube side negative and the exterior side positive at every one. *** AND PR #261's MONOPOLE TRUNCATION WAS NEVER A STABILITY LIMITATION. *** A one-channel tube presents a single number at each mouth, so it drives only l = 0; the higher sectors are the exterior's own modes with positive DtN, at least 1.62x stiffer than the monopole, and can neither be driven nor go unstable. PR #250's (a/d)^l screening bounds missed AMPLITUDE, not the answer. WHAT THE FIXED AMBIENT COST, PRICED: PR #261's f(a)G(a) against this round's 1/N_0 agrees to 3.8e-03 at the working radius and departs to 10.7% at a = 0.35, lam = -4, growing with the fraction of the sphere wrongly left in -- so that round was right where it looked, and now by a measured margin rather than a caveat. THE ONE APPROXIMATION LEFT IN THE REDUCED 2x2 IS ALSO PRICED: its cross term f(a)^2 G(d) is SINGLE-SCATTERING, and the neglected series expands in cross/self = 0.8 (a/d), which is 9.5e-04 at the working point and at worst 0.16 anywhere sampled -- too small to flip the sign the reduced model is asked to decide, and irrelevant to the theorem, which does not go through it. THE SOFT MODE IS FORCED, NOT FOUND: the same style of argument that kills the growing mode produces this one from the two ends of the gap, F_sym -> +inf as lam -> 0+ and -> -inf as lam -> 1-, the latter because the exterior stiffness VANISHES at the free ESU threshold, N_0 -> 2 pi (pi - a + sin a cos a)(1 - lam), reproduced to 3.1e-05 with first-order convergence (error ratio 0.1000 per decade). *** ONE CORRECTION TO PR #261: *** its 'exactly one state below the gap' is a statement about L < pi, not a structural one. The channel functions have POLES at lam = (2 pi j/L)^2 and (pi(2j-1)/L)^2; above L = pi these enter the gap and each brings a genuine extra state just above it -- three states at L = 8. A pole is a sign change with NO ZERO, so crossing-counting alone reports states that are not there (4 of them across the lengths swept here); classifying by the residual separates roots from poles by fifteen orders of magnitude. None of it touches the stability answer -- every one of those states is positive, as the form requires. AND THE SOFT MODE SURVIVES THE REMOVAL: lam_0 moves by 7.0e-04 at the working radius and still tends to 8 pi a/(A L), matched to 1.0%. WHAT IS STILL PUT IN: the tube has ONE transverse channel, so A is a coupling rather than a solved neck cross-section; the ambient metric is FIXED; and there is NO BACKREACTION. PR #260's gate is now closed from both sides -- resolved mouths and removed balls -- so stationary action and backreaction can proceed on a background whose positivity is proved rather than sampled.
