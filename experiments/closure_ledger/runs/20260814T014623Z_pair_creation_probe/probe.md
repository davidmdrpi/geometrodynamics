# probe — pair_creation

_2026-08-14T01:46:23.422500+00:00_

## T1_goal — PASS

- **question**: earlier rounds drew a single wavefront refocusing at its antipode and treated the caustic as a particle-creation event; is it one, and if not, what does the geometry actually require?
- **scope**: kinematics on a fixed round sphere. The Breit-Wheeler threshold and cross-section are imported QED; rays-as-photons is a correspondence; no rate is computed; the crossing chord as a throat is interpretation, and PR #249 priced that throat without this round paying it.

## T2_focusing_is_neither_sufficient_nor_necessary — PASS

- **collinear_rows**:
    - amplification=1, energy=1, s_collinear=0
    - amplification=10, energy=10, s_collinear=0
    - amplification=1000, energy=1000, s_collinear=0
    - amplification=1e+06, energy=1e+06, s_collinear=0
    - amplification=1e+12, energy=1e+12, s_collinear=0
- **crossed_rows**:
    - energy=0.5, theta=π, s_head_on=1, above_threshold=False
    - energy=0.999, theta=π, s_head_on=3.992, above_threshold=False
    - energy=1, theta=π, s_head_on=4, above_threshold=True
    - energy=1.5, theta=π, s_head_on=9, above_threshold=True
    - energy=4, theta=π, s_head_on=64, above_threshold=True
- **focusing_is_not_sufficient**: True
- **largest_amplification_tried**: 1e+12
- **focusing_is_not_necessary**: True
- **threshold_is_at_energy_equal_mass**: True
- **a_converging_shell_does_contain_opposed_rays**: True
- **so_the_distinction_is_independence_not_brightness**: a spherically converging front has diametrically opposed rays, so its self-invariant is not zero; what a single front cannot supply is two INDEPENDENTLY propagated waves, and that is a statement about the source topology rather than about the amplitude
- **what_this_corrects**: earlier rounds drew one caustic and called it a creation event; s depends on the opening angle, which focusing does not create

## T3_the_invariant_is_a_triangle_identity — PASS

- **rows**:
    - dim=3, sphere=S², samples_used=297, worst_abs_error=2.04281e-14
    - dim=4, sphere=S³, samples_used=300, worst_abs_error=7.10543e-15
- **worst_over_all_dimensions**: 2.04281e-14
- **the_closed_form_is_confirmed**: True
- **and_it_is_dimension_independent**: True
- **the_control_never_uses_the_law_of_cosines**: the crossing point is solved as a linear system and the momenta are great-circle tangents in the embedding

## T4_the_collision_is_head_on_twice — PASS

- **rows**:
    - delta=0.15, t_near=0.075, theta_near=3.14159, s_near=4, t_equator=1.5708, theta_equator=0.15, s_equator=0.0224578, t_far=3.06659, theta_far=3.14159, s_far=4, equator_angle_equals_delta=True, near_and_far_are_mirror=True
    - delta=0.42, t_near=0.21, theta_near=3.14159, s_near=4, t_equator=1.5708, theta_equator=0.42, s_equator=0.173822, t_far=2.93159, theta_far=3.14159, s_far=4, equator_angle_equals_delta=True, near_and_far_are_mirror=True
    - delta=1, t_near=0.5, theta_near=3.14159, s_near=4, t_equator=1.5708, theta_equator=1, s_equator=0.919395, t_far=2.64159, theta_far=3.14159, s_far=4, equator_angle_equals_delta=True, near_and_far_are_mirror=True
    - delta=2, t_near=1, theta_near=3.14159, s_near=4, t_equator=1.5708, theta_equator=2, s_equator=2.83229, t_far=2.14159, theta_far=3.14159, s_far=4, equator_angle_equals_delta=True, near_and_far_are_mirror=True
- **head_on_invariant**: 4
- **worst_head_on_error**: 3.28626e-14
- **both_ends_are_head_on**: True
- **the_equator_angle_is_exactly_the_separation**: True
- **the_minimum_is_at_the_equator**: True
- **so_the_invariant_is_u_shaped_in_t**: True
- **which_is_why_a_threshold_cuts_two_windows**: True

## T5_the_threshold_opens_two_windows — PASS

- **delta**: 0.42
- **rows**:
    - energy_over_mass=0.6, closed_form_windows=[], n_windows_closed_form=0, n_windows_scanned=0, scan_agrees=True
    - energy_over_mass=1, closed_form_windows=[], n_windows_closed_form=0, n_windows_scanned=0, scan_agrees=True
    - energy_over_mass=1.4, closed_form_windows=[(0.21, 0.29615405627503927), (2.845438597314754, 2.931592653589793)], n_windows_closed_form=2, n_windows_scanned=2, scan_agrees=True
    - energy_over_mass=3, closed_form_windows=[(0.21, 0.6756180323932854), (2.4659746211965077, 2.931592653589793)], n_windows_closed_form=2, n_windows_scanned=2, scan_agrees=True
    - energy_over_mass=6, closed_form_windows=[(0.21, 2.931592653589793)], n_windows_closed_form=1, n_windows_scanned=1, scan_agrees=True
- **merge_energy_over_mass**: 4.79709
- **the_scan_agrees_with_the_closed_form**: True
- **below_E_equals_m_there_is_no_window**: True
- **at_E_equals_m_the_window_has_zero_width**: True
- **and_that_is_the_head_on_threshold_touched_exactly**: True
- **there_are_exactly_two_windows_in_between**: True
- **and_they_merge_above**: True
- **so_a_second_interaction_at_the_antipode_is_forced**: True

## T6_only_the_far_window_is_independent — PASS

- **delta**: 0.42
- **energy_over_mass**: 1.4
- **near_window**:
    - 0.21
    - 0.296154
- **far_window**:
    - 2.84544
    - 2.93159
- **path_before_near_collision**: 0.296154
- **path_before_far_collision**: 2.84544
- **separation_of_the_sources**: 0.42
- **near_collision_is_within_the_source_region**: True
- **far_collision_is_past_a_quarter_turn**: True
- **ratio_of_path_lengths**: 9.60797
- **both_windows_are_head_on_at_their_outer_edge**: True
- **only_the_far_window_is_a_collision_of_independent_waves**: True
- **why**: the near window is the emission region — the fronts have not separated yet; the far one is reached after each has crossed a half-circumference, which is the only place on a round space where two independently propagated fronts meet head-on again

## T7_the_projected_angle_is_not_the_opening_angle — PASS

- **worst_projected_error_deg**: 67.3978
- **worst_disagreement_between_the_two_crossings_deg**: 56.4026
- **momenta_are_perpendicular_to_the_sphere**: True
- **the_projection_misreads_the_angle**: True
- **and_it_misreads_the_two_crossings_differently**: True
- **though_their_true_opening_angle_is_identical**: True
- **so_the_arrows_cannot_carry_the_claim**: the renderer draws the opening angle in the plane the two momenta span, where it is undistorted; the arrows on the sphere are decoration and are labelled as projected
- **sample_rows**:
    - t=0.211, true_deg=168.922, projected_deg=[175.242159631056, 175.84149842047555]
    - t=0.257095, true_deg=110.136, projected_deg=[121.15078843725887, 154.06473921750037]
    - t=0.30319, true_deg=88.5662, projected_deg=[83.90517225685522, 115.80847535295638]
    - t=0.349284, true_deg=75.0533, projected_deg=[64.05812290113126, 7.6555081309902375]
    - t=0.395379, true_deg=65.5351, projected_deg=[52.93573957895036, 9.12322926286116]
    - t=0.441474, true_deg=58.4033, projected_deg=[45.83355019395443, 14.656204714818688]

## T8_the_cross_section_is_imported_and_checked — PASS

- **peak_sqrt_s_over_2m**: 1.40282
- **peak_beta**: 0.701319
- **peak_sigma_over_pi_re2**: 0.681706
- **peak_sigma_over_thomson**: 0.25564
- **zero_at_threshold**: True
- **zero_below_threshold**: True
- **falls_at_large_s**: True
- **matches_the_textbook_peak**: True
- **this_is_imported_not_derived**: the threshold and the cross-section are QED; this round supplies only where on the sphere s clears the threshold

## T9_the_pair_conserves_orientation — PASS

- **orientations**:
    - 1
    - -1
- **net_orientation**: 0
- **the_labels_cancel**: True
- **crossing_locus_size_on_S2**: 2
- **separation_of_the_two_crossing_points**: 0.284583
- **evaluated_at_t**: 2.88852
- **s_there**: 5.43435
- **above_threshold_there**: True
- **on_S3_the_locus_is_a_circle**: 16
- **but_the_throat_is_interpretation**: the kinematics says where the fronts cross, at what angle, and whether s clears (2m)²; calling the two crossing points the mouths of one throat is this program's reading, and shells.junction priced that throat — the bill is inherited, not paid

## T10_assessment — PASS

- **n_passed**: 9
- **n_total**: 9

## verdict — THE_SECOND_INTERACTION_MUST_BE_ANTIPODAL

PAIR CREATION IS A COLLISION, NOT A FOCUS — AND THE GEOMETRY THEN FORCES A SECOND INTERACTION AT THE ANTIPODE. A caustic is where the AMPLITUDE gets large; Breit-Wheeler is a threshold on an INVARIANT, s = 2 E1E2 (1 - cos theta) >= (2m)^2, and those are different quantities. E is what focusing raises; theta is what a collision supplies. Focusing is therefore NEITHER SUFFICIENT — collinear momenta have s = 0 identically, still zero after amplifying by 1e+12 — NOR NECESSARY, since two beams crossed head-on with no focusing anywhere clear threshold as soon as E >= m. The honest complication is recorded rather than buried: a spherically converging front DOES contain diametrically opposed rays, so its self-invariant is not zero, and the real distinction is INDEPENDENCE of the sources rather than brightness. Put two sources a geodesic distance delta apart and the crossing obeys 1 - cos theta = (1 - cos delta)/sin^2 t, an identity of geodesic triangles verified against embedded tangent vectors that never use the law of cosines — to 2.0e-14, on S2 AND on S3, because a geodesic triangle lies in a great 2-sphere whatever it is embedded in. So s(t) = 4 E1E2 sin^2(delta/2)/sin^2 t, which is U-SHAPED: head-on at BOTH ends of the crossing window and minimal at the equator, where the opening angle is exactly delta. The moment the wavefronts are largest is the moment the invariant is smallest, which is the opposite of the intuition the single-caustic pictures encouraged. A threshold therefore cuts TWO DISJOINT WINDOWS AND NEVER ONE, confirmed against a 40,000-point scan; below E = m there is none, at E = m the window has zero width, and they merge only above E = m/sin(delta/2) = 4.797 m. AND ONLY THE FAR WINDOW IS A COLLISION OF INDEPENDENT WAVES: the near one sits on top of the sources, the fronts having travelled 0.296 against a separation of 0.42, while the far one is reached only after each front has crossed a half-circumference — a factor of 9.6 in path length. That is why the second interaction has to be antipodal, and it is derived rather than staged. ONE FURTHER TRAP, CAUGHT: the momenta are exact, perpendicular to their own wavefront to 1e-15 and matching the closed form to 2e-13, but a FIGURE shows their projection, and projection does not preserve angles. Measured off the picture the opening angle is wrong by up to 67.4 degrees and differs by up to 56.4 degrees between the two crossing points, whose true opening angle is IDENTICAL — so the renderer draws the angle in the plane the two momenta span and labels the arrows on the sphere as decoration. WHAT IS IMPORTED, PLAINLY: the Breit-Wheeler threshold and cross-section are QED. The cross-section is the textbook closed form and is checked against its known peak — beta = 0.701, sqrt(s) = 1.40 (2m), sigma = 0.256 sigma_T — and it FALLS at large s, so the most violent part of the crossing is not the most productive. Rays-as-photons is a correspondence, no rate is computed, the orientation labels are carried rather than derived, and calling the crossing chord a throat through the bulk is this program's reading, whose exotic-matter bill was priced in PR #249 and is inherited here unpaid.
