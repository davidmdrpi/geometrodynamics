# probe — branch_coupling

_2026-08-15T03:51:37.455049+00:00_

## T1_goal — PASS

- **question**: is the throat identification part of the field problem, or a shift applied to the free branches after they are computed? and what is the primitive object -- is it indexed by a pair of branches?
- **scope**: a linear scalar field on a fixed Einstein static universe. The throat is still an identification map with a coupling kappa put in by hand, and PR #249's exotic-matter bill is inherited and unpaid. No backreaction, no stress tensor, no topology change, no rate. The two-throat cross term in T7 is a throat-throat interference, not roadmap step 3's two-source invariant.

## T2_the_closed_form_transfer_is_the_branch_sum — PASS

- **rows**:
    - chi=0.5, omega=0.4, series=0.151597-0.00355667j, closed=0.151597-0.00355667j, abs_error=6.44012e-17
    - chi=0.5, omega=1.7, series=0.196248-0.028117j, closed=0.196248-0.028117j, abs_error=5.97561e-16
    - chi=0.5, omega=3.3, series=-0.128911-0.0349941j, closed=-0.128911-0.0349941j, abs_error=6.93369e-16
    - chi=0.5, omega=6.1, series=-0.196938-0.0298564j, closed=-0.196938-0.0298564j, abs_error=9.16354e-16
    - chi=1.3, omega=0.4, series=0.0580173-0.00292181j, closed=0.0580173-0.00292181j, abs_error=9.07554e-17
    - chi=1.3, omega=1.7, series=-0.00212744-0.00905676j, closed=-0.00212744-0.00905676j, abs_error=1.97471e-16
    - chi=1.3, omega=3.3, series=0.0194271+0.0113007j, closed=0.0194271+0.0113007j, abs_error=4.08187e-16
    - chi=1.3, omega=6.1, series=-0.207231-0.10511j, closed=-0.207231-0.10511j, abs_error=1.01905e-15
    - chi=2.2, omega=0.4, series=0.0377705-0.00256512j, closed=0.0377705-0.00256512j, abs_error=4.82478e-17
    - chi=2.2, omega=1.7, series=-0.118753+0.0132731j, closed=-0.118753+0.0132731j, abs_error=1.92077e-16
    - chi=2.2, omega=3.3, series=-0.00344579-0.00604665j, closed=-0.00344579-0.00604665j, abs_error=8.21871e-16
    - chi=2.2, omega=6.1, series=-0.126593-0.0734155j, closed=-0.126593-0.0734155j, abs_error=2.72037e-15
- **worst_abs_error**: 2.72037e-15
- **the_series_is_the_closed_form**: True
- **poles**:
    - omega=1, scaling=[{'damping': 0.01, 'abs_T': 2.533075544560594, 'ratio_to_previous': None}, {'damping': 0.001, 'abs_T': 25.330300507723564, 'ratio_to_previous': 9.99982040098119}, {'damping': 0.0001, 'abs_T': 253.3029595655644, 'ratio_to_previous': 9.999998203271563}], residue=-0.0253302-1.29502e-08j, residue_real=-0.0253302, predicted_mode_over_2omega=-0.0253303, matches=True
    - omega=2, scaling=[{'damping': 0.01, 'abs_T': 1.3558064348636607, 'ratio_to_previous': None}, {'damping': 0.001, 'abs_T': 13.551713152542943, 'ratio_to_previous': 9.995315558378875}, {'damping': 0.0001, 'abs_T': 135.5164961142369, 'ratio_to_previous': 9.999953112113179}], residue=-0.0135516+4.14838e-08j, residue_real=-0.0135516, predicted_mode_over_2omega=-0.0135516, matches=True
    - omega=3, scaling=[{'damping': 0.01, 'abs_T': 1.8083701781508013, 'ratio_to_previous': None}, {'damping': 0.001, 'abs_T': 18.08023053461427, 'ratio_to_previous': 9.998080455575035}, {'damping': 0.0001, 'abs_T': 180.80195810737732, 'ratio_to_previous': 9.999980794560958}], residue=0.0180801+3.5144e-08j, residue_real=0.0180801, predicted_mode_over_2omega=0.0180802, matches=True
- **the_poles_are_the_esu_spectrum**: True
- **the_residues_are_the_mode_functions**: True
- **what_this_shows**: the image sum and the mode sum are the same function; the branch labels are a representation, not an approximation

## T3_solving_the_throat_resums_every_traversal — PASS

- **rows**:
    - kappa=0.05, omega=0.4, loop_gain=0.00290454, solved=0.000230884-0.00012233j, walked_series=0.000230884-0.00012233j, abs_error=0, post_processed_one_traversal=0.000230435-0.000121718j, relative_miss_of_one_traversal=0.00290454
    - kappa=0.05, omega=1.1, loop_gain=0.0107476, solved=0.00192977-1.94626e-07j, walked_series=0.00192977-1.94626e-07j, abs_error=2.1684e-19, post_processed_one_traversal=0.00194679-1.20438e-05j, relative_miss_of_one_traversal=0.0107476
    - kappa=0.05, omega=2.7, loop_gain=0.00486356, solved=-0.000106923-9.3356e-05j, walked_series=-0.000106923-9.3356e-05j, abs_error=3.83323e-20, post_processed_one_traversal=-0.000106592-9.275e-05j, relative_miss_of_one_traversal=0.00486356
    - kappa=0.05, omega=5.3, loop_gain=0.00171366, solved=9.55447e-05+0.000300796j, walked_series=9.55447e-05+0.000300796j, abs_error=5.58785e-20, post_processed_one_traversal=9.58186e-05+0.00030033j, relative_miss_of_one_traversal=0.00171366
    - kappa=0.2, omega=0.4, loop_gain=0.0116182, solved=0.000928946-0.000496766j, walked_series=0.000928946-0.000496766j, abs_error=0, post_processed_one_traversal=0.000921742-0.000486872j, relative_miss_of_one_traversal=0.0116182
    - kappa=0.2, omega=1.1, loop_gain=0.0429905, solved=0.00751863+0.00013181j, walked_series=0.00751863+0.00013181j, abs_error=5.42101e-20, post_processed_one_traversal=0.00778718-4.81752e-05j, relative_miss_of_one_traversal=0.0429905
    - kappa=0.2, omega=2.7, loop_gain=0.0194542, solved=-0.000431685-0.000380857j, walked_series=-0.000431685-0.000380857j, abs_error=2.29994e-19, post_processed_one_traversal=-0.000426369-0.000371j, relative_miss_of_one_traversal=0.0194542
    - kappa=0.2, omega=5.3, loop_gain=0.00685463, solved=0.000378849+0.00120879j, walked_series=0.000378849+0.00120879j, abs_error=0, post_processed_one_traversal=0.000383274+0.00120132j, relative_miss_of_one_traversal=0.00685463
    - kappa=0.45, omega=0.4, loop_gain=0.0261409, solved=0.00211055-0.00114646j, walked_series=0.00211055-0.00114646j, abs_error=8.94056e-19, post_processed_one_traversal=0.00207392-0.00109546j, relative_miss_of_one_traversal=0.0261409
    - kappa=0.45, omega=1.1, loop_gain=0.0967287, solved=0.0161951+0.000728599j, walked_series=0.0161951+0.000728599j, abs_error=3.47114e-18, post_processed_one_traversal=0.0175211-0.000108394j, relative_miss_of_one_traversal=0.0967287
    - kappa=0.45, omega=2.7, loop_gain=0.0437721, solved=-0.000986423-0.000886079j, walked_series=-0.000986423-0.000886079j, abs_error=0, post_processed_one_traversal=-0.000959331-0.00083475j, relative_miss_of_one_traversal=0.0437721
    - kappa=0.45, omega=5.3, loop_gain=0.0154229, solved=0.000839588+0.00274086j, walked_series=0.000839588+0.00274086j, abs_error=4.47028e-19, post_processed_one_traversal=0.000862368+0.00270297j, relative_miss_of_one_traversal=0.0154229
- **worst_abs_error**: 3.47114e-18
- **the_resolvent_is_the_traversal_sum**: True
- **worst_relative_miss_of_post_processing**: 0.0967287
- **post_processing_is_the_n_equals_zero_term**: True
- **what_this_shows**: writing the identification into the field problem is not a rearrangement of PR #254 — it adds every history that returns to the mouth

## T4_the_solve_adds_arrivals_post_processing_cannot — PASS

- **n_words_one_traversal**: 10
- **n_words_multi**: 25
- **reconstruction_relative_error**: 5.36486e-06
- **the_waveform_is_the_sum_over_history_words**: True
- **isolated_echoes**:
    - t=5.4, traversals=2, isolation=1.58319, word_sign=1, solved_level=0.00190197, control_level=3.91252e-16, contrast=4.86123e+12, field_sign=1, sign_agrees=True
    - t=9.88319, traversals=2, isolation=0.5, word_sign=-1, solved_level=0.00173424, control_level=5.24836e-16, contrast=3.30435e+12, field_sign=-1, sign_agrees=True
- **worst_echo_contrast**: 3.30435e+12
- **every_isolated_echo_stands_above_the_control**: True
- **every_echo_carries_its_maslov_word_sign**: True
- **amplitude_ladder**:
    - traversals=1, largest_amplitude=0.00499015, over_first=1
    - traversals=2, largest_amplitude=0.000240927, over_first=0.0482805
    - traversals=3, largest_amplitude=1.16321e-05, over_first=0.00233101
    - traversals=4, largest_amplitude=5.61603e-07, over_first=0.000112542
- **kappa**: 0.6
- **what_this_shows**: the extra arrivals are not a correction to the amplitudes of PR #254's ledger; they are events at times that ledger does not contain

## T5_closure_is_broadband_coherence — PASS

- **delay**: -14.6664
- **n_pairs**: 64
- **n_closed**: 3
- **worst_closed_coherence**: 1
- **best_other_coherence**: 0.0908531
- **closed_pairs_are_broadband_coherent**: True
- **every_other_pair_dephases**: True
- **contrast**: 11.0068
- **pairs_a_single_index_rule_would_select**: 9
- **the_condition_does_not_factorize**: True
- **the_amplitude_does_factorize**: True
- **what_this_shows**: the amplitude is an outer product over the pair index and the closure condition is not; that is why the pair is the primitive

## T6_the_condition_does_not_factorize — PASS

- **n_closed**: 3
- **pairs_a_single_index_rule_would_select**: 9
- **the_condition_does_not_factorize**: True
- **the_amplitude_does_factorize**: True
- **delay**: -14.6664

## T7_one_throat_is_rank_one_and_two_are_rank_two — PASS

- **cross_visibility_max**: 0.999973
- **cross_visibility_min**: -1
- **the_cross_term_is_a_full_fringe**: True
- **singular_values_throat_A**:
    - 1
    - 1.31982e-16
    - 4.63439e-17
    - 3.3065e-17
- **singular_values_throat_B**:
    - 1
    - 9.44176e-17
    - 5.18723e-17
    - 3.43986e-17
- **singular_values_both**:
    - 1
    - 0.541745
    - 1.04711e-16
    - 7.13296e-17
- **rank_one_throat**: 1
- **rank_two_throats**: 2
- **one_throat_is_rank_one**: True
- **two_throats_are_rank_two**: True
- **cross_term**: 4.72287e-06
- **cross_term_direct**: 4.72287e-06
- **cross_term_agrees**: True
- **why_it_vanishes_without_a_second_throat**: it is bilinear — one factor from each throat's pair sum — so removing either sets it to zero identically, not merely underdetermined; that is structural, not measured
- **scope**: a throat–throat interference term, not the two-source invariant of roadmap step 3

## T8_the_expansion_fails_at_the_eigenfrequencies — PASS

- **rows**:
    - kappa_critical=3.0465, omega_of_the_peak=5.9992, peak_transfer=0.328245, damping=0.08, kappa_c_ratio_to_previous=None, damping_ratio=None, exponent=None
    - kappa_critical=1.52368, omega_of_the_peak=6, peak_transfer=0.656308, damping=0.04, kappa_c_ratio_to_previous=1.99944, damping_ratio=2, exponent=0.999597
    - kappa_critical=0.761888, omega_of_the_peak=6, peak_transfer=1.31253, damping=0.02, kappa_c_ratio_to_previous=1.99987, damping_ratio=2, exponent=0.999905
    - kappa_critical=0.38095, omega_of_the_peak=6, peak_transfer=2.62501, damping=0.01, kappa_c_ratio_to_previous=1.99997, damping_ratio=2, exponent=0.999976
- **exponents**:
    - 0.999597
    - 0.999905
    - 0.999976
- **kappa_critical_scales_like_damping**: True
- **mean_exponent**: 0.999826
- **the_peak_sits_on_an_esu_eigenfrequency**: True
- **the_relative_miss_is_the_loop_gain**: |1/(1−L) − 1| / |1/(1−L)| = |L| exactly, so the error of post-processing IS the round-trip gain
- **on_and_off_resonance**:
    - omega=6, loop_gain=0.5, relative_miss_of_one_traversal=0.5
    - omega=6.5, loop_gain=0.0176638, relative_miss_of_one_traversal=0.0176638
- **resonance_is_where_post_processing_is_worst**: True
- **kappa_used**: 0.380944
- **what_this_shows**: the free-branch answer is an expansion in the loop gain, and the loop gain has no bound as the regulator is removed

## T9_assessment — PASS

- **n_passed**: 8
- **n_total**: 8

## verdict — THE_THROAT_IS_SOLVED_FOR_AND_THE_PRIMITIVE_IS_A_PAIR

YES -- AND IT WAS NOT A REARRANGEMENT OF PR #254. That round applied the identification phi(M+,t) = eta phi(M-,t+Delta) to the free branches after computing them, which gives one traversal by construction: a post-processing step cannot notice that what it re-emits will come back. Written into the field problem instead, the amplitude re-emitted at M- is driven by everything reaching M+ including its own return, and the solution carries a resolvent 1/(1 - L). FIRST, THE BRANCH SERIES SUMS IN CLOSED FORM: the short-way images all carry the Maslov factor +1 and the long-way images all carry -1, so the winding sum is two geometric series and equals (e^{-u chi} - e^{-u(2pi-chi)})/(1 - e^{-2pi u}), verified against the term-by-term sum to 2.7e-15. Its poles, as the regulator is removed, sit at omega = 1, 2, 3, ... -- the conformal ESU eigenfrequencies -- with residues equal to the mode functions over 2 omega. The image representation and the mode representation are the same function seen from two sides, which is the strongest statement so far that the branch labels are a REPRESENTATION rather than an approximation. SECOND, THE SOLVE IS THE SUM OVER EVERY TRAVERSAL: the resolvent agrees with an explicit walk over 400 traversals to 3.5e-18, and PR #254's answer is the n = 0 term, whose relative error is exactly the round-trip gain |L|. THIRD -- AND THIS IS THE SHARPEST PART -- THE SOLVE ADDS EVENTS, NOT AMPLITUDES. The solved waveform is the sum over history words to 5.4e-06, so the word enumeration is checked rather than asserted; at echo times ell_a + Delta + n(ell_c + Delta) + ell_b that no one-traversal word can reach, the solved field stands 3.3e+12 above the control, with amplitudes on the kappa^n ladder and signs equal to the product of every Maslov factor in the word. Those are arrivals at times PR #254's ledger does not contain. FOURTH, THE PRIMITIVE IS INDEXED BY A PAIR OF BRANCHES, AND FOR A GOOD REASON: K_ab carries the phase e^{-i omega (l_a + Delta + l_b)}, so PR #253's closure condition is EXACTLY the statement that K_ab is independent of omega -- closed pairs have band coherence 1.000 while every other pair dephases below 0.091. The amplitude factorizes over that index (K is rank one) and the CONDITION DOES NOT: 3 pairs close where any rule phrased on a alone and b alone would admit 9. FIFTH, the rank counts histories -- one throat rank one, two throats rank two -- and the throat-throat cross term is a full fringe, visibility -1.000 to 1.000, bilinear and therefore identically zero without either throat. That is the same shape as roadmap step 3's invariant and is explicitly NOT that invariant: these are throats, not sources. SIXTH, THERE IS A REGIME WHERE POST-PROCESSING IS NOT AN APPROXIMATION AT ALL: |T_d| peaks where 1 - e^{-2pi u} = 2 pi gamma, so the critical coupling falls linearly in the regulator -- measured exponent 0.9998 -- and the peak sits exactly on an ESU eigenfrequency. As the regulator is removed every coupling is critical somewhere, and there the one-traversal answer is the first term of a divergent series. WHAT IS STILL PUT IN: the throat remains an identification map with kappa by hand, the background is fixed, the field is linear, and when Delta + ell_c < 0 the loop is closed in time and 1/(1 - L) is a self-consistency condition rather than a history sum -- it has a unique solution exactly when the branch-resolved loop gain is subcritical, which is a bound on kappa and not a derivation of it. No backreaction, no stress tensor, no topology change, no rate, and no two-source invariant.
