# probe — field_solve

_2026-08-14T04:16:10.548809+00:00_

## T1_goal — PASS

- **question**: does a solved wave field, with the throat imposed as the same identification, reproduce the exact branch ledger of PR #253 -- arrival times, amplitudes, and phases the ray picture never had?
- **scope**: a linear scalar field on a fixed Einstein static universe. The throat is still an identification map, not a solution, with PR #249's exotic-matter bill inherited and unpaid. No backreaction, no stress tensor, no topology change, no rate, and no two-source invariant yet. Repeated throat traversals are PR #251's fixed point and are not redone here.

## T2_the_solve_matches_the_image_sum — PASS

- **rows**:
    - chi=0.4, peak_scale=1.35873, worst_abs_error=8.29125e-13, relative=6.10221e-13
    - chi=1.1, peak_scale=0.593704, worst_abs_error=1.24844e-15, relative=2.10279e-15
    - chi=2, peak_scale=0.581893, worst_abs_error=3.27516e-15, relative=5.62846e-15
    - chi=2.7, peak_scale=1.23804, worst_abs_error=1.76525e-14, relative=1.42585e-14
- **worst_abs_error**: 8.29125e-13
- **the_two_constructions_agree**: True
- **what_this_shows**: the solved field IS the image sum; the branch structure is exact, not an approximation

## T3_the_field_support_is_the_branch_ledger — PASS

- **rows**:
    - chi=0.7, n_peaks=4, n_branches=4, matches=[{'t_peak': 0.7002978395810671, 't_branch': 0.7, 'error': 0.0002978395810671053, 'branch': (0, 0)}, {'t_peak': 5.583445017679947, 't_branch': 5.583185307179586, 'error': 0.0002597105003605549, 'branch': (1, 0)}, {'t_peak': 6.983297296733512, 't_branch': 6.983185307179586, 'error': 0.00011198955392544008, 'branch': (0, 1)}, {'t_peak': 11.86644447483239, 't_branch': 11.866370614359173, 'error': 7.386047321666922e-05, 'branch': (1, 1)}]
    - chi=1.1, n_peaks=4, n_branches=4, matches=[{'t_peak': 1.0998966050559624, 't_branch': 1.1, 'error': 0.00010339494403766025, 'branch': (0, 0)}, {'t_peak': 5.183217952259336, 't_branch': 5.183185307179587, 'error': 3.264507974964914e-05, 'branch': (1, 0)}, {'t_peak': 7.382896062208407, 't_branch': 7.383185307179586, 'error': 0.00028924497117888137, 'branch': (0, 1)}, {'t_peak': 11.46621740941178, 't_branch': 11.466370614359173, 'error': 0.00015320494739334833, 'branch': (1, 1)}]
    - chi=1.9, n_peaks=4, n_branches=4, matches=[{'t_peak': 1.8997224359514686, 't_branch': 1.9, 'error': 0.000277564048531298, 'branch': (0, 0)}, {'t_peak': 4.38339212136383, 't_branch': 4.383185307179586, 'error': 0.00020681418424395304, 'branch': (1, 0)}, {'t_peak': 8.183350193049627, 't_branch': 8.183185307179587, 'error': 0.00016488587004026556, 'branch': (0, 1)}, {'t_peak': 10.666391578516274, 't_branch': 10.666370614359172, 'error': 2.096415710184374e-05, 'branch': (1, 1)}]
- **worst_time_error**: 0.00029784
- **every_branch_has_a_peak_and_no_peak_is_spurious**: True
- **peaks_land_on_the_ledger**: True
- **so_the_field_reproduces_the_ray_ledger**: True

## T4_the_phases_are_the_maslov_index — PASS

- **rows**:
    - chi=0.7, t=0.700298, crossings=0, field_sign=1, maslov_sign=1, agrees=True
    - chi=0.7, t=5.58345, crossings=1, field_sign=-1, maslov_sign=-1, agrees=True
    - chi=0.7, t=6.9833, crossings=2, field_sign=1, maslov_sign=1, agrees=True
    - chi=0.7, t=11.8664, crossings=3, field_sign=-1, maslov_sign=-1, agrees=True
    - chi=1.1, t=1.0999, crossings=0, field_sign=1, maslov_sign=1, agrees=True
    - chi=1.1, t=5.18322, crossings=1, field_sign=-1, maslov_sign=-1, agrees=True
    - chi=1.1, t=7.3829, crossings=2, field_sign=1, maslov_sign=1, agrees=True
    - chi=1.1, t=11.4662, crossings=3, field_sign=-1, maslov_sign=-1, agrees=True
    - chi=1.9, t=1.89972, crossings=0, field_sign=1, maslov_sign=1, agrees=True
    - chi=1.9, t=4.38339, crossings=1, field_sign=-1, maslov_sign=-1, agrees=True
    - chi=1.9, t=8.18335, crossings=2, field_sign=1, maslov_sign=1, agrees=True
    - chi=1.9, t=10.6664, crossings=3, field_sign=-1, maslov_sign=-1, agrees=True
- **every_sign_is_the_maslov_factor**: True
- **the_rule**: sign = (−1)^(number of focal crossings)
- **the_ray_ledger_could_not_supply_this**: path lengths give arrival times; the sign comes from the focusing history, which only a field carries

## T5_the_amplitude_is_the_shell_law — PASS

- **rows**:
    - chi=0.35, peak=1.85168, peak_times_sin_chi=0.634936
    - chi=0.7, peak=0.985593, peak_times_sin_chi=0.634936
    - chi=1.1, peak=0.712445, peak_times_sin_chi=0.634936
    - chi=1.5708, peak=0.634936, peak_times_sin_chi=0.634936
    - chi=2.1, peak=0.735553, peak_times_sin_chi=0.634936
    - chi=2.6, peak=1.23169, peak_times_sin_chi=0.634936
- **relative_spread**: 6.99423e-16
- **the_product_is_constant**: True
- **constant**: 0.634936
- **predicted_constant**: 0.634936
- **matches_one_over_four_pi_sin_chi**: True
- **so_the_shell_law_is_derived_here_not_imposed**: True

## T6_the_ledger_belongs_to_the_conformal_field — PASS

- **conformal**: peak=0.593698, worst_between_branches=2.36449e-08, ratio=3.98265e-08
- **minimal**: peak=0.61893, worst_between_branches=0.387897, ratio=0.626722
- **conformal_is_sharp**: True
- **minimal_is_not**: True
- **the_ledger_belongs_to_the_conformal_field**: True
- **why**: ω = n+1 is integer only with the conformal term ξR = 1; without it ω = √(n(n+2)) is irrational and the Green function has no images

## T7_the_throat_reproduces_the_closure_condition — PASS

- **rows**:
    - leg_in=1.2, leg_out=0.9, delay=-2.1, closes_at_zero=True, sign_of_the_closing_arrival=1, expected_sign=1, sign_agrees=True
    - leg_in=1.2, leg_out=5.38319, delay=-6.58319, closes_at_zero=True, sign_of_the_closing_arrival=-1, expected_sign=-1, sign_agrees=True
    - leg_in=1.2, leg_out=7.18319, delay=-8.38319, closes_at_zero=True, sign_of_the_closing_arrival=1, expected_sign=1, sign_agrees=True
    - leg_in=5.08319, leg_out=0.9, delay=-5.98319, closes_at_zero=True, sign_of_the_closing_arrival=-1, expected_sign=-1, sign_agrees=True
    - leg_in=5.08319, leg_out=5.38319, delay=-10.4664, closes_at_zero=True, sign_of_the_closing_arrival=1, expected_sign=1, sign_agrees=True
    - leg_in=5.08319, leg_out=7.18319, delay=-12.2664, closes_at_zero=True, sign_of_the_closing_arrival=-1, expected_sign=-1, sign_agrees=True
    - leg_in=7.48319, leg_out=0.9, delay=-8.38319, closes_at_zero=True, sign_of_the_closing_arrival=1, expected_sign=1, sign_agrees=True
    - leg_in=7.48319, leg_out=5.38319, delay=-12.8664, closes_at_zero=True, sign_of_the_closing_arrival=-1, expected_sign=-1, sign_agrees=True
    - leg_in=7.48319, leg_out=7.18319, delay=-14.6664, closes_at_zero=True, sign_of_the_closing_arrival=1, expected_sign=1, sign_agrees=True
- **closure_puts_an_arrival_on_the_emission_event**: True
- **every_closing_sign_is_eta_times_the_maslov_factors**: True
- **the_closure_condition_is_a_field_statement**: ℓ₁ + Δ + ℓ₂ = 0 is where a through-throat arrival lands back on its own emission event
- **and_the_field_adds_the_sign**: η times the two Maslov factors — the ray ledger carried no phase
- **what_is_still_put_in**: the throat is still an identification map, not a solution, and the self-consistent sum over repeated traversals is PR #251's fixed point and is not redone here

## T8_assessment — PASS

- **n_passed**: 7
- **n_total**: 7

## verdict — THE_FIELD_REPRODUCES_THE_LEDGER_AND_ADDS_ITS_PHASES

YES -- AND THE BRANCHES ARE EXACT SUPPORT, NOT STATIONARY-PHASE CONTRIBUTIONS. On the Einstein static universe the scalar Laplacian has eigenvalues -n(n+2) and R = 6, so the CONFORMALLY coupled massless field has omega^2 = n(n+2) + 1 = (n+1)^2: the frequencies are exactly the integers. The retarded Green function is therefore exactly periodic, and in closed form it is a sum of images, G = 1/(4 pi sin chi) [sum_k delta(t - chi - 2pi k) - sum_k delta(t + chi - 2pi k)]. The geometric-optics branches are the EXACT SUPPORT of the solved field; nothing has to be argued asymptotically. A truncated eigenmode sum -- which never sees an image -- agrees with the closed-form image sum -- which never sees a mode -- to 8.3e-13. THE FIELD'S SUPPORT IS PR #253's RAY LEDGER: peaks read off the SOLVED field land on the branch times to 3.0e-04, which is half a grid cell and therefore grid-limited rather than a disagreement, with every branch matched by a peak and no peak spurious. AND THE FIELD SUPPLIES WHAT THE RAY LEDGER STRUCTURALLY COULD NOT: every arrival carries a sign, and it is (-1)^m with m the number of focal crossings -- the antipode at t = pi, the source point again at t = 2pi, and so on. That is the MASLOV INDEX, and 12 of 12 signs agree. Path-length bookkeeping gives arrival times and has no way to produce a phase; this is the first quantity in the arc that the ray picture could not in principle have carried. THE AMPLITUDE IS PR #251's SHELL LAW, NOW DERIVED RATHER THAN IMPOSED: that round set A ~ 1/sin chi by conserving energy across a shell of area 4 pi sin^2 chi, and here peak * sin(chi) is the same constant at every chi to 7.0e-16, equal to 1/(4 pi sin chi). AND THE LEDGER BELONGS TO THE CONFORMAL FIELD SPECIFICALLY, which PR #253 had no way to notice: the MINIMALLY coupled field has omega = sqrt(n(n+2)), irrational, so no images and no sharp branches -- 63% of the peak amplitude sits BETWEEN the arrivals, against 4.0e-08 for the conformal field. Rays cannot tell the two apart. FINALLY THE THROAT REPRODUCES THE CLOSURE CONDITION: a through-throat contribution lands at l_1 + Delta + l_2, and setting Delta to minus a branch-pair sum -- exactly PR #253's closure condition -- puts an arrival back on the emission event, with the field adding the sign as eta times the two Maslov factors. WHAT IS STILL PUT IN: the throat remains an identification map rather than a solution, the background is fixed, the field is linear, repeated traversals are PR #251's fixed point and are not redone, and the two-source invariant that vanishes when a source is removed is NOT built here -- it is the next step, and the branch structure established above is why it will need a branch-resolved definition.
