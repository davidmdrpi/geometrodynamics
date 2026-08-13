# multipole_coupling probe

_generated 2026-08-13T02:26:37.639716+00:00_

## T1_goal — PASS
- **question**: shell_junction found no throat-shell coupling at l = 0 in GR; in the static Newtonian analogue, where does the coupling start and on what terms?
- **scope**: static Newtonian (Laplace) two-shell model; establishes the shell-theorem/multipole structure. Birkhoff is a GR result and remains what shell_junction relies on — this is its Newtonian analogue, not that theorem. Not Regge-Wheeler/Zerilli, no quasinormal frequencies.

## T2_the_closed_form_against_brute_force — PASS
- **rows**:
    - ell=0, closed_form=0, brute_force=1.73472e-12, relative_error=1.73472e-12, laplacian_eigenvalue=0
    - ell=1, closed_form=0.0177778, brute_force=0.0177779, relative_error=4.79918e-06, laplacian_eigenvalue=2
    - ell=2, closed_form=0.00768, brute_force=0.00768004, relative_error=5.71494e-06, laplacian_eigenvalue=6
    - ell=3, closed_form=0.00313469, brute_force=0.00313472, relative_error=7.19198e-06, laplacian_eigenvalue=12
    - ell=4, closed_form=0.0012642, brute_force=0.00126421, relative_error=8.96083e-06, laplacian_eigenvalue=20
- **b**: 2
- **a**: 5
- **worst_relative_error**: 8.96083e-06
- **the_closed_form_is_confirmed**: True
- **the_prefactor_is_ell_times_ell_plus_one**: True
- **formula**: G m_b m_a ℓ(ℓ+1) (b/a)^ℓ / (a (2ℓ+1)²)

## T3_the_ell_zero_decoupling_is_a_zero_eigenvalue — PASS
- **rows**:
    - a=5, ell_0=0, ell_1=0.0177778, ell_2=0.00768, ell_3=0.00313469
    - a=20, ell_0=0, ell_1=0.00111111, ell_2=0.00012, ell_3=1.22449e-05
    - a=100, ell_0=0, ell_1=4.44444e-05, ell_2=9.6e-07, ell_3=1.95918e-08
- **b**: 2
- **ell_zero_is_exactly_zero_at_every_separation**: True
- **every_other_ell_couples**: True
- **ell_zero_is_separation_independent**: True
- **why**: the mutual stiffness carries the angular-Laplacian eigenvalue ℓ(ℓ+1), which vanishes only for the constant mode; in this Newtonian model that zero IS the ℓ = 0 decoupling
- **this_is_not_birkhoff**: Birkhoff is a GR theorem and remains what shells.junction relies on; this is its static Newtonian analogue

## T4_the_translation_mode_does_not_couple — PASS
- **rows**:
    - displacement=0.1, energy=-0.2, relative_drift=0, inner_stays_inside=True
    - displacement=0.8, energy=-0.2, relative_drift=2.77556e-16, inner_stays_inside=True
    - displacement=2.5, energy=-0.2, relative_drift=1.66533e-15, inner_stays_inside=True
- **reference_energy**: -0.2
- **newtonian_prediction**: -0.2
- **shell_theorem_holds**: True
- **rigid_translation_cross_derivative**: 8.32667e-13
- **pure_P1_shape_coupling**: 0.0177778
- **the_translation_mode_does_not_couple**: True
- **the_shape_mode_does**: True
- **so_coupling_starts_at_ell_two**: True
- **why_the_area_test_was_not_enough**: translation invariance of the area is a statement about the surface, not about the mutual gravitational energy; the energy control is the one that decides whether ℓ = 1 couples, and it says it does not

## T5_the_pure_dipole_is_not_a_translation — PASS
- **rows**:
    - ell=0, numeric=2, closed_form=2, relative_error=2.86736e-10
    - ell=1, numeric=1.33333, closed_form=1.33333, relative_error=2.86736e-10
    - ell=2, numeric=1.6, closed_form=1.6, relative_error=3.56262e-08
    - ell=3, numeric=2, closed_form=2, relative_error=7.35464e-09
    - ell=4, numeric=2.44444, closed_form=2.44444, relative_error=2.14904e-08
    - ell=5, numeric=2.90909, closed_form=2.90909, relative_error=6.47115e-09
- **worst_relative_error**: 3.56262e-08
- **the_closed_form_is_confirmed**: True
- **the_dipole_area_cost_is_not_zero**: True
- **dipole_second_variation**: 1.33333
- **so_a_pure_P1_test_would_be_wrong**: True

## T6_translation_invariance_is_exact — PASS
- **rows**:
    - displacement=0.02, exact_sphere=4.24074e-16, second_order_family=2.13309e-08, pure_P1=0.000266656, family_beats_pure_by=12500.9
    - displacement=0.05, exact_sphere=4.24074e-16, second_order_family=8.32738e-07, pure_P1=0.00166625, family_beats_pure_by=2000.93
    - displacement=0.1, exact_sphere=2.02142e-14, second_order_family=1.32952e-05, pure_P1=0.00666, family_beats_pure_by=500.931
- **the_exact_sphere_is_area_preserving**: True
- **the_truncated_family_is_preserving_to_order_d4**: True
- **the_pure_dipole_is_not**: True
- **translation_invariance_holds**: True
- **the_naive_test_does_not_test_it**: True

## T7_the_area_cost_is_the_exact_rational — PASS
- **rows**:
    - ell=0, value=2, rational=2, matches=True
    - ell=1, value=1.33333, rational=1.33333, matches=True
    - ell=2, value=1.6, rational=1.6, matches=True
    - ell=3, value=2, rational=2, matches=True
    - ell=4, value=2.44444, rational=2.44444, matches=True
    - ell=5, value=2.90909, rational=2.90909, matches=True
- **every_value_is_the_exact_rational**: True
- **it_grows_without_bound**: True

## T8_the_coupling_is_screened_by_separation — PASS
- **rows**:
    - ell=1, mutual_stiffness=0.0177778, geometric_suppression=0.4
    - ell=2, mutual_stiffness=0.00768, geometric_suppression=0.16
    - ell=3, mutual_stiffness=0.00313469, geometric_suppression=0.064
    - ell=4, mutual_stiffness=0.0012642, geometric_suppression=0.0256
    - ell=6, mutual_stiffness=0.000203588, geometric_suppression=0.004096
    - ell=8, mutual_stiffness=3.26546e-05, geometric_suppression=0.00065536
- **b_over_a**: 0.4
- **suppression_from_ell_1_to_ell_8**: 544.419
- **every_shape_mode_couples**: True
- **but_the_ell_1_entry_is_a_shape_mode_not_a_translation**: the rigid translation has zero mutual coupling by the shell theorem — see measure_the_translation_mode_does_not_couple
- **but_it_falls_geometrically**: True
- **so_the_answer_has_two_halves**: genuine coupling starts at ℓ = 2, and the same formula screens it by (b/a)^ℓ — the modes that carry a spin-2 signal are the ones separation suppresses hardest

## T9_the_shear_response_is_an_input — PASS
- **rows**:
    - ell=2, fluid_shear_stiffness=0, elastic_shear_stiffness=10.0531, tension_stiffness=20.1062
    - ell=3, fluid_shear_stiffness=0, elastic_shear_stiffness=17.952, tension_stiffness=25.1327
    - ell=4, fluid_shear_stiffness=0, elastic_shear_stiffness=25.1327, tension_stiffness=30.7178
- **a_fluid_shell_has_no_shear_response**: True
- **an_elastic_one_needs_an_extra_input**: True
- **the_shear_modulus_is_never_fitted**: True
- **what_this_costs**: the previous round's conclusion that ℓ ≥ 2 is where the coupling lives is only half an answer: the coupling is there, but what the shell does with it depends on a constitutive choice spherical symmetry never had to make

## T10_assessment — PASS
- **n_passed**: 9
- **n_total**: 9

## verdict — COUPLING_STARTS_AT_ELL_TWO

COUPLING STARTS AT L = 2. In the static Newtonian two-shell model the monopole mutual stiffness vanishes while higher angular multipoles couple, with the coupling suppressed geometrically by separation. The closed form is G m_b m_a l(l+1) (b/a)^l / (a (2l+1)^2), agreeing with brute-force double integration over both deformed surfaces — which never expands in multipoles — to 9.0e-06, and exactly zero at l = 0. The prefactor is the angular-Laplacian eigenvalue, so the l = 0 decoupling IS that zero eigenvalue. This is the Newtonian analogue of what shell_junction established in GR; BIRKHOFF IS A GR THEOREM AND REMAINS WHAT THAT ROUND RELIES ON, and nothing here replaces it. Where the coupling starts is the part an earlier draft got wrong, and the correction matters. That draft concluded 'everything l >= 1 couples', having checked translation invariance of the AREA — the wrong quantity. Run on the mutual ENERGY the two disagree: a rigidly translated inner sphere leaves the mutual energy at exactly -G m_b m_a / a, Newton's shell theorem, held to 1e-15 out to d = 2.5, so the cross-derivative with respect to rigid displacements is 8.3e-13. THE TRANSLATION MODE DOES NOT COUPLE. The pure P1 SHAPE mode is a different object and does, at 1.78e-02. So the honest ordering is that l = 0 decouples by the vanishing eigenvalue, the l = 1 translation mode decouples by the shell theorem, and genuine coupling starts at l = 2 — which is what PR #249 guessed and this establishes. Both of this round's errors are the same mistake: a pure P1 deformation is not a translation past linear order, its area second variation being 1.3333 and not zero, and a zero-mode test has to be run on the quantity the claim is about. The second half of the answer is that the same formula screens the coupling as (b/a)^l — a factor of 544 from l = 1 to l = 8 at b/a = 0.4 — so the modes carrying a spin-2 signal are the ones separation suppresses hardest. And the shear response is an INPUT: a perfect-fluid shell resists area change and nothing else, so resisting shape change at fixed area needs an elastic modulus no equation of state supplies. It is carried explicitly and never fitted.
