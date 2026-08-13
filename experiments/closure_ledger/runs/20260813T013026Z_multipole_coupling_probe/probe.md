# multipole_coupling probe

_generated 2026-08-13T01:30:26.308859+00:00_

## T1_goal — PASS
- **question**: shell_junction found no throat-shell coupling at l = 0 and imported Birkhoff to explain it; does l >= 2 supply the coupling, and on what terms?
- **scope**: static Newtonian (Laplace) multipole problem, the weak-field limit of the junction problem; not Regge-Wheeler/Zerilli

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

## T3_birkhoff_is_the_ell_zero_case — PASS
- **rows**:
    - a=5, ell_0=0, ell_1=0.0177778, ell_2=0.00768, ell_3=0.00313469
    - a=20, ell_0=0, ell_1=0.00111111, ell_2=0.00012, ell_3=1.22449e-05
    - a=100, ell_0=0, ell_1=4.44444e-05, ell_2=9.6e-07, ell_3=1.95918e-08
- **b**: 2
- **ell_zero_is_exactly_zero_at_every_separation**: True
- **every_other_ell_couples**: True
- **ell_zero_is_separation_independent**: True
- **why**: the mutual stiffness carries the Laplacian eigenvalue ℓ(ℓ+1), which vanishes only for the constant mode — the decoupling of the spherical model is that zero, not a separate theorem about spheres

## T4_the_pure_dipole_is_not_a_translation — PASS
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

## T5_translation_invariance_is_exact — PASS
- **rows**:
    - displacement=0.02, exact_sphere=4.24074e-16, second_order_family=2.13309e-08, pure_P1=0.000266656, family_beats_pure_by=12500.9
    - displacement=0.05, exact_sphere=4.24074e-16, second_order_family=8.32738e-07, pure_P1=0.00166625, family_beats_pure_by=2000.93
    - displacement=0.1, exact_sphere=2.02142e-14, second_order_family=1.32952e-05, pure_P1=0.00666, family_beats_pure_by=500.931
- **the_exact_sphere_is_area_preserving**: True
- **the_truncated_family_is_preserving_to_order_d4**: True
- **the_pure_dipole_is_not**: True
- **translation_invariance_holds**: True
- **the_naive_test_does_not_test_it**: True

## T6_the_area_cost_is_the_exact_rational — PASS
- **rows**:
    - ell=0, value=2, rational=2, matches=True
    - ell=1, value=1.33333, rational=1.33333, matches=True
    - ell=2, value=1.6, rational=1.6, matches=True
    - ell=3, value=2, rational=2, matches=True
    - ell=4, value=2.44444, rational=2.44444, matches=True
    - ell=5, value=2.90909, rational=2.90909, matches=True
- **every_value_is_the_exact_rational**: True
- **it_grows_without_bound**: True

## T7_the_coupling_is_screened_by_separation — PASS
- **rows**:
    - ell=1, mutual_stiffness=0.0177778, geometric_suppression=0.4
    - ell=2, mutual_stiffness=0.00768, geometric_suppression=0.16
    - ell=3, mutual_stiffness=0.00313469, geometric_suppression=0.064
    - ell=4, mutual_stiffness=0.0012642, geometric_suppression=0.0256
    - ell=6, mutual_stiffness=0.000203588, geometric_suppression=0.004096
    - ell=8, mutual_stiffness=3.26546e-05, geometric_suppression=0.00065536
- **b_over_a**: 0.4
- **suppression_from_ell_1_to_ell_8**: 544.419
- **the_coupling_exists_for_every_ell_at_least_one**: True
- **but_it_falls_geometrically**: True
- **so_the_answer_has_two_halves**: ℓ ≥ 2 restores the coupling ℓ = 0 forbade, and screens it by (b/a)^ℓ — the modes that carry a spin-2 signal are the ones separation suppresses hardest

## T8_the_shear_response_is_an_input — PASS
- **rows**:
    - ell=2, fluid_shear_stiffness=0, elastic_shear_stiffness=10.0531, tension_stiffness=20.1062
    - ell=3, fluid_shear_stiffness=0, elastic_shear_stiffness=17.952, tension_stiffness=25.1327
    - ell=4, fluid_shear_stiffness=0, elastic_shear_stiffness=25.1327, tension_stiffness=30.7178
- **a_fluid_shell_has_no_shear_response**: True
- **an_elastic_one_needs_an_extra_input**: True
- **the_shear_modulus_is_never_fitted**: True
- **what_this_costs**: the previous round's conclusion that ℓ ≥ 2 is where the coupling lives is only half an answer: the coupling is there, but what the shell does with it depends on a constitutive choice spherical symmetry never had to make

## T9_assessment — PASS
- **n_passed**: 8
- **n_total**: 8

## verdict — BIRKHOFF_IS_A_VANISHING_EIGENVALUE

BIRKHOFF IS A VANISHING EIGENVALUE. The previous round measured that two concentric surfaces in a vacuum spherical model cannot talk, and imported Birkhoff's theorem to say why. The multipole form shows that import was not needed: the mutual stiffness of two concentric deformed shells is G m_b m_a l(l+1) (b/a)^l / (a (2l+1)^2), with the LAPLACIAN EIGENVALUE l(l+1) as its prefactor, so it vanishes identically at l = 0 and nowhere else. The decoupling is a property of the constant mode, not of spheres. The closed form agrees with brute-force double integration over both deformed surfaces — which never expands in multipoles — to 9.0e-06, and gives exactly zero at l = 0. So l >= 2 does supply the coupling the earlier round found forbidden, and the answer has a second half that matters as much as the first: the same formula screens it as (b/a)^l, a factor of 544 from l = 1 to l = 8 at b/a = 0.4. The multipoles that carry a spin-2 signal are precisely the ones separation suppresses hardest, so having the channel and being able to use it are different statements. Two things had to be got right before any of this could be trusted. A pure P1 deformation is NOT a translation past linear order: the area second variation is [2+l(l+1)]/(2l+1) exactly, which at l = 1 is 1.3333 and not zero, so a naive translation-invariance check would have reported a stiffness that is not there. Translation invariance does hold exactly — the rigid displacement carries induced l = 0 and l = 2 pieces, and along that family the exact displaced sphere is area-preserving to 4e-16 at every displacement. And the shear response is an INPUT: a perfect-fluid shell resists area change and nothing else, so resisting shape change at fixed area needs an elastic modulus no equation of state supplies. It is carried explicitly and never fitted. The coupling l >= 2 restores is real; what a shell would do with it is a constitutive choice spherical symmetry never had to make.
