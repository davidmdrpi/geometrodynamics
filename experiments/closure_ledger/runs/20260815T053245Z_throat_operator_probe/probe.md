# probe — throat_operator

_2026-08-15T05:32:45.478656+00:00_

## T1_goal — PASS

- **question**: does a genuinely flux-conserving two-mouth throat -- a self-adjoint extension, with reflection, transmission and a unitary boundary operator -- behave differently from PR #255's directional mouth relation? and what is its spectrum?
- **scope**: a linear scalar field on a fixed Einstein static universe. The throat is POINT-SUPPORTED: no interior, no proper length, and no delay -- Delta is not a parameter of a point extension. The boundary matrix A is a choice of four real parameters, not a derivation; PR #249 is what would fix it. No backreaction, no stress tensor, no topology change, no rate, no two-source invariant.

## T2_the_green_function_has_a_finite_part — PASS

- **rows**:
    - omega=0.37, chi=0.4, closed_form=0.189076, branch_series_limit=0.189076, residual_imaginary=6.53559e-08, abs_error=1.33643e-13
    - omega=0.37, chi=1.3, closed_form=0.0566811, branch_series_limit=0.0566811, residual_imaginary=5.16584e-08, abs_error=1.13382e-13
    - omega=0.37, chi=2.6, closed_form=0.0334809, branch_series_limit=0.0334809, residual_imaginary=4.37576e-08, abs_error=1.00815e-13
    - omega=1.63, chi=0.4, closed_form=0.216089, branch_series_limit=0.216089, residual_imaginary=4.41006e-07, abs_error=8.539e-13
    - omega=1.63, chi=1.3, closed_form=-0.0125391, branch_series_limit=-0.0125391, residual_imaginary=1.47058e-07, abs_error=1.59303e-13
    - omega=1.63, chi=2.6, closed_form=-0.12994, branch_series_limit=-0.12994, residual_imaginary=2.34498e-07, abs_error=9.40942e-13
    - omega=2.4, chi=0.4, closed_form=0.0628065, branch_series_limit=0.0628065, residual_imaginary=4.99236e-07, abs_error=4.35624e-13
    - omega=2.4, chi=1.3, closed_form=-0.0831472, branch_series_limit=-0.0831472, residual_imaginary=3.87515e-08, abs_error=3.0885e-13
    - omega=2.4, chi=2.6, closed_form=0.156391, branch_series_limit=0.156391, residual_imaginary=1.36108e-07, abs_error=8.87734e-13
    - omega=5.21, chi=0.4, closed_form=-0.329837, branch_series_limit=-0.329837, residual_imaginary=1.46934e-06, abs_error=6.33893e-12
    - omega=5.21, chi=1.3, closed_form=0.0227861, branch_series_limit=0.0227861, residual_imaginary=1.52288e-07, abs_error=5.4292e-13
    - omega=5.21, chi=2.6, closed_form=-0.0792028, branch_series_limit=-0.0792028, residual_imaginary=4.50268e-07, abs_error=2.20288e-12
- **worst_abs_error**: 6.33893e-12
- **the_closed_form_is_the_branch_series**: True
- **finite_part**:
    - omega=0.37, chi=0.01, finite_part=-0.0126634, g=-0.0127414, error=7.79758e-05
    - omega=0.37, chi=0.001, finite_part=-0.0127336, g=-0.0127414, error=7.814e-06
    - omega=0.37, chi=0.0001, finite_part=-0.0127406, g=-0.0127414, error=7.81565e-07
    - omega=1.63, chi=0.01, finite_part=0.055205, g=0.0561311, error=0.000926061
    - omega=1.63, chi=0.001, finite_part=0.0560386, g=0.0561311, error=9.24673e-05
    - omega=1.63, chi=0.0001, finite_part=0.0561218, g=0.0561311, error=9.24533e-06
    - omega=2.4, chi=0.01, finite_part=-0.0642093, g=-0.0620551, error=0.00215421
    - omega=2.4, chi=0.001, finite_part=-0.062271, g=-0.0620551, error=0.000215871
    - omega=2.4, chi=0.0001, finite_part=-0.0620767, g=-0.0620551, error=2.15915e-05
- **convergence_ratios_per_decade**:
    - 9.97898
    - 9.99789
    - 10.015
    - 10.0015
    - 9.97914
    - 9.99794
- **the_remainder_is_first_order_in_chi**: True
- **what_this_shows**: the coincidence divergence is the universal 1/(4πχ), so the finite part g(ω) exists and a point interaction is definable

## T3_the_boundary_operator_is_unitary_with_both_channels — PASS

- **rows**:
    - alpha1=0.05, alpha2=-0.011, beta=0.03+0.02j, scale=0.05, reflection=0.654173, transmission=0.756345, sum_of_squares=1, unitarity_defect=2.23367e-16, cayley_round_trip=6.93889e-18, reciprocal=False
    - alpha1=0.05, alpha2=-0.011, beta=0.03+0.02j, scale=0.1, reflection=0.816014, transmission=0.578032, sum_of_squares=1, unitarity_defect=2.2205e-16, cayley_round_trip=7.75792e-18, reciprocal=False
    - alpha1=0.05, alpha2=-0.011, beta=0.03+0.02j, scale=0.2, reflection=0.940865, transmission=0.338783, sum_of_squares=1, unitarity_defect=1.14054e-16, cayley_round_trip=8.39866e-18, reciprocal=False
    - alpha1=0.2, alpha2=0.2, beta=0.15, scale=0.05, reflection=0.8, transmission=0.6, sum_of_squares=1, unitarity_defect=4.44208e-16, cayley_round_trip=7.47342e-17, reciprocal=True
    - alpha1=0.2, alpha2=0.2, beta=0.15, scale=0.1, reflection=0.675725, transmission=0.737154, sum_of_squares=1, unitarity_defect=2.22091e-16, cayley_round_trip=7.78923e-17, reciprocal=True
    - alpha1=0.2, alpha2=0.2, beta=0.15, scale=0.2, reflection=0.691905, transmission=0.721988, sum_of_squares=1, unitarity_defect=4.47634e-17, cayley_round_trip=1.38778e-17, reciprocal=True
    - alpha1=-0.4, alpha2=0.07, beta=-0.09+0.31j, scale=0.05, reflection=0.971298, transmission=0.237866, sum_of_squares=1, unitarity_defect=1.13724e-16, cayley_round_trip=2.00148e-16, reciprocal=False
    - alpha1=-0.4, alpha2=0.07, beta=-0.09+0.31j, scale=0.1, reflection=0.896889, transmission=0.442256, sum_of_squares=1, unitarity_defect=4.44089e-16, cayley_round_trip=1.75542e-16, reciprocal=False
    - alpha1=-0.4, alpha2=0.07, beta=-0.09+0.31j, scale=0.2, reflection=0.713985, transmission=0.700161, sum_of_squares=1, unitarity_defect=2.22784e-16, cayley_round_trip=1.46948e-16, reciprocal=False
- **worst_unitarity_defect**: 4.44208e-16
- **worst_sum_of_squares_defect**: 4.44089e-16
- **the_cayley_transform_is_unitary**: True
- **every_mouth_conserves**: True
- **both_channels_are_present**: True
- **the_family_is_U2**: four real parameters — two self-energies and a complex mouth-to-mouth amplitude
- **pr255_in_the_same_language**:
    - kappa=0.3, reflection=0, transmission=0.3, sum_of_squares=0.09, in_U2=False
    - kappa=0.6, reflection=0, transmission=0.6, sum_of_squares=0.36, in_U2=False
    - kappa=1, reflection=0, transmission=1, sum_of_squares=1, in_U2=True
- **pr255_is_outside_U2_unless_kappa_is_one**: True

## T4_flux_conservation_is_exactly_hermiticity — PASS

- **n_draws**: 200
- **worst_relative_net_flux**: 1.83305e-16
- **flux_is_conserved_identically**: True
- **worst_pairwise_imbalance**: 1.72379e-16
- **what_one_mouth_absorbs_the_other_emits**: True
- **control_directional_worst_net_flux**: 0.998555
- **control_median_net_flux**: 0.538982
- **the_control_does_not_conserve**: True
- **control_anti_hermitian_part**: 1.66667
- **the_identity**: Im(q† A q) = 0 for all q  ⟺  A = A†

## T5_self_adjointness_makes_the_spectrum_real — PASS

- **rows**:
    - alpha1=0.05, alpha2=0.05, beta=0.03, hermitian=True, worst_relative_imaginary_part_of_det=0, n_real_roots=12, n_complex_seeds_converged=12, n_off_axis=0, worst_abs_imaginary=6.35338e-24, nothing_off_the_axis=True
    - alpha1=0.2, alpha2=-0.13, beta=0.15+0.07j, hermitian=True, worst_relative_imaginary_part_of_det=6.10623e-16, n_real_roots=11, n_complex_seeds_converged=11, n_off_axis=0, worst_abs_imaginary=3.88633e-18, nothing_off_the_axis=True
    - alpha1=-0.4, alpha2=0.07, beta=-0.09+0.31j, hermitian=True, worst_relative_imaginary_part_of_det=1.48525e-15, n_real_roots=11, n_complex_seeds_converged=8, n_off_axis=0, worst_abs_imaginary=4.48791e-18, nothing_off_the_axis=True
- **worst_relative_imaginary_part**: 1.48525e-15
- **the_secular_function_is_real_on_the_axis**: True
- **every_root_found_is_real**: True
- **nothing_off_the_axis**: True
- **worst_abs_imaginary_over_all_seeds**: 4.48791e-18
- **control_directional**: anti_hermitian_part=1.66667, worst_relative_imaginary_part=1, secular_is_real_on_the_axis=False, n_roots=9, n_off_axis=9, n_growing=2, worst_abs_imaginary=0.683863
- **the_control_fails_both_tests**: True
- **and_the_control_is_unstable_even_at_unit_transmission**: True
- **what_this_retires**: PR #255 found its poles off the real axis and had to separate three thresholds to say what that meant; a flux-conserving throat has a real spectrum for every coupling, so the instability was the model's non-conservation, not a throat

## T6_the_coupled_spectrum_interlaces_the_free_one — PASS

- **by_channel**:
    - gap=1, symmetric=1.67978, antisymmetric=1.45153
    - gap=2, symmetric=2.7234, antisymmetric=2.39818
    - gap=3, symmetric=3.61906, antisymmetric=3.50647
    - gap=4, symmetric=4.50474, antisymmetric=4.5782
    - gap=5, symmetric=5.5177, antisymmetric=5.55752
    - gap=6, symmetric=6.57067, antisymmetric=6.48235
    - gap=7, symmetric=7.58463, antisymmetric=7.46956
    - gap=8, symmetric=8.5379, antisymmetric=8.51074
- **n_roots**: 16
- **roots_per_gap**: 1=2, 2=2, 3=2, 4=2, 5=2, 6=2, 7=2, 8=2
- **exactly_two_per_gap**: True
- **every_root_strictly_between_free_ones**: True
- **both_channel_functions_are_monotone**: True
- **decoupling**:
    - strength=1, inverse_strength=1, worst_shift=0.49526, exponent=None
    - strength=10, inverse_strength=0.1, worst_shift=0.348831, exponent=0.152219
    - strength=100, inverse_strength=0.01, worst_shift=0.0602253, exponent=0.762835
    - strength=1000, inverse_strength=0.001, worst_shift=0.00621173, exponent=0.986567
    - strength=10000, inverse_strength=0.0001, worst_shift=0.000622578, exponent=0.999018
- **shift_exponents**:
    - 0.152219
    - 0.762835
    - 0.986567
    - 0.999018
- **asymptotic_exponents**:
    - 0.986567
    - 0.999018
- **the_shift_goes_like_one_over_the_boundary_norm**: True
- **free_spectrum_recovered**: True
- **off_is_large_A_not_small_A**: the diagonal of A is an inverse scattering length, so α → ∞ decouples and α = 0 is a resonant throat
- **what_this_shows**: the throat interlaces the ESU spectrum two per gap and returns it when switched off; the coupled problem is a rank-two perturbation of a self-adjoint operator, which is why

## T7_the_free_spectrum_returns_on_decoupling — PASS

- **decoupling**:
    - strength=1, inverse_strength=1, worst_shift=0.49526, exponent=None
    - strength=10, inverse_strength=0.1, worst_shift=0.348831, exponent=0.152219
    - strength=100, inverse_strength=0.01, worst_shift=0.0602253, exponent=0.762835
    - strength=1000, inverse_strength=0.001, worst_shift=0.00621173, exponent=0.986567
    - strength=10000, inverse_strength=0.0001, worst_shift=0.000622578, exponent=0.999018
- **asymptotic_exponents**:
    - 0.986567
    - 0.999018
- **the_shift_goes_like_one_over_the_boundary_norm**: True
- **free_spectrum_recovered**: True
- **off_is_large_A_not_small_A**: the diagonal of A is an inverse scattering length, so α → ∞ decouples and α = 0 is a resonant throat

## T8_where_pr255_sits — PASS

- **boundary_entry_by_frequency**:
    - 0.891663+3.21186j
    - -3.10583+1.21033j
    - -1.88871-2.74661j
    - 2.58522-2.10422j
- **the_boundary_matrix_depends_on_frequency**: True
- **only_through_its_phase**: True
- **rows**:
    - omega=1.3, reflection_entries=[np.float64(0.0), np.float64(0.0)], anti_hermitian_norm=1.66667, hermitian_norm=1.66667, secular_relative_imaginary=0.964481
    - omega=2.77, reflection_entries=[np.float64(0.0), np.float64(0.0)], anti_hermitian_norm=1.66667, hermitian_norm=1.66667, secular_relative_imaginary=0.324235
    - omega=4.11, reflection_entries=[np.float64(0.0), np.float64(0.0)], anti_hermitian_norm=1.66667, hermitian_norm=1.66667, secular_relative_imaginary=0.890329
    - omega=5.6, reflection_entries=[np.float64(0.0), np.float64(0.0)], anti_hermitian_norm=1.66667, hermitian_norm=1.66667, secular_relative_imaginary=0.596812
- **no_reflection_channel**: True
- **not_hermitian_at_any_frequency**: True
- **anti_hermitian_part_is_the_whole_matrix**: True
- **what_is_not_wrong_with_it**: the resolvent of PR #255 is exact for the model it posed; these are statements about which model
- **what_a_self_adjoint_extension_forbids**: a frequency-dependent boundary matrix — A is a boundary condition, not a dynamical response, so the delay Δ has no place in a point extension and does not survive into it

## T9_assessment — PASS

- **n_passed**: 8
- **n_total**: 8

## verdict — A_CONSERVING_THROAT_HAS_A_REAL_SPECTRUM

YES -- AND THE DIFFERENCE IS THE ONE THAT MATTERED. PR #255 owed a boundary operator and this is it: a point-supported throat is a SELF-ADJOINT EXTENSION of the Laplacian on S3 minus two points, and von Neumann's theorem parametrizes those by a unitary between the deficiency spaces -- U(2) -- equivalently, by Krein's formula, a Hermitian 2x2 boundary matrix A. Everything below follows from that one substitution. FIRST, IT IS DEFINABLE AT ALL: the free Green function has the closed form sin(omega(pi-chi))/(4 pi sin chi sin(pi omega)), real on the real axis with poles exactly at the free spectrum omega = n+1, agreeing with PR #255's branch series to 6.3e-12; and its short-distance split is 1/(4 pi chi) + g(omega) + O(chi) with g = -(omega/4 pi) cot(pi omega), the remainder first order in chi. The divergence is the universal Coulomb one, so the subtraction is not a choice. SECOND, THE BOUNDARY OPERATOR IS UNITARY AND HAS BOTH CHANNELS: the Cayley transform of any Hermitian A is unitary to 4.4e-16 and inverts back, with reflection on the diagonal, transmission off it, and |r|^2 + |t|^2 = 1 at each mouth to 4.4e-16. PR #255's model in the same language has r = 0 and |t| = kappa: outside U(2) unless kappa = 1. THIRD, FLUX CONSERVATION IS EXACTLY HERMITICITY. The current through a small sphere at mouth j is Im(q_j* phi_j^reg), independent of the sphere, so the total absorbed is Im(q^dag A q) -- zero for every q when A = A^dag, measured at 1.8e-16 over 200 random draws, with a purely off-diagonal throat moving flux from one mouth to the other exactly. FOURTH -- AND THIS IS WHAT THE ROUND EXISTS FOR -- THE SPECTRUM IS REAL FOR EVERY COUPLING. Gamma is real symmetric on the axis, so M = A - Gamma is Hermitian there and det M is a real function; Newton from a grid of complex seeds, the same method PR #255 used to find its poles, converges only onto the axis: 0 off-axis roots, worst |Im omega| = 4.5e-18. The directional control gives 9 roots of which 9 are off-axis and 2 growing, worst |Im omega| = 0.684 -- and it is unstable even at unit transmission, so it is the DIRECTIONALITY and not the loss. PR #255's instability was its own non-conservation, and the three thresholds that round had to separate collapse here to one statement: a conserving throat cannot ring up. FIFTH, THE COUPLED SPECTRUM INTERLACES THE FREE ONE. For an exchange-symmetric pair the secular equation splits into g + G_d = alpha + beta and g - G_d = alpha - beta, both monotone from -infinity to +infinity across every unit gap, so there are EXACTLY TWO coupled frequencies strictly between consecutive free ones -- verified over 8 gaps. SIXTH, THE FREE SPECTRUM RETURNS WHEN THE THROAT IS SWITCHED OFF, and off is ||A|| -> infinity rather than A -> 0, because the diagonal of A is an INVERSE scattering length and alpha = 0 is a resonant throat rather than no throat. The shift to the nearest free frequency falls like 1/||A||, exponent 0.9990. WHAT IS STILL PUT IN: the boundary matrix. A is four real numbers chosen, not derived, and PR #249 is what would fix them from a matter model; nothing here computes the exotic-matter bill. The throat is POINT-SUPPORTED, so it has no interior and no proper length, and the delay Delta of PRs #251-#255 is not a parameter of a point extension and does not survive into one -- a real loss of structure relative to those rounds, stated rather than hidden. No backreaction, no stress tensor, no topology change, no rate, and no two-source invariant, which is the next step.
