# probe — finite_mouth

_2026-08-19T04:46:54.899359+00:00_

## T1_goal — PASS

- **question**: PR #260 gated the roadmap on one question: does its negative mode survive a finite-radius mouth? Replace the point mouths by spheres of radius a, coupled through one channel each, and settle it.
- **scope**: a linear conformally coupled scalar on a fixed Einstein static universe. The mouths are spheres in a FIXED ambient with MONOPOLE coupling -- not a solved neck geometry -- so the l >= 1 content on each sphere is dropped, quantified in T8. No backreaction.

## T2_the_mean_value_identities — PASS

- **cross**:
    - radius=0.05, separation=0.8, lambda=0.16, quadrature=0.0939949, prediction=0.0939949, relative_error=2.06393e-11
    - radius=0.05, separation=1.3, lambda=0.16, quadrature=0.0583577, prediction=0.0583577, relative_error=2.05863e-11
    - radius=0.05, separation=2.2, lambda=0.16, quadrature=0.0380771, prediction=0.0380771, relative_error=2.05651e-11
    - radius=0.05, separation=0.8, lambda=2.89, quadrature=0.101943, prediction=0.101943, relative_error=2.06547e-11
    - radius=0.05, separation=1.3, lambda=2.89, quadrature=-0.00111029, prediction=-0.00111029, relative_error=1.82832e-11
    - radius=0.05, separation=2.2, lambda=2.89, quadrature=-0.121512, prediction=-0.121512, relative_error=2.05585e-11
    - radius=0.15, separation=0.8, lambda=0.16, quadrature=0.0942587, prediction=0.0942587, relative_error=2.12853e-11
    - radius=0.15, separation=1.3, lambda=0.16, quadrature=0.0585215, prediction=0.0585215, relative_error=2.07865e-11
    - radius=0.15, separation=2.2, lambda=0.16, quadrature=0.038184, prediction=0.038184, relative_error=2.05957e-11
    - radius=0.15, separation=0.8, lambda=2.89, quadrature=0.101301, prediction=0.101301, relative_error=2.14245e-11
    - radius=0.15, separation=1.3, lambda=2.89, quadrature=-0.0011033, prediction=-0.0011033, relative_error=1.20675e-13
    - radius=0.15, separation=2.2, lambda=2.89, quadrature=-0.120746, prediction=-0.120746, relative_error=2.05321e-11
    - radius=0.35, separation=0.8, lambda=0.16, quadrature=0.0955952, prediction=0.0955952, relative_error=2.52622e-11
    - radius=0.35, separation=1.3, lambda=0.16, quadrature=0.0593513, prediction=0.0593513, relative_error=2.18641e-11
    - radius=0.35, separation=2.2, lambda=0.16, quadrature=0.0387254, prediction=0.0387254, relative_error=2.07523e-11
    - radius=0.35, separation=0.8, lambda=2.89, quadrature=0.0981002, prediction=0.0981002, relative_error=2.62085e-11
    - radius=0.35, separation=1.3, lambda=2.89, quadrature=-0.00106844, prediction=-0.00106844, relative_error=1.00992e-10
    - radius=0.35, separation=2.2, lambda=2.89, quadrature=-0.116931, prediction=-0.116931, relative_error=2.03933e-11
- **self**:
    - radius=0.05, lambda=0.16, quadrature=1.58148, prediction=1.5821, relative_error=0.00039522
    - radius=0.05, lambda=2.89, quadrature=1.68272, prediction=1.68335, relative_error=0.000371463
    - radius=0.15, lambda=0.16, quadrature=0.522719, prediction=0.522823, relative_error=0.0002
    - radius=0.15, lambda=2.89, quadrature=0.608438, prediction=0.608542, relative_error=0.000171836
    - radius=0.35, lambda=0.16, quadrature=0.223001, prediction=0.223092, relative_error=0.000408513
    - radius=0.35, lambda=2.89, quadrature=0.275583, prediction=0.275674, relative_error=0.000330666
- **worst_cross_error**: 1.00992e-10
- **worst_self_error**: 0.000408513
- **the_cross_identity_is_exact**: True
- **the_self_identity_is_grid_limited**: True
- **what_it_justifies**: the smeared mouth's ambient block, and with it every result in this module

## T3_the_negative_mode_does_not_survive — PASS

- **samples**: 3078
- **positives**: 0
- **worst_approach_to_zero**: -0.000510941
- **the_tube_side_is_always_negative**: True
- **the_ambient_side_is_always_positive**: True
- **there_is_no_growing_mode**: True
- **the_answer**: no: a finite-radius mouth removes the negative mode, and the sign structure says no parameter choice can bring it back

## T4_the_instability_was_the_linearization — PASS

- **rows**:
    - radius=0.02, linearized_root=50.02, root_times_radius=1.0004, exact_has_a_root=False, agreement=[kappa_a=0.02, exact=-4.11117, linearized=-4.1096, relative_gap=0.000383437, kappa_a=0.1, exact=-3.62312, linearized=-3.59738, relative_gap=0.00710222, kappa_a=1, exact=-1.72202, linearized=-0.00159155, relative_gap=0.999076, kappa_a=3, exact=-0.662121, linearized=7.95722, relative_gap=13.0178]
    - radius=0.05, linearized_root=20.0499, root_times_radius=1.00249, exact_has_a_root=False, agreement=[kappa_a=0.02, exact=-2.7142, linearized=-2.71244, relative_gap=0.000646737, kappa_a=0.1, exact=-1.5054, linearized=-1.49407, relative_gap=0.00752532, kappa_a=1, exact=-0.692631, linearized=-0.00397887, relative_gap=0.994255, kappa_a=3, exact=-0.266148, linearized=3.18177, relative_gap=12.9549]
    - radius=0.15, linearized_root=6.8142, root_times_radius=1.02213, exact_has_a_root=False, agreement=[kappa_a=0.02, exact=-10.5146, linearized=-10.5103, relative_gap=0.000408778, kappa_a=0.1, exact=-0.925131, linearized=-0.917809, relative_gap=0.00791525, kappa_a=1, exact=-0.243103, linearized=-0.0120102, relative_gap=0.950596, kappa_a=3, exact=-0.0928436, linearized=1.05705, relative_gap=12.3853]
    - radius=0.35, linearized_root=3.21976, root_times_radius=1.12692, exact_has_a_root=False, agreement=[kappa_a=0.02, exact=-54.4294, linearized=-54.4188, relative_gap=0.00019346, kappa_a=0.1, exact=-2.4299, linearized=-2.41823, relative_gap=0.00480374, kappa_a=1, exact=-0.137769, linearized=-0.0344745, relative_gap=0.749766, kappa_a=3, exact=-0.0486884, linearized=0.445435, relative_gap=10.1487]
- **worst_gap_below_kappa_a_of_0p1**: 0.00791525
- **worst_gap_at_kappa_a_of_3**: 10.1487
- **linearized_roots_times_radius**:
    - 1.0004
    - 1.00249
    - 1.02213
    - 1.12692
- **every_linearized_model_has_a_root**: True
- **no_exact_model_has_one**: True
- **the_root_sits_at_kappa_a_of_order_one**: True
- **the_diagnosis**: the linearization is excellent for κa ≪ 1 and wrong in SIGN beyond it; the mode it produces sits at κa ≈ 1, the edge of its own validity

## T5_the_mode_became_soft_and_positive — PASS

- **rows**:
    - radius=0.4, lambda_0=0.496462, omega_0=0.704601, closed_form=0.888889, ratio=0.55852, antisymmetric=None
    - radius=0.2, lambda_0=0.351588, omega_0=0.592949, closed_form=0.444444, ratio=0.791074, antisymmetric=None
    - radius=0.1, lambda_0=0.204256, omega_0=0.451947, closed_form=0.222222, ratio=0.919154, antisymmetric=None
    - radius=0.05, lambda_0=0.107591, omega_0=0.328011, closed_form=0.111111, ratio=0.968321, antisymmetric=None
    - radius=0.02, lambda_0=0.0439836, omega_0=0.209723, closed_form=0.0444444, ratio=0.98963, antisymmetric=None
    - radius=0.01, lambda_0=0.0221157, omega_0=0.148713, closed_form=0.0222222, ratio=0.995207, antisymmetric=None
    - radius=0.005, lambda_0=0.0110856, omega_0=0.105288, closed_form=0.0111111, ratio=0.997701, antisymmetric=None
- **every_one_is_positive**: True
- **every_one_is_below_the_gap**: True
- **no_antisymmetric_bound_state**: True
- **worst_closed_form_error**: 0.0103697
- **the_capacitance_formula_holds**: True
- **the_point_limit_approaches_zero_from_above**: True
- **what_it_means**: the mode is soft, not unstable; PR #260 put it on the wrong side of zero

## T6_the_delay_survives — PASS

- **rows**:
    - length=0.4, onset_same_mouth=0, onset_opposite=0.208282, prediction=0.36
    - length=0.6, onset_same_mouth=0, onset_opposite=0.409698, prediction=0.56
    - length=0.9, onset_same_mouth=0, onset_opposite=0.709534, prediction=0.86
    - length=1.2, onset_same_mouth=0, onset_opposite=1.00937, prediction=1.16
    - length=2, onset_same_mouth=0, onset_opposite=1.06659, prediction=1.26
    - length=3, onset_same_mouth=0, onset_opposite=1.06659, prediction=1.26
- **by_radius**:
    - radius=0.01, onset=0.7164
    - radius=0.02, onset=0.709534
    - radius=0.03, onset=0.704956
    - radius=0.04, onset=0.704956
- **slope_in_length**: 1.00101
- **slope_in_radius**: -0.389099
- **onset_spread_above_d**: 0
- **contour**: 0.4
- **the_traversal_time_survives**: True
- **the_ambient_path_still_takes_over**: True
- **the_mouth_shift_is_subleading**: True
- **the_correction**: the mouth shifts the onset earlier by O(a), well below the leading min(L,d); a first draft predicted −2a from an ambient block missing the shell form factor

## T7_the_static_results_survive — PASS

- **rows**:
    - lambda=1e-08, det_S=-3.63685e-08, det_S_over_lambda=-3.63685, antisymmetry=1.75853e-07, defect=-8.84194e+06, minus_beta=-8.84194e+06, defect_error=8.30423e-13
    - lambda=1e-06, det_S=-3.63688e-06, det_S_over_lambda=-3.63688, antisymmetry=1.75856e-05, defect=-88419.4, minus_beta=-88419.4, defect_error=3.55506e-12
    - lambda=0.0001, det_S=-0.000364019, det_S_over_lambda=-3.64019, antisymmetry=0.0017617, defect=-884.206, minus_beta=-884.206, defect_error=2.76436e-14
- **linear_coefficient**: -3.63797
- **det_S_is_linear_in_lambda**: True
- **worst_antisymmetry**: 1.75856e-05
- **the_static_response_is_still_rank_one**: True
- **worst_defect_error**: 3.55506e-12
- **the_defect_is_still_minus_beta**: True
- **why**: both came from the tube's zero mode, which the mouth does not touch

## T8_monopole_matching_is_what_is_left — PASS

- **lambda**: 0.16
- **separation**: 1.3
- **rows**:
    - radius=0.02, monopole=0.0583406, dipole=0.000838565, dipole_over_monopole=0.0143736, ratio_over_a_over_d=0.934285, dropped_power_fraction=6.88716e-05
    - radius=0.05, monopole=0.0583577, dipole=0.00209737, dipole_over_monopole=0.0359399, ratio_over_a_over_d=0.934437, dropped_power_fraction=0.000430662
    - radius=0.15, monopole=0.0585215, dipole=0.00631955, dipole_over_monopole=0.107987, ratio_over_a_over_d=0.935885, dropped_power_fraction=0.00389543
    - radius=0.35, monopole=0.0593513, dipole=0.0150718, dipole_over_monopole=0.253942, ratio_over_a_over_d=0.943213, dropped_power_fraction=0.0217584
- **dipole_scales_like_a_over_d**: True
- **worst_dropped_fraction**: 0.0217584
- **smallest_dropped_fraction**: 6.88716e-05
- **the_scope**: one channel per mouth couples only ℓ = 0; the dropped multipoles are screened as (a/d)^ℓ, PR #250's law, so the leading omission is the dipole

## T9_assessment — PASS

- **n_passed**: 8
- **n_total**: 8

## verdict — THE_NEGATIVE_MODE_DOES_NOT_SURVIVE_A_FINITE_MOUTH

*** THE ANSWER IS NO. *** PR #260 gated the roadmap on whether its growing mode survives a finite-radius mouth. It does not, and the statement is STRUCTURAL rather than parametric. Smearing the coupling over a sphere of radius a -- the same operator on both sides, so the composite stays manifestly self-adjoint -- replaces the renormalized self-energy g(lam), which is negative, by the UNSUBTRACTED Green function G(a,lam), which is positive. The mean-value identities that give it, G_self = f(a)G(a) and G_cross = f(a)^2 G(d), are checked against direct quadrature on S^3 to 1.0e-10 for the cross term and 4.1e-04 for the self term, the latter grid-limited by the integrable singularity at coincidence. THEN THE SIGNS DECIDE IT: at lam < 0 the tube's channel functions are strictly NEGATIVE -- a passive interior supplies no restoring force -- and the ambient's are strictly POSITIVE, because f and G are positive on the imaginary axis and G(a) > f(a)G(d) whenever a < d/2, which disjoint mouths require anyway. A difference of a negative and a positive number has no zero, so no parameter choice can produce a growing mode: 3078 samples over (a, d, L, A, m, sigma) give 0 positives, worst approach -5.1e-04. *** AND PR #260's MODE WAS THE LINEARIZATION. *** That round modelled the mouth as a CONSTANT shift 1/(4 pi a) -- the leading term of G(a,lam) = 1/(4 pi a) + g(lam) + O(a). The exact G(a,-k^2) is SCREENED, ~ e^{-k a}/(4 pi a), and dies; the constant does not, so it eventually beats the tube's -1/(A k) and crosses. The two agree to 0.8% for k a <= 0.1 and differ by 1015% at k a = 3, and the linearized root sits at k a = 1.000, 1.002, 1.022, 1.127 -- AT THE EDGE OF ITS OWN VALIDITY. A mode that lives at the scale where its derivation fails is an artifact, and that is now shown rather than suspected. WHERE THE MODE WENT: SOFT AND POSITIVE. The composite has exactly one state below the free gap lam = 1, in the symmetric channel, and as a -> 0 it goes to 8 pi a/(A L) -- two mouth capacitances 4 pi a restoring a tube of volume A L -- matched to 1.0%. So the point limit drives the mode to zero FROM ABOVE. PR #260 did not get a rate slightly wrong; it took a mode approaching 0+ and put it on the other side of zero, at lam ~ -1/a^2. THE GOOD RESULTS SURVIVE: the traversal delay has slope 1.0010 in L against a predicted 1, saturating at the ambient path to 0.0e+00, with the mouth adding only a sub-leading O(a) shift (slope -0.39; a first draft predicted -2a from an ambient block missing the shell form factor, and that prediction is recorded as wrong). The static response is still rank one, det S linear in lam, and PR #258's W = -beta(lam) holds to 3.6e-12 -- all three came from the TUBE's zero mode, which the mouth does not touch. The contour is easier too: with nothing to clear, eps = 0.4 suffices where PR #260 needed eps > 2. WHAT IS STILL PUT IN: one channel per mouth, so only l = 0 couples; the dropped multipoles obey PR #250's screening law, dipole/monopole = 0.934 (a/d) across a decade in a, and the dropped power fraction is 6.9e-05 at a = 0.02. THE MOUTHS ARE SPHERES IN A FIXED AMBIENT, NOT A SOLVED NECK, AND THERE IS NO BACKREACTION. But the gate PR #260 set is answered: with the mouth resolved there is no growing mode, so the stationary-action and backreaction rounds can proceed on this background and measure the physics rather than an artifact.
