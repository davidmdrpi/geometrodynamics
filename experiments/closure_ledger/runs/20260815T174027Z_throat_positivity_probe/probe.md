# probe — throat_positivity

_2026-08-15T17:40:27.230390+00:00_

## T1_goal — PASS

- **question**: for which Hermitian A is the self-adjoint point-throat operator non-negative? PR #256 showed flux conservation does not give stability and mapped a two-parameter slice by scanning; what is the answer on all four parameters?
- **scope**: a linear scalar field on a fixed Einstein static universe. The throat is point-supported: no interior, no proper length, no delay. The boundary data is four real numbers chosen, not derived; PR #249 is what would fix them. No backreaction, no stress tensor, no topology change, no rate, no two-source invariant.

## T2_the_positive_sector_is_a_shifted_psd_cone — PASS

- **n_draws**: 200
- **n_mismatches**: 0
- **the_criterion_is_exact**: True
- **n_stable**: 19
- **both_verdicts_occur**: True
- **n_with_complex_beta**: 200
- **the_light_cone_form_agrees**: True
- **rows**:
    - alpha1=0.0417328, alpha2=-0.047118, beta=0.19251+0.211478j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.0226377, norm_x=0.259732, inside_the_cone=False
    - alpha1=-0.0571451, alpha2=0.228211, beta=0.0404291-0.0721901j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.110863, norm_x=0.1601, inside_the_cone=False
    - alpha1=-0.120968, alpha2=-0.147441, beta=0.0476877+0.062449j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=-0.108874, norm_x=0.0638404, inside_the_cone=False
    - alpha1=0.200835, alpha2=0.1804, beta=0.0764064-0.124528j, predicted_non_negative=True, scan_says_stable=True, agrees=True, x0=0.215948, norm_x=0.128044, inside_the_cone=True
    - alpha1=0.175882, alpha2=-0.110468, beta=-0.0994049-0.0711357j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.0580374, norm_x=0.217737, inside_the_cone=False
    - alpha1=-0.0775596, alpha2=0.0860073, beta=-0.120761-0.0406378j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.0295541, norm_x=0.192249, inside_the_cone=False
    - alpha1=-0.0113279, alpha2=0.0472096, beta=-0.127049+0.203885j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.0432712, norm_x=0.270578, inside_the_cone=False
    - alpha1=-0.179791, alpha2=-0.0875452, beta=0.0470248+0.0845535j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=-0.108338, norm_x=0.0963252, inside_the_cone=False
    - alpha1=0.00580716, alpha2=-0.110248, beta=-0.30185-0.302138j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=-0.0268903, norm_x=0.466196, inside_the_cone=False
    - alpha1=0.0417614, alpha2=-0.127233, beta=0.072181+0.136931j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=-0.0174057, norm_x=0.16265, inside_the_cone=False
    - alpha1=-0.124597, alpha2=0.12079, beta=-0.155596+0.00611296j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.0234271, norm_x=0.238139, inside_the_cone=False
    - alpha1=0.00624145, alpha2=0.0024526, beta=-0.229326+0.07671j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.0296773, norm_x=0.288143, inside_the_cone=False
    - alpha1=0.20415, alpha2=-0.0767168, beta=0.0839499+0.343871j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.089047, norm_x=0.373137, inside_the_cone=False
    - alpha1=-0.0397358, alpha2=0.30337, beta=0.0122995-0.307378j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.157147, norm_x=0.353858, inside_the_cone=False
    - alpha1=0.229564, alpha2=0.100382, beta=0.136767-0.252937j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.190303, norm_x=0.275601, inside_the_cone=False
    - alpha1=-0.010751, alpha2=-0.0524034, beta=0.191724+0.0800495j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=-0.00624686, norm_x=0.165469, inside_the_cone=False
    - alpha1=-0.153364, alpha2=0.0803491, beta=0.0248667-0.0597569j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=-0.0111771, norm_x=0.133344, inside_the_cone=False
    - alpha1=-0.0224785, alpha2=-0.133783, beta=0.147176+0.00614801j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=-0.0528004, norm_x=0.11353, inside_the_cone=False
    - alpha1=0.170968, alpha2=-0.0100906, beta=0.0698971-0.0548269j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.105769, norm_x=0.107996, inside_the_cone=False
    - alpha1=-0.117027, alpha2=0.112191, beta=0.126235+0.211448j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.0229124, norm_x=0.252788, inside_the_cone=False
    - alpha1=-0.156436, alpha2=-0.0486781, beta=-0.127865+0.0956607j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=-0.0772268, norm_x=0.207672, inside_the_cone=False
    - alpha1=-0.153976, alpha2=0.100251, beta=-0.0893099+0.142505j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=-0.00153243, norm_x=0.235442, inside_the_cone=False
    - alpha1=0.0821237, alpha2=0.000463962, beta=-0.0702033+0.190705j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.0666241, norm_x=0.228266, inside_the_cone=False
    - alpha1=-0.0810118, alpha2=0.0932199, beta=0.190234+0.0634504j, predicted_non_negative=False, scan_says_stable=False, agrees=True, x0=0.0314344, norm_x=0.178125, inside_the_cone=False
- **the_criterion**: A ⪰ Γ(0) in the Löwner order
- **the_geometry**: A − Γ(0) = x₀I + x·σ is PSD iff x₀ ≥ |x| — the stable set is a forward light cone in the four dimensions of Hermitian boundary data
- **why**: dΓ/dλ ≻ 0 below threshold, so every eigenvalue of A − Γ(λ) is decreasing in λ and both run to +∞ as λ → −∞; one crosses zero below λ* iff it is already negative at λ*

## T3_the_inertia_theorem — PASS

- **rows**:
    - lambda_star=-2, n_draws=40, mismatches=0
    - lambda_star=0, n_draws=40, mismatches=0
    - lambda_star=0.5, n_draws=40, mismatches=0
    - lambda_star=0.9, n_draws=40, mismatches=0
- **n_tested**: 160
- **n_mismatches**: 0
- **the_inertia_theorem_holds**: True
- **gamma_derivative**:
    - lmbda=-100, eigenvalues=[0.00397886, 0.00397889]
    - lmbda=-9, eigenvalues=[0.0129007, 0.0136251]
    - lmbda=-1, eigenvalues=[0.0256265, 0.0523743]
    - lmbda=-0.05, eigenvalues=[0.0307173, 0.125708]
    - lmbda=-0.001, eigenvalues=[0.031056, 0.135392]
- **d_gamma_d_lambda_is_positive_definite**: True
- **the_theorem**: #{eigenvalues < λ*} = #{negative eigenvalues of A − Γ(λ*)} for every λ* below the free ground state
- **stability_is_the_lambda_star_zero_case**: True

## T4_the_boundary_is_a_zero_mode — PASS

- **rows**:
    - direction=[1, 0, 0], x0=0.02, lightlike_defect=-3.46945e-18, smallest_eigenvalue=-3.46945e-18, is_a_zero_mode=True, secular_at_zero=-1.38778e-19, marginal_lambda_found_by_root_finding=-1.43051e-14, n_strictly_growing=0, q=[-0.70711+0.j  0.70711+0.j]
    - direction=[1, 0, 0], x0=0.1, lightlike_defect=0, smallest_eigenvalue=0, is_a_zero_mode=True, secular_at_zero=0, marginal_lambda_found_by_root_finding=-1.43051e-14, n_strictly_growing=0, q=[-0.70711+0.j  0.70711+0.j]
    - direction=[1, 0, 0], x0=0.4, lightlike_defect=0, smallest_eigenvalue=0, is_a_zero_mode=True, secular_at_zero=0, marginal_lambda_found_by_root_finding=-1.43051e-14, n_strictly_growing=0, q=[-0.70711+0.j  0.70711+0.j]
    - direction=[0, 1, 0], x0=0.02, lightlike_defect=0, smallest_eigenvalue=0, is_a_zero_mode=True, secular_at_zero=0, marginal_lambda_found_by_root_finding=-1.43051e-14, n_strictly_growing=0, q=[-0.70711+0.j       0.     +0.70711j]
    - direction=[0, 1, 0], x0=0.1, lightlike_defect=0, smallest_eigenvalue=0, is_a_zero_mode=True, secular_at_zero=0, marginal_lambda_found_by_root_finding=-1.43051e-14, n_strictly_growing=0, q=[-0.70711+0.j       0.     +0.70711j]
    - direction=[0, 1, 0], x0=0.4, lightlike_defect=0, smallest_eigenvalue=0, is_a_zero_mode=True, secular_at_zero=0, marginal_lambda_found_by_root_finding=-1.43051e-14, n_strictly_growing=0, q=[-0.70711+0.j       0.     +0.70711j]
    - direction=[0, 0, 1], x0=0.02, lightlike_defect=0, smallest_eigenvalue=0, is_a_zero_mode=True, secular_at_zero=0, marginal_lambda_found_by_root_finding=-1.43051e-14, n_strictly_growing=0, q=[0.+0.j 1.+0.j]
    - direction=[0, 0, 1], x0=0.1, lightlike_defect=0, smallest_eigenvalue=0, is_a_zero_mode=True, secular_at_zero=0, marginal_lambda_found_by_root_finding=-1.43051e-14, n_strictly_growing=0, q=[0.+0.j 1.+0.j]
    - direction=[0, 0, 1], x0=0.4, lightlike_defect=0, smallest_eigenvalue=0, is_a_zero_mode=True, secular_at_zero=0, marginal_lambda_found_by_root_finding=-1.43051e-14, n_strictly_growing=0, q=[0.+0.j 1.+0.j]
    - direction=[0.6, -0.5, 0.62], x0=0.02, lightlike_defect=3.46945e-18, smallest_eigenvalue=2.60209e-18, is_a_zero_mode=True, secular_at_zero=1.12531e-19, marginal_lambda_found_by_root_finding=0, n_strictly_growing=0, q=[-0.43489+0.j       0.69177-0.57648j]
    - direction=[0.6, -0.5, 0.62], x0=0.1, lightlike_defect=0, smallest_eigenvalue=0, is_a_zero_mode=True, secular_at_zero=1.12531e-18, marginal_lambda_found_by_root_finding=0, n_strictly_growing=0, q=[-0.43489+0.j       0.69177-0.57648j]
    - direction=[0.6, -0.5, 0.62], x0=0.4, lightlike_defect=0, smallest_eigenvalue=0, is_a_zero_mode=True, secular_at_zero=1.8005e-17, marginal_lambda_found_by_root_finding=0, n_strictly_growing=0, q=[-0.43489+0.j       0.69177-0.57648j]
- **every_boundary_point_has_a_zero_mode**: True
- **worst_secular_at_zero**: 1.8005e-17
- **the_secular_function_vanishes_there**: True
- **the_marginal_mode_sits_at_lambda_zero**: True
- **worst_marginal_lambda**: 1.43051e-14
- **no_boundary_point_is_strictly_unstable**: True
- **apex_degeneracy**: 2
- **the_apex_carries_two_zero_modes**: True
- **what_this_shows**: the cone's boundary is where λ = 0 enters the spectrum — located independently by root-finding, not read off the eigenvalue — so it is detectable rather than conventional

## T5_the_growth_rate_turns_on_with_a_square_root — PASS

- **rows**:
    - epsilon=0.01, lmbda=-0.078196, sigma=0.279635, lambda_over_epsilon=-7.8196, sigma_over_sqrt_epsilon=2.79635
    - epsilon=0.001, lmbda=-0.00741719, sigma=0.0861231, lambda_over_epsilon=-7.41719, sigma_over_sqrt_epsilon=2.72345
    - epsilon=0.0001, lmbda=-0.000737869, sigma=0.0271637, lambda_over_epsilon=-7.37869, sigma_over_sqrt_epsilon=2.71637
    - epsilon=1e-05, lmbda=-7.37486e-05, sigma=0.0085877, lambda_over_epsilon=-7.37486, sigma_over_sqrt_epsilon=2.71567
    - epsilon=1e-06, lmbda=-7.37448e-06, sigma=0.0027156, lambda_over_epsilon=-7.37448, sigma_over_sqrt_epsilon=2.7156
- **exponents**:
    - 0.511472
    - 0.50113
    - 0.500113
    - 0.500011
- **asymptotic_exponent**: 0.500011
- **lambda_over_epsilon_limit**: -7.37448
- **predicted_from_the_eigenvalue_slope**: -7.37443
- **mu_prime**: 0.135604
- **exponent_is_one_half**: True
- **lambda_is_linear_in_epsilon**: True
- **the_coefficient_matches_the_eigenvalue_slope**: True
- **what_this_shows**: the boundary is a genuine threshold, not a numerical artefact: the growth rate is continuous and rises like √ε

## T6_the_wedge_is_a_slice — PASS

- **n_slice_points**: 143
- **slice_mismatches**: 0
- **the_wedge_is_exactly_the_slice**: True
- **worst_off_slice_coordinate**: 0
- **the_slice_really_is_x2_equals_x3_equals_zero**: True
- **n_general_draws**: 400
- **n_the_wedge_rule_gets_wrong**: 65
- **the_slice_rule_is_not_enough_in_general**: True
- **why**: Im β and the mouth asymmetry are the two dimensions the wedge does not see; both push A out of the cone without changing α ± Re β

## T7_where_the_apex_sits — PASS

- **rows**:
    - separation=0.2, trace=-0.0506606, zero_A_margin=-0.349722, determinant=-0.140023, eigenvalues=[-0.400383, 0.349722], indefinite=True, negative_semidefinite=False, zero_matrix_is_stable=False, channel_thresholds=[-0.400383, 0.349722]
    - separation=0.5, trace=-0.0506606, zero_A_margin=-0.114237, determinant=-0.0188375, eigenvalues=[-0.164898, 0.114237], indefinite=True, negative_semidefinite=False, zero_matrix_is_stable=False, channel_thresholds=[-0.164898, 0.114237]
    - separation=0.8, trace=-0.0506606, zero_A_margin=-0.0573528, determinant=-0.00619487, eigenvalues=[-0.108013, 0.0573528], indefinite=True, negative_semidefinite=False, zero_matrix_is_stable=False, channel_thresholds=[-0.108013, 0.0573528]
    - separation=1.3, trace=-0.0506606, zero_A_margin=-0.023082, determinant=-0.00170213, eigenvalues=[-0.0737426, 0.023082], indefinite=True, negative_semidefinite=False, zero_matrix_is_stable=False, channel_thresholds=[-0.0737426, 0.023082]
    - separation=2, trace=-0.0506606, zero_A_margin=-0.00647105, determinant=-0.000369702, eigenvalues=[-0.0571316, 0.00647105], indefinite=True, negative_semidefinite=False, zero_matrix_is_stable=False, channel_thresholds=[-0.0571316, 0.00647105]
    - separation=2.8, trace=-0.0506606, zero_A_margin=-0.000499403, determinant=-2.55494e-05, eigenvalues=[-0.05116, 0.000499403], indefinite=True, negative_semidefinite=False, zero_matrix_is_stable=False, channel_thresholds=[-0.05116, 0.000499403]
    - separation=3, trace=-0.0506606, zero_A_margin=-8.48374e-05, determinant=-4.30511e-06, eigenvalues=[-0.0507454, 8.48374e-05], indefinite=True, negative_semidefinite=False, zero_matrix_is_stable=False, channel_thresholds=[-0.0507454, 8.48374e-05]
    - separation=3.14159, trace=-0.0506606, zero_A_margin=0, determinant=0, eigenvalues=[-0.0506606, 0], indefinite=False, negative_semidefinite=True, zero_matrix_is_stable=True, channel_thresholds=[-0.0506606, 0]
- **trace_is_separation_independent**: True
- **trace_value**: -0.0506606
- **predicted_trace**: -0.0506606
- **trace_matches_minus_one_over_two_pi_squared**: True
- **the_apex_is_indefinite_away_from_the_antipode**: True
- **the_apex_is_negative_semidefinite_at_the_antipode**: True
- **the_zero_matrix_is_unstable_away_from_the_antipode**: True
- **antipodal**: Gamma_0=[[-0.02533  0.02533]
 [ 0.02533 -0.02533]], eigenvalues=[-0.0506606, 0], trace=-0.0506606, determinant=0, trace_is_two_g0=-0.0506606, channel_thresholds=[-0.0506606, 0], indefinite=False, negative_semidefinite=True, zero_matrix_is_stable=True, zero_A_margin=0, G_pi_at_zero=0.0253303, minus_g0=0.0253303, marginal_channel=symmetric
- **the_antipodal_endpoint_is_marginal**: True
- **at_the_antipode_A_zero_sits_on_the_boundary**: True
- **the_marginal_channel_is_symmetric**: symmetric
- **eigenvalues_are_the_channel_thresholds**: True
- **the_symmetric_threshold_closes_at_the_antipode**: True

## T8_how_big_is_it — PASS

- **half_width**: 0.2
- **n_draws**: 4000
- **n_stable**: 332
- **fraction**: 0.083
- **the_box**: |α_j| ≤ w and |Re β|, |Im β| ≤ w
- **caveat**: a cone is unbounded; the fraction is box-dependent

## T9_the_monotonicity_is_a_gram_matrix — PASS

- **rows**:
    - separation=1.3, lmbda=-9, abs_error=6.70833e-13, eigenvalues=[0.0129007, 0.0136251], positive_definite=True
    - separation=1.3, lmbda=-1, abs_error=3.66391e-12, eigenvalues=[0.0256265, 0.0523743], positive_definite=True
    - separation=1.3, lmbda=-0.05, abs_error=3.56e-12, eigenvalues=[0.0307173, 0.125708], positive_definite=True
    - separation=1.3, lmbda=0, abs_error=2.86288e-12, eigenvalues=[0.031063, 0.135604], positive_definite=True
    - separation=1.3, lmbda=0.5, abs_error=8.11287e-12, eigenvalues=[0.0351769, 0.445356], positive_definite=True
    - separation=2.6, lmbda=-9, abs_error=6.70833e-13, eigenvalues=[0.013237, 0.0132888], positive_definite=True
    - separation=2.6, lmbda=-1, abs_error=2.63625e-12, eigenvalues=[0.0311837, 0.0468171], positive_definite=True
    - separation=2.6, lmbda=-0.05, abs_error=3.56e-12, eigenvalues=[0.0402153, 0.11621], positive_definite=True
    - separation=2.6, lmbda=0, abs_error=2.93818e-12, eigenvalues=[0.0408588, 0.125808], positive_definite=True
    - separation=2.6, lmbda=0.5, abs_error=6.81355e-12, eigenvalues=[0.0487623, 0.43177], positive_definite=True
    - separation=3.14159, lmbda=-9, abs_error=6.70833e-13, eigenvalues=[0.0132449, 0.0132809], positive_definite=True
    - separation=3.14159, lmbda=-1, abs_error=2.63625e-12, eigenvalues=[0.0315815, 0.0464193], positive_definite=True
    - separation=3.14159, lmbda=-0.05, abs_error=3.56e-12, eigenvalues=[0.0409932, 0.115432], positive_definite=True
    - separation=3.14159, lmbda=0, abs_error=2.98788e-12, eigenvalues=[0.0416667, 0.125], positive_definite=True
    - separation=3.14159, lmbda=0.5, abs_error=6.81355e-12, eigenvalues=[0.0499636, 0.430569], positive_definite=True
- **worst_abs_error**: 8.11287e-12
- **the_gram_sum_is_the_closed_form_derivative**: True
- **positive_definite_everywhere**: True
- **including_at_the_antipode**: True
- **the_identity**: dΓ_ij/dλ = ⟨δ_i, (H₀−λ)^{-2} δ_j⟩ — a Gram matrix
- **what_this_upgrades**: Löwner monotonicity of Γ from a fact sampled at a few λ to an analytic consequence; the root scans elsewhere are regression checks

## T10_the_criterion_extends_beyond_the_chart — PASS

- **chart_draws**: 60
- **chart_mismatches**: 0
- **the_general_form_agrees_with_the_cone_on_the_chart**: True
- **stratum_draws**: 60
- **stratum_mismatches**: 0
- **the_general_form_agrees_with_the_scan_on_the_strata**: True
- **every_stratum_has_one_dirichlet_direction**: True
- **worst_hermitian_defect**: 1.76736e-12
- **worst_row_space_defect**: 9.60933e-16
- **the_reduction_is_legitimate**: True
- **free_stratum**: k=0, non_negative=True, stratum=free (q = 0): no mouth-active spectrum
- **the_free_stratum_is_non_negative**: True
- **the_scope**: A ⪰ Γ(0) is the criterion in the finite-A chart; the general one is A_eff ⪰ P†Γ(0)P on the allowed-charge subspace, and the chart is its k = 2 case

## T11_assessment — PASS

- **n_passed**: 10
- **n_total**: 10

## verdict — THE_POSITIVE_SECTOR_IS_A_LIGHT_CONE_AT_GAMMA_ZERO

NON-NEGATIVE IF AND ONLY IF A >= Gamma(0), IN THE LOEWNER ORDER -- in the finite-A chart, with the general form on the whole U(2) family got by restricting to the allowed-charge subspace. PR #256 left the question open on three of its four parameters and answered the fourth by scanning; the general answer is a single inequality, and the reason is one line. Gamma(lambda) has dGamma/dlambda POSITIVE DEFINITE below threshold, so every eigenvalue of M(lambda) = A - Gamma(lambda) is strictly decreasing in lambda, while as lambda -> -infinity Gamma -> -(sigma/4pi) I and both eigenvalues run to +infinity. An eigenvalue therefore crosses zero somewhere below threshold if and only if it is already negative AT threshold. Checked against an actual negative-lambda root scan on 200 random Hermitian A -- every one with complex beta and unequal mouths, so all four parameters are exercised -- with 0 mismatches and 19 of them stable, so both verdicts occur. THE GEOMETRY IS A LIGHT CONE: Hermitian 2x2 matrices are R^4 under A - Gamma(0) = x0 I + x.sigma, and positive semidefiniteness is x0 >= |x|, so the stable set is the forward light cone with apex at A = Gamma(0) -- convex, closed under positive scaling from the apex, and four-dimensional. AND THE SAME ARGUMENT COUNTS: #{mouth-active eigenvalues < lambda*} = #{negative eigenvalues of A - Gamma(lambda*)} for any lambda* below the free ground state, 0 mismatches in 160 tests at lambda* = -2, 0, 0.5 and 0.9. That is a Krein-type inertia theorem, and stability is its lambda* = 0 case. THE BOUNDARY IS DETECTABLE, NOT CONVENTIONAL: on the null surface A - Gamma(0) is rank one, so lambda = 0 enters the spectrum as a genuine ZERO MODE -- a static solution supported by the throat, below the free ground state -- located independently by root-finding at 1.4e-14 with the secular function vanishing to 1.8e-17. At the apex there are TWO. AND THE INSTABILITY OUTSIDE TURNS ON LIKE A SQUARE ROOT: lambda is linear in the distance eps past the boundary (-7.374476, against -7.374433 predicted from the eigenvalue slope rather than fitted), so the growth rate rises with exponent 0.50001. PR #256's WEDGE IS THE x2 = x3 = 0 SLICE, reproduced exactly at all 143 sampled points -- and applied to general boundary data by averaging the mouths and dropping Im beta it gets 65 of 400 draws wrong, which is why the general form was needed rather than a wider scan. FINALLY, WHERE THE APEX SITS: tr Gamma(0) = 2 g0 = -1/(2 pi^2) at EVERY mouth separation, its eigenvalues are exactly PR #256's two channel thresholds, and and its determinant g0^2 - G0^2 is negative for 0 < d < pi -- so Gamma(0) is indefinite there and A = 0 is unstable. THE EXACT ANTIPODE IS A DIFFERENT STATEMENT, and for this geometry it is the load-bearing one: G_d has a REMOVABLE singularity at d = pi, with G_pi(0) = +1/(4 pi^2) = -g0, so Gamma(0) = g0[[1,-1],[-1,1]] has eigenvalues (2 g0, 0) -- negative SEMIdefinite, not indefinite -- and A = 0 is MARGINALLY non-negative there, sitting on the cone's boundary with a zero mode in the symmetric channel. TWO FURTHER TIGHTENINGS. First, the monotonicity everything rests on is not a sampled fact: Gamma_ij(lambda) = <delta_i,(H0-lambda)^-1 delta_j> up to a lambda-independent subtraction, so dGamma/dlambda is the GRAM MATRIX of (H0-lambda)^-1 delta_j -- PSD for free, positive definite for distinct mouths -- computed independently from the S3 addition theorem and agreeing with the closed form to 8.1e-12, antipode included. Second, A >= Gamma(0) is the criterion IN A CHART: phi_reg = A q needs B invertible, and the strata it misses are Dirichlet directions, reached as ||A|| -> infinity and not represented by any finite Hermitian A. The general criterion is A_eff >= P^dag Gamma(0) P on the allowed-charge subspace, agreeing with the cone on 60 chart draws and with a root scan on 60 stratum draws, with the reduction's own assumptions checked rather than assumed. HOW BIG THE REGION IS depends on the box, and the box is stated: 0.083 of a uniform draw over |alpha_j|, |Re beta|, |Im beta| <= 0.2. WHAT IS STILL PUT IN: the boundary data itself, four real numbers chosen and not derived, with PR #249 still the thing that would fix them from matter. The throat is point-supported -- no interior, no proper length, no delay. No backreaction, no stress tensor, no topology change, no rate, and no two-source invariant; what this round buys the next one is a stated region to work inside and a count of what goes wrong outside it.
