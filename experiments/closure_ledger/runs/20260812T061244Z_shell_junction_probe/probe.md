# shell_junction probe

_generated 2026-08-12T06:12:44.229708+00:00_

## T1_goal — PASS
- **question**: does a detached oppositely-glued bulk shell change the throat's effective force and stability while itself remaining non-exotic?
- **observables**: ["1 shell's own Israel surface stress sigma", '2 shell-induced force on the throat, -dV_shell/db', '3 stiffness d2V/db2 and the coupled normal modes']
- **why_three**: a positive stiffness alone would mean restoring and would not establish that the shell supplied it

## T2_the_machinery_reproduces_published_shells — PASS
- **rows**:
    - case=bubble, mass=0.001, sigma=7.95815e-07, rest_mass=0.00100005, rest_mass_error=5.0005e-05
    - case=bubble, mass=0.01, sigma=7.96173e-06, rest_mass=0.010005, rest_mass_error=0.000500501
    - case=bubble, mass=0.1, sigma=7.99794e-05, rest_mass=0.100505, rest_mass_error=0.00505063
    - case=visser_throat, mass=0.1, radius=3, sigma=-0.0512528, reference=-0.0512528
    - case=visser_throat, mass=0.1, radius=5, sigma=-0.0311879, reference=-0.0311879
    - case=visser_throat, mass=0.5, radius=3, sigma=-0.0433165, reference=-0.0433165
    - case=visser_throat, mass=0.5, radius=5, sigma=-0.0284705, reference=-0.0284705
- **worst_visser_error**: 6.93889e-18
- **the_bubble_is_ordinary**: True
- **its_rest_mass_is_the_bulk_mass**: True
- **visser_is_reproduced**: True

## T3_every_minimal_surface_is_exotic — PASS
- **samples**: 200000
- **anti_aligned_positive_sigma**: 0
- **worst_anti_aligned_sigma**: -0.000670072
- **aligned_positive_sigma**: 100170
- **aligned_positive_fraction**: 0.50085
- **every_minimal_surface_is_exotic**: True
- **the_sweep_can_find_positive_sigma**: True
- **it_is_an_identity_not_a_statistic**: σ = −(β₊ + β₋)/4πGR with β± ≥ 0 for any timelike shell

## T4_an_aligned_shell_can_be_ordinary — PASS
- **rows**:
    - radius=8, aligned_sigma=0.000420465, aligned_rest_mass=0.338158, anti_aligned_sigma=-0.0176495, anti_aligned_rest_mass=-14.1946
    - radius=20, aligned_sigma=6.23958e-05, aligned_rest_mass=0.313636, anti_aligned_sigma=-0.00761178, anti_aligned_rest_mass=-38.261
    - radius=60, aligned_sigma=6.72747e-06, aligned_rest_mass=0.304344, anti_aligned_sigma=-0.00261473, anti_aligned_rest_mass=-118.287
- **screened_mass**: 0.3
- **the_aligned_shell_is_ordinary**: True
- **the_anti_aligned_shell_is_exotic**: True
- **so_the_oppositely_glued_shell_moves_the_exotic_matter**: True

## T5_the_shell_supports_the_throat — PASS
- **rows**:
    - screened_mass=0, potential_shift_at_b0=0, force=-0, analytic_force=0
    - screened_mass=0.1, potential_shift_at_b0=0.04, force=0.008, analytic_force=0.008
    - screened_mass=0.3, potential_shift_at_b0=0.12, force=0.024, analytic_force=0.024
    - screened_mass=0.6, potential_shift_at_b0=0.24, force=0.048, analytic_force=0.048
- **the_force_opposes_closure**: True
- **it_grows_with_the_screened_mass**: True
- **it_matches_2_G_dM_over_b_squared**: True
- **zero_shell_gives_zero_force**: True

## T6_the_stability_window_and_what_screening_does — PASS
- **rows**:
    - interior_mass=1, sigma=-0.0246562, is_exotic=True, beta2_critical=-1.08333, stable_below_critical=True
    - interior_mass=0.9, sigma=-0.0254648, is_exotic=True, beta2_critical=-0.946332, stable_below_critical=True
    - interior_mass=0.8, sigma=-0.0262485, is_exotic=True, beta2_critical=-0.843891, stable_below_critical=True
    - interior_mass=0.7, sigma=-0.0270095, is_exotic=True, beta2_critical=-0.764847, stable_below_critical=True
    - interior_mass=0.5, sigma=-0.0284705, is_exotic=True, beta2_critical=-0.651786, stable_below_critical=True
- **throat_radius**: 5
- **screening_raises_the_critical_beta2**: True
- **the_window_always_needs_negative_beta2**: True
- **the_throat_is_always_exotic**: True
- **so_screening_helps_but_does_not_rescue**: True

## T7_the_three_observables_disagree — PASS
- **observable_1_shell_sigma**: 6.23958e-05
- **observable_1_shell_is_exotic**: False
- **observable_2_force_on_throat**: 0.024
- **observable_2_supports_the_throat**: True
- **observable_3_beta2_critical**: -0.764847
- **observable_3_needs_negative_beta2**: True
- **throat_sigma**: -0.0270095
- **the_throat_is_still_exotic**: True
- **they_do_not_agree**: True
- **verdict**: the detached shell is ordinary and does support the throat, and the throat's own exotic requirement is untouched — three answers, three different signs

## T8_birkhoff_decouples_them — PASS
- **rows**:
    - radius=8, shell_rest_mass=0.338158, shell_sigma=0.000420465, shell_pressure=2.85163e-05, throat_sigma=-0.0270095, off_diagonal=0
    - radius=20, shell_rest_mass=0.313636, shell_sigma=6.23958e-05, shell_pressure=1.45137e-06, throat_sigma=-0.0270095, off_diagonal=0
    - radius=60, shell_rest_mass=0.304344, shell_sigma=6.72747e-06, shell_pressure=4.90654e-08, throat_sigma=-0.0270095, off_diagonal=0
    - radius=200, shell_rest_mass=0.301283, shell_sigma=5.99384e-07, shell_pressure=1.28478e-09, throat_sigma=-0.0270095, off_diagonal=0
- **shell_rest_mass_range**: (0.3012832723257919, 0.338157619558471)
- **shell_rest_mass_varies_by**: 1.12239
- **shell_sigma_varies_by**: 701.494
- **throat_sigma_spread**: 0
- **worst_off_diagonal**: 0
- **the_shells_are_genuinely_different**: True
- **the_throat_never_notices**: True
- **the_off_diagonal_vanishes**: True
- **but_that_is_structural_not_measured**: Birkhoff is assumed when the region between is written as Schwarzschild with a constant mass; the zero confirms the implementation is consistent with it and proves nothing more
- **so_there_is_no_two_mode_coupling**: True
- **why_it_matters**: spherical symmetry has no radiative channel — the same ℓ = 0 fact wave_constraints found for the scalar — so ℓ ≥ 2 internal modes are not a later refinement but the only place a genuine throat-shell coupling could live

## T9_no_flat_direction — PASS
- **rows**:
    - scale=1, eigenvalues=[0.008693087975390134, 0.11462222222222218], off_diagonal=0, scaled_eigenvalues=[0.008693087975390134, 0.11462222222222218]
    - scale=2, eigenvalues=[0.0021732719938475334, 0.028655555555555545], off_diagonal=0, scaled_eigenvalues=[0.008693087975390134, 0.11462222222222218]
    - scale=4, eigenvalues=[0.0005433179984618834, 0.007163888888888886], off_diagonal=0, scaled_eigenvalues=[0.008693087975390134, 0.11462222222222218]
- **eigenvalues_scale_as_inverse_length_squared**: True
- **worst_scaling_drift**: 0
- **smallest_absolute_eigenvalue**: 0.000543318
- **no_flat_direction**: True
- **the_scale_is_set_by_the_masses_not_by_a_boundary**: True

## T10_assessment — PASS
- **n_passed**: 9
- **n_total**: 9

## verdict — THE_EXOTIC_MATTER_MOVES_BUT_DOES_NOT_LEAVE

THE EXOTIC MATTER MOVES BUT DOES NOT LEAVE. Of the three observables the hope needed, two come back positive and the decisive one comes back negative. An ANTI-ALIGNED detached shell — the oppositely-glued one the proposal is about — is necessarily exotic: its two extrinsic-curvature terms add rather than cancel, so sigma = -(beta+ + beta-)/4piGR with both roots non-negative. Over 200000 random Schwarzschild, de Sitter and Reissner-Nordstrom pairs there is not one counterexample, and there cannot be, because it is an identity and not a statistic. The same identity applies to the throat, which IS a minimal surface, so no arrangement of bulk content can relieve it. What an ORDINARY aligned shell can do is real but different: it screens mass, pushing the throat outward with F = 2G dM/b^2, and it raises the critical beta2 monotonically as it screens more, enlarging the throat's stability window. Both normal modes can then be positive at once. But the window never reaches beta2 >= 0, and in exactly the configuration where the shell is ordinary and supporting, the throat's own sigma is -0.027 — still exotic. Three observables, three different signs, on one system: which is why they had to be reported separately rather than collapsed into a single stability verdict. The last finding is the one that shapes what comes next. Birkhoff decouples the two surfaces exactly: the shell's surface density varies 701x as it is moved at fixed screened mass, and the throat's sigma does not change in its last bit. Spherical symmetry has no radiative channel — the same l = 0 fact wave_constraints found for the scalar — so a genuine two-mode trapped resonator cannot exist here at all, and the l >= 2 internal modes are not a later refinement but the only place such a coupling could live.
