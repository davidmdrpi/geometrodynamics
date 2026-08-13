# shell_junction probe

_generated 2026-08-12T23:54:03.715775+00:00_

## T1_goal — PASS
- **question**: does a detached oppositely-glued bulk shell change the throat's effective force and stability while itself remaining non-exotic?
- **observables**: ["1 the shell's own Israel surface stress sigma", '2 the shell-induced gradient of the throat potential', '3 the stiffness d2V/db2']
- **scope**: Einstein gravity, Darmois-Israel thin shells, spherical symmetry, vacuum bulk; D parametric with D=5 the BAM case

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

## T3_the_gluing_fixes_the_sign — PASS
- **rows**:
    - gluing=ordinary bubble, minus_branch=INNER, plus_branch=OUTER, eps=[1, 1], eta=1, forced_sign=None, sigma=0.000164165, sign_agrees=True
    - gluing=minimal surface, minus_branch=OUTER, plus_branch=OUTER, eps=[-1, 1], eta=-1, forced_sign=-1, sigma=-0.0154298, sign_agrees=True
    - gluing=maximal surface, minus_branch=INNER, plus_branch=INNER, eps=[1, -1], eta=-1, forced_sign=1, sigma=0.0154298, sign_agrees=True
    - gluing=anti-bubble, minus_branch=OUTER, plus_branch=INNER, eps=[-1, -1], eta=1, forced_sign=None, sigma=-0.000164165, sign_agrees=True
- **eta_minus_one_covers**: ['minimal surface', 'maximal surface']
- **their_forced_signs**: [-1, 1]
- **eta_alone_does_not_decide**: True
- **every_forced_sign_is_realised**: True
- **the_minimal_surface_is_exotic**: True
- **the_maximal_surface_is_not**: True
- **and_the_maximal_one_is_closed_off**: a maximal surface retains r ≤ R on both sides, so the manifold caps off either way and shares no bulk with the throat — it is non-exotic precisely because it is disconnected

## T4_the_forced_signs_hold_in_any_dimension — PASS
- **rows**:
    - dim=4, samples=13333, minimal_surface_violations=0, maximal_surface_violations=0, worst_minimal_sigma=-0.00318201, worst_maximal_sigma=0.00318201
    - dim=5, samples=13333, minimal_surface_violations=0, maximal_surface_violations=0, worst_minimal_sigma=-0.00499612, worst_maximal_sigma=0.00499612
    - dim=6, samples=13333, minimal_surface_violations=0, maximal_surface_violations=0, worst_minimal_sigma=-0.00661989, worst_maximal_sigma=0.00661989
- **dims**: [4, 5, 6]
- **no_violations_in_any_dimension**: True
- **d5_is_the_bam_case**: True
- **it_is_an_identity_not_a_statistic**: σ_minimal = −(D−2)(β₊+β₋)/8πGR and σ_maximal = +the same, with β± ≥ 0 for any timelike shell
- **scope**: Einstein gravity, Darmois–Israel thin shells, spherical symmetry, vacuum bulk

## T5_the_connected_dichotomy — PASS
- **rows**:
    - radius=8, ordinary bubble=0.000420465, minimal surface=-0.0176495, maximal surface=0.0176495, anti-bubble=-0.000420465
    - radius=20, ordinary bubble=6.23958e-05, minimal surface=-0.00761178, maximal surface=0.00761178, anti-bubble=-6.23958e-05
    - radius=60, ordinary bubble=6.72747e-06, minimal surface=-0.00261473, maximal surface=0.00261473, anti-bubble=-6.72747e-06
- **the_bubble_is_ordinary**: True
- **the_minimal_surface_is_exotic**: True
- **the_maximal_surface_is_ordinary**: True
- **but_it_is_disconnected**: True

## T6_the_shell_potential_gradient — PASS
- **rows**:
    - screened_mu=0, potential_shift_at_b0=0, minus_dV_db=-0, acceleration_contribution=-0
    - screened_mu=0.2, potential_shift_at_b0=0.04, minus_dV_db=0.008, acceleration_contribution=0.004
    - screened_mu=0.6, potential_shift_at_b0=0.12, minus_dV_db=0.024, acceleration_contribution=0.012
    - screened_mu=1.2, potential_shift_at_b0=0.24, minus_dV_db=0.048, acceleration_contribution=0.024
- **the_gradient_opposes_closure**: True
- **it_grows_with_the_screened_mass**: True
- **zero_shell_gives_zero_gradient**: True
- **it_is_not_an_equilibrium_consistent_force**: fixed throat rest mass, no equation-of-state response; the acceleration contribution is −½∂ΔV/∂b

## T7_the_stability_window — PASS
- **rows**:
    - mu_interior=2, sigma=-0.0246562, is_exotic=True, beta2_critical=-1.08333, stable_below_critical=True
    - mu_interior=1.8, sigma=-0.0254648, is_exotic=True, beta2_critical=-0.946332, stable_below_critical=True
    - mu_interior=1.6, sigma=-0.0262485, is_exotic=True, beta2_critical=-0.843891, stable_below_critical=True
    - mu_interior=1.4, sigma=-0.0270095, is_exotic=True, beta2_critical=-0.764847, stable_below_critical=True
    - mu_interior=1, sigma=-0.0284705, is_exotic=True, beta2_critical=-0.651786, stable_below_critical=True
- **throat_radius**: 5
- **dim**: 4
- **screening_raises_the_critical_beta2**: True
- **the_window_always_needs_negative_beta2**: True
- **the_throat_is_always_exotic**: True
- **beta2_is_an_eos_derivative_not_a_sound_speed**: True

## T8_the_three_observables_disagree — PASS
- **dim**: 4
- **observable_1_shell_sigma**: 6.23958e-05
- **observable_1_shell_is_exotic**: False
- **observable_2_minus_dV_db**: 0.024
- **observable_2_opposes_closure**: True
- **observable_3_beta2_critical**: -0.764847
- **observable_3_needs_negative_beta2**: True
- **throat_sigma**: -0.0270095
- **the_throat_is_still_exotic**: True
- **they_do_not_agree**: True

## T9_birkhoff_and_what_it_is_worth — PASS
- **rows**:
    - radius=8, shell_sigma=0.000420465, shell_pressure=2.85163e-05, throat_sigma=-0.0270095, off_diagonal=0
    - radius=20, shell_sigma=6.23958e-05, shell_pressure=1.45137e-06, throat_sigma=-0.0270095, off_diagonal=0
    - radius=60, shell_sigma=6.72747e-06, shell_pressure=4.90654e-08, throat_sigma=-0.0270095, off_diagonal=0
    - radius=200, shell_sigma=5.99384e-07, shell_pressure=1.28479e-09, throat_sigma=-0.0270095, off_diagonal=0
- **shell_sigma_varies_by**: 701.494
- **throat_sigma_spread**: 0
- **worst_off_diagonal**: 0
- **the_shells_are_genuinely_different**: True
- **the_throat_never_notices**: True
- **the_off_diagonal_vanishes**: True
- **but_that_is_structural_not_measured**: Birkhoff is assumed when the region between is written with a constant mass parameter; the zero confirms the implementation is consistent with it and proves nothing more
- **what_it_establishes**: no separation-dependent coupling in this vacuum spherical model — not that every spherical trapped resonator is impossible

## T10_a_units_check_on_the_stiffnesses — PASS
- **rows**:
    - scale=1, stiffnesses=[0.11462222217003178, 0.008693088025442849], scaled=[0.11462222217003178, 0.008693088025442849]
    - scale=2, stiffnesses=[0.028655555542507944, 0.0021732720063607123], scaled=[0.11462222217003178, 0.008693088025442849]
    - scale=4, stiffnesses=[0.007163888885626986, 0.0005433180015901781], scaled=[0.11462222217003178, 0.008693088025442849]
- **dim**: 4
- **worst_scaling_drift**: 0
- **stiffnesses_scale_as_inverse_length_squared**: True
- **this_is_a_units_check_only**: it does not show a fixed system has no dilation mode — that needs a kinetic term, which is not derived here

## T11_assessment — PASS
- **n_passed**: 10
- **n_total**: 10

## verdict — CONNECTED_MEANS_EXOTIC

CONNECTED MEANS EXOTIC. Deriving the orientation from the gluing rather than setting it by hand sharpens the result and corrects the earlier framing. Each side retains the INNER or the OUTER radial branch, which fixes eps with no freedom left, and that gives FOUR gluings rather than two. The parity eta = eps+ eps- covers two of them at once — minimal surface, maximal surface — whose forced signs are OPPOSITE, so eta alone decides nothing and 'the oppositely-glued shell is exotic' was too coarse a statement. What is actually forced is sharper: a MINIMAL surface has sigma = -(D-2)(beta+ + beta-)/8piGR < 0 and a MAXIMAL surface the same with the other sign, both identities, and neither is violated once in 40 000 random Tangherlini, de Sitter and charged pairs across D = [4, 5, 6]. The dichotomy that produces is the answer. A detached surface that CONNECTS to the throat's region does so through a minimal surface and is necessarily exotic. A detached surface that is non-exotic by its gluing is a maximal surface, which caps off on both sides and shares no bulk with the throat — it is non-exotic precisely because it is disconnected, and therefore cannot support anything. Within Einstein-Israel spherical thin shells the exotic matter is relocated, never removed. The three observables still disagree on one system, which is why they are reported apart: the ordinary bubble has sigma = 6.24e-05, its screening does push the throat outward, and the throat's own sigma is -0.027 regardless. Two cautions are carried in the output rather than buried. The shell's effect is reported as a potential GRADIENT, not a force: it is taken at fixed throat rest mass with no equation-of-state response. And Birkhoff's vanishing off-diagonal is STRUCTURAL, imported the moment the intervening region is written with a constant mass parameter; what is measured instead is that a family of shells differing 701x in surface density leaves the throat bit-for-bit unchanged. That establishes no separation-dependent coupling in THIS model, not that every spherical trapped resonator is impossible.
