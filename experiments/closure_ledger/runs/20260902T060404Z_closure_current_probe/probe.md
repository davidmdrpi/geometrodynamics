# Closure-current fork probe — 7/7

Pre-registration: `docs/closure_current_prereg.md @ f954e3d`. Verdict: **`FORK_UNRESOLVED_BY_CURRENT_STRUCTURES`**

| candidate | E(1) | S_max |
|--|--|--|
| `POSITIVE_COAREA` | `+0.39850` | `2.1423` |
| `HOLONOMY_WEIGHTED_COAREA` | `+0.54030` | `2.8284` |

R1: `q5 = -q0 G^-1` to `4.4e-16`; frame holonomy = Ω(x;u,−v) to `2.2e-15`; lift = cos(Ω/2)+sin(Ω/2)x to `9.3e-16`.
C1: e^{iΩ/2} = sgn D on the closure set to `1.2e-16`.
R2: holonomy-weighted current gives P_like = (1+cos γ)/4, P_unlike = (1−cos γ)/4: deviations γ=0.3: `5.6e-17`, γ=1.0: `5.6e-17`, γ=2.0: `5.6e-17`.
R3: sector prior ratio 0.5/1/2 → E(1) = `0.0751`, `0.3985`, `0.6460`; marginals 1/2 throughout.
R4: the repository's condition is *stationarity: the history has extremal action*, implemented: False. The phase-stationarity **proxy** is analytically disjoint from sharp closure: `∇Ω|_{N=0} = 2(u×v)/D`, minimum norm `0.5107` on the closure set (finite-difference residual `3.4e-07`); the only points with `D = 0` are `x = −u, −v`, where the `arg` chart is singular, not stationary.
Oriented-current sector integrals: `∫_Γ D_s dσ = 2π(1+u·v) ≥ 0` (residual `0.0e+00`) — the cancellation is internal, so the holonomy weighting is wave interference, not negative event probabilities.

## Checks

* PASS — R1_pin_loop_reduces_to_triangle_with_partner_sign
* PASS — C1_branch_holonomy_is_sign_D
* PASS — R2_holonomy_weighted_current_is_born_without_projectors
* PASS — R3_sector_prior_chosen_marginals_fixed_E_moves
* PASS — R4_stationarity_unimplemented_and_proxy_is_disjoint_from_closure
* PASS — sector_integrals_of_the_oriented_current_are_nonnegative
* PASS — R5_pin_supplies_label_not_weight

## Where the gap is

Not the spin structure (derived), not the loop (derived), not Bell (evaded by measurement dependence either way), not the relative outcome law (the oriented current gives it analytically). It is whether observed event frequencies are the positive count of closed histories or their oriented, holonomy-weighted sum. Nothing classical in the repository decides; the rule forbids deciding it by the answer. Sharper form, recorded as an open audit item and not a criterion: if the observable is an integral of a section of the local system the Pin/Hopf data define over the closure locus, the sign is geometrically mandatory; if the physical object is a measure on histories, positivity forces |D|.

## Dependency ledger

* `loop` = loop( source pair [x, -x via J], geodesic realignment at A, B [chosen], J = L_{-j} [derived] ) -> triangle x -> u -> -v -> x [derived]
* `closure holonomy` = -sgn D(x, u, -v) [derived from J^2 = -1 and the lift G]
* `positive coarea` = |D|/(2|u x v|) [derived conditioning; counting measure on sectors: chosen]
* `oriented current` = e^{i Phi} x coarea [derived label; adopting it as weight: the open step]
* `sector prior` = counting [chosen]
