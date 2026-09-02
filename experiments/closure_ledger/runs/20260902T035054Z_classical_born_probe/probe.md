# Classical Born-rule probe — 10/10

Pre-registration: `docs/classical_born_prereg.md @ 7ff6e41`. Verdict: **`BORN_REQUIRES_AN_IMPORTED_MEASURE_OR_DETECTOR_LAW`**

| candidate | coupling | induced f(θ) under fibre Haar | max miss from Born |
|--|--|--|--|
| C1 | sign of a linear functional of the frame | arccos family with plateaus | `0.109` (best, α/ρ = 0.80) |
| C2 | classical Malus intensities | step | `0.500` |
| C3 | spinor + frame harmonics | irregular | `0.500` (best of four) |
| C4 | repository winding Stern–Gerlach | step | `0.500` |
| C5 | sign(a·(x + κ y)), y Haar on S² | clip((1 + cos θ/κ)/2) | κ=0.5: `0.250`, κ=0.9: `0.050`, κ=1.0: `0.000`, κ=1.1: `0.045`, κ=2.0: `0.250` |

Monte Carlo at κ = 1: θ=0.3: `0.9776`, θ=1.0: `0.7701`, θ=2.0: `0.2913`

## Checks

* PASS — H1_symmetry_plus_haar_do_not_select_born
* PASS — C1_linear_family_misses_by_pre_registered_amount
* PASS — C2_intensity_is_the_step
* PASS — C3_two_harmonic_misses_by_at_least_half
* PASS — C4_repository_detector_is_the_step
* PASS — C5_archimedes_born_iff_kappa_one
* PASS — C5_monte_carlo
* PASS — basin_control_tuned_reproduces_born
* PASS — measure_control_fires
* PASS — reversal_control

## Typology

C: deterministic hidden outcome. The only classical route to Born found (C5) is Bell 1964 / Kochen-Specker 1967: outcomes fixed by (x, y), probabilities from ignorance of y. Outcome D (setting-dependent ensemble from a derived global boundary problem) was not found and nothing here supplies it.

## Dependency ledger

* `f under fibre Haar` = f( rotational covariance [derived], Haar on S^1 [derived], basin shape [from the coupling: derived for C1-C4; tuned for the control] )
* `C5 Born` = ( Haar on S^2 [chosen: the detector mouth is unprepared], kappa = 1 [chosen], D = sign(a.(x + kappa y)) [chosen] )

Nothing about composition or CHSH. Nothing about the field sector beyond the spin-frame degree of freedom.
