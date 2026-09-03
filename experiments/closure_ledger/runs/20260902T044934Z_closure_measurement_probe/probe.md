# Closure measurement-dependence probe — 14/14

Pre-registration: `docs/closure_measurement_dependence_prereg.md @ 1b0144e`. Verdict: **`CLOSURE_INDUCES_SETTING_DEPENDENT_SOURCE_MEASURE_NO_SIGNALLING_NOT_BORN`**

| γ | E (quadrature) | E (closed form) | cos γ | P(A=+) |
|--|--|--|--|--|
| `0.3000` | `+0.82095` | `+0.82095` | `+0.95534` | `0.5000` |
| `0.7854` | `+0.53557` | `+0.53557` | `+0.70711` | `0.5000` |
| `1.0000` | `+0.39850` | `+0.39850` | `+0.54030` | `0.5000` |
| `1.5708` | `+0.00000` | `+0.00000` | `+0.00000` | `0.5000` |
| `2.0000` | `-0.30350` | `-0.30350` | `-0.41615` | `0.5000` |

CHSH at (0, π/2, π/4, −π/4): `2.1423`; maximum over settings: `2.1423` at `(0.0, -1.571, -0.785, 0.785)` (below 2√2 = 2.8284: True).
Signed variant: `E = cos γ` (deviation `0.0e+00`), `D/4 = Re Tr(P_x P_u P_v)` to `2.2e-16`; negative-weight fraction of the closure circle for like outcomes: γ=0.5: `0.0795`, γ=1.0: `0.1592`, γ=2.0: `0.3183`.
Window Monte Carlo at γ = 1: ε=0.4: `0.3729`, ε=0.2: `0.3923`, ε=0.1: `0.3973`, ε=0.05: `0.3972` → coarea `0.3985`.
Measurement dependence: non-coplanar settings TV `1.000`; coplanar settings share the circle, TV of the densities `0.0182` (standard pairs `0.0631`).
Local-detector control: CHSH `0.9397`; Gaussian-width control: σ=0.6: `0.1700`, σ=0.3: `0.2925`, σ=0.1: `0.3804`, σ=0.03: `0.3947`.

## Checks

* PASS — P1_measurement_dependence_non_coplanar_maximal
* PASS — P1_measurement_dependence_coplanar_in_the_density
* PASS — P2_no_signalling
* PASS — P3_quadrature_equals_closed_form
* PASS — P3_pre_registered_values
* PASS — P3_not_cos
* PASS — P4_bell_violation_standard_angles
* PASS — P4_below_tsirelson
* PASS — P5_signed_is_bargmann_and_cos
* PASS — P6_strict_zero_variant
* PASS — P7_window_converges
* PASS — control_two_leg_loop_no_conditioning
* PASS — control_local_detectors_respect_bell
* PASS — control_gaussian_width_matters

## Typology

D: the admissible source ensemble depends on both future settings through the closure constraint (support = the great circle through a and b), with exact no-signalling. It is not a deterministic detector: the outcome is which closed history is realised. Bell's theorem is evaded by measurement dependence, S = 2.1423 at the standard angles, but the correlation is not the quantum one; the quantum one is the signed closure measure.

## Dependency ledger

* `rho(x|a,b)` = rho( Haar on S^2 [invariant prior; pair direction physical: chosen], closure axiom Omega = 0 mod 2pi [repository axiom], geodesic-realignment detection model [chosen], coarea conditioning [derived; window limit] )
* `E(gamma)` = E( rho(x|a,b) [above], outcome signs as history boundary data [chosen: D-type] )
* `Born` = signed closure measure [identity: Re Bargmann] -- not a probability; imported if used
