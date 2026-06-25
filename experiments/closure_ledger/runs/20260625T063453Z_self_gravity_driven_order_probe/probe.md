# Self-gravity-driven throat-order instability: bind, or drive a new order parameter? (PR #178)

**Run:** 2026-06-25T06:34:53+00:00

Couples the #176/#177 weak-field self-gravity solver to the #178 throat-order field and asks: does the gravitational concentration merely BIND the wave, or DRIVE the order parameter? *(QFT on the classical throat, not quantum gravity.)*

- **Coupling**: V(q;ρ) = ½(a₀ − gρ)|q|² + (λ/4)|q|⁴; a(ρ)=a₀−gρ; ρ_c=a₀/g
- **Merely bound**: sub-threshold ρ_peak < ρ_c → q = 0 (no geometric order)
- **Drives order**: super-threshold collapse ρ_peak > ρ_c → q nucleates at the core
- **Gravitational**: G=0 never crosses ρ_c → no order; restoring G nucleates it

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | couple #176/#177 self-gravity to the #178 order field | **PASS** |
| T2 | `T2_density_coupled_landau_potential` | density-coupled Landau potential V(q;ρ); the critical ρ_c | **PASS** |
| T3 | `T3_merely_bound_no_order` | merely bound: sub-threshold ρ<ρ_c → q=0 (no order) | **PASS** |
| T4 | `T4_drives_order_nucleation` | drives order: super-threshold ρ>ρ_c → q nucleates at the core | **PASS** |
| T5 | `T5_ordering_is_gravitational` | gravitational: G=0 never crosses ρ_c → no order | **PASS** |
| T6 | `T6_dynamical_order_front` | dynamical: q switches on as ρ_peak(t) crosses ρ_c | **PASS** |
| T7 | `T7_honest_scope` | honest scope (one-way coupling; effective constants) | **PASS** |
| T8 | `T8_assessment` | SELF_GRAVITY_DRIVES_THROAT_ORDER | **PASS** |

## Concentration drives the order parameter

| regime | M | G | ρ_peak | vs ρ_c | max \|q\| | outcome |
|---|---:|---:|---:|---|---:|---|
| sub-threshold | 1 | 1 | 0.0617 | < 0.3 | 1.1e-10 | merely bound (no order) |
| super-threshold | 3 | 1 | 0.8957 | > 0.3 | 0.6813 | order nucleates |
| gravity off | 3 | 0 | 0.1839 | < 0.3 | 1.5e-08 | no order (it is gravity) |

## Verdict

**SELF_GRAVITY_DRIVES_THROAT_ORDER_ABOVE_A_CRITICAL_CONCENTRATION_NOT_MERELY_BINDING.** DRIVES ORDER — NOT MERELY BINDING. Coupling the #176/#177 self-gravity solver to the #178 throat-order field, the gravitational concentration drives the order parameter off zero.

THE COUPLING. The matter density ρ = |ψ|² is the control field of a density-dependent Landau potential V(q; ρ) = ½(a₀ − gρ)|q|² + (λ/4)|q|⁴; the order field's effective mass² a(ρ) = a₀ − gρ changes sign at ρ_c = 0.300.

MERELY BOUND (sub-threshold). A sub-threshold packet (M = 1) reaches only ρ_peak = 0.062 < ρ_c, and the order field relaxes to zero (max |q| = 1.1e-10) — bound, but no geometric order.

DRIVES ORDER (super-threshold). Above the mass threshold (M = 3) the collapse drives ρ_peak = 0.896 > ρ_c and the order field NUCLEATES (max |q| = 0.681) — a localized symmetry-broken domain at the density peak (the throat core).

GRAVITATIONAL. With gravity off (G = 0) the same mass never crosses ρ_c (ρ_peak = 0.184) and no order nucleates; restoring gravity it does — the ordering inherits the M_c ∝ 1/G gravity of #176/#177.

DYNAMICAL. Driving q by the time-dependent ρ_peak(t) of the collapse, the order parameter switches on (step 380) after the density crosses ρ_c (step 210) — a moving order front following the gravitational concentration.

SCOPE. One-way coupling (q's metric back-reaction pending); the constants a₀, g, λ — and so ρ_c — are effective (the existence of a gravitationally-crossed concentration threshold is the result, not its microscopic value); weak-field. The collapse of #176/#177 carries the matter density across the order transition and nucleates the throat-order field of #178.
