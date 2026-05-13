# Throat boundary-condition probe

**Run:** 2026-05-13T06:11:44+00:00

Opens the throat-dynamics thread (`docs/throat_dynamics_research_plan.md`) by testing whether the inner-boundary hard wall (Dirichlet) is uniquely required for a discrete Tangherlini radial spectrum, or whether other local BCs (Neumann, Robin) give ε-converged spectra that remove the regularization dependence.

Outer endpoint always Dirichlet at r = R_OUTER − ε. R\* = 1.262636 (closure-quantum cross-species fixed point). Reference points: ε at which Dirichlet ω = 1 is ε* = 3.5087e-04; closure-quantum reading is `7π/(100·5⁴)` = 3.5186e-04.

## ω(1, 0; R*, ε) per boundary condition

Each row gives the lowest eigenfrequency at the labeled BC and ε. Spread is `max(ω) − min(ω)` across the ε sweep.

| BC | ε = 1e-02 | ε = 5e-03 | ε = 1e-03 | ε = 5e-04 | ε = 1e-04 | ε = 5e-05 | ε = 1e-05 | spread | converged? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Dirichlet (hard wall) | 1.8953 | 1.5978 | 1.1754 | 1.0535 | 0.8454 | 0.7778 | 0.6544 | 1.2409 | no |
| Neumann (u′ = 0) | 1.0467 | 0.8766 | 0.6283 | 0.5571 | 0.4385 | 0.4010 | 0.3340 | 0.7127 | no |
| Robin κ = -10.0 | 1.9972 | 1.6697 | 1.2152 | 1.0860 | 0.8669 | 0.7962 | 0.6677 | 1.3295 | no |
| Robin κ = -1.0 | 2.5211 | 2.0891 | 1.4950 | 1.3274 | 1.0434 | 0.9518 | 0.7856 | 1.7355 | no |
| Robin κ = -0.1 | 0.9863 | 0.8146 | 0.5618 | 0.4887 | 0.3659 | 0.3267 | 0.2556 | 0.7307 | no |
| Robin κ = +0.1 | 1.1000 | 0.9301 | 0.6823 | 0.6113 | 0.4923 | 0.4546 | 0.3867 | 0.7133 | no |
| Robin κ = +1.0 | 1.3977 | 1.2078 | 0.9242 | 0.8397 | 0.6924 | 0.6437 | 0.5531 | 0.8446 | no |
| Robin κ = +10.0 | 1.8032 | 1.5316 | 1.1379 | 1.0228 | 0.8248 | 0.7602 | 0.6417 | 1.1615 | no |

## Compton-bridge ε per BC

For each BC, search for the ε at which ω(1, 0; R*, ε) = 1 (the Compton-bridge condition). A natural BC for the BAM framework would give this ε at a closure-quantum value.

| BC | ε at which ω = 1 |
|---|---:|
| Dirichlet (hard wall) | 3.5087e-04 |
| Neumann (u′ = 0) | 8.4538e-03 |
| Robin κ = -10.0 | 2.8754e-04 |
| Robin κ = -1.0 | 7.2996e-05 |
| Robin κ = -0.1 | (no bracket) |
| Robin κ = +0.1 | 6.8437e-03 |
| Robin κ = +1.0 | 1.6851e-03 |
| Robin κ = +10.0 | 4.2814e-04 |

## ω at the two reference ε values per BC

Reference 1: ε = ε*_Dirichlet = 3.5087e-04 (where Dirichlet hits ω = 1 exactly).
Reference 2: ε = `7π/(100·5⁴)` = 3.5186e-04 (closure-quantum reading).

| BC | ω at ε*_Dirichlet | ω at ε_closure-quantum |
|---|---:|---:|
| Dirichlet (hard wall) | 0.999998 | 1.000403 |
| Neumann (u′ = 0) | 0.526239 | 0.526472 |
| Robin κ = -10.0 | 1.029427 | 1.029854 |
| Robin κ = -1.0 | 1.254125 | 1.254679 |
| Robin κ = -0.1 | 0.456890 | 0.457130 |
| Robin κ = +0.1 | 0.580342 | 0.580575 |
| Robin κ = +1.0 | 0.802198 | 0.802482 |
| Robin κ = +10.0 | 0.972077 | 0.972462 |

## Verdict

**Converged BCs (ω spread < 0.05 over 3 orders of magnitude in ε):** 0 of 8.

**No simple local BC at the inner endpoint removes the ε-dependence of ω(1, 0).** Every BC tested produces a spectrum that drifts with the regularization at the same order as Dirichlet — the spread is O(1.74).

This is the expected outcome: in tortoise coordinates the throat is asymptotically free (V → 0 as r* → −∞), so any local BC at finite ε produces only a *reflection phase* shift, not a change in the asymptotic spectrum. The discrete spectrum at finite ε is set by the box-width L(ε) which is log-divergent.

**Implication.** The physical inner boundary cannot be a purely local BC at the throat. The closure-quantum identification `ε = resistance / k_5⁴` is therefore either (i) the correct effective description with the hard-wall scheme acting as a mean-field stand-in for non-local throat physics, or (ii) waiting on a finite throat-thickness model (sub-target 2 of the research plan) where a soft confining potential replaces the wall.

**One important caveat.** The closure-quantum ε = `7π/(100·5⁴)` closes the Compton bridge *only* under Dirichlet. The table above shows ω at this same ε is 0.526 under Neumann, 1.030 under Robin κ = −10, and 0.580 under Robin κ = +0.1. The Compton-bridge ε varies by two orders of magnitude across BCs (from 7.3×10⁻⁵ to 8.5×10⁻³). The closure-quantum reading established in PR #18 is therefore *Dirichlet-specific* — not a BC-independent geometric identity. This is consistent with the closure-quantum scaffolding of PRs #15–18 having been derived entirely within the Dirichlet scheme; transport = 8π, resistance = 7π/100, γ = Σ V_max, etc. would all shift under a different BC.

## What this leaves open

Per `docs/throat_dynamics_research_plan.md`, three routes remain in the thread:

1. **Throat-thickness model (sub-target 2).** Replace the hard wall with a smooth confining potential `V_throat(r)` of finite thickness δ. Scan (V_0, δ) and check whether a natural closure-quantum value gives the locked spectrum without an ε regularization.
2. **Non-orientable identification (sub-target 1 cont'd).** The BAM throat has a T = iσ_y action with T² = −I. The natural BC may not be a single local condition but a *matching* condition across the throat that identifies the wavefunction with its image under T at the antipode. This is non-local in the radial coordinate.
3. **Reflection-phase analysis (sub-target 3).** The throat imposes a specific reflection phase φ(ω) on the asymptotic free waves. If φ has a natural form derived from the Tangherlini geometry or the throat T-action, the discrete spectrum emerges from matching φ to the outer BC.