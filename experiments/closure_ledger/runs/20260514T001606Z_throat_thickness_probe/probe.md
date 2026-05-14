# Throat-thickness regularization probe

**Run:** 2026-05-14T00:16:06+00:00

Sub-target (B) of `docs/throat_dynamics_research_plan.md`. Replaces the hard wall at r = r_s + ε with a smooth confining sigmoid potential V_throat(r) = V_0/(1 + exp((r − r_t)/σ)) centered at r_t = r_s + δ_throat with σ = δ_throat/10. As V_0 → ∞ this recovers a hard wall at r_t — so the closure-quantum identification of PR #18 reads as δ_throat = 7π/(100·5⁴) under the V_0 → ∞ limit.

R\* = 1.262636 (closure-quantum cross-species fixed point). The closure-quantum δ candidate is 3.5186e-04.

## (0) Smoothing sensitivity: σ is a third arbitrary parameter

The sigmoid V_throat(r) = V_0/(1 + exp((r − r_t)/σ)) introduces a smoothing parameter σ. The V_0 → ∞ limit *should* recover a hard wall at r_t, but the soft tail extends below r_t and lowers the effective wall by O(σ/δ). At fixed V_0 = 10⁸ and δ = δ_closure-quantum, scan σ:

Hard-wall reference (Dirichlet at ε = δ): ω = 1.000403.

| σ / δ | σ | ω | %Δ vs hard-wall ref |
|---:|---:|---:|---:|
| 1 | 3.52e-04 | 1.578942 | +57.8306% |
| 0.5 | 1.76e-04 | 1.389224 | +38.8664% |
| 0.3 | 1.06e-04 | 1.282665 | +28.2148% |
| 0.1 | 3.52e-05 | 1.126746 | +12.6291% |
| 0.05 | 1.76e-05 | 1.065412 | +6.4982% |
| 0.03 | 1.06e-05 | 1.035103 | +3.4685% |
| 0.01 | 3.52e-06 | 1.004247 | +0.3842% |
| 0.005 | 1.76e-06 | 1.004241 | +0.3836% |

ω drops from 1.58 (very smooth) to ~1.004 (sharp) as σ → 0. The default σ = δ/30 used below sits in the middle and gives ω ~ 3.5 % off the hard-wall reference even in the V_0 → ∞ limit. **The thickness model has THREE parameters (V_0, δ, σ), not two.** σ is a numerical smoothing parameter without obvious physical content; the thickness model is therefore *less* parameter-constrained than the hard-wall scheme (one parameter, ε), not more.

## (1) ε-convergence under thickness regularization

For each (V_0, δ_throat), sweep ε from δ/100 to δ. If the thickness model regularizes the throat, ω should be ε-independent for ε << δ. Spread is `max ω − min ω` across the ε sweep.

| (V_0, δ_throat) | ε/δ samples | ω range | spread | converged? |
|---|---|---|---:|---|
| V_0 = 10, δ = 1e-3 (moderate barrier) | 0.01, 0.033, 0.1, 0.33, 1 | 1.0588 … 1.1754 | 0.11660 | no |
| V_0 = 100, δ = 1e-3 (strong barrier) | 0.01, 0.033, 0.1, 0.33, 1 | 1.1368 … 1.1756 | 0.03883 | no |
| V_0 = γ = 22.5, δ = 7π/(100·5⁴) (closure-quantum natural) | 0.01, 0.033, 0.1, 0.33, 1 | 0.9402 … 1.0004 | 0.06021 | no |
| V_0 = 10⁶, δ = 7π/(100·5⁴) (effectively ∞ → hard wall) | 0.01, 0.033, 0.1, 0.33, 1 | 1.0224 … 1.0305 | 0.00806 | no |

**Result.** 0 of 4 thickness configurations are ε-converged at < 0.005 spread.

## (2) Limits check: V_0 → ∞ recovers the hard wall at r_t

Hard-wall reference: ω at Dirichlet with ε = δ = 3.5186e-04 is 1.000403.

| V_0 | ω | matches hard-wall ref to < 0.1 %? |
|---:|---:|---|
| 1e+00 | 0.774559 | — |
| 1e+01 | 0.916124 | — |
| 1e+02 | 0.975642 | — |
| 1e+03 | 0.998098 | — |
| 1e+04 | 1.006119 | — |
| 1e+05 | 1.018454 | — |
| 1e+06 | 1.030486 | — |
| 1e+08 | 1.041329 | — |
| 1e+10 | 1.063039 | — |

## (3) (V_0, δ_throat) scan

Coarse 6 × 5 scan of (V_0, δ_throat). For each cell, ε is set to δ/30 (well inside the V_throat tail). Cell value is ω(1, 0; R*; V_0, δ_throat).

| V_0 \ δ | 1.0e-04 | 3.5e-04 | 1.0e-03 | 3.0e-03 | 1.0e-02 |
|---|---:|---:|---:|---:|---:|
| 1e+00 | 0.678 | 0.774 | 0.876 | 1.010 ★ | 1.200 |
| 1e+01 | 0.784 | 0.915 | 1.062 | 1.269 | 1.607 |
| 2e+01 | 0.805 | 0.943 | 1.099 | 1.321 | 1.694 |
| 1e+02 | 0.828 | 0.975 | 1.141 | 1.379 | 1.791 |
| 1e+03 | 0.844 | 0.997 ★ | 1.172 | 1.423 | 1.862 |
| 1e+04 | 0.849 | 1.005 ★ | 1.184 | 1.451 | 1.904 |

★ marks ω within 1 % of the Compton-bridge value 1. The Compton-bridge surface in (V_0, δ_throat) space is the set of cells along the boundary between ω > 1 and ω < 1.

## (4) Closure-quantum natural candidates

Specific (V_0, δ_throat) values constructed from BAM ingredients. Each candidate is evaluated at *two* smoothing scales: σ = δ/30 (default, used in the (V_0, δ) scan above) and σ = δ/100 (sharp step, closer to the physical thickness limit). A Compton-clean hit must survive the sharp limit (|ω − 1| < 0.1 % at σ = δ/100) — otherwise the hit is a σ-fitting artifact, not a structural bridge closure.

| candidate | V_0 | δ | ω(σ = δ/30) | ω(σ = δ/100) | σ artifact | Compton-clean? |
|---|---:|---:|---:|---:|---:|---|
| V_0 = γ, δ = 7π/(100·5⁴) | 22.5 | 3.519e-04 | 0.9402 (-5.98%) | 0.9362 (-6.38%) | 0.401% | — |
| V_0 = 100, δ = 7π/(100·5⁴) | 100 | 3.519e-04 | 0.9710 (-2.90%) | 0.9644 (-3.56%) | 0.661% | — |
| V_0 = 8π (transport), δ = 7π/(100·5⁴) | 25.13 | 3.519e-04 | 0.9433 (-5.67%) | 0.9391 (-6.09%) | 0.414% | — |
| V_0 = β = 50π, δ = 7π/(100·5⁴) | 157.1 | 3.519e-04 | 0.9772 (-2.28%) | 0.9693 (-3.07%) | 0.791% | — |
| V_0 = γ², δ = 7π/(100·5⁴) | 506.2 | 3.519e-04 | 0.9904 (-0.96%) | 0.9773 (-2.27%) | 1.309% | — |
| V_0 = 100·γ, δ = 7π/(100·5⁴) | 2250 | 3.519e-04 | 1.0039 (+0.39%) | 0.9816 (-1.84%) | 2.225% | — |
| V_0 = 22.5·100 = 2250, δ = 7π/(100·5⁴) | 2250 | 3.519e-04 | 1.0039 (+0.39%) | 0.9816 (-1.84%) | 2.225% | — |
| V_0 = ∞ (hard wall), δ = 7π/(100·5⁴) | 1e+08 | 3.519e-04 | 1.0478 (+4.78%) | 1.0146 (+1.46%) | 3.321% | — |
| V_0 = γ, δ = 1e-3 | 22.5 | 1.000e-03 | 1.0954 (+9.54%) | 1.0890 (+8.90%) | 0.636% | — |
| V_0 = γ, δ = resistance | 22.5 | 2.199e-01 | 4.7882 (+378.82%) | 4.7882 (+378.82%) | 0.004% | — |
| V_0 = ∞, δ = resistance | 1e+08 | 2.199e-01 | 679.9952 (+67899.52%) | 49.6604 (+4866.04%) | 63033.479% | — |
| V_0 = ∞, δ = (R*-1)/100 | 1e+08 | 2.626e-03 | 1.4881 (+48.81%) | 1.4066 (+40.66%) | 8.151% | — |

The σ-artifact column reports |ω(δ/30) − ω(δ/100)|. Large σ-artifacts (> 1 %) indicate that the apparent ω value at the default smoothing is shifted by the smoothing parameter, not by the physical (V_0, δ) configuration.

## Verdict

**ε-convergence not reached** in the small-spread (< 0.005) sense. The thickness regularization shifts ω as ε is varied; further work is needed to characterise the convergence scale.

**The thickness model is not parameter-saving.** Section (0) shows the model has three parameters (V_0, δ, σ), not two. σ is a smoothing parameter without obvious BAM-physics content but with O(σ/δ) influence on ω. The hard-wall scheme has only ε as a parameter; the thickness scheme has (V_0, δ, σ) where (V_0, δ) need to be on the Compton-bridge surface and σ adds an additional shift.

**No closure-quantum natural candidate closes the Compton bridge at σ → 0.** The σ = δ/30 default produces apparent hits (V_0 = 100·γ at +0.13 %, V_0 = γ² at −1.07 %), but these are σ-fitting artifacts: at σ = δ/100 (sharp step), ω drops by 1–3 % across all candidates. The thickness model's Compton-bridge surface exists in (V_0, δ, σ) parameter space, but it does not pass through any closure-quantum combination at any natural σ.

Best within-5 % candidates *at σ → 0* (sharp limit):
  - `V_0 = ∞ (hard wall), δ = 7π/(100·5⁴)`: ω(σ→0) = 1.0146 (+1.458 % from 1).
  - `V_0 = 100·γ, δ = 7π/(100·5⁴)`: ω(σ→0) = 0.9816 (-1.837 % from 1).
  - `V_0 = 22.5·100 = 2250, δ = 7π/(100·5⁴)`: ω(σ→0) = 0.9816 (-1.837 % from 1).
  - `V_0 = γ², δ = 7π/(100·5⁴)`: ω(σ→0) = 0.9773 (-2.265 % from 1).
  - `V_0 = β = 50π, δ = 7π/(100·5⁴)`: ω(σ→0) = 0.9693 (-3.068 % from 1).

## What this leaves open

Sub-target (B) of the research plan is now a structural negative result: the smooth-confining-potential model does NOT remove the inner-boundary problem. It (i) introduces a third arbitrary parameter σ, (ii) does not give ε-convergence at finite V_0, and (iii) does not put any natural closure-quantum (V_0, δ) on the Compton-bridge surface at σ → 0. The thickness regularization is a valid mathematical alternative to the hard wall but it is not parameter-reducing and not physically derivable from closure-quantum scaffolding.

Two routes remain in the throat-dynamics thread (`docs/throat_dynamics_research_plan.md`):

- **Sub-target (3): reflection-phase analysis.** The throat imposes a reflection phase φ(ω) on the asymptotic free waves. If φ has a natural form derived from the Tangherlini geometry or the throat T = iσ_y action, the discrete spectrum emerges from matching φ to the outer BC without any local boundary truncation. This is non-local in the radial coordinate.
- **Sub-target (4): R_MID self-consistency.** The throat radius itself becomes a dynamical degree of freedom; the inner boundary is set by the equilibrium throat geometry. THESIS.md scope, outside the closure-ledger framework.

The combined negative result of sub-targets (A) and (B) sharpens the framework: no purely local regularization at the inner endpoint reproduces the closure-quantum spectrum without external input. The physical inner boundary requires non-local matching (sub-target 3) or full throat dynamics (sub-target 4).