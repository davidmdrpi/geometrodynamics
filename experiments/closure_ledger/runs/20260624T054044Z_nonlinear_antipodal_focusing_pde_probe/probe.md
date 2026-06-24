# The nonlinear antipodal focusing PDE sandbox (PR #175)

**Run:** 2026-06-24T05:40:44+00:00

**Can a continuous time-dependent geometry actually evolve into the discrete sector?** A nonlinear antipodal-focusing PDE sandbox (a focusing NLS on the antipodal ring; the discrete sector = the winding Q). *(QFT on the classical throat, not quantum gravity.)*

- **Smooth conserves**: winding frozen under smooth evolution (discrete sector locked out; #174)
- **The gate**: Q changes only at an amplitude-zero node, forced at the antipode (the focus)
- **The threshold**: focusing reaches the core only above a critical mass (#58/#166)
- **The jump**: quantized ±1 (discrete response to a smooth drive)
- **Answer**: yes, but only through the antipodal focusing singularity — never smoothly

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the question and the sandbox | **PASS** |
| T2 | `T2_sandbox_model` | the model (winding Q = the discrete sector) | **PASS** |
| T3 | `T3_winding_conserved_smooth_evolution` | winding conserved under smooth evolution (#174 dynamical) | **PASS** |
| T4 | `T4_topological_gate_at_antipode` | the topological gate: a node forced at the antipode | **PASS** |
| T5 | `T5_focusing_threshold` | the focusing threshold: disperse vs concentrate-to-core | **PASS** |
| T6 | `T6_quantized_jump` | the jump is quantized ±1 (discrete response to a smooth drive) | **PASS** |
| T7 | `T7_synthesis_and_scope` | synthesis + honest scope (reduced 1D model) | **PASS** |
| T8 | `T8_assessment` | CONTINUOUS_REACHES_DISCRETE_ONLY_THROUGH_CAUSTIC | **PASS** |

## The focusing threshold (critical-NLS mass scan)

| mass | peak growth | outcome |
|---:|---:|---|
| 0.99 | ×1.0 | disperse |
| 1.588 | ×1.0 | disperse |
| 1.94 | ×1.51 | concentrate |
| 2.534 | ×5.1 | concentrate |

## Verdict

**CONTINUOUS_REACHES_DISCRETE_SECTOR_ONLY_THROUGH_ANTIPODAL_FOCUSING_SINGULARITY.** YES — BUT ONLY THROUGH THE CAUSTIC. A continuous time-dependent geometry can enter the discrete sector, but never by smooth deformation: only by developing a focusing singularity at the antipode.

SMOOTH EVOLUTION CONSERVES THE WINDING. A smooth Q = 1 field (|ψ| > 0) keeps its winding exactly (1.00000 → 1.00000) under the nonlinear evolution — the discrete sector is dynamically locked out, the dynamical confirmation of #174.

THE GATE IS AT THE ANTIPODE. Changing Q forces an amplitude-zero node (min|ψ| → 6e-17 along the interpolation), located EXACTLY at the antipode (χ/π = 1.0) — the focus. The discrete sector is gated by a singular core at the antipodal caustic, and the antipodal focusing is precisely what drives the field there.

THE THRESHOLD. Whether the nonlinear focusing reaches that core depends on the mass: below ~1.588 it disperses (continuous, Q frozen), above ~1.94 it concentrates toward the core (nucleation) — the #58/#166 threshold, now simulated nonlinearly.

THE JUMP IS QUANTIZED. At the core the winding changes by exactly ±1 (ΔQ = +1) — a discrete response to the smooth focusing drive. Continuous driver, discrete response, mediated by the node.

SCOPE. A reduced 1D ring model (Q proxies the discrete k, the collapse core proxies throat nucleation, the critical-NLS collapse is marginal): the conceptual answer is robust, the numbers model-dependent.
