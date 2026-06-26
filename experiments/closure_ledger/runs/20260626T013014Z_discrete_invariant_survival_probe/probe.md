# Discrete invariant survival on the ψ–Φ–q throat-soliton (PR #181)

**Run:** 2026-06-26T01:30:14+00:00

Shows the throat's discrete winding invariant `Q=(1/2π)∮∇φ` survives on the #180 continuous self-consistent ψ–Φ–q soliton — conserved to machine precision under continuous evolution while `|q|>0`, changing only through `|q|=0` (the #182 phase slip). *(QFT on the classical throat, not quantum gravity.)*

- **Quantized**: Q = k exactly on an ordered loop of the soliton
- **Survival**: machine-precision under continuous evolution while |q|>0 (wave: all k; dissipative: sustained k=1,3)
- **Criterion**: Q changes ONLY through |q|=0 (the unsustained k=5 slips — the #182 event)
- **Rigidity**: unchanged under all |q|>0-preserving homotopies

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the discrete winding invariant survives on the #180 soliton | **PASS** |
| T2 | `T2_quantized_invariant_on_the_soliton` | quantized invariant Q=k on an ordered loop; sustainability A² | **PASS** |
| T3 | `T3_survival_continuous_wave_evolution` | survival under continuous wave evolution (all k; ΔW~machine) | **PASS** |
| T4 | `T4_survival_dissipative_relaxation` | survival under dissipative relaxation (sustained k=1,3) | **PASS** |
| T5 | `T5_invariant_changes_only_through_q_zero` | criterion: Q changes only through |q|=0 (k=5 slip → #182) | **PASS** |
| T6 | `T6_rigidity_under_homotopy` | rigidity under |q|>0-preserving homotopies (all k) | **PASS** |
| T7 | `T7_honest_scope` | honest scope + odd-k ladder / #183 bridge | **PASS** |
| T8 | `T8_assessment` | DISCRETE_INVARIANT_SURVIVES | **PASS** |

## The invariant on the soliton loop

| k | winding Q | sustain A² = (gρ−a₀)−(κ/R²)k² |
|---:|---:|---:|
| 1 | 1.0 | 0.1526 |
| 3 | 3.0 | 0.0815 |
| 5 | 5.0 | -0.0607 |

Loop radius `R=0.75` in the ordered core. The soliton sustains k=1,3 (A²>0); k=5 exceeds the loop's capacity (A²<0) and phase-slips under dissipation (the #182 event).

## Verdict

**DISCRETE_WINDING_INVARIANT_SURVIVES_CONTINUOUS_EVOLUTION_ON_THE_SOLITON_WHILE_Q_NONZERO.** SURVIVES — THE CONTINUOUS SOLITON CARRIES THE DISCRETE CHARGE. On the #180 self-consistent ψ–Φ–q throat-soliton, the winding invariant rides the continuous dynamics untouched.

QUANTIZED. On an ordered equatorial loop (R = 0.75, |q| > 0) the charge Q = (1/2π)∮∇φ = k exactly for k ∈ {1, 3, 5}; the soliton sustains k = 1, 3 (A² = {'1': 0.1526, '3': 0.0815, '5': -0.0607}).

SURVIVAL — WAVE. Continuous norm-conserving evolution conserves Q to machine precision for ALL k ∈ {1, 3, 5} (ΔW ~ 10⁻¹⁶, min|q| > 0).

SURVIVAL — DISSIPATIVE. The order field's own gradient flow preserves the sustained windings k = 1, 3 (ΔW ~ 10⁻¹⁵, the perturbed vortex relaxes back, min|q| > 0).

THE CRITERION (→ #182). Q changes ONLY through |q| = 0: the unsustained k = 5 (A² < 0) is driven to a zero (min|q| = 0.0001) and slips (5 → 2) — survival ⟺ |q| > 0, exactly; that slip is the PR #182 topology-change event.

RIGID. Under random |q|>0-preserving homotopies the charge is unchanged in every case ({'1': '40/40', '3': '40/40', '5': '40/40'}) — a superselection charge outside the continuous moduli (the #173/#174 rigidity on the dynamical soliton).

So the #174/#178 winding ladder rides the #179/#180 soliton untouched, except at the amplitude zeros where topology changes. The odd-k ladder's survival under a deformed bulk is PR #183.
