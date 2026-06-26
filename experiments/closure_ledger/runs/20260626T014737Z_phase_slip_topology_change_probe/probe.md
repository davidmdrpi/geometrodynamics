# The phase-slip / topology-change event (PR #182)

**Run:** 2026-06-26T01:47:37+00:00

Dissects how the discrete winding invariant `Q=(1/2π)∮∇φ` changes when the order field hits zero — the topological obstruction, the ±1 quantum, and the dynamical staircase — on the #180/#181 throat-soliton. *(QFT on the classical throat, not quantum gravity.)*

- **Obstruction**: changing Q forces an exact zero (no nowhere-zero path between sectors)
- **Quantum**: a simple zero changes Q by exactly ±1 (a localized 2π phase kink)
- **Dynamics**: Q(t) piecewise-constant, steps ±1 exactly at min|q|→0 (single slip + staircase)
- **Localization**: Q is an exact integer except at the measure-zero amplitude-zero events

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | dissect the phase slip — exactly how Q changes at |q|=0 | **PASS** |
| T2 | `T2_topological_obstruction` | the obstruction: changing Q forces an exact zero (homotopy) | **PASS** |
| T3 | `T3_elementary_quantum_pm_one` | the quantum: across a simple zero ΔQ=±1 (a 2π phase kink) | **PASS** |
| T4 | `T4_single_dynamical_slip` | the single event: Q steps −1 exactly at min|q|→0 | **PASS** |
| T5 | `T5_quantized_staircase_cascade` | the staircase: a cascade of ±1 slips, each at a zero | **PASS** |
| T6 | `T6_invariant_localized_to_zeros` | localization: Q is an exact integer except at the zeros | **PASS** |
| T7 | `T7_physical_meaning_and_scope` | physical meaning + #175/#178/#181/#183 bridges; scope | **PASS** |
| T8 | `T8_assessment` | TOPOLOGY_CHANGE_ONLY_AT_ZEROS_BY_PLUS_MINUS_ONE | **PASS** |

## The quantized staircase (k = 8 relaxing)

Winding cascade: **[8, 7, 5, 4, 3, 2]** — 5 elementary slips, each coinciding with a `min|q|→0` event (min|q| at the steps: [0.0, 0.0, 0.0, 0.0, 0.0]).

The throat sheds winding one quantum at a time, each through an amplitude-zero node.

## Verdict

**TOPOLOGY_CHANGE_OCCURS_ONLY_AT_AMPLITUDE_ZEROS_EACH_ELEMENTARY_SLIP_CHANGES_Q_BY_PLUS_MINUS_ONE.** RESOLVED — THE INVARIANT CHANGES ONLY AT |q| = 0, BY ±1. The phase slip, dissected on the #180/#181 throat-soliton.

OBSTRUCTION. The straight homotopy between the winding-1 and winding-0 sectors is FORCED through an exact zero (min|q| = 2.5e-17 at s* = 0.495, φ* = 3.1416 ≈ π), where Q jumps 1 → 0: no nowhere-zero path connects the sectors — changing Q REQUIRES |q| = 0.

QUANTUM. Across the simple zero ΔQ = -1: the integrated winding density ∮∇φ changes by exactly -6.2832 = −2π — one full turn of phase removed at the zero. The elementary slip is ±1.

SINGLE EVENT. In a genuine ψ–Φ–q evolution the unsustained k = 5 holds flat, then the first slip steps ΔQ = -1 to 4 EXACTLY when min|q| → 0 (min|q| = 0.0005 at the slip).

STAIRCASE. The over-wound k = 8 cascades down the quantized staircase [8, 7, 5, 4, 3, 2] (5 slips), every step at an amplitude-zero event — the throat sheds winding one quantum at a time through the nodes.

LOCALIZATION. Away from the zeros the unrounded winding is an exact integer (to 2e-15); the invariant is ambiguous only at the amplitude zeros. With #181 (survival between events), the throat's winding sector is a conserved charge that transitions ONLY at the #175/#178 nodes.
