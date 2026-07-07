# The GHZ sector: multipartite entanglement is bridge valence - companion probe (PR #208)

**Run:** 2026-07-07T18:23:03+00:00

The deliverable is `docs/multi_mouth_bridge_ghz.md` - the three-mouth bridge and the GHZ sector, with the charged no-go. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the GHZ open: the trousers junction | **PASS** |
| T2 | `T2_charged_no_go` | charged GHZ superselection-forbidden (derived) | **PASS** |
| T3 | `T3_y_junction_identification` | the Y-junction live; leg-cut -> the #206 pair | **PASS** |
| T4 | `T4_ghz_emergence` | GHZ at fidelity 0.999+; holonomy law; 3-tangle 1 | **PASS** |
| T5 | `T5_mermin` | Mermin 4 vs local bound 2; pairwise marginals empty | **PASS** |
| T6 | `T6_valence_ledger` | bridge valence is the entanglement class (the ledger) | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | one open left; two successor questions | **PASS** |

## The multipartite holonomy orbit

| legs (s_B, s_C) | phase measured | predicted | fidelity |
|---|---:|---:|---:|
| (2,2) | 0.0 | 0.0 | 1.0 |
| (0,2) | 3.1416 | 3.1416 | 1.0 |
| (2,1) | 1.5708 | 1.5708 | 1.0 |

3-tangle = 1.0; Mermin = 4.0 (local bound 2); pairwise negativities [-0.0, -0.0, -0.0], pairwise CHSH [2.0, 2.0, 2.0].

## Verdict

**GHZ_IS_A_THREE_MOUTH_BRIDGE_MERMIN_4_MEASURED_CHARGED_GHZ_SUPERSELECTION_FORBIDDEN_PAIRWISE_MARGINALS_EMPTY_BRIDGE_VALENCE_IS_THE_ENTANGLEMENT_CLASS.** THE GHZ SECTOR IS BRIDGE VALENCE (the argument is in docs/multi_mouth_bridge_ghz.md; this probe measures every claim).

THE NO-GO. For a flux (charge) label the Y-junction has 0 conserving channels in the doublet, and GHZ components straddle charge sectors [-1, 1]: charged GHZ is superselection-forbidden, exactly as in QM. Multipartite entanglement belongs to the transported-frame (spin/Pin) label - a derived charge/spin split, invisible in the pairwise sector where (k, -k) is zero-sum.

THE JUNCTION. One bulk fiber, three reading frames: exactly symmetric distribution, per-leg deck phases as transport demands, leg-cut collapsing to the #206 pair with the composed holonomy (consistency across three PRs).

THE STATE. The three-frame embedding is an isometry; the extracted state is GHZ at fidelity 1.0 obeying the multipartite holonomy law phi = -pi(s_B + s_C)/2; 3-tangle = 1.0. Mermin: local bound 2 (64 strategies enumerated), the Y-state reaches 4.0 - the algebraic maximum - while every pairwise marginal is empty (negativities [-0.0, -0.0, -0.0], pairwise CHSH exactly 2). Valence 2: pairs carry everything (#207). Valence 3: the triple carries everything. Bridge valence is the entanglement class; one open remains (spatial/measurement), two successor questions named (the 5D pants; W-class reachability).
