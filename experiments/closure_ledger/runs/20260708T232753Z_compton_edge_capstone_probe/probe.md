# The Compton-edge capstone: the second release - companion probe (PR #211)

**Run:** 2026-07-08T23:27:53+00:00

The deliverable is `docs/compton_edge_capstone.md` - the #204-#210 arc consolidated, with the Compton-edge law. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the closing PR: consolidation + one new law | **PASS** |
| T2 | `T2_compton_edge_law` | the Compton-edge law: the O(1) confined (new) | **PASS** |
| T3 | `T3_ledger_commitment_chain` | ledger 1 green: the commitment chain (#204/#205) | **PASS** |
| T4 | `T4_ledger_topology_chain` | ledger 2 green: the topology chain (#206-#208) | **PASS** |
| T5 | `T5_ledger_measurement_mass` | ledger 3 green: measurement + mass (#209/#210) | **PASS** |
| T6 | `T6_claim_graph_and_register` | the claim graph + the updated register | **PASS** |
| T7 | `T7_honest_scope` | honest scope (the arc's stated conditions) | **PASS** |
| T8 | `T8_assessment` | Release II | **PASS** |

## The Compton-edge law

| x = sigma/lambda_C | S(x) | D(x) |
|---:|---:|---:|
| 0.25 | 1.0105 | 0.9947 |
| 0.5 | 1.0415 | 0.9793 |
| 0.647 | 1.0692 | 0.9657 |
| 1.0 | 1.1634 | 0.9207 |
| 1.509 | 1.3628 | 0.8306 |
| 2.0 | 1.6176 | 0.7256 |
| 3.0 | 2.2824 | 0.501 |
| 4.0 | 3.0786 | 0.3114 |
| 5.0 | 3.9508 | 0.1785 |
| 6.0 | 4.8665 | 0.0962 |
| 8.0 | 6.765 | 0.0244 |

(natural windows: x <= 2.58 at S <= 2, x <= 5.58 at S <= 4.48; corrected self-consistent band [0.626, 1.307])

## Verdict

**RELEASE_II_GREEN_THE_ENTANGLED_SECTOR_CLOSED_OPERATIONALLY_THE_MASS_LADDER_ON_THE_ALPHA_ANCHOR_THE_O1_CONFINED_TO_THE_COMPTON_EDGE_OF_THE_NATURAL_DOMAIN.** RELEASE II - THE ARC CLOSES GREEN (the argument is in docs/compton_edge_capstone.md).

THE NEW LAW. The massive bridge equation deforms #202's exact suppression law by a UNIVERSAL function of sigma_mode/lambda_C alone (m r_s collapse 3e-04; massless limit exact to 0e+00): the sensitivity leaves 1 at the Compton scale and naturalness confines the #210 O(1) to x <= 2.58 (S <= 2) or x <= 5.58 (the ladder's own 4.48), with both convention anchors inside and the self-consistent band tightened to [0.626, 1.307] - still bracketing 1.

THE LEDGER. Commitment chain: retarded-real norm 6e-14, continuity 2e-03; classical channel entropy -0e+00 vs pairwise 0.0256; SN/BMV nulls armed. Topology chain: singlet fidelity 1.0, swapping fidelity 1.0, LHV bounds CHSH = 2 / Mermin = 2 enumerated, GHZ Mermin = 4.0 with empty pairs. Measurement + mass: the k-odd identity 2e-16, live Born within 0.0024, live Kaup 0.6318, the alpha bands [0.647, 1.509] and [0.00313, 0.0073] bracketing their targets.

THE REGISTER. Six items, all named and bounded; the falsification card armed with two near-term nulls and a x1.5-level mass prediction. The program's deepest machinery is now theorems, measurements, and windows - not imports.
