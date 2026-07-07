# Entanglement swapping is bridge surgery - companion probe (PR #207)

**Run:** 2026-07-07T17:12:04+00:00

The deliverable is `docs/bridge_surgery_entanglement_swapping.md` - local pair annihilation as the bridge surgery that links never-co-nucleated throats. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the dynamical path: annihilation relinks bridges | **PASS** |
| T2 | `T2_composition_is_swapping` | the composition law == the QM swapping law (a+b+c) | **PASS** |
| T3 | `T3_composite_bridge_orbit` | the composite bridge: purity 0.9998; the 4-outcome orbit | **PASS** |
| T4 | `T4_annihilation_event` | the event: local pinch, distant linkage, radiated middle | **PASS** |
| T5 | `T5_swapped_bell_pair` | the swapped pair saturates Tsirelson; repeater composes | **PASS** |
| T6 | `T6_monogamy_is_the_matching` | monogamy is the matching (1/2 vs 0, machine) | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | the register moves; two named opens | **PASS** |

## The four-outcome holonomy orbit

| sj | phase measured | predicted | fidelity |
|---:|---:|---:|---:|
| 0 | 0.0 | 0.0 | 1.0 |
| 1 | 1.5874 | 1.5708 | 0.99991 |
| 2 | 3.165 | 3.1416 | 0.99983 |
| 3 | 4.7288 | 4.7124 | 0.99992 |

## The annihilation event

- (1,4) response to Alice's phase: 3e-06 before -> 1e-02 after (x4652); swapped-channel power ratio 288
- middle population: 0.5017 -> 0.2853 (control 0.5431); radiated fraction 0.0 -> 0.0448

## Verdict

**ENTANGLEMENT_SWAPPING_IS_BRIDGE_SURGERY_THE_JUNCTION_HOLONOMY_IS_THE_BELL_OUTCOME_LOCAL_ANNIHILATION_LINKS_NEVER_CONUCLEATED_THROATS_MONOGAMY_IS_THE_MATCHING.** THE DYNAMICAL PATH IS OPEN - LOCAL ANNIHILATION LINKS NEVER-CO-NUCLEATED THROATS (the argument is in docs/bridge_surgery_entanglement_swapping.md; this probe measures both sides).

THE LAW. Bridge surgery composes transports, and the composition law IS the quantum swapping law: phi_14 = phi_a + phi_b + phi_c, with the Bell outcome phi_c supplied by the junction holonomy of the pinch (QM side machine-checked at fidelity 1.0, probabilities 1/4; bulk side measured across the four-outcome orbit to 0.0234 rad, fidelities >= 0.99983). Winding superselection confines outcomes to the charge-zero sector, and the unconditioned mixture is separable (negativity 2e-17) - no outcome record, no usable correlation: swapping's classical-communication requirement and no-signaling in one identity.

THE EVENT. Two populated disjoint bridges; the proximity pair (2,3) annihilates - locally. The distant (1,4) response rises x4652, lands in the swapped channel (power ratio 288), carries Alice's phase linearly, and the middle population (0.5017 -> 0.2853, control 0.5431) radiates to wave fronts - the pair returns to radiation, the linkage persists.

THE NETWORK. The swapped pair saturates Tsirelson (worst CHSH 2.8283); monogamy is the perfect matching (partner negativity 1/2, non-partner exactly 0); surgeries compose associatively (the repeater chain) - pairwise entanglement distribution over the whole throat network is mechanism-complete. Remaining, named: the spatial/measurement sector; multipartite (>= 3-mouth) junctions.
