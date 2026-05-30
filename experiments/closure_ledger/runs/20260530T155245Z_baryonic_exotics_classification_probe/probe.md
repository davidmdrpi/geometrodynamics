# BAM-specific baryonic exotics: classification + experimental constraints (PR #102)

**Run:** 2026-05-30T15:52:45+00:00

Classifies BAM's non-orientable (Möbius / hybrid) BARYONIC exotics and ranks the channels by experimental constraint. **Key subtlety:** unlike mesons (where `1-+` is a smoking-gun exotic via `C`), baryons have **no forbidden `J^P`** — so BAM's Möbius/hybrid baryons are *supernumerary ordinary-`J^P`* states, identifiable only by counting. They land in the densely-measured light N*/Δ* region (~1.8–2.1 GeV), making them the **most** experimentally constrained corner of BAM's non-orientable predictions — the opposite extreme from the unobserved glueballs (PR #100).

- **Identification**: BAM-specific baryonic exotics (Möbius / hybrid baryon) have no exotic-J^P smoking gun (supernumerary ordinary-J^P states), land in the densely-measured light N*/Δ* region (~1.8–2.1 GeV, the 2√σ gap), and are the MOST experimentally constrained corner of BAM's non-orientable predictions
- **No smoking gun**: no exotic J^P for baryons (any J^P ordinary for qqq)
- **Masses**: hybrid N ≈ 1.79 GeV, hybrid Δ ≈ 2.08 GeV (nucleon/Δ + 2√σ)
- **Constraint ranking**: light N*/Δ* (strongest) > strange hyperons > charm/bottom (weakest)
- **Test**: Möbius doubling must coincide with observed states or decouple; else excluded
- **Contrast**: opposite extreme from the unobserved glueballs (PR #100)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_baryonic_exotics_setup` | baryonic exotics from non-orientable topology (Möbius/hybrid) | **PASS** |
| T2 | `T2_no_exotic_jp_for_baryons` | no exotic J^P for baryons ⟹ supernumerary, no smoking gun | **PASS** |
| T3 | `T3_baryonic_exotic_masses` | hybrid N ≈ 1.79 GeV, hybrid Δ ≈ 2.08 GeV (base + 2√σ) | **PASS** |
| T4 | `T4_experimental_constraint_ranking` | constraint ranking: light N*/Δ* > hyperons > heavy baryons | **PASS** |
| T5 | `T5_lightest_exotic_in_dense_region` | lightest exotic in the dense (most-constrained) N*/Δ* region | **PASS** |
| T6 | `T6_test_and_tension` | Möbius doubling: coincide or decouple, else excluded | **PASS** |
| T7 | `T7_opposite_extreme_from_glueballs` | opposite extreme from glueballs (most vs least constrained) | **PASS** |
| T8 | `T8_assessment` | BARYONIC_EXOTICS_LACK_EXOTIC_JP_LIGHT_SECTOR_MOST_CONSTRAINED | **PASS** |

## T2: Every baryon J^P is ordinary (no exotic)

| J^P | reachable by qqq? |
|---|:---:|
| 1/2+ | ✓ |
| 1/2− | ✓ |
| 3/2+ | ✓ |
| 3/2− | ✓ |
| 5/2+ | ✓ |
| 5/2− | ✓ |
| 7/2+ | ✓ |

Contrast with mesons: `1-+` is exotic (C-forbidden). No baryon `J^P` is forbidden, so a Möbius/hybrid baryon is a supernumerary ordinary-`J^P` state — countable, not manifestly exotic.

## T4: Experimental-constraint ranking

| channel | data density | constraint |
|---|---|---|
| light N*/Δ* (<2.5 GeV) | ~40 PDG states; missing-resonance problem | strongest |
| strange hyperons (Λ*, Σ*) | well measured | strong |
| charm/bottom baryons (Ξ_c*, Ω_c*, Λ_b*) | sparse (LHCb handful) | weakest (BAM freest) |

The lightest BAM baryonic exotic (`1.79 GeV` N-based, `2.08 GeV` Δ-based) sits in the densely-measured light sector — the strongest constraint.

## Verdict

**BARYONIC_EXOTICS_LACK_EXOTIC_JP_LIGHT_SECTOR_MOST_CONSTRAINED.** BAM-SPECIFIC BARYONIC EXOTICS HAVE NO SMOKING GUN AND LIVE IN THE MOST CONSTRAINED CHANNEL. PR #101 matched BAM's non-orientable topology to the OBSERVED mesonic exotics and flagged the Möbius baryon as a BAM-specific prediction. This probe classifies the baryonic exotics and identifies which channels are experimentally constrained.

NO EXOTIC J^P FOR BARYONS. For mesons the Möbius flux tube gave a smoking-gun exotic quantum number, J^PC = 1-+, forbidden to qq̄ by the C=(−1)^{L+S} constraint. Baryons have no such constraint: a qqq baryon has P=(−1)^L, S∈{½,3/2}, and no good C, so every half-integer J^P is reachable by an ordinary baryon — there is NO forbidden (exotic) J^P. So a BAM Möbius / hybrid baryon carries ORDINARY quantum numbers; it is a SUPERNUMERARY state, identifiable only by counting (an extra resonance beyond the qqq quark model), not by a smoking-gun J^P.

THE BAM BARYONIC EXOTICS. The hybrid baryon (excited Y-junction glue) sits ≈ nucleon + 2√σ ≈ 1.79 GeV; the Δ-based partner ≈ 2.08 GeV — the same 2√σ flux-tube quantum as the mesonic hybrids (PR #101). The Möbius baryon (make_mobius_baryon) is the Z₂-twisted partner. Both land in the LIGHT N*/Δ* region.

EXPERIMENTAL-CONSTRAINT RANKING. By data density: light N*/Δ* (< 2.5 GeV) — ~40 PDG states, the missing-resonance problem active — is the STRONGEST constraint; strange hyperons (Λ*, Σ*) strong; charm/bottom baryons (Ξ_c*, Ω_c*, Λ_b*) sparse — the WEAKEST (BAM freest). The lightest BAM baryonic exotics (~1.8–2.1 GeV) sit squarely in the MOST-constrained light sector (N(1710), N(1875), N(1900), …; Δ(1900), Δ(1950), Δ(2000), …).

THE TEST / TENSION. BAM's Möbius topology doubles the spectrum (a Z₂-twisted partner per state). The dense light sector cannot absorb arbitrary extras, so the Möbius/hybrid baryons must either COINCIDE with observed-but-unexplained resonances (filling missing resonances the quark model under-predicts) or DECOUPLE from the dominant πN production/decay (the standard missing-resonance mechanism); over-prediction without decoupling would be excluded. So this is the opposite extreme from the unobserved glueballs (PR #100): the light baryonic exotics are BAM's MOST exposed topological prediction.

HONEST SCOPE. ESTABLISHED (BAM-native): baryonic exotics have no exotic-J^P smoking gun (any J^P is ordinary for qqq), so BAM's Möbius/hybrid baryons are supernumerary ordinary-J^P states identifiable only by counting; they land in the light N*/Δ* region (~1.8–2.1 GeV, the 2√σ gap), the MOST experimentally constrained channel; the constraint ranking is light N*/Δ* > strange hyperons > heavy-quark baryons (the freest). NOT established: a specific Möbius-baryon state confirmed in data — without a smoking-gun J^P the prediction is a counting one, testable only by matching the dense spectrum or demonstrating decoupling.

## What this leaves open

- **A confirmed Möbius-baryon state** — without a smoking-gun `J^P`, the prediction is a counting one; testing it means matching the dense N*/Δ* spectrum or demonstrating decoupling (the missing-resonance mechanism).
- **Heavy-quark baryon exotics** — the least-constrained channel, where BAM's Möbius predictions are freest and a clean new state is most likely findable.
