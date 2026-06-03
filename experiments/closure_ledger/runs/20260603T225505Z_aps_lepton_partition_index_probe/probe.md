# The APS lepton partition index from the factorized sector sum (PR #124)

**Run:** 2026-06-03T22:55:05+00:00

Runs the same APS audit as PR #123, now on the **lepton** sector. Because the lepton partition descends entirely from **derived geometry** (`k₅ = 5` = bulk dimension), the index fixes **both the structure AND the value** — no residual, in contrast to the quark `n_part`.

- **Lepton partition**: N_lepton = 4·k₅² = 100 (k₅ = 5 = bulk dim, derived; β_lepton = 50π)
- **APS index**: I = Tr(−1)^k topological; ξ(a) = 1/2 − a; spectral flow = 1 (universal)
- **Fully topological**: k₅ derived ⟹ structure AND value fixed, NO residual
- **Contrast**: quark N_q = 2·n_part (n_part residual, drifts) vs lepton N_lepton = 4·k₅² (k₅ derived, fixed)
- **Conclusion**: leptons the clean APS case; quark n_part the program's lone compensator residual

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | apply the APS audit (PR #123) to the lepton sector | **PASS** |
| T2 | `T2_lepton_partition_derived` | N_lepton = 4·k₅² = 100 (k₅ = bulk dim derived); 3 generations | **PASS** |
| T3 | `T3_aps_index_applies_identically` | APS index identical: ξ(a) = 1/2 − a, spectral flow = 1 | **PASS** |
| T4 | `T4_lepton_partition_fully_topological` | fully topological: k₅ derived ⟹ structure AND value fixed | **PASS** |
| T5 | `T5_contrast_with_quarks` | contrast: quark N_q = 2·n_part (residual) vs lepton 4·k₅² (derived) | **PASS** |
| T6 | `T6_s8_stability` | §8: N_lepton = 100 fixed; n_part drifts 216–255 | **PASS** |
| T7 | `T7_scope` | scope: lepton index fully determined; quark n_part lone residual | **PASS** |
| T8 | `T8_assessment` | APS_LEPTON_PARTITION_INDEX_FULLY_TOPOLOGICAL_DERIVED_NO_RESIDUAL | **PASS** |

## Quark vs lepton: the same index, opposite outcomes

| sector | partition | feeding integer | APS index fixes | residual |
|---|---|---|---|---|
| quark | `N_q = 2·n_part = 466` | `n_part` (drifts 216–255) | structure only (doubling) | **n_part value** |
| lepton | `N_lepton = 4·k₅² = 100` | `k₅ = 5` (bulk dim, derived) | structure **AND** value | **none** |

The universal spectral-flow index is `1` and `ξ(a) = 1/2 − a` for both. The difference is entirely in *what feeds the partition*: a derived integer (`k₅`) for leptons, a compensator residual (`n_part`) for quarks.

## Verdict

**APS_LEPTON_PARTITION_INDEX_FULLY_TOPOLOGICAL_DERIVED_NO_RESIDUAL.** THE APS LEPTON PARTITION INDEX IS FULLY DETERMINED — STRUCTURE AND VALUE — BY DERIVED GEOMETRY; THERE IS NO RESIDUAL, IN CONTRAST TO THE QUARK n_part. PR #123 ran the APS audit on the quark sector; this probe runs the same audit on the leptons.

THE LEPTON PARTITION. N_lepton = 4·k₅² = 100, with k₅ = 5 the bulk dimension dim(S³)+2 (DERIVED, PR #73) and β_lepton = k₅²·2π = 50π (PR #71). The factor of 4 is the closure structure 4β_lepton/(2π) = 4k₅². The three generations are the odd-k closures k ∈ {1,3,5}, with #gen = (k₅+1)/2 = 3.

THE APS INDEX APPLIES IDENTICALLY. The same Z₂-graded Witten/APS index of PR #123: the orientation sign (−1)^k makes I = Tr(−1)^k topological; the APS ξ-invariant is ξ(a) = (η+h)/2 = 1/2 − a; and the spectral flow over one holonomy cycle is ξ(0⁺) − ξ(1⁻) = 1, an integer — universal, the same for the lepton and quark sectors.

THE LEPTON PARTITION IS FULLY TOPOLOGICAL. Here is the contrast with the quark sector. For the quarks, N_q = 2·n_part and the feeding integer n_part is the phenomenological compensator residual (it drifts 216–255 across the quark_axioms §8 ablations), so the APS index fixes only the topological doubling (N_q even). For the leptons, the feeding integer k₅ = 5 is a fixed DERIVED structural integer (the bulk dimension), so there is no §8 ablation that moves it: N_lepton = 4·k₅² = 100 is fixed in BOTH its structure (the 4k₅² closure form) AND its value. The lepton partition index is therefore FULLY topological — there is no residual.

THE CONCLUSION. Leptons are the clean APS case; the quark n_part is the program's lone compensator residual. The same index machinery, applied to both, sharpens the program ledger: the lepton partition is fully determined by derived geometry (k₅, the bulk dimension), while the quark partition carries the one undetermined dimensionless residual (n_part). This is exactly why the lepton sector has been the clean, predictive one throughout.
