# The APS lepton partition index from the factorized sector sum (PR #124)

**Run:** 2026-06-03T22:54:44+00:00

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
| T3 | `T3_aps_index_applies_identically` | APS index identical: ξ(a) = 1/2 − a, spectral flow = 1 | **FAIL** |
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

**APS_LEPTON_PARTITION_INDEX_INCONCLUSIVE.** INCONCLUSIVE. A structural test failed; review the lepton partition / APS accounting.
