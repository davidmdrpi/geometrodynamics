# The APS lepton partition index from the factorized sector sum (PR #124)

PR #123 applied the Witten / Atiyah–Patodi–Singer index of the factorized
sector sum (PR #122) to the **quark** sector, finding that the index fixes
the §8-stable topological structure (the doubling `N_q = 2·n_part`) but
**not** the bare value `n_part` (the compensator, which drifts 216–255).
This PR runs the **same** APS audit on the **lepton** sector — and finds the
opposite: because the lepton partition descends entirely from **derived
geometry**, the index fixes **both the structure AND the value**, with no
residual.

## The lepton partition (derived geometry)

```
N_lepton = 4·k₅² = 100,      β_lepton = k₅²·(2π) = 50π,
```

with `k₅ = 5` the bulk dimension `dim(S³) + 2` (DERIVED, PR #73) and
`β_lepton` the lepton closure-winding mass parameter (DERIVED, PR #71). The
factor of 4 is the closure structure `4β_lepton/(2π) = 4k₅²`. The three
generations are the odd-k closures `k ∈ {1, 3, 5}`, with `#gen = (k₅+1)/2 =
3` (PR #73).

## The APS index applies identically

The same Z₂-graded Witten/APS index of PR #123: the orientation sign
`(−1)^k` makes `I = Tr(−1)^k` topological; the APS ξ-invariant is `ξ(a) =
(η+h)/2 = 1/2 − a`; and the spectral flow over one holonomy cycle is
`ξ(0⁺) − ξ(1⁻) = 1`, an integer — **universal**, the same for the lepton and
quark sectors.

## The lepton partition is FULLY topological — no residual

| sector | partition | feeding integer | APS index fixes | residual |
|---|---|---|---|---|
| quark | `N_q = 2·n_part = 466` | `n_part` (RESIDUAL, drifts 216–255) | structure only (doubling) | **n_part value** |
| lepton | `N_lepton = 4·k₅² = 100` | `k₅ = 5` (DERIVED bulk dim) | structure **AND** value | **none** |

For the quarks the index fixes only the topological doubling (`N_q` even),
because the feeding integer `n_part` is the compensator residual. For the
leptons the feeding integer `k₅ = 5` is a fixed **derived** structural
integer (the bulk dimension), so there is no §8 ablation that moves it:
`N_lepton = 4·k₅² = 100` is fixed in **both** its structure (the `4k₅²`
closure form) and its value. The lepton partition index is therefore
**fully topological** — no residual.

## The conclusion

Leptons are the clean APS case; the quark `n_part` is the program's lone
compensator residual. The same index machinery, applied to both, sharpens
the program ledger: the lepton partition is fully determined by derived
geometry (`k₅`, the bulk dimension), while the quark partition carries the
one undetermined dimensionless residual (`n_part`). This is exactly why the
lepton sector has been the clean, predictive one throughout.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | apply the APS audit (PR #123) to the lepton sector |
| T2 | lepton partition | `N_lepton = 4·k₅² = 100` (k₅ = bulk dim derived); 3 generations |
| T3 | APS index | identical: `ξ(a) = 1/2 − a`, spectral flow = 1 (universal) |
| T4 | fully topological | `k₅` derived ⟹ structure AND value fixed, no residual |
| T5 | contrast | quark `N_q = 2·n_part` (residual) vs lepton `4·k₅²` (derived) |
| T6 | §8 stability | `N_lepton = 100` fixed; `n_part` drifts 216–255 |
| T7 | scope | lepton index fully determined; quark `n_part` lone residual |
| T8 | assessment | `APS_LEPTON_PARTITION_INDEX_FULLY_TOPOLOGICAL_DERIVED_NO_RESIDUAL` |

## Established and open

  - **Established (BAM-native):** the same APS audit, applied to the lepton
    sector, finds the lepton partition `N_lepton = 4·k₅² = 100` fully
    determined — structure **and** value — because `k₅ = 5` is the derived
    bulk dimension (PR #73), not a compensator. The universal spectral-flow
    index (`= 1`) and `ξ(a) = 1/2 − a` are as in PR #123. The lepton
    partition has **no residual**: leptons are the clean APS case.

  - **Open / residual (quark sector, unchanged):** the quark `N_q = 2·n_part`
    fixes only the structure; the value `n_part` remains the §8-drifting
    compensator residual (PR #97/#107/#123) — the program's lone undetermined
    dimensionless integer.

## Cross-references

  - `docs/aps_quark_partition_index_research_plan.md` — PR #123, the quark
    APS audit this contrasts with.
  - `docs/bam_factorized_sector_sum_research_plan.md` — PR #122, the
    factorized sector sum the index is read off.
  - `docs/k5_origin_research_plan.md` / `docs/beta_lepton_derivation_research_plan.md`
    — PR #73/#71, the derived `k₅` and `β_lepton`.

## Run

```
python -m experiments.closure_ledger.aps_lepton_partition_index_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_aps_lepton_partition_index_probe/`.
Expected verdict:
`APS_LEPTON_PARTITION_INDEX_FULLY_TOPOLOGICAL_DERIVED_NO_RESIDUAL`, 8/8 PASS.
