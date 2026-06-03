# The combined matter-sector APS ledger (PR #125)

PRs #123 and #124 ran the Atiyah–Patodi–Singer index audit of the factorized
sector sum (PR #122) on the quark and lepton sectors separately. This PR
**combines** them into one matter-sector APS ledger and ties it to the
program's full input budget (PRs #104–#108, #112): the same index machinery,
applied uniformly, classifies every matter sector's partition into a
**derived topological part** and (at most) **one residual feeding integer** —
certifying that the matter sectors contribute exactly **one** dimensionless
partition residual (the quark `n_part`), with **leptons fully derived**.

## The universal APS structure

Every matter-sector partition has the form

```
N_sector = (structural factor) × (feeding integer),
```

and the Z₂-graded Witten/APS index of the factorized sum is universal: the
spectral flow is the integer `1`, the APS ξ-invariant is `ξ(a) = (η+h)/2 =
1/2 − a` (PRs #119–#124). The **topological content** — the structural factor
and the integer spectral flow — is derived in every sector; only the
**feeding integer** can be a residual.

## The matter-sector partition ledger

| sector | partition | feeding integer | residual |
|---|---|---|---|
| lepton | `N_lepton = 4·k₅² = 100` | `k₅ = 5` (DERIVED bulk dim, #73) | **none** |
| quark | `N_q = 2·n_part = 466` | `n_part = 233` (RESIDUAL, drifts 216–255) | **n_part** |
| neutrino | (ε compliance / healing length) | `ε ~ R_c³` (order-of-mag DERIVED, #112) | **ε value** |

So the unique matter-**partition** residual is the quark `n_part`: leptons
contribute 0, quarks exactly 1. (The neutrino's `ε` is a compliance/healing
length, not a closure-partition count; its order of magnitude is derived,
its value residual, #112.)

## Tie to the full input budget

Combining with PRs #104–#108 and #112, the BAM input budget is:

  - **One dimensionful anchor: `G`** — the bulk-gravity scale, from which
    both `m_e` and `√σ` descend (#105/#106); `ℏ` geometric, `c` units.
  - **Four dimensionless residuals:**
      - `n_part` (quark partition — the APS-confirmed lone matter-partition
        residual, #123);
      - `√σ/m_e ≈ 830` (lepton/QCD ratio — irreducible, #108);
      - `ε` (neutrino compliance value — order-of-mag derived, #112);
      - `α` (universal coupling — #105).
  - **The universal flavor puzzle** (Yukawa hierarchy — not BAM-specific).

The APS audit's specific contribution: of the matter-sector **partition
counts**, exactly one (`n_part`) is residual — leptons are fully derived
from `k₅` (the bulk dimension). The other dimensionless residuals (`√σ/m_e`,
`ε`, `α`) are a cross-sector ratio, a compliance, and a coupling — not
partition counts — so APS isolates `n_part` as the unique matter-partition
residual.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | combine #123/#124 into a matter-sector ledger + input budget |
| T2 | universal structure | `N = factor × feeding integer`; spectral flow = 1; `ξ(a) = 1/2 − a` |
| T3 | partition ledger | lepton/quark/neutrino |
| T4 | residual count | leptons 0, quarks 1 (`n_part`) |
| T5 | input budget | 1 anchor `G` + 4 residuals + flavor puzzle |
| T6 | APS sharpening | isolates `n_part` as the unique matter-partition residual |
| T7 | scope | classification established; residuals not removed |
| T8 | assessment | `COMBINED_MATTER_SECTOR_APS_LEDGER_LEPTON_DERIVED_QUARK_ONE_RESIDUAL` |

## Established and open

  - **Established (BAM-native):** a uniform topological classification of the
    matter sectors — partition = (derived structural factor) × (feeding
    integer), with the topology (factor + spectral flow) derived everywhere;
    the matter sectors carry exactly **one** partition residual (`n_part`);
    leptons fully derived; the full input budget assembled (one anchor `G` +
    four dimensionless residuals + the flavor puzzle).

  - **Does not / open:** the residuals (`n_part`, `√σ/m_e`, `ε`, `α`) stand —
    the APS index organizes and isolates them, it does not remove them.

## Cross-references

  - `docs/aps_quark_partition_index_research_plan.md` /
    `docs/aps_lepton_partition_index_research_plan.md` — PR #123/#124, the
    sector audits combined here.
  - `docs/program_synthesis_research_plan.md` — PR #104, the 5-tier epistemic
    ledger / input budget.
  - `docs/scale_count_anchors_research_plan.md` /
    `docs/lepton_qcd_ratio_legitimate_search_research_plan.md` — PR #106/#108,
    the one-scale-`G` reduction and the `√σ/m_e` residual.
  - `docs/epsilon_bulk_compliance_research_plan.md` — PR #112, the neutrino
    `ε` (order-of-mag derived, value residual).

## Run

```
python -m experiments.closure_ledger.combined_matter_sector_aps_ledger_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_combined_matter_sector_aps_ledger_probe/`.
Expected verdict:
`COMBINED_MATTER_SECTOR_APS_LEDGER_LEPTON_DERIVED_QUARK_ONE_RESIDUAL`, 8/8 PASS.
