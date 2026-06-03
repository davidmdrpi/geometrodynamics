# The combined matter-sector APS ledger (PR #125)

**Run:** 2026-06-03T23:10:21+00:00

Combines the quark (#123) and lepton (#124) APS index audits into one matter-sector ledger and ties it to the program's full input budget (PRs #104–#108, #112). The uniform result: every matter partition is **(derived topological factor) × (feeding integer)**, so the matter sectors carry **exactly one partition residual** (the quark `n_part`), with **leptons fully derived**.

- **Universal structure**: N = (derived factor) × (feeding integer); spectral flow = 1; ξ(a) = 1/2 − a
- **Partition ledger**: lepton 4·k₅² (k₅ derived, no residual); quark 2·n_part (n_part residual); neutrino ε (value residual)
- **Unique partition residual**: n_part (leptons 0, quarks 1)
- **Input budget**: 1 anchor G + 4 residuals {n_part, √σ/m_e, ε, α} + flavor puzzle
- **APS role**: organizes and isolates the residuals; does not remove them

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | combine #123/#124 into a matter-sector ledger + input budget | **PASS** |
| T2 | `T2_universal_aps_structure` | universal: N = factor × feeding integer; spectral flow = 1 | **PASS** |
| T3 | `T3_matter_sector_partition_ledger` | partition ledger (lepton/quark/neutrino) | **PASS** |
| T4 | `T4_matter_partition_residual_count` | residual count: leptons 0, quarks 1 (n_part) | **PASS** |
| T5 | `T5_full_input_budget` | input budget: 1 anchor G + 4 residuals + flavor puzzle | **PASS** |
| T6 | `T6_aps_sharpening` | APS isolates n_part as the unique matter-partition residual | **PASS** |
| T7 | `T7_scope` | scope: classification established; residuals not removed | **PASS** |
| T8 | `T8_assessment` | COMBINED_MATTER_SECTOR_APS_LEDGER_LEPTON_DERIVED_QUARK_ONE_RESIDUAL | **PASS** |

## The matter-sector partition ledger

| sector | partition | feeding integer | residual |
|---|---|---|---|
| lepton | `N_lepton = 4·k₅² = 100` | k₅ = 5 (DERIVED bulk dim, #73) | none |
| quark | `N_q = 2·n_part = 466` | n_part = 233 (RESIDUAL, drifts 216–255) | n_part |
| neutrino | `(ε compliance / healing length)` | ε ~ R_c³ (order-of-mag DERIVED, #112) | ε value |

The topology (the structural factor `4`/`2` and the integer spectral flow `1`) is **derived in every sector**; only the feeding integer can be residual. Leptons: `k₅` derived ⟹ no residual. Quarks: `n_part` the lone partition residual.

## The full input budget

- **One dimensionful anchor:** `G` (bulk-gravity scale; `m_e`, `√σ` descend, #105/#106; `ℏ` geometric, `c` units).
- **Four dimensionless residuals:** `n_part` (quark partition, APS-confirmed lone matter-partition residual), `√σ/m_e ≈ 830` (lepton/QCD ratio, irreducible #108), `ε` (neutrino compliance value, order-of-mag derived #112), `α` (universal coupling #105).
- **The universal flavor puzzle** (Yukawa hierarchy — not BAM-specific).

APS isolates `n_part` as the unique matter-**partition** residual; the other three are a ratio, a compliance, and a coupling.

## Verdict

**COMBINED_MATTER_SECTOR_APS_LEDGER_LEPTON_DERIVED_QUARK_ONE_RESIDUAL.** THE COMBINED MATTER-SECTOR APS LEDGER: EVERY PARTITION IS (DERIVED TOPOLOGICAL FACTOR) × (FEEDING INTEGER), SO THE MATTER SECTORS CARRY EXACTLY ONE PARTITION RESIDUAL — THE QUARK n_part — WITH LEPTONS FULLY DERIVED. PRs #123/#124 ran the APS index on the quark and lepton sectors; this probe combines them and ties the result to the input budget.

THE UNIVERSAL APS STRUCTURE. Every matter-sector partition is N_sector = (structural factor) × (feeding integer), and the Z₂-graded Witten/APS index of the factorized sum (PR #122) is universal: the spectral flow is the integer 1, the APS ξ-invariant is ξ(a) = (η+h)/2 = 1/2 − a. The topological content — the structural factor and the integer spectral flow — is derived in every sector; only the feeding integer can be a residual.

THE PARTITION LEDGER. Lepton: N_lepton = 4·k₅² = 100, feeding k₅ = 5 (the derived bulk dimension, #73) — NO residual. Quark: N_q = 2·n_part = 466, feeding n_part = 233 (the compensator residual, drifts 216–255) — residual n_part. Neutrino: ε (compliance/healing length), order-of-magnitude derived (~R_c³, #112) — residual the ε value. So the unique matter-PARTITION residual is the quark n_part: leptons contribute 0, quarks exactly 1.

THE FULL INPUT BUDGET. Combining with PRs #104–#108 and #112: one dimensionful anchor G (the bulk-gravity scale, from which both m_e and √σ descend, #105/#106; ℏ geometric, c units); four dimensionless residuals — n_part (quark partition, the APS-confirmed lone matter-partition residual), √σ/m_e ≈ 830 (the lepton/QCD ratio, irreducible, #108), ε (the neutrino compliance value, order-of-mag derived, #112), α (the universal coupling, #105); and the universal flavor puzzle. The APS audit's specific contribution is to isolate n_part as the unique matter-PARTITION residual — the other three residuals are a ratio, a compliance, and a coupling, not partition counts.

SCOPE. The ledger ESTABLISHES the uniform topological classification (partition = derived factor × feeding integer; topology derived everywhere; exactly one partition residual; leptons fully derived) and assembles the input budget. It does NOT derive n_part, √σ/m_e, ε, or α — APS organizes and isolates the residuals, it does not remove them.
