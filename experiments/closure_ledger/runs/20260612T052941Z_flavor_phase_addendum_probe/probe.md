# Flavor phase addendum: Hopf CP derivation and full CKM realization (PR #162)

**Run:** 2026-06-12T05:29:41+00:00

The consolidating addendum for the CP-phase arc #156→#161: one keystone from every step re-verified in a single run, and the full account written into docs/THESIS.md — replacing the interim #158–#159 postscript and updating the #157 card's quark-CP row. *(QFT on the classical throat, not quantum gravity.)*

- **Arc**: #156 corrected → #158 relocation → #159/#160 φ_h = π/k₅ derived → #161 dataset realized
- **Keystones**: unitarity/quartet exact; k·π/k₅ exact; Weyl exact; nine observables ≤ 1%
- **Bookkeeping**: #149–#161: net inputs 0; knobs −1
- **Thesis**: postscript → full addendum; #157 card quark-CP row updated
- **Remaining**: knob-level v3+CP re-lock; lepton anarchic draw

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the addendum: consolidate #156→#161 into the thesis | **PASS** |
| T2 | `T2_relocation_keystone_158` | #158: unitarity/quartet exact; partition-mixing J ≈ 0 | **PASS** |
| T3 | `T3_transport_keystone_159` | #159: sector-arc transport k·π/k₅ exact | **PASS** |
| T4 | `T4_weyl_keystone_160` | #160: the Weyl commutator exact | **PASS** |
| T5 | `T5_realization_keystone_161` | #161: nine observables ≤ 1%, masses to 1e-14 (re-solved) | **PASS** |
| T6 | `T6_chain_and_bookkeeping` | chain sourced (3 derived + 1 locked); net inputs 0, knobs −1 | **PASS** |
| T7 | `T7_scope` | thesis edits; remaining items unchanged | **PASS** |
| T8 | `T8_assessment` | FLAVOR_PHASE_ADDENDUM_ALL_KEYSTONES_REVERIFIED | **PASS** |

## The derivation chain (consolidated)

| ingredient | value | status |
|---|---|---|
| rate ½ | A_φ(χ=0) — the spin-½ factor | derived (connection module) |
| sign ± | Z₂ partition orientation | derived (#63 C-swap) |
| winding dk | max(k, k′) | locked (the v3 mass calibration) |
| arc 2π/k₅ | the Weyl commutator quantum | derived (#160 algebra) |
| φ_h = π/k₅ | 0.6283 | derived (#159; alternatives excluded) |

## The #161 realization, re-solved

Residual 0.00519; masses preserved to 8e-15; predicted observables: V_td ×1.008, V_ts ×1.0021, J ×1.0044, α = 91.8106°, sin δ = 0.8888.

## Verdict

**FLAVOR_PHASE_ADDENDUM_ALL_KEYSTONES_REVERIFIED_THESIS_CONSOLIDATED.** EVERY STEP OF THE #156→#161 CP ARC RE-VERIFIES IN ONE RUN AND THE CONSOLIDATED ACCOUNT IS WRITTEN INTO THE THESIS. The arc: the #156 partition-mixing calibration was corrected away (#158: 16% non-unitarity, quartet-inconsistent J, unitarized core J ≈ 0, first-row unitarity ×40); the phase relocated to the Hopf-fiber transport of the same-partition shell couplings; its scale derived as φ_h = π/k₅ (#159: connection ½ × C-swap sign × mass-locked dk × sector arc, alternatives excluded; #160: the arc = the Weyl commutator quantum, algebra); and the soft V_us direction resolved through the mass-preserving family until the complete nine-observable flavor-CP dataset landed at ≤ 1% (#161).

THE KEYSTONES, RE-RUN. #158: exact unitarity and quartet consistency at the Hopf phases; the excluded route's unitarized J still ≈ 0. #159: sector-arc transport k·π/k₅ exact. #160: the Weyl commutator exact. #161: re-solved — residual 0.00519, masses to 8e-15, the four predicted observables landing (V_td ×1.01, V_ts ×1.00, J ×1.00, α = 91.8°, sin δ = 0.889).

THE BOOKKEEPING. Thirteen probes #149–#161: net ZERO inputs (the #156 input consumed, then returned by the #159 derivation) and one modelling knob retired. THE THESIS now carries the full account — the derivation chain, the realized dataset, the re-lock targets — with the interim postscript replaced and the #157 card's quark-CP row updated. Remaining: the knob-level v3+CP re-lock and the lepton anarchic draw.
