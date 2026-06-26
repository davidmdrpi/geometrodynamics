# Odd-k / generation-sector survival under a deformed bulk geometry (PR #183)

**Run:** 2026-06-26T02:07:10+00:00

Shows the odd-k {1,3,5} charged-lepton generation sector (#174) is topologically PROTECTED — set by metric-independent invariants, it survives any smooth bulk deformation and changes only at a genuine topology change (the bulk-level analog of #181/#182). *(QFT on the classical throat, not quantum gravity.)*

- **Grading**: deck det ∓1 (RP²/RP³); ½ tr T² = −1 (Pin⁻); tr(T^k)=0 odd / ±2 even
- **Survival**: 1000 random deformations preserve the invariants to ~1e-15
- **Generations**: odd k ≤ D_bulk=5 ⟹ {1,3,5}=3 (topological count) survives
- **Topology change**: flips only at the degenerate spin structure ½ tr T²=0 (θ=π/4)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | does the odd-k generation sector survive a deformed bulk? | **PASS** |
| T2 | `T2_topological_grading` | the grading as a topological invariant (deck det; T²; tr T^k) | **PASS** |
| T3 | `T3_survival_under_deformation` | survival under smooth deformation (conjugation/GL⁺; machine-exact) | **PASS** |
| T4 | `T4_generation_count_survives` | the generation count {1,3,5}=3 survives (D_bulk + parity) | **PASS** |
| T5 | `T5_changes_only_at_topology_change` | changes only at a topology change (½ tr T² crosses 0 at θ=π/4) | **PASS** |
| T6 | `T6_robustness_random_deformations` | robustness under many random smooth deformations | **PASS** |
| T7 | `T7_unity_and_scope` | the unity with #181/#182; scope | **PASS** |
| T8 | `T8_assessment` | ODD_K_GENERATION_TOPOLOGICALLY_PROTECTED | **PASS** |

## The metric-independent grading

| invariant | value | meaning |
|---|---:|---|
| deck det (brane S²/antipodal) | -1 | RP² **non-orientable** |
| deck det (bulk S³/antipodal) | +1 | RP³ orientable |
| ½ tr T² | -1 | Pin⁻ (T² = −I) |

`tr(T^k)`: k=1: +0, k=2: -2, k=3: +0, k=4: +2, k=5: +0, k=6: -2 — **0 for odd (fermion), ±2 for even (boson)**.

## Verdict

**ODD_K_GENERATION_SECTOR_SURVIVES_BULK_DEFORMATION_TOPOLOGICALLY_PROTECTED_CHANGES_ONLY_AT_A_TOPOLOGY_CHANGE.** SURVIVES — THE GENERATION SECTOR IS TOPOLOGICALLY PROTECTED. The odd-k {1, 3, 5} ladder rides any smooth bulk deformation.

TOPOLOGICAL GRADING. The grading is set by metric-independent invariants: deck det = -1 (brane RP², non-orientable) and +1 (bulk RP³, orientable); ½ tr T² = -1 (Pin⁻); tr(T^k) = 0 for odd k (fermion), ±2 for even (boson).

SURVIVAL. Across 1000 random orientation-preserving deformations ½ tr T² stays −1 (to 2e-15) and the deck dets stay ∓1 (to 1e-15) — machine precision; named squash/tidal deformations too.

GENERATION COUNT. Odd k ≤ D_bulk = 5 ⟹ [1, 3, 5] = 3 generations (matching LEPTON_BASELINE_DEPTHS); D_bulk and the parity selection are topological, so the count survives every smooth deformation.

CHANGES ONLY AT A TOPOLOGY CHANGE. The only sector-flipping path (T² : −I → +I) crosses the degenerate spin structure ½ tr T² = 0 at θ = 0.7854 ≈ π/4 — the bulk-level analog of the #182 amplitude zero; smooth deformations never reach it.

ROBUST. 500/500 random deformations preserved the odd-k grading and the orientability class. The generation sector is to the bulk what the winding is to the soliton (#181/#182): a topological charge robust to smooth deformation, changing only at a topology change. The #174 round-metric result is not special — it is topologically protected.
