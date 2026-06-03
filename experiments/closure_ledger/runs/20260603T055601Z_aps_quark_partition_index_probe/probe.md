# The APS quark partition index from the factorized sector sum (PR #123)

**Run:** 2026-06-03T05:56:01+00:00

Reads the Witten / Atiyah–Patodi–Singer **index** off the factorized sector sum (PR #122). The Z₂-graded sum has a topological index, computed via the APS η-invariant; for the quark sector it fixes the **§8-stable topological structure** (the graded doubling `N_q = 2·n_part`) — but **not** the residual value of `n_part`.

- **Index**: I = Tr(−1)^k (topological, β-independent); spectral flow = 1 (integer)
- **APS ξ**: ξ(a) = (η+h)/2 = 1/2 − a = ζ_H(0,a) (η-boundary correction)
- **Quark partition**: N_q = 2·n_part = 466 (the Z₂-graded doubling, §8-stable)
- **Topological vs residual**: §8-stable: doubling N_q=2·n_part (even) + integer index; residual: n_part value (drifts 216–255, PR #97/#107)
- **Derives**: the topological structure (graded doubling, spectral-flow integer), NOT the n_part value

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | derive APS quark partition index from the factorized sum | **PASS** |
| T2 | `T2_witten_aps_index_from_z2_grading` | index I = Tr(−1)^k from the Z₂ grading (topological) | **PASS** |
| T3 | `T3_aps_xi_invariant` | APS ξ-invariant ξ(a) = (η+h)/2 = 1/2 − a | **PASS** |
| T4 | `T4_spectral_flow_integer_index` | integer index = spectral flow = ξ(0⁺)−ξ(1⁻) = 1 | **PASS** |
| T5 | `T5_apply_to_quark_partition` | quark partition N_q = 2·n_part (the Z₂-graded doubling) | **PASS** |
| T6 | `T6_topological_stable_vs_value_residual` | topological (doubling/index) §8-stable vs n_part value residual | **PASS** |
| T7 | `T7_scope` | scope: index fixes structure, not n_part value | **PASS** |
| T8 | `T8_assessment` | APS_QUARK_PARTITION_INDEX_..._NPART_VALUE_RESIDUAL | **PASS** |

## The APS ξ-invariant `ξ(a) = (η+h)/2 = 1/2 − a`

| a | η(a) = 1−2a | ξ(a) = 1/2 − a |
|---:|---:|---:|
| 0.0 | 1.0 | 0.5 |
| 0.25 | 0.5 | 0.25 |
| 0.3333 | 0.3333 | 0.1667 |
| 0.5 | 0.0 | 0.0 |
| 0.6667 | -0.3333 | -0.1667 |
| 0.75 | -0.5 | -0.25 |
| 1.0 | -1.0 | -0.5 |

As `a : 0 → 1`, `ξ` runs `1/2 → −1/2` (the `n=0` eigenvalue crosses zero): the **spectral flow is `1`, an integer** — the topological APS index. The fractional `ξ` is the η-boundary correction.

## Topological vs residual (the honest split)

- **§8-stable (topological):** the index — the graded doubling `N_q = 2·n_part` (even across all the `quark_axioms` §8 ablations) and the integer spectral flow — the Z₂ mod-2 invariant, protected by APS.
- **Residual (non-topological):** the bare value `n_part` — the continuous, ξ-type (η-boundary) content — drifting `216–255` across §8 (the compensator, PR #97/#107).

So the APS quark partition index **derives the structure** (the graded doubling, the integer index) but **not the value** of `n_part`, exactly as the compensator status requires.

## Verdict

**APS_QUARK_PARTITION_INDEX_TOPOLOGICAL_DOUBLING_STABLE_NPART_VALUE_RESIDUAL.** THE APS QUARK PARTITION INDEX IS A TOPOLOGICAL SPECTRAL-FLOW INTEGER — IT FIXES THE §8-STABLE STRUCTURE OF THE QUARK PARTITION (THE GRADED DOUBLING N_q = 2·n_part), NOT THE VALUE OF n_part. PR #122 assembled the factorized sector sum Z = Σ_{k odd, c₁, n_part} (−1)^k ∫(dL/L) det^{−1/2}_matter e^{i(π/2)(1−2a)} e^{−S_BAM}; this probe reads off its index.

THE INDEX FROM THE Z₂-GRADING. The orientation sign (−1)^k is a Z₂ grading, so I = Tr(−1)^k e^{−βH} is a topological (β-independent) invariant — the orientable (+) and Möbius (−) nonzero modes pair, leaving the net graded content. This is the Witten/APS index of the factorized sum.

THE APS η-BOUNDARY TERM. For the first-order closure operator ∂_τ with holonomy a (eigenvalues 2πi(n+a)/L), the Atiyah–Patodi–Singer ξ-invariant is ξ(a) = (η_A(0)+h)/2 = 1/2 − a = ζ_H(0,a) (twisted, h = 0), continuous in a — the η-boundary correction of PRs #119–#121.

THE INDEX IS AN INTEGER. As the holonomy winds once, a : 0 → 1, exactly one eigenvalue 2π(n+a)/L crosses zero (the n = 0 mode), so the APS index — the spectral flow — is ξ(0⁺) − ξ(1⁻) = 1/2 − (−1/2) = 1, an integer. The fractional ξ is the boundary η-correction; the spectral flow is the integer topological index (the discrete Z₂ content of the sector-phase ledger, PR #121, made into an index).

THE QUARK PARTITION. The quark closure count is N_q = 2·n_part = 466 (n_part = 233, N_lepton = 4 k₅² = 100). The factor of 2 — the EVEN doubling — is precisely the Z₂-graded structure: the orientation index pairs the modes, doubling the count. So the APS index extracts the TOPOLOGICAL (mod-2 / doubling) content of the quark partition.

TOPOLOGICAL vs RESIDUAL. The APS index formalises exactly the split PRs #97/#107 found empirically: the §8-STABLE part is the topological doubling N_q = 2·n_part (even across all the quark_axioms §8 ablations) and the integer spectral flow — the Z₂ mod-2 invariant, protected by APS; the RESIDUAL part is the bare value n_part (the continuous, ξ-type / η-boundary content), which drifts 216–255 across §8, the phenomenological compensator. So the APS quark partition index DERIVES the topological structure (the graded doubling, the integer index) but NOT the value of n_part — exactly as the compensator status requires.
